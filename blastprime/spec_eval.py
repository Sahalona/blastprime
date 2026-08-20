"""引物特异性评估引擎(guide_sup2.md §26-§40)。

主流程(引物设计)使用新引擎,结合位点发现自 2026-08-14 起由 blastn-short
全库批处理给出(k-mer 索引只建模板+侧翼,用户确认):

    blastn-short 批处理(§26 快速淘汰:命中 > prefilter_max_hits → truncated)
    → HSP 3' 端锚定 → BindingSite(§31-§32,含 3' 端错配计数/可延伸判定)
    → in-silico PCR(§33-§35,O(HF log HR))→ 成对特异性最终判定(§38-§40)。

成对判定必须先跑 simulate_pcr,不能单独按 F/R 命中数判定(§44)。
off-target amplicon 具有最高优先级:可扩增脱靶产物直接淘汰(§39)。

sgRNA / 引物分析保留旧 blastn-short 三级函数(运行于 _design_sgrna),
与新引擎并存;两者互不依赖。

Reverse-complement specificity assessment for primer design (rewritten per
guide_sup2.md §26-§40): batched blastn-short binding-site discovery →
HSP 3'-end anchoring → BindingSite index → in-silico PCR → pair-level
classification. The legacy blastn-short three-level functions are retained
for sgRNA mode.
"""

from __future__ import annotations

import tempfile
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from . import blast
from .blast import CancelFlag

# ---------------------------------------------------------------- 等级与评分
# ---------------------------------------------------------------- Pair levels

PAIR_UNIQUE = "unique"
PAIR_SAFE_3P = "safe_with_3p_mismatch"
PAIR_NO_PRODUCT = "no_offtarget_product"
PAIR_OFFTARGET = "offtarget_amplification"
PAIR_KMER_ONLY = "kmer_depth_only"   # skip_spec_eval:未做 blastn-short 逆向验证,仅 k-mer 深度剖面
PAIR_TRUNCATED = "truncated_hits"    # 命中数超上限(§26 快速淘汰):位点未完整评估

PAIR_LABELS = {
    PAIR_UNIQUE: "基因组唯一",
    PAIR_SAFE_3P: "3' 端错配豁免",
    PAIR_NO_PRODUCT: "无脱靶产物",
    PAIR_OFFTARGET: "可扩增脱靶产物(淘汰)",
    PAIR_KMER_ONLY: "基于 k-mer 深度(未做 blastn-short 逆向验证)",
    PAIR_TRUNCATED: "命中数超上限(无法完整评估,淘汰)",
}
PAIR_SCORES = {
    PAIR_UNIQUE: 100,
    PAIR_SAFE_3P: 80,
    PAIR_NO_PRODUCT: 60,
    PAIR_OFFTARGET: 0,
    PAIR_KMER_ONLY: 100,   # 每对实际分 = 3' 端 k-mer 深度剖面折算(见 primer_design)
    # per-pair score is derived from the 3'-end k-mer depth profile (see primer_design)
    PAIR_TRUNCATED: 0,
}

# ---------------------------------------------------------------- BindingSite

@dataclass(frozen=True)
class BindingSite:
    """引物在基因组上的真实可结合位点(§32)。

    strand "+": 引物以自身序列结合正链,其 3' 端在结合窗口下游端(朝右延伸);
    strand "-": rc(引物) 以正链序列出现,引物结合负链,其 3' 端在上游端
    (朝左延伸)。primer_3p 是基因组坐标上的实际 3' 端位置(1-based),
    对 in-silico PCR 极其重要(§32)。

    A real binding site of the primer on the genome (§32). strand "+" means
    the primer itself matches the plus strand with its 3' end at the
    downstream end of the footprint; strand "-" means the reverse complement
    of the primer matches the plus strand, so the primer binds the minus
    strand with its 3' end at the upstream end. primer_3p is the actual
    3'-end position in genomic coordinates (1-based), essential for
    in-silico PCR.
    """
    seq_id: str
    strand: str            # "+" / "-"
    start: int             # 1-based 结合窗口起点(+ 链坐标)
    end: int               # 1-based 结合窗口终点(+ 链坐标)
    primer_3p: int         # 1-based 引物实际 3' 端位置(+ 链坐标)
    extension_capable: bool
    mismatch_total: int
    mismatch_3p: int       # 引物 3' 端最后 3 bp 的错配数(结合取向)
    alignment_score: float  # 百分身份比 (percent identity)

    def to_dict(self) -> dict:
        return {
            "seq_id": self.seq_id, "strand": self.strand,
            "start": self.start, "end": self.end,
            "primer_3p": self.primer_3p,
            "extension_capable": self.extension_capable,
            "mismatch_total": self.mismatch_total,
            "mismatch_3p": self.mismatch_3p,
            "alignment_score": self.alignment_score,
        }


def is_extension_capable(site: BindingSite) -> bool:
    """独立策略函数(§30):3' 端最后 3 bp 完全配对 ⟹ 可延伸。

    Policy function (§30): the last 3 bp at the 3' end fully paired
    ⟹ extension-capable. 修改此策略只改这一处。
    """
    return site.mismatch_3p == 0


# ---------------------------------------------------------------- 局部验证
# ---------------------------------------------------------------- Local validation

def evaluate_binding_site(
    primer: str,
    genome_seq: str,
    pos0: int,
    strand: str,
    seed_k: int,
    seq_id: str,
    min_identity: float = 0.5,
    anchor: str = "3p",
) -> BindingSite | None:
    """在 seed 命中位点执行完整引物局部验证(§27)。

    pos0: seed 起点(0-based)。anchor 决定窗口锚定:
      "3p"  seed 是引物结合取向上 3' 端 seed(默认)——
            strand "+" 时 = 引物自身末 L 碱基(strand "-" 时 = rc(引物) 末 L
            碱基,即 revcomp(引物前 L 碱基));窗口 = [pos0-(n-L), pos0+L)。
      "5p"  seed 是 5' 端 seed(近误位点检测)——窗口 = [pos0, pos0+n)。
            3' 端锚定 seed 匹配 ⟹ 引物末 L 碱基全配 ⟹ mismatch_3p 恒为 0,
            "3' 端错配豁免"级(§38 SAFE_WITH_3P_MISMATCH)在数学上不可达;
            因此另查 5' 端 seed 检测"5' 端配对、3' 端错配"的近误位点
            (旧 blastn-short word_size 7 能发现这些位点)。

    窗口与 footprint(引物 / rc(引物))逐碱基比较;身份比 < min_identity 的
    命中是内部短种子命中(§28),不是结合位点 → None。

    Full-primer local validation at a seed-hit locus (§27). pos0 is the
    0-based start of the seed; anchor decides the window: "3p" (default)
    anchors at the 3'-end seed with window [pos0-(n-L), pos0+L); "5p"
    anchors at the 5'-end seed with window [pos0, pos0+n). The 5'-anchored
    query detects near-miss loci whose 3' end carries mismatches — a
    3'-anchored seed match would always imply mismatch_3p == 0, making the
    SAFE_WITH_3P_MISMATCH level (§38) unreachable; the legacy word_size-7
    blastn found these loci. Hits below the identity threshold are internal
    short-seed hits (§28), not binding sites → None.
    """
    P = primer.upper()
    n = len(P)
    if n < seed_k:
        return None
    footprint = P if strand == "+" else blast.revcomp(P)
    if anchor == "5p":
        start0 = pos0
    else:
        start0 = pos0 - (n - seed_k)
    if start0 < 0 or start0 + n > len(genome_seq):
        return None
    win = genome_seq[start0:start0 + n]
    if len(win) != n:
        return None
    total = 0
    mm3 = 0
    for i in range(n):
        if win[i] != footprint[i]:
            total += 1
            # 3' 端最后 3 bp(结合取向):"+" 在窗口右端,"-" 在窗口左端
            if (strand == "+" and i >= n - 3) or (strand == "-" and i < 3):
                mm3 += 1
    ident = (n - total) / n
    if ident < min_identity:
        return None  # 内部短种子命中 (§28 internal short-seed hit)
    site = BindingSite(
        seq_id=seq_id,
        strand=strand,
        start=start0 + 1,
        end=start0 + n,
        primer_3p=start0 + n if strand == "+" else start0 + 1,
        extension_capable=mm3 == 0,
        mismatch_total=total,
        mismatch_3p=mm3,
        alignment_score=round(100.0 * ident, 1),
    )
    return site


def blast_binding_sites(
    primers: list[tuple[str, str]],
    db_prefix: str,
    params: dict,
    on_log: Callable[[str], None] | None = None,
    cancel: CancelFlag | None = None,
) -> dict[str, list[dict]]:
    """批 blastn-short 全库结合位点发现(用户 2026-08-14 确认:索引只建
    模板+侧翼,结合位点改由 blastn-short 全库给出)。

    全部引物合并为一条 FASTA 单次调用:库只加载/扫描一次,逐引物分发 HSP。
    word_size 7 / E-value 1000 / DUST off(与 sgRNA 通道一致,CHANGE_REPORT
    Q10 文档设计);max_target_seqs 截断 → truncated 语义由有效位点
    (身份过滤后)计数判定,见 hits_to_binding_sites。

    primers = [(引物序列, label)];label 仅限字母数字下划线(写入 FASTA 头,
    经 qseqid 回读)。返回 {label: [HSP dict]}(parse_tabular 结构)。

    Batch blastn-short binding-site discovery against the whole database
    (user-confirmed 2026-08-14: the k-mer index covers only the template +
    flanks, so binding sites are now found by blastn-short). All primers are
    merged into one FASTA for a single call: the database is loaded/scanned
    once and HSPs are distributed per primer.
    """
    if not primers:
        return {}
    tool = blast.require_tool("blastn")
    evalue = params.get("stage4_evalue", 1000)
    max_targets = int(params.get("stage4_max_targets", 1000))
    with tempfile.TemporaryDirectory(prefix="blastprime_bind_") as td:
        qfa = Path(td) / "primers.fa"
        qfa.write_text("".join(
            f">{i}|{label}\n{seq.strip()}\n"
            for i, (seq, label) in enumerate(primers)), encoding="utf-8")
        cmd = [tool, "-query", str(qfa), "-db", db_prefix,
               "-task", "blastn-short", "-word_size", "7",
               # 禁用 DUST:特异性评估必须看到全部脱靶命中(同 sgRNA 通道)
               "-dust", "no",
               "-evalue", str(evalue),
               "-max_target_seqs", str(max_targets),
               "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore qseq sseq",
               "-out", "-"]
        if on_log:
            on_log(f"结合位点 blastn-short: {len(primers)} 条引物 vs {db_prefix}")
        r = blast.run_proc(cmd, on_log=None, cancel=cancel, timeout=600)
        out: dict[str, list[dict]] = {}
        for h in blast.parse_tabular(r.stdout):
            qid = str(h.get("qseqid", ""))
            label = qid.split("|", 1)[1] if "|" in qid else qid
            out.setdefault(label, []).append(h)
    return out


def hits_to_binding_sites(
    primer: str,
    hits: list[dict],
    params: dict,
) -> dict:
    """blastn-short HSP 列表 → BindingSite 集合(§27 语义的 blastn-short 实现)。

    - HSP 覆盖引物 3' 端(plus: qend ≥ n-2;minus: qstart ≤ 3)→ 3p 锚定位点,
      extension_capable = 3' 端最后 3 bp 完全配对(逐位比较 qseq/sseq);
    - HSP 仅覆盖 5' 端 → 5p 近误位点(3' 端不可延伸,§38 SAFE 级可达);
    - 其余为内部短命中(§28 internal),不计位点;
    - identity < binding_min_identity 的 HSP 丢弃(内部短种子命中);
    - 端点修剪(2026-08-14):blastn 在 3' 端错配处停止延伸并把错配碱基修剪
      出比对(qend-qstart+1 > length)→ 修剪掉的碱基视为 3' 端错配,
      extension_capable 随之失效 —— 否则末位错配的脱靶位点会被误判为完美;
    - 次优短 HSP 剔除:引物内部重复会让 blastn 在同一区域反向报 8-12 bp
      次优 HSP(实验确认,blastn-short 对短 query 的已知行为);按 subject
      窗口长度降序处理,新位点窗口被已放置位点窗口包含时剔除(跨链)——
      同一区域的正/反两链位点是同一结合位点的重复报告;
    - hits 总数 > prefilter_max_hits → truncated(§26 快速淘汰,不逐位验证)。

    blastn-short HSP → BindingSite set. HSPs covering the primer 3' end
    anchor 3p sites (extension_capable = last 3 bp fully paired, compared
    position by position on qseq/sseq); HSPs covering only the 5' end are
    near-miss 5p sites; the rest are internal short-seed hits (§28); hits
    below the identity threshold are dropped; 3'-end trimmed HSPs
    (qend-qstart+1 > length, blastn stops at a terminal mismatch and clips
    the mismatched base out of the alignment) count the trimmed bases as
    3'-end mismatches so the site is never misjudged extension-capable;
    overlapping shorter HSPs caused by internal primer repeats are dropped
    (sorted by subject-window length, cross-strand containment test); total
    hits above prefilter_max_hits → truncated (§26 fast elimination).
    """
    P = primer.upper()
    n = len(P)
    max_hits = int(params.get("prefilter_max_hits", 200))
    min_ident = float(params.get("binding_min_identity", 80)) / 100.0
    sites: list[BindingSite] = []
    internal = 0
    # 注意:truncated 判定在身份过滤**之后**(有效位点计数)——原始 HSP 数
    # 在真实大库动辄数百(3' 端 7-mer seed 平均出现数万次,blastn-short
    # 报出的多数是全长身份 <80% 的低 Tm 无效短命中),先截断会把有效
    # 位点寥寥的引物误淘汰(R24 实测:436 个原始 HSP 仅 3 个 ≥80% 身份)。
    # Truncation is judged on the filtered site count, not the raw HSP
    # count — raw hits run into the hundreds in real genomes (3'-end 7-mer
    # seeds occur tens of thousands of times, most HSPs are low-Tm invalid
    # short hits below 80% identity); judging before filtering eliminated
    # primers with only a handful of valid sites.
    # 按 subject 窗口长度降序:先放长位点,后放的短位点被包含时剔除
    hits = sorted(hits, key=lambda h: -(
        abs(int(h.get("send") or 0) - int(h.get("sstart") or 0)) + 1))
    placed: dict[str, list[tuple[int, int]]] = {}   # seq_id → 窗口(跨链)
    for h in hits:
        qs, qe = int(h.get("qstart") or 0), int(h.get("qend") or 0)
        ss, se = int(h.get("sstart") or 0), int(h.get("send") or 0)
        plus = h.get("sstrand", "plus") != "minus"
        if not (qs and qe and ss and se):
            continue
        # 端点修剪:blastn-short 在 3' 端错配处停止延伸并剪掉错配碱基,
        # HSP 的 qend 只覆盖到修剪后的位置(qseq 与 qstart..qend 一致) →
        # 修剪量 = 引物长度 − 实际比对碱基数。
        aln_len = min(abs(qe - qs) + 1, len(str(h.get("qseq", ""))))
        if aln_len == 0:
            internal += 1
            continue
        s_lo, s_hi = min(ss, se), max(ss, se)
        # 3' 端锚定判定与 3' 端错配计数(qseq/sseq 逐位,结合取向下最后 3 bp)
        if plus:
            cov_3p = qe >= n - 2
            cov_5p = qs <= 3
            mm3 = _mm3_tail(h, "plus")
        else:
            cov_3p = qs <= 3
            cov_5p = qe >= n - 2
            mm3 = _mm3_tail(h, "minus")
        # 引物全长身份比:blastn 只延伸种子区,未覆盖端/修剪端都是错配,
        # blastn 报的 pident 只对 HSP 自身长度 → 全长身份比 = 匹配数/引物长
        # (7-mer 随机命中 7/20=35% < 50% → 内部短种子命中 §28,不是位点)。
        # **身份过滤保留**:真实基因组上 3' 端 7-mer seed 平均出现数万次
        # (4^7≈1.6 万 → 460 Mb 库 ≈2.8 万次),max_target_seqs 截断后命中
        # 远超 prefilter_max_hits → truncated 淘汰;中等库中全长身份 <50%
        # 的 7-8 bp 短命中是低 Tm 的无效结合(20 bp 引物只配 8 bp,
        # 退火温度 ~16 °C,体内不会有效延伸),不是真实脱靶位点——保留
        # 它们会把几乎所有引物压到 80 分档,100 分(UNIQUE)不可达。
        # 身份 ≥binding_min_identity(默认 80%)的 3' 端完全配对位点是
        # 有效结合,由修剪端修正(见下)正确标为可延伸 → 60/淘汰。
        # Identity filtering is kept: in real genomes a 3'-end 7-mer seed
        # occurs tens of thousands of times (4^7≈16k → ≈28k in a 460 Mb db),
        # so max_target_seqs truncation pushes hits past prefilter_max_hits
        # → elimination; in medium dbs, full-length identity <50% short hits
        # are low-Tm invalid bindings (an 8 bp footprint on a 20 bp primer,
        # ~16 °C annealing — no effective in-vivo extension), not real
        # off-target sites — keeping them would push nearly every primer to
        # the 80 tier and make 100 unreachable. Fully-paired 3' ends with
        # identity ≥binding_min_identity (default 80%) are valid bindings and
        # are correctly marked extension-capable by the trimming-end fix.
        if (aln_len - int(h.get("mismatch", 0))) / n < min_ident:
            internal += 1
            continue
        if cov_3p:
            anchor = "3p"
        elif cov_5p:
            anchor = "5p"
        else:
            internal += 1   # 纯内部命中(§28)
            continue
        # 修剪掉的碱基计入 3' 端错配,否则末位错配的脱靶位点会被误判
        # 为完美可延伸(2026-08-14 实测:blastn 修剪后 qseq 不含错配位)。
        # **修剪端判定(R23)**:修剪发生在哪端取决于 HSP 起点——
        #   plus 链:qstart==1 时 HSP 从引物起点开始,修剪在 3' 端(blastn
        #     在 3' 端错配处停止延伸,计入 3' 错配);qstart>1 时 seed 在
        #     3' 端、修剪在 5' 端(5' 端错配被剪),不影响 3' 端延伸——
        #     这类"3' 端完全一致"的命中此前被误标为不可延伸,只能判
        #     80(safe_3p)而非其真实档位(可延伸 → 60/淘汰),与用户
        #     "3' 端完全一致被判高分的脱靶位点被漏判"的报告吻合;
        #   minus 链对称:qe==n 时修剪在 3' 端,qe<n 时在 5' 端。
        # Trimming-end resolution: the trimmed bases count as 3'-end
        # mismatches only when the HSP starts at the primer's 3'-end edge —
        # plus: qstart==1 (blastn stopped at a terminal mismatch); minus:
        # qend==n. When the 3'-end seed is intact and the trim is on the
        # 5' end (plus qstart>1 / minus qend<n), the site's 3' end is fully
        # paired and must stay extension-capable — previously mis-marked
        # non-extensible, capping it at 80 instead of its true tier.
        trimmed = max(0, n - aln_len)
        if trimmed:
            trimmed_3p = trimmed if (plus and qs == 1) or (not plus and qe == n) else 0
            mm3 = max(mm3, min(3, trimmed_3p))
        # 3' 端未覆盖(5p/内部)→ 3' 端不可能配对延伸
        site = BindingSite(
            seq_id=str(h.get("sseqid", "")),
            strand="+" if plus else "-",
            start=s_lo, end=s_hi,
            primer_3p=(se if plus else ss),
            extension_capable=cov_3p and mm3 == 0,
            mismatch_total=int(h.get("mismatch", 0)) + trimmed,
            mismatch_3p=mm3 if cov_3p else 3,
            alignment_score=round(float(h.get("pident", 0)), 1),
        )
        # 次优短 HSP 剔除:窗口被已放置(更长)位点包含 → 同一结合位点的
        # 重复报告(引物内部重复触发的 blastn 反向次优 HSP,跨链判定)
        wl = placed.get(site.seq_id)
        if wl and any(s <= site.start and site.end <= e for s, e in wl):
            continue
        placed.setdefault(site.seq_id, []).append((site.start, site.end))
        sites.append(site)
        # 有效位点超上限才截断(身份过滤后仍大量高身份位点 = 真重复引物)
        # Truncate only when valid sites exceed the cap (high-identity sites
        # beyond it mean a genuinely repetitive primer)
        if len(sites) > max_hits:
            return {"sites": [], "internal_hits": internal, "truncated": True}
    return {"sites": sites, "internal_hits": internal, "truncated": False}


def _mm3_tail(h: dict, plus: bool) -> int:
    """HSP 结合取向末 3 bp 错配数:plus = qseq/sseq 末 3 位,minus = 头 3 位。"""
    q = str(h.get("qseq", "")).upper()
    s = str(h.get("sseq", "")).upper()
    k = min(3, len(q), len(s))
    if k <= 0:
        return 3
    rng = range(len(q) - k, len(q)) if plus else range(0, k)
    return sum(1 for i in rng
               if q[i] != s[i] and q[i] != "N" and s[i] != "N")


def find_candidate_binding_sites(
    primer: str,
    db_prefix: str,
    params: dict,
    on_log: Callable[[str], None] | None = None,
    cancel: CancelFlag | None = None,
) -> dict:
    """blastn-short 全库结合位点发现(替换旧的 k-mer seed 索引查询)。

    单引物入口(单引物模式/单元测试):批处理入口在 primer_design 管线里
    按级合并调用,见 blast_binding_sites。

    Full-database binding-site discovery via blastn-short (replacing the old
    k-mer seed index lookups). Single-primer entry (single-primer mode /
    unit tests); the design pipeline batches per level, see
    blast_binding_sites.
    """
    hits = blast_binding_sites([(primer, "p0")], db_prefix, params,
                               on_log, cancel).get("p0", [])
    return hits_to_binding_sites(primer, hits, params)


def merge_binding_sites(*lists: list[BindingSite]) -> list[BindingSite]:
    """合并多个位点列表并去重(键 = seq_id/strand/start/end)。

    目标位点由 blastn-short 直接给出(不需要单独构造),此处只做防御性
    去重:同一引物可能被 3p/5p 两种锚定重复发现,不去重会虚增
    binding_site_count / off_target_sites。
    """
    seen: set[tuple] = set()
    out: list[BindingSite] = []
    for lst in lists:
        for s in lst:
            key = (s.seq_id, s.strand, s.start, s.end)
            if key not in seen:
                seen.add(key)
                out.append(s)
    return out


# ---------------------------------------------------------------- in-silico PCR

def _site_on_target(site: BindingSite, loci: list[dict] | None,
                    buffer: int = 50) -> bool:
    """位点是否落在目标位点(±buffer)内。loci 由 target_loci_from_hsps 产出;
    None(无目标位点信息)→ 一律视为非目标。locus 的 sstart/send 为负链时
    可能 sstart > send,须归一化。

    Whether the site falls within the target loci (±buffer). loci come from
    target_loci_from_hsps; None (no target information) → never on target.
    Negative-strand loci may report sstart > send; normalized here.
    """
    if not loci:
        return False
    return any(l.get("sseqid") == site.seq_id
               and _overlaps((site.start, site.end),
                             (min(l.get("sstart", 0), l.get("send", 0)),
                              max(l.get("sstart", 0), l.get("send", 0))),
                             buffer)
               for l in loci)


def simulate_pcr(
    f_sites: list[BindingSite],
    r_sites: list[BindingSite],
    product_min: int,
    product_max: int,
    target_loci: list[dict] | None = None,
    buffer: int = 50,
    cancel: CancelFlag | None = None,
) -> list[dict]:
    """in-silico PCR(§33-§35):F 正链位点 × R 负链位点,F 在上游、产物长度
    在 [min, max] 内、双方均可延伸 → 构成产物。

    按 seq_id 分组 + bisect 区间查询(§34):O(HF log HR),禁止 HF×HR 双重循环。
    产物 = [F.primer_3p, R.primer_3p](含两端),is_on_target = 两侧位点
    都在目标位点内(预期产物,§36 排除)。

    In-silico PCR (§33-§35): F plus-strand sites × R minus-strand sites,
    F upstream, product length within range, both extension-capable.
    Grouped by seq_id with bisect interval queries (§34): O(HF log HR).
    The product spans [F.primer_3p, R.primer_3p]; is_on_target = both sites
    fall within the target loci (expected product, §36 exclusion).
    """
    if not f_sites or not r_sites:
        return []
    # 契约:第一个参数 = 上游位点(3' 端在右端、向右延伸,strand "+");
    # 第二个参数 = 下游位点(3' 端在左端、向左延伸,strand "-")。
    # evaluate_primer_pair 已按此过滤;此处再校验一次,防错误调用产生假产物。
    f_sites = [s for s in f_sites if s.strand == "+"]
    r_sites = [s for s in r_sites if s.strand == "-"]
    if not f_sites or not r_sites:
        return []
    by_seq: dict[str, list[BindingSite]] = {}
    for s in r_sites:
        by_seq.setdefault(s.seq_id, []).append(s)
    amps: list[dict] = []
    for f in f_sites:
        if cancel is not None and cancel.cancelled:
            raise blast.BlastError("任务已取消")
        rs = by_seq.get(f.seq_id)
        if not rs:
            continue
        rs_sorted = sorted(rs, key=lambda s: s.primer_3p)
        keys = [s.primer_3p for s in rs_sorted]
        lo = f.primer_3p + product_min - 1
        hi = f.primer_3p + product_max - 1
        for r in rs_sorted[bisect_left(keys, lo):bisect_right(keys, hi)]:
            # F 上游 R 已由 lo ≥ f.3p 保证;两 3' 端须可延伸(§35)
            if not (f.extension_capable and r.extension_capable):
                continue
            amps.append({
                "seq_id": f.seq_id,
                "f": f, "r": r,
                "product_len": r.primer_3p - f.primer_3p + 1,
                "is_on_target": (_site_on_target(f, target_loci, buffer)
                                 and _site_on_target(r, target_loci, buffer)),
            })
    return amps


def find_offtarget_pairs(
    f_sites: list[BindingSite],
    r_sites: list[BindingSite],
    product_min: int,
    product_max: int,
    target_loci: list[dict] | None = None,
    buffer: int = 50,
    cancel: CancelFlag | None = None,
) -> list[dict]:
    """simulate_pcr + 过滤掉预期产物。返回脱靶扩增产物列表。"""
    return [a for a in simulate_pcr(f_sites, r_sites, product_min, product_max,
                                    target_loci, buffer, cancel)
            if not a["is_on_target"]]


# ---------------------------------------------------------------- 成对判定
# ---------------------------------------------------------------- Pair classification

def evaluate_primer_pair(
    f_primer: str,
    r_primer: str,
    f_sites: list[BindingSite],
    r_sites: list[BindingSite],
    target_loci: list[dict],
    params: dict,
    cancel: CancelFlag | None = None,
) -> dict:
    """成对特异性最终判定(§38-§40,§44)。

    必须先 simulate_pcr;任何可扩增脱靶产物 → OFFTARGET_AMPLIFICATION(淘汰)。
    两个方向组合都模拟(正常取向 F+ × R-;反向嵌合取向 R+ × F-:
    目标位点反向嵌入基因组时预期产物是后者,两个组合都可能产生真实产物)。

    分级:
      UNIQUE               无非目标 binding site(任何类型)
      SAFE_WITH_3P_MISMATCH 存在非目标 binding,但全部 3' 端不满足延伸条件
      NO_OFFTARGET_PRODUCT  存在可延伸位点,但无法组成可扩增产物对
      OFFTARGET_AMPLIFICATION 存在可延伸、方向正确、长度合理的非目标产物

    Pair-level final classification (§38-§40, §44): must run simulate_pcr;
    any amplifiable off-target product → OFFTARGET_AMPLIFICATION (rejected).
    Both orientation combinations are simulated (normal F+ × R-; inverted
    R+ × F-, which is the expected product when the target locus is embedded
    on the reverse strand) — either combination can yield real products.
    """
    buffer = int(params.get("target_buffer", 50))
    omin = int(params.get("offtarget_product_min", 50))
    omax = int(params.get("offtarget_product_max", 4000))
    f_sites = merge_binding_sites(f_sites)  # 防调用方传入重复位点(防御)
    r_sites = merge_binding_sites(r_sites)
    f_plus = [s for s in f_sites if s.strand == "+"]
    f_minus = [s for s in f_sites if s.strand == "-"]
    r_plus = [s for s in r_sites if s.strand == "+"]
    r_minus = [s for s in r_sites if s.strand == "-"]
    amps = []
    amps += simulate_pcr(f_plus, r_minus, omin, omax, target_loci, buffer, cancel)
    amps += simulate_pcr(r_plus, f_minus, omin, omax, target_loci, buffer, cancel)
    off_amps = [a for a in amps if not a["is_on_target"]]

    # 所有位点(含反向角色类型)中不在目标位点内的部分
    # Sites (both role types) that fall outside the target loci
    off_sites = [s for s in f_sites + r_sites
                 if not _site_on_target(s, target_loci, buffer)]
    out = {
        "level": PAIR_OFFTARGET,
        "label": PAIR_LABELS[PAIR_OFFTARGET],
        "spec_score": PAIR_SCORES[PAIR_OFFTARGET],
        "note": "",
        "off_target_sites": len(off_sites),
        "binding_site_count": len(f_sites) + len(r_sites),
        "offtarget_amplicons": len(off_amps),
        "amplifiable_pairs": [_amp_to_dict(a) for a in off_amps],
    }
    if off_amps:
        # 结构化 note_key/note_params:前端语言切换时按当前语言渲染
        # Structured note_key/note_params: the frontend re-renders per language
        out["note"] = (f"{len(off_amps)} 个可扩增脱靶产物"
                       f"(产物 {min(a['product_len'] for a in off_amps)}~"
                       f"{max(a['product_len'] for a in off_amps)} bp),已淘汰")
        out["note_key"] = "offtarget_amps"
        out["note_params"] = {"n": len(off_amps),
                              "lo": min(a["product_len"] for a in off_amps),
                              "hi": max(a["product_len"] for a in off_amps)}
        return out

    if not off_sites:
        level = PAIR_UNIQUE
        note = "两条引物无非目标结合位点"
        note_key, note_params = "pair_unique", None
    elif all(not s.extension_capable for s in off_sites):
        level = PAIR_SAFE_3P
        note = f"存在 {len(off_sites)} 个非目标结合位点,但 3' 端均不可延伸"
        note_key, note_params = "safe_3p", {"n": len(off_sites)}
    else:
        level = PAIR_NO_PRODUCT
        note = (f"存在 {len(off_sites)} 个非目标结合位点,"
                f"但无法构成可扩增产物对(共 {len(amps)} 个预期/无效组合)")
        note_key, note_params = "no_product", {"n": len(off_sites), "m": len(amps)}
    out.update({
        "level": level, "label": PAIR_LABELS[level],
        "spec_score": PAIR_SCORES[level], "note": note,
        "note_key": note_key, "note_params": note_params,
    })
    return out


def classify_single_primer(
    sites: list[BindingSite],
    target_loci: list[dict],
    params: dict,
) -> dict:
    """单引物模式特异性(新引擎):不进行成对扩增判定。

    UNIQUE = 无非目标位点;SAFE_WITH_3P_MISMATCH = 全部非目标位点 3' 端
    不可延伸;否则为存在可延伸脱靶位点(不可用)。

    Single-primer mode specificity (new engine): no pair amplification check.
    """
    buffer = int(params.get("target_buffer", 50))
    off = [s for s in sites if not _site_on_target(s, target_loci, buffer)]
    if not off:
        return {"level": PAIR_UNIQUE, "label": PAIR_LABELS[PAIR_UNIQUE],
                "spec_score": PAIR_SCORES[PAIR_UNIQUE], "off_target_sites": 0,
                "note": "无非目标结合位点",
                "note_key": "single_unique", "note_params": None}
    if all(not s.extension_capable for s in off):
        return {"level": PAIR_SAFE_3P, "label": PAIR_LABELS[PAIR_SAFE_3P],
                "spec_score": PAIR_SCORES[PAIR_SAFE_3P],
                "off_target_sites": len(off),
                "note": "存在非目标结合位点,但 3' 端均不可延伸",
                "note_key": "single_safe_3p", "note_params": None}
    ext = [s for s in off if s.extension_capable]
    return {"level": PAIR_NO_PRODUCT, "label": PAIR_LABELS[PAIR_NO_PRODUCT],
            "spec_score": PAIR_SCORES[PAIR_NO_PRODUCT],
            "off_target_sites": len(off),
            "note": f"存在 {len(ext)} 个可延伸脱靶位点",
            "note_key": "single_ext", "note_params": {"n": len(ext)}}


def _amp_to_dict(a: dict) -> dict:
    return {
        "seq_id": a["seq_id"],
        "product_len": a["product_len"],
        "f": a["f"].to_dict(),
        "r": a["r"].to_dict(),
    }


def _overlaps(a: tuple[int, int], b: tuple[int, int], buffer: int = 0) -> bool:
    return not (a[1] + buffer < b[0] or b[1] + buffer < a[0])


# ================================================================ sgRNA 保留区
# ================================================================ Legacy (sgRNA)

L1_UNIQUE = "unique"
L2_EXEMPT = "exempt_3p"
L3_UNPAIRED = "unpaired"

SPEC_LABELS = {
    L1_UNIQUE: "基因组唯一匹配",
    L2_EXEMPT: "3' 端错配豁免",
    L3_UNPAIRED: "不可成对扩增",
    "eliminated": "已淘汰",
}

SPEC_SCORES = {L1_UNIQUE: 100, L2_EXEMPT: 80, L3_UNPAIRED: 60}


def _three_prime_mismatch(qseq: str, sseq: str, window: int = 3) -> int | None:
    """检查比对序列 3' 端(window bp)内是否存在错配。返回第一个错配位置偏移,无错配返回 None。

    Check whether the alignment has mismatches within the last window bp at
    the 3' end. Returns the offset of the first mismatch, or None if none.
    """
    q, s = qseq.upper(), sseq.upper()
    n = min(len(q), len(s), window)
    for i in range(len(q) - n, len(q)):
        if i < len(s) and q[i] != s[i] and q[i] != "N" and s[i] != "N":
            return len(q) - i
    return None


def _binding_hits(hits: list[dict], primer_len: int | None) -> list[dict]:
    """仅保留比对延伸至引物 3' 端(qend == 引物全长)的命中。

    短种子命中(word_size 7 触发、只对齐引物内部片段)的引物 3' 端并未与
    基因组配对,不可能作为引物结合位点,不参与特异性判定(仅留作详情展示)。
    """
    if not primer_len:
        return hits
    return [h for h in hits if h.get("qend") == primer_len]


def _extension_capable(h: dict, primer_len: int | None = None, window: int = 3) -> bool:
    """该命中是否可延伸:比对须延伸至引物 3' 端,且 3' 端最后 window bp 完全配对。"""
    if primer_len is not None and h.get("qend") != primer_len:
        return False
    return _three_prime_mismatch(h.get("qseq", ""), h.get("sseq", ""), window) is None


def target_loci_from_hsps(hsps: list[dict], template_len: int,
                          target: tuple[int, int] | None,
                          buffer: int = 50) -> list[dict]:
    """从第一步全序列比对中,找出模板"主位点"上覆盖目标区域(±buffer)的库位点。

    主位点 = 按 subject 区间重叠分组后 bitscore 总和最高的组(模板的基因组家园);
    同源重复拷贝(间隔不相邻的其他位点)不并入目标位点,其命中作为脱靶位点
    参与三级排序与淘汰 —— 纯重复模板因此整体淘汰并给出失败诊断。
    相邻串联重复(全长 HSP 覆盖多个拷贝)属于同一主位点,全部命中视为目标命中。
    """
    if target:
        t1, t2 = target
    else:
        t1, t2 = 1, template_len
    cand: list[dict] = []
    for h in hsps:
        qs, qe = h.get("qstart"), h.get("qend")
        if qs is None:
            continue
        if qs > qe:
            qs, qe = qe, qs
        if qe < t1 - buffer or qs > t2 + buffer:
            continue
        cand.append(h)
    if not cand:
        return []
    # 按 (sseqid, 库区间重叠) 分组:不同库条目绝不合并 —— 否则 chrC 的命中会
    # 并入 chrB 的串联重复组,合并组 bitscore 总和虚高,把真正的基因组家园
    # (chrA)挤出主位点,目标位点缺失导致预期结合位点全部漏报。
    # Group by (sseqid, subject-interval overlap): HSPs from different
    # subjects must never merge — otherwise a chrC hit folds into chrB's
    # tandem group, the merged group's bitscore sum inflates, the true
    # genome home (chrA) is dropped from the target loci, and all expected
    # binding sites are lost.
    groups: list[dict] = []
    for h in sorted(cand, key=lambda x: (x.get("sseqid", ""), x.get("sstart", 0), x.get("send", 0))):
        ss, se = h.get("sstart"), h.get("send")
        if ss > se:
            ss, se = se, ss
        sid = h.get("sseqid", "")
        placed = False
        for g in groups:
            if g["seq_id"] != sid:
                continue
            if g["se"] >= ss and se >= g["ss"]:
                g["ss"], g["se"] = min(g["ss"], ss), max(g["se"], se)
                g["bits"] += float(h.get("bitscore", 0) or 0)
                g["hsps"].append(h)
                placed = True
                break
        if not placed:
            groups.append({"seq_id": sid, "ss": ss, "se": se,
                           "bits": float(h.get("bitscore", 0) or 0), "hsps": [h]})
    best = max(groups, key=lambda g: g["bits"])
    # 保留 qstart/qend/sstrand:模板↔数据库映射(§37 预期位点构造)需要
    return [{"sseqid": h.get("sseqid"),
             "sstart": int(h.get("sstart", 0)), "send": int(h.get("send", 0)),
             "qstart": int(h.get("qstart", 0)), "qend": int(h.get("qend", 0)),
             "sstrand": h.get("sstrand", "plus")}
            for h in best["hsps"]]


def _is_on_target(h: dict, loci: list[dict], buffer: int = 50) -> bool:
    ss, se = h.get("sstart"), h.get("send")
    if ss is None:
        return False
    if ss > se:
        ss, se = se, ss
    return any(l["sseqid"] == h.get("sseqid")
               and _overlaps((ss, se), (l["sstart"], l["send"]), buffer)
               for l in loci)


def _can_form_product(f_hit: dict, r_hit: dict, min_len: int, max_len: int) -> bool:
    """两条引物命中是否可形成 PCR 产物(sgRNA 用旧实现)。"""
    if f_hit.get("sseqid") != r_hit.get("sseqid"):
        return False
    f_s, f_e = f_hit.get("sstart"), f_hit.get("send")
    r_s, r_e = r_hit.get("sstart"), r_hit.get("send")
    if f_s is None or r_s is None:
        return False
    f_plus = f_hit.get("sstrand", "plus") != "minus"
    r_minus = r_hit.get("sstrand", "plus") == "minus"
    if not (f_plus and r_minus):
        return False
    f_3p = max(f_s, f_e)
    r_3p = min(r_s, r_e)
    if r_3p <= f_3p:
        return False
    prod = r_3p - f_3p + 1
    return min_len <= prod <= max_len


def evaluate_single_primer(
    hits: list[dict], target_loci: list[dict], params: dict,
    primer_len: int | None = None,
) -> dict:
    """单引物(sgRNA / 单引物模式)3' 端法则评估(旧实现)。

    L1: 无脱靶命中;L2: 全部脱靶命中 3' 端 1~3 bp 存在错配;
    其余: 存在可延伸(3' 端完全配对)的脱靶命中,记录脱靶位点数。
    """
    buffer = int(params.get("target_buffer", 50))
    bind = _binding_hits(hits, primer_len)
    off = [h for h in bind if not _is_on_target(h, target_loci, buffer)]
    if not off:
        return {"level": L1_UNIQUE, "off_target_sites": 0}
    if all(_three_prime_mismatch(h.get("qseq", ""), h.get("sseq", "")) is not None
           for h in off):
        return {"level": L2_EXEMPT, "off_target_sites": 0}
    sites = len({h["sseqid"] for h in off if _extension_capable(h, primer_len)})
    return {"level": "else", "off_target_sites": sites}


def run_specificity_blast(
    primers: list[dict],
    db_prefix: str,
    params: dict,
    on_log: Callable[[str], None] | None = None,
    cancel: CancelFlag | None = None,
) -> None:
    """对引物列表执行 blastn-short 并回填 hits(sgRNA 用旧实现)。"""
    tool = blast.require_tool("blastn")
    evalue = params.get("stage4_evalue", 1000)
    with tempfile.TemporaryDirectory(prefix="blastprime_spec_") as td:
        td = Path(td)
        idx = 0
        for pr in primers:
            for side in ("left", "right"):
                if side not in pr or not pr[side].get("seq"):
                    continue
                idx += 1
                seq = pr[side]["seq"]
                pr[side]["primer_len"] = len(seq)
                qfa = td / f"p{idx}.fa"
                qfa.write_text(f">p{idx}\n{seq}\n", encoding="utf-8")
                cmd = [tool, "-query", str(qfa), "-db", db_prefix,
                       "-task", "blastn-short", "-word_size", "7",
                       # 禁用 DUST:特异性评估必须看到全部脱靶命中,
                       # DUST 会屏蔽简单重复引物的真实脱靶,导致特异性被高估
                       "-dust", "no",
                       "-evalue", str(evalue),
                       "-max_target_seqs", str(params.get("stage4_max_targets", 1000)),
                       "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore qseq sseq",
                       "-out", "-"]
                if on_log:
                    on_log(f"特异性 BLAST: {seq} (blastn-short)")
                r = blast.run_proc(cmd, on_log=None, cancel=cancel, timeout=300)
                pr[side]["hits"] = blast.parse_tabular(r.stdout)
