"""引物设计管线(guide_sup2.md 新算法:k-mer 计数深度 + binding-site + in-silico PCR)。

流程(每条模板查询):
  1. 整段模板 blastn + 纯 Python k-mer 计数窗口查询(第一步,8/10/12/15 bp
     四种 k-mer 长度,逐碱基窗口)→ step1_hsps + target_loci(模板↔库映射,§37)
     + 逐碱基匹配深度
  2. 逐碱基特异性剖面 global / 3'(§16-§19):深度 d → 1/(1+log2(d)),
     与模板内自重复 k-mer 分(补充串联重复等)逐碱基取 min
  3. 分级设计(§23):specificity_levels 阈值逐级放宽
     - 允许区 = global ≥ g_th 且 three_prime ≥ t_th(补集作 primer3 excluded)
     - primer3 物理设计(§24-§25,职责边界 §45)
     - 3' seed 深度预筛(§26,快速淘汰明显不特异候选)
     - binding-site 局部验证(§27-§32:索引位点 + 预期目标位点合并,§37)
     - pair-level in-silico PCR 判定(§33-§40,evaluate_primer_pair 必须模拟 PCR)
     - 综合评分(§39 可配置权重)+ 排序(§40 tuple)→ 成功即停
  4. 全级别失败 → 失败诊断(§42,failure_stage + 各级候选统计 + 低特异区列表)

附加模式:sgRNA(§58 四段式:参与 k-mer 评分,按剖面阈值逐级过滤 guide 3' 端位置)、单引物(不做成对扩增判定,§59)。

Primer design pipeline (guide_sup2.md: k-mer-count depth specificity +
binding-site analysis + in-silico PCR).
Per template record:
  1. whole-template blastn + native k-mer counting (step 1) → step1_hsps +
     per-base depth
  2. specificity profiles (global/3') from depth, min'ed with the in-template
     self-repeat k-mer score
  3. tiered design (§23): specificity_levels thresholds, relaxed progressively
  4. failure diagnosis (§42) when all levels fail
Additional modes: sgRNA (legacy kept, §58), single primer (no pair PCR, §59).
"""

from __future__ import annotations

import hashlib
import io
from types import SimpleNamespace
from typing import Callable

from Bio import Seq, SeqIO

from . import blast
from . import kmer_count
from .blast import CancelFlag
from .masking import compute_depth, compute_windowed_depth_compact
from .primer_index import KmerIndex, compute_profiles, kmer_occurrence_score
from .primer_metrics import (
    composite_score, design_primers, dimer_stats, gc, gc_clamp_3p, physical_score, tm,
)
from .spec_eval import (
    PAIR_KMER_ONLY, PAIR_LABELS, PAIR_NO_PRODUCT, PAIR_OFFTARGET,
    PAIR_SCORES, PAIR_TRUNCATED, SPEC_LABELS, SPEC_SCORES,
    blast_binding_sites, classify_single_primer, evaluate_primer_pair,
    evaluate_single_primer, hits_to_binding_sites, run_specificity_blast,
    target_loci_from_hsps,
)

# count 分档(R33):阈值按得分函数 count^(-1/3) 换算——count=1 → 1.0、
# count=2/3 → ≥3^(-1/3)=0.6934、count=4/5/6 → ≥6^(-1/3)=0.5503、
# count≥7 → L4 全放行。L1/L2 在 global/3' 双维度按最长尺度(15-mer)
# 的 count 分档(global/three_prime 均为 max 合并,分档取其"最佳证据")。
# L2/L3 的 global/3' 阈值可由高级参数"阶段二/三放行阈值"覆盖
# (level2_global_th/level3_global_th,见 _levels_from_params)。
DEFAULT_SPECIFICITY_LEVELS = [[1.0, 1.0], [0.6934, 0.6934], [0.5503, 0.5503], None]
LEVEL_NAMES = {
    1: "Level 1(count=1,唯一):global≥1.0,3'≥1.0",
    2: "Level 2(count 2-3):global≥0.6934,3'≥0.6934",
    3: "Level 3(count 4-6):global≥0.5503,3'≥0.5503",
    4: "Level 4(count≥7 或无预屏蔽)",
}


def _levels_from_params(params: dict) -> list:
    """按 count 分档 + 可配置的 L2/L3 放行阈值生成四段式 levels。

    L1 固定 (1.0, 1.0)(count=1 唯一);L2/L3 的 global 与 3' 阈值取同一
    参数值(高级参数"阶段二/三放行阈值",默认 3^(-1/3)/6^(-1/3))——
    参数化后可视化色带阈值随结果 levels.thresholds 下发,与设计一致。
    """
    t2 = float(params.get("level2_global_th", 3 ** (-1.0 / 3.0)))
    t3 = float(params.get("level3_global_th", 6 ** (-1.0 / 3.0)))
    return [[1.0, 1.0], [t2, t2], [t3, t3], None]
# 滑动窗口尺度(R11/R12):四种 k-mer 长度逐碱基窗口(步长 = 1)。R11 起
# 每尺度只统计"全长精确命中";R12 起窗口查询由纯 Python k-mer 计数
# (kmer_count)直接统计,构造上即全长精确匹配。
# global = 各尺度分数 max 合并(combine_kmer_scores 语义:最长 k-mer 占优,
# 0.95 阈值在真实基因组上仍可达);three_prime = 各尺度以该位置结尾的
# k-mer 分数取 max(与自重复分量 compute_profiles 的 3' 窗口约定一致)。
# Sliding-window scales (R11/R12): four k-mer lengths, per-base windows
# (step = 1). Since R11 each scale counts only full-length EXACT hits;
# since R12 the window queries are counted natively (kmer_count), which is
# exact by construction.
# global = max-combined scores across scales (combine_kmer_scores
# semantics: the longest k-mer dominates, keeping the 0.95 thresholds
# reachable on real genomes); three_prime = max over the scales of the
# k-mer ending at the position (same 3'-window convention as the
# self-repeat component's compute_profiles).
WINDOW_KMERS = (8, 10, 12, 15)
# 失败阶段(§42)(Failure stages)
FAIL_NO_HIGH_SPEC = "TARGET_NO_HIGH_SPECIFICITY_REGION"
FAIL_PRIMER3 = "PRIMER3_NO_CANDIDATE"
FAIL_SINGLE_SPEC = "PRIMER_SINGLE_SPECIFICITY_FAILED"
FAIL_PAIR_OFFTARGET = "PRIMER_PAIR_OFFTARGET"
FAIL_NO_PAIR = "NO_ACCEPTABLE_PAIR"

# 查询内进度映射(§48 建议,可按级细分;总和 1.0)
# Query-local progress mapping (per §48 recommendations, sub-divided per level; sums to 1.0)
_F_PROFILE = 0.25           # 剖面(10-30 段)
_F_STEP1_BLAST = 0.05       # 第一步整段 blastn 完成(剖面内 0-5)
_F_KMER_S = 0.05            # k-mer 计数起点(剖面内 5)
_F_KMER_SPAN = 0.18         # k-mer 计数跨度(剖面内 5-23,按尺度推进)
_F_SELF_S = 0.92            # 自重复剖面起点(剖面内 23-25 段)
_F_LEVEL_START = 0.25       # 分级设计起点
_F_LEVEL_SPAN = 0.70        # 4 级共占
_F_RANK = 0.95


# ---------------------------------------------------------------- 主入口 (Main entry)

def design_pipeline(
    template_text: str,
    db_prefix: str,
    target: tuple[int, int] | None,
    params: dict,
    on_log: Callable[[str], None] | None = None,
    cancel: CancelFlag | None = None,
    locate_ctx: dict | None = None,
    spec_db: str | None = None,
    on_progress: Callable[[float], None] | None = None,
    kmer_caches: list[dict] | None = None,
) -> dict:
    """对模板 FASTA(可多条)执行设计,返回按查询分组的结果。

    locate_ctx: 定位模式上下文(guide_sup1.md §5,由 /api/design/extract 产出),
    非空时对每条结果做基因组绝对坐标标注与产物覆盖标记。
    spec_db: 特异性比对库(名称型 targetbase= 解析结果);None = 使用
    db_prefix 作为特异性比对库。
    on_progress: 任务级进度 0~100(索引 0-10,各查询平分 10-100)。

    Run design on the template FASTA (possibly multiple records), returning
    results grouped by query. spec_db: the specificity-comparison db
    (resolved name-style targetbase=); None = db_prefix. on_progress:
    job-level progress 0-100 (index 0-10; queries share 10-100).
    """
    mode = params.get("mode", "standard")
    # `#` 注释行自动忽略 + 重名标题去重(提交/解析时改写,输入框原文不动)
    # Lines starting with `#` are auto-ignored; duplicate headers are deduped
    # (submit/parse time only, the input box keeps the original text)
    template_text = blast.dedupe_fasta_headers(blast.strip_comment_lines(template_text))
    stripped = template_text.lstrip()
    if stripped.startswith(">"):
        recs = list(SeqIO.parse(io.StringIO(template_text), "fasta"))
    else:
        # 无 > 头的纯序列
        # Bare sequence without a ">" header
        seq = "".join(template_text.split())
        if not seq:
            raise ValueError("模板序列为空")
        recs = [SeqIO.SeqRecord(seq=Seq.Seq(seq), id="query", description="query")]

    # 自重复 k-mer 索引(仅模板,构建极快):窗口 blastn 深度把"同一条目内
    # subject 区间重叠的相邻拷贝"合并为同一位点 → 串联重复/微卫星显示为
    # 单拷贝;此索引只为捕获模板内部的自重复,与深度剖面逐碱基取 min
    # (见 _design_with_profile)。基因组特异性的主信号是第一步的窗口
    # blastn 深度(全库),不再是索引本身——索引构建失败只降级、不阻断。
    # This template-only k-mer index serves only as a self-repeat component:
    # windowed blastn depth merges adjacent in-entry copies (overlapping
    # subject ranges) into one locus, hiding tandem repeats; the index
    # recovers those. Genome-wide specificity comes from the step-1 window
    # depth instead, so index failure only degrades (never aborts) a design.
    index = None
    if mode in ("standard", "single"):
        try:
            index = KmerIndex(
                sequences=[str(r.seq).upper() for r in recs],
                seq_ids=[r.id for r in recs],
                use_cache=False,
                on_log=on_log,
                cancel=cancel,
            )
        except Exception as e:
            index = None
            if on_log:
                on_log(f"自重复 k-mer 索引构建失败(忽略,深度剖面不受影响): {e}")

    results = []
    total = len(recs)
    for qi, rec in enumerate(recs):
        if cancel is not None and cancel.cancelled:
            break
        if on_log:
            on_log(f"[{qi + 1}/{total}] 设计查询: {rec.description[:80]}")
        q_progress = (lambda f: on_progress(100.0 * (qi + f) / total)
                      if on_progress else None)
        try:
            res = _design_one(str(rec.seq).upper(), rec.description, db_prefix,
                              target, params, on_log, cancel, locate_ctx,
                              index=index, spec_db=spec_db, on_progress=q_progress,
                              kmer_caches=kmer_caches)
            results.append(res)
        except Exception as e:
            if cancel is not None and cancel.cancelled:
                raise
            results.append({
                "query": rec.description, "success": False, "error": str(e),
                "stage_reached": 0, "template_len": len(rec.seq),
            })
            if on_log:
                on_log(f"查询 {rec.description} 失败: {e}")
    ok = sum(1 for r in results if r.get("success"))
    return {"results": results, "total": total, "succeeded": ok,
            "mode": mode, "params": params}


def _kmer_cache_key(seq: str, depth_db: str, params: dict) -> str:
    """k-mer 结果缓存键:模板序列 + 特异性库 + 尺度 + 模式。
    F/R 设计范围不再入键——范围变化只改变"检测范围",命中后由
    `_find_kmer_cache` 按缓存记录的 detected_ranges 检查当前设计范围是否
    被覆盖:范围外位置从未计数(R15),未覆盖则视为未命中重算,覆盖则复用
    (同序列改范围可免重算,前提是新范围不超出缓存检测范围)。模式归一化
    为 "sgrna" vs "pair"(standard/single 的剖面结构完全相同,共享缓存;
    R34:标准→单引物切换应复用)。sgRNA 独立是因为其设计路径特殊。

    Cache key for k-mer results: template sequence + specificity db + scales
    + mode. F/R design ranges are deliberately not in the key — a range
    change only alters the detected region, so `_find_kmer_cache` checks the
    cached detected_ranges against the current design ranges on a hit.
    Mode is normalised to "sgrna" vs "pair": standard and single share
    identical profile structures and reuse each other's caches (R34);
    sgRNA stays separate due to its distinct design path.
    """
    mode = "sgrna" if params.get("mode", "standard") == "sgrna" else "pair"
    raw = f"{seq}\0{depth_db}\0{mode}\0{WINDOW_KMERS!r}"
    return hashlib.sha256(raw.encode()).hexdigest()


def _detected_ranges(designable: set[int] | None, n: int) -> list[list[int]]:
    """k-mer 实际检测位置(0-based designable 集合)→ 1-based 连续区间列表;
    designable=None(无设计范围,全检测)→ [[1, n]]。缓存随结果携带该元数据,
    供复用前判断"当前设计范围是否已被检测覆盖"。

    Converts the detected positions (0-based designable set) into 1-based
    contiguous ranges; None (no design ranges, full detection) -> [[1, n]].
    The cache carries this metadata so reuse can verify that the current
    design ranges were actually covered by the detection.
    """
    if designable is None:
        return [[1, n]]
    pos = sorted(designable)
    out: list[list[int]] = []
    s = prev = pos[0]
    for p in pos[1:]:
        if p == prev + 1:
            prev = p
        else:
            out.append([s + 1, prev + 1])
            s = prev = p
    out.append([s + 1, prev + 1])
    return out


def _covered_by(ranges: list[list[int]], designable: set[int]) -> bool:
    """当前设计的设计位置(0-based)是否全部落在缓存检测区间(1-based)内。"""
    for p in designable:
        if not any(s - 1 <= p <= e - 1 for s, e in ranges):
            return False
    return True


def _merge_ranges(a: list[list[int]], b: list[list[int]]) -> list[list[int]]:
    """两个 1-based 区间列表([start, end] 闭区间)合并为并集(排序去重)。"""
    seg = sorted((s, e) for s, e in a + b if s <= e)
    out: list[list[int]] = []
    for s, e in seg:
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def _uncovered_positions(ranges: list[list[int]],
                         designable: set[int]) -> set[int]:
    """designable(0-based)中未被缓存检测区间(1-based)覆盖的位置集合。
    全命中 → 空集;部分重叠 → 未覆盖位置(需重新计数);designable=None
    (无设计范围,全检测)时以模板全长计,返回全 0..n-1(缓存区间按 [[1,n]]
    语义比较)。

    Positions in designable (0-based) not covered by the cached detected
    ranges (1-based). Empty set = full hit; non-empty = partial overlap
    (those positions must be re-counted); designable=None (no ranges, full
    detection) compares against the full template.
    """
    if designable is None:
        # 全检测:未覆盖 = 缓存检测区间补集(按 0-based 位置)
        n_max = max((e for _, e in ranges), default=0)
        return {i for i in range(n_max) if not any(s - 1 <= i <= e - 1 for s, e in ranges)}
    return {p for p in designable
            if not any(s - 1 <= p <= e - 1 for s, e in ranges)}


def _find_kmer_cache(kmer_caches: list[dict] | None, key: str,
                     n: int, designable: set[int] | None = None,
                     seq: str = "", depth_db: str = "",
                     mode: str = "standard"
                     ) -> tuple[dict | None, set[int] | None, int]:
    """按键 + 模板长度取缓存,返回 (cache, uncovered, offset):
    - cache 为 None:未命中(键/长度不符),整体重算;
    - uncovered 为空集:缓存检测范围覆盖当前设计范围,直接复用;
    - uncovered 非空:部分重叠——重叠区复用缓存 depth/profiles,未覆盖
      位置重新计数后合并(R28:范围变化不再整体重算);
    - offset:当前模板在缓存模板序列中的 0-based 位置(子串命中时 >0,
      精确命中为 0)——缓存数组/检测范围按 offset 映射到当前坐标。
    子串命中(R29):侧翼长度变化导致模板重新提取(嵌套关系,当前模板是
    缓存模板的子串)时,键不匹配但重叠区仍可复用——按 offset 映射
    检测范围与数据。缺省/损坏条目安全忽略(重算兜底)。

    Pick the cache entry by key + template length, returning (cache,
    uncovered, offset): None cache = miss; empty uncovered = full coverage
    (reuse as-is); non-empty uncovered = partial overlap — reuse the cached
    depth/profiles on the overlap and re-count the uncovered positions.
    offset locates the current template within the cached template's
    sequence (0 for an exact-key hit, >0 for a substring hit) so arrays and
    detected ranges map onto current coordinates. Substring hits (R29)
    cover flank-length changes, which re-extract a nested template whose
    key no longer matches but whose overlap is still reusable.
    """
    if not kmer_caches:
        return None, None, 0
    for c in kmer_caches:
        if c.get("key") != key or c.get("template_len") != n:
            continue
        ranges = c.get("detected_ranges")
        if designable is not None and ranges:
            uncovered = _uncovered_positions(ranges, designable)
            return c, uncovered, 0
        return c, set(), 0
    # 子串命中:当前模板是缓存模板的子串(flank 变化),要求同库同模式
    # (模式归一化:standard/single 同为 "pair")
    # Substring hit: the current template is a substring of a cached one
    # (flank change); requires the same db and mode (normalised: standard
    # and single both read "pair")
    if seq and depth_db and mode != "sgrna":
        mode_norm = "pair"
        for c in kmer_caches:
            if c.get("db") != depth_db or c.get("mode") != mode_norm:
                continue
            cseq = c.get("seq") or ""
            off = cseq.find(seq)
            if off < 0:
                continue
            ranges = [[max(1, s - off), e - off]
                      for s, e in (c.get("detected_ranges") or [])]
            ranges = [r for r in ranges if r[0] <= r[1]]
            if designable is not None and ranges:
                uncovered = _uncovered_positions(ranges, designable)
                return c, uncovered, off
            return c, set(), off
    return None, None, 0


def _design_one(
    seq: str, query_name: str, db_prefix: str,
    target: tuple[int, int] | None, params: dict,
    on_log: Callable[[str], None] | None, cancel: CancelFlag | None,
    locate_ctx: dict | None = None,
    index: KmerIndex | None = None,
    spec_db: str | None = None,
    on_progress: Callable[[float], None] | None = None,
    kmer_caches: list[dict] | None = None,
) -> dict:
    n = len(seq)
    if n < 40:
        raise ValueError("模板序列过短(<40 bp)")
    mode = params.get("mode", "standard")

    # 第一步:整段模板 blastn(R12 起窗口查询不再走 blastn——窗口 k-mer 的
    # "出现次数"由 kmer_count 纯 Python 精确计数,blastn 只查整段)
    # R8:整段 blastn 对每个库条目只报最优 HSP 并沿查询全长延伸,模板内部
    # 与跨位点的重复区会被单条主 HSP 掩盖 → 逐碱基深度恒为 1,剖面全绿
    # ("全部可用+单拷贝")。R11 起用滑动窗口各自独立查询使每个拷贝都能命中,
    # 窗口位点计数即该处的基因组拷贝数;R12 把窗口查询替换为直接统计
    # k-mer 出现次数(blastn 对窗口查询的 seed 延伸/对齐产物被 R11 的全长
    # 精确过滤全部丢弃,计数程序构造上即"全长精确匹配",语义逐位一致)。
    # Step 1: whole-template blastn (since R12 the window queries no longer
    # go through blastn — k-mer occurrence counts come from the native
    # kmer_count module; blastn only handles the whole template).
    # R8: a whole-template query only reports the best HSP per entry (extended
    # across the full query), hiding in-template and cross-locus repeats, so
    # the per-base depth stays 1 ("all usable + single-copy"). R11 added
    # independent sliding windows so each copy hits; R12 replaces the window
    # blastn with direct k-mer counting (the seed-extension/alignment outputs
    # blastn produced for window queries were discarded entirely by R11's
    # full-length exact filter, so the counter is exactly that semantics).
    # F/R 设计范围外区域(R15):正、反向引物 3' 端均不可放置的区域
    # (left_range ∪ right_range 之外)不进行 k-mer 计数/评分 —— 不可用
    # 区域的重复区计数纯属浪费,也不把"未评分"伪装成任何分数。边界窗口
    # (起点或终点落在范围内)仍计数:范围内位置(尤其范围边界)的 3'
    # 剖面要读"以 j 结尾的 k-mer",其起点可能在范围外。范围外位置的
    # 剖面值在下方统一置 0,可视化以专属灰色标注(R13 图例)。
    # R15: positions outside the F/R design ranges (the complement of
    # left_range ∪ right_range) are not k-mer scored — counting repeats in
    # regions where no primer can land is wasted work, and "unscored" is not
    # dressed up as a score. Boundary windows (start or end inside a range)
    # are still counted so in-range positions (especially at range edges)
    # read true 3' profiles; out-of-range profile values are zeroed below
    # and the visualization marks them gray (R13 legend).
    designable = None
    if mode != "sgrna" and (params.get("left_range") or params.get("right_range")):
        designable = set()
        for rng in (params.get("left_range"), params.get("right_range")):
            if not rng:
                continue
            s, e = max(1, int(rng[0])), min(n, int(rng[1]))
            if e >= s:
                designable.update(range(s - 1, e))

    skip_kmer = bool(params.get("skip_kmer_scoring"))
    # 深度/特异性比对库(R13):spec_db 是 targetbase= 解析出的基因组,缺省
    # = db_prefix。深度剖面语义是"基因组拷贝数",必须以特异性库为基准:
    # 模板来源库若是 CDS/转录本集合,其中每条 k-mer 天然单拷贝 → 剖面
    # 全 1.0("全部单拷贝"假象),设计被当作全可用 —— R12 把 spec_db 只透传
    # 给 blast_binding_sites,深度统计仍留在 db_prefix,名称型 targetbase=
    # 流程即触发。第一步 blastn 与 k-mer 计数统一走 depth_db(整段查询还
    # 驱动 target_loci 的基因组定位,同样必须以基因组为基准)。
    # Depth/specificity db (R13): spec_db is the genome resolved from
    # targetbase=; None falls back to db_prefix. The depth profile means
    # "genome copy number" and must be counted against the specificity db:
    # a template sourced from a CDS/transcript set is naturally single-copy
    # there, yielding an all-1.0 profile (an "all single-copy" illusion) and
    # treating every region as usable. R12 threaded spec_db only into
    # blast_binding_sites, leaving the counting on db_prefix — the name-style
    # targetbase= flow hit exactly that. Step-1 blastn and k-mer counting
    # both use depth_db (the whole-query blastn also drives genomic
    # target_loci mapping, which must likewise be genome-based).
    depth_db = spec_db or db_prefix
    if depth_db != db_prefix and on_log:
        on_log(f"深度统计/特异性比对库: {depth_db} "
               f"(模板来源库: {db_prefix})")
    # k-mer 结果缓存(R17):同一项目内对同一序列重复设计时,第一步 blastn +
    # 计数(整条管线最耗时的部分,大库可达分钟级)由前端 sessionStorage /
    # 项目文件带回的缓存直接替代。键含序列/库/尺度/范围,命中即逐位一致。
    # 缓存查找前置于 windows 构建:部分命中(R28)时窗口只保留未覆盖区域。
    # K-mer result cache (R17): re-designing the same sequence within a
    # project reuses the cached step-1 blastn + counting (the most expensive
    # phase, minutes on large DBs) that the frontend carried back from
    # sessionStorage / the project file. The key covers sequence/db/scales/
    # ranges, so a hit is bit-for-bit identical. The lookup runs before the
    # window build so a partial hit (R28) can keep only the uncovered windows.
    # 跳过 k-mer 评分:不查缓存(占位剖面/全零深度无复用价值,还会把
    # 后续正常设计污染成全零),计数与剖面构建一并跳过
    # Skip k-mer scoring: no cache lookup (placeholder profiles/zero depths
    # are worthless to reuse and would poison later normal designs), and the
    # counting/profile construction are skipped below
    cache_key = cache_hit = None
    cache_uncovered = None      # None=无缓存/整体重算;空集=全命中;非空=部分命中
    cache_offset = 0            # 子串命中:当前模板在缓存模板中的偏移
    # Uncovered positions beyond the cached detected ranges (empty set = full
    # coverage; non-empty = partial overlap that must be re-counted)
    if not skip_kmer:
        cache_key = _kmer_cache_key(seq, depth_db, params)
        cache_hit, cache_uncovered, cache_offset = _find_kmer_cache(
            kmer_caches, cache_key, n, designable, seq, depth_db, mode)
        if cache_hit is None:
            cache_uncovered = None
        elif cache_offset > 0 and on_log:
            on_log(f"k-mer 缓存子串命中:当前模板为缓存模板的子串"
                   f"(偏移 {cache_offset} bp),重叠区复用缓存")
    windows: dict[str, tuple[int, int]] = {}   # qseqid → (0 基偏移, 窗长)
    if not skip_kmer:
        # 步长 = 1:逐碱基 k-mer 窗口。步长 = 窗长的平铺窗口会把位置 i 的
        # 深度取成"覆盖它的那个窗口"的命中数 —— 跨重复边界的窗口(如
        # A1 起点 100 的 96 号 12-mer 窗)因读到唯一侧翼而报"唯一",
        # max 合并后把重复区边界也放成 1.0("全绿"复发)。逐碱基窗口
        # 保证 depth[i] = 模板 i 处 k-mer 的基因组出现次数,3' 剖面也能
        # 精确读"以 j 结尾的 k-mer"。真实基因组上短尺度因 evalue 过滤
        # 自然静默,查询量(≈4n)可接受。
        # Step = 1: per-base k-mer windows. Tiling (step = window length)
        # gives position i the hit count of "the window covering i", so a
        # window straddling a repeat boundary reads the unique flank and
        # reports "unique" — max-combining then lifts repeat-edge bases
        # back to 1.0 ("all green" recurs). Per-base windows make
        # depth[i] = the genome occurrence count of the k-mer at i, and let
        # the 3' profile read the k-mer ending at j exactly. On real
        # genomes the short scales go silent via the evalue filter, so the
        # query count (~4n) stays acceptable.
        for k in WINDOW_KMERS:
            for s, e in _sliding_windows(n, k, 1):
                if designable is None or s in designable or (s + k - 1) in designable:
                    # 部分缓存命中:只保留未覆盖区域的窗口(重叠区用缓存值)
                    # Partial cache hit: keep only windows in the uncovered
                    # region (the overlap reuses the cached values)
                    if cache_hit is not None and cache_uncovered:
                        if s not in cache_uncovered and (s + k - 1) not in cache_uncovered:
                            continue
                    windows[f"__bpx{k}_{s}"] = (s, k)
    # 深度/特异性比对库(R13):spec_db 是 targetbase= 解析出的基因组,缺省
    # = db_prefix。深度剖面语义是"基因组拷贝数",必须以特异性库为基准:
    # 模板来源库若是 CDS/转录本集合,其中每条 k-mer 天然单拷贝 → 剖面
    # 全 1.0("全部单拷贝"假象),设计被当作全可用 —— R12 把 spec_db 只透传
    # 给 blast_binding_sites,深度统计仍留在 db_prefix,名称型 targetbase=
    # 流程即触发。第一步 blastn 与 k-mer 计数统一走 depth_db(整段查询还
    # 驱动 target_loci 的基因组定位,同样必须以基因组为基准)。
    # Depth/specificity db (R13): spec_db is the genome resolved from
    # targetbase=; None falls back to db_prefix. The depth profile means
    # "genome copy number" and must be counted against the specificity db:
    # a template sourced from a CDS/transcript set is naturally single-copy
    # there, yielding an all-1.0 profile (an "all single-copy" illusion) and
    # (see the depth_db/cache blocks above, which now run before the window
    # build so partial hits can filter the windows)
    if skip_kmer and on_log:
        on_log("已跳过 k-mer 评分过程(按全放行设计)")
    if cache_hit is not None:
        if cache_uncovered:
            if on_log:
                on_log(f"k-mer 缓存部分命中:检测范围重叠 "
                       f"{n - len(cache_uncovered)} bp,未覆盖 "
                       f"{len(cache_uncovered)} bp 重新计数")
        elif on_log:
            on_log("命中 k-mer 结果缓存,跳过第一步 blastn 与 k-mer 计数")
        # 新缓存携带 target_loci(模板↔库坐标映射,仅数 KB);R17 旧缓存只带
        # step1_hsps 全文(大命中可达数十 MB),此处兼容读取并在下方重新推导
        # New caches carry target_loci (the template↔db coordinate mapping,
        # a few KB); R17 caches only held the full step1_hsps (tens of MB on
        # heavy hits), read compatibly below and re-derived once
        full_hsps = cache_hit.get("step1_hsps") or []
        if on_progress:
            on_progress(_F_KMER_S + _F_KMER_SPAN)   # 跳过计数 → 直接到自重复段
    else:
        if on_log:
            if windows:
                rng = ",设计范围外不评分" if designable is not None else ""
                on_log(f"第一步:整段模板 blastn + k-mer 计数窗口查询"
                       f"({WINDOW_KMERS} bp,步长=1{rng},共 {1 + len(windows)} 个窗口)...")
            else:
                on_log("第一步:全序列 blastn 基因组比对 ...")
        if on_progress:
            on_progress(_F_STEP1_BLAST * 0.5)   # 整段 blastn 开始(子进程无内部进度)
        with _temp_query(seq, query_name) as qfa:
            hsps = blast.run_blast_tabular(
                depth_db, qfa, program="blastn",
                evalue=float(params.get("blast_evalue", 10)),
                max_targets=int(params.get("blast_max_targets", 5000)),
                # 必须禁用 DUST:特异性检测的目标正是发现重复区,
                # DUST 会把低复杂度串联重复整段屏蔽,导致重复区显示为"无匹配"
                # DUST must be disabled: specificity detection aims precisely to
                # find repeat regions, but DUST masks low-complexity tandem
                # repeats entirely, making repeat regions appear as "no match"
                # word_size 7(与 blastn-short 的 seed 约定一致)保持整段查询
                # 行为与 R11 完全一致(target_loci / step1_hsps 不变);窗口查询
                # 已由 kmer_count 纯 Python 计数替代,不再依赖 word_size。
                # word_size 7 (the blastn-short seed convention) keeps the
                # whole-query behavior identical to R11 (target_loci / step1_hsps
                # unchanged); window queries are counted natively by kmer_count
                # and no longer depend on the word size.
                extra_args=["-dust", "no", "-task", "blastn", "-word_size", "7"],
                on_log=on_log, cancel=cancel,
                timeout=int(params.get("timeout_sec", 600)),
            )
        if cancel is not None and cancel.cancelled:
            raise blast.BlastError("任务已取消")
        full_hsps = [h for h in hsps if h["qseqid"] not in windows]
        if on_progress:
            on_progress(_F_STEP1_BLAST)   # 整段 blastn 完成

        # 窗口 k-mer 出现次数(R12):纯 Python 计数程序替代窗口 blastn 查询,
        # 返回 (qids, per_scale) 紧凑结构化数组(12 B/命中,R15.1:list[dict]
        # 在 115 Mb 库实测 3.4 GB 撑爆内存),由 compute_windowed_depth_compact
        # 做与 R11 相同的位点合并。
        # Window k-mer occurrence counts (R12): the native counter replaces the
        # windowed blastn queries; it returns compact structured arrays
        # (qids, per_scale) at 12 B/hit (R15.1: the list-of-dicts form measured
        # 3.4 GB on a 115 Mb db), consumed by compute_windowed_depth_compact
        # for the same locus merging as R11.
    # 计数在 blastn 分支之外:全命中(缓存全覆盖)时无窗口要数跳过;
    # 部分命中时 windows 已按未覆盖区域过滤,只数未覆盖区(重叠区用缓存值)。
    # The counting runs outside the blastn branch: a full cache hit has no
    # windows to count; a partial hit counts only the uncovered windows (the
    # windows dict was filtered above), merging with the cached overlap.
    if skip_kmer or (cache_hit is not None and not cache_uncovered):
        qids, counted = [], {}
    else:
        qids, counted = kmer_count.count_window_occurrences(
            depth_db, seq, windows, on_log=on_log, cancel=cancel,
            on_progress=(lambda f: on_progress(_F_KMER_S + _F_KMER_SPAN * f)
                         if on_progress else None))

    target_range = None
    if target:
        t1, t2 = max(1, target[0]), min(n, target[1])
        if t2 >= t1:
            target_range = (t1, t2)

    # 目标位点(§37):主位点上覆盖目标区(±buffer)的库位点,含 qstart/qend/
    # sstrand 供模板↔库映射(预期位点构造)。缓存命中时直接取坐标映射
    # (新缓存);R17 旧缓存无 target_loci 时从 HSP 推导一次作为兼容降级
    # Target loci (§37): db loci on the primary site covering the target range
    # (±buffer), with qstart/qend/sstrand for template↔db mapping. On a cache
    # hit the mapping comes straight from the cache (new caches); R17 caches
    # without target_loci fall back to one HSP-based derivation.
    if cache_hit is not None and cache_hit.get("target_loci") is not None:
        tl = cache_hit["target_loci"]
        if cache_offset > 0:
            # 子串命中:缓存模板坐标映射到当前模板(q 坐标减 offset)
            # Substring hit: map the cached template coordinates onto the
            # current template (q coordinates minus the offset)
            target_loci = [{
                **l,
                "qstart": (int(l.get("qstart") or 0) - cache_offset)
                          if l.get("qstart") else None,
                "qend": (int(l.get("qend") or 0) - cache_offset)
                        if l.get("qend") else None,
            } for l in tl]
        else:
            target_loci = tl
    else:
        target_loci = target_loci_from_hsps(
            full_hsps, n, target_range, buffer=int(params.get("target_buffer", 50)))

    # 逐碱基深度:各 k-mer 尺度窗口位点计数(模板短于窗长时无该尺度窗口,
    # 直接用整段 HSP 的深度)。深度统计同时作旧 UI 展示兼容
    # (§43:主数据结构已是 specificity profile)。
    # Per-base depth from per-scale window locus counts; templates shorter
    # than a window size fall back to whole-query depth for that scale.
    if skip_kmer:
        # 跳过 k-mer 评分:深度/3' 深度全零占位(可视化显示为未评分,下游
        # 预筛/评分按 0 处理,四段式只有 Level 4 全放行能产出)
        # Skipped k-mer scoring: all-zero placeholder depth (the visualization
        # shows unscored; downstream prefilters/scoring read 0 and only the
        # Level-4 free pass yields candidates)
        depth_profile = SimpleNamespace(depth=[0] * n, matched=0.0,
                                        repeat_frac=0.0, histogram={})
        depth_3p = [0] * n
    elif cache_hit is not None:
        # 缓存重建:深度/剖面直接取自缓存(键含范围,设计范围外置零逻辑
        # 在下方 profile 段按同一 designable 幂等重放)。子串命中时按
        # offset 切片(当前模板 = 缓存模板的子区间)。
        # Cache rebuild: depth/profiles come straight from the cache (the key
        # covers ranges; the out-of-range zeroing replays idempotently below).
        # Substring hits slice by offset (the current template is a sub-range
        # of the cached one).
        if cache_offset > 0:
            depth_profile = SimpleNamespace(
                depth=list(cache_hit["depth"][cache_offset:cache_offset + n]),
                **cache_hit["depth_stats"])
            depth_3p = list(cache_hit["depth_3p"][cache_offset:cache_offset + n])
            # 子串切片继承了大模板的边界窗口值(大模板 designable 边界
            # 窗口落在小模板坐标上可能越出当前设计范围)——仅在"当前
            # 模板下该位置无任何保留窗口"时置 0:窗口保留条件 = 起点或
            # 终点 ∈ designable,故位置 i 有值当且仅当 ∃k: i 或 i+k-1
            # ∈ designable;边界窗口终点落在 designable 内的位置(i 本身
            # 在外)保留切片值,与全新计数一致
            # Substring slices inherit the parent template's boundary-window
            # values, which can fall outside the current design ranges. Zero
            # a position only when no window starts there under the current
            # designable (a window is kept when its start OR end falls in
            # designable, so i keeps its value iff some k has i or i+k-1 in
            # designable) — matching a fresh computation
            if designable is not None:
                for i in range(n):
                    if i not in designable and \
                            not any(i + k - 1 in designable for k in WINDOW_KMERS):
                        depth_profile.depth[i] = 0
                        depth_3p[i] = 0
        else:
            depth_profile = SimpleNamespace(depth=list(cache_hit["depth"]),
                                            **cache_hit["depth_stats"])
            depth_3p = list(cache_hit["depth_3p"])
        # 部分命中:未覆盖区域重新计数并合并进缓存数组(重叠区保持缓存值)
        # Partial hit: re-count the uncovered region and merge it into the
        # cached arrays (the overlap keeps the cached values)
        if cache_uncovered:
            off_by_k = {k: {t: off for t, (off, win) in windows.items() if win == k}
                        for k in WINDOW_KMERS}
            new_depths = {}
            for k in WINDOW_KMERS:
                if off_by_k[k]:
                    new_depths[k] = compute_windowed_depth_compact(
                        counted.get(k), qids, n, off_by_k[k])
            if new_depths:
                # 合并位置 = uncovered ∪ 窗口终点落在 uncovered 的窗口起点:
                # 全新计算的窗口保留规则是"起点或终点 ∈ designable",终点
                # 在未覆盖区的窗口(起点在重叠区)会被计数,其起点位置的
                # 深度/剖面也受新计数影响——只覆盖 uncovered 会让这些
                # 边界位置残留缓存的 0,与全新计算不一致
                # Merge positions = uncovered ∪ window starts whose window
                # ends in uncovered: the fresh counting keeps windows whose
                # end falls in the designable region, so a window starting in
                # the overlap but ending in uncovered also updates its start
                # position — covering only uncovered would leave those
                # boundary positions at the cached 0
                merge_pos = set(cache_uncovered)
                for i in cache_uncovered:
                    for k in new_depths:
                        for s in range(max(0, i - k + 1), i):
                            merge_pos.add(s)
                rep_k = max(new_depths)          # 代表性尺度(最大有窗口者)
                for i in merge_pos:
                    depth_profile.depth[i] = new_depths[rep_k].depth[i]
                if 12 in new_depths:
                    for i in merge_pos:
                        depth_3p[i] = new_depths[12].depth[i]
    else:
        off_by_k = {k: {t: off for t, (off, win) in windows.items() if win == k}
                    for k in WINDOW_KMERS}
        depths = {}
        for k in WINDOW_KMERS:
            if not off_by_k[k]:
                depths[k] = compute_depth(full_hsps, n)
                continue
            # 计数程序构造上即"全长精确匹配"(R11 的精确过滤被原生计数取代):
            # 无错配、无缺口、无 seed 级短 HSP;per_scale 已按尺度分好,无需过滤。
            # 尺度无命中 → 全零深度(与旧空列表行为一致;整段 blastn 的偶发
            # 短命中会把唯一模板抬到 2~3)。
            # The counter is exact by construction (R11's full-length exact
            # filter is subsumed): no mismatches, gaps, or seed-level short
            # HSPs exist; per_scale is already grouped by scale. A scale with no
            # hits yields an all-zero profile (matching the old empty-list
            # behavior; the whole-query blastn's chance short HSPs would lift a
            # unique template to 2~3).
            depths[k] = compute_windowed_depth_compact(
                counted.get(k), qids, n, off_by_k[k])
        depth_profile = depths[max(WINDOW_KMERS)]   # 代表性深度(旧 UI 展示/诊断)
        depth_3p = depths[12].depth                  # 3' 预筛:12-mer seed 深度
    if on_log:
        if designable is not None:
            on_log(f"设计范围外 {n - len(designable)} bp 未评分,"
                   f"深度统计仅覆盖可设计区")
        on_log(f"深度分析: 匹配率 {depth_profile.matched * 100:.1f}%, "
               f"重复区占比 {depth_profile.repeat_frac * 100:.1f}%, "
               f"平均深度 {sum(depth_profile.depth) / max(1, n):.1f}")

    # 特异性剖面:深度 → 1/(1+log2(深度))(§19 等价单调递减函数);
    # 0(全库无匹配)→ 1.0。global = 四尺度 max 合并("最强特异性证据",
    # 长 k-mer 占优,阈值可达);three_prime = 各尺度以该位置结尾的 k-mer
    # 分数取 max(与自重复分量 compute_profiles 的 3' 窗口约定一致;模板
    # 前 7 bp 无任何尺度窗口结尾 → 0,保守)。sgRNA 模式不建剖面(沿用
    # 旧通道)。
    # Specificity profiles: depth mapped through 1/(1+log2(depth)) (§19).
    # global = max-combined across the four scales (the longest k-mer
    # dominates, thresholds stay reachable); three_prime = max over the
    # scales of the k-mer ending at the position (same 3'-window convention
    # as the self-repeat component's compute_profiles; the first 7 bp have
    # no k-mer ending there → 0, conservative).
    # 特异性剖面:所有模式(sgRNA 同样参与 k-mer 评分与四段式,用户需求)
    # 构建;skip_kmer 时全零占位剖面(四段式自然只有 Level 4 全放行能产出,
    # 特异性评分按 0 处理)。缓存命中直接取缓存剖面。
    # Specificity profiles: built for every mode (sgRNA also takes part in
    # k-mer scoring and the staged design, per user request); skip_kmer yields
    # an all-zero placeholder (only the Level-4 free pass produces candidates,
    # specificity scores read 0). Cache hits reuse the cached profiles.
    if skip_kmer:
        # skip_kmer:全零占位剖面——四段式各级阈值永不可达(仅 Level 4
        # 全放行有产出),特异性评分读 0;深度占位全零供可视化显示"未评分"
        # skip_kmer: all-zero placeholder profiles — no level threshold is
        # reachable (only the Level-4 free pass yields candidates) and
        # specificity scores read 0; the all-zero depth placeholder drives
        # the "unscored" visualization
        profiles = {"global": [0.0] * n, "three_prime": [0.0] * n}
    elif cache_hit is not None:
        cp = cache_hit["profiles"]
        if cache_offset > 0:
            profiles = {"global": list(cp["global"][cache_offset:cache_offset + n]),
                        "three_prime": list(cp["three_prime"][cache_offset:cache_offset + n])}
        else:
            profiles = dict(cp)
        # 部分命中:未覆盖区域的剖面按新计数深度重算并合并
        # Partial hit: recompute the uncovered region's profiles from the
        # freshly counted depths and merge
        if cache_uncovered and new_depths:
            scores = {k: [round(kmer_occurrence_score(v), 4)
                          for v in new_depths[k].depth] for k in new_depths}
            # merge_pos 与深度合并一致:uncovered ∪ 边界窗口起点
            merge_pos = set(cache_uncovered)
            for i in cache_uncovered:
                for k in new_depths:
                    for s in range(max(0, i - k + 1), i):
                        merge_pos.add(s)
            for i in merge_pos:
                g = 0.0
                t = 0.0
                for k, sc in scores.items():
                    if sc[i] > g:
                        g = sc[i]
                    j = i - k + 1
                    if j >= 0 and sc[j] > t:
                        t = sc[j]
                profiles["global"][i] = g
                profiles["three_prime"][i] = t
    else:
        scores = {k: [round(kmer_occurrence_score(v), 4)
                      for v in depths[k].depth] for k in WINDOW_KMERS}
        profiles = {
            "global": [max(scores[k][i] for k in WINDOW_KMERS)
                       for i in range(n)],
            "three_prime": [
                max((scores[k][i - k + 1] for k in WINDOW_KMERS
                     if i >= k - 1), default=0.0)
                for i in range(n)],
        }
    if designable is not None:
            # 范围外位置无 k-mer 分数 → 置 0(未评分):色带灰标(R13)、
            # 曲线归零、概览统计只算范围内。范围内位置不受影响(边界
            # 窗口按上述保留规则仍计数)。
            # Out-of-range positions have no k-mer score → 0 ("unscored"):
            # the band shows gray (R13), the curve drops to zero, and the
            # overview stats cover only in-range positions.
            for i in range(n):
                if i not in designable:
                    profiles["global"][i] = 0.0
                    profiles["three_prime"][i] = 0.0

    base = {
        "query": query_name, "template_len": n, "success": False,
        "depth": depth_profile.depth, "depth_stats": {
            "matched": depth_profile.matched, "repeat_frac": depth_profile.repeat_frac,
            "histogram": depth_profile.histogram,
        },
        # 注意:profile_stats 不在此处设置——标准模式需在自重复分量合并
        # 之后记录(见 _design_with_profile),sgRNA 在 _design_sgrna 设置
        # (profile_stats is set after the self-repeat merge in
        # _design_with_profile; _design_sgrna sets it for sgRNA)
        "target": {"start": target_range[0] if target_range else 1,
                   "end": target_range[1] if target_range else n},
        # step1_hsps 只保留坐标等轻量字段:qseq/sseq 全长序列在大命中时可把
        # 结果撑到数百 MB(实测 29.6 万 HSP → 156 MB),前端渲染与深度统计
        # 均不需要它们;缓存命中的新缓存不含 HSP 全文,此处即为空列表
        # step1_hsps keeps only lightweight fields: full-length qseq/sseq
        # bloated results to hundreds of MB on heavy hits (measured 296k HSPs
        # -> 156 MB) and neither rendering nor depth statistics need them;
        # caches without the full HSP list hit this as an empty list
        "step1_hsps": [{k: h.get(k) for k in
                        ("qstart", "qend", "sseqid", "sstart", "send", "pident", "evalue", "bitscore", "sstrand", "length")}
                       for h in full_hsps],
        "stages": [],
        "levels": [],
        # 完整四段式阈值配置随结果下发(含未实际跑到的级别——设计可能
        # Level 1 即成功,levels 只含已跑级别,可视化色带/参考线/图例需要
        # 完整阈值才能按参数绘制;与子函数 levels 解析同源)。
        # The full four-level threshold config rides with the result (levels
        # only lists stages actually reached — a Level-1 success would
        # otherwise leave L2/L3 thresholds out, breaking the parameterised
        # band/reference lines).
        "level_thresholds": (params.get("specificity_levels")
                              or _levels_from_params(params)),
    }
    # k-mer 结果缓存随结果回传:前端存入 sessionStorage 并在项目保存时写入
    # 文件(同一序列再次设计/换机器加载项目 → 免重算)。命中缓存时原样回传,
    # 未命中时打包本次计算结果(幂等,前端重复存储无害)。
    # The k-mer cache travels with the result so the frontend can store it in
    # sessionStorage and into project files (re-designing the same sequence,
    # or loading the project elsewhere, skips the recomputation). On a cache
    # hit the entry is echoed back; otherwise the computed analysis is packed
    # (idempotent — re-storing is harmless).
    if cache_hit is not None:
        if cache_uncovered:
            # 部分命中:合并检测范围与数据写回新缓存——重叠区保留原值,
            # 未覆盖区为新计数,检测范围取并集(下次更大范围可继续复用);
            # 子串命中时缓存范围先按 offset 映射到当前模板坐标
            # Partial hit: merge the detected ranges and data into a new cache
            # — the overlap keeps its values, the uncovered region is freshly
            # counted, and the ranges are unioned (larger future ranges can
            # keep reusing); substring hits map the cached ranges onto the
            # current template first
            c_ranges = cache_hit.get("detected_ranges") or []
            if cache_offset > 0:
                # 子串命中:切片已把当前设计范围外的旧值置 0(边界窗口例外
                # 但保守视为未检测),检测范围 = 当前设计范围本身——并集会
                # 引入未检测区间(负坐标/范围外),下次复用会拿到错误的全 0
                # Substring hits zero the out-of-range old values in the slice
                # (boundary windows are conservatively treated as undetected),
                # so the detected ranges equal the current design ranges —
                # unioning would pull in undetected intervals (negative
                # coordinates / out-of-range) that reuse would read as bogus 0
                merged_detected = _detected_ranges(designable, n)
            else:
                merged_detected = _merge_ranges(
                    c_ranges, _detected_ranges(designable, n))
            merged = dict(cache_hit)
            merged["key"] = cache_key
            merged["template_len"] = n
            merged["seq"] = seq
            merged["mode"] = "sgrna" if mode == "sgrna" else "pair"
            merged["depth"] = depth_profile.depth
            merged["depth_3p"] = depth_3p
            merged["profiles"] = profiles
            merged["detected_ranges"] = merged_detected
            base["kmer_cache"] = merged
        else:
            base["kmer_cache"] = cache_hit
    elif not skip_kmer:
        base["kmer_cache"] = {
            "key": cache_key,
            "template_len": n,
            # 缓存自描述元数据:所属特异性库、基因(模板名)、检测范围
            # (排除 F/R 均不设计而未检测的区域)——复用前据此判断覆盖;
            # seq/mode 供子串命中(flank 变化导致模板嵌套)映射坐标
            # Self-describing metadata: specificity db, gene (query name) and
            # detected ranges (excluding regions where neither primer can be
            # designed) — reuse checks coverage against these; seq/mode feed
            # substring hits (nested templates from flank changes)
            "db": depth_db,
            "gene": query_name,
            "seq": seq,
            "mode": "sgrna" if mode == "sgrna" else "pair",
            "detected_ranges": _detected_ranges(designable, n),
            "depth": depth_profile.depth,
            "depth_3p": depth_3p,
            "profiles": profiles,
            "depth_stats": {
                "matched": depth_profile.matched,
                "repeat_frac": depth_profile.repeat_frac,
                "histogram": depth_profile.histogram,
            },
            # target_loci 替代 step1_hsps 全文:复用只需坐标推导预期位点,
            # HSP 全文使缓存膨胀到数十 MB、进不了 sessionStorage,复用机制
            # 形同虚设。结果里的 step1_hsps 已瘦身为坐标,缓存不再重复携带
            # target_loci replaces the full step1_hsps: reuse only needs the
            # coordinates to derive expected loci; the full HSP list bloated
            # the cache to tens of MB, never fitting sessionStorage, and made
            # reuse pointless. The slimmed step1_hsps already rides in the
            # result, so the cache does not duplicate it.
            "target_loci": target_loci,
        }

    if mode == "sgrna":
        base = _design_sgrna(seq, base, db_prefix, params, on_log, cancel,
                             spec_db, profiles)
    else:
        base = _design_with_profile(seq, base, target_range, target_loci,
                                    index, profiles, depth_3p,
                                    params, db_prefix, on_log,
                                    cancel, spec_db, on_progress,
                                    designable)

    if locate_ctx:
        _annotate_locate(base, locate_ctx)
    return base


# ---------------------------------------------------------------- 新管线 (New pipeline)

def _design_with_profile(
    seq: str, base: dict, target_range: tuple[int, int] | None,
    target_loci: list[dict], index: KmerIndex | None,
    profiles: dict | None, depth_3p: list[int], params: dict,
    db_prefix: str, on_log: Callable[[str], None] | None,
    cancel: CancelFlag | None, spec_db: str | None,
    on_progress: Callable[[float], None] | None,
    designable: set[int] | None = None,
) -> dict:
    n = base["template_len"]
    mode = params.get("mode", "standard")
    single_mode = mode == "single"

    if profiles is None:
        # skip_kmer:无剖面 → 单级 Level 4 全放行设计(不再视为失败;
        # 特异性评分路径各自处理无剖面场景)
        # skip_kmer: no profiles -> a single Level-4 free pass (no longer a
        # failure; each scoring path handles the profile-less case)
        levels = [None]
        if on_log:
            on_log("已跳过 k-mer 评分过程:按 Level 4 全放行设计")
    else:
        levels = params.get("specificity_levels") or _levels_from_params(params)

    if profiles is not None and index is not None:
        # 补充模板内自重复 k-mer 分(§16-§19,仅模板索引):窗口深度把同一条目
        # 内相邻拷贝(subject 区间重叠)合并为同一位点,串联重复/微卫星因此
        # 显示为"单拷贝";自重复分量捕获这些模板内部拷贝,与深度分逐碱基
        # 取 min —— 只加强掩蔽、从不放宽,深度剖面仍是主信号。
        # Self-repeat k-mer component (template-only index): window depth
        # merges adjacent in-entry copies (overlapping subject ranges) into
        # one locus, so tandem repeats would read as "single-copy"; this
        # component recovers in-template copies and is min'ed elementwise
        # into the depth profile (it can only mask more, never relax).
        if on_log:
            on_log("补充模板内自重复 k-mer 剖面(8/10/12/15-mer)...")
        self_profiles = compute_profiles(
            index, seq,
            kmer_set=tuple(int(k) for k in params.get("specificity_kmers", [8, 10, 12, 15])),
            three_prime_windows=tuple(int(w) for w in params.get("three_prime_windows", [8, 10, 12, 15])),
            on_progress=(lambda f: on_progress(_F_PROFILE * (_F_SELF_S + 0.08 * f))
                         if on_progress else None),
            cancel=cancel,
        )
        if cancel is not None and cancel.cancelled:
            raise blast.BlastError("任务已取消")
        profiles = {
            "global": [min(a, b) for a, b in
                       zip(profiles["global"], self_profiles["global"])],
            "three_prime": [min(a, b) for a, b in
                            zip(profiles["three_prime"], self_profiles["three_prime"])],
        }
    # 最终剖面(含自重复合并)随结果下发,前端可视化/测试读取
    # The final profile (after the self-repeat merge) rides with the result
    base["profile_stats"] = {
        "global": profiles["global"] if profiles is not None else [],
        "three_prime": profiles["three_prime"] if profiles is not None else [],
    }
    # F/R 设计范围外区域(R13):正、反向引物 3' 端都无法放置的区间补集
    # (left_range ∪ right_range 之外)。级内循环把它并入 primer3 excluded
    # (primer3 从范围内出候选),但深度剖面仍照常计数——这些位置不是重复
    # (可能显示"单拷贝绿"),只是设计上不可用;随结果下发 range_excluded,
    # 供可视化以专属颜色标注并加图例。
    # R13: the F/R design-range complement (no primer of either strand can
    # land there) is folded into primer3's excluded regions inside the level
    # loop below, but the depth profile still counts normally — these
    # positions are not repeats (they may read "single-copy green"), they
    # are simply not designable. Expose them with the result so the
    # visualization can mark them in their own color with a legend entry.
    range_ex = _range_excluded(n, params.get("left_range"),
                               params.get("right_range"))
    if range_ex:
        base["range_excluded"] = range_ex

    buffer_len = int(params.get("buffer_len", 8))
    max_hits = int(params.get("prefilter_max_hits", 200))
    candidate_count = int(params.get("candidate_count", 50))

    accepted: list[dict] = []
    level_used = None
    counts = {"primer3": 0, "postfilter": 0, "prefilter": 0, "local": 0, "pair": 0}
    for li, th in enumerate(levels):
        level_no = li + 1
        ls = _F_LEVEL_START + _F_LEVEL_SPAN * li / len(levels)
        le = _F_LEVEL_START + _F_LEVEL_SPAN * (li + 1) / len(levels)
        if on_log:
            on_log(f"{LEVEL_NAMES.get(level_no, 'Level')} 开始")
        if th is None:
            g_th = t_th = 0.0
            allowed = None          # 全放行(§23 Level 4)
            allowed_runs = None
            excluded = None
        else:
            g_th, t_th = float(th[0]), float(th[1])
            allowed = [profiles["global"][j] >= g_th
                       and profiles["three_prime"][j] >= t_th
                       for j in range(n)]
            if not any(allowed):
                if on_log:
                    on_log(f"{LEVEL_NAMES.get(level_no)}: 无满足阈值的允许区,下一级")
                base["levels"].append(_level_info(level_no, th, 0, 0, 0, 0, 0))
                continue
            # 可用区两侧缓冲:允许区 = 每个连续可用区上下游各外扩 buffer_len;
            # 原始可用区区间另行交给级内后置过滤,保证引物 3' 端落在特异区
            # Both-sides buffer: the designable region is each available run
            # extended by buffer_len on both sides; the original runs go to the
            # level postfilter so a primer's 3' end always lands in specific
            # territory while its 5' end may sit in the buffer band
            allowed_runs = _allowed_runs(allowed)
            excluded = _excluded_with_both_sides_buffer(allowed, buffer_len,
                                                        designable)
            if on_log:
                on_log(f"{LEVEL_NAMES.get(level_no)}: 允许区(含两侧缓冲) "
                       f"{n - sum(l for _, l in excluded)} bp / 屏蔽 "
                       f"{sum(l for _, l in excluded)} bp")
        # F/R 设计范围外区域并入 excluded:primer3 不知道范围,只做后置
        # 过滤会把范围外候选全灭 → 让 primer3 直接从范围内出候选(§6)
        # Fold the F/R design-range complement into excluded so primer3 picks
        # inside the ranges instead of only the post-filter rejecting them (§6)
        if range_ex:
            excluded = _merge_excluded((excluded or []) + range_ex)

        # --- primer3 物理设计(§24-§25,§45)---
        # 不把 target_range 传给 primer3:SEQUENCE_TARGET 强制产物覆盖整个
        # 目标区,目标触及模板起点或超过产物上限时 primer3 必然 0 候选
        # (定位模式曾因此全灭,浏览器任务 candidate_count=0);目标区已由
        # target_loci 负责命中判定与评分,primer3 的位置约束只有 excluded。
        # target_range is deliberately not passed to primer3: SEQUENCE_TARGET
        # forces the product to span the whole target region, which yields
        # zero candidates whenever the target touches the template start or
        # exceeds the product range (this broke locate mode entirely); the
        # target region is handled by target_loci for hit/scoring, and
        # primer3's placement constraints are excluded_regions alone.
        if on_progress:
            on_progress(ls + (le - ls) * 0.2)
        pool = design_primers(seq, excluded_regions=excluded,
                              params=params, num_return=candidate_count,
                              single=single_mode)
        n_p3 = len(pool)
        if on_log:
            on_log(f"primer3 候选 {n_p3} 对")
        counts["primer3"] += n_p3
        if not pool:
            base["levels"].append(_level_info(level_no, th, n_p3, 0, 0, 0, 0))
            continue
        pool = _diversity_pool(pool, max_pairs=candidate_count)

        # --- 级内后置过滤:3' 端必须落在原始可用区(方向性缓冲)+ 设计范围 + 二聚体 + GC 夹子 ---
        pool = _level_postfilter(pool, profiles, g_th, t_th, params,
                                 single=single_mode, allowed_runs=allowed_runs)
        n_post = len(pool)
        if on_log:
            on_log(f"级内过滤后 {n_post} 对")
        counts["postfilter"] += n_post
        if not pool:
            base["levels"].append(_level_info(level_no, th, n_p3, n_post, 0, 0, 0))
            continue

        # --- 3' seed 深度预筛(§26:快速淘汰明显不特异候选)---
        if on_progress:
            on_progress(ls + (le - ls) * 0.5)
        kept = []
        for p in pool:
            f3p = p["left"]["start"] + p["left"]["len"] - 1
            if _kmer_prefilter(depth_3p, f3p, max_hits):
                continue
            if not single_mode:
                r3p = p["right"]["start"] - p["right"]["len"] + 1
                if _kmer_prefilter(depth_3p, r3p, max_hits):
                    continue
            kept.append(p)
        pool = kept
        n_pre = len(pool)
        if on_log:
            on_log(f"k-mer 预筛后 {n_pre} 对")
        counts["prefilter"] += n_pre
        if not pool:
            base["levels"].append(_level_info(level_no, th, n_p3, n_post, n_pre, 0, 0))
            continue

        # --- binding-site 验证(§27-§32)+ pair-level PCR(§33-§40)---
        # 批处理:本级全部候选引物合并一次 blastn-short 全库(库只加载一次)
        # skip_spec_eval 开启时跳过逆向验证:特异性分直接取 3' 端 k-mer
        # 深度剖面(基因组拷贝数代理),速度大幅提升、但无错配豁免/脱靶
        # 产物判定 —— 3' 端深度已由 k-mer 预筛(≤ max_hits)把关。
        # With skip_spec_eval the reverse check is skipped entirely: the
        # specificity score comes straight from the 3'-end k-mer depth
        # profile (a genome copy-number proxy) — much faster, but no
        # mismatch exemption / off-target-product judgement; the 3'-end
        # depth is still bounded by the k-mer prefilter (≤ max_hits).
        skip_spec = bool(params.get("skip_spec_eval"))
        if on_progress:
            on_progress(ls + (le - ls) * 0.7)
        if skip_spec:
            if on_log:
                on_log(f"{LEVEL_NAMES.get(level_no)}: 跳过 blastn-short 逆向验证"
                       f"(特异性分按 k-mer 深度剖面)")
            for p in pool:
                _score_pair_kmer_only(p, profiles, params, single_mode)
        else:
            bf_primers = []
            for pi, p in enumerate(pool):
                bf_primers.append((p["left"]["seq"], f"f{pi}"))
                if p.get("right", {}).get("seq"):
                    bf_primers.append((p["right"]["seq"], f"r{pi}"))
            hits_map = blast_binding_sites(bf_primers, spec_db or db_prefix, params,
                                           on_log, cancel)
            for pi, p in enumerate(pool):
                f_res = hits_to_binding_sites(p["left"]["seq"],
                                              hits_map.get(f"f{pi}", []), params)
                f_sites = f_res["sites"]
                p["seed_hits"] = {"forward": f_res["internal_hits"], "reverse": 0}
                r_res = None
                if p.get("right", {}).get("seq"):
                    r_res = hits_to_binding_sites(p["right"]["seq"],
                                                  hits_map.get(f"r{pi}", []), params)
                    r_sites = r_res["sites"]
                    p["seed_hits"]["reverse"] = r_res["internal_hits"]
                else:
                    r_sites = []
                if f_res["truncated"] or (r_res is not None and r_res["truncated"]):
                    # 命中数 > prefilter_max_hits(§26 快速淘汰):引物在库中
                    # 大量出现(重复序列),sites 为空只是截断假象——绝不能
                    # 当"无非目标位点"判 UNIQUE(100)。独立 level 淘汰,
                    # 与"可扩增脱靶产物"区分开(UI 误读为后者会误导用户)。
                    # Hits beyond prefilter_max_hits (§26 fast elimination):
                    # the primer matches the db massively (repetitive); the
                    # empty sites are a truncation artefact and must never be
                    # read as "no off-target site" (UNIQUE 100). Rejected
                    # under its own level, distinct from off-target
                    # amplification (users misread the shared label).
                    n_trunc = max(f_res["internal_hits"],
                                  (r_res or {}).get("internal_hits", 0))
                    ev = {
                        "level": PAIR_TRUNCATED,
                        "label": PAIR_LABELS[PAIR_TRUNCATED],
                        "spec_score": PAIR_SCORES[PAIR_TRUNCATED],
                        "off_target_sites": n_trunc,
                        "note": f"命中数超过上限({n_trunc}),无法完整评估脱靶位点,已淘汰",
                        "note_key": "spec_truncated",
                        "note_params": {"n": n_trunc},
                        "amplifiable_pairs": [],
                    }
                    _score_pair(p, ev, [], [], profiles, params, single_mode)
                    continue
                if single_mode:
                    ev = classify_single_primer(
                        f_sites, target_loci, params)
                else:
                    ev = evaluate_primer_pair(p["left"]["seq"], p["right"]["seq"],
                                              f_sites, r_sites, target_loci,
                                              params, cancel)
                _score_pair(p, ev, f_sites, r_sites, profiles, params, single_mode)
        if on_progress:
            on_progress(ls + (le - ls) * 0.9)

        survivors = [p for p in pool if not p["specificity"]["reject"]]
        eliminated = [p for p in pool if p["specificity"]["reject"]]
        counts["local"] += n_pre
        counts["pair"] += len(survivors)
        base["levels"].append(_level_info(level_no, th, n_p3, n_post, n_pre,
                                          len(pool), len(survivors)))
        for p in survivors:
            p["level"] = level_no
        if survivors:
            accepted = survivors
            level_used = level_no
            base["stages"].append({
                "stage": level_no, "name": LEVEL_NAMES.get(level_no, f"Level {level_no}"),
                "pairs": survivors, "eliminated": eliminated, "success": True,
                "thresholds": th,
            })
            break
        # 本级别全部被 pair 判定淘汰 → 下一级(更宽松)
        if on_log:
            on_log(f"{LEVEL_NAMES.get(level_no)}: 全部候选被特异性判定淘汰,下一级")
        base["stages"].append({
            "stage": level_no, "name": LEVEL_NAMES.get(level_no, f"Level {level_no}"),
            "pairs": [], "eliminated": eliminated, "success": False,
            "thresholds": th,
        })

    if accepted:
        accepted.sort(key=_rank_key)
        base["success"] = True
        base["level_used"] = level_used
        base["pairs"] = accepted
        base["stage_reached"] = level_used
        if on_log:
            on_log(f"设计成功(Level {level_used}): {len(accepted)} 对引物")
    else:
        base["failure"] = _failure_diagnosis(base, params, counts,
                                             level_used=level_used,
                                             profile=profiles["global"],
                                             db=db_prefix or "")
    if on_progress:
        on_progress(1.0)
    return base


def _score_pair_kmer_only(p: dict, profiles: dict | None, params: dict,
                          single_mode: bool) -> None:
    """skip_spec_eval 模式下的评分:不做 blastn-short 逆向验证,特异性分
    取两条引物 3' 端 k-mer 深度剖面(three_prime 剖面 = 以该位置结尾的
    k-mer 基因组出现次数的分数)的较小值 —— 深度 1 → 100,深度越高越接近
    0;不淘汰任何候选(无脱靶产物判定信息)。profiles=None(skip_kmer)
    时无剖面可评分,spec_score 置 0 并附说明。

    Scoring under skip_spec_eval: no reverse blastn-short check; the
    specificity score is the smaller of the two primers' 3'-end k-mer depth
    profile values (the three_prime profile = the occurrence score of the
    k-mer ending at that position) — depth 1 → 100, decreasing with depth.
    No candidate is rejected (no off-target-product information exists).
    profiles=None (skip_kmer) leaves no profile to score: spec_score is 0
    with an explanatory note.
    """
    if profiles is not None:
        n = len(profiles["global"])
        l3p = int(p["left"]["start"]) + int(p["left"]["len"]) - 1
        r3p = (int(p.get("right", {}).get("start", l3p))
               - int(p.get("right", {}).get("len", 1)) + 1)
        t = (min(profiles["three_prime"][l3p - 1],
                 profiles["three_prime"][r3p - 1])
             if 1 <= l3p <= n and 1 <= r3p <= n else 0.0)
        note = "未进行 blastn-short 逆向验证,特异性分仅来自 k-mer 深度剖面"
        note_key = "spec_skip_note"
    else:
        t = 0.0
        note = "已跳过 k-mer 评分与 blastn-short 逆向验证,特异性分未评估"
        note_key = "spec_skipped_all"
    ev = {
        "level": PAIR_KMER_ONLY,
        "spec_score": round(t * 100, 1),
        "off_target_sites": 0,
        "note": note,
        "note_key": note_key,
        "note_params": None,
        "amplifiable_pairs": [],
    }
    _score_pair(p, ev, [], [], profiles, params, single_mode)


def _score_pair(p: dict, ev: dict, f_sites: list, r_sites: list,
                profiles: dict | None, params: dict, single_mode: bool) -> None:
    """候选对评分(§39):物理分 + 特异性拆分 + 综合分;reject 标记(§40)。
    profiles=None(skip_kmer)时 global/3' 分记 0。"""
    n = len(profiles["global"]) if profiles is not None else 0
    l3p = int(p["left"]["start"]) + int(p["left"]["len"]) - 1
    r3p = (int(p.get("right", {}).get("start", l3p))
           - int(p.get("right", {}).get("len", 1)) + 1)
    if profiles is not None and 1 <= l3p <= n and 1 <= r3p <= n:
        g = min(profiles["global"][l3p - 1], profiles["global"][r3p - 1])
        t = min(profiles["three_prime"][l3p - 1], profiles["three_prime"][r3p - 1])
    else:
        g = t = 0.0
    p["physical_score"] = round(physical_score(p, params), 1)
    p["global_score"] = round(g * 100, 1)
    p["three_prime_score"] = round(t * 100, 1)
    p["pair_spec_score"] = float(ev["spec_score"])
    level = ev["level"]
    reject = (level == PAIR_OFFTARGET or level == PAIR_TRUNCATED
              or (single_mode and level == PAIR_NO_PRODUCT))
    p["specificity"] = {
        "level": level,
        "label": PAIR_LABELS.get(level, SPEC_LABELS.get(level, level)),
        "spec_score": ev["spec_score"],
        "off_target_sites": ev.get("off_target_sites", 0),
        "note": ev.get("note", ""),
        # 结构化本地化键:spec_eval 各分支携带 note_key/note_params,
        # 前端语言切换时按当前语言渲染;旧数据缺失则降级用中文 note
        # Structured i18n keys: spec_eval branches carry note_key/note_params,
        # rendered per language by the frontend; legacy data falls back to `note`
        "note_key": ev.get("note_key"),
        "note_params": ev.get("note_params"),
        "amplifiable_pairs": ev.get("amplifiable_pairs", []),
        "reject": reject,
        "seed_hits": p.get("seed_hits", {"forward": 0, "reverse": 0}),
    }
    p["binding_sites"] = {
        "forward": [s.to_dict() for s in f_sites],
        "reverse": [s.to_dict() for s in r_sites],
    }
    p["composite_score"] = composite_score(
        p["physical_score"], p["pair_spec_score"], ev.get("off_target_sites", 0), params)


def _rank_key(p: dict) -> tuple:
    """排序优先级(§40):reject 最高优先,然后 pair spec / 3' / global / 物理 / 位置。"""
    return (1 if p["specificity"]["reject"] else 0,
            -p["pair_spec_score"],
            -p["three_prime_score"],
            -p["global_score"],
            p["penalty"],
            p["left"]["start"])


def _kmer_prefilter(depth_3p: list[int], pos: int, max_hits: int) -> bool:
    """3' seed 深度预筛(§26 的窗口深度实现):引物 3' 端所在位置的
    12-mer 窗口位点数(全基因组拷贝数)> max_hits → 明显不可能特异,
    快速淘汰。O(1)。

    窗口深度 ≥ 其内任意 3' seed 的命中数(blastn 以 seed 扩展并报告
    HSP,窗口命中即包含 seed 命中),因此深度超阈值 ⇒ seed 必超阈值,
    预筛只承担粗剪枝;精确判定交给 blastn-short 与 3' 端延伸。

    Depth-based 3' seed prefilter (§26): if the 60 bp window at the
    primer's 3' end matches more than max_hits genome loci, the primer is
    obviously non-specific and is dropped early.
    """
    if not (1 <= pos <= len(depth_3p)):
        return False
    return depth_3p[pos - 1] > max_hits


def _in_allowed_runs(pos: int, runs: list[tuple[int, int]] | None) -> bool:
    """1-based 坐标是否落在任一原始可用区(0-based 闭区间)内;runs=None 全放行。"""
    if runs is None:
        return True
    for s, e in runs:
        if s <= pos - 1 <= e:
            return True
    return False


def _level_postfilter(pairs: list[dict], profiles: dict, g_th: float, t_th: float,
                      params: dict, single: bool = False,
                      allowed_runs: list[tuple[int, int]] | None = None) -> list[dict]:
    """级内后置过滤:

    - 两条引物 3' 端必须落在原始可用区内(可用区两侧缓冲的方向性约束:
      primer3 的排除区允许整条引物落在两侧外扩的缓冲带,此处把 3' 端
      拉回特异区——正向引物 3' 端在右、反向引物 3' 端在左,因此缓冲带
      只可能被对应引物的 5' 端占用,与"5' 端可落入非特异区、3' 端在
      特异区内"的方向性语义一致);
    - F/R 设计范围(定位模式 §6);二聚体 ≤ max_dimer;3' GC 夹子限制。
    """
    max_dimer = int(params.get("max_dimer", 5))
    max_clamp = int(params.get("max_gc_clamp_3p", 3))
    n = len(profiles["global"])
    out = []
    for p in pairs:
        left = p["left"]
        l3p = left["start"] + left["len"] - 1
        if not (1 <= l3p <= n):
            continue
        if not _in_range(l3p, params.get("left_range")):
            continue
        if not _in_allowed_runs(l3p, allowed_runs):
            continue
        if not single:
            right = p["right"]
            r3p = right["start"] - right["len"] + 1  # 右引物 3' 端 = 覆盖区左端 (right primer 3' end = coverage left edge)
            if not (1 <= r3p <= n):
                continue
            if not _in_range(right["start"], params.get("right_range")):
                continue
            if not _in_allowed_runs(r3p, allowed_runs):
                continue
            if gc_clamp_3p(left["seq"]) > max_clamp or gc_clamp_3p(right["seq"]) > max_clamp:
                continue
            d1 = dimer_stats(left["seq"], right["seq"])
            d2 = dimer_stats(left["seq"], blast.revcomp(right["seq"]))
            if max(d1["max_consec"], d2["max_consec"]) >= max_dimer:
                continue
            p["dimer"] = {"max_consec": max(d1["max_consec"], d2["max_consec"]),
                          "max_total": max(d1["max_total"], d2["max_total"])}
        else:
            if gc_clamp_3p(left["seq"]) > max_clamp:
                continue
            p["dimer"] = {"max_consec": 0, "max_total": 0}
        out.append(p)
    return out


def _allowed_to_excluded(allowed: list[bool]) -> list[tuple[int, int]]:
    """允许位点布尔数组 → primer3 excluded 区列表 (start, length) 1-based。"""
    n = len(allowed)
    ex: list[tuple[int, int]] = []
    i = 0
    while i < n:
        if allowed[i]:
            i += 1
            continue
        j = i
        while j < n and not allowed[j]:
            j += 1
        ex.append((i + 1, j - i))
        i = j
    return ex


def _allowed_runs(allowed: list[bool]) -> list[tuple[int, int]]:
    """允许位点布尔数组 → 连续可用区区间列表 (0-based 闭区间)。"""
    runs: list[tuple[int, int]] = []
    n = len(allowed)
    i = 0
    while i < n:
        if not allowed[i]:
            i += 1
            continue
        s = i
        while i < n and allowed[i]:
            i += 1
        runs.append((s, i - 1))
    return runs


def _excluded_with_both_sides_buffer(allowed: list[bool],
                                     buffer_len: int,
                                     designable: set[int] | None = None
                                     ) -> list[tuple[int, int]]:
    """允许位点布尔数组 → primer3 excluded 区列表:"可用区两侧缓冲"——
    每个连续可用区上下游各外扩 buffer_len 作为设计允许区,之外全部排除。
    designable 非空时缓冲外扩**裁剪到设计范围内**:范围边界处的可用区
    不得向范围外扩展(否则 primer3 会在范围外的缓冲带出候选,引物 5' 端
    超出用户规定的 F/R 区域)。

    方向性(引物 5' 端可落入非特异缓冲带,3' 端必须在特异区内):
    - 正向引物 5' 端在左 → 左扩缓冲带供其 5' 端落入;
    - 反向引物 5' 端在右 → 右扩缓冲带供其 5' 端落入;
    - 两条引物的 3' 端位置由级内后置过滤约束在**原始可用区**内
      (`_level_postfilter` 用 `_allowed_runs` 检查),缓冲带自身
      profile < 阈值,primer3 排除区又保证整条引物不越过扩展边界。
    buffer_len=0 时退化为 `_allowed_to_excluded`。

    Both-sides available-region buffer: every contiguous available run is
    extended by buffer_len on each side to form the designable region;
    everything outside is excluded from primer3. When designable is given
    the extension is clipped to the design ranges — an available run at a
    range boundary must not extend outward, or primer3 would place
    candidates in the out-of-range buffer band (5' end beyond the user's
    F/R region). The 3' end of either primer is constrained to the original
    available runs by the level postfilter (`_allowed_runs`), so the 5' end
    may fall into the non-specific buffer band while the 3' end always stays
    in specific territory. buffer_len=0 reproduces `_allowed_to_excluded`.
    """
    n = len(allowed)
    runs = _allowed_runs(allowed)
    # designable(0-based 位置集合)→ 区间列表,供裁剪
    d_runs: list[tuple[int, int]] = []
    if designable:
        dp = sorted(designable)
        ds = de = dp[0]
        for p in dp[1:]:
            if p == de + 1:
                de = p
            else:
                d_runs.append((ds, de))
                ds = de = p
        d_runs.append((ds, de))
    ext: list[tuple[int, int]] = []
    for s, e in runs:
        cs, ce = max(0, s - buffer_len), min(n - 1, e + buffer_len)
        # 缓冲外扩裁剪到设计范围内:范围边界处不向范围外扩展
        # Clip the extension to the design ranges: no outward expansion at
        # range boundaries
        if d_runs:
            for ds, de in d_runs:
                a, b = max(cs, ds), min(ce, de)
                if a > b:
                    continue
                if ext and a <= ext[-1][1] + 1:
                    ext[-1] = (ext[-1][0], max(ext[-1][1], b))
                else:
                    ext.append((a, b))
        else:
            if ext and cs <= ext[-1][1] + 1:    # 扩展后相邻/重叠 → 合并
                ext[-1] = (ext[-1][0], max(ext[-1][1], ce))
            else:
                ext.append((cs, ce))
    ex: list[tuple[int, int]] = []
    pos = 0
    for s, e in ext:                            # 扩展区补集 → excluded
        if s > pos:
            ex.append((pos + 1, s - pos))
        pos = e + 1
    if pos < n:
        ex.append((pos + 1, n - pos))
    return ex


def _merge_excluded(regions: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """相邻/重叠 excluded 区间合并 (start, length) 1-based,primer3 安全。"""
    if not regions:
        return []
    seg = sorted((s, s + l - 1) for s, l in regions)  # (start, end) 1-based 闭区间
    merged: list[tuple[int, int]] = []
    for s, e in seg:
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return [(s, e - s + 1) for s, e in merged]


def _range_product_conflict(params: dict, n: int) -> str | None:
    """F/R 设计范围与产物长度范围的几何矛盾诊断。

    双侧范围存在时:F 3' 端限于 left_range、R 3' 端限于 right_range,产物
    长度 = R3' − F3' 被范围钳制在 [right0 − left1, right1 − left0] 内;
    该区间与产物 [min,max] 无交 → primer3 物理上无解(默认滑块范围与
    默认产物 150-300 在长目标区间时即触发此情形)。返回中文建议文案。

    Geometric contradiction between F/R design ranges and the product-length
    range. With both ranges set, the product length is constrained to
    [right0 - left1, right1 - left0]; if that interval does not intersect the
    product [min, max], primer3 has no physical solution.
    """
    info = _range_product_conflict_info(params, n)
    return info["reason"] if info else None


def _range_product_conflict_info(params: dict, n: int) -> dict | None:
    """几何矛盾诊断的结构化版本:文案 + 结构化键(供前端 i18n)。

    Structured variant of the geometric-conflict diagnostic: message plus
    i18n keys/params for the frontend (the Chinese `reason` stays for
    backward compatibility and tests).
    """
    left = params.get("left_range")
    right = params.get("right_range")
    if not left or not right:
        return None
    l0, l1 = max(1, int(left[0])), min(n, int(left[1]))
    r0, r1 = max(1, int(right[0])), min(n, int(right[1]))
    if l1 >= r0:  # 范围重叠 → 产物长度无下限约束,不构成矛盾
        return None
    pmin = float(params.get("product_len_min") or 0)
    pmax = float(params.get("product_len_max") or 0)
    min_possible, max_possible = r0 - l1, r1 - l0
    if pmax and min_possible > pmax:
        return {
            "reason": (f"正向/反向设计范围要求产物至少 {min_possible} bp"
                       f"(F 范围终点 {l1}、R 范围起点 {r0}),"
                       f"超过产物最大长度 {pmax:.0f} bp:"
                       f"primer3 无解。请增大产物长度范围,或将滑块向目标区间内侧收拢"),
            "reason_key": "range_product_conflict_min",
            "reason_params": {"min_possible": min_possible, "l1": l1,
                              "r0": r0, "pmax": pmax},
        }
    if pmin and max_possible < pmin:
        return {
            "reason": (f"正向/反向设计范围允许产物至多 {max_possible} bp"
                       f"(F 范围起点 {l0}、R 范围终点 {r1}),"
                       f"低于产物最小长度 {pmin:.0f} bp:"
                       f"primer3 无解。请缩小产物长度范围,或将滑块向外侧展开"),
            "reason_key": "range_product_conflict_max",
            "reason_params": {"max_possible": max_possible, "l0": l0,
                              "r1": r1, "pmin": pmin},
        }
    return None


def _range_excluded(n: int, left_range=None, right_range=None) -> list[tuple[int, int]]:
    """设计范围外区域 → primer3 excluded 列表 (start, length) 1-based。

    primer3 本身不知道 F/R 设计范围(级内过滤只是后置淘汰,范围外候选
    直接全灭);把范围外的区域并入 excluded,primer3 从范围内出候选。
    范围端点钳制到 [1, n];无范围 → 空列表。

    Design-range complement as primer3 excluded regions (start, length),
    1-based. primer3 has no knowledge of the F/R design ranges (the in-level
    filter only rejects candidates afterwards, killing every out-of-range
    pair); excluding the complement lets primer3 pick inside the ranges.
    """
    if not left_range and not right_range:
        return []
    ranges = [r for r in (left_range, right_range) if r]
    covered = set()
    for r in ranges:
        s, e = max(1, int(r[0])), min(n, int(r[1]))
        if e >= s:
            covered.update(range(s, e + 1))
    ex: list[tuple[int, int]] = []
    i = 1
    while i <= n:
        if i in covered:
            i += 1
            continue
        j = i
        while j <= n and j not in covered:
            j += 1
        ex.append((i, j - i))
        i = j
    return ex


def _level_info(level_no: int, th, cand: int, after_post: int, after_pre: int,
                after_local: int, after_pair: int) -> dict:
    return {
        "level": level_no,
        "name": LEVEL_NAMES.get(level_no, f"Level {level_no}"),
        "thresholds": th,
        "candidates": cand, "after_postfilter": after_post,
        "after_prefilter": after_pre, "after_local": after_local,
        "after_pair_check": after_pair,
    }


# ---------------------------------------------------------------- 失败诊断 (Failure diagnosis)

def _failure_diagnosis(base: dict, params: dict, counts: dict,
                       level_used: int | None,
                       profile: list[float], db: str = "") -> dict:
    """失败诊断(§42):failure_stage + 各级候选统计 + 低特异区列表 + 建议。"""
    n = base["template_len"]
    # 设计范围外区域(R15)未评分(剖面为 0),不参与低特异区判定 ——
    # 否则范围外区域会被当作"最低特异区"报告,误导用户
    # Out-of-range positions (R15) are unscored (profile 0); they must not
    # be reported as low-specificity regions
    excl = base.get("range_excluded") or []
    xmask = None
    if excl:
        xmask = [False] * n
        for st, ln in excl:
            for j in range(st - 1, min(n, st - 1 + ln)):
                xmask[j] = True
    # 低特异性区:global < 0.5 的连续段(按最低分排序)
    lows = []
    start = None
    min_s = 1.0
    for i in range(n):
        if xmask and xmask[i]:
            if start is not None:
                lows.append((start + 1, i, round(min_s, 3)))
                start = None
            continue
        s = profile[i]
        if s < 0.5:
            if start is None:
                start = i
                min_s = s
            min_s = min(min_s, s)
        elif start is not None:
            lows.append((start + 1, i, round(min_s, 3)))
            start = None
    if start is not None:
        lows.append((start + 1, n, round(min_s, 3)))
    lows.sort(key=lambda r: r[2])
    if not lows:
        # 全模板高特异仍失败:取最低分单点(如 primer3 无候选);
        # 跳过未评分的范围外位置
        js = [j for j in range(n) if not (xmask and xmask[j])] or [0]
        i = min(js, key=lambda j: profile[j])
        lows = [(i + 1, i + 1, round(profile[i], 3))]

    single_mode = params.get("mode") == "single"
    # 按最具体情况选择阶段(§42):primer3 无产出 → 物理失败优先;
    # 全部候选被 k-mer 预筛淘汰 → 无特异 3' 位置;其余 → 特异性判定失败。
    # Pick the most specific stage (§42): no primer3 output → physical
    # failure first; all candidates killed by the k-mer prefilter → no
    # specific 3' positions; otherwise → specificity failure.
    if counts["primer3"] == 0 or counts["postfilter"] == 0:
        stage = FAIL_PRIMER3
        reason = "primer3 未能在允许区设计出合格引物(可放宽 Tm/GC/产物长度范围)"
        reason_key, reason_params = "physical_no_candidate", {}
        conflict = _range_product_conflict_info(params, n)
        if conflict:
            reason = conflict["reason"]
            reason_key = conflict["reason_key"]
            reason_params = conflict["reason_params"]
    elif counts["prefilter"] == 0:
        stage = FAIL_NO_HIGH_SPEC
        reason = "所有候选引物 3' 端 seed 均过于重复,模板内无特异位置可用"
        reason_key, reason_params = "no_specific_3p", {}
    elif counts["pair"] == 0:
        stage = FAIL_SINGLE_SPEC if single_mode else FAIL_PAIR_OFFTARGET
        if single_mode:
            reason = "全部候选引物存在可延伸的脱靶结合位点"
            reason_key = "offtarget_single"
        else:
            reason = "全部候选引物对存在可扩增脱靶产物"
            reason_key = "offtarget_pair"
        reason_params = {}
    else:
        stage = FAIL_NO_PAIR
        reason = "存在合格候选但无引物对通过全部特异性判定"
        reason_key, reason_params = "no_pair", {}

    return {
        "reason": reason,
        # 结构化 i18n 键:前端据此本地化(旧 `reason` 中文文案保留作降级)
        # Structured i18n key: frontend localizes from this (Chinese `reason` kept as fallback)
        "reason_key": reason_key,
        "reason_params": reason_params,
        "stage": level_used or 0,
        "failure_stage": stage,
        "template_length": n,
        "database": db,
        "candidate_count": counts["primer3"],
        "candidate_after_prefilter": counts["prefilter"],
        "candidate_after_local_check": counts["local"],
        "candidate_after_pair_check": counts["pair"],
        # 旧 UI 兼容键(§42 新键见上)
        "repeat_stats": {
            "matched": base["depth_stats"]["matched"],
            "repeat_frac": base["depth_stats"]["repeat_frac"],
            "top_repeat_regions": [
                {"start": r[0], "end": r[1], "max_depth": 2 if r[2] < 0.5 else 1}
                for r in lows[:10]],
        },
        "top_low_specificity_regions": [
            {"start": r[0], "end": r[1], "min_score": r[2]} for r in lows[:10]],
        "suggestion": ("建议放宽: 特异性等级阈值 / 候选数 / Tm 与 GC 范围 / 产物长度范围;"
                       "或换用更长、重复度更低的模板区域"),
        "suggestion_key": "relax_params",
    }


def _failure(base: dict, params: dict, stage: str, reason: str,
             n: int, db: str = "", reason_key: str = "other",
             reason_params: dict | None = None) -> dict:
    return {
        "reason": reason, "stage": 0, "failure_stage": stage,
        "reason_key": reason_key,
        "reason_params": reason_params or {},
        "template_length": n, "database": db,
        "candidate_count": 0, "candidate_after_prefilter": 0,
        "candidate_after_local_check": 0, "candidate_after_pair_check": 0,
        "repeat_stats": {"matched": 0.0, "repeat_frac": 0.0,
                         "top_repeat_regions": []},
        "top_low_specificity_regions": [],
        "suggestion": "检查数据库索引与序列文件",
        "suggestion_key": "check_db_index",
    }


# ---------------------------------------------------------------- 附加模式 (Additional modes)

def _sgrna_3p_index(g: dict) -> int:
    """sgRNA 引导序列 3' 端在模板上的 0-based 位置:正链 = start+len-2;
    负链引导序列为模板上游 PAM 的 revcomp,3' 端对应模板 start-1。

    0-based position of a guide's 3' end on the template: plus strand =
    start+len-2; minus-strand guides are the revcomp of the upstream PAM, so
    the 3' end sits at template start-1.
    """
    if g["strand"] == "+":
        return g["left"]["start"] + g["left"]["len"] - 2
    return g["left"]["start"] - 1


def _sgrna_score(g: dict, base: dict, target_loci: list[dict] | None,
                 params: dict, profiles: dict | None,
                 skip_spec: bool) -> None:
    """单个 sgRNA 候选的特异性/物理/综合评分。skip_spec → k-mer 剖面;
    否则 blastn-short 评估(evaluate_single_primer)。skip_kmer 后无剖面
    可评分时 spec_score 置 0 并附 note(诚实:未做任何特异性评估)。

    Score one guide: skip_spec -> k-mer profile; otherwise blastn-short
    (evaluate_single_primer). With skip_kmer there is no profile to score
    against, so spec_score drops to 0 with a note (honest: no specificity
    assessment was made).
    """
    g["right"] = {}
    if skip_spec:
        # 3' 剖面评分(与标准模式 _score_pair_kmer_only 同语义):guide 3' 端
        # 位置的 three_prime 剖面分数 ×100;skip_kmer 占位剖面全 0 → 0 分
        # Score via the 3' profile (same semantics as standard
        # _score_pair_kmer_only): the three_prime profile at the guide's 3'
        # end ×100; the skip_kmer placeholder profile reads 0
        dep_idx = _sgrna_3p_index(g)
        prof = profiles["three_prime"][dep_idx] if (
            profiles is not None and 0 <= dep_idx < len(profiles["three_prime"])) else 0.0
        spec_score = round(prof * 100, 1)
        spec = {
            "level": PAIR_KMER_ONLY,
            "label": PAIR_LABELS.get(PAIR_KMER_ONLY, PAIR_KMER_ONLY),
            "spec_score": spec_score,
            "off_target_sites": 0,
            "reject": False,
        }
        if profiles is None:
            spec.update({
                "note": "已跳过 k-mer 评分与 blastn-short 逆向验证,特异性分未评估",
                "note_key": "spec_skipped_all",
                "note_params": None,
            })
        g["specificity"] = spec
    else:
        ev = evaluate_single_primer(g["left"]["hits"], target_loci, params,
                                    g["left"].get("primer_len"))
        level = ev["level"]
        g["specificity"] = {
            "level": level,
            "label": SPEC_LABELS.get(level, "有可延伸脱靶命中"),
            "spec_score": SPEC_SCORES.get(level, 30),
            "off_target_sites": ev.get("off_target_sites", 0),
            "reject": False,
        }
    gcval = g["left"]["gc"]
    physical = max(0.0, 100.0 - abs(gcval - 50) * 2)
    g["physical_score"] = round(physical, 1)
    g["composite_score"] = composite_score(
        physical, g["specificity"]["spec_score"],
        g["specificity"].get("off_target_sites", 0), params)


def _design_sgrna(seq, base, db_prefix, params, on_log, cancel,
                  spec_db: str | None = None,
                  profiles: dict | None = None) -> dict:
    """sgRNA:20 bp 引导序列 + PAM(NGG),两条链扫描,GC 40%~60%。

    四段式(k-mer 剖面分级,与标准引物对同 levels 阈值):guide 的 3' 端
    位置必须落在该级允许区(global/three_prime ≥ 阈值)内,逐级放宽,
    有产出即成功;Level 4 全放行。skip_kmer(占位全零剖面)自然只有
    Level 4 能产出。skip_spec → 按 k-mer 剖面评分,否则 blastn-short。

    sgRNA: 20 bp guide + PAM (NGG), both strands scanned, GC 40%~60%.
    Staged by the k-mer profiles (same level thresholds as standard pairs):
    a guide's 3' end must sit in the level's allowed region, thresholds relax
    level by level, the first level with output wins; Level 4 passes all.
    skip_kmer (all-zero placeholder profiles) yields output only at Level 4.
    skip_spec -> k-mer profile scoring, otherwise blastn-short.
    """
    if on_log:
        on_log("sgRNA 模式:双链扫描 20 bp + PAM(NGG) ...")
    # sgRNA 的剖面随结果下发(前端可视化/测试读取;sgRNA 无自重复合并)
    # The sgRNA profile rides with the result (visualization/tests; sgRNA
    # does not run the self-repeat merge)
    if profiles is not None:
        base["profile_stats"] = {
            "global": profiles["global"],
            "three_prime": profiles["three_prime"],
        }
    sgrna_len = int(params.get("sgrna_len", 20))
    pam = params.get("sgrna_pam", "NGG").upper()
    if len(pam) != 3 or pam[1:3] != "GG":
        raise ValueError("PAM 需为 NGG 形式")

    guides: list[dict] = []
    n = len(seq)
    rc = blast.revcomp(seq)
    for strand, s, side in (("+", seq, "sense"), ("-", rc, "antisense")):
        for i in range(n - sgrna_len + 1):
            g = s[i:i + sgrna_len]
            if any(c not in "ACGT" for c in g):
                continue
            if not (40 <= gc(g) <= 60):
                continue
            if strand == "+":
                pam_seq = seq[i + sgrna_len:i + sgrna_len + 3]
                ok = len(pam_seq) == 3 and pam_seq[1:3] == "GG"
                start = i + 1
                guide_seq = g
            else:
                # 负链:引导序列绑定模板 [i, i+len),PAM 为模板上游 CCN(展示为 NGG 形式)
                # Minus strand: the guide binds the template at [i, i+len);
                # PAM is the upstream CCN on the template (shown in NGG form)
                pam_seq = seq[n - i - sgrna_len - 3: n - i - sgrna_len]
                ok = len(pam_seq) == 3 and pam_seq[0:2] == "CC"
                start = n - i - sgrna_len + 1
                guide_seq = g
                pam_seq = blast.revcomp(pam_seq)   # CCN → NGG
            if not ok:
                continue
            guides.append({
                "left": {
                    "seq": guide_seq, "start": start, "len": sgrna_len,
                    "tm": round(tm(guide_seq)[0], 2), "gc": gc(guide_seq),
                    "penalty": 0.0, "self_any": 0.0, "self_end": 0.0, "hairpin": 0.0,
                },
                "strand": strand, "pam": pam_seq, "side": side,
                "penalty": 0.0,
            })

    # 定位模式(或指定目标区时):"引导序列必须位于目标区段内"选项(默认开启,§7)
    # Locate mode (or when a target region is specified): the "guide must lie
    # within the target region" option (enabled by default, §7)
    if params.get("sgrna_target_only", True):
        t1, t2 = base["target"]["start"], base["target"]["end"]
        guides = [g for g in guides
                  if t1 <= g["left"]["start"] <= t2 - g["left"]["len"] + 1]
        if not guides:
            base["stage_reached"] = 4
            base["failure"] = _failure(
                base, params, FAIL_NO_HIGH_SPEC,
                "目标区段内未找到 GC 40%~60% 的 NGG PAM 位点"
                "(「引导序列必须位于目标区段内」已开启)", n,
                reason_key="no_ngg_pam_in_target")
            return base
    skip_spec = bool(params.get("skip_spec_eval"))
    if not skip_spec:
        if on_log:
            on_log(f"扫描到 {len(guides)} 个 sgRNA 候选,"
                   f"执行 blastn-short 特异性评估 ...")
        for g in guides:
            g["right"] = {}
        run_specificity_blast(guides, spec_db or db_prefix, params,
                              on_log=on_log, cancel=cancel)
    target_loci = None
    if not skip_spec:
        target_loci = target_loci_from_hsps(
            base["step1_hsps"], base["template_len"],
            (base["target"]["start"], base["target"]["end"]),
            buffer=int(params.get("target_buffer", 50)))

    # 四段式:逐级放宽 k-mer 剖面阈值,guide 3' 端位置必须落在该级允许区;
    # 有产出即成功(sgRNA 候选已过 GC/PAM 过滤,分级只约束位置)。
    # Staged design: thresholds relax level by level; a guide's 3' end must
    # sit in the level's allowed region. The first level with output wins
    # (candidates already passed GC/PAM filters; staging constrains position).
    levels = params.get("specificity_levels") or _levels_from_params(params)
    accepted: list[dict] = []
    level_used = None
    for li, th in enumerate(levels):
        level_no = li + 1
        if on_log:
            on_log(f"{LEVEL_NAMES.get(level_no, 'Level')} 开始 (sgRNA)")
        if th is None or profiles is None:
            pool = guides
        else:
            g_th, t_th = float(th[0]), float(th[1])
            pool = []
            for g in guides:
                dep_idx = _sgrna_3p_index(g)
                if 0 <= dep_idx < len(profiles["global"]) and \
                        profiles["global"][dep_idx] >= g_th and \
                        profiles["three_prime"][dep_idx] >= t_th:
                    pool.append(g)
        if not pool:
            if on_log:
                on_log(f"{LEVEL_NAMES.get(level_no)}: 无满足阈值的 sgRNA 候选,下一级")
            continue
        if on_log:
            on_log(f"{LEVEL_NAMES.get(level_no)}: {len(pool)} 个 sgRNA 候选"
                   + ("(跳过 blastn-short,按 k-mer 剖面评分)" if skip_spec else ""))
        for g in pool:
            _sgrna_score(g, base, target_loci, params, profiles, skip_spec)
        accepted = pool
        level_used = level_no
        break
    if not accepted:
        base["stage_reached"] = 4
        base["failure"] = _failure(
            base, params, FAIL_NO_HIGH_SPEC,
            "四段式各级均无满足 k-mer 剖面阈值的 sgRNA 候选"
            "(GC 40%~60% 且 3' 端位于允许区)", n,
            reason_key="no_sgrna_candidates")
        return base
    accepted.sort(key=lambda g: (-g["composite_score"], g["left"]["start"]))
    base["stage_reached"] = 4
    base["success"] = True
    base["pairs"] = accepted
    base["stages"].append({"stage": level_used or 4,
                           "name": f"sgRNA 扫描 + 特异性评估(Level {level_used or 4})",
                           "pairs": accepted, "success": True})
    return base


# ---------------------------------------------------------------- 辅助 (Helpers)

def _temp_query(seq: str, name: str):
    """整段模板查询 FASTA 临时文件(第一步 blastn 只查整段;窗口 k-mer
    由 kmer_count 纯 Python 计数,不再进 blastn)。

    Temp FASTA for the step-1 whole-template blastn query; window k-mers are
    counted natively (kmer_count) and never go through blastn.
    """
    import tempfile
    import contextlib

    @contextlib.contextmanager
    def _cm():
        lines = [f">{name}\n" + _wrap_seq(seq)]
        with tempfile.NamedTemporaryFile(
                "w", suffix=".fa", prefix="blastprime_q_", delete=False) as f:
            f.write("\n".join(lines) + "\n")
            path = f.name
        try:
            yield path
        finally:
            try:
                import os
                os.unlink(path)
            except OSError:
                pass
    return _cm()


def _sliding_windows(n: int, win: int, step: int) -> list[tuple[int, int]]:
    """[0, n) 上的滑动窗口 [(s, s+win)];末尾补一窗保证覆盖到 n。
    模板短于窗长 → 无窗口(直接用整段 HSP 的深度)。

    Sliding windows [(s, s+win)] over [0, n); a final window is appended so
    the template tail is covered. Templates shorter than the window size
    yield no windows (whole-query depth is used instead).
    """
    if n <= win:
        return []
    starts = list(range(0, n - win + 1, step))
    if starts[-1] + win < n:
        starts.append(n - win)
    return [(s, s + win) for s in starts]


def _wrap_seq(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def _in_range(pos: int | None, rng) -> bool:
    """模板 1-based 坐标是否落在设计范围 [start, end] 内(§6,前端换算后传入)。

    Whether a 1-based template coordinate falls within the design range
    [start, end] (§6; coordinates are converted by the frontend before passing).
    """
    if pos is None or not rng:
        return True
    return int(rng[0]) <= int(pos) <= int(rng[1])


def _diversity_pool(pairs: list[dict], max_pairs: int = 50) -> list[dict]:
    """按 penalty 排序后,以坐标多样性均衡选取(同一 50 bp 窗口至多保留 1 对)。

    After sorting by penalty, pick evenly by coordinate diversity
    (at most 1 pair kept per 50 bp window).
    """
    if not pairs:
        return []
    pairs = sorted(pairs, key=lambda p: (p["penalty"], p["left"]["start"]))
    picked = []
    seen_windows: set[int] = set()
    for p in pairs:
        if len(picked) >= max_pairs:
            break
        w = p["left"]["start"] // 50
        if w in seen_windows:
            continue
        seen_windows.add(w)
        picked.append(p)
    if len(picked) < max_pairs:
        for p in pairs:
            if len(picked) >= max_pairs:
                break
            if p not in picked:
                picked.append(p)
    return picked


# ---------------------------------------------------------------- 定位模式标注 (Locate-mode annotation)

def _annotate_locate(base: dict, ctx: dict) -> None:
    """定位模式标注(guide_sup1.md §5/§7/§10):
    引物坐标换算为基因组绝对坐标,逐对标记"产物是否覆盖目标区"。

    换算规则固定为:正链 abs = extract_start + t;负链 abs = extract_end − t
    (t 为模板内 0 基偏移,§5.3)。产物区间 = 正向引物 3' 端 → 反向引物 3' 端
    (与 _postfilter 的 3' 端锚定约定一致:right.start 即反向引物 3' 端)。

    Locate-mode annotation (guide_sup1.md §5/§7/§10):
    converts primer coordinates to absolute genomic coordinates and marks,
    per pair, whether the product covers the target region.

    Fixed conversion rules: plus strand abs = extract_start + t; minus strand
    abs = extract_end − t (t is the 0-based offset within the template, §5.3).
    Product interval = forward primer 3' end → reverse primer 3' end
    (consistent with the 3'-end anchoring convention: right.start is the
    reverse primer 3' end).
    """
    minus = ctx.get("strand") in ("minus", "-")
    ext_s = int(ctx.get("extract_start", 0))
    ext_e = int(ctx.get("extract_end", 0))

    def abs_of(p1: int) -> int:
        return ext_e - p1 + 1 if minus else ext_s + p1 - 1

    tg = base.get("target") or {}
    t1 = int(tg.get("start", 1))
    t2 = int(tg.get("end", base["template_len"]))
    for p in base.get("pairs", []):
        l = p.get("left") or {}
        if l.get("start"):
            p["left"]["abs_start"] = abs_of(int(l["start"]))
            p["left"]["abs_end"] = abs_of(int(l["start"]) + int(l["len"]) - 1)
        r = p.get("right") or {}
        if r.get("seq"):
            p["right"]["abs_start"] = abs_of(int(r["start"]))
            p["right"]["abs_end"] = abs_of(int(r["start"]) + int(r["len"]) - 1)
            f3p = int(l["start"]) + int(l["len"]) - 1
            r3p = int(r["start"]) - int(r["len"]) + 1  # 右引物 3' 端 = 覆盖区左端
            prod_s, prod_e = min(f3p, r3p), max(f3p, r3p)
            p["covers_target"] = prod_s <= t1 and prod_e >= t2
            if not p["covers_target"]:
                p["coverage_note"] = (
                    f"产物 {prod_s}-{prod_e} 未完整覆盖目标区 {t1}-{t2}"
                    "(受设计范围或重复区屏蔽影响)")
                # 结构化 i18n 键(前端本地化;旧字段保留作降级)
                p["coverage_note_key"] = "coverage_note"
                p["coverage_note_params"] = {"prod_s": prod_s, "prod_e": prod_e,
                                             "t1": t1, "t2": t2}
        else:
            p["covers_target"] = None  # 单引物/sgRNA 不受产物覆盖约束(§7) (Single primer/sgRNA are not subject to product-coverage constraints (§7))

    # 覆盖失败的候选降级但保留展示:排在不覆盖候选之后(§7)
    # Candidates failing coverage are demoted but kept for display:
    # sorted after non-covering candidates (§7)
    if base.get("pairs"):
        base["pairs"].sort(key=lambda p: (
            1 if p.get("covers_target") is False else 0,
            -(p.get("composite_score") or 0),
            (p.get("left") or {}).get("start", 0)))

    # 目标 HSP 的基因组坐标(与 extract 端点的换算,用于摘要与导出)
    # Genomic coordinates of the target HSP (converted from the extract
    # endpoints, used in the summary and export)
    if minus:
        g1, g2 = ext_e - t2 + 1, ext_e - t1 + 1
    else:
        g1, g2 = ext_s + t1 - 1, ext_s + t2 - 1
    base["locate"] = {
        "entry": ctx.get("entry", ""),
        # name= 显示名(前端摘要头优先展示;缺省回退 entry)
        # The name= display name (the frontend summary shows it first;
        # falls back to entry when absent)
        "display_name": ctx.get("display_name") or ctx.get("entry", ""),
        "strand": "minus" if minus else "plus",
        "extract_start": ext_s, "extract_end": ext_e,
        "target_start": t1, "target_end": t2,
        "target_genomic": [g1, g2],
    }
