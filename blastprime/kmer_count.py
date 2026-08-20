"""纯 Python + numpy 精确 k-mer 出现次数计数(R12,替代引物设计第一步的
窗口 blastn 查询)。

窗口 blastn 查询的本质是精确 k-mer 出现次数统计:blastn 的 seed 延伸/
对齐产物被 R11 的全长精确过滤(length ≥ k 且 mismatch/gapopen = 0)全部
丢弃——本模块直接实现该统计:每尺度一趟**numpy 滚动哈希**(4 进制精确
编码,模板 k-mer ↔ 整数一一对应,无碰撞)扫描库序列(正链 + 反向互补链),
输出与管线兼容的合成 HSP({qseqid, sseqid, sstart, send},1-based subject
坐标,正反链统一到正链坐标空间),由 masking.compute_windowed_depth 完成
位点合并(重叠/单碱基间隙 → 1 位点)。

为什么不用正则/逐 k-mer 搜索(实测 26 Mb 库,400 bp 模板):
  * sre 交替 `a|b|c|...` 在每位置逐分支尝试,26 Mb × ~400 分支 ≈ 20 s/尺度,
    且 finditer 是左起非重叠语义——起点落在上一命中区间内的重叠命中被吞
    (实测 k=8 丢失 ~27%),blastn 的每个窗口查询是独立的,必须全部报告;
  * 逐 k-mer str.find 每 k-mer 一趟全库扫描,≈ 11 s/尺度;
  * numpy 滚动哈希一趟 O(库长),四尺度 × 双链对 26 Mb 库 < 1 s。
  numpy 为 Biopython 既有依赖,不新增第三方包。

与 R11 窗口 blastn + 精确过滤语义逐位一致:
  * 出现次数 = 全长精确匹配(构造上无错配、无缺口、无 seed 级短 HSP);
  * 双链:blastn 按 sstrand 报告正反链命中,本模块扫两链、坐标归一,
    同一 subject 坐标上的正反链命中可正常合并;
  * 模板含 N/IUPAC 的 k-mer 无法精确匹配 → 无命中(深度 0),与 blastn 一致
    (含非 ACGT 的窗口哈希为垃圾值,命中后字符串校验剔除,绝不误报);
  * 模板内部重复的 k-mer 字符串 → 每个对应窗口各得一份命中(与 blastn
    每个窗口查询独立报告一致);
  * 重叠出现(如 poly-A 内逐位命中)**不合并**——由 compute_windowed_depth
    合并,与 R11 位点语义一致。

Pure-Python+numpy k-mer occurrence counting (R12): replaces the windowed
blastn queries of design step 1. Per-scale rolling-hash (exact base-4
encoding, collision-free) scans of the database (forward +
reverse-complement strands, chunked to bound memory) emit synthetic HSPs
consumed by masking.compute_windowed_depth, preserving the R11 locus
semantics exactly. numpy comes via Biopython — no new dependency.
"""

from __future__ import annotations

import io
from pathlib import Path
from typing import Callable

import numpy as np
from Bio import SeqIO

from . import blast
from .blast import BlastError, CancelFlag

_COMPLEMENT = str.maketrans("ACGTUN", "TGCAAN")
_ACGT = frozenset("ACGT")

# 碱基 → 0-3;非 ACGT(IUPAC/N)→ -1:含它的窗口哈希变垃圾值,命中后
# 字符串校验剔除(与 blastn 无法精确匹配含 N 查询一致)。int8 使每块
# arr 数组只有 8 MB(int64 则为 64 MB)。
# Base → 0-3; non-ACGT (IUPAC/N) → -1: windows containing one get a
# garbage hash, filtered out by the post-hit string check. int8 keeps the
# per-chunk arr array at 8 MB (64 MB as int64).
_TABLE = np.full(256, -1, dtype=np.int8)
for _i, _c in enumerate("ACGT"):
    _TABLE[ord(_c)] = _i

# 每块 ≤ 8 Mb:arr 8 MB(int8) + h 64 MB(int64) + 命中掩码,~90 MB 峰值
# (染色体级条目安全);块间重叠 k-1 碱基保证窗口完整、起点不重不漏。
# Chunk ≤ 8 Mb: arr 8 MB (int8) + h 64 MB (int64) + hit mask, ~90 MB
# peak (chromosome-scale entries stay safe); chunks overlap by k-1 bases
# for complete windows.
_CHUNK = 8 * 1024 * 1024


def _revcomp(seq: str) -> str:
    """反向互补(输入须为大写)。

    Reverse complement (input must be uppercase).
    """
    return seq.translate(_COMPLEMENT)[::-1]


def _kmer_hash(kmer: str) -> int:
    h = 0
    for ch in kmer:
        h = h * 4 + int(_TABLE[ord(ch)])
    return h


def _window_hashes(arr: np.ndarray, k: int, pow4: np.ndarray) -> np.ndarray:
    """每位置窗口的 4 进制精确编码(0..4^k-1,一一对应,无碰撞)。

    向量化滚动递推 h[i+1] = 4·h[i] + δ[i] 有顺序依赖(右侧 h[:-1] 含尚未
    计算出的 h[1],直接写会读到 np.empty 垃圾值),不能整段向量化;改卷积
    h[i] = Σ_j arr[i+j]·4^(k-1-j) = np.convolve(arr, pow4[::-1], 'valid')
    (C 实现,O(L·k),k ≤ 15 时 ~0.1 s/8 Mb)。

    Exact base-4 encodings of every window (0..4^k-1, bijective — no
    collisions). The sequential rolling recurrence cannot be vectorized
    directly (its RHS reads h[1] before it is written); a convolution with
    the reversed weights is C-implemented and O(L·k) (k <= 15 → ~0.1 s per
    8 Mb chunk).
    """
    if len(arr) < k:
        return np.empty(0, dtype=np.int64)
    return np.convolve(arr, pow4[::-1], mode="valid")


def _membership_lookup(k: int, vals: np.ndarray) -> tuple[np.ndarray, np.ndarray | None]:
    """构建成员判定查找表(每尺度一次,全部块复用)。

    小尺度(4^k ≤ 2^24,k ≤ 12):直接布尔表,`tbl[h]` O(L) 判定,最大
    16.7 MB(k=12)。大尺度(k=15,哈希范围 2^30):两级——高 16 位布尔表
    (64 KB)预筛 + 命中候选的小数组精确校验。

    替代 `np.isin(h, vals)`:哈希范围 ≤ 2^24 时 numpy 内部还能走 table
    快速路径,但 k=15 的 2^30 范围会落入 unique/argsort/searchsorted 的
    sort 路径——实测每 8 Mb 块 1.26 s(小尺度 0.08~0.10 s 的 14 倍),且
    每次块扫描产生 ~400 MB 排序临时数组;真实基因组库(115 Mb)上该尺度
    累计数分钟,内存压力触发换页后膨胀到数十分钟,表现为"卡死在尺度
    15"(R15.1 修复)。

    Returns (tbl, None) for k ≤ 12 — membership is `tbl[h]`; or
    (tbl_hi, sorted_vals) for k ≥ 13 — prefilter `tbl_hi[h >> 14]` then
    verify the exact hash at candidates via sorted_vals.
    """
    if k <= 12:                        # 直接布尔表:4^12 = 16.7 MB
        tbl = np.zeros(1 << (2 * k), dtype=bool)
        tbl[vals] = True
        return tbl, None
    hi = vals >> 14                    # [0, 2^16) 高 16 位预筛表
    tbl_hi = np.zeros(1 << 16, dtype=bool)
    tbl_hi[hi] = True
    return tbl_hi, np.sort(vals)


def _scan_strand(text: str, L: int, k: int,
                 pow4: np.ndarray, hash2kmer: dict[int, str],
                 kmer2qids: dict[str, list[str]], lookup: tuple,
                 out: list[tuple[int, int, int]], qid_index: dict[str, int],
                 forward: bool, cancel: CancelFlag | None) -> None:
    """扫描一条链(正链或反链)的所有块,追加 (qidx, sstart, send) 到 out。

    反链命中(rc 位置 p,0 基)→ 正链 1-based [L-(p+k)+1, L-p];正链 →
    [p+1, p+k]。

    **命中预合并为"位点"**:同一条目同一链上、同一 k-mer 的命中若起点差
    ≤ k(区间在正链坐标中重叠或相邻),下游 compute_windowed_depth 必然把
    它们并成 1 个位点 —— 本层提前合并成一条区间,把同聚物/微卫星/重复
    区产生的逐位置命中压成 1 条 HSP,避免 11M 级 dict 列表撑爆内存
    (115 Mb 库实测 3.4 GB)。合并是幂等的:下游 _merge_loci 对已合并的
    区间再合并得到相同位点集合(最终位点只取决于区间集合)。

    Hits are pre-merged into loci: hits of the same k-mer on the same entry
    and strand whose start positions differ by ≤ k (their ranges overlap or
    touch in forward coordinates) would necessarily collapse into one locus
    downstream in compute_windowed_depth — merging here turns homopolymer /
    microsatellite / repeat runs into a single HSP instead of one dict per
    position (a 115 Mb db measured 3.4 GB of dicts). The merge is idempotent:
    _merge_loci re-merging these intervals yields the same loci set.
    """
    n = len(text)
    start = 0
    runs: dict[str, list[int]] = {}    # kmer → [run 首命中, 尾命中] (0 基)

    def flush(kmer: str, first: int, last: int) -> None:
        if forward:
            ss, se = first + 1, last + k
        else:
            ss, se = L - (last + k) + 1, L - first
        for qid in kmer2qids[kmer]:
            out.append((qid_index[qid], ss, se))

    while start < n - k + 1:
        if cancel is not None and cancel.cancelled:
            raise blast.BlastError("任务已取消")
        end = min(start + _CHUNK, n - k + 1)
        chunk = text[start:end + k - 1]  # 末尾多取 k-1 保证窗口完整
        arr = _TABLE[np.frombuffer(chunk.encode("ascii"), dtype=np.uint8)]
        h = _window_hashes(arr, k, pow4)   # len = end - start
        tbl, svals = lookup
        # 含非 ACGT 的窗口哈希为负数,不能直接索引表(负数索引会回绕到表
        # 内任意槽位,产生假命中 → 逐命中 KeyError);先掩掉再查表。
        # Windows containing non-ACGT bases hash negative; negative indices
        # would wrap around inside the table (spurious hits → KeyError in
        # the hit loop), so mask them out before the table lookup.
        ok = h >= 0
        hs = np.where(ok, h, 0)
        if svals is None:
            cand = np.flatnonzero(ok & tbl[hs])
        else:                          # 大尺度:高 16 位预筛 + 精确校验
            cand = np.flatnonzero(ok & tbl[hs >> 14])
            if cand.size:
                cand = cand[np.isin(h[cand], svals)]
        for p in cand:
            pos = start + int(p)
            kmer = hash2kmer.get(int(h[p]))
            if kmer is None:
                continue
            if chunk[p:p + k] != kmer:
                continue  # 垃圾哈希(含 N/IUPAC 的窗口)误报校验
            run = runs.get(kmer)
            if run is not None and pos - run[1] <= k:
                run[1] = pos             # 延伸当前位点
                continue
            if run is not None:
                flush(kmer, run[0], run[1])
            runs[kmer] = [pos, pos]
        start = end
    for kmer, run in runs.items():
        flush(kmer, run[0], run[1])


def _scan_strand_all(text: str, L: int, scales: list[int],
                     pow4s: dict[int, np.ndarray],
                     hash2kmers: dict[int, dict[int, str]],
                     kmer2qids_all: dict[int, dict[str, list[str]]],
                     lookups: dict[int, tuple],
                     out_by_k: dict[int, list[tuple[int, int, int]]],
                     qid_index: dict[str, int],
                     forward: bool, cancel: CancelFlag | None) -> None:
    """单遍扫描一条链,同时提取全部尺度的 k-mer 命中(R33)。

    原实现每尺度独立 `_scan_strand` 一趟(4 尺度 × 2 链 = 8 趟全库扫描,
    每趟重新编码/查表);此处一遍遍历同时维护各尺度的 numpy 滚动哈希:
    块编码只做一次、revcomp 只做一次、块循环与 run 合并只走一遍 ——
    大库上 4~8 趟的编码/循环开销降到 2 趟(正反链),卷积与查表次数
    不变(每尺度每链各一次,数学上必要)。

    Per-scale rolling hashes are computed from the same encoded block in a
    single pass (the old code scanned the database once per scale × strand,
    re-encoding each time). Encoding, revcomp, chunking and run merging now
    happen once per strand; the per-scale convolutions/table lookups are
    unchanged (mathematically necessary). 4-8 passes become 2.
    """
    n = len(text)
    k_max = max(scales)
    start = 0
    runs: dict[str, list[int]] = {}    # kmer → [run 首命中, 尾命中] (0 基)

    def flush(k: int, kmer: str, first: int, last: int) -> None:
        if forward:
            ss, se = first + 1, last + k
        else:
            ss, se = L - (last + k) + 1, L - first
        for qid in kmer2qids_all[k][kmer]:
            out_by_k[k].append((qid_index[qid], ss, se))

    while start < n - k_max + 1:
        if cancel is not None and cancel.cancelled:
            raise blast.BlastError("任务已取消")
        end = min(start + _CHUNK, n - k_max + 1)
        chunk = text[start:end + k_max - 1]  # 末尾多取 k_max-1 保证最长窗口完整
        arr = _TABLE[np.frombuffer(chunk.encode("ascii"), dtype=np.uint8)]
        w = end - start                      # 有效窗口起点数(每尺度同)
        for k in scales:
            h = _window_hashes(arr, k, pow4s[k])[:w]
            tbl, svals = lookups[k]
            ok = h >= 0
            hs = np.where(ok, h, 0)
            if svals is None:
                cand = np.flatnonzero(ok & tbl[hs])
            else:                          # 大尺度:高 16 位预筛 + 精确校验
                cand = np.flatnonzero(ok & tbl[hs >> 14])
                if cand.size:
                    cand = cand[np.isin(h[cand], svals)]
            hash2kmer = hash2kmers[k]
            for p in cand:
                pos = start + int(p)
                kmer = hash2kmer.get(int(h[p]))
                if kmer is None:
                    continue
                if chunk[p:p + k] != kmer:
                    continue  # 垃圾哈希(含 N/IUPAC 的窗口)误报校验
                run = runs.get(kmer)
                if run is not None and pos - run[1] <= k:
                    run[1] = pos             # 延伸当前位点
                    continue
                if run is not None:
                    flush(k, kmer, run[0], run[1])
                runs[kmer] = [pos, pos]
        start = end
    for kmer, (first, last) in runs.items():
        flush(len(kmer), kmer, first, last)   # 尺度 = k-mer 长度


def count_window_occurrences(
    db_prefix: str,
    template_seq: str,
    windows: dict[str, tuple[int, int]],
    on_log: Callable[[str], None] | None = None,
    cancel: CancelFlag | None = None,
    on_progress: Callable[[float], None] | None = None,
) -> tuple[list[str], dict[int, np.ndarray]]:
    """统计窗口 k-mer 在库中的出现次数。

    `windows`: qseqid → (0 基偏移, 窗长 k)(与管线第一步同构)。返回
    **(qids, per_scale)**:qids = 窗口 qseqid 列表(qidx 下标);per_scale:
    尺度 k → 结构化数组 [("qidx", i4), ("sstart", i4), ("send", i4)],
    每命中仅 12 字节 —— list[dict] 旧格式在 115 Mb 库上 11.3M 命中实测
    3.4 GB(加上下游 per_query 副本峰值 >6 GB,触发换页即"卡死"),数组
    版同量级 ≈ 135 MB(R15.1 修复)。坐标 1-based、正反链统一到正链坐标
    空间;重叠出现按 R11 位点语义**预合并**(同链同 k-mer 起点差 ≤ k 的
    命中并成一条区间,幂等:compute_windowed_depth_compact 再合并得到
    相同位点集合)。无窗口 → ([], {})(不访问库)。

    Count window k-mer occurrences in the database, returning
    **(qids, per_scale)**: qids lists window qseqids (qidx indices);
    per_scale maps scale k to a structured [(qidx, sstart, send)] array at
    12 B per hit — the old list-of-dicts measured 3.4 GB for 11.3M hits on a
    115 Mb db (peaking >6 GB with the downstream per_query copy, swapping
    into a "hang"), the array form is ~135 MB (R15.1 fix). Coordinates are
    1-based with both strands mapped to forward space; occurrences are
    pre-merged into loci with the R11 semantics (same-strand same-k-mer
    hits with start difference <= k collapse into one interval; idempotent —
    compute_windowed_depth_compact re-merging yields the identical locus
    set). An empty window dict returns ([], {}) without touching the db.
    """
    if not windows:
        return [], {}
    if not Path(str(db_prefix) + ".nin").exists():
        raise BlastError("k-mer 计数需要核酸数据库")
    seq = template_seq.upper()

    # 按尺度分组并去重 k-mer → 对应窗口 qseqid 列表(模板内部重复字符串
    # 共享同一命中列表);含 N/IUPAC 的 k-mer 直接丢弃(深度 0)
    # Group by scale, dedupe k-mers to their window qseqids (an in-template
    # repeat shares one hit list); N/IUPAC k-mers are dropped (depth 0).
    off_by_k: dict[int, dict[str, int]] = {}
    for qid, (off, k) in windows.items():
        off_by_k.setdefault(k, {})[qid] = off
    by_k: dict[int, dict[str, list[str]]] = {}
    for k, offs in off_by_k.items():
        kmer2qids: dict[str, list[str]] = {}
        for qid, off in offs.items():
            kmer = seq[off:off + k]
            if len(kmer) < k or any(c not in _ACGT for c in kmer):
                continue
            kmer2qids.setdefault(kmer, []).append(qid)
        if kmer2qids:
            by_k[k] = kmer2qids
    if not by_k:
        return []

    # 一次取全库(blastdbcmd -entry all,与 primer_index 同法);解析为
    # (id, 大写序列) 后立即释放 FASTA 文本与 SeqRecord 容器,序列字符串
    # 只建一次、各尺度复用(此前每尺度对每条序列重复 upper(),还额外持有
    # fa_text 全文)
    # Fetch the whole db once (blastdbcmd -entry all); parse into
    # (id, uppercased seq) and drop the FASTA text and SeqRecord wrappers —
    # the strings are built once and reused across all scales.
    fa_text = blast.fetch_entry(db_prefix, "all")
    entries = [(rec.id, str(rec.seq).upper())
               for rec in SeqIO.parse(io.StringIO(fa_text), "fasta")]
    del fa_text
    if not entries:
        raise BlastError("数据库为空或无法读取条目序列")
    log_step = max(1, len(entries) // 20)

    qids = list(windows)
    qid_index = {q: i for i, q in enumerate(qids)}
    scales = sorted(by_k.keys())
    # 预构建各尺度查找表(单遍扫描共用)
    # Prebuild per-scale lookup tables (shared across the single-pass scans)
    pow4s: dict[int, np.ndarray] = {}
    hash2kmers: dict[int, dict[int, str]] = {}
    lookups: dict[int, tuple] = {}
    for k, kmer2qids in by_k.items():
        pow4s[k] = np.array([4 ** (k - 1 - j) for j in range(k)], dtype=np.int64)
        hash2kmers[k] = {_kmer_hash(km): km for km in kmer2qids}
        vals = np.fromiter(hash2kmers[k].keys(), dtype=np.int64)
        lookups[k] = _membership_lookup(k, vals)
    if on_log:
        on_log(f"k-mer 计数: 单遍扫描 {len(scales)} 个尺度"
               f"({len(entries)} 条序列 ×2 链)...")
    out_by_k: dict[int, list[tuple[int, int, int]]] = {k: [] for k in by_k}
    k_max = max(scales)
    for i, (rec_id, subj) in enumerate(entries, 1):
        if cancel is not None and cancel.cancelled:
            raise blast.BlastError("任务已取消")
        L = len(subj)
        if L >= k_max:
            _scan_strand_all(subj, L, scales, pow4s, hash2kmers, by_k,
                             lookups, out_by_k, qid_index,
                             forward=True, cancel=cancel)
            _scan_strand_all(_revcomp(subj), L, scales, pow4s, hash2kmers,
                             by_k, lookups, out_by_k, qid_index,
                             forward=False, cancel=cancel)
        else:
            # 短条目(不足最长尺度):仅对长度足够的尺度单独扫描
            # Short entries (shorter than the longest scale): scan only the
            # scales their length can host
            for k in scales:
                if L < k:
                    continue
                _scan_strand(subj, L, k, pow4s[k], hash2kmers[k], by_k[k],
                             lookups[k], out_by_k[k], qid_index,
                             forward=True, cancel=cancel)
                _scan_strand(_revcomp(subj), L, k, pow4s[k], hash2kmers[k],
                             by_k[k], lookups[k], out_by_k[k], qid_index,
                             forward=False, cancel=cancel)
        # 进度按条目推进(单遍,每 5% 一次)
        # Progress advances per entry (single pass, once every 5%)
        if on_progress and (i % log_step == 0 or i == len(entries)):
            on_progress(i / len(entries))
        if on_log and i % log_step == 0:
            on_log(f"k-mer 计数: 已扫描 {i}/{len(entries)} 条序列(×2 链)...")
    if cancel is not None and cancel.cancelled:
        raise blast.BlastError("任务已取消")
    if on_log:
        on_log("k-mer 计数: 全部尺度完成")

    # 汇成紧凑结构化数组(每尺度一块;元组暂存峰值 ~1 GB,数组 12 B/命中)
    # Assemble compact structured arrays (one per scale; the tuple staging
    # peaks ~1 GB while the arrays cost 12 B per hit).
    dtype = [("qidx", "<i4"), ("sstart", "<i4"), ("send", "<i4")]
    per_scale = {k: (np.array(v, dtype=dtype) if v else
                     np.empty(0, dtype=dtype))
                 for k, v in out_by_k.items()}
    return qids, per_scale
