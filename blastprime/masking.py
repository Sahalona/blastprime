"""匹配深度分析与渐进式屏蔽掩码(旧实现,保留兼容;新管线不使用——
四段式分级与"可用区两侧缓冲"见 primer_design.py / guide_sup2.md §57)。

深度规则:同一库条目内的重叠 HSP 合并后计 1 次,跨条目累计。
状态:available(单拷贝/已放行)、buffer(重复区边缘 15 bp)、core(不可用)。

Match-depth analysis and progressive masking masks (legacy, kept for
compatibility; the new pipeline uses the profile-threshold staging and the
both-sides available-region buffer instead — see primer_design.py).

Depth rule: overlapping HSPs within the same database entry are merged and
counted once; they accumulate across entries. States: available (single-copy /
released), buffer (15 bp at the edge of repeat regions), core (unusable).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class DepthProfile:
    depth: list[int]          # 逐碱基匹配深度 d(i) (per-base match depth d(i))
    matched: float            # 全序列匹配率 0~1 (overall sequence match fraction, 0~1)
    repeat_frac: float        # 重复区占比 (fraction of repeat regions)
    histogram: dict[int, int] # 深度直方图 (depth histogram)


def compute_depth(hsps: list[dict], template_len: int) -> DepthProfile:
    """由 tabular HSP 列表计算逐碱基匹配深度。

    同一库条目内 subject 区间重叠的 HSP 视为同一位点(同一拷贝),合并后计 1 次;
    不同位点(不同拷贝)跨条目/跨位点累计。

    Compute per-base match depth from a tabular HSP list.

    HSPs with overlapping subject ranges within the same database entry are
    treated as the same locus (same copy), merged and counted once; different
    loci (different copies) accumulate across entries/loci.
    """
    return _depth_from_hsps(hsps, template_len)


def compute_windowed_depth(hsps: list[dict], template_len: int,
                           offsets: dict[str, int]) -> DepthProfile:
    """k-mer 滑动窗口 blastn 的逐碱基出现次数(R11,guide.md §10)。

    `offsets`:窗口 qseqid → 该窗口在模板中的 0 基偏移;只处理出现在映射中
    的 HSP。**每个位置 j 的深度 = 从 j 起始的 k-mer 的出现次数** = 该窗口
    查询的位点数:同窗口内同 subject 重叠(或单碱基间隙)的命中合并为
    1 个位点(内部错配拆 HSP 的情形),跨 subject 累计。调用方必须先过滤
    出全长精确命中(length ≥ k 且 mismatch/gapopen = 0):word_size 7 下
    blastn 的 seed 级/错配延伸命中会把唯一区抬成"重复"。位置 j 只取自己的
    窗口计数 —— 不做跨窗口 q 区间覆盖累计:每个窗口查询的自身位点(模板
    在库中的拷贝)都会覆盖到相邻窗口的位置,若按覆盖累计,深度会被"覆盖
    该位置的窗口数"抬高(模板自身位点存在时 ≈ 窗长数),唯一区也会被相邻
    窗口的偶然命中抬到 2。跨窗口的 subject 位点合并(整段查询语义)则把
    串联重复阵列(单元 60 bp,边界窗口命中跨到相邻单元)并成单个位点,
    深度塌回 1 —— 两者都错,正确语义就是"每个位置的 k-mer 各自数出现
    次数"。

    Per-base k-mer occurrence count from windowed blastn HSPs (R11).

    `offsets` maps each window qseqid to its 0-based template offset; only
    HSPs listed there are used. **Depth at position j = occurrence count of
    the k-mer starting at j** = the locus count of that window query:
    overlapping (or single-gap) subject ranges within one window merge into
    one locus (internal-mismatch HSP splits), distinct subject loci
    accumulate. Position j takes only its own window's count — no cross-
    window q-range coverage accumulation: every window's self-hit (the
    template's own copy in the DB) would cover its neighbours' positions and
    inflate depth by the number of covering windows; in unique regions,
    adjacent chance hits of different windows would likewise raise depth to
    2. Subject-wide merging (whole-query semantics) would collapse a tandem
    array (60 bp unit, boundary windows hitting into the next unit) into a
    single locus with depth 1. The correct semantics is per-position k-mer
    occurrence counting.
    """
    if not offsets:
        return compute_depth(hsps, template_len)
    # 按 (qseqid, sseqid) 分组:每个窗口查询独立数位点(k-mer 出现次数)
    # Group by (qseqid, sseqid): each window query counts its own loci
    per_query: dict[tuple, list[dict]] = {}
    for h in hsps:
        off = offsets.get(h.get("qseqid"))
        if off is None:
            continue
        ss, se = h.get("sstart"), h.get("send")
        if ss is None:
            continue
        ss, se = min(ss, se), max(ss, se)
        per_query.setdefault(h["qseqid"], []).append({"ss": ss, "se": se})
    # 点增量直接是绝对值,不能再做差分数组前缀和(会把每个窗口的计数
    # 泄漏到所有后续位置 —— 深度变成累计值而非"该位置 k-mer 的出现次数")
    # Point increments are absolute; do NOT prefix-sum them (a difference
    # array would leak each window's count into every later position).
    depth = [0] * (template_len + 2)  # 1-based 点位数组 (1-based position array)
    for qseqid, hits in per_query.items():
        loci = _merge_loci(hits)
        # 该窗口 = 从 off 起始的 k-mer;出现次数计入 off 自己的位置
        # This window is the k-mer starting at off; its occurrence count
        # lands on off itself.
        depth[offsets[qseqid] + 1] += len(loci)
    d = depth[1:template_len + 1]

    matched = sum(1 for v in d if v >= 1) / template_len
    repeat = sum(1 for v in d if v >= 2) / template_len
    hist: dict[int, int] = {}
    for v in d:
        hist[v] = hist.get(v, 0) + 1
    return DepthProfile(depth=d, matched=round(matched, 4),
                        repeat_frac=round(repeat, 4), histogram=hist)


def compute_windowed_depth_compact(
    hits: np.ndarray, qids: list[str], template_len: int,
    offsets: dict[str, int],
) -> DepthProfile:
    """`compute_windowed_depth` 的紧凑数组版:kmer_count 计数程序直接返回
    结构化数组 [(qidx, sstart, send)](每命中 12 字节),本函数在 numpy 内
    完成分组 + 线性扫掠合并,与字典版语义逐位一致(位点集合相同)。

    输入 hits 必须只含本尺度窗口的命中(调用方按尺度切分);qidx 是
    qids 的下标,offsets 仍按 qseqid 字符串映射窗口偏移。空数组/None(该
    尺度无命中,如模板全部为含 N 的 k-mer)→ 全零深度剖面 —— 与旧字典版
    空列表 → 全零一致,不要退化为整段 blastn 深度:偶发短命中会把唯一
    模板抬到 2~3。

    Compact-array counterpart of compute_windowed_depth: the counter returns
    a structured [(qidx, sstart, send)] array (12 B per hit) and this
    function groups + merges loci vectorized in numpy, with semantics
    identical to the dict version (identical locus sets). `hits` must
    contain only this scale's windows; qidx indexes `qids`; `offsets` still
    maps qseqid strings to window offsets. An empty array / None (no hits at
    this scale, e.g. all-N k-mers) yields an all-zero profile — matching the
    dict version's empty-list → all-zero behavior; do NOT fall back to the
    whole-query blastn depth, whose chance short HSPs would lift a unique
    template to 2~3.
    """
    depth = [0] * (template_len + 2)
    if hits is not None and len(hits):
        # 按 (qidx, sstart) 排序 → 每个窗口的命中连续,组内 ss 递增
        # Sort by (qidx, sstart): each window's hits become contiguous with
        # ascending ss inside the group.
        order = np.lexsort((hits["sstart"], hits["qidx"]))
        ss = hits["sstart"][order]
        se = hits["send"][order]
        qi = hits["qidx"][order]
        n = len(ss)
        bounds = np.flatnonzero(qi[1:] != qi[:-1]) + 1
        starts = np.concatenate(([0], bounds))
        ends = np.concatenate((bounds, [n]))
        for s, e in zip(starts, ends):
            off = offsets.get(qids[int(qi[s])])
            if off is None:
                continue
            gss, gse = ss[s:e], se[s:e]
            # 线性扫掠(同 _merge_loci):组内每个元素只与**之前**的最大 se
            # 比较 —— maxse[i] = max(se[:i])(gse[0]-2 哨兵只保证首元素
            # 打开新位点,不参与 i ≥ 1 的最大值);新位点计数 = 首元素恒开
            # 1 个 + 其后 ss > 之前最大 se + 1 的元素数。
            # Linear sweep (same rule as _merge_loci): each element compares
            # only against the max se of its predecessors — maxse[i] =
            # max(se[:i]) (the gse[0]-2 sentinel merely opens the first
            # locus; it never participates in i >= 1 maxima). Locus count =
            # 1 for the leading element + elements whose ss > max previous
            # se + 1.
            maxse = np.maximum.accumulate(
                np.concatenate(([gse[0] - 2], gse[:-1])))
            depth[off + 1] += 1 + int(np.count_nonzero(
                gss[1:] > maxse[1:] + 1))
    d = depth[1:template_len + 1]
    matched = sum(1 for v in d if v >= 1) / template_len
    repeat = sum(1 for v in d if v >= 2) / template_len
    hist: dict[int, int] = {}
    for v in d:
        hist[v] = hist.get(v, 0) + 1
    return DepthProfile(depth=d, matched=round(matched, 4),
                        repeat_frac=round(repeat, 4), histogram=hist)


def _merge_loci(hits: list[dict]) -> list[dict]:
    """按 subject 区间重叠(或单碱基间隙)合并为位点(拷贝)。

    线性扫掠(线性扫掠替代原 O(n²) 两两比较):区间按 ss 排序后,合并条件
    (se+1 ≥ ss 且 h.se ≥ loc.ss)退化为 `h.ss ≤ 末位点.se+1`(h.se ≥ h.ss ≥
    loc.ss 恒成立;位点间 se 严格递增,首个匹配位点必是末位点)——这是
    标准区间闭包,结果与逐对合并完全一致。k=8 窗口在真实基因组上可有
    ~2 万孤立命中,旧实现每个命中扫描全部位点 → O(n²) 每窗口,整页任务
    累计数小时(R15.1 修复)。

    Linear sweep (replaces the old O(n²) pairwise scan): after sorting by ss,
    the merge condition (se+1 >= ss and h.se >= loc.ss) reduces to
    `h.ss <= last.se+1` (h.se >= h.ss >= loc.ss always; locus se's strictly
    increase, so the first matching locus is the last) — the standard
    interval closure, identical to pairwise merging. A k=8 window can have
    ~20k isolated hits in a real genome; the old scan was O(n²) per window,
    accumulating to hours per design (R15.1 fix).
    """
    loci: list[dict] = []
    for h in sorted(hits, key=lambda x: (x["ss"], x["se"])):
        if loci and loci[-1]["se"] + 1 >= h["ss"]:
            if h["se"] > loci[-1]["se"]:
                loci[-1]["se"] = h["se"]
        else:
            loci.append({"ss": h["ss"], "se": h["se"]})
    return loci


def _depth_from_hsps(hsps: list[dict], template_len: int) -> DepthProfile:
    # 整段查询语义:按 subject 分组,同条目 subject 重叠 → 同一位点。
    # Whole-query semantics: group by subject; overlapping subject ranges
    # on one entry = one locus.
    per_subject: dict[tuple, list[dict]] = {}
    for h in hsps:
        qs, qe, ss, se = (h.get("qstart"), h.get("qend"),
                          h.get("sstart"), h.get("send"))
        if qs is None or ss is None:
            continue
        qs, qe = max(1, qs), min(template_len, qe)
        if qe < qs:
            qs, qe = qe, qs
        ss, se = min(ss, se), max(ss, se)
        per_subject.setdefault(h.get("sseqid", "?"), []).append(
            {"qs": qs, "qe": qe, "ss": ss, "se": se})

    depth = [0] * (template_len + 2)  # 1-based 差分数组 (1-based difference array)
    for hits in per_subject.values():
        # 按 subject 区间重叠分组 → 位点(拷贝)
        # Group by overlapping subject ranges → loci (copies)
        loci = _merge_loci([{"ss": h["ss"], "se": h["se"]} for h in hits])
        for loc in loci:
            # 组内查询区间按邻接合并,合并后的覆盖区段深度 +1
            # Merge query intervals within the group when adjacent; depth of merged coverage +1
            qints = sorted((h["qs"], h["qe"]) for h in hits
                           if h["ss"] >= loc["ss"] and h["se"] <= loc["se"]
                           and h["ss"] <= loc["se"] and h["se"] >= loc["ss"])
            merged: list[tuple[int, int]] = []
            for a, b in qints:
                if merged and a <= merged[-1][1] + 1:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], b))
                else:
                    merged.append((a, b))
            for a, b in merged:
                depth[a] += 1
                depth[b + 1] -= 1
    for i in range(1, template_len + 1):
        depth[i] += depth[i - 1]
    d = depth[1:template_len + 1]

    matched = sum(1 for v in d if v >= 1) / template_len
    repeat = sum(1 for v in d if v >= 2) / template_len
    hist: dict[int, int] = {}
    for v in d:
        hist[v] = hist.get(v, 0) + 1
    return DepthProfile(depth=d, matched=round(matched, 4),
                        repeat_frac=round(repeat, 4), histogram=hist)


# 掩码状态
# Mask states
AVAILABLE = "available"
BUFFER = "buffer"
CORE = "core"


def build_mask(depth: list[int], released: set[int] | None = None,
               buffer_len: int = 15) -> list[str]:
    """构建设计可用性掩码。

    - available:深度==1 或已放行(released)位置;
    - buffer:位于屏蔽区(深度≥2 且未放行)内、距任一 available 位置 ≤buffer_len;
    - core:其余屏蔽区(不可用)。

    Build the design availability mask.

    - available: positions with depth == 1 or released;
    - buffer: masked positions (depth >= 2 and not released) within buffer_len
      of any available position;
    - core: remaining masked regions (unusable).
    """
    n = len(depth)
    state = ["?"] * n
    released = released or set()
    for i in range(n):
        if depth[i] <= 1 or i in released:
            state[i] = AVAILABLE
        else:
            state[i] = CORE
    # buffer:被屏蔽位置中距最近 available ≤ buffer_len 者
    # buffer: masked positions within buffer_len of the nearest available position
    last_avail = -10 ** 9
    dist_left = [10 ** 9] * n
    for i in range(n):
        if state[i] == AVAILABLE:
            last_avail = i
        dist_left[i] = i - last_avail
    next_avail = 10 ** 9
    for i in range(n - 1, -1, -1):
        if state[i] == AVAILABLE:
            next_avail = i
        if state[i] != AVAILABLE:
            d = min(dist_left[i], next_avail - i)
            if d <= buffer_len:
                state[i] = BUFFER
    return state


def mask_regions(state: list[str], which: str) -> list[tuple[int, int]]:
    """提取状态区间的连续区间 [(start,end)] 1-based 闭区间。

    Extract contiguous runs of the given state as [(start, end)] 1-based closed intervals.
    """
    regions: list[tuple[int, int]] = []
    start = None
    for i, st in enumerate(state):
        if st == which and start is None:
            start = i
        elif st != which and start is not None:
            regions.append((start + 1, i))  # 1-based
            start = None
    if start is not None:
        regions.append((start + 1, len(state)))
    return regions


def release_lowest_depth(depth: list[int], percent: float,
                         already: set[int] | None = None) -> set[int]:
    """按深度升序放行最低的 percent% 位置(覆盖模板长度占比),与已放行累积。

    Release the lowest percent% of positions by ascending depth (as a fraction
    of template length), accumulating with already-released positions.
    """
    n = len(depth)
    target = int(n * percent / 100.0)
    already = set(already or [])
    # 只需放行被屏蔽的(深度≥2 且未放行)
    # Only masked positions (depth >= 2 and not released) need releasing.
    masked = [(depth[i], i) for i in range(n) if depth[i] >= 2 and i not in already]
    masked.sort()
    released = set(already)
    take = min(target, len(masked))
    for _, i in masked[:take]:
        released.add(i)
    return released


def regions_to_excluded(state: list[str]) -> list[tuple[int, int]]:
    """core 区域转 primer3 SEQUENCE_EXCLUDED_REGION 参数 [(start,length)] 1-based。

    primer3 要求 start + length ≤ 序列长度(end 排他),故钳制长度。

    Convert core regions to primer3 SEQUENCE_EXCLUDED_REGION parameters
    [(start, length)] 1-based.

    primer3 requires start + length <= sequence length (end is exclusive),
    so the length is clamped.
    """
    n = len(state)
    out = []
    for s, e in mask_regions(state, CORE):
        length = min(e - s + 1, n - s)
        if length > 0:
            out.append((s, length))
    return out


def target_single_copy_fraction(depth: list[int], target: tuple[int, int] | None) -> float:
    """目标区内单拷贝程度(用于阶段一~三成功时的特异性分折算)。

    Single-copy fraction within the target region (used to scale the
    specificity score when stages 1-3 succeed).
    """
    if target:
        s, e = target
        s, e = max(1, s), min(len(depth), e)
    else:
        s, e = 1, len(depth)
    if e < s:
        return 0.0
    return round(sum(1 for i in range(s - 1, e) if depth[i] <= 1) / (e - s + 1), 4)
