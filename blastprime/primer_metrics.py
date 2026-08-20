"""引物物理指标:primer3-py 封装 + 模块三报告 + 引物对设计任务。

Primer physical metrics: primer3-py wrapper + Module 3 report + primer pair design tasks.
"""

from __future__ import annotations

from typing import Callable

# 从 primer3.bindings 导入:primer3-py 1.x 只在 bindings 下暴露这些函数,
# 2.x 才在顶层再导出 —— 顶层导入在 1.x 环境报
# "module 'primer3' has no attribute 'design_primers'"(exe 建库后设计全灭)
# Import from primer3.bindings: primer3-py 1.x exposes these functions only
# under bindings; 2.x re-exports them at the top level. A top-level import
# fails on 1.x with "module 'primer3' has no attribute 'design_primers'",
# killing every design on such an exe.
# p3_ 前缀别名:本模块自带 design_primers 包装(下方),裸名导入会被覆盖。
# 必须从 primer3.bindings 直接导入(不走顶层):primer3-py 2.3.0 的顶层导出
# 包在 `try: ... except BaseException: pass` 里 —— exe 内 thermoanalysis 因
# p3helpers 缺失初始化失败时,顶层导出被静默吞掉,`import primer3` 成功但
# 没有 design_primers(报 "no attribute");bindings 是底层模块,无此包装。
# The p3_ aliases avoid shadowing: this module defines its own design_primers
# wrapper below, which would override a bare import. Import from
# primer3.bindings directly (not the top level): primer3-py 2.3.0 wraps its
# top-level exports in `try: ... except BaseException: pass` — when
# thermoanalysis fails to init (missing p3helpers) inside the exe, the
# top-level exports are silently dropped, so `import primer3` succeeds
# without design_primers ("no attribute"). bindings has no such wrapper.
from primer3.bindings import (calc_hairpin as p3_calc_hairpin,
                              calc_tm as p3_calc_tm,
                              design_primers as p3_design_primers)

# PyInstaller 打包注意:thermoanalysis(.pyx 编译产物)初始化时 import
# primer3.p3helpers,而 PyInstaller 静态分析读不到 .pyx 源码 → p3helpers
# 漏打包,exe 内 thermoanalysis 加载失败。**不能**用显式 import 让分析器
# 收集它:实测 p3helpers 以扩展模块注册进包会使 PyInstaller 引导器
# longjmp(1.txt/2.txt 归档对比确认唯一差异)。正确做法是打包命令用
# --add-data 把 p3helpers 的 .pyd 以数据文件放进包(见 README 打包命令),
# 运行时由 thermoanalysis 内部的 import 经标准路径搜索加载。
# PyInstaller note: thermoanalysis (a compiled .pyx) imports primer3.p3helpers
# at init, but the analyzer cannot read .pyx sources, so p3helpers is missed.
# Do NOT collect it via an explicit import: bundling p3helpers as a
# registered extension module crashes the PyInstaller bootloader with
# "longjmp" (confirmed by archive diffs). Instead the build command ships
# p3helpers' .pyd as a data file via --add-data (see the README build
# command); thermoanalysis' internal import then finds it through the
# standard path search.

from .config import DEFAULT_PRIMER_PARAMS


# ---------------------------------------------------------------- 基础指标
# ---------------------------------------------------------------- Basic metrics

def tm(seq: str, method: str = "primer3") -> tuple[float, str]:
    """Tm 计算。
    - 短序列(<14 bp):2×(A+T)+4×(G+C) 法则;
    - 较长序列:primer3 盐校正物理模型(注明来源)。

    Tm calculation.
    - Short sequences (<14 bp): the 2×(A+T)+4×(G+C) rule;
    - Longer sequences: primer3 salt-corrected physical model (source noted).
    """
    s = seq.upper()
    if len(s) < 14:
        val = 2.0 * (s.count("A") + s.count("T")) + 4.0 * (s.count("G") + s.count("C"))
        return val, "2×(A+T)+4×(G+C)"
    return round(p3_calc_tm(s), 2), "primer3 盐校正模型"


def gc(seq: str) -> float:
    s = seq.upper()
    if not s:
        return 0.0
    return round(100.0 * (s.count("G") + s.count("C")) / len(s), 1)


def gc_clamp_3p(seq: str, window: int = 5) -> int:
    """3' 端末 window bp 内 G/C 数(GC 夹子)。

    Number of G/C bases within the last window bp at the 3' end (GC clamp).
    """
    tail = seq.upper()[-window:]
    return tail.count("G") + tail.count("C")


def dimer_stats(seq1: str, seq2: str | None = None) -> dict:
    """自互补/异源二聚体:所有位移下最大连续配对长度与最大总配对长度。

    Self-complementary / heterologous dimer: the maximum consecutive paired
    length and the maximum total paired length across all shifts.
    """
    a = seq1.upper()
    b = (seq2 if seq2 is not None else _revcomp(a)).upper()
    best_consec = 0
    best_total = 0
    for shift in range(-len(a) + 1, len(b)):
        start = max(0, -shift)
        end = min(len(a), len(b) - shift)
        if start >= end:
            continue
        run = 0
        total = 0
        for i in range(start, end):
            if _complement(a[i]) == b[i + shift]:
                run += 1
                total += 1
                if run > best_consec:
                    best_consec = run
            else:
                run = 0
    return {"max_consec": best_consec, "max_total": best_total}


def hairpin_stats(seq: str) -> dict:
    """发夹结构:primer3 结果或回文检测。

    Hairpin structure: primer3 result or palindrome detection.
    """
    s = seq.upper()
    try:
        r = p3_calc_hairpin(s)
        delta = float(r.dg)
        stem = int(r.stem_length) if hasattr(r, "stem_length") else 0
        return {"stem_len": stem, "dg": delta}
    except Exception:
        return {"stem_len": 0, "dg": 0.0}


def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def _complement(c: str) -> str:
    return _revcomp(c)


# ---------------------------------------------------------------- 模块三报告
# ---------------------------------------------------------------- Module 3 report

def analyze_short_sequence(seq: str, hsps: list[dict] | None = None) -> dict:
    """模块三:短序列/引物物理指标 + 特异性评估(依赖当前比对 HSP 列表)。

    Module 3: physical metrics for short sequences/primers + specificity
    assessment (depends on the current alignment HSP list).
    """
    s = seq.strip().upper()
    if not s or any(c not in "ACGTUN" for c in s):
        raise ValueError("序列含非法字符")
    tmv, tmethod = tm(s)
    gc_val = gc(s)
    clamp = gc_clamp_3p(s)
    dimer = dimer_stats(s)
    hairpin = hairpin_stats(s)

    items: list[dict] = []
    items.append(_level_item("GC", gc_val, [40, 60], [30, 40, 60, 70], unit="%",
                             good=lambda v: 40 <= v <= 60,
                             ok=lambda v: (30 <= v < 40) or (60 < v <= 70)))
    items.append(_level_item("Tm", tmv, [50, 65], [45, 50, 65, 70], unit="°C",
                             good=lambda v: 50 <= v <= 65,
                             ok=lambda v: (45 <= v < 50) or (65 < v <= 70)))
    items.append(_level_item("3'端GC夹子", clamp, [1, 3], [0, 4], unit=" bp",
                             good=lambda v: 1 <= v <= 3,
                             ok=lambda v: v == 0 or v == 4,
                             bad=lambda v: v >= 5))
    items.append(_level_item("自互补二聚体(最大连续配对)", dimer["max_consec"], [0, 3], [4],
                             unit=" bp",
                             good=lambda v: v <= 3,
                             ok=lambda v: v == 4,
                             bad=lambda v: v >= 5))

    spec = _specificity_assessment(s, hsps or [])
    items.append({
        "key": "specificity", "name": "BLAST 特异性", "value": spec["label"],
        "level": spec["level"], "detail": spec["detail"], "hits": spec["hits"],
    })

    levels = [it["level"] for it in items]
    if "red" in levels:
        verdict = {"level": "red", "label": "不可用", "reason": "存在红色(高风险)指标"}
    elif "yellow" in levels:
        verdict = {"level": "yellow", "label": "风险/待定", "reason": "无红色但有黄色指标"}
    else:
        verdict = {"level": "green", "label": "可用", "reason": "全部指标合格"}
    reasons = [
        ("达标" if it["level"] == "green" else "违规" if it["level"] == "red" else "待定",
         it["name"] + " = " + str(it["value"]))
        for it in items
    ]
    return {
        "seq": s, "tm": tmv, "tm_method": tmethod, "gc": gc_val, "gc_clamp": clamp,
        "dimer": dimer, "hairpin": hairpin,
        "items": items, "verdict": verdict, "reasons": reasons,
    }


def _level_item(key, value, good_range, ok_bounds, unit="", good=None, ok=None, bad=None) -> dict:
    """按 guide.md 6.2 表格给出绿/黄/红等级。

    Assign green/yellow/red levels according to the guide.md 6.2 table.
    """
    if bad and bad(value):
        level = "red"
    elif ok and ok(value):
        level = "yellow"
    elif good(value):
        level = "green"
    else:
        level = "red"
    name_map = {"GC": "GC 含量", "Tm": "Tm 值", "3'端GC夹子": "3' 端 G/C 数(GC 夹子,末 5 bp)",
                "自互补二聚体(最大连续配对)": "自互补二聚体(最大连续配对)"}
    return {"key": key, "name": name_map.get(key, key), "value": value,
            "level": level, "unit": unit}


def _specificity_assessment(seq: str, hsps: list[dict]) -> dict:
    """模块三特异性:首个 HSP 外的其余 HSP 按 L ≥ 8+2k+3×errors 判高危,
    L ≥ max(8,6)+2k+3×errors 判潜在(两者均要求 E-value ≤10)。

    Module 3 specificity: among HSPs other than the first, L ≥ 8+2k+3×errors
    is judged high-risk, and L ≥ max(8,6)+2k+3×errors is judged potential
    (both require E-value ≤ 10).
    """
    if not hsps:
        return {"level": "red", "label": "库中未找到匹配靶标",
                "detail": "当前比对结果中无该查询的命中", "hits": []}
    sorted_h = sorted(hsps, key=lambda h: h.get("evalue", 1e300))
    first = sorted_h[0]
    rest = sorted_h[1:]
    qlen = len(seq)
    details: list[dict] = []
    has_high = has_potential = False
    for h in rest:
        ev = h.get("evalue", 0)
        if ev is None or ev > 10:
            continue
        qend = h.get("qend") or 0
        k = max(0, qlen - qend)           # 3' 端偏移 (3' end offset)
        errors = h.get("mismatch")
        if errors is None:
            ident = str(h.get("identity", "0/0"))
            try:
                n = int(ident.split("/")[0])
                errors = max(0, (h.get("length") or n) - n)
            except Exception:
                errors = 0
        L = h.get("length") or 0
        high_thr = 8 + 2 * k + 3 * errors
        pot_thr = max(8, 6) + 2 * k + 3 * errors
        is_high = L >= high_thr
        is_pot = L >= pot_thr
        if is_high:
            has_high = True
        elif is_pot:
            has_potential = True
        details.append({
            "subject": h.get("sseqid") or h.get("sseqid", ""),
            "evalue": ev, "length": L, "offset_3p": k, "errors": errors,
            "level": "red" if is_high else ("yellow" if is_pot else "green"),
        })
    target = {
        "subject": first.get("sseqid", ""),
        "evalue": first.get("evalue", 0), "length": first.get("length", 0),
        "qstart": first.get("qstart"), "qend": first.get("qend"),
    }
    if len(sorted_h) == 1:
        return {"level": "green", "label": "特异性高",
                "detail": f"仅命中 1 个靶标({target['subject']})", "hits": details, "target": target}
    if has_high:
        return {"level": "red", "label": "高非特异性风险",
                "detail": "存在 3' 端可能结合的脱靶命中(高危)", "hits": details, "target": target}
    if has_potential or any(d["level"] == "yellow" for d in details):
        return {"level": "yellow", "label": "潜在非特异风险",
                "detail": "存在潜在非特异命中", "hits": details, "target": target}
    return {"level": "yellow", "label": "匹配单个基因但存在多处结合区",
            "detail": "命中同一库条目内多处(同一目标内多 HSP)", "hits": details, "target": target}


# ---------------------------------------------------------------- 引物对设计(primer3)
# ---------------------------------------------------------------- Primer pair design (primer3)

def design_primers(
    template: str,
    excluded_regions: list[tuple[int, int]] | None = None,   # (start,length) 1-based
    params: dict | None = None,
    num_return: int = 20,
    single: bool = False,
) -> list[dict]:
    """调用 primer3 设计引物对(或单引物)。返回结构化结果列表。

    Call primer3 to design primer pairs (or a single primer); returns a list
    of structured results.

    注意:不传 SEQUENCE_TARGET —— 该约束强制产物覆盖整个目标区,目标触及模板
    起点或产物上限不足时 primer3 必然 0 候选(定位模式曾因此全灭,见
    primer_design.py 级内循环);目标区由调用方负责命中判定与评分,primer3
    的位置约束只有 SEQUENCE_EXCLUDED_REGION。

    Note: SEQUENCE_TARGET is deliberately not used — it forces the product to
    span the entire target region, yielding zero candidates whenever the
    target touches the template start or exceeds the product range (this once
    broke locate mode; see the level loop in primer_design.py). The target is
    the caller's concern (hit determination and scoring); primer3 receives
    SEQUENCE_EXCLUDED_REGION as its only placement constraint.
    """
    p = {**DEFAULT_PRIMER_PARAMS, **(params or {})}
    seq_args = {
        "SEQUENCE_ID": "template",
        "SEQUENCE_TEMPLATE": template.upper(),
    }
    if excluded_regions:
        # boulder-append 格式:每区一项,不带方括号。起点+长度不得越过模板
        # 末端(primer3 对越界区抛 OSError "EXCLUDED_REGION beyond end of
        # sequence",被下层 catch 吞掉后整体 0 候选);末端被钳掉的位置仍由
        # 级内 3' 后置过滤兜底(剖面=0 → 拒绝)。
        # boulder-append format: one item per region, without brackets.
        # start+length must not exceed the template end (primer3 raises an
        # OSError "EXCLUDED_REGION beyond end of sequence" for out-of-range
        # regions, silently yielding zero candidates); positions dropped by
        # the clamp are still guarded by the per-level 3' post-filter
        # (profile=0 → rejected).
        n = len(template)
        clamped = []
        for s, l in excluded_regions:
            if l <= 0 or s > n:
                continue
            l = min(l, n - s)
            if l > 0:
                clamped.append(f"{s},{l}")
        seq_args["SEQUENCE_EXCLUDED_REGION"] = clamped

    if p.get("product_len_mode") == "relative":
        # 相对产物长度:基准=模板长度 + 偏移(防御层;前端已换算成 absolute,
        # 语义与 config.py product_offset1/2 一致)
        # Relative product length: baseline = template length + offset (defense
        # layer; the frontend has already converted to absolute, consistent with
        # the product_offset1/2 semantics in config.py)
        tlen = len(template)
        o1 = float(p.get("product_offset1", 0.0))
        o2 = float(p.get("product_offset2", 300.0))
        lo = max(50, int(tlen + min(o1, o2)))
        hi = max(tlen, int(tlen + max(o1, o2)))
        prod_range = [[lo, hi]]
    elif p.get("product_len_mode") == "unlimited" or single:
        prod_range = [[50, 100000]]
    else:
        prod_range = [[int(p.get("product_len_min", 150)), int(p.get("product_len_max", 300))]]

    # 产物范围上限钳制到模板长度:模板短于产物下限时,primer3 直接报
    # SEQUENCE_INCLUDED_REGION length < min PRIMER_PRODUCT_SIZE_RANGE 的
    # OSError,钳制后让它自然返回 0 候选 → 管线走 PRIMER3_NO_CANDIDATE
    # Clamp the product range upper bound to the template length: when the
    # template is shorter than the product minimum, primer3 raises an
    # OSError; clamping lets it return 0 candidates naturally and the
    # pipeline falls into PRIMER3_NO_CANDIDATE.
    tlen = len(template)
    prod_range = [[max(1, min(lo, tlen)), max(1, min(hi, tlen))]
                  for lo, hi in prod_range]
    prod_range = [[lo, hi] for lo, hi in prod_range if lo <= hi]

    global_args = {
        "PRIMER_TASK": "generic" if single else "pick_pcr_primers",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        # 单引物模式同样设计反向引物(R35:正向与反向各作为独立单引物项)
        # Single-primer mode also designs reverse primers (R35: forward and
        # reverse each become an independent single-primer item)
        "PRIMER_PICK_RIGHT_PRIMER": 1 if single else 1,
        "PRIMER_NUM_RETURN": max(num_return, 1),
        "PRIMER_PRODUCT_SIZE_RANGE": prod_range or [[1, tlen]],
        "PRIMER_MIN_TM": float(p.get("tm_min", 55)),
        "PRIMER_OPT_TM": float(p.get("tm_opt", 60)),
        "PRIMER_MAX_TM": float(p.get("tm_max", 65)),
        "PRIMER_MIN_GC": float(p.get("gc_min", 30)),
        "PRIMER_MAX_GC": float(p.get("gc_max", 70)),
        "PRIMER_MIN_LENGTH": int(p.get("primer_len_min", 18)),
        "PRIMER_MAX_LENGTH": int(p.get("primer_len_max", 25)),
        "PRIMER_MAX_DIFF_TM": float(p.get("max_tm_diff", 2)),
    }

    try:
        res = p3_design_primers(seq_args, global_args)
    except OSError:
        # 输入校验拒绝(排除区把可用区压缩到产物下限以下等)→ 无候选
        # Input validation rejected (excluded regions shrink the usable
        # region below the product minimum, etc.) → no candidates
        return []
    except Exception as e:
        raise RuntimeError(f"primer3 设计失败: {e}")

    n_left = int(res.get("PRIMER_LEFT_NUM_RETURNED", 0))
    n_right = int(res.get("PRIMER_RIGHT_NUM_RETURNED", 0))
    n_pair = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))
    pairs: list[dict] = []
    for i in range(max(n_left, n_right)):
        def _side(kind: str, idx: int) -> dict:
            pos = res.get(f"PRIMER_{kind}_{idx}")
            if pos is None:
                return {}
            # primer3 报告 0-based 5' 端;统一换算为模板 1-based 坐标
            # (下游 profile 索引、binding-site 映射、定位换算全部按 1-based)
            # primer3 reports the 0-based 5' end; normalize to 1-based
            # template coordinates for all downstream consumers
            return {
                "seq": str(res.get(f"PRIMER_{kind}_{idx}_SEQUENCE", "")),
                "start": int(pos[0]) + 1,
                "len": int(pos[1]),
                "tm": round(float(res.get(f"PRIMER_{kind}_{idx}_TM", 0)), 2),
                "gc": round(float(res.get(f"PRIMER_{kind}_{idx}_GC_PERCENT", 0)), 1),
                "penalty": round(float(res.get(f"PRIMER_{kind}_{idx}_PENALTY", 0)), 2),
                "self_any": float(res.get(f"PRIMER_{kind}_{idx}_SELF_ANY_TH", 0)),
                "self_end": float(res.get(f"PRIMER_{kind}_{idx}_SELF_END_TH", 0)),
                "hairpin": float(res.get(f"PRIMER_{kind}_{idx}_HAIRPIN_TH", 0)),
            }
        left = _side("LEFT", i)
        item = {"left": left,
                "penalty": round(float(res.get(f"PRIMER_PAIR_{i}_PENALTY",
                                               left.get("penalty", 0))), 2),
                "product_len": int(res.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0))}
        if not single:
            item["right"] = _side("RIGHT", i)
        pairs.append(item)
        if single and i < n_right:
            # 单引物模式:反向引物(结合负链)作为独立单引物项,与正向项
            # 同结构(left 承载引物数据),side="reverse" 供前端/导出区分
            # Single-primer mode: the reverse primer (binds the minus strand)
            # becomes an independent single-primer item with the same shape
            # (data in `left`), side="reverse" for the UI/exports
            rev = _side("RIGHT", i)
            pairs.append({
                "left": rev, "side": "reverse",
                "penalty": round(float(res.get(f"PRIMER_RIGHT_{i}_PENALTY",
                                               rev.get("penalty", 0))), 2),
                "product_len": 0,
            })
    return pairs


# ---------------------------------------------------------------- 物理分与评分
# ---------------------------------------------------------------- Physical score and rating

def physical_score(pair: dict, params: dict) -> float:
    """物理分(0~100):primer3 惩罚分归一 + Tm/GC/发夹二聚体加减分。

    Physical score (0-100): normalized primer3 penalty plus Tm/GC/hairpin and
    dimer additions/deductions.
    """
    penalty = pair.get("penalty", 100)
    base = max(0.0, 100.0 - 5.0 * penalty)
    bonus = 0.0
    tm_opt = float(params.get("tm_opt", 60))
    ltm = pair["left"]["tm"]
    rtm = pair.get("right", {}).get("tm", ltm)
    if abs((ltm + rtm) / 2 - tm_opt) <= 2:
        bonus += 5.0
    if all(40 <= pair[k]["gc"] <= 60 for k in ("left", "right") if k in pair):
        bonus += 5.0
    if all(pair[k]["hairpin"] <= -3.0 for k in ("left", "right") if k in pair) or \
       all(pair[k]["self_end"] > -4.0 for k in ("left", "right") if k in pair):
        bonus += 5.0
    return max(0.0, min(100.0, base + bonus))


def composite_score(physical: float, specificity: float, off_target_sites: int = 0,
                    params: dict | None = None) -> float:
    """综合评分 = 物理分×w_phys + 特异性分×w_spec,脱靶位点每个扣分(§39)。

    权重可配置(score_physical_weight / score_specificity_weight,默认各 0.5),
    不再机械使用旧 0.6/0.4。缺省调用保持向后兼容。

    Composite score = physical × w_phys + specificity × w_spec, with points
    deducted per off-target site (§39). Weights are configurable
    (score_physical_weight / score_specificity_weight, default 0.5 each) —
    no longer mechanically 0.6/0.4. Defaults keep backward compatibility.
    """
    w_phys = float((params or {}).get("score_physical_weight", 0.5))
    w_spec = float((params or {}).get("score_specificity_weight", 0.5))
    spec = max(0.0, specificity - 10.0 * off_target_sites)
    wsum = w_phys + w_spec or 1.0
    return round((w_phys * physical + w_spec * spec) / wsum, 1)
