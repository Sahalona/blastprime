/* 引物设计页:输入/参数/四阶段任务/深度图/结果表/详情/导出
   名称型模板(>[gene],range=,database=):后端解析 → 同条目侧翼提取 →
   F/R 设计范围滑块(模板坐标 *N 记法)→ 产物覆盖 → 四阶段设计 → 基因组绝对坐标呈现/导出 */
// Primer design page: input / parameters / four-stage tasks / depth plot / results table / detail / export
// Named template (>[gene], range=, database=): parsed by the backend -> same-entry flank extraction ->
// F/R design-range sliders (template coordinates in *N notation) -> product coverage -> four-stage design -> genomic absolute-coordinate presentation/export
"use strict";

/* 与 config.py DEFAULT_PRIMER_PARAMS 一致的内置默认(恢复默认用) */
// Built-in defaults consistent with config.py DEFAULT_PRIMER_PARAMS (used to restore defaults)
const DEFAULT_PARAMS = {
  mode: "standard", tm_min: 55, tm_opt: 60, tm_max: 65,
  gc_min: 30, gc_max: 70, primer_len_min: 18, primer_len_max: 25,
  product_len_mode: "absolute", product_len_min: 150, product_len_max: 300,
  relative_product_min: 0, relative_product_max: 1, flank_extension: 150,
  product_offset1: 0, product_offset2: 300,
  max_dimer: 5, max_tm_diff: 2, max_gc_clamp_3p: 3, high_risk_threshold: 8,
  level2_global_th: 0.6934, level3_global_th: 0.5503, buffer_len: 8,
  stage4_pool_min: 20, stage4_pool_max: 50, candidate_count: 50,
  offtarget_product_min: 50, offtarget_product_max: 4000,
  target_buffer: 50, blast_evalue: 10, blast_max_targets: 5000,
  timeout_sec: 600, sgrna_len: 20, sgrna_pam: "NGG",
  sgrna_target_only: true, skip_spec_eval: false, skip_kmer_scoring: false,
};

const PARAM_KEYS = [
  "tm_min", "tm_opt", "tm_max", "gc_min", "gc_max", "primer_len_min", "primer_len_max",
  "product_len_mode", "product_len_min", "product_len_max",
  "product_offset1", "product_offset2",
  "max_dimer", "max_gc_clamp_3p",
  "buffer_len", "level2_global_th", "level3_global_th", "target_buffer",
  "offtarget_product_min", "offtarget_product_max", "blast_evalue", "blast_max_targets",
  "timeout_sec", "sgrna_len", "sgrna_pam", "stage4_pool_max", "candidate_count", "sgrna_target_only",
];

/* CSS 变量 → 实际颜色值。canvas 2D 不支持 "var(--xxx)" 字符串(直接赋值会解析失败
   全部渲染为默认黑色),必须经 getComputedStyle 取值;jsdom 下取不到时回退默认色 */
// CSS variables -> actual color values. Canvas 2D does not support "var(--xxx)" strings (direct assignment fails to parse
// and everything renders as the default black); values must be read via getComputedStyle; fall back to the default color when unavailable under jsdom
function cssVar(name, fb) {
  try {
    const v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
    return v || fb;
  } catch (e) { return fb; }
}

// hex 颜色加深(负 factor 变浅):供"选中范围"在恒长轨道上作对比色
// Darken a hex color (negative factor lightens it): used as the contrast color for the selected range over the constant-length track
function shadeColor(hex, factor) {
  try {
    let h = String(hex).trim();
    if (h.startsWith("#")) h = h.slice(1);
    if (h.length === 3) h = h.split("").map((c) => c + c).join("");
    if (!/^[0-9a-fA-F]{6}$/.test(h)) return hex;
    const f = Math.max(-100, Math.min(100, factor));
    const n = parseInt(h, 16);
    const ch = (v) => {
      const c = v + Math.round(f * 2.55);
      return Math.max(0, Math.min(255, c)).toString(16).padStart(2, "0");
    };
    return "#" + ch((n >> 16) & 255) + ch((n >> 8) & 255) + ch(n & 255);
  } catch (e) { return hex; }
}

/* 双端滑块:轴坐标 → 模板坐标;F 轴 0=目标区间起点,R 轴 seq_len=目标区间终点(guide_sup1.md §6)。
   显示记法:*N 表示距端点 N bp(上游/下游),*0=恰在端点,正数=模板内坐标。
   两手柄最小间距 12px:左右值相等在滑块层不可达(用户要求);最近拖动者画上层+加粗描边,
   重叠时也能区分拖的是哪个 */
// Dual-handle slider: axis coordinates -> template coordinates; F axis 0 = target-region start, R axis seq_len = target-region end (guide_sup1.md §6).
// Display notation: *N means N bp from the endpoint (upstream/downstream), *0 = exactly at the endpoint, positive numbers = in-template coordinates.
// Minimum handle gap of 12px: equal left/right values are unreachable at the slider level (user requirement); the most recently dragged handle is drawn on top with a thicker outline,
// so it remains distinguishable which handle is being dragged when they overlap
class DualSlider {
  constructor(canvas, opts) {
    this.cv = canvas;
    this.ctx = canvas.getContext("2d");
    this.minV = opts.minV; this.maxV = opts.maxV;
    this.axis = opts.axis;      // "f" | "r"(决定 *N 记法)
    // "f" | "r" (determines the *N notation)
    this.seqLen = opts.seqLen;  // 目标区间长度(R 轴端点基准)
    // Target-region length (endpoint reference for the R axis)
    this.onSet = opts.onSet || (() => {});
    this.start = opts.start; this.end = opts.end;
    // 构造路径同样执行最小间距钳制(与 set() 一致):轴被目标区间边界压缩时
    // 默认位置可能退化为两值相等(F 轴起点贴条目起点)甚至越序(R 轴终点贴条目末尾)
    // The constructor path applies the same minimum-gap clamp (consistent with set()): when the axis is compressed by the target-region bounds,
    // the default positions may degenerate to equal values (F-axis start stuck at the entry start) or even out of order (R-axis end stuck at the entry end)
    const gap = this.gapUnits();
    if (this.end - this.start < gap) {
      this.end = Math.min(this.maxV, this.start + gap);
      if (this.end - this.start < gap) this.start = Math.max(this.minV, this.end - gap);
    }
    this._drag = null;     // {which: "start"|"end", px: x}
    this._lastMoved = null; // 最近拖动的手柄(重叠时视觉区分)
    // The most recently dragged handle (for visual distinction when overlapping)
    this._bind();
    this.render();
  }
  _xOf(v) { const W = this.cv.width; return 6 + (v - this.minV) / (this.maxV - this.minV || 1) * (W - 12); }
  _vOf(x) { const W = this.cv.width; return Math.round(this.minV + (x - 6) / (W - 12) * (this.maxV - this.minV || 1)); }
  /* 最小间距(轴单位):对应 12px 屏幕距离 */
  // Minimum gap (in axis units): corresponds to a 12px screen distance
  gapUnits() { const W = this.cv.width; return Math.max(1, Math.round(12 * (this.maxV - this.minV) / Math.max(1, W - 12))); }
  _bind() {
    const pos = (e) => {
      const r = this.cv.getBoundingClientRect();
      return (e.clientX - r.left) * (this.cv.width / r.width);
    };
    this._onDown = (e) => {
      const x = pos(e);
      const ds = Math.abs(x - this._xOf(this.start)), de = Math.abs(x - this._xOf(this.end));
      const which = ds <= de ? "start" : "end";
      this._drag = { which, px: x };
      this._lastMoved = which;
      e.preventDefault();
    };
    this._onMove = (e) => {
      if (!this._drag) return;
      const v = this._vOf(pos(e));
      this.set(this._drag.which === "start" ? v : this.start, this._drag.which === "end" ? v : this.end, this._drag.which);
    };
    this._onUp = () => { this._drag = null; };
    this.cv.addEventListener("mousedown", this._onDown);
    window.addEventListener("mousemove", this._onMove);
    window.addEventListener("mouseup", this._onUp);
  }
  /* 解绑全部事件(rebuildSliders 重建时调用):旧实例的 canvas mousedown
     与 window mousemove 若不移除,拖动会同时驱动新旧两个实例——旧实例
     (越界值)最后渲染覆盖新实例,画布横跳且渲染与数值框分离 */
  // Unbind all listeners (called when rebuildSliders recreates the sliders):
  // without this, dragging drives both the old and the new instance — the
  // old one (stale values) renders last, causing the canvas to jump between
  // positions while the number boxes stay correct
  unbind() {
    if (!this._onDown) return;
    this.cv.removeEventListener("mousedown", this._onDown);
    window.removeEventListener("mousemove", this._onMove);
    window.removeEventListener("mouseup", this._onUp);
    this._drag = null;
  }
  set(start, end, which) {
    start = Math.max(this.minV, Math.min(this.maxV, start));
    end = Math.max(this.minV, Math.min(this.maxV, end));
    // 最小间距钳制:相等或越序都不允许(which 指定拖动哪端,另一端保持)
    // Minimum-gap clamp: equal or out-of-order values are not allowed (which specifies the dragged end; the other end stays)
    const gap = this.gapUnits();
    if (end - start < gap) {
      if (which === "start") start = end - gap;
      else if (which === "end") end = start + gap;
      else { end = start + gap; if (end > this.maxV) { end = this.maxV; start = end - gap; } }
    }
    if (start === this.start && end === this.end) return;
    if (which) this._lastMoved = which;
    this.start = start; this.end = end;
    this.render();
    this.onSet(start, end);
  }
  /* 显示记法:轴值 → *N 形式;F 轴 0 为端点,R 轴 seqLen 为端点 */
  // Display notation: axis value -> *N form; on the F axis 0 is the endpoint, on the R axis seqLen is the endpoint
  static fmt(v, axis, seqLen) {
    if (axis === "f") {
      if (v <= 0) return v === 0 ? "*0" : `*${-v}`;
      return String(v);
    }
    if (v >= seqLen) return v === seqLen ? "*0" : `*${v - seqLen}`;
    return String(v);
  }
  render() {
    const ctx = this.ctx, W = this.cv.width, H = this.cv.height;
    ctx.clearRect(0, 0, W, H);
    const y = H / 2;
    const accent = cssVar("--accent", "#2563eb");
    const outside = cssVar("--slider-outside", "#4b5563");   // 范围之外:深灰(深浅主题各异)
    const between = cssVar("--slider-between", "#93c5fd");   // 范围之间:浅蓝
    const fg = cssVar("--fg", "#1c2330");
    const bg = cssVar("--bg", "#ffffff");
    // 轨道:全段 = 范围之外,深灰绘制,进度条恒等长(用户要求:不随可选总数
    // 忽长忽短);手柄之间用浅蓝加粗覆盖,标出当前可设计范围
    // Track: the full length is drawn in dark gray (out-of-range), so the bar
    // has a constant length (user request: no more fluctuating with the
    // selectable range); the region between the handles is overlaid with a
    // thicker light-blue band to mark the currently selectable range
    ctx.lineCap = "round";
    ctx.strokeStyle = outside;
    ctx.lineWidth = 4;
    ctx.beginPath(); ctx.moveTo(this._xOf(this.minV), y); ctx.lineTo(this._xOf(this.maxV), y); ctx.stroke();
    ctx.strokeStyle = between;
    ctx.lineWidth = 6;
    ctx.beginPath(); ctx.moveTo(this._xOf(this.start), y); ctx.lineTo(this._xOf(this.end), y); ctx.stroke();
    // 端点刻度值
    // Endpoint tick values
    ctx.fillStyle = fg;
    ctx.font = "9px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText(String(this.minV), Math.max(10, this._xOf(this.minV)), H - 3);
    ctx.fillText(String(this.maxV), Math.min(W - 10, this._xOf(this.maxV)), H - 3);
    // 值标签:向两端外侧偏移,避免两手柄接近时文字重叠
    // Value labels: offset outward at both ends to avoid text overlap when the two handles are close
    ctx.font = "10px monospace";
    ctx.fillText(DualSlider.fmt(this.start, this.axis, this.seqLen),
      Math.max(14, this._xOf(this.start) - 14), 8);
    ctx.fillText(DualSlider.fmt(this.end, this.axis, this.seqLen),
      Math.min(W - 14, this._xOf(this.end) + 14), 8);
    // 手柄:未动的先画(下层),最近拖动的后画(上层+加粗描边);accent 填充 + fg 描边 + bg 内点
    // Handles: the unmoved one is drawn first (bottom layer), the most recently dragged one last (top layer + thicker outline); accent fill + fg stroke + bg inner dot
    const order = this._lastMoved === "start" ? ["end", "start"] : ["start", "end"];
    order.forEach((w) => {
      const v = w === "start" ? this.start : this.end;
      const x = this._xOf(v);
      ctx.fillStyle = accent;
      ctx.strokeStyle = fg;
      ctx.lineWidth = w === this._lastMoved ? 2.5 : 1.5;
      ctx.beginPath(); ctx.arc(x, y, 5, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
      ctx.fillStyle = bg;
      ctx.beginPath(); ctx.arc(x, y, 2, 0, Math.PI * 2); ctx.fill();
    });
  }
}

window.PAGE = {
  taskId: null,
  result: null,
  resultIdx: 0,
  selPair: null,
  canvasRects: [],
  sliders: {},
  /* 名称型模板状态:解析成功后 tmpl 就绪,滑块/范围输入框可用 */
  // Named-template state: once parsing succeeds tmpl is ready and the sliders/range inputs are usable
  named: { query: "", flankLen: 150, tmpl: null, f: null, r: null },

  async init() {
    await this.loadDbs();
    this.fillParams(App.settings?.primer_params || DEFAULT_PARAMS);
    this.updateProductFields();
    this.restoreState();
    this.updateParamsCaret();   // restoreState 提前返回时也要同步 caret
    this.restoreResult();
    // 兜底持久化完成事件:切页回来时任务在后台完成、init 的同步 restoreResult
    // 早于异步写入 → 事件到达时再恢复渲染(仅当前无结果时,避免与 run()
    // 的结果渲染重复)
    // Fallback-persist completion event: when a task finished while the page
    // was away, init's synchronous restoreResult ran before the async write
    // — the event triggers the recovery render (only when no result exists,
    // so run()'s render never duplicates)
    window.addEventListener("bp:result-persisted", (e) => {
      const d = e && e.detail;
      if (d && d.kind === "design" && !this.result) this.restoreResult();
    });
    this.bindStateSavers();
    document.getElementById("btn-design").addEventListener("click", () => this.run());
    // 模板序列:从 FASTA 文件导入(追加模式,多次导入会累积)
    // Template sequence: imported from a FASTA file (append mode; repeated imports accumulate)
    document.getElementById("dp-file").addEventListener("change", async (e) => {
      const { added, failed } = await appendFilesToTextarea(e.target, document.getElementById("dp-template"));
      if (added > 0) toast(t("design.import_added").replace("{n}", added));
      if (failed > 0) toast(t("common.read_failed"));
    });
    document.getElementById("btn-restore").addEventListener("click", async () => {
      this.fillParams(DEFAULT_PARAMS);
      this.updateProductFields();
      try {
        await api("/api/settings", { method: "POST", json: true, body: { primer_params: this.collectParams() } });
        toast(t("settings.saved"));
      } catch (e) { toast(e.message); }
    });
    document.getElementById("params-toggle").addEventListener("click", () => {
      const body = document.getElementById("params-body");
      body.style.display = body.style.display === "none" ? "" : "none";
      this.updateParamsCaret();
    });
    document.getElementById("btn-export-csv").addEventListener("click", () => {
      if (this.taskId) downloadURL(`/api/export/design/${this.taskId}/stats.csv`, "primers.csv");
    });
    document.getElementById("btn-export-fasta").addEventListener("click", () => {
      if (this.taskId) downloadURL(`/api/export/design/${this.taskId}/pairs.fasta`, "primers.fasta");
    });
    /* 复制 .fasta 导出内容到剪贴板(与导出按钮同一后端端点,内容一致) */
    // Copy the .fasta export content to the clipboard (same backend endpoint as the export button)
    document.getElementById("btn-copy-fasta").addEventListener("click", async () => {
      if (!this.taskId) { toast(t("design.export_no_task")); return; }
      try {
        const resp = await fetch(`/api/export/design/${this.taskId}/pairs.fasta`);
        if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
        await copyText(await resp.text());
      } catch (e) {
        toast(t("common.copy_failed"));
      }
    });
    window.addEventListener("resize", () => {
      if (this.result) this.renderDepth();
      // 窗口宽度变化 → 滑块轨道按新容器宽度重排(防抖 120ms)
      // Window width change -> re-fit slider tracks to the new container width (debounced 120ms)
      clearTimeout(this._rszT);
      this._rszT = setTimeout(() => this.fitSliders(), 120);
    });
    /* 产物长度模式/设计模式 → 条件显示产物最短/最长/偏移 */
    // Product-length mode / design mode -> conditionally show the product min/max/offset fields
    document.getElementById("p-product_len_mode").addEventListener("change", () => {
      this.updateProductFields();
      this.saveState();
    });
    document.getElementById("dp-mode").addEventListener("change", () => {
      this.updateProductFields();
      this.saveState();
    });
    /* k-mer 缓存:清空按钮 + 条数(语言切换时计数文案重绘) */
    // K-mer cache: clear button + count (re-rendered on language switch)
    document.getElementById("btn-clear-kmer-cache").addEventListener("click", () => {
      this.clearKmerCaches();
    });
    this.updateKmerCacheCount();
    /* 名称型模板:模板文本/库/侧翼长度变化 → 重新解析(防抖 500ms) */
    // Named template: changes to template text / database / flank length -> re-parse (debounced 500ms)
    let namedTimer = null;
    document.getElementById("dp-template").addEventListener("input", () => {
      clearTimeout(namedTimer);
      namedTimer = setTimeout(() => this.loadNamedTemplate(), 500);
    });
    document.getElementById("dp-flank").addEventListener("change", () => {
      this.named.flankLen = Math.max(0, parseInt(document.getElementById("dp-flank").value) || 0);
      if (this.named.tmpl) this.loadNamedTemplate();   // 侧翼变化 → 重取模板
      // Flank change -> re-fetch the template
      this.saveState();
    });
    document.getElementById("dp-db").addEventListener("change", () => {
      if (this.named.tmpl || this.isNameQuery(document.getElementById("dp-template").value)) {
        this.loadNamedTemplate();   // 换库 → 名称型解析需重跑
        // Database change -> the named-query parse needs to re-run
      }
    });
    /* 范围输入框 ←→ 滑块双向同步 */
    // Two-way sync between the range input fields and the sliders
    ["p-f-start", "p-f-end", "p-r-start", "p-r-end"].forEach((id) => {
      document.getElementById(id).addEventListener("change", (e) => this.onRangeInput(e.target));
    });
    this.consumeLoadedProject();
  },

  /* ---------------- 页面状态临时保存(切页复原,见 app.js savePageState) ---------------- */
  // ---------------- Page state temporary save (restored on page switch-back, see app.js savePageState) ----------------

  saveState() {
    const params = {};
    PARAM_KEYS.forEach((k) => {
      const el = document.getElementById("p-" + k);
      if (el) params[k] = el.value;
    });
    savePageState("design", {
      template: document.getElementById("dp-template").value,
      db: document.getElementById("dp-db").value,
      mode: document.getElementById("dp-mode").value,
      flank: document.getElementById("dp-flank").value,
      skipSpec: document.getElementById("dp-skip-spec").checked,
      skipKmer: document.getElementById("dp-skip-kmer").checked,
      params,
      named: this.serializeNamed(),
    });
  },

  restoreState() {
    const s = loadPageState("design");
    if (!s) return;
    const set = (id, v) => { if (v != null) document.getElementById(id).value = v; };
    set("dp-template", s.template);
    set("dp-mode", s.mode);
    set("dp-flank", s.flank);
    document.getElementById("dp-skip-spec").checked = !!s.skipSpec;
    document.getElementById("dp-skip-kmer").checked = !!s.skipKmer;
    // 库下拉仅恢复仍存在的库
    // Restore only databases that still exist in the dropdown
    const sel = document.getElementById("dp-db");
    if (s.db && [...sel.options].some((o) => o.value === s.db)) sel.value = s.db;
    if (s.params) {
      PARAM_KEYS.forEach((k) => {
        const el = document.getElementById("p-" + k);
        if (el && s.params[k] != null) el.value = s.params[k];
      });
    }
    // 设计参数区默认展开,不恢复折叠状态
    // Design params section always starts expanded
    this.updateParamsCaret();
    this.updateProductFields();
    if (s.named) this.restoreNamed(s.named);
  },

  bindStateSavers() {
    const save = () => this.saveState();
    document.getElementById("dp-template").addEventListener("input", save);
    ["dp-db", "dp-mode", "dp-flank", "dp-skip-spec", "dp-skip-kmer"].forEach((id) => document.getElementById(id).addEventListener("change", save));
    // 参数区:委托监听全部 p-* 输入/选择
    // Params section: delegated listener for all p-* inputs/selects
    const body = document.getElementById("params-body");
    body.addEventListener("input", save);
    body.addEventListener("change", save);
    document.getElementById("params-toggle").addEventListener("click", save);
  },

  async loadDbs() {
    const sel = document.getElementById("dp-db");
    try {
      const { records } = await api("/api/db/records");
      if (!records.length) return;
      // 优先展示核酸库
      // Prioritize showing nucleotide databases
      const sorted = records.sort((a, b) => {
        const na = /\.nsq|\.nin|\.n/.test(a.prefix) ? 0 : 1;
        const nb = /\.nsq|\.nin|\.n/.test(b.prefix) ? 0 : 1;
        return na - nb;
      });
      sel.innerHTML = sorted.map((r) => {
        const name = r.prefix.split("/").pop();
        const label = r.note ? `${name}（${r.note}）` : name;
        return `<option value="${escapeHtml(r.prefix)}" title="${escapeHtml(r.prefix)}">${escapeHtml(label)}</option>`;
      }).join("");
    } catch (e) { /* 无库 */ }
    // No databases
  },

  /* ---------------- 参数表单 ---------------- */
  // ---------------- Parameter form ----------------

  fillParams(p) {
    PARAM_KEYS.forEach((k) => {
      const el = document.getElementById("p-" + k);
      if (!el) return;
      if (k === "sgrna_target_only") el.checked = !!p[k];
      else if (p[k] !== undefined) el.value = p[k];
    });
    const cb = document.getElementById("dp-skip-spec");
    if (cb && p.skip_spec_eval !== undefined) cb.checked = !!p.skip_spec_eval;
    const kb = document.getElementById("dp-skip-kmer");
    if (kb && p.skip_kmer_scoring !== undefined) kb.checked = !!p.skip_kmer_scoring;
    this.updateProductFields();
  },

  collectParams() {
    const p = { ...DEFAULT_PARAMS };
    PARAM_KEYS.forEach((k) => {
      const el = document.getElementById("p-" + k);
      if (!el) return;
      if (k === "sgrna_target_only") p[k] = el.checked;
      else if (k === "sgrna_pam" || k === "product_len_mode") p[k] = el.value;
      else p[k] = el.value === "" ? p[k] : Number(el.value);
    });
    p.mode = document.getElementById("dp-mode").value;
    if (p.mode === "single") p.product_len_mode = "unlimited";
    p.skip_spec_eval = document.getElementById("dp-skip-spec").checked;
    p.skip_kmer_scoring = document.getElementById("dp-skip-kmer").checked;
    return p;
  },

  /* 产物条件显示:绝对值→最短/最长(同一行,各占一半宽);相对→偏移两框;
     不限制/single→三者全隐 */
  // Conditional product display: absolute -> min/max on one row (half width
  // each); relative -> the two offset fields; unlimited/single -> all hidden
  updateProductFields() {
    const mode = document.getElementById("dp-mode").value;
    const eff = mode === "single" ? "unlimited" : document.getElementById("p-product_len_mode").value;
    const pw = document.getElementById("p-product_len-wrap");
    const ow = document.getElementById("p-product-offset-wrap");
    if (pw) pw.style.display = eff === "absolute" ? "" : "none";
    if (ow) ow.style.display = eff === "relative" ? "" : "none";
  },

  /* ================= 名称型模板(guide_sup1.md §5/§6 语义) ================= */
  // ================= Named template (guide_sup1.md §5/§6 semantics) =================

  /* 名称型格式:首行 >[gene],range=...,database=...,name=...,targetbase=...
     除 gene 外全部可省略。判定与后端同 gate:首行含逗号(key=value 语法),
     或首行后无任何序列字符(裸 >gene 即名称查询);普通 FASTA 不受影响。
     后端 named:false 兜底(解析失败时按普通模板处理)。 */
  // Named-query format: first line >[gene],range=...,database=...,name=...,
  // targetbase=... — everything but the gene is optional. Same gate as the
  // backend: the first line has a comma (key=value syntax), or no sequence
  // characters follow it (a bare >gene is a name query). Plain FASTA is
  // unaffected; the backend's named:false is the fallback.
  isNameQuery(text) {
    const lines = (text || "").split(/\r?\n/);
    const first = (lines[0] || "").trim();
    if (!first.startsWith(">")) return false;
    if (first.includes(",")) return true;
    return !/[ACGTUNacgtun]/.test(lines.slice(1).join(""));
  },

  /* 可序列化摘要:query + flank + 提取模板 + 当前 F/R 范围 */
  // Serializable summary: query + flank + extracted template + current F/R ranges
  serializeNamed() {
    const n = this.named;
    if (!n.tmpl) return null;
    // v:2 = R 轴坐标帧 v2(0 = 目标片段起点,与 F 轴同原点);旧帧(v 缺省)的
    // R 数值语义完全不同,restoreNamed 见版本不符即丢弃范围回默认
    // v:2 = R-axis frame v2 (0 = target-region start, same origin as F); old-frame
    // (v absent) R numbers mean something else entirely — restoreNamed drops them
    return { v: 2, query: n.query, flankLen: n.flankLen, tmpl: n.tmpl, f: n.f, r: n.r };
  },

  restoreNamed(s) {
    if (!s || !s.tmpl) return;
    // 旧版/损坏的命名模板状态缺少核心字段:忽略而非抛异常,避免 init 中断导致事件绑定全部丢失
    // Old/corrupted named-template state missing core fields: ignore instead of throwing, to avoid init interruption that would lose all event bindings
    const t = s.tmpl;
    if (!t.entry || !t.target || !(t.target.start > 0) || !(t.target.end > 0) || !(t.template_len > 0)) return;
    this.named.query = s.query || "";
    this.named.flankLen = s.flankLen != null ? s.flankLen : 150;
    document.getElementById("dp-flank").value = this.named.flankLen;
    this.named.tmpl = s.tmpl;
    // 坐标帧版本:旧版(v 缺省)R 轴以 template_len 为基准,数值含义与现帧
    // 完全不同——版本不符即丢弃保存的范围回默认(*0~*flank),查询/侧翼/
    // 模板照常恢复
    // Coordinate-frame version: old states (v absent) used the template_len-based
    // R axis, whose numbers mean something entirely different now — drop the saved
    // ranges back to defaults; query/flank/template still restore
    const frameOk = s.v === 2;
    this.named.f = frameOk && Array.isArray(s.f) ? s.f.slice() : null;
    this.named.r = frameOk && Array.isArray(s.r) ? s.r.slice() : null;
    // 先显示再建滑块(位图宽度依赖可见布局)
    // Show first, then build the sliders (bitmap width depends on visible layout)
    document.getElementById("flank-view").style.display = "";
    this.rebuildSliders(s.tmpl);
    // 回填用户保存的范围(程序化设置,不记 _lastMoved;越界即回默认——
    // 模板/范围随会话变化时旧值可能出界)
    // Backfill the ranges saved by the user (programmatic set, no _lastMoved;
    // out-of-range values fall back to defaults — the template/ranges may
    // change between sessions)
    if (this.named.f && this.named.f[0] >= this.sliders.f.minV && this.named.f[1] <= this.sliders.f.maxV)
      this.sliders.f.set(this.named.f[0], this.named.f[1], null);
    if (this.named.r && this.named.r[0] >= this.sliders.r.minV && this.named.r[1] <= this.sliders.r.maxV)
      this.sliders.r.set(this.named.r[0], this.named.r[1], null);
    this.syncRangeInputs();
    this.renderTmplInfo();
  },

  /* 名称型模板解析:后端按 [range−flank, range+flank] 同条目截取。
     防抖去重:同一 (文本,库,侧翼) 已解析成功则跳过 */
  // Named-template parsing: the backend extracts the same entry over [range-flank, range+flank].
  // Debounce/dedupe: skip if the same (text, database, flank) has already been parsed successfully
  async loadNamedTemplate() {
    const text = document.getElementById("dp-template").value.trim();
    if (!this.isNameQuery(text)) {
      this._namedKey = null;
      this.named = { query: "", flankLen: this.named.flankLen, tmpl: null, f: null, r: null };
      document.getElementById("flank-view").style.display = "none";
      ["p-f-start", "p-f-end", "p-r-start", "p-r-end"].forEach((id) => {
        const el = document.getElementById(id);
        el.value = ""; el.disabled = true;
      });
      return;
    }
    const db = document.getElementById("dp-db").value;
    // 名称型自带 database=/targetbase= 时页面未选库也能解析(后端从历史
    // 记录解析);否则库下拉可能尚未填充完成,跳过
    // With database=/targetbase= in the query the page db is optional (the
    // backend resolves from history records); otherwise the dropdown may not
    // be populated yet — skip.
    if (!db && !/,\s*(database|targetbase)=/i.test(text)) return;
    const flank = Math.max(0, parseInt(document.getElementById("dp-flank").value) || 0);
    const key = `${text}|${db}|${flank}`;
    if (this._namedKey === key && this.named.tmpl) return;
    this._namedKey = key;
    // 解析中提示(名称型解析要查库,大库上可达数秒)
    // Parse-in-progress hint (named resolution queries the DB; can take seconds on large DBs)
    const ns = document.getElementById("named-status");
    if (ns) {
      ns.style.display = "";
      ns.innerHTML = `<span class="spinner"></span>${t("design.named_parsing")}`;
    }
    try {
      const resp = await api("/api/design/template", {
        method: "POST", json: true,
        body: { template_text: text, db_prefix: db, flank_len: flank },
      });
      if (!resp.named) {
        this._namedKey = null;
        document.getElementById("flank-view").style.display = "none";
        return;
      }
      const prevF = this.named.f, prevR = this.named.r;   // 修改前的 F/R(可能 null)
      this.named.query = text;
      this.named.flankLen = flank;
      this.named.tmpl = resp;
      // 注意:不得在此置 this.named.f/r 为 null——set() 的 onSet → syncRangeInputs
      // 会读 named.f/r(此时可能为 null 而抛 TypeError,异常被 catch 后 R 的
      // set 从未执行,新解析的 F/R 全部回默认且模板卡被误隐藏,R22 已踩过)
      // Never null out named.f/r here: set()'s onSet → syncRangeInputs reads
      // them (a null throws, the exception is swallowed by catch, the R-axis
      // set never runs, and the fresh F/R all fall back to defaults)
      // 先显示卡片再建滑块:位图宽度 = 容器可用宽度,display:none 时量不到
      // Show the card first, then build the sliders: bitmap width = container width, which is 0 while display:none
      document.getElementById("flank-view").style.display = "";
      this.rebuildSliders(resp);
      // 修改名称型输入时保留用户之前设置的 F/R:新模板滑块范围内则恢复,
      // 越界(目标区间/模板变化)回默认——避免每次修改输入都被重置成默认
      // Preserve the user's previous F/R when the named input changes: restore
      // if still within the new slider bounds, fall back to defaults when out
      // of range (the target/template changed)
      if (prevF && prevF[0] >= this.sliders.f.minV && prevF[1] <= this.sliders.f.maxV)
        this.sliders.f.set(prevF[0], prevF[1], null);
      if (prevR && prevR[0] >= this.sliders.r.minV && prevR[1] <= this.sliders.r.maxV)
        this.sliders.r.set(prevR[0], prevR[1], null);
      // 同步最终滑块位置到 named.f/r(rebuildSliders 不再代做)
      this.named.f = [this.sliders.f.start, this.sliders.f.end];
      this.named.r = [this.sliders.r.start, this.sliders.r.end];
      this.syncRangeInputs();
      this.renderTmplInfo();
      this.saveState();
    } catch (e) {
      this._namedKey = null;
      document.getElementById("flank-view").style.display = "none";
      toast(e.message || t("design.named_failed"));
    } finally {
      const ns = document.getElementById("named-status");
      if (ns) ns.style.display = "none";
    }
  },

  /* 滑块:参照目标区间(§6)。F 轴 0=区间起点、R 轴 seq_len=区间终点。
     F [start,end] → 模板 1-based [a+start+1, a+end+1](a=target.start−1)
     R [start,end] → 模板 1-based [a+start+1, a+end+1](两轴同原点,a=target.start−1) */
  // Sliders: relative to the target region (§6). F axis 0 = region start, R axis 0 = region start,
  // R seqLen = region end, R tail = region end + outer flank.
  // F [start,end] -> template 1-based [a+start+1, a+end+1] (a=target.start-1)
  // R [start,end] -> template 1-based [a+start+1, a+end+1] (same origin as F)
  rebuildSliders(tmpl) {
    const T = tmpl.template_len;
    const a = tmpl.target.start - 1;
    const L = tmpl.target.end - tmpl.target.start + 1;   // seq_len = 目标片段长度(规范 §6.1)
    const flank = tmpl.flank_len;
    const fMin = Math.max(-flank, -a);
    const fMax = Math.min(T, T - a - 1);
    // R 轴与 F 轴同原点:0 = 目标片段起点,seqLen(L-1) = 片段终点(*0),
    // 末端 = 片段终点 + 外侧翼(*flank),截断到模板末端。修复:起点恒为 0,
    // 不再随侧翼长度漂移(旧实现以 template_len 为基准,起点值 = T-e+1
    // 随提取范围伸缩,且允许反向范围伸入片段上游侧翼)
    // R axis shares the F origin: 0 = target-region start, seqLen (L-1) = region end (*0),
    // tail = region end + outer flank (*flank), truncated at the template end. Fix: the
    // origin is constant 0 — it no longer drifts with the flank length (the old
    // implementation was template_len-based: the origin value T-e+1 moved with the
    // extraction span and let the reverse range reach into the upstream flank)
    const rMin = 0;
    const rMax = Math.min(L - 1 + flank, T - a - 1);
    document.getElementById("loc-f-min").textContent = String(fMin);
    document.getElementById("loc-r-min").textContent = String(rMin);
    const mk = (id, axis, seqLen) => {
      const cv = document.getElementById(id);
      // 位图宽度 = 容器当前可用宽度(flex:1 占满剩余空间):轨道长度按窗口当前宽度的
      // 百分比设计,窗口越宽轨道越长;显示后才测量(display:none 时宽度为 0,兜底 300)
      // Bitmap width = the container's current available width (flex:1 fills the remaining space):
      // the track length follows a percentage of the current window width; measured after the
      // card is shown (width is 0 while display:none, fall back to 300)
      cv.width = Math.max(120, Math.round(cv.getBoundingClientRect().width || 300));
      const sl = new DualSlider(cv, {
        axis, seqLen,
        minV: axis === "f" ? fMin : rMin,
        maxV: axis === "f" ? fMax : rMax,
        start: axis === "f" ? Math.max(fMin, -flank) : Math.max(rMin, L - 1),
        end: axis === "f" ? 0 : Math.min(rMax, L - 1 + flank),
        onSet: () => {
          // 只同步当前轴:restoreNamed 逐轴 set 时,另一轴若尚未恢复,
          // 全量同步会把其默认值写进 named 并落盘(F/R 参数切页丢失的
          // 根因之一)——拖动时另一轴值本就不变,单轴同步语义等价
          // Sync only the current axis: during restoreNamed's per-axis set,
          // a full sync would write the other (not-yet-restored) axis's
          // default into named and persist it. Dragging never changes the
          // other axis anyway, so per-axis sync is equivalent.
          this.named[axis === "f" ? "f" : "r"] =
            [this.sliders[axis].start, this.sliders[axis].end];
          this.syncRangeInputs();
          this.saveState();
        },
      });
      this.sliders[axis] = sl;
      return sl;
    };
    // 先解绑旧实例事件(否则拖动会同时驱动新旧两个滑块,画布横跳)
    // Unbind old instances first (otherwise dragging drives both sliders)
    if (this.sliders.f) this.sliders.f.unbind();
    if (this.sliders.r) this.sliders.r.unbind();
    this.sliders = {};
    mk("loc-f-slider", "f", T);
    mk("loc-r-slider", "r", L - 1);   // R 轴 seqLen = 目标片段长度(0 基终点,*N 基准)
    // 注意:不再在此处覆盖 this.named.f/r —— rebuildSliders 被 restoreNamed
    // 调用时,调用方先恢复的保存值会被这里的默认值吞掉(F/R 参数切页丢失)。
    // named.f/r 由调用方在 rebuildSliders 之后设置。
    // Deliberately not overwriting this.named.f/r here: when restoreNamed calls
    // rebuildSliders, this would clobber the saved ranges it just restored
    // (the F/R range loss on page switches). Callers set named.f/r after.
  },

  /* 窗口 resize → 滑块重适配:位图宽度 = 容器当前可用宽度(轨道长度 = 窗口宽度百分比),
     值不变只重渲染(刻度/手柄位置随新宽度自动换算) */
  // Window resize -> re-fit sliders: bitmap width = current container width (track = a percentage of the window width);
  // values are unchanged, only re-rendered (ticks/handles are recomputed for the new width)
  fitSliders() {
    if (!this.sliders || !this.sliders.f) return;
    ["f", "r"].forEach((ax) => {
      const sl = this.sliders[ax];
      const w = Math.max(120, Math.round(sl.cv.getBoundingClientRect().width || 300));
      if (w !== sl.cv.width) {
        sl.cv.width = w;
        sl.render();
      }
    });
  },

  /* 参数卡范围输入框 ← 滑块 同步;*N 记法(§6) */
  // Params-card range inputs <-- slider sync; *N notation (§6)
  syncRangeInputs() {
    const tmpl = this.named.tmpl;
    if (!tmpl || !this.sliders.f) return;
    const T = tmpl.template_len;
    const L = tmpl.target.end - tmpl.target.start + 1;   // seq_len = 目标片段长度(R 轴 *N 基准)
    document.getElementById("p-f-start").value = DualSlider.fmt(this.named.f[0], "f", T);
    document.getElementById("p-f-end").value = DualSlider.fmt(this.named.f[1], "f", T);
    document.getElementById("p-r-start").value = DualSlider.fmt(this.named.r[0], "r", L - 1);
    document.getElementById("p-r-end").value = DualSlider.fmt(this.named.r[1], "r", L - 1);
    ["p-f-start", "p-f-end", "p-r-start", "p-r-end"].forEach((id) => {
      const el = document.getElementById(id);
      el.disabled = false;
      el.classList.remove("invalid");
    });
  },

  /* 手动输入:解析 *N → 轴值;非法/相等/越序红框提示,不静默交换 */
  // Manual input: parse *N -> axis value; invalid/equal/out-of-order values get a red-border hint, never silently swapped
  onRangeInput(el) {
    const tmpl = this.named.tmpl;
    if (!tmpl || !this.sliders.f) return;
    const id = el.id;
    const axis = id === "p-f-start" || id === "p-f-end" ? "f" : "r";
    const which = id.endsWith("-start") ? "start" : "end";
    const T = tmpl.template_len;
    const L = tmpl.target.end - tmpl.target.start + 1;   // seq_len = 目标片段长度(R 轴 *N 基准)
    const raw = el.value.trim();
    if (raw === "") { el.classList.remove("invalid"); return; }
    let v = null;
    if (/^\*\d+$/.test(raw)) {
      const n = parseInt(raw.slice(1), 10);
      v = axis === "f" ? -n : (L - 1) + n;
    } else if (/^-?\d+$/.test(raw)) {
      v = parseInt(raw, 10);
    }
    const sl = this.sliders[axis];
    if (v === null || v < sl.minV || v > sl.maxV) {
      el.classList.add("invalid");
      return;
    }
    el.classList.remove("invalid");
    const other = which === "start" ? sl.end : sl.start;
    if (v === other) {
      el.classList.add("invalid");
      toast(t("design.range_equal"));
      return;
    }
    if ((which === "start" && v > other) || (which === "end" && v < other)) {
      el.classList.add("invalid");
      toast(t("design.range_order"));
      return;
    }
    sl.set(which === "start" ? v : sl.start, which === "end" ? v : sl.end, which);
  },

  renderTmplInfo() {
    const el = document.getElementById("loc-tmpl-info");
    const info = this.named.tmpl;
    if (!info) { el.textContent = "—"; return; }
    // range 为条目上的基因组坐标(旧会话缓存可能缺失,兜底用目标区间)
    // range is the genomic coordinates on the entry (may be missing in old session caches; falls back to the target region)
    const rg = info.range && info.range.length === 2 ? `${info.range[0]}-${info.range[1]}`
      : `${info.target.start}-${info.target.end}`;
    const extract = t("design.tmpl_extract", { s: info.extract.start, e: info.extract.end });
    const trunc = info.truncated
      ? " " + t("design.tmpl_trunc", { s: info.requested.start, e: info.requested.end }) : "";
    // 信息卡显示基因/染色体名(entry):name= 是用户显示标签,只用于
    // 结果摘要头与导出(R36),不替换信息卡中的真实条目名
    // The info card shows the gene/chromosome name (entry): name= is a
    // user display label used only in the result summary/exports (R36),
    // not a replacement for the real entry name here
    el.innerHTML =
      `${t("design.tmpl_summary", { entry: `<b>${escapeHtml(info.entry)}</b>`, len: info.entry_len, strand: "+" })}<br>
       ${t("design.tmpl_target", { rg })} · ${extract}${trunc}<br>
       ${t("design.tmpl_stats", { len: info.template_len, flank: info.flank_len, s: info.target.start, e: info.target.end })}`;
  },

  /* ---------------- 运行 ---------------- */
  // ---------------- Run ----------------

  /* 普通模式相对产物基准=模板序列长度;名称型=目标区间长度 */
  // Relative-product baseline: plain mode = template sequence length; named mode = target-region length
  plainSeqLen(text) {
    // `#` 注释行不计入长度(与后端 strip_comment_lines 一致)
    // `#` comment lines are excluded (matches backend strip_comment_lines)
    return (text || "").split(/\r?\n/)
      .filter((l) => !l.startsWith(">") && !l.trim().startsWith("#"))
      .join("").replace(/\s+/g, "").length;
  },

  async run() {
    const db = document.getElementById("dp-db").value;
    const params = this.collectParams();
    const body = { db_prefix: db, params };
    if (this.named.tmpl) {
      // 名称型:单任务,模板=提取序列,产物换算基准=目标区间长度(§6/§7)。
      // 条目库与特异性比对库可能来自名称型参数(database=/targetbase=),
      // 页面未选库也能运行——缺省回落页面库。
      // Named mode: single job; template = extracted sequence; product
      // conversion baseline = target-region length (§6/§7). The entry db and
      // specificity db may come from the named query itself (database=/
      // targetbase=), so the page db is optional here — defaults apply.
      const tmpl = this.named.tmpl;
      if (tmpl.db_prefix) body.db_prefix = tmpl.db_prefix;
      if (tmpl.spec_db) body.spec_db = tmpl.spec_db;
      const T = tmpl.template_len;
      const a = tmpl.target.start - 1;
      const p = { ...params };
      // targetbase=False:该查询不进行特异性比对(不做 blastn-short 逆向验证,
      // 深度统计照常走页面库)——与页面上勾选"跳过特异性查询"同语义,查询级
      // 参数优先(覆盖页面复选框)
      // targetbase=False: this query skips specificity (no blastn-short reverse
      // validation; depth counting still runs on the page db) — same semantics
      // as the "skip specificity" checkbox, with per-query precedence
      if (tmpl.spec_skip) p.skip_spec_eval = true;
      p.left_range = [a + this.named.f[0] + 1, a + this.named.f[1] + 1];
      // R 轴与 F 轴同原点(0 = 目标片段起点):模板 1-based = a + v + 1
      // R axis shares the F origin (0 = target-region start): template 1-based = a + v + 1
      p.right_range = [a + this.named.r[0] + 1, a + this.named.r[1] + 1];
      if (p.product_len_mode === "relative" && p.mode !== "single") {
        const T0 = tmpl.target.end - tmpl.target.start + 1;
        const o1 = Number(document.getElementById("p-product_offset1").value) || 0;
        const o2 = Number(document.getElementById("p-product_offset2").value) || 0;
        p.product_len_min = Math.max(1, T0 + Math.min(o1, o2));
        p.product_len_max = Math.max(p.product_len_min, T0 + Math.max(o1, o2));
        p.product_len_mode = "absolute";
      }
      body.locate_jobs = [{
        template: tmpl.template_fasta,
        target: [tmpl.target.start, tmpl.target.end],
        params: p,
        ctx: {
          entry: tmpl.entry, strand: tmpl.strand,
          display_name: tmpl.display_name || tmpl.entry,
          extract_start: tmpl.extract.start, extract_end: tmpl.extract.end,
          requested: tmpl.requested, flank_len: tmpl.flank_len,
          template_len: T,
        },
      }];
      // k-mer 缓存(sessionStorage 中键匹配本模板/库/模式的既往分析;
      // 检测范围覆盖由后端按缓存 detected_ranges 把关)
      // K-mer caches from sessionStorage whose keys match this template/db/mode
      // (detected-range coverage is enforced by the backend via detected_ranges)
      body.kmer_caches = await this.collectKmerCaches(
        tmpl.template_fasta, body.db_prefix, body.spec_db, p.mode);
    } else {
      if (!db) { toast(t("design.db_empty")); return; }
      // 重名标题加后缀(仅提交/解析时改写,输入框保持原文)
      // Duplicate headers get " (N)" suffixes at submit/parse time only (the
      // input box keeps the original text)
      const template = dedupeFastaHeaders(document.getElementById("dp-template").value).trim();
      if (!template) { toast(t("design.template_empty")); return; }
      const p = { ...params };
      if (p.product_len_mode === "relative" && p.mode !== "single") {
        // 普通模式:基准=模板总长
        // Plain mode: baseline = total template length
        const T0 = this.plainSeqLen(template);
        const o1 = Number(document.getElementById("p-product_offset1").value) || 0;
        const o2 = Number(document.getElementById("p-product_offset2").value) || 0;
        p.product_len_min = Math.max(1, T0 + Math.min(o1, o2));
        p.product_len_max = Math.max(p.product_len_min, T0 + Math.max(o1, o2));
        p.product_len_mode = "absolute";
      }
      body.template = template;
      body.target = null;
      // k-mer 缓存(标准模式无设计范围,全检测)
      // K-mer caches (plain mode has no design ranges: full detection)
      body.kmer_caches = await this.collectKmerCaches(
        template, body.db_prefix, undefined, params.mode);
    }
    const btn = document.getElementById("btn-design");
    btn.disabled = true;
    btn.textContent = t("design.running");
    try {
      const { task_id } = await api("/api/design/run", { method: "POST", json: true, body });
      this.taskId = task_id;
      this.result = null;
      // 运行状态行:清掉上次的终态,显示任务号;spinner/秒表由 watchTask 驱动
      // Run status line: clear the previous terminal state, show the task id;
      // the spinner/stopwatch are driven by watchTask
      document.getElementById("run-status").textContent = "";
      document.getElementById("run-hint").textContent = t("design.run_title") + " … " + task_id;
      clearResult("bp_design_result");   // 新旧结果都清(localStorage + IndexedDB)
      // Clear any old result from both stores (localStorage + IndexedDB)
      document.getElementById("design-result").style.display = "none";
      await this.watchTask(task_id);
      const { status, result } = await api(`/api/tasks/${task_id}/result`);
      if (status === "succeeded" && result) {
        this.result = result;
        this.resultIdx = 0;
        this.selPair = null;
        this.storeKmerCaches(result);   // k-mer 分析入 sessionStorage,同序列再设计即复用
        // K-mer analyses go into sessionStorage; re-designing the same sequence reuses them
        this.renderResult();
        await this.saveResult();   // 结果就绪即存会话(切页/刷新回来仍显示设计结果;await 确保写完再允许切页)
        // Save to session as soon as the result is ready (the design result still shows after page switch or refresh; await so the write completes before any navigation)
      } else if (status === "failed") {
        toast(result?.error || t("design.failed"));
      }
    } catch (err) {
      toast(err.message);
    } finally {
      btn.disabled = false;
      btn.textContent = t("design.start");
    }
  },

  /* 等待任务结束(SSE 已在 App 维护任务状态,轮询终态)。
     运行感知与比对/建库页一致:按钮禁用 + spinner/秒表状态行 + "可继续
     其他页面"提示,大库设计耗时较长,秒表让"仍在运行"可见可感;
     步骤条(updateRunSteps)随任务日志在页面上推进管线阶段。 */
  // Wait for the task to finish (SSE already tracks task state in App; poll
  // for the terminal state). Run feedback mirrors the blast/db pages:
  // disabled button + spinner/elapsed status line + "you may continue
  // elsewhere" hint; long designs need the timer to make "still running"
  // visible and felt; the step bar (updateRunSteps) advances the pipeline
  // stages on the page from the task logs.
  watchTask(taskId) {
    const status = document.getElementById("run-status");
    const hint = document.getElementById("run-hint-extra");
    this._runTaskId = taskId;
    this._runStarted = Date.now();
    this._runElapsed = null;
    const started = Date.now();
    const setStatus = (cls, html) => {
      if (!status) return;
      status.className = "build-status" + (cls ? " " + cls : "");
      status.innerHTML = html;
    };
    const fmt = (s) => (s < 60 ? s + "s" : Math.floor(s / 60) + "m" + (s % 60) + "s");
    if (hint) hint.textContent = t("design.run_elsewhere");
    this.updateRunSteps(taskId);
    return new Promise((resolve) => {
      const tick = setInterval(() => {
        const s = App.tasks.get(taskId);
        if (!s) return;
        if (s.status !== "running" && s.status !== "pending") {
          clearInterval(tick);
          if (hint) hint.textContent = "";
          // 终态:隐藏步骤条,运行状态行接管(✓ 完成 / ✗ 失败)
          // Terminal state: hide the step bar; the run status line takes over
          const steps = document.getElementById("run-steps");
          if (steps) steps.style.display = "none";
          const secs = Math.floor((Date.now() - started) / 1000);
          this._runElapsed = secs;   // 固定耗时:语言切换重绘不再随时间增长
          if (s.status === "succeeded") setStatus("done", `✓ ${t("design.run_done", { s: fmt(secs) })}`);
          else if (s.status === "cancelled") setStatus("err", `✗ ${t("design.run_canceled")}`);
          else setStatus("err", `✗ ${t("design.run_failed")}`);
          resolve(s);
        } else {
          const secs = Math.floor((Date.now() - started) / 1000);
          setStatus("", `<span class="spinner"></span>${t("design.run_running", { s: fmt(secs) })}`);
          this.updateRunSteps(taskId);
        }
      }, 400);
      setTimeout(() => clearInterval(tick), 40 * 60 * 1000);
    });
  },

  /* 设计管线步骤条(任务 #9):从任务日志识别当前阶段并渲染胶囊序列 ——
     已完成打 ✓、当前步 accent 高亮 + spinner、未到达置灰。日志消息与后端
     primer_design.py 的 on_log 文案对应;命中 k-mer 缓存时第一步被跳过,
     缓存日志排在计数之后 → 第一步/计数同时标记完成。 */
  // Design pipeline step bar: the current stage is identified from the task
  // logs and rendered as a capsule sequence — finished steps get a check
  // mark, the active step is accent-highlighted with a spinner, later steps
  // stay greyed out. Log messages correspond to the on_log texts in the
  // backend primer_design.py; on a k-mer cache hit step 1 is skipped, and
  // the cache log follows the counting one, so both read as done.
  RUN_STEPS: [
    { key: "step1",   re: /第一步|设计查询|定位模板设计/ },
    { key: "kmer",    re: /k-mer 计数/ },
    { key: "profile", re: /命中 k-mer 结果缓存|深度分析|自重复 k-mer 剖面|sgRNA 模式|设计范围外 [0-9]+ bp/ },
    { key: "level",   re: /Level \d|分级设计/ },
    { key: "primer3", re: /primer3 候选|引物生成/ },
    { key: "spec",    re: /blastn-short|逆向验证|结合位点|k-mer 预筛/ },
  ],

  updateRunSteps(taskId) {
    const wrap = document.getElementById("run-steps");
    if (!wrap) return;
    const s = App.tasks.get(taskId);
    const logs = (s && s.logs || []).map((l) => l.msg || "");
    let cur = -1;
    for (let i = 0; i < this.RUN_STEPS.length; i++) {
      if (logs.some((m) => this.RUN_STEPS[i].re.test(m))) cur = i;
    }
    if (cur < 0) { wrap.style.display = "none"; return; }
    wrap.style.display = "";
    wrap.innerHTML = this.RUN_STEPS.map((st, i) => {
      const done = i < cur;
      const active = i === cur;
      const mark = done ? "✓ " : (active ? "<span class='spinner'></span>" : "");
      return `<span class="step${done ? " done" : ""}${active ? " active" : ""}">${mark}${t("design.run_step_" + st.key)}</span>`;
    }).join("");
  },

  /* ---------------- k-mer 结果缓存(sessionStorage + 项目文件, R17) ----------------
     同一项目内对同一序列重复设计时,第一步 blastn + k-mer 计数(最耗时环节)
     由缓存直接跳过:结果携带 kmer_cache(键 = 序列/库/模式/尺度/范围哈希,
     与后端 _kmer_cache_key 逐字一致),前端存入 sessionStorage;项目保存时
     kmer_cache 随结果一并写入,加载项目时再回灌 sessionStorage。 */
  // ---------------- K-mer result cache (sessionStorage + project file, R17) ----------------
  // Re-designing the same sequence in a project skips step-1 blastn + counting
  // (the most expensive phase): results carry kmer_cache (key = hash of
  // seq/db/mode/scales/ranges, byte-identical to the backend _kmer_cache_key),
  // stored in sessionStorage; project save writes it with the result, and
  // loading a project pours it back into sessionStorage.

  /* 与 blastprime/primer_design.py _kmer_cache_key 一致的键:
     sha256(seq \0 depth_db \0 mode \0 "(8, 10, 12, 15)")
     F/R 范围不入键:复用前由后端按缓存记录的检测范围检查覆盖。
     mode 归一化 "sgrna" vs "pair"(standard/single 共享,R34)。 */
  // Key identical to blastprime/primer_design.py _kmer_cache_key:
  // sha256(seq \0 depth_db \0 mode \0 "(8, 10, 12, 15)")
  // F/R ranges are not in the key; the backend checks the cached
  // detected_ranges for coverage before reuse. Mode is normalised to
  // "sgrna" vs "pair" (standard and single share caches, R34).
  async kmerCacheKey(seq, db, mode) {
    const m = mode === "sgrna" ? "sgrna" : "pair";
    const raw = `${seq}\0${db}\0${m}\0(8, 10, 12, 15)`;
    const buf = await crypto.subtle.digest("SHA-256", new TextEncoder().encode(raw));
    return [...new Uint8Array(buf)].map((b) => b.toString(16).padStart(2, "0")).join("");
  },

  /* 模板 FASTA → 每条记录的纯序列(去空白、大写,与后端 SeqIO 解析一致) */
  // Template FASTA -> per-record raw sequence (stripped and uppercased,
  // matching the backend SeqIO parse)
  fastaRecords(text) {
    const recs = [];
    let cur = null;
    for (const line of text.split("\n")) {
      const s = line.trim();
      if (!s || s.startsWith("#")) continue;          // 空行 / # 注释(后端同样忽略)
      if (s.startsWith(">")) { if (cur !== null) recs.push(cur); cur = ""; continue; }
      if (cur === null) cur = "";                     // 无 > 头的纯序列 → 单条
      cur += s.toUpperCase();
    }
    if (cur !== null) recs.push(cur);                 // 收尾:字符串不可变,必须累积完再入列
    return recs.filter((r) => r.length >= 40);        // 短于后端下限的无需缓存
  },

  /* 收集与当前模板匹配的 k-mer 缓存(sessionStorage,带 bp_kmer_cache: 前缀)。
     纯序列入键(去头),返回 [cacheObj...] 供请求体 kmer_caches。 */
  // Collect k-mer caches matching the current template(s) from sessionStorage
  // (prefix bp_kmer_cache:); keys are computed on the raw sequence (headers
  // stripped). Returns [cacheObj...] for the request's kmer_caches.
  async collectKmerCaches(template, db, specDb, mode) {
    try {
      const seqs = this.fastaRecords(template);
      const keys = [];
      for (const seq of seqs) {
        keys.push(await this.kmerCacheKey(seq, specDb || db, mode));
      }
      const keySet = new Set(keys);
      const dbRef = specDb || db;
      const out = [];
      // 遍历 sessionStorage 全部缓存:精确键 或 子串命中(flank 变化导致
      // 模板嵌套,当前模板是缓存模板的子串——后端按 offset 映射复用)。
      // 旧逻辑只按精确键收集,子串候选(不同键的大模板缓存)从未回送,
      // 后端 R29 的子串命中形同虚设。
      // Iterate all sessionStorage caches: exact-key match, or substring
      // match (flank changes nest the template — the current template is a
      // substring of the cached one, reused with an offset by the backend).
      // The old key-only collection never sent substring candidates, making
      // the backend R29 substring reuse dead code.
      for (let i = 0; i < sessionStorage.length; i++) {
        const k = sessionStorage.key(i);
        if (!k || !k.startsWith("bp_kmer_cache:")) continue;
        let c = null;
        try { c = JSON.parse(sessionStorage.getItem(k)); }
        catch (e) { sessionStorage.removeItem(k); continue; }
        if (!c) continue;
        const sub = typeof c.seq === "string" && c.db === dbRef &&
          seqs.some((s) => c.seq.includes(s));
        if (!keySet.has(k.slice(14)) && !sub) continue;
        // 回送前剔除旧结构(R17)里的 step1_hsps 全文:大命中时可达数十 MB,
        // 塞进请求体拖垮整个设计流程;后端只读 key/template_len/depth/
        // profiles/target_loci,不会用到它
        // Strip the old (R17) step1_hsps from the cache before sending:
        // on heavy hits it reached tens of MB and bloated the request;
        // the backend only reads key/template_len/depth/profiles/target_loci
        if ("step1_hsps" in c) { const { step1_hsps, ...rest } = c; out.push(rest); }
        else out.push(c);
      }
      return out;
    } catch (e) { return []; }   // crypto.subtle 不可用等 → 无缓存,正常重算
    // crypto.subtle unavailable etc. -> no cache, recompute normally
  },

  /* 设计结果 → 逐查询写入 sessionStorage(配额溢出静默放弃,旧条目让位) */
  // Design result -> per-query sessionStorage writes (quota overflow is
  // silently dropped, making room for newer entries)
  storeKmerCaches(result) {
    if (!result || !result.results) return;
    for (const q of result.results) {
      const c = q && q.kmer_cache;
      if (!c || !c.key) continue;
      try { sessionStorage.setItem("bp_kmer_cache:" + c.key, JSON.stringify(c)); }
      catch (e) { /* quota: 放弃该条缓存,重算兜底 */ }
    }
    this.updateKmerCacheCount();
  },

  /* 高级参数区:k-mer 缓存条数 + 清空按钮 */
  // Advanced params: k-mer cache count + clear button
  updateKmerCacheCount() {
    const el = document.getElementById("kmer-cache-count");
    if (!el) return;
    let n = 0;
    for (let i = 0; i < sessionStorage.length; i++) {
      if ((sessionStorage.key(i) || "").startsWith("bp_kmer_cache:")) n++;
    }
    el.textContent = t("design.kmer_cache_count", { n });
  },
  clearKmerCaches() {
    let n = 0;
    for (let i = sessionStorage.length - 1; i >= 0; i--) {
      const k = sessionStorage.key(i);
      if (k && k.startsWith("bp_kmer_cache:")) { sessionStorage.removeItem(k); n++; }
    }
    this.updateKmerCacheCount();
    toast(n ? t("design.kmer_cache_cleared", { n }) : t("design.kmer_cache_empty"));
  },

  /* ---------------- 结果会话保存(切页/刷新回来恢复设计结果,仿 blast.js) ---------------- */
  // ---------------- Result session save (restores the design result after page switch/refresh, modeled after blast.js) ----------------

  async saveResult() {
    if (!this.result) return;
    // 瘦身结果(剔除渲染不读的 step1_hsps):大命中时该字段可达数十 MB,
    // 超 5 MB 配额静默失败 → 切页/刷新回来结果丢失(见 app.js slimDesignResult);
    // persistResult 在 localStorage 超限时回退 IndexedDB
    // Slim the result (strip step1_hsps, which rendering never reads): on
    // heavy hits it reaches tens of MB, silently blowing the 5 MB quota and
    // losing the result on page switch/refresh (see app.js slimDesignResult);
    // persistResult falls back to IndexedDB when localStorage overflows
    await persistResult("bp_design_result", {
      _v: 1, taskId: this.taskId, result: slimDesignResult(this.result),
      resultIdx: this.resultIdx, selPair: this.selPair,
    });
  },

  async restoreResult() {
    const data = await loadResult("bp_design_result");
    if (!data || data._v !== 1 || !data.result) {
      await clearResult("bp_design_result");
      return;
    }
    this.result = data.result;
    this.taskId = data.taskId || null;
    this.resultIdx = data.resultIdx || 0;
    this.selPair = data.selPair || null;
    this.storeKmerCaches(data.result);   // 刷新后回灌 k-mer 缓存,再设计即复用
    // Pour k-mer caches back after a refresh so re-designing reuses them
    this.renderResult();
  },

  /* ---------------- 结果渲染 ---------------- */
  // ---------------- Result rendering ----------------

  renderResult() {
    document.getElementById("design-result").style.display = "";
    const r = this.result;
    const locateMode = !!r.locate_mode;
    // 多查询:结果选择器(名称型单任务也经 locate_jobs,保留)
    // Multi-query: result selector (the named-mode single job also goes through locate_jobs; kept)
    let selWrap = document.getElementById("result-query-sel");
    if (locateMode) {
      if (!selWrap) {
        const toolbar = document.querySelector("#design-result .toolbar");
        const div = document.createElement("span");
        div.id = "result-query-sel";
        div.className = "row";
        div.style.cssText = "gap:6px; margin-left:8px";
        toolbar.insertBefore(div, toolbar.querySelector(".spacer"));
        selWrap = div;  // 回填引用,否则下方 selWrap.innerHTML 对 null 报错
        // Backfill the reference, otherwise selWrap.innerHTML below would throw on null
      }
      const sel = document.createElement("select");
      sel.id = "result-query";
      sel.innerHTML = (r.results || []).map((res, i) => {
        const loc = res.locate || {};
        const ok = res.success;
        return `<option value="${i}">${escapeHtml(loc.display_name || loc.entry || res.query || i + 1)} ${ok ? "✓" : "✗"}</option>`;
      }).join("");
      sel.value = this.resultIdx;
      sel.addEventListener("change", () => {
        this.resultIdx = +sel.value;
        this.selPair = null;
        this.renderResult();
      });
      selWrap.innerHTML = "";
      selWrap.appendChild(sel);
    } else if (selWrap) {
      selWrap.innerHTML = "";
    }
    const res0 = (r.results || [])[this.resultIdx];
    if (!res0) return;
    const mode = r.mode || "standard";
    // 概览
    // Overview
    const ov = document.getElementById("design-overview");
    const stats = res0.depth_stats || {};
    const fail = res0.failure;
    // 新引擎:profile_stats 逐碱基特异性剖面统计(global 0~1,3' 0~1);
    // 旧字段 depth_stats 保留作降级
    // New engine: per-base specificity profile stats; depth_stats kept for fallback
    const ps = res0.profile_stats || {};
    const gprof = ps.global || [];
    const tprof = ps.three_prime || [];
    // 范围外位置(R15)未评分(值为 0),概览统计只算可设计区
    // Out-of-range positions (R15) are unscored (0); overview stats cover
    // only the designable region
    const rex = res0.range_excluded || null;
    const statOf = (arr) => {
      if (!arr.length) return null;
      let mn = 1, mx = 0, sum = 0, cnt = 0;
      for (let i = 0; i < arr.length; i++) {
        if (rex && rex.some((rng) => i >= rng[0] - 1 && i < rng[0] - 1 + rng[1])) continue;
        const v = arr[i]; cnt++;
        if (v < mn) mn = v; if (v > mx) mx = v; sum += v;
      }
      if (!cnt) return null;
      return { mn, mx, mean: sum / cnt };
    };
    const gs = statOf(gprof), ts = statOf(tprof);
    const specTags = [];
    if (res0.success) {
      const levels = {};
      (res0.pairs || []).forEach((p) => { levels[p.specificity?.label] = (levels[p.specificity?.label] || 0) + 1; });
      specTags.push(`<span class="tag ok">${t("design.stage_reached")}: ${res0.stage_reached}</span>`);
      specTags.push(`<span class="tag info">${t("design.pairs_count")}: ${(res0.pairs || []).length}</span>`);
      Object.entries(levels).forEach(([k, v]) => specTags.push(`<span class="tag">${escapeHtml(this.specLabel(k))} × ${v}</span>`));
    } else {
      specTags.push(`<span class="tag err">${t("design.failed")}</span>`);
    }
    let title = escapeHtml(res0.query || "query");
    if (res0.locate) {
      const L = res0.locate;
      title = `${escapeHtml(L.display_name || L.entry)} ${t("loc.summary")}` +
        ` (${t("loc.genomic")}: ${L.target_genomic[0]}-${L.target_genomic[1]}, ` +
        `${t("loc.strand")}: ${L.strand === "minus" ? "−" : "+"})`;
    }
    // 四阶段摘要:每级候选 → 通过数
    // Four-stage summary: candidates -> accepted per level
    const lvTags = (res0.levels || []).map((lv) =>
      `<span class="tag" title="${escapeHtml(this.levelName(lv.name))}">L${lv.level}: ${lv.candidates ?? "-"}→${lv.after_pair_check ?? "-"}</span>`).join(" ");
    ov.innerHTML = `
      <div class="row">
        <b class="mono">${title}</b>
        <span class="muted">len ${res0.template_len} bp</span>
        ${specTags.join(" ")}
      </div>
      <div class="hint" style="margin-top:6px">
        ${gs ? t("design.profile_stats", {
          mn: (gs.mn * 100).toFixed(0) + "%",
          mx: (gs.mx * 100).toFixed(0) + "%",
          mean: (gs.mean * 100).toFixed(1) + "%",
          t3: ts ? (ts.mean * 100).toFixed(1) + "%" : "-",
        }) : t("design.overview_stats", {
          m: stats.matched !== undefined ? (stats.matched * 100).toFixed(1) + "%" : "-",
          r: stats.repeat_frac !== undefined ? (stats.repeat_frac * 100).toFixed(1) + "%" : "-",
          d: stats.histogram ? Math.max(...Object.keys(stats.histogram).map(Number)) : "-",
        })}
        ${mode === "sgrna" ? ` · ${t("design.mode_sgrna")}` : mode === "single" ? ` · ${t("design.mode_single")}` : ""}
      </div>
      ${lvTags ? `<div class="row" style="margin-top:4px;gap:4px">${lvTags}</div>` : ""}
      ${fail ? this.renderFailure(fail) : res0.error ? `<div class="tag err" style="margin-top:8px">${escapeHtml(res0.error)}</div>` : ""}
      <div class="row" style="margin-top:8px">
        ${(res0.pairs || []).length ? `
          <span class="muted">${t("design.spec_levels")}</span>` : ""}
      </div>`;
    this.renderDepth();
    this.renderPairs(mode);
  },

  /* 特异性等级标签本地化:后端结果内嵌中文标签(PAIR_LABELS/SPEC_LABELS),
     语言切换时按当前语言映射;未知值原样返回 */
  // Localize backend-embedded Chinese spec-level labels (PAIR_LABELS / SPEC_LABELS)
  specLabel(l) {
    if (!l) return "";
    const M = {
      "基因组唯一": t("design.spec_unique"),
      "基因组唯一匹配": t("design.spec_unique"),
      "3' 端错配豁免": t("design.spec_3p_exempt"),
      "无脱靶产物": t("design.spec_no_product"),
      "可扩增脱靶产物(淘汰)": t("design.spec_offtarget"),
      "不可成对扩增": t("design.spec_unpaired"),
      "已淘汰": t("design.spec_eliminated"),
      "基于 k-mer 深度(未做 blastn-short 逆向验证)": t("design.spec_kmer_only"),
      "命中数超上限(无法完整评估,淘汰)": t("design.spec_truncated"),
    };
    return M[l] || l;
  },

  /* 失败诊断文案本地化:后端带 reason_key 时按当前语言渲染并回填参数;
     旧结果/未知键 → 原样显示后端中文 reason(降级) */
  // Localize the failure-diagnosis text: render via reason_key when present;
  // legacy results / unknown keys fall back to the backend Chinese `reason`
  failureText(f) {
    if (!f.reason_key) return escapeHtml(f.reason || "");
    const key = "design.fail_reason_" + f.reason_key;
    const s = t(key, f.reason_params || {});
    return s === key ? escapeHtml(f.reason || "") : s;
  },

  /* 覆盖失败备注本地化:coverage_note_key 存在时按当前语言渲染(带参数),
     旧结果降级为后端中文原文 */
  // Localize the coverage note: render via coverage_note_key when present,
  // legacy results fall back to the backend Chinese text
  coverageNote(p) {
    const raw = p.coverage_note || "";
    if (!p.coverage_note_key) return escapeHtml(raw);
    const key = "design.coverage_note";
    const s = t(key, p.coverage_note_params || {});
    return s === key ? escapeHtml(raw) : s;
  },

  /* 特异性 note 本地化:specificity.note 内嵌中文,携带 note_key/note_params
     时按当前语言渲染(带参数);旧数据/未知键降级为中文原文 */
  // Localize the specificity note (backend embeds Chinese): render via
  // note_key/note_params when present; legacy data falls back to the raw text
  noteText(sp) {
    const raw = sp.note || "";
    if (!sp.note_key) return escapeHtml(raw);
    const key = "design.note_" + sp.note_key;
    const s = t(key, sp.note_params || {});
    return s === key ? escapeHtml(raw) : s;
  },

  /* 四级名称本地化:后端 LEVEL_NAMES 内嵌中文,按精确串映射当前语言 */
  // Localize the four level names (backend LEVEL_NAMES are Chinese literals)
  levelName(name) {
    const M = {
      "Level 1(count=1,唯一):global≥1.0,3'≥1.0": t("design.level1_name"),
      "Level 2(count 2-3):global≥0.6934,3'≥0.6934": t("design.level2_name"),
      "Level 3(count 4-6):global≥0.5503,3'≥0.5503": t("design.level3_name"),
      "Level 4(count≥7 或无预屏蔽)": t("design.level4_name"),
    };
    return M[name] || name;
  },

  renderFailure(f) {
    // 新引擎失败诊断:失败阶段 + 候选数漏斗 + 低特异区 + 建议
    // New engine failure diagnosis: failing stage + candidate funnel + low-specificity regions + suggestion
    const stgName = (s) => {
      const map = {
        "PRIMER3_NO_CANDIDATE": t("design.fail_stage_primer3"),
        "TARGET_NO_HIGH_SPECIFICITY_REGION": t("design.fail_stage_nospec"),
        "PRIMER_PAIR_OFFTARGET": t("design.fail_stage_offtarget"),
        "NO_ACCEPTABLE_PAIR": t("design.fail_stage_nopair"),
        "DB_INDEX_FAILED": t("design.fail_stage_dbindex"),
      };
      return map[s] || s;
    };
    const funnel = (f.candidate_count != null) ? `
      <div class="hint" style="margin-top:6px">${t("design.fail_funnel", {
        cand: f.candidate_count,
        pre: f.candidate_after_prefilter,
        local: f.candidate_after_local_check,
        pair: f.candidate_after_pair_check,
      })}</div>` : "";
    // 低特异性区(新)或重复区(旧降级)
    // Low-specificity regions (new) or repeat regions (legacy fallback)
    const lows = (f.top_low_specificity_regions || []).map((x) =>
      `<span class="tag err">${x.start}-${x.end} (${(x.min_score * 100).toFixed(0)}%)</span>`).join(" ");
    const regs = lows || (f.repeat_stats?.top_repeat_regions || []).map((x) =>
      `<span class="tag err">${x.start}-${x.end} (depth ${x.max_depth})</span>`).join(" ");
    return `
      <div class="row" style="margin-top:8px;gap:6px">
        <span class="tag err">${escapeHtml(stgName(f.failure_stage))}</span>
        <span class="muted">${this.failureText(f)}</span>
      </div>
      ${funnel}
      ${regs ? `<div class="hint" style="margin-top:6px">${t("design.fail_repeat", { regs })}</div>` : ""}
      ${f.suggestion_key
        ? `<div class="hint" style="margin-top:4px">${t("design.fail_suggest_" + f.suggestion_key)}</div>`
        : f.suggestion ? `<div class="hint" style="margin-top:4px">${t("design.fail_suggest", { sug: escapeHtml(f.suggestion) })}</div>` : ""}`;
  },

  /* ---------------- 特异性剖面图 ---------------- */
  // ---------------- Specificity profile plot ----------------

  renderDepth() {
    const cv = document.getElementById("depth-canvas");
    const r = this.result;
    const res0 = (r.results || [])[this.resultIdx];
    if (!res0) return;
    const n = res0.template_len;
    const depth = res0.depth || [];
    // 新引擎:逐碱基特异性分数(0~1);旧引擎:匹配深度 0/1/≥2(降级)
    // New engine: per-base specificity score (0~1); legacy depth (fallback)
    const ps = res0.profile_stats || {};
    const gprof = ps.global || null;
    const tprof = ps.three_prime || null;
    const hasProfile = !!gprof && gprof.length >= n;
    const dpr = window.devicePixelRatio || 1;
    const W = cv.parentElement.clientWidth - 4;
    const BAND = 26, DEPTH = 96, PAD = 12, TOP = 6;
    const H = TOP + BAND + DEPTH + 26;
    cv.width = W * dpr; cv.height = H * dpr;
    cv.style.height = H + "px";
    const ctx = cv.getContext("2d");
    ctx.scale(dpr, dpr);
    ctx.clearRect(0, 0, W, H);
    this.canvasRects = [];
    const colors = getComputedStyle(document.documentElement);
    const xOf = (pos) => PAD + (pos - 1) * ((W - 2 * PAD) / Math.max(1, n));
    const cvColors = {
      avail: colors.getPropertyValue("--ok"),
      buffer: colors.getPropertyValue("--warn"),
      // L3(count 4-6)中间档:主题无关的固定橙(与 --warn 黄/--err 红区分)
      // L3 (count 4-6) mid tier: fixed orange, distinct from warn-yellow and err-red
      l3: "#f59e0b",
      core: colors.getPropertyValue("--fg-faint"),
      err: colors.getPropertyValue("--err") || "#dc2626",
      fg: colors.getPropertyValue("--fg"),
    };
    // 取设计时的屏蔽掩码:最后带 mask 的 stage(旧引擎降级用)
    // Take the masking mask from design time: the last stage that carries a mask (legacy fallback)
    let mask = null;
    const stages = res0.stages || [];
    for (let i = stages.length - 1; i >= 0; i--) {
      if (stages[i].mask) { mask = stages[i].mask; break; }
    }
    // 设计范围外区域(R13/R15):正、反向引物 3' 端都无法放置的区间补集,
    // 后端 _range_excluded 下发为 (start, len) 1-based 列表(仅定位模式
    // 存在)。R15 起该区域不进行 k-mer 评分(剖面为 0),不再有可展示的
    // 重复语义 → 灰标无条件覆盖。
    // R13/R15: F/R design-range complement (no primer of either strand can
    // land here) from the backend's _range_excluded as 1-based (start, len)
    // pairs (locate mode only). Since R15 these positions are unscored
    // (profile 0) and carry no repeat semantics, so the gray mark wins.
    const rex = res0.range_excluded || null;
    const con = rex ? new Array(n).fill(false) : null;
    if (con) {
      for (const rng of rex) {
        const s = rng[0] - 1, e = Math.min(n, s + rng[1]);
        for (let p = s; p < e; p++) con[p] = true;
      }
    }
    // 色带:global 特异性分按 count 分档阈值分级(R33/R43/R48)——L1 固定
    // ≥1.0(count=1)绿;L2/L3 的阈值优先从结果 level_thresholds 读取
    // (后端下发的完整四级配置,含未实际跑到的级别——Level 1 即成功时
    // levels 不含 L2/L3,仅读 levels 会回退默认),再回退 levels,最后
    // 回退 count 分档默认;旧引擎按深度/掩码
    // Color band: global score tiered by the count-band thresholds
    // (R33/R43/R48) — L1 fixed at ≥1.0 (count=1) green; L2/L3 thresholds
    // come from the result's level_thresholds (the full four-level config
    // the backend sends, including levels never reached — reading only
    // levels would fall back to defaults when Level 1 succeeds), then
    // levels, then the count-band defaults.
    const TH_DEFAULT = [1.0, 0.6934, 0.5503];
    const thSrc = Array.isArray(res0.level_thresholds) ? res0.level_thresholds
      : res0.levels || [];
    const pick = (idx, fallback) => {
      const th = thSrc[idx];
      return (th && Array.isArray(th) && th[0] != null) ? Number(th[0]) : fallback;
    };
    const L1_TH = TH_DEFAULT[0];
    let L2_TH = pick(1, TH_DEFAULT[1]);
    let L3_TH = pick(2, TH_DEFAULT[2]);
    for (let i = 0; i < n; i++) {
      let col;
      if (hasProfile) {
        const s = gprof[i] ?? 0;
        col = s >= L1_TH ? cvColors.avail
          : s >= L2_TH ? cvColors.buffer
          : s >= L3_TH ? cvColors.l3
          : cvColors.err;
      } else if (mask) {
        col = mask[i] === "available" ? cvColors.avail
          : mask[i] === "buffer" ? cvColors.buffer : cvColors.core;
      } else {
        col = (depth[i] || 0) < 1 ? cvColors.avail : cvColors.core;
      }
      // 范围外位置剖面恒为 0(未评分)→ 会落到 err 红;灰标无条件覆盖
      // Out-of-range positions have profile 0 (unscored) -> err red; the
      // gray mark overrides unconditionally
      if (con && con[i]) col = cvColors.core;
      ctx.fillStyle = col;
      const x = xOf(i + 1);
      ctx.fillRect(x, TOP, Math.max(1, (W - 2 * PAD) / n + 0.4), BAND);
    }
    // 目标区域
    // Target region
    const tg = res0.target || {};
    if (tg.start && tg.end && tg.end > tg.start) {
      ctx.fillStyle = "rgba(37,99,235,.22)";
      ctx.fillRect(xOf(tg.start), TOP, Math.max(2, xOf(tg.end + 1) - xOf(tg.start)), BAND);
      // 目标区蓝框加粗:用户反馈 1px 不够明显
      // Thicker target-region outline: the 1px frame was too subtle
      ctx.strokeStyle = colors.getPropertyValue("--accent");
      ctx.lineWidth = 2.5;
      ctx.strokeRect(xOf(tg.start) + .5, TOP + .5, Math.max(2, xOf(tg.end + 1) - xOf(tg.start)) - 1, BAND - 1);
    }
    if (hasProfile) {
      // 曲线区:global 实线 + 3' 端虚线,分数 0~1
      // Curve area: global solid + 3'-end dashed, score 0~1
      const yOf = (s) => TOP + BAND + DEPTH - Math.max(0, Math.min(1, s)) * (DEPTH - 8);
      const line = (arr, style, dash) => {
        ctx.strokeStyle = style;
        ctx.lineWidth = dash ? 1 : 1.4;
        ctx.setLineDash(dash || []);
        ctx.beginPath();
        for (let i = 0; i < n; i++) {
          const x = xOf(i + 1);
          const y = yOf(arr[i] ?? 0);
          if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        }
        ctx.stroke();
        ctx.setLineDash([]);
      };
      // 3' 端剖面用 accent 色:深色主题下 --border 与背景同色阶,虚线会不可见
      // 3'-end profile in accent: --border is nearly invisible on dark themes
      line(tprof || [], colors.getPropertyValue("--accent"), [3, 3]);
      line(gprof, cvColors.fg, null);
      // 参考线:count 分档阈值 1.0 / L2 / L3(R33/R43,L2/L3 随参数)
      // Reference lines at the count-band thresholds (R33/R43; L2/L3
      // follow the parameterised stage thresholds)
      ctx.strokeStyle = colors.getPropertyValue("--border");
      ctx.setLineDash([2, 4]);
      [L1_TH, L2_TH, L3_TH].forEach((th) => {
        const y = yOf(th);
        ctx.beginPath(); ctx.moveTo(PAD, y); ctx.lineTo(W - PAD, y); ctx.stroke();
        ctx.fillStyle = colors.getPropertyValue("--fg-faint");
        ctx.font = "9px sans-serif";
        ctx.fillText(String(th), 2, y + 3);
      });
      ctx.setLineDash([]);
      ctx.fillStyle = colors.getPropertyValue("--fg-soft");
      ctx.font = "11px sans-serif";
      ctx.fillText(`${t("design.depth_legend")}  ·  ${n} bp`, PAD, H - 6);
    } else {
      // 旧引擎降级:深度台阶 0/1/≥2
      // Legacy fallback: depth steps 0/1/>=2
      const maxD = Math.max(2, ...(depth || []).slice(0, 5000));
      const yOf = (d) => TOP + BAND + DEPTH - Math.min(1, d / maxD) * (DEPTH - 8);
      ctx.strokeStyle = colors.getPropertyValue("--fg-soft");
      ctx.lineWidth = 1;
      ctx.beginPath();
      for (let i = 0; i < n; i++) {
        const x = xOf(i + 1);
        const y = yOf(depth[i] || 0);
        if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.stroke();
      ctx.strokeStyle = colors.getPropertyValue("--border");
      ctx.setLineDash([3, 3]);
      [0, 1, 2].forEach((d) => {
        if (d <= maxD) {
          const y = yOf(d);
          ctx.beginPath(); ctx.moveTo(PAD, y); ctx.lineTo(W - PAD, y); ctx.stroke();
          ctx.setLineDash([]);
          ctx.fillStyle = colors.getPropertyValue("--fg-faint");
          ctx.font = "10px sans-serif";
          ctx.fillText(String(d), 2, y + 3);
          ctx.setLineDash([3, 3]);
        }
      });
      ctx.setLineDash([]);
      ctx.fillStyle = colors.getPropertyValue("--fg-soft");
      ctx.font = "11px sans-serif";
      ctx.fillText(`depth (0/1/≥2)  ·  ${n} bp`, PAD, H - 6);
    }
    // 引物对标记(相邻标签交替两行,避免短模板 + 多引物对时文字重叠)
    // Primer-pair markers (adjacent labels alternate between two rows to
    // avoid text overlap on short templates with many primer pairs)
    const pairs = res0.pairs || [];
    let lastX = -1e9, lastRow = 0;
    const mk = (x, col, label) => {
      ctx.fillStyle = col;
      ctx.beginPath();
      ctx.moveTo(x, TOP + BAND + 2);
      ctx.lineTo(x + 7, TOP + BAND + 2);
      ctx.lineTo(x + 3.5, TOP + BAND + 10);
      ctx.closePath();
      ctx.fill();
      let row = 0;
      if (Math.abs(x - lastX) < 13) row = lastRow ^ 1;   // 相邻 → 换行
      lastX = x; lastRow = row;
      ctx.fillStyle = colors.getPropertyValue("--fg-soft");
      ctx.font = "9px sans-serif";
      ctx.fillText(label, x - 2, TOP + BAND + 20 + row * 10);
      ctx.font = "11px sans-serif";
    };
    pairs.forEach((p, i) => {
      const l = p.left || {};
      const rp = p.right || {};
      if (l.start) mk(xOf(l.start), cvColors.err, `F${i + 1}`);
      if (rp.start) mk(xOf(rp.start + rp.len - 1), "#c026d3", `R${i + 1}`);
      this.canvasRects.push({ pairIdx: i, x: xOf(l.start || 1) });
    });
    // 图例随参数化阈值动态更新(R48):L2/L3 区间文字填入实际阈值
    // Legend follows the parameterised thresholds (R48): the L2/L3 band
    // labels show the actual threshold values
    const l2Text = `${L2_TH.toFixed(4)}~1.0`;
    const l3Text = `${L3_TH.toFixed(4)}~${L2_TH.toFixed(4)}`;
    const l2El = document.querySelector('[data-i18n="design.legend_buffer"]');
    const l3El = document.querySelector('[data-i18n="design.legend_l3"]');
    const lmaskEl = document.querySelector('[data-i18n="design.legend_masked"]');
    if (l2El) l2El.textContent = t("design.legend_buffer", { th: l2Text });
    if (l3El) l3El.textContent = t("design.legend_l3", { th: l3Text });
    if (lmaskEl) lmaskEl.textContent = t("design.legend_masked", { th: L3_TH.toFixed(4) });
    // 点击画布 → 定位表格行
    // Clicking the canvas -> locate the table row
    cv.onclick = (e) => {
      const rect = cv.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const hit = this.canvasRects.find((r) => Math.abs(mx - r.x) < 8);
      if (!hit) return;
      const row = document.querySelector(`#pr-table-wrap tr[data-idx="${hit.pairIdx}"]`);
      if (row) { row.scrollIntoView({ block: "center" }); row.click(); }
    };
  },

  /* ---------------- 引物表 ---------------- */
  // ---------------- Primer table ----------------

  renderPairs(mode) {
    const r = this.result;
    const res0 = (r.results || [])[this.resultIdx];
    const wrap = document.getElementById("pr-table-wrap");
    const pairs = res0.pairs || [];
    if (!pairs.length) {
      wrap.innerHTML = `<div class="empty">${res0.failure ? this.failureText(res0.failure) : "—"}</div>`;
      return;
    }
    const locate = !!res0.locate;
    const specTag = (p) => {
      const sp = p.specificity || {};
      const lv = sp.spec_score >= 100 ? "ok" : sp.spec_score >= 60 ? "warn" : "err";
      return `<span class="tag ${lv}" title="${escapeHtml(this.specLabel(sp.label))}">${sp.spec_score ?? "-"}</span>`;
    };
    const covTag = (p) => {
      if (p.covers_target === undefined || p.covers_target === null) return "";
      return p.covers_target
        ? `<span class="tag ok" title="${escapeHtml(p.coverage_note || "")}">${t("loc.covers")}</span>`
        : `<span class="tag warn" title="${escapeHtml(p.coverage_note || "")}">${t("loc.not_covers")}</span>`;
    };
    const absTxt = (pr) => (pr.abs_start ? ` <span class="muted">abs ${pr.abs_start}-${pr.abs_end}</span>` : "");
    const rows = pairs.map((p, i) => {
      const l = p.left || {}, rp = p.right || {};
      if (mode === "sgrna") {
        return `<tr class="pr-row" data-idx="${i}">
          <td class="num">${i + 1}</td>
          <td class="mono" style="font-size:12px">${escapeHtml(l.seq || "")}</td>
          <td>${escapeHtml(p.pam || "")} ${escapeHtml(p.strand || "")}</td>
          <td class="num">${l.gc ?? "-"}%</td>
          <td class="num">${p.specificity?.off_target_sites ?? 0}</td>
          <td>${specTag(p)}</td>
          <td class="num"><b>${p.composite_score ?? "-"}</b></td>
        </tr>`;
      }
      if (mode === "single") {
        // R35:单引物模式含正向(F)与反向(R)两套独立引物,行首标记方向
        // R35: single-primer mode lists forward (F) and reverse (R) primers
        const sideTag = p.side === "reverse" ? "<span class='tag'>R</span> " : "<span class='tag'>F</span> ";
        return `<tr class="pr-row" data-idx="${i}">
          <td class="num">${i + 1}</td>
          <td class="mono" style="font-size:12px">${sideTag}${escapeHtml(l.seq || "")}${absTxt(l)}</td>
          <td class="num">${l.tm ?? "-"} / ${l.gc ?? "-"}%</td>
          <td class="num">${l.len ?? "-"}</td>
          <td class="num">${p.specificity?.off_target_sites ?? 0}</td>
          <td>${specTag(p)}</td>
          <td class="num"><b>${p.composite_score ?? "-"}</b></td>
        </tr>`;
      }
      return `<tr class="pr-row" data-idx="${i}">
        <td class="num">${i + 1}</td>
        <td class="mono" style="font-size:12px">${escapeHtml(l.seq || "")}<br><span class="muted">${l.start}-${(l.start || 0) + (l.len || 0) - 1}</span>${absTxt(l)}</td>
        <td class="num">${l.tm ?? "-"} / ${l.gc ?? "-"}%</td>
        <td class="mono" style="font-size:12px">${escapeHtml(rp.seq || "")}<br><span class="muted">${rp.start}-${(rp.start || 0) + (rp.len || 0) - 1}</span>${absTxt(rp)}</td>
        <td class="num">${l.tm && rp.tm ? Math.abs(l.tm - rp.tm).toFixed(1) : "-"}</td>
        <td class="num">${p.product_len ?? "-"}</td>
        <td class="num">${p.dimer?.max_consec ?? "-"}</td>
        <td class="num">${p.specificity?.off_target_sites ?? 0}</td>
        <td>${specTag(p)}</td>
        <td class="num"><b>${p.composite_score ?? "-"}</b></td>
        ${locate ? `<td>${covTag(p)}</td>` : ""}
      </tr>`;
    });
    const head = mode === "sgrna"
      ? `<tr><th>#</th><th>Guide (20 bp)</th><th>PAM / strand</th><th class="num">GC</th><th class="num">${t("design.col_offtarget")}</th><th>${t("design.col_spec")}</th><th class="num">${t("design.col_score")}</th></tr>`
      : mode === "single"
        ? `<tr><th>#</th><th>Primer</th><th class="num">Tm / GC</th><th class="num">Len</th><th class="num">${t("design.col_offtarget")}</th><th>${t("design.col_spec")}</th><th class="num">${t("design.col_score")}</th></tr>`
        : `<tr><th>#</th><th>${t("design.col_fwd")}</th><th class="num">Tm / GC</th><th>${t("design.col_rev")}</th><th class="num">ΔTm</th><th class="num">${t("design.col_product")}</th><th class="num">${t("design.col_dimer")}</th><th class="num">${t("design.col_offtarget")}</th><th>${t("design.col_spec")}</th><th class="num">${t("design.col_score")}</th>${locate ? `<th>${t("loc.covers")}</th>` : ""}</tr>`;
    wrap.innerHTML = `<table class="data"><thead>${head}</thead><tbody>${rows.join("")}</tbody></table>`;
    wrap.querySelectorAll(".pr-row").forEach((tr) => {
      // 语言切换/重绘后恢复选中高亮
      // restore the selected-row highlight after a language-switch re-render
      if (this.selPair != null && +tr.dataset.idx === this.selPair) tr.classList.add("sel");
      tr.addEventListener("click", () => {
        wrap.querySelectorAll(".pr-row.sel").forEach((x) => x.classList.remove("sel"));
        tr.classList.add("sel");
        this.selPair = +tr.dataset.idx;  // 记录选中行,语言切换时重绘详情
        // remember the selected row so the detail panel re-renders on language switch
        this.renderDetail(+tr.dataset.idx);
      });
    });
  },

  /* ---------------- 详情 ---------------- */
  // ---------------- Detail ----------------

  renderDetail(idx) {
    const res0 = (this.result.results || [])[this.resultIdx];
    const p = res0.pairs?.[idx];
    const el = document.getElementById("pr-detail");
    if (!p) return;
    const sp = p.specificity || {};
    const l = p.left || {}, rp = p.right || {};
    // 新引擎 amplifiable_pairs:{seq_id, product_len, f, r} 位点 dict
    // New engine amplifiable_pairs: {seq_id, product_len, f, r} site dicts
    let ofPairs = (sp.amplifiable_pairs || []).map((s) =>
      `<div class="spec-hit">✕ ${escapeHtml(s.seq_id)}: ${s.f?.start}-${s.r?.end} (${s.product_len} bp, F ${s.f?.strand === "-" ? "−" : "+"} R ${s.r?.strand === "-" ? "−" : "+"})</div>`).join("");
    // 旧引擎降级字段(带 qstart/sseqid 的列表形式)
    // Legacy fallback fields (list form with qstart/sseqid)
    if (!ofPairs && sp.amplifiable_pairs?.length) {
      ofPairs = sp.amplifiable_pairs.map((s) =>
        `<div class="spec-hit">✕ ${escapeHtml(s.subject)} ${s.start}-${s.end}</div>`).join("");
    }
    // 内部短种子命中:新引擎为 {forward, reverse} 计数;旧引擎为列表
    // Internal short-seed hits: new engine {forward, reverse} counts; legacy list
    let seedTxt = "";
    const sh = p.seed_hits;
    if (sh && typeof sh === "object" && !Array.isArray(sh) && (sh.forward !== undefined || sh.reverse !== undefined)) {
      seedTxt = `${t("design.col_fwd")} ${sh.forward ?? 0} / ${t("design.col_rev")} ${sh.reverse ?? 0}`;
    }
    // 新引擎 pair 级剖面分(0~100)
    // New engine pair-level profile scores (0-100)
    const profRow = (p.global_score !== undefined || p.three_prime_score !== undefined) ? `
        <dt>${t("design.col_profile")}</dt><dd>${p.global_score ?? "-"} / ${p.three_prime_score ?? "-"}</dd>` : "";
    const cov = p.covers_target === undefined ? ""
      : `<dt>${t("loc.covers")}</dt><dd>${p.covers_target === null ? t("loc.covers_na") : p.covers_target ? t("loc.covers_yes") : t("loc.covers_no")}${p.coverage_note ? ` — ${this.coverageNote(p)}` : ""}</dd>`;
    const absRow = (pr) => pr.abs_start ? ` · abs ${pr.abs_start}-${pr.abs_end}` : "";
    el.innerHTML = `
      <dl class="kv">
        <dt>${t("design.col_spec")}</dt><dd>${escapeHtml(this.specLabel(sp.label))} <span class="tag info">${sp.spec_score ?? "-"}</span>${p.level ? ` <span class="tag">${t("design.pair_level")} ${p.level}</span>` : ""}</dd>
        <dt>${t("design.col_offtarget")}</dt><dd>${sp.off_target_sites ?? 0}${p.binding_sites?.forward || p.binding_sites?.reverse ? ` · binding F ${p.binding_sites.forward.length} / R ${p.binding_sites.reverse.length}` : ""}</dd>
        <dt>${t("design.col_score")}</dt><dd><b>${p.composite_score ?? "-"}</b> = ${t("design.score_physical")} ${p.physical_score ?? "-"} × ${Math.round((this.result?.params?.score_physical_weight ?? 0.6) * 100)}% + ${t("design.score_spec")} ${sp.spec_score ?? "-"} × ${Math.round((this.result?.params?.score_specificity_weight ?? 0.4) * 100)}% ${(sp.off_target_sites || 0) ? ` − 10×${sp.off_target_sites}` : ""}</dd>
        ${profRow}
        <dt>${t("design.col_product")}</dt><dd>${p.product_len ?? "-"} bp</dd>
        ${cov}
        <dt>F</dt><dd class="mono">${escapeHtml(l.seq || "")} · start ${l.start}${absRow(l)} · Tm ${l.tm} · GC ${l.gc}% · ${t("design.hairpin", { h: l.hairpin ?? "-" })}</dd>
        ${rp.seq ? `<dt>R</dt><dd class="mono">${escapeHtml(rp.seq || "")} · start ${rp.start}${absRow(rp)} · Tm ${rp.tm} · GC ${rp.gc}% · ${t("design.hairpin", { h: rp.hairpin ?? "-" })}</dd>` : ""}
      </dl>
      ${sp.note ? `<div class="hint">${this.noteText(sp)}</div>` : ""}
      ${ofPairs ? `<div style="margin-top:6px"><b class="hint">${t("design.offtarget_sites_label")}</b>${ofPairs}</div>` : ""}
      ${seedTxt ? `<div style="margin-top:6px"><b class="hint">${t("design.seed_hits_label")}</b>${seedTxt}</div>` : ""}`;
  },

  /* ---------------- 项目保存/加载(guide 9.4,顶栏全局按钮) ---------------- */
  // ---------------- Project save/load (guide 9.4, global buttons in the top bar) ----------------

  async saveProject() {
    const options = {
      template: document.getElementById("dp-template").value,
      db: document.getElementById("dp-db").value,
      mode: document.getElementById("dp-mode").value,
      flank_len: parseInt(document.getElementById("dp-flank").value) || 0,
      params: this.collectParams(),
      named: this.serializeNamed(),
    };
    try {
      const payload = await api("/api/project/save", {
        method: "POST", json: true,
        body: {
          kind: "primer_design",
          db_prefix: document.getElementById("dp-db").value,
          options,
          // 项目文件同样瘦身(step1_hsps 渲染不读,大命中可达数十 MB);
          // kmer_cache 保留 —— 加载项目时回灌 sessionStorage 供再设计复用
          // Project files slim too (step1_hsps is never rendered and can hit
          // tens of MB); kmer_cache stays — loading pours it back into
          // sessionStorage so re-designing reuses it
          primer_design: { named: this.serializeNamed(), results: slimDesignResult(this.result), task_id: this.taskId },
        },
      });
      const blob = new Blob([JSON.stringify(payload, null, 2)], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `blastprime_design_${new Date().toISOString().slice(0, 10)}.json`;
      a.click();
      URL.revokeObjectURL(url);
      toast(t("design.project_saved"));
    } catch (e) { toast(e.message); }
  },

  /* 加载项目(顶栏 load-project 经 app.js dispatchLoadedProject 分发):
     回填模板/库/模式/参数/侧翼卡/设计结果;旧版定位项目已不支持 */
  // Load project (the top-bar load-project is dispatched via app.js dispatchLoadedProject):
  // backfills template/database/mode/params/flank card/design result; old locate projects are no longer supported
  async applyLoadedProject(data) {
    if (!data || data.kind !== "primer_design") {
      toast(t("top.unknown_project"));
      return;
    }
    const opts = data.options || {};
    const set = (id, v) => { if (v != null) document.getElementById(id).value = v; };
    if (opts.template != null) document.getElementById("dp-template").value = opts.template;
    set("dp-mode", opts.mode);
    set("dp-flank", opts.flank_len != null ? opts.flank_len : opts.flank);
    if (opts.db) {
      const sel = document.getElementById("dp-db");
      if ([...sel.options].some((o) => o.value === opts.db)) sel.value = opts.db;
      else if (sel.options.length) sel.value = sel.options[0].value;
    }
    if (opts.params) this.fillParams(opts.params);
    this.updateProductFields();
    const pd = data.primer_design || {};
    if (pd.locate) toast(t("design.old_project_locate"));
    const n = pd.named || opts.named || null;
    if (n) this.restoreNamed(n);
    // 恢复设计结果
    // Restore the design result
    this.result = pd.results || null;
    this.taskId = pd.task_id || null;
    this.resultIdx = 0;
    this.selPair = null;
    this.storeKmerCaches(this.result);   // 项目文件内的 k-mer 分析回灌 sessionStorage
    // Pour the project file's k-mer analyses back into sessionStorage
    if (this.result && (this.result.results || []).length) {
      document.getElementById("design-result").style.display = "";
      this.renderResult();
      await this.saveResult();   // 项目加载的结果也进入会话,切页不丢
      // Results loaded from a project also go into the session, so they are not lost on page switch
    } else {
      clearResult("bp_design_result");   // 新项目无结果 → 清掉旧会话结果(两处存储)
      // New project has no result -> clear the old session result (both stores)
    }
    this.saveState();
    toast(t("design.project_loaded"));
  },

  /* 跨页加载通道:app.js 跳转前写入 bp_design_load,本页消费后即删 */
  // Cross-page load channel: app.js writes bp_design_load before navigation; this page deletes it after consuming
  consumeLoadedProject() {
    let data = null;
    try { data = JSON.parse(localStorage.getItem("bp_design_load") || "null"); } catch (e) {}
    if (!data) return;
    localStorage.removeItem("bp_design_load");
    this.applyLoadedProject(data);
  },

  /* 主题切换 → 滑块 + 剖面图即时重绘(app.js applyTheme 钩子) */
  // Theme switch -> instant slider + profile plot re-render (app.js hook)
  onThemeChanged() {
    if (this.sliders.f) this.sliders.f.render();
    if (this.sliders.r) this.sliders.r.render();
    if (this.result) this.renderDepth();
  },

  /* 语言切换 → 动态渲染内容即时重绘(app.js setLang 钩子;静态 data-i18n 由 app.js 处理) */
  // Language switch -> re-render dynamic content (app.js setLang hook; static data-i18n elements are handled by app.js)
  /* 语言切换后重绘运行状态行:spinner/秒表/✓ 完成等是动态 innerHTML
     快照,setLang 只扫静态 data-i18n,切换语言后这里按当前任务状态重画 */
  // Re-render the run status line after a language switch: the spinner/
  // stopwatch/✓ done text are dynamic innerHTML snapshots that setLang
  // (static data-i18n only) cannot refresh, so redraw from task state here
  _renderRunStatus() {
    const s = this._runTaskId ? App.tasks.get(this._runTaskId) : null;
    if (!s) return;
    const status = document.getElementById("run-status");
    const hint = document.getElementById("run-hint-extra");
    const title = document.getElementById("run-hint");
    const fmt = (v) => (v < 60 ? v + "s" : Math.floor(v / 60) + "m" + (v % 60) + "s");
    const running = s.status === "running" || s.status === "pending";
    if (title) title.textContent = t("design.run_title") + " … " + this._runTaskId;
    if (!status) return;
    if (running) {
      if (hint) hint.textContent = t("design.run_elsewhere");
      const secs = Math.floor((Date.now() - (this._runStarted || Date.now())) / 1000);
      status.className = "build-status";
      status.innerHTML = `<span class="spinner"></span>${t("design.run_running", { s: fmt(secs) })}`;
    } else {
      if (hint) hint.textContent = "";
      const secs = this._runElapsed != null ? this._runElapsed : 0;
      status.className = "build-status" + (s.status === "succeeded" ? " done" : " err");
      status.innerHTML = s.status === "succeeded" ? `✓ ${t("design.run_done", { s: fmt(secs) })}`
        : s.status === "cancelled" ? `✗ ${t("design.run_canceled")}` : `✗ ${t("design.run_failed")}`;
    }
  },

  onLangChanged() {
    if (this.named.tmpl) this.renderTmplInfo();
    if (this.result) {
      this.renderResult();
      this.renderDepth();
      // 详情面板内容(特异性 note/覆盖备注等)依赖语言,选中行时一并重绘
      // The detail panel (specificity note / coverage note) is language-dependent:
      // re-render it when a row is selected
      if (this.selPair != null) this.renderDetail(this.selPair);
    }
    this._renderRunStatus();
    this.updateParamsCaret();
    this.updateKmerCacheCount();
  },

  /* 参数区收起提示(用户需求:可收起但无提示):▾ 展开 / ▸ 收起,
     悬停 title 说明可点击。setLang 只替换 data-i18n 文本节点,
     caret span 无 data-i18n,方向指示不被语言切换抹掉 */
  // Params collapse affordance (user request: collapsible but no hint):
  // ▾ expanded / ▸ collapsed, hover title explains clickability. setLang only
  // replaces data-i18n text nodes, so the caret span (no data-i18n) survives.
  updateParamsCaret() {
    const body = document.getElementById("params-body");
    const caret = document.getElementById("params-caret");
    if (caret) caret.textContent = body.style.display === "none" ? " ▸" : " ▾";
    document.getElementById("params-toggle").title = t("design.params_collapse");
  },
};
