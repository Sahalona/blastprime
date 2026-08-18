/* BLAST 比对页:运行/图形摘要/结果树/详情/原始输出懒加载/导入/导出/引物分析 */
// BLAST comparison page: run / graphical summary / result tree / details / lazy-loaded raw output / import / export / primer analysis
"use strict";

window.PAGE = {
  taskId: null,
  result: null,           // {raw_output, parsed, options}
  selHsp: null,           // {query_idx, subject_idx, hit_idx} — 单个 HSP
  // Single HSP
  selAgg: null,           // {type: "query"|"subject", qi, si?} — 聚合展示
  // Aggregated display
  hlMode: "mismatch",     // 高亮模式: mismatch(错配/缺口,默认) | match(完全一致)
  // Highlight mode: mismatch (mismatch/gap, default) | match (identical)
  vizScale: 1,
  vizQueryIdx: 0,         // 图形摘要/命中列表当前显示的查询(下拉联动)
  // Query currently shown in the graphical summary / hit list (linked to dropdown)
  vizRects: [],
  rawShown: 0,
  _anaToken: 0,           // 引物分析请求令牌(丢弃过期响应)
  // Primer analysis request token (stale responses are discarded)

  async init() {
    await this.loadDbs();
    this.restoreState();
    this.restoreResult();
    // 兜底持久化完成事件:切页回来时任务在后台完成、init 的同步 restoreResult
    // 早于异步写入 → 事件到达时再恢复渲染(仅当前无结果时)
    // Fallback-persist completion event: when a task finished while the page
    // was away, init's synchronous restoreResult ran before the async write
    // — the event triggers the recovery render (only when no result exists)
    window.addEventListener("bp:result-persisted", (e) => {
      const d = e && e.detail;
      if (d && d.kind === "blast" && !this.result) this.restoreResult();
    });
    this.bindStateSavers();
    this.bindRun();
    this.bindImport();
    this.bindExport();
    this.bindAnalysis();
    this.bindRemoteToggle();
    this.bindSeqTools();
    this.bindDetail();
    document.getElementById("q-file").addEventListener("change", async (e) => {
      // 大文件(>800KB)直接走文件上传通道:预读进文本框会让 options 字段撑爆 multipart 1MB 上限
      // Large files (>800 KB) go straight through the file-upload channel:
      // preloading into the text box would blow the 1 MB multipart field limit
      if (e.target.files.length && [...e.target.files].some((f) => f.size > 800 * 1024)) {
        toast(t("blast.big_file_note"));
        return; // 保留文件在 input 中,提交时随 multipart 上传
      }
      const { added, failed } = await appendFilesToTextarea(e.target, document.getElementById("q-seq"));
      if (added > 0) toast(t("blast.import_added").replace("{n}", added));
      if (failed > 0) toast(t("common.read_failed"));
    });
    this.initViz();
  },

  /* 详情控件:高亮模式切换 + 复制 Sbjct 序列(guide 5.4.3) */
  // Detail controls: highlight-mode toggle + copy Sbjct sequence (guide 5.4.3)
  bindDetail() {
    document.getElementById("hl-mode").addEventListener("change", (e) => {
      this.hlMode = e.target.value;
      this.renderDetail();
    });
    document.getElementById("btn-copy-sbjct").addEventListener("click", () => this.copySbjct());
    document.getElementById("btn-copy-name").addEventListener("click", () => this.copyNameQuery());
  },

  bindRun() {
    document.getElementById("btn-run").addEventListener("click", () => this.runBlast());
    document.getElementById("btn-restore-defaults")?.addEventListener("click", () => this.restoreDefaults());
    document.getElementById("btn-clear-result")?.addEventListener("click", () => this.clearResult());
  },

  /* ---------------- 页面状态临时保存(切页复原,见 app.js savePageState) ---------------- */
  // Temporary page-state persistence (restored on page switch-back; see app.js savePageState)

  saveState() {
    savePageState("blast", {
      seq: document.getElementById("q-seq").value,
      db: document.getElementById("q-db").value,
      program: document.getElementById("q-program").value,
      evalue: document.getElementById("q-evalue").value,
      maxtargets: document.getElementById("q-maxtargets").value,
      short: document.getElementById("q-short").checked,
      remote: document.getElementById("q-remote").checked,
      remoteDb: document.getElementById("q-remote-db").value,
    });
  },

  restoreState() {
    const s = loadPageState("blast");
    if (!s) return;
    const set = (id, v) => { if (v != null) document.getElementById(id).value = v; };
    set("q-seq", s.seq);
    set("q-program", s.program);
    set("q-evalue", s.evalue);
    set("q-maxtargets", s.maxtargets);
    document.getElementById("q-short").checked = !!s.short;
    document.getElementById("q-remote").checked = !!s.remote;
    set("q-remote-db", s.remoteDb);
    // 库下拉仅恢复仍存在的库
    // Restore only databases that still exist in the database dropdown
    const sel = document.getElementById("q-db");
    if (s.db && [...sel.options].some((o) => o.value === s.db)) sel.value = s.db;
    document.getElementById("q-remote-db").style.display = s.remote ? "" : "none";
  },

  bindStateSavers() {
    const save = () => this.saveState();
    document.getElementById("q-seq").addEventListener("input", save);
    ["q-db", "q-program", "q-evalue", "q-maxtargets", "q-short", "q-remote", "q-remote-db"]
      .forEach((id) => document.getElementById(id).addEventListener("change", save));
  },

  /* ---------------- 序列工具(guide.md 5.3) ---------------- */
  // Sequence tools (guide.md 5.3)

  bindSeqTools() {
    const apply = async (action) => {
      const ta = document.getElementById("q-seq");
      if (!ta.value.trim()) { toast(t("blast.seq_empty")); return; }
      try {
        const r = await api("/api/seq/" + action, { method: "POST", json: true, body: { text: ta.value } });
        ta.value = r.text;
        ta.focus();
      } catch (e) { toast(e.message); }
    };
    document.getElementById("btn-revcomp").addEventListener("click", () => apply("revcomp"));
    document.getElementById("btn-clean").addEventListener("click", () => apply("clean"));
    document.getElementById("btn-translate6").addEventListener("click", () => this.openTranslate6());
    document.querySelector("#translate-modal .btn-tr-close").addEventListener("click", () =>
      modalHide(document.getElementById("translate-modal")));
    document.getElementById("btn-tr-copy").addEventListener("click", () => this.copyTranslateAll());
    document.getElementById("tr-table").addEventListener("change", () => this.loadTranslate());
    document.getElementById("tr-seq").addEventListener("change", () => this.renderTranslate());
  },

  async openTranslate6() {
    const text = document.getElementById("q-seq").value.trim();
    if (!text) { toast(t("blast.seq_empty")); return; }
    this._trText = text;
    this._trSeqs = [];
    const tableSel = document.getElementById("tr-table");
    if (!tableSel.options.length) {
      try {
        const { tables } = await api("/api/seq/tables");
        tableSel.innerHTML = "";
        tables.forEach((tb) => {
          const opt = document.createElement("option");
          opt.value = tb.id;
          opt.textContent = `${tb.id}: ${tb.name}`;
          tableSel.appendChild(opt);
        });
      } catch (e) { toast(e.message); return; }
    }
    modalShow(document.getElementById("translate-modal"));
    await this.loadTranslate();
  },

  async loadTranslate() {
    const table = document.getElementById("tr-table").value;
    try {
      const { seqs } = await api("/api/seq/translate6", {
        method: "POST", json: true,
        body: { text: this._trText, table: parseInt(table) || 1 },
      });
      this._trSeqs = seqs;
      const seqSel = document.getElementById("tr-seq");
      const prev = seqSel.value;
      seqSel.innerHTML = "";
      seqs.forEach((s, i) => {
        const opt = document.createElement("option");
        opt.value = i;
        opt.textContent = `${s.name || "seq" + (i + 1)} (${s.length} bp)`;
        seqSel.appendChild(opt);
      });
      if (prev && seqSel.querySelector(`option[value="${prev}"]`)) seqSel.value = prev;
      this.renderTranslate();
    } catch (e) {
      document.getElementById("tr-out").innerHTML = `<div class="empty">${escapeHtml(e.message)}</div>`;
    }
  },

  renderTranslate() {
    const idx = parseInt(document.getElementById("tr-seq").value) || 0;
    const seq = this._trSeqs[idx];
    const out = document.getElementById("tr-out");
    if (!seq) { out.innerHTML = `<div class="empty">${t("blast.seq_empty")}</div>`; return; }
    const labels = { f1: "+1", f2: "+2", f3: "+3", r1: "-1", r2: "-2", r3: "-3" };
    out.innerHTML = Object.entries(labels).map(([k, lab]) => `
      <div style="margin-bottom:8px">
        <b class="muted">${escapeHtml(seq.name)} ${lab}</b>
        <pre class="code" style="white-space:pre-wrap; min-height:auto; max-height:110px; overflow:auto; margin-top:2px">${escapeHtml(seq.frames[k])}</pre>
      </div>`).join("");
  },

  copyTranslateAll() {
    const idx = parseInt(document.getElementById("tr-seq").value) || 0;
    const seq = this._trSeqs[idx];
    if (!seq) return;
    const labels = { f1: "+1", f2: "+2", f3: "+3", r1: "-1", r2: "-2", r3: "-3" };
    const lines = [];
    Object.entries(labels).forEach(([k, lab]) => lines.push(`>${seq.name} ${lab}`, seq.frames[k]));
    copyText(lines.join("\n"));
  },

  async initViz() {
    // 可视化(guide 5.4.2):Ctrl+滚轮(含 Linux Button4/5)横向缩放 + 拖拽/滚动条平移 + 悬停 + 点击
    // Visualization (guide 5.4.2): Ctrl+wheel (incl. Linux Button4/5) horizontal zoom + drag/scrollbar pan + hover + click
    const cv = document.getElementById("viz-canvas");
    const scroll = document.getElementById("viz-scroll");
    // 查询下拉:图形摘要与命中列表只显示所选查询(用户需求:多查询不上下堆叠)
    // Query dropdown: the graphical summary and hit list show only the selected query (user request: multiple queries not stacked vertically)
    document.getElementById("viz-query").addEventListener("change", (e) => {
      this.setVizQuery(parseInt(e.target.value) || 0);
    });
    // 上一个/下一个:不用展开下拉直接切换查询序列(用户需求)
    // Prev / Next: switch query sequence without opening the dropdown (user request)
    document.getElementById("btn-viz-prev").addEventListener("click", () => {
      this.setVizQuery(this.vizQueryIdx - 1);
    });
    document.getElementById("btn-viz-next").addEventListener("click", () => {
      this.setVizQuery(this.vizQueryIdx + 1);
    });
    // R12 分级切换:identity|score 双态按钮组(普通/短序列模式分开存)
    // Grade toggle (R12): identity|score two-state button group (normal /
    // short-seq modes stored separately)
    document.getElementById("btn-viz-grade-identity").addEventListener("click", () => {
      this.setVizGrade("identity");
    });
    document.getElementById("btn-viz-grade-score").addEventListener("click", () => {
      this.setVizGrade("score");
    });
    this.updateGradeButtons();
    this.updateVizNav();
    cv.addEventListener("wheel", (e) => {
      if (!e.ctrlKey || !this.result) return;
      e.preventDefault();
      this.vizScale = Math.min(8, Math.max(0.5, this.vizScale * (e.deltaY < 0 ? 1.2 : 1 / 1.2)));
      this.renderViz();
    }, { passive: false });
    // 纵向滚动时标尺固定在头部,仅同步横向偏移
    // On vertical scroll the ruler stays fixed in the header; only horizontal offset is synced
    scroll.addEventListener("scroll", () => this.syncRuler(scroll));
    cv.addEventListener("mousemove", (e) => this.vizHover(e));
    cv.addEventListener("mouseleave", () => {
      document.getElementById("viz-tip").style.display = "none";
      this._vizDrag = null;
      cv.style.cursor = "crosshair";
    });
    cv.addEventListener("click", (e) => this.vizClick(e));
    // 鼠标拖拽横向平移
    // Mouse drag pans horizontally
    cv.addEventListener("mousedown", (e) => {
      this._vizDrag = { startX: e.clientX, startLeft: scroll.scrollLeft, moved: false };
      cv.style.cursor = "grabbing";
    });
    window.addEventListener("mousemove", (e) => {
      if (!this._vizDrag) return;
      const dx = e.clientX - this._vizDrag.startX;
      if (Math.abs(dx) > 3) this._vizDrag.moved = true;
      scroll.scrollLeft = this._vizDrag.startLeft - dx;
    });
    window.addEventListener("mouseup", () => {
      if (this._vizDrag) { this._vizDrag = null; cv.style.cursor = "crosshair"; }
    });
    window.addEventListener("resize", () => { if (this.result) this.renderViz(); });
  },

  /* ---------------- 库与选项 ---------------- */
  // Databases and options

  async loadDbs() {
    const sel = document.getElementById("q-db");
    try {
      const { records } = await api("/api/db/records");
      // R12:不再有"选择本地数据库…"占位项;无库时放一个 disabled 提示项
      // R12: no "select a local database..." placeholder option; when there
      // are no records, a disabled hint option is shown instead
      sel.innerHTML = records.map((r) => {
        const name = r.prefix.split("/").pop();
        const label = r.note ? `${name}（${r.note}）` : name;
        return `<option value="${escapeHtml(r.prefix)}" title="${escapeHtml(r.prefix)}">${escapeHtml(label)}</option>`;
      }).join("");
      if (!records.length) {
        sel.innerHTML = `<option value="" disabled>${escapeHtml(t("blast.db_none"))}</option>`;
      }
    } catch (e) { /* 无库 */ }
    // No databases
  },

  bindRemoteToggle() {
    const remote = document.getElementById("q-remote");
    remote.addEventListener("change", () => {
      const rdb = document.getElementById("q-remote-db");
      rdb.style.display = remote.checked ? "" : "none";
    });
  },

  /* ---------------- 运行 ---------------- */
  // Run

  async runBlast() {
    // 重名标题加后缀(仅提交/解析时改写,输入框保持原文)
    // Duplicate headers get " (N)" suffixes at submit/parse time only (the
    // input box keeps the original text)
    const queryText = dedupeFastaHeaders(document.getElementById("q-seq").value).trim();
    const remote = document.getElementById("q-remote").checked;
    const db = document.getElementById("q-db").value;
    const remoteDb = document.getElementById("q-remote-db").value.trim();
    if (remote && !remoteDb) { toast(t("blast.remote_db_empty")); return; }
    // 名称型 query 自带 targetbase=(运行库由后端解析)时,页面未选库也能跑
    // With targetbase= in a name-style query the run db comes from the
    // backend, so the page db is optional there
    const namedDb = /,\s*targetbase=/i.test(queryText.split("\n")[0] || "");
    if (!remote && !db && !namedDb) { toast(t("blast.db_empty")); return; }
    if (!queryText && !document.getElementById("q-file").files.length) {
      toast(t("blast.seq_empty")); return;
    }
    // 非法字符检测(guide 5.3):序列行仅允许字母、*、-,否则提示先执行序列清洗
    // Invalid-character check (guide 5.3): sequence lines may contain only letters, *, -; otherwise prompt to run sequence cleaning first
    for (const line of queryText.split("\n")) {
      const s = line.trim();
      if (!s || s.startsWith(">") || s.startsWith("#")) continue;
      if (/[^A-Za-z*\-]/.test(s)) { toast(t("blast.warn_invalid")); return; }
    }
    // 短序列警告
    // Short-sequence warning
    const firstLine = queryText.split("\n")[0].trim();
    const seqOnly = queryText && !firstLine.startsWith(">");
    const seqLen = seqOnly ? queryText.replace(/[^ACGTUNacgtun]/g, "").length : 0;
    if (seqLen && seqLen <= 30 && !document.getElementById("q-short").checked) {
      if (!confirm(t("blast.warn_short"))) return;
    }
    // 远程警告
    // Remote-search warning
    if (remote && !confirm(t("blast.warn_remote"))) return;

    const options = {
      query_text: queryText,
      db, remote, remote_db: remoteDb,
      program: document.getElementById("q-program").value,
      evalue: document.getElementById("q-evalue").value,
      max_targets: parseInt(document.getElementById("q-maxtargets").value) || 500,
      short_seq_mode: document.getElementById("q-short").checked,
    };
    const fd = buildForm({ options: JSON.stringify(options) }, document.getElementById("q-file"));
    const btn = document.getElementById("btn-run");
    btn.disabled = true;
    try {
      const { task_id } = await api("/api/blast/run", { method: "POST", body: fd });
      this.taskId = task_id;
      document.getElementById("run-hint").textContent = t("blast.run_title") + " … " + task_id;
      await this.watchTask(task_id, () => {
        btn.disabled = false;
        this.taskId = task_id;
        this.fetchResult(task_id);
      });
    } catch (e) {
      btn.disabled = false;
      toast(e.message);
    }
  },

  /* 等待任务结束(SSE 已在 App 维护任务状态,轮询终态)。
     运行感知与建库页一致:按钮禁用 + spinner/秒表状态行 + "可继续其他页面"提示,
     大库/远程比对耗时较长,秒表让"仍在运行"可见可感。 */
  // Wait for the task to finish (SSE already tracks task state in App; poll for the terminal state).
  // Run feedback mirrors the database-build page: disabled button + spinner/elapsed status line +
  // "you may continue elsewhere" hint; long alignments (large DBs / remote) need the timer to
  // make "still running" visible and felt.
  async watchTask(taskId, onDone) {
    const status = document.getElementById("run-status");
    const hint = document.getElementById("run-hint-extra");
    const btn = document.getElementById("btn-run");
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
    if (hint) hint.textContent = t("blast.run_elsewhere");
    return new Promise((resolve) => {
      const tick = setInterval(async () => {
        const s = App.tasks.get(taskId);
        if (!s) return;
        if (s.status !== "running" && s.status !== "pending") {
          clearInterval(tick);
          if (btn) btn.disabled = false;
          if (hint) hint.textContent = "";
          const secs = Math.floor((Date.now() - started) / 1000);
          this._runElapsed = secs;   // 固定耗时:语言切换重绘不再随时间增长
          if (s.status === "succeeded") setStatus("done", `✓ ${t("blast.run_done", { s: fmt(secs) })}`);
          else if (s.status === "cancelled") setStatus("err", `✗ ${t("blast.run_canceled")}`);
          else setStatus("err", `✗ ${t("blast.run_failed")}`);
          onDone(s);
          resolve(s);
        } else {
          const secs = Math.floor((Date.now() - started) / 1000);
          setStatus("", `<span class="spinner"></span>${t("blast.run_running", { s: fmt(secs) })}`);
        }
      }, 400);
      setTimeout(() => clearInterval(tick), 30 * 60 * 1000);
    });
  },

  /* 语言切换后重绘运行状态行(动态 innerHTML 快照,setLang 不刷新) */
  // Re-render the run status line after a language switch (dynamic innerHTML
  // snapshots that setLang cannot refresh)
  _renderRunStatus() {
    const s = this._runTaskId ? App.tasks.get(this._runTaskId) : null;
    if (!s) return;
    const status = document.getElementById("run-status");
    const hint = document.getElementById("run-hint-extra");
    const title = document.getElementById("run-hint");
    const fmt = (v) => (v < 60 ? v + "s" : Math.floor(v / 60) + "m" + (v % 60) + "s");
    const running = s.status === "running" || s.status === "pending";
    if (title) title.textContent = t("blast.run_title") + " … " + this._runTaskId;
    if (!status) return;
    if (running) {
      if (hint) hint.textContent = t("blast.run_elsewhere");
      const secs = Math.floor((Date.now() - (this._runStarted || Date.now())) / 1000);
      status.className = "build-status";
      status.innerHTML = `<span class="spinner"></span>${t("blast.run_running", { s: fmt(secs) })}`;
    } else {
      if (hint) hint.textContent = "";
      const secs = this._runElapsed != null ? this._runElapsed : 0;
      status.className = "build-status" + (s.status === "succeeded" ? " done" : " err");
      status.innerHTML = s.status === "succeeded" ? `✓ ${t("blast.run_done", { s: fmt(secs) })}`
        : s.status === "cancelled" ? `✗ ${t("blast.run_canceled")}` : `✗ ${t("blast.run_failed")}`;
    }
  },

  async fetchResult(taskId) {
    try {
      // 大结果 JSON 在 Web Worker 中解析,避免主线程阻塞数秒(页面无响应)
      // Large result JSON is parsed in a Web Worker to avoid blocking the main thread for seconds (page unresponsive)
      let res = typeof Worker === "function"
        ? await this.parseResultWorker(`/api/tasks/${taskId}/result`)
        : await api(`/api/tasks/${taskId}/result`);
      if (typeof Worker === "function") {
        // Worker 返回 {ok, data}/{ok:false, error} 信封,成功时解开为接口响应本体
        // Worker returns an {ok, data}/{ok:false, error} envelope; unwrap to the API response body on success
        if (!res.ok) { toast(res.error || t("common.task_failed")); return; }
        res = res.data;
      }
      const { status, result } = res;
      if (status === "failed") { toast(result?.error || t("common.task_failed")); return; }
      if (status !== "succeeded" || !result) return;
      this.result = result;
      this.rawShown = 0;
      this.selHsp = null;
      this.selAgg = null;
      this.vizQueryIdx = 0;
      // 先显示结果区再渲染:display:none 时容器宽为 0,canvas 会按 0 宽绘制导致内容水平拉伸;
      // 且首帧容器宽不含垂直滚动条(cw 在设 canvas 高度前读取),推迟一帧等布局稳定再渲染
      // Show the result area before rendering: with display:none the container width is 0 and
      // the canvas drawn at 0 width stretches content horizontally; the first frame also lacks
      // the vertical scrollbar (cw is read before the canvas height is set) — defer one frame
      // until the layout stabilizes, then render
      document.getElementById("result-area").style.display = "";
      requestAnimationFrame(() => {
        this.renderAll();
        this.selectFirstHsp();   // R12:默认选中第一个 HSP
        this.updateAnalysisPanel();
      });
      await this.saveResult();   // 结果就绪即存会话(切页/刷新回来仍显示比对结果;await 确保写完再允许切页)
      // Save to session as soon as the result is ready (the comparison result survives page switches/refreshes; await so the write completes before any navigation)
    } catch (e) { toast(e.message); }
  },

  /* ---------------- 结果会话保存(切页/刷新回来恢复比对结果) ---------------- */
  // Result session persistence (restores the comparison result after page switch/refresh)

  async saveResult() {
    if (!this.result) return;
    // raw_output 全文可达数十 MB,超 localStorage 5 MB 配额时静默失败 →
    // 切页回来结果丢失;persistResult 超限回退 IndexedDB(无配额)
    // raw_output runs to tens of MB; silently blowing the 5 MB localStorage
    // quota used to lose the result on page switch; persistResult falls back
    // to IndexedDB (no quota) when localStorage overflows
    await persistResult("bp_blast_result", {
      _v: 1, taskId: this.taskId,
      result: this.result,              // {raw_output, parsed, options}
      rawShown: this.rawShown,
    });
  },

  async restoreResult() {
    const data = await loadResult("bp_blast_result");
    if (!data || data._v !== 1 || !data.result || !Array.isArray(data.result.parsed)) {
      await clearResult("bp_blast_result");
      return;
    }
    this.result = data.result;
    this.taskId = data.taskId || null;
    this.rawShown = data.rawShown || 0;
    this.selHsp = null;
    this.selAgg = null;
    this.vizQueryIdx = 0;
    // 任务状态标签:任务仍在会话内则显示终态,否则仅显示 id
    // Task status label: show the terminal state if the task is still in this session, otherwise just the id
    const st = data.taskId && App.tasks.get(data.taskId);
    const hint = st && st.status === "succeeded" ? t("blast.done") + " · " + data.taskId
      : (st && st.status === "failed" ? t("blast.failed") + " · " + data.taskId
        : t("blast.run_title") + " … " + (data.taskId || ""));
    document.getElementById("run-hint").textContent = hint;
    document.getElementById("result-area").style.display = "";
    requestAnimationFrame(() => {
      this.renderAll();
      this.selectFirstHsp();   // R12:默认选中第一个 HSP
      this.updateAnalysisPanel();
    });
  },

  /* Worker 解析结果 JSON;失败时回退主线程解析 */
  // Parse result JSON in a Worker; fall back to main-thread parsing on failure
  parseResultWorker(url) {
    return new Promise((resolve, reject) => {
      try {
        const w = new Worker("/static/js/parse-worker.js");
        w.onmessage = (ev) => { w.terminate(); resolve(ev.data); };
        w.onerror = (e) => { w.terminate(); reject(e); };
        w.postMessage({ url });
      } catch (e) { reject(e); }
    });
  },

  /* ---------------- 结果渲染 ---------------- */
  // Result rendering

  renderAll() {
    this.fillVizQuery();
    this.renderViz();
    this.renderTree();
    this.renderDetail();
    this.renderRaw();
  },

  /* 查询下拉填充(与图形摘要/命中列表联动,见 initViz change) */
  // Fill the query dropdown (linked to the graphical summary / hit list; see initViz change)
  fillVizQuery() {
    const sel = document.getElementById("viz-query");
    const parsed = this.result?.parsed || [];
    sel.innerHTML = "";
    parsed.forEach((q, i) => {
      const opt = document.createElement("option");
      opt.value = String(i);
      const nm = q.name || `Query ${i + 1}`;
      opt.textContent = nm.length > 42 ? nm.slice(0, 42) + "…" : nm;
      opt.title = nm;
      sel.appendChild(opt);
    });
    if (this.vizQueryIdx >= parsed.length) this.vizQueryIdx = 0;
    sel.value = String(this.vizQueryIdx);
    this.updateVizNav();
  },

  /* 切换显示查询(下拉 change / 上一个 / 下一个共用):钳制边界、下拉同步、
     清选择后重渲染图形摘要/命中列表/详情三视图 */
  // Switch the displayed query (dropdown change / prev / next share this):
  // clamp to bounds, sync the dropdown, clear the selection, re-render all three panes
  setVizQuery(i) {
    const parsed = this.result?.parsed || [];
    if (!parsed.length) return;
    const idx = Math.max(0, Math.min(parsed.length - 1, i));
    if (idx === this.vizQueryIdx) { this.updateVizNav(); return; }
    this.vizQueryIdx = idx;
    const sel = document.getElementById("viz-query");
    if (sel) sel.value = String(idx);
    this.selHsp = null;
    this.selAgg = null;
    this.renderViz();
    this.renderTree();
    this.renderDetail();
    this.updateVizNav();
  },

  /* 上一个/下一个按钮可用态:无结果或单查询禁用;边界处禁用同向按钮 */
  // Prev/Next enabled state: disabled with no result or a single query;
  // the boundary-direction button is disabled at either end
  updateVizNav() {
    const parsed = this.result?.parsed || [];
    const prev = document.getElementById("btn-viz-prev");
    const next = document.getElementById("btn-viz-next");
    if (!prev || !next) return;
    prev.disabled = parsed.length <= 1 || this.vizQueryIdx <= 0;
    next.disabled = parsed.length <= 1 || this.vizQueryIdx >= parsed.length - 1;
  },

  /* R12 默认选中第一个 HSP:parsed 首个含 subjects 的查询的
     {query_idx, subject_idx:0, hit_idx:0} → 树高亮(scrollIntoView)+ 详情
     + 图形摘要选中环;无命中则无操作 */
  // Select the first HSP by default (R12): the first query with subjects,
  // subject 0, hit 0 -> tree highlight (scrollIntoView) + details + selection
  // ring in the graphic summary; no-op when there are no hits
  selectFirstHsp() {
    const parsed = this.result?.parsed || [];
    for (let qi = 0; qi < parsed.length; qi++) {
      const subs = parsed[qi].subjects || [];
      if (!subs.length) continue;
      this.vizQueryIdx = qi;
      this.selAgg = null;
      this.selHsp = { query_idx: qi, subject_idx: 0, hit_idx: 0 };
      const sel = document.getElementById("viz-query");
      if (sel) sel.value = String(qi);
      this.renderViz();
      this.renderTree(true);
      this.renderDetail();
      this.updateVizNav();
      return;
    }
  },

  /* score 保留一位小数显示:后端 float 解析后 JSON 会把 50.0 变 50,末尾 0 不能省略 */
  // Display score with one decimal place: after backend float parsing, JSON turns 50.0 into 50; the trailing 0 must not be dropped
  fmtScore(s) {
    return s == null || s === "" ? "" : Number(s).toFixed(1);
  },

  /* 条带颜色按 identity 分级(用户 2026-08-14 确认,替换 guide.md 的 score 分级):
     ≥95% 红 / ≥80% 品红 / ≥50% 绿 / ≥40% 蓝 / 其余黑。
     tblastn/blastp 用 BLOSUM62 打分,保守替换也给正分,identity 只有 ~50% 的命中
     bit score 也可能 ≥200 —— 按 score 分级会把半错配条带标红,误导用户。 */
  // Band color graded by identity (user-confirmed 2026-08-14, replacing the guide.md score grading):
  // >=95% red / >=80% magenta / >=50% green / >=40% blue / else black.
  // tblastn/blastp score via BLOSUM62, where conservative substitutions still score positively,
  // so ~50%-identity hits can reach bit score >=200 -- score grading paints half-mismatched bands red.
  /* 条带颜色按 identity 分级(用户 2026-08-14 确认,替换 guide.md 的 score 分级):
     ≥95% 红 / ≥80% 品红 / ≥50% 绿 / ≥40% 蓝 / 其余黑。 */
  // Band color graded by identity (user-confirmed 2026-08-14, replacing the guide.md score grading)
  identityColor(f) {
    const v = Number(f) || 0;
    if (v >= 0.95) return "#dc2626";
    if (v >= 0.80) return "#c026d3";
    if (v >= 0.50) return "#16a34a";
    if (v >= 0.40) return "#2563eb";
    return "#3f3f46";
  },

  /* 旧 score 分级兜底(数据缺少 identity 时;阈值同 guide.md L299) */
  // Legacy score grading fallback (when identity is missing; thresholds per guide.md L299)
  scoreColor(score) {
    const s = Number(score) || 0;
    if (s >= 200) return "#dc2626";
    if (s >= 80) return "#c026d3";
    if (s >= 50) return "#16a34a";
    if (s >= 40) return "#2563eb";
    return "#3f3f46";
  },

  /* HSP identity 分数(0~1):outfmt0 解析有 ident_frac,tabular 有 pident(0~100),都没有返回 null */
  // HSP identity fraction (0~1): ident_frac from outfmt0 parsing, pident (0~100) from tabular; null when neither exists
  hspIdent(h) {
    if (h && h.ident_frac != null) return Number(h.ident_frac);
    if (h && h.pident != null) return Number(h.pident) / 100;
    return null;
  },

  /* 条带色统一入口:分级模式 = identity 时 identity 优先、缺失回退 score;
     score 模式直走 score 分级(R12 切换按钮) */
  // Unified band-color entry: identity mode grades by identity with score
  // grading as fallback when identity is missing; score mode grades by score
  hspColor(h) {
    const idf = this.hspIdent(h);
    if (this.vizGrade() === "score") return this.scoreColor(h.score ?? h.bitscore);
    return idf != null ? this.identityColor(idf) : this.scoreColor(h.score ?? h.bitscore);
  },

  /* R12 分级模式:普通模式默认 identity,短序列模式默认 score;两键分开
     存储(经 savePageState 助手,即 localStorage) */
  // Grading mode (R12): identity by default in normal mode, score by default
  // in short-seq mode; two separate keys (localStorage via savePageState)
  vizGrade() {
    const short = document.getElementById("q-short")?.checked;
    const saved = short ? loadPageState("viz_grade_short") : loadPageState("viz_grade");
    return saved || (short ? "score" : "identity");
  },

  /* 切换分级模式:写入当前生效的键(普通/短序列分开存),刷新按钮态与提示,
     重绘图形摘要与命中列表(条带色随模式变化) */
  // Switch the grading mode: persist to the currently active key (normal /
  // short-seq stored separately), refresh the buttons and hint, re-render the
  // graphic summary and hit list (band colors change with the mode)
  setVizGrade(grade) {
    const short = document.getElementById("q-short")?.checked;
    if (short) savePageState("viz_grade_short", grade);
    else savePageState("viz_grade", grade);
    this.updateGradeButtons();
    this.renderViz();
    this.renderTree();
  },

  /* 分级按钮组 + 提示行状态(init / 切换 / 语言切换 / 短序列模式切换共用) */
  // Grade toggle group + hint state (shared by init / toggle / lang switch /
  // short-mode toggle)
  updateGradeButtons() {
    const idBtn = document.getElementById("btn-viz-grade-identity");
    const scBtn = document.getElementById("btn-viz-grade-score");
    const hint = document.getElementById("viz-hint");
    if (idBtn) idBtn.classList.toggle("active", this.vizGrade() === "identity");
    if (scBtn) scBtn.classList.toggle("active", this.vizGrade() === "score");
    if (hint) hint.textContent = this.vizGrade() === "identity"
      ? t("blast.viz_hint_identity") : t("blast.viz_hint_score");
  },

  /* R12 恢复默认参数:程序自动/e-value 10/最大命中 500/短序列关/远程关。
     仅参数,序列与库保留(用户 2026-08-15 确认);关短序列后分级模式随 q-short
     回到普通模式的默认(identity) */
  // Restore default parameters (R12): program auto / e-value 10 / max targets 500 /
  // short-seq off / remote off. Parameters only — sequence and database are kept
  // (user-confirmed 2026-08-15); with short-seq off the grading mode returns to the
  // normal-mode default (identity)
  restoreDefaults() {
    document.getElementById("q-program").value = "auto";
    document.getElementById("q-evalue").value = "10";
    document.getElementById("q-maxtargets").value = "500";
    document.getElementById("q-short").checked = false;
    document.getElementById("q-remote").checked = false;
    document.getElementById("q-remote-db").value = "";
    document.getElementById("q-remote-db").style.display = "none";
    this.saveState();
    this.updateGradeButtons();
    toast(t("blast.params_restored"));
  },

  /* R12 清空结果:结果状态与全部结果面板置空,隐藏结果区,移除持久化结果,
     复位运行状态行;序列/库/参数不受影响 */
  // Clear the result (R12): null the result state and all result panes, hide the
  // result area, remove the persisted result, reset the run-status row;
  // sequence / database / parameters are untouched
  clearResult() {
    this.result = null;
    this.taskId = null;
    this.rawShown = 0;
    this.selHsp = null;
    this.selAgg = null;
    this.vizQueryIdx = 0;
    this.vizRects = [];
    this._report = null;
    clearResult("bp_blast_result");   // 新旧结果都清(localStorage + IndexedDB)
    ["viz-query", "result-tree", "detail-view", "raw-view"].forEach((id) => {
      const el = document.getElementById(id);
      if (el) el.innerHTML = "";
    });
    ["viz-tip", "viz-empty", "result-tabs", "tab-analysis"].forEach((id) => {
      const el = document.getElementById(id);
      if (el) el.style.display = "none";
    });
    const tabs = document.getElementById("tab-results");
    if (tabs) tabs.style.display = "";
    const area = document.getElementById("result-area");
    if (area) area.style.display = "none";
    // 图形摘要画布/头部与无命中占位一起复位
    // Reset the graphic-summary canvas/head along with the no-hits placeholder
    const cv = document.getElementById("viz-canvas");
    if (cv) cv.style.display = "none";
    const head = document.getElementById("viz-head");
    if (head) head.style.display = "none";
    ["ana-cards", "ana-verdict", "ana-hits"].forEach((id) => {
      const el = document.getElementById(id);
      if (el) el.innerHTML = "";
    });
    const notes = document.getElementById("ana-notes");
    if (notes) notes.textContent = "";
    const rs = document.getElementById("run-status");
    if (rs) { rs.className = "build-status"; rs.innerHTML = ""; }
    const hint = document.getElementById("run-hint");
    if (hint) hint.textContent = "";
    const extra = document.getElementById("run-hint-extra");
    if (extra) extra.textContent = "";
    this.updateVizNav();
    toast(t("blast.result_cleared"));
  },

  /* --- 图形摘要 (guide 5.4.2) ---
     结构:固定头部(#viz-head:查询信息 + 坐标轴标尺,不随纵向滚动;标尺画布
     随 scrollLeft 同步平移)+ 可滚动条带区。横向缩放扩展内容宽度,滚动条平移。
     行布局(用户 2026-08-12 澄清):覆盖度条带位于画面层级最底(先画,被其它
     所有遮挡),垂直居中与基因名重叠;名称下方垫 65% 透明底框(同画布颜色,
     左缘对齐 0% 即 PAD,宽=文字宽+10)。条带右缘 x2 = PAD + qe*unit,
     全长命中恰抵 100% 参考线;选中改 --bg 外环 + --fg 内环自适应双描边。 */
  // Graphical summary (guide 5.4.2): fixed header (#viz-head: query info + axis ruler, not
  // scrolling vertically; the ruler canvas translates with scrollLeft) + scrollable band area.
  // Horizontal zoom widens the content; the scrollbar pans. Row layout (clarified by the user
  // 2026-08-12): the coverage band sits at the bottom z-order (drawn first, covered by all
  // others), vertically centered and overlapping the gene name; a 65% transparent backing box
  // (canvas color, left edge at 0% i.e. PAD, width = text width + 10) sits under the name.
  // Band right edge x2 = PAD + qe*unit, so a full-length hit exactly reaches the 100% line;
  // selection uses an adaptive double stroke: --bg outer ring + --fg inner ring.
  renderViz() {
    const cv = document.getElementById("viz-canvas");
    const scroll = document.getElementById("viz-scroll");
    const empty = document.getElementById("viz-empty");
    const headInfo = document.getElementById("viz-head-info");
    const parsed = this.result?.parsed || [];
    if (!parsed.length) {           // 无命中占位
    // No-hits placeholder
      cv.style.display = "none";
      document.getElementById("viz-head").style.display = "none";
      empty.style.display = "";
      return;
    }
    cv.style.display = "";
    document.getElementById("viz-head").style.display = "";
    empty.style.display = "none";
    const dpr = window.devicePixelRatio || 1;
    const cw = scroll.clientWidth || cv.parentElement.clientWidth;
    if (!cw) return;                // 容器不可见(宽 0)时延后渲染,避免 0 宽绘制拉伸
    // Container invisible (width 0): defer rendering to avoid 0-width drawing that stretches
    const BAND = 20, PAD = 14, RULER_H = 14;   // 行高:条带间距=BAND-14(名称框底 bandY+15 须 < BAND)
    // Row height: band spacing = BAND-14 (name-box bottom bandY+15 must be < BAND)
    const BAND_H = 14, LANE_GAP = 2;           // 条带高度 + 同基因多 HSP 分道间距
    // Band height + lane spacing for multiple HSPs of the same gene
    const usable = cw - 2 * PAD;
    const scale = this.vizScale;
    // 内容宽度随缩放扩展,外层滚动容器横向滚动(guide:缩放后标尺与条带同步)
    // Content width grows with zoom; the outer scroll container pans horizontally (guide: ruler and bands stay in sync after zoom)
    const contentW = Math.max(cw, Math.ceil(PAD * 2 + usable * scale));
    // 仅绘制下拉选中的查询(用户需求:多查询图形摘要不上下堆叠)
    // Draw only the query selected in the dropdown (user request: no vertical stacking of multiple-query summaries)
    const qi = this.vizQueryIdx < parsed.length ? this.vizQueryIdx : 0;
    const q = parsed[qi];
    // 同基因 HSP 分道:按查询坐标起点排序,与前一 HSP 重叠(查询坐标)则开新道,
    // 否则复用最先结束的道 → 重叠 HSP 错开不堆叠。行高随道数扩展(BAND + 额外道×(高+间距))
    // Same-gene HSP lane allocation: iterated in PARSED order (score desc, the same order as
    // the result tree) so lane 0 = the tree's first HSP — previously sorted by query-coordinate
    // start, which made the topmost band disagree with the tree whenever score order and
    // q-coordinate order diverged. An overlap (in query coords) with an earlier HSP opens a new
    // lane, otherwise the earliest-ending lane is reused, so overlapping HSPs stagger without
    // stacking. Row height grows with lane count (BAND + extra lanes x (height + gap))
    const layout = q.subjects.map((s) => {
      const items = s.hits.map((h, hi) => ({ h, hi,
        qs: Math.min(h.qstart, h.qend), qe: Math.max(h.qstart, h.qend) }));
      const lanes = [];
      const laneEnd = [];            // 每条道的实际末端(max qe)——输入不再按 qs 排序,
      // 不能只看末条插入项;The real end (max qe) per lane — input is no longer qs-sorted,
      // so the last-inserted item's end is not the lane's end
      items.forEach((it) => {
        // 优先复用末端 < 新条起点的最早空闲道;否则开新道
        // Reuse the first lane whose end is before the new item's start, else open a new lane
        let li = -1;
        for (let i = 0; i < lanes.length; i++) {
          if (laneEnd[i] < it.qs) { li = i; break; }
        }
        if (li < 0) { li = lanes.length; lanes.push([]); laneEnd.push(0); }
        lanes[li].push(it);
        laneEnd[li] = Math.max(laneEnd[li], it.qe);
      });
      return { rowH: BAND + Math.max(0, lanes.length - 1) * (BAND_H + LANE_GAP), lanes };
    });
    const H = layout.reduce((a, r) => a + r.rowH, 0) + 6;
    cv.width = contentW * dpr; cv.height = H * dpr;
    cv.style.width = contentW + "px"; cv.style.height = H + "px";
    const ctx = cv.getContext("2d");
    ctx.scale(dpr, dpr);
    ctx.clearRect(0, 0, contentW, H);
    this.vizRects = [];
    const colors = getComputedStyle(document.documentElement);
    const qlen = q.length || 1;
    const unit = usable * scale / qlen;
    const span = usable * scale;                       // 缩放后 0~100% 的绘制宽度
    // Draw width for 0-100% after zooming
    // 固定头部:查询信息(文字,不走滚动)
    // Fixed header: query info (text, not scrolled)
    headInfo.textContent = `${q.name}  (${qlen} nt) · ${Math.round(scale * 100)}%`;
    // 坐标轴标尺画布(与条带共享坐标空间,随横向 scrollLeft 平移,纵向冻结)
    // Axis ruler canvas (shares the coordinate space with the bands; translates with horizontal scrollLeft, frozen vertically)
    const rulerCv = document.getElementById("viz-ruler-canvas");
    rulerCv.width = contentW * dpr; rulerCv.height = RULER_H * dpr;
    rulerCv.style.width = contentW + "px"; rulerCv.style.height = RULER_H + "px";
    const rctx = rulerCv.getContext("2d");
    rctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    rctx.clearRect(0, 0, contentW, RULER_H);
    // 0~100% 等分刻度线
    // Evenly spaced tick marks from 0 to 100%
    rctx.strokeStyle = colors.getPropertyValue("--border-strong");
    rctx.beginPath(); rctx.moveTo(PAD, 3); rctx.lineTo(PAD + span, 3); rctx.stroke();
    rctx.fillStyle = colors.getPropertyValue("--fg-faint");
    rctx.font = "9px sans-serif";
    [0, 0.25, 0.5, 0.75, 1].forEach((frac) => {
      const x = PAD + span * frac;
      rctx.beginPath(); rctx.moveTo(x, 3); rctx.lineTo(x, 6); rctx.stroke();
      rctx.textAlign = "center";
      rctx.fillText(Math.round(frac * 100) + "%", x, RULER_H - 1);
    });
    this.syncRuler(scroll, rulerCv);
    // 第一遍:全部条带(层级最底,被名称底框与文字遮挡)。同基因多道:**整块条带块**垂直
    // 居中,第 li 道 = 块顶 + li×(高+间距);单 HSP 行与旧版一致(14px 居中)。
    // 关键:此前只居中第 0 道,多道块的底边超出行高(2 道超 5px、3 道超 13px),
    // 下道条带溢出伸进下一行,与下一行条带粘连成"堆叠"视觉效果
    // Pass 1: all bands (bottom z-order, covered by the name box and text). For multiple
    // lanes of the same gene the whole BLOCK is vertically centered; lane li = block top +
    // li x (height + gap); single-HSP rows match the old version (14px, centered).
    // Key: only lane 0 used to be centered, so the block bottom exceeded the row height
    // (5px for 2 lanes, 13px for 3), letting the lower lane spill into the next row and
    // visually glue the bands together as "stacking"
    let yy = 0;
    q.subjects.forEach((s, si) => {
      const L = layout[si];
      const center = yy + L.rowH / 2;
      // 块高 = 14 + (n-1)×16;块顶使整块居中,行内上/下余量各 3px(与单道一致)
      // Block height = 14 + (n-1)x16; the block top centers the whole block, leaving a
      // 3px margin above and below inside the row (same as single-lane rows)
      const blockTop = center - (BAND_H + (L.lanes.length - 1) * (BAND_H + LANE_GAP)) / 2;
      L.lanes.forEach((lane, li) => {
        const ly = blockTop + li * (BAND_H + LANE_GAP);
        lane.forEach((it) => {
          const x1 = PAD + (it.qs - 1) * unit;
          const x2 = PAD + it.qe * unit;               // 右缘在 qe 处,PAD+span 恰为 100% 线
          // Right edge at qe; PAD+span is exactly the 100% line
          const col = this.hspColor(it.h);
          ctx.fillStyle = col;
          ctx.fillRect(Math.min(x1, x2), ly, Math.max(2, Math.abs(x2 - x1)), BAND_H);
          this.vizRects.push({ qi, si, hi: it.hi, h: it.h, x: Math.min(x1, x2), y: ly, w: Math.max(2, Math.abs(x2 - x1)), hh: BAND_H });
        });
      });
      yy += L.rowH;
    });
    // 第二遍:选中条带高亮环——第四遍(最顶层)绘制,见下方 pass 4
    // Pass 2: the selected-band highlight ring is drawn in pass 4 (topmost), below
    // 第三遍:基因名 + 65% 透明底框(压在条带与高亮环之上,重合时名称仍可看清)
    // 框与文字中点对齐第 0 道条带中点(多道时放在最上面那条带上,下方条带不被遮挡;
    // 单道时第 0 道即整块,与旧版一致)。随行高变化与条带一起移动
    // 基因名是序列名,用等宽字体(--mono)显示
    // Pass 3: gene name + 65% transparent backing box (on top of the bands and highlight ring;
    // the name stays legible where they overlap). Box and text centers align with lane 0 of the
    // band block (for multiple lanes the box sits on the topmost band so the lower bands stay
    // unobstructed; for a single lane lane 0 is the whole block, matching the old version).
    // Box and text move with the bands as the row height changes. The gene name is a sequence
    // name, shown in monospace (--mono)
    ctx.font = `10px ${colors.getPropertyValue("--mono")}`;
    ctx.textBaseline = "middle";
    let yy3 = 0;
    q.subjects.forEach((s, si) => {
      const L = layout[si];
      const center = yy3 + L.rowH / 2;               // 与第一遍条带块相同
      // Same as the band block in pass 1
      // 第 0 道中点 = 块顶 + 道高/2;单道时 = 行中点,多道时上移至第一条带
      // Lane-0 center = block top + band height/2; = row center for a single lane,
      // shifted up onto the first band when there are multiple lanes
      const blockTop = center - (BAND_H + (L.lanes.length - 1) * (BAND_H + LANE_GAP)) / 2;
      const boxC = blockTop + BAND_H / 2;
      const name = s.name || "";
      const tw = ctx.measureText(name).width;
      ctx.globalAlpha = 0.65;
      ctx.fillStyle = colors.getPropertyValue("--bg-soft");
      ctx.fillRect(PAD, boxC - 8, tw + 10, 16);      // [boxC-8, boxC+8],中点=第 0 道中点
      // [boxC-8, boxC+8]; center = lane-0 center
      ctx.globalAlpha = 1;
      ctx.fillStyle = colors.getPropertyValue("--fg");
      ctx.textAlign = "left";
      ctx.fillText(name, PAD + 5, boxC);             // 文字中点=第 0 道中点
      // Text center = lane-0 center
      yy3 += L.rowH;
    });
    // 第四遍(最顶层):选中条带高亮环——画在所有条带与基因名框之上,保证高亮框
    // 完整可见,不被下方相邻条带或它旁边的基因名框遮挡(此前画在名称框之下,
    // 环顶/环底与相邻条带重叠、重叠部分被名称框底衬淡化,环显得被"挡住")。
    // 环中间透明,基因名文字仍清晰;环线压过相邻条带顶部或名称框边缘是环完整
    // 的必然代价(环优先级最高)。
    // Pass 4 (topmost): highlight ring for the selected band — drawn above all bands and
    // the gene-name box so it stays fully visible, never obscured by a band below it or
    // by the name box beside it (previously drawn under the box, the ring's top/bottom
    // strokes overlapped the adjacent band and the box backing faded them, so the ring
    // looked "blocked"). The ring interior is transparent so the name text stays legible;
    // ring lines crossing a neighboring band's top edge or the box edge are the necessary
    // price of ring completeness (the ring has topmost priority)
    if (this.selHsp && this.selHsp.query_idx === qi) {
      const hit = this.vizRects.find((r) => r.qi === qi
        && r.si === this.selHsp.subject_idx && r.hi === this.selHsp.hit_idx);
      if (hit) {
        ctx.strokeStyle = "#000000";          // 外层:黑色高亮框(先画,被白环盖住内侧)
        // Outer: black highlight box (drawn first; its inner side is covered by the white ring)
        ctx.lineWidth = 2;
        ctx.strokeRect(hit.x - 3, hit.y - 3, hit.w + 6, hit.hh + 6);
        ctx.strokeStyle = "#FFFFFF";          // 内层:白色过渡间隙 1px,与条带本体隔开
        // Inner: 1px white transition gap separating it from the band body
        ctx.lineWidth = 1;
        ctx.strokeRect(hit.x - 1, hit.y - 1, hit.w + 2, hit.hh + 2);
        ctx.lineWidth = 1;
      }
    }
    ctx.textBaseline = "alphabetic";
    ctx.font = "11px sans-serif";
  },

  /* 标尺画布随横向滚动同步平移(纵向冻结在固定头部,不随 #viz-scroll 滚动) */
  // The ruler canvas translates with horizontal scrolling (vertically frozen in the fixed header, not scrolling with #viz-scroll)
  syncRuler(scroll, rulerCv) {
    if (!rulerCv) rulerCv = document.getElementById("viz-ruler-canvas");
    rulerCv.style.transform = `translateX(${-scroll.scrollLeft}px)`;
  },

  /* 画布内容坐标:getBoundingClientRect 已含滚动偏移,
     再加 scrollTop/scrollLeft 会双倍计数导致滚动后点击/悬停错位 */
  // Canvas content coordinates: getBoundingClientRect already includes scroll offsets;
  // adding scrollTop/scrollLeft again double-counts them and misaligns click/hover after scrolling
  vizPos(e) {
    const rect = e.target.getBoundingClientRect();
    return {
      x: e.clientX - rect.left,
      y: e.clientY - rect.top,
    };
  },

  vizHover(e) {
    const tip = document.getElementById("viz-tip");
    const { x: mx, y: my } = this.vizPos(e);
    const hit = this.vizRects.find((r) => mx >= r.x && mx <= r.x + r.w && my >= r.y && my <= r.y + r.hh);
    if (!hit) { tip.style.display = "none"; if (!this._vizDrag) e.target.style.cursor = "crosshair"; return; }
    const h = hit.h;
    if (!this._vizDrag) {
      tip.style.display = "block";  // 样式表默认 display:none,空串会回落为隐藏
      // Stylesheet default is display:none; an empty string would fall back to hidden
      // 提示浮动保持在可视区内(相对固定包装层定位,不随滚动移走)
      // Tooltip floats within the viewport (positioned relative to the fixed wrapper, not carried away by scrolling)
      const wRect = document.getElementById("viz-wrap").getBoundingClientRect();
      // 提示浮动保持在可视区内,置于鼠标右下方(用户需求)
      // Tooltip stays inside the viewport, placed to the bottom-right of the mouse (user request)
      tip.style.left = Math.min(Math.max(e.clientX - wRect.left + 14, 4), wRect.width - 240) + "px";
      tip.style.top = Math.max(e.clientY - wRect.top + 12, 4) + "px";
      tip.textContent =
        `${hit.h.sseqid} — score ${this.fmtScore(h.score ?? h.bitscore)} · E-value ${h.evalue} · ` +
        `q ${h.qstart}-${h.qend} · s ${h.sstart}-${h.send}`;
      e.target.style.cursor = "pointer";
    }
  },

  vizClick(e) {
    if (this._vizDrag?.moved) return;   // 拖拽后不触发选中
    // Do not trigger selection after a drag
    const { x: mx, y: my } = this.vizPos(e);
    const hit = this.vizRects.find((r) => mx >= r.x && mx <= r.x + r.w && my >= r.y && my <= r.y + r.hh);
    if (!hit) return;
    this.selAgg = null;
    this.selHsp = { query_idx: hit.qi, subject_idx: hit.si, hit_idx: hit.hi };
    this.renderViz();
    this.renderTree(true);
    this.renderDetail();
  },

  /* --- 结果树 --- */
  // --- Result tree ---
  renderTree(scrollToHsp) {
    const root = document.getElementById("result-tree");
    const parsed = this.result?.parsed || [];
    if (!parsed.length) { root.innerHTML = `<div class="empty">${t("blast.result_empty")}</div>`; return; }
    // 与图形摘要联动:只显示下拉选中的查询(用户需求:命中列表不全部堆叠)
    // Linked with the graphical summary: only the dropdown-selected query is shown (user request: hit list not fully stacked)
    const qi = this.vizQueryIdx < parsed.length ? this.vizQueryIdx : 0;
    const q = parsed[qi];
    let html = "";
    const qAgg = this.selAgg?.type === "query" && this.selAgg.qi === qi;
    html += `<div class="tree-node q${qAgg ? " sel" : ""}" data-qi="${qi}"><span class="caret">▼</span>${escapeHtml(q.name)}
      <span class="muted"> · ${(q.subjects || []).length} ${t("blast.hits_count")}</span></div>`;
    html += `<div class="q-body" data-qi="${qi}">`;
    (q.subjects || []).forEach((s, si) => {
      const sAgg = this.selAgg?.type === "subject" && this.selAgg.qi === qi && this.selAgg.si === si;
      html += `<div class="tree-node s${sAgg ? " sel" : ""}" data-qi="${qi}" data-si="${si}"><span class="caret">▼</span>${escapeHtml(s.name)}
        <span class="muted"> · ${(s.hits || []).length} HSP · len ${s.length}</span></div>`;
      (s.hits || []).forEach((h, hi) => {
        const sel = this.selHsp && this.selHsp.query_idx === qi
          && this.selHsp.subject_idx === si && this.selHsp.hit_idx === hi;
        html += `<div class="tree-node h${sel ? " sel" : ""}" data-qi="${qi}" data-si="${si}" data-hi="${hi}">
          <span class="swatch" style="background:${this.hspColor(h)}"></span>
          score ${this.fmtScore(h.score ?? h.bitscore)} · E ${h.evalue} · q:${h.qstart}-${h.qend} s:${h.sstart}-${h.send}
          <span class="muted"> · id ${h.identity}</span></div>`;
      });
    });
    html += "</div>";
    root.innerHTML = html;

    // 折叠/展开 + 点击节点聚合展示(guide 5.4.3:选中 Query/Subject 节点聚合其全部 HSP)
    // Collapse/expand + click a node for aggregated display (guide 5.4.3: selecting a Query/Subject node aggregates all its HSPs)
    root.querySelectorAll(".tree-node.q").forEach((n) => {
      n.querySelector(".caret").addEventListener("click", (e) => {
        e.stopPropagation();
        const body = root.querySelector(`.q-body[data-qi="${n.dataset.qi}"]`);
        const open = body.style.display !== "none";
        body.style.display = open ? "none" : "";
        n.querySelector(".caret").textContent = open ? "▶" : "▼";
      });
      n.addEventListener("click", () => {
        this.selAgg = { type: "query", qi: +n.dataset.qi };
        this.selHsp = null;
        this.renderTree();
        this.renderViz();
        this.renderDetail();
      });
    });
    root.querySelectorAll(".tree-node.s").forEach((n) => {
      n.querySelector(".caret").addEventListener("click", (e) => {
        e.stopPropagation();
        // 简单实现:折叠该 subject 下的所有 h
        // Simple implementation: collapse all h nodes under this subject
        const subs = root.querySelectorAll(`.tree-node.h[data-qi="${n.dataset.qi}"][data-si="${n.dataset.si}"]`);
        const open = subs[0]?.style.display !== "none";
        subs.forEach((h) => { h.style.display = open ? "none" : ""; });
        n.querySelector(".caret").textContent = open ? "▶" : "▼";
      });
      n.addEventListener("click", () => {
        this.selAgg = { type: "subject", qi: +n.dataset.qi, si: +n.dataset.si };
        this.selHsp = null;
        this.renderTree();
        this.renderViz();
        this.renderDetail();
      });
    });
    // HSP 点击 → 详情
    // HSP click -> details
    root.querySelectorAll(".tree-node.h").forEach((n) => {
      n.addEventListener("click", () => {
        // 清除全部行(q/s/h)的选中态——只清 .h.sel 会让基因名行背景残留
        // Clear the selection of all rows (q/s/h) — clearing only .h.sel would leave the gene-name row background stuck
        root.querySelectorAll(".tree-node.sel").forEach((x) => x.classList.remove("sel"));
        n.classList.add("sel");
        this.selAgg = null;
        this.selHsp = { query_idx: +n.dataset.qi, subject_idx: +n.dataset.si, hit_idx: +n.dataset.hi };
        this.renderViz();          // 对应条带高亮
        // Highlight the corresponding band
        this.renderDetail();
      });
    });
    if (scrollToHsp && this.selHsp) {
      const n = root.querySelector(`.tree-node.h[data-qi="${this.selHsp.query_idx}"][data-si="${this.selHsp.subject_idx}"][data-hi="${this.selHsp.hit_idx}"]`);
      if (n) n.scrollIntoView({ block: "nearest" });
    }
  },

  /* --- 详情 (guide 5.4.3:三行式带行号 + 高亮切换 + 聚合 + 复制) --- */
  // --- Details (guide 5.4.3: three-line layout with row numbers + highlight toggle + aggregation + copy) ---

  /* 当前展示的全部 HSP(聚合模式或单 HSP),供详情渲染与复制共用 */
  // All HSPs currently displayed (aggregation mode or single HSP), shared by detail rendering and copy
  currentHsps() {
    const parsed = this.result?.parsed || [];
    const agg = this.selAgg;
    if (agg?.type === "query" || agg?.type === "subject") {
      const q = parsed[agg.qi];
      if (!q) return [];
      const subs = agg.type === "subject" ? (q.subjects || []).slice(agg.si, agg.si + 1) : (q.subjects || []);
      const out = [];
      subs.forEach((s) => (s.hits || []).forEach((h) => out.push({ q, s, h })));
      return out;
    }
    if (!this.selHsp) return [];
    const q = parsed[this.selHsp.query_idx];
    const s = q?.subjects?.[this.selHsp.subject_idx];
    const h = s?.hits?.[this.selHsp.hit_idx];
    return h ? [{ q, s, h }] : [];
  },

  renderDetail() {
    const view = document.getElementById("detail-view");
    const parsed = this.result?.parsed || [];
    const agg = this.selAgg;
    let html = "";
    if (agg?.type === "query" || agg?.type === "subject") {
      const q = parsed[agg.qi];
      if (!q) { view.innerHTML = ""; return; }
      const subs = agg.type === "subject" ? (q.subjects || []).slice(agg.si, agg.si + 1) : (q.subjects || []);
      const n = subs.reduce((m, s) => m + (s.hits || []).length, 0);
      html += `<div class="hint" style="margin-bottom:8px">${
        agg.type === "query"
          ? t("blast.agg_query").replace("{q}", escapeHtml(q.name)).replace("{n}", n)
          : t("blast.agg_subject").replace("{q}", escapeHtml(q.name)).replace("{s}", escapeHtml(subs[0]?.name || "")).replace("{n}", n)
      }</div>`;
      subs.forEach((s) => (s.hits || []).forEach((h) => { html += this.hspBlock(q, s, h); }));
    } else {
      if (!this.selHsp || !parsed.length) {
        view.innerHTML = `<div class="empty">${t("blast.detail_empty")}</div>`;
        return;
      }
      const { query_idx, subject_idx, hit_idx } = this.selHsp;
      const q = parsed[query_idx];
      const s = q?.subjects?.[subject_idx];
      const h = s?.hits?.[hit_idx];
      if (!h) { view.innerHTML = ""; return; }
      html = this.hspBlock(q, s, h);
    }
    view.innerHTML = html;
  },

  /* 单个 HSP 块:元数据 + 三行式比对 */
  // Single HSP block: metadata + three-line alignment
  hspBlock(q, s, h) {
    const score = this.fmtScore(h.score ?? h.bitscore);
    return `
      <dl class="kv">
        <dt>Subject</dt><dd class="mono">${escapeHtml(s.name)} (len ${s.length})</dd>
        <dt>Score / E-value</dt><dd>${score} / ${h.evalue ?? ""}</dd>
        <dt>Identity</dt><dd>${h.identity ?? ""}${h.gaps ? " · Gaps " + h.gaps : ""}</dd>
        <dt>Query range</dt><dd>${h.qstart}-${h.qend} (${Math.abs((h.qend ?? 1) - (h.qstart ?? 1)) + 1} bp)</dd>
        <dt>Subject range</dt><dd>${h.sstart}-${h.send} ${h.strand ? "· " + h.strand : ""}${h.frame ? " · frame " + h.frame : ""}</dd>
      </dl>
      ${this.alnBlock(h)}`;
  },

  /* NCBI 三行式:Query/中间行/Sbjct,左右行号,等宽逐字符对齐
     minus 链 Sbjct 坐标递减(blast 原生行为,左号大右号小) */
  // NCBI three-line layout: Query / middle line / Sbjct, row numbers on both sides,
  // monospace per-character alignment. On the minus strand Sbjct coordinates decrease
  // (native BLAST behavior: left number larger than right)
  alnBlock(h) {
    const qseq = h.qseq || "", sseq = h.sseq || "", mid = h.mid || "";
    const n = Math.max(qseq.length, sseq.length, mid.length);
    if (!n) return "";
    const minus = (h.strand || "").includes("Minus");
    const qstart = Number(h.qstart) || 1;
    const sstart = Number(h.sstart) || 1;
    const send = Number(h.send) || 1;
    const matchMode = this.hlMode === "match";
    let out = "";
    for (let i = 0; i < n; i += 60) {
      const qc = qseq.slice(i, i + 60);
      const sc = sseq.slice(i, i + 60);
      const mc = mid.slice(i, i + 60);
      const L = Math.max(qc.length, sc.length);
      const qh = [], sh = [], mh = [];
      for (let j = 0; j < L; j++) {
        const a = qc[j] || " ", b = sc[j] || " ";
        const aUp = a.toUpperCase(), bUp = b.toUpperCase();
        const same = a !== "-" && b !== "-" && aUp === bUp;
        const isGap = a === "-" || b === "-";
        const isMismatch = !isGap && !same;
        let qcls = "", scls = "", mcls = "";
        if (matchMode) {
          if (same) { qcls = scls = "hl-id"; mcls = "hl-idm"; }
        } else if (isMismatch || isGap) { qcls = scls = "hl-mm"; }
        const qch = a === " " ? "&nbsp;" : escapeHtml(a);
        const sch = b === " " ? "&nbsp;" : escapeHtml(b);
        qh.push(qcls ? `<b class="${qcls}">${qch}</b>` : qch);
        sh.push(scls ? `<b class="${scls}">${sch}</b>` : sch);
        const mch = mc[j] || " ";
        mh.push(mcls ? `<b class="${mcls}">${mch === " " ? "&nbsp;" : escapeHtml(mch)}</b>` : (mch === " " ? "&nbsp;" : escapeHtml(mch)));
      }
      const qLeft = qstart + i, qRight = qstart + i + L - 1;
      const sLeft = minus ? send - i : sstart + i;
      const sRight = minus ? send - i - L + 1 : sstart + i + L - 1;
      out += `<div class="aln-line"><span class="aln-lab">Query</span><span class="aln-no">${qLeft}</span><span class="aln-seq">${qh.join("")}</span><span class="aln-no">${qRight}</span></div>`;
      out += `<div class="aln-line"><span class="aln-lab"></span><span class="aln-no"></span><span class="aln-seq aln-mid">${mh.join("")}</span><span class="aln-no"></span></div>`;
      out += `<div class="aln-line"><span class="aln-lab">Sbjct</span><span class="aln-no">${sLeft}</span><span class="aln-seq">${sh.join("")}</span><span class="aln-no">${sRight}</span></div>`;
    }
    return `<div class="aln-block">${out}</div>`;
  },

  /* 复制 Sbjct 序列(guide 5.4.3):纯净片段去缺口,反向链给基因组方向;
     多 HSP(聚合)时输出 FASTA,单 HSP 输出纯序列 */
  // Copy the Sbjct sequence (guide 5.4.3): clean fragment with gaps removed, reverse strand
  // given in genomic orientation; outputs FASTA for multiple HSPs (aggregation), plain
  // sequence for a single HSP
  copySbjct() {
    const hsps = this.currentHsps();
    if (!hsps.length) { toast(t("blast.detail_empty")); return; }
    const rc = (seq) => seq.split("").reverse().map((b) => {
      const m = { A: "T", T: "A", G: "C", C: "G", N: "N", U: "A", R: "R", Y: "Y", K: "K", M: "M", B: "B", D: "D", H: "H", V: "V" };
      const u = b.toUpperCase();
      return m[u] || b;
    }).join("");
    const parts = hsps.map(({ s, h }) => {
      const minus = (h.strand || "").includes("Minus");
      let seq = (h.sseq || "").replace(/-/g, "");
      if (minus) seq = rc(seq);
      if (hsps.length > 1) return `>${s.name} [${h.sstart}-${h.send}]\n${seq}`;
      return seq;
    });
    copyText(parts.join(hsps.length > 1 ? "\n" : " "));
  },

  /* 复制名称型查询(guide 5.4.3 旁新增,引物设计页侧翼模板输入格式):
     `>[gene],range=起-止,database=库basename`;仅单个 HSP 可复制(聚合模式多
     HSP 无单一坐标),远程比对无本地库不可用 */
  // Copy a name-style query (added next to guide 5.4.3; input format for flanking templates
  // on the primer-design page): `>[gene],range=start-end,database=db_basename`; only a single
  // HSP can be copied (aggregation mode has no single coordinate), and it is unavailable for
  // remote searches without a local database
  copyNameQuery() {
    const hsps = this.currentHsps();
    if (!hsps.length) { toast(t("blast.detail_empty")); return; }
    if (hsps.length !== 1) { toast(t("blast.copy_name_single")); return; }
    const opts = this.result?.options || {};
    const db = opts.db || "";
    if (opts.remote || !db) { toast(t("blast.copy_name_remote")); return; }
    const { s, h } = hsps[0];
    const name = (s.name || "").split(/[\s,]/)[0];
    if (!name) { toast(t("blast.copy_name_single")); return; }
    const s1 = Math.min(h.sstart, h.send);
    const s2 = Math.max(h.sstart, h.send);
    copyText(`>${name},range=${s1}-${s2},database=${db.split("/").pop()}`);
  },

  /* ---------------- 项目保存/加载(guide 9.4,顶栏全局按钮) ---------------- */
  // Project save/load (guide 9.4, global top-bar buttons)

  async saveProject() {
    const opts = this.result?.options || {
      query_text: document.getElementById("q-seq").value,
      db: document.getElementById("q-db").value,
      program: document.getElementById("q-program").value,
      evalue: document.getElementById("q-evalue").value,
      max_targets: document.getElementById("q-maxtargets").value,
      short_seq_mode: document.getElementById("q-short").checked,
      remote: document.getElementById("q-remote").checked,
      remote_db: document.getElementById("q-remote-db").value,
    };
    try {
      const payload = await api("/api/project/save", {
        method: "POST", json: true,
        body: {
          kind: "blast",
          raw_output: this.result?.raw_output || "",
          parsed: this.result?.parsed || [],
          db_prefix: document.getElementById("q-db").value,
          options: opts,
        },
      });
      const blob = new Blob([JSON.stringify(payload, null, 2)], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `blastprime_blast_${new Date().toISOString().slice(0, 10)}.json`;
      a.click();
      URL.revokeObjectURL(url);
      toast(t("blast.project_saved"));
    } catch (e) { toast(e.message); }
  },

  /* 加载项目(顶栏 load-project 经 app.js dispatchLoadedProject 分发):
     回填表单 + 恢复结果区;taskId=null(无原始任务,导出不可用) */
  // Load a project (top-bar load-project dispatched via app.js dispatchLoadedProject):
  // refill the form + restore the result area; taskId=null (no original task, export disabled)
  async applyLoadedProject(data) {
    if (!data || data.kind !== "blast") {
      toast(t("top.unknown_project"));
      return;
    }
    const opts = data.options || {};
    const set = (id, v) => { if (v != null) document.getElementById(id).value = v; };
    set("q-seq", opts.query_text);
    set("q-program", opts.program);
    set("q-evalue", opts.evalue);
    set("q-maxtargets", opts.max_targets);
    document.getElementById("q-short").checked = !!opts.short_seq_mode;
    document.getElementById("q-remote").checked = !!opts.remote;
    set("q-remote-db", opts.remote_db);
    if (opts.db) {
      const sel = document.getElementById("q-db");
      if ([...sel.options].some((o) => o.value === opts.db)) sel.value = opts.db;
    }
    document.getElementById("q-remote-db").style.display = opts.remote ? "" : "none";
    if (data.parsed && data.parsed.length) {
      this.result = { raw_output: data.raw_output || "", parsed: data.parsed, options: opts };
      this.taskId = null;
      this.rawShown = 0;
      this.selHsp = null;
      this.selAgg = null;
      this.vizQueryIdx = 0;
      document.getElementById("run-hint").textContent = t("blast.done") + " · " + t("blast.project_loaded_short");
      document.getElementById("result-area").style.display = "";
      requestAnimationFrame(() => {
        this.renderAll();
        this.selectFirstHsp();   // R12:默认选中第一个 HSP
        this.updateAnalysisPanel();
      });
      await this.saveResult();
    }
    this.saveState();
    toast(t("blast.project_loaded"));
  },

  /* --- 原始输出懒加载 --- */
  // --- Lazy loading of the raw output ---
  renderRaw() {
    const view = document.getElementById("raw-view");
    const raw = this.result?.raw_output || "";
    if (!raw) { view.innerHTML = `<div class="empty">${t("blast.raw_empty")}</div>`; return; }
    const lines = raw.split("\n");
    const chunk = lines.slice(0, this.rawShown + 1000).join("\n");
    const more = this.rawShown + 1000 < lines.length;
    view.innerHTML = `<pre class="code" style="max-height:none">${escapeHtml(chunk)}</pre>` +
      (more ? `<div class="row" style="margin-top:6px"><button class="small" id="btn-more-raw">${escapeHtml(t("common.load_more"))} (${lines.length - this.rawShown - 1000})</button></div>` : "");
    const btn = document.getElementById("btn-more-raw");
    if (btn) btn.addEventListener("click", () => { this.rawShown += 1000; this.renderRaw(); });
  },

  /* ---------------- 导出 ---------------- */
  // Export

  bindExport() {
    const exp = (path) => () => {
      if (!this.taskId) { toast(t("blast.export_no_task")); return; }
      downloadURL(`/api/export/blast/${this.taskId}/${path}`);
    };
    document.getElementById("btn-export-raw").addEventListener("click", exp("raw"));
    document.getElementById("btn-export-aln").addEventListener("click", exp("aln"));
    document.getElementById("btn-export-csv").addEventListener("click", exp("stats.csv"));
    document.getElementById("btn-export-fasta").addEventListener("click", exp("hits.fasta"));
  },

  /* ---------------- 导入 ---------------- */
  // Import

  bindImport() {
    document.getElementById("btn-import").addEventListener("click", () => {
      modalShow(document.getElementById("import-modal"));
    });
    // 文件导入(用户需求:复用 q-file 的追加导入模式,清空 fileInput 允许重复选同一文件)
    // File import (user request: reuse q-file's append-import mode; clear fileInput so the same file can be selected again)
    document.getElementById("import-file").addEventListener("change", async (e) => {
      const { added, failed } = await appendFilesToTextarea(e.target, document.getElementById("import-text"));
      if (added > 0) toast(t("blast.import_added").replace("{n}", added));
      if (failed > 0) toast(t("common.read_failed"));
    });
    document.querySelector("#import-modal .btn-cancel").addEventListener("click", () =>
      modalHide(document.getElementById("import-modal")));
    document.querySelector("#import-modal .btn-do-import").addEventListener("click", async () => {
      const content = document.getElementById("import-text").value;
      try {
        const res = await api("/api/blast/import", { method: "POST", json: true, body: { content } });
        this.result = res;
        this.taskId = null;
        this.rawShown = 0;
        this.selHsp = null;
        this.selAgg = null;
        this.vizQueryIdx = 0;
        // 与 fetchResult 一致:先显示再渲染,避免容器宽 0 导致画布拉伸
        // Same as fetchResult: show before rendering to avoid canvas stretching at container width 0
        document.getElementById("result-area").style.display = "";
        requestAnimationFrame(() => {
          this.renderAll();
          this.selectFirstHsp();   // R12:默认选中第一个 HSP
        });
        // 导入后短序列模式开启且查询为核酸时,仅显示分析标签,默认停在比对结果页
        // (用户需求:导入应显示比对结果;分析可手动切换)
        // After import, when short-sequence mode is on and the query is nucleotide, only the
        // analysis tab is shown while staying on the comparison-result page by default
        // (user request: import should show the comparison result; analysis can be switched manually)
        this.updateAnalysisPanel();
        modalHide(document.getElementById("import-modal"));
        toast(t("blast.import_ok"));
      } catch (e) { toast(e.message); }
    });
  },

  /* ---------------- 模块三:引物分析(短序列模式,guide 5.2/6) ----------------
     触发:勾选"短序列模式"且比对结果中所有查询 ≤100 bp(guide 6.1)。
     复用当前 BLAST 结果(查询序列 + HSP),不重新跑比对。 */
  // Module 3: primer analysis (short-sequence mode, guide 5.2/6):
  // Triggered when "short-sequence mode" is checked and all queries in the comparison result
  // are <= 100 bp (guide 6.1). Reuses the current BLAST result (query sequence + HSPs)
  // without rerunning the search.

  bindAnalysis() {
    document.querySelectorAll("#result-tabs .tab").forEach((b) => {
      b.addEventListener("click", () => this.switchTab(b.dataset.tab));
    });
    document.getElementById("ana-query").addEventListener("change", () => this.runAnalysis());
    document.getElementById("btn-copy-report").addEventListener("click", () => this.copyReport());
    // 勾选状态变化时联动标签可见性(用当前结果即时判断)
    // Update tab visibility when the checkbox state changes (judged against the current result)
    document.getElementById("q-short").addEventListener("change", () => {
      this.updateAnalysisPanel();
      this.updateGradeButtons();   // R12:分级模式随短序列模式切换生效键
    });
  },

  /* 触发条件(guide 6.1):短序列模式 + 有结果 + 全部查询长度 ≤100 bp 且为核酸 */
  // Trigger conditions (guide 6.1): short-sequence mode + result present + all query lengths <= 100 bp and nucleotide
  analysisAvailable() {
    if (!document.getElementById("q-short").checked) return false;
    const parsed = this.result?.parsed || [];
    if (!parsed.length) return false;
    if (!parsed.every((q) => (q.length || 0) > 0 && (q.length || 0) <= 100)) return false;
    // 比对结果带查询序列时校验核酸;导入场景无序列信息,交给后端校验
    // Validate nucleotide when the result carries query sequences; imported scenarios lack sequence info and defer to the backend
    const qs = this.result?.queries || [];
    if (qs.length) return qs.every((q) => /^[ACGTUN]+$/.test(q.seq || ""));
    return true;
  },

  switchTab(tab) {
    document.querySelectorAll("#result-tabs .tab").forEach((b) =>
      b.classList.toggle("active", b.dataset.tab === tab));
    document.getElementById("tab-results").style.display = tab === "results" ? "" : "none";
    document.getElementById("tab-analysis").style.display = tab === "analysis" ? "" : "none";
    if (tab === "analysis") this.runAnalysis();
  },

  /* 结果加载后刷新标签可见性并预计算分析(autoOpen: 兼容参数,已无调用方——导入也停在比对结果) */
  // After results load, refresh tab visibility and precompute the analysis (autoOpen: legacy
  // parameter with no remaining callers — import also stays on the comparison result)
  updateAnalysisPanel(autoOpen) {
    const tabs = document.getElementById("result-tabs");
    const ok = this.analysisAvailable();
    tabs.style.display = ok ? "" : "none";
    if (!ok) {
      document.getElementById("tab-results").style.display = "";
      document.getElementById("tab-analysis").style.display = "none";
      document.querySelectorAll("#result-tabs .tab").forEach((b) =>
        b.classList.toggle("active", b.dataset.tab === "results"));
      return;
    }
    // 填充查询下拉(多查询可切换)
    // Fill the query dropdown (switchable among multiple queries)
    const sel = document.getElementById("ana-query");
    const prev = sel.value;
    sel.innerHTML = (this.result?.parsed || []).map((q, i) =>
      `<option value="${i}">${escapeHtml(q.name)} (${q.length} bp)</option>`).join("");
    if (prev && sel.querySelector(`option[value="${prev}"]`)) sel.value = prev;
    this.runAnalysis();
    if (autoOpen) this.switchTab("analysis");
  },

  /* 复用当前比对结果分析选中查询:物理指标 + 特异性(不重新跑 BLAST) */
  // Analyze the selected query with the current comparison result: physical metrics + specificity (no BLAST rerun)
  async runAnalysis() {
    const idx = parseInt(document.getElementById("ana-query").value) || 0;
    const parsed = this.result?.parsed || [];
    const q = parsed[idx];
    if (!q) return;
    const token = ++this._anaToken;
    const cards = document.getElementById("ana-cards");
    cards.innerHTML = `<div class="muted">${t("blast.ana_run")}</div>`;
    const body = {
      query_name: q.name,
      db_prefix: document.getElementById("q-db").value,
      parsed: parsed.slice(idx, idx + 1),
      queries: this.result?.queries || [],
    };
    try {
      const report = await api("/api/analysis/primer", { method: "POST", json: true, body });
      if (token !== this._anaToken) return;   // 过期响应丢弃
      // Discard stale responses
      this._report = report;
      this.renderAnalysis(report);
    } catch (e) {
      if (token !== this._anaToken) return;
      this._report = null;
      cards.innerHTML = `<div class="tag err">${escapeHtml(e.message)}</div>`;
      document.getElementById("ana-verdict").innerHTML = "";
      document.getElementById("ana-hits").innerHTML = "";
      document.getElementById("ana-notes").textContent = "";
    }
  },

  /* 指标卡片 + 判定大字 + 特异性命中列表(可展开)+ 备注(guide 6.3) */
  // Metric cards + verdict headline + specificity hit list (expandable) + notes (guide 6.3)
  renderAnalysis(r) {
    const dot = (lv) => `<span class="ana-dot ${lv === "red" ? "red" : lv === "yellow" ? "yellow" : "green"}"></span>`;
    const items = r.items || [];
    const numeric = items.filter((it) => it.key !== "specificity");
    const spec = items.find((it) => it.key === "specificity");
    document.getElementById("ana-cards").innerHTML = numeric.map((it) => `
      <div class="ana-card">
        <div class="name">${dot(it.level)}${escapeHtml(this.anaTr(it.name))}</div>
        <div class="val">${escapeHtml(String(it.value ?? ""))}${it.unit ? `<span class="unit"> ${escapeHtml(it.unit)}</span>` : ""}</div>
        <div class="detail">${escapeHtml(this.anaTr(it.detail || ""))}</div>
      </div>`).join("");
    // 可用性判定大字
    // Verdict headline
    const vl = r.verdict?.level || "red";
    const vcol = vl === "red" ? "var(--err)" : vl === "yellow" ? "var(--warn)" : "var(--ok)";
    document.getElementById("ana-verdict").innerHTML = `
      <div class="ana-verdict" style="color:${vcol}">
        ${dot(vl)}${escapeHtml(this.anaTr(r.verdict?.label || ""))}
        <span class="sub">${escapeHtml(this.anaTr(r.verdict?.reason || ""))}</span>
      </div>`;
    // 特异性:目标 + 可展开的脱靶命中列表(3' 端结合情况)
    // Specificity: target + expandable off-target hit list (3' end binding)
    const hits = spec?.hits || [];
    const target = spec?.target;
    let hhtml = "";
    if (!spec) {
      hhtml = `<div class="empty">—</div>`;
    } else if (!target && !hits.length) {
      hhtml = `<div class="muted">${escapeHtml(this.anaTr(spec.detail || ""))}</div>`;
    } else {
      if (target) {
        hhtml += `<div class="row" style="margin-bottom:8px"><span class="tag ok">${t("blast.ana_target")}</span>
          <span class="mono">${escapeHtml(target.subject || "")}</span>
          <span class="muted">E ${target.evalue} · L ${target.length} · q ${target.qstart}-${target.qend}</span></div>`;
      }
      if (hits.length) {
        hhtml += hits.map((h) => {
          const hl = h.level || "green";
          const k = h.offset_3p || 0, errs = h.errors || 0;
          const thr = 8 + 2 * k + 3 * errs;   // 高危阈值(guide 6.2)
          // High-risk threshold (guide 6.2)
          const bound = hl === "red" ? t("blast.ana_high_risk")
            : hl === "yellow" ? t("blast.ana_potential") : t("blast.ana_safe");
          return `<div class="ana-hit">
            <span class="mono">${dot(hl)}${escapeHtml(h.subject || "")}</span>
            <span class="muted"> · E ${h.evalue} · L ${h.length} · ${t("blast.ana_off3p")} ${k} · ${t("blast.ana_errors")} ${errs}</span>
            <div class="ana-hit-body">${t("blast.ana_l_vs")}: L = ${h.length} vs ${t("blast.ana_threshold")} ${thr} (8 + 2×${k} + 3×${errs}) → ${bound}</div>
          </div>`;
        }).join("");
      }
    }
    document.getElementById("ana-hits").innerHTML = hhtml;
    document.getElementById("ana-hits").querySelectorAll(".ana-hit").forEach((n) =>
      n.addEventListener("click", () => n.classList.toggle("open")));
    // 备注:Tm 方法 / 二聚体 / 发夹
    // Notes: Tm method / dimer / hairpin
    document.getElementById("ana-notes").textContent =
      `${t("blast.ana_note_tm")}: ${this.anaTr(r.tm_method || "")} · ${t("blast.ana_note_dimer")}: ` +
      `${t("blast.ana_max_consec")} ${r.dimer?.max_consec ?? "-"} bp / ${t("blast.ana_max_total")} ` +
      `${r.dimer?.max_total ?? "-"} bp · ${t("blast.ana_note_hairpin")}: ${r.hairpin?.stem_len ?? "-"} bp`;
  },

  /* 后端固定中文文本 → 当前语言(引物分析;未命中的原样返回,数字不变) */
  // Backend fixed Chinese text -> current language (primer analysis; unmatched text is returned as-is, numbers unchanged)
  anaTr(text) {
    if (!text || typeof text !== "string") return text;
    // 带变量的特殊文本:仅命中 1 个靶标(chr1)
    // Special text with a variable: e.g. only 1 target hit (chr1)
    if (text.startsWith("仅命中 1 个靶标(")) {
      const subj = text.slice(text.indexOf("(") + 1, text.lastIndexOf(")"));
      return t("blast.ana_spec_good_d").replace("{subject}", subj);
    }
    const M = {
      "GC 含量": "blast.ana_item_gc", "Tm 值": "blast.ana_item_tm",
      "3' 端 G/C 数(GC 夹子,末 5 bp)": "blast.ana_item_clamp",
      "自互补二聚体(最大连续配对)": "blast.ana_item_dimer",
      "BLAST 特异性": "blast.ana_item_spec",
      "不可用": "blast.ana_v_red", "风险/待定": "blast.ana_v_yellow", "可用": "blast.ana_v_green",
      "存在红色(高风险)指标": "blast.ana_r_red", "无红色但有黄色指标": "blast.ana_r_yellow",
      "全部指标合格": "blast.ana_r_green",
      "库中未找到匹配靶标": "blast.ana_spec_nofound",
      "当前比对结果中无该查询的命中": "blast.ana_spec_nofound_d",
      "特异性高": "blast.ana_spec_good",
      "高非特异性风险": "blast.ana_spec_highrisk",
      "存在 3' 端可能结合的脱靶命中(高危)": "blast.ana_spec_highrisk_d",
      "潜在非特异风险": "blast.ana_spec_potential",
      "存在潜在非特异命中": "blast.ana_spec_potential_d",
      "匹配单个基因但存在多处结合区": "blast.ana_spec_multihit",
      "命中同一库条目内多处(同一目标内多 HSP)": "blast.ana_spec_multihit_d",
      "primer3 盐校正模型": "blast.ana_tm_method",
      "达标": "blast.ana_ok_word", "违规": "blast.ana_bad_word", "待定": "blast.ana_pending_word",
    };
    return M[text] ? t(M[text]) : text;
  },

  /* 复制完整报告(纯文本,guide 6.3) */
  // Copy the full report (plain text, guide 6.3)
  copyReport() {
    const r = this._report;
    if (!r) { toast(t("blast.analysis_empty")); return; }
    const lv = (l) => l === "green" ? t("blast.ana_lv_good") : l === "yellow" ? t("blast.ana_lv_ok") : t("blast.ana_lv_bad");
    const lines = [t("blast.ana_report_title"), "──────────────────",
      `${t("blast.analysis_seq")}: ${r.seq} (${r.seq.length} bp)`];
    (r.items || []).forEach((it) => {
      const val = typeof it.value === "string" ? this.anaTr(it.value) : it.value;
      lines.push(`${this.anaTr(it.name)}: ${val ?? ""}${it.unit || ""} — ${lv(it.level)}`);
    });
    lines.push("──────────────────",
      `${t("blast.ana_verdict")}: ${this.anaTr(r.verdict?.label || "")} (${lv(r.verdict?.level || "red")}) — ${this.anaTr(r.verdict?.reason || "")}`);
    (r.reasons || []).forEach(([st, txt]) => lines.push(`  - ${this.anaTr(st)}: ${this.anaTr(txt)}`));
    copyText(lines.join("\n"));
  },

  /* 语言切换后重绘动态渲染的本地化内容(静态 data-i18n 由 app.js 处理) */
  // Redraw dynamically rendered localized content after a language switch (static data-i18n is handled by app.js)
  onLangChanged() {
    this.updateGradeButtons();   // R12:分级按钮/提示行随语言重绘
    this._renderRunStatus();
    if (!this.result) return;
    this.renderTree();
    this.renderDetail();
    if (this._report) this.renderAnalysis(this._report);
    this.renderRaw();
  },

  /* 主题切换后重绘 canvas 绘制内容(基因名底框 --bg-soft、名称文字 --fg、
     坐标轴标尺等随主题变量变化;结果树/详情为 DOM 自动跟随) */
  // Redraw canvas content after a theme switch (gene-name backing box --bg-soft, name text
  // --fg, axis ruler, etc. follow theme variables; the result tree/details are DOM and follow automatically)
  onThemeChanged() {
    if (!this.result || !document.getElementById("viz-canvas")?.clientWidth) return;
    this.renderViz();
  },
};

