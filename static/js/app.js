/* BlastPrime Studio — 共享骨架:任务/SSE/状态栏/日志抽屉/设置/主题/i18n */
// Shared skeleton: tasks / SSE / status bar / log drawer / settings / theme / i18n
"use strict";

const App = {
  settings: null,
  env: null,
  tasks: new Map(),           // id -> snapshot
  activeTask: null,           // 当前状态栏展示的任务 id
  // id of the task currently shown in the status bar
  logs: [],                   // [{task_id, t, level, msg}]
  sse: null,
  pageInit: null,

  /* ================================ 启动 ================================ */
  // Startup

  async init() {
    // 全站文本框禁用拼写检查(序列/日志/备注等不应出现红波浪线)
    // Disable spell-check on all textboxes site-wide (sequences/logs/notes must not show red underlines)
    document.querySelectorAll("textarea, input[type='text'], input[type='search'], input:not([type])")
      .forEach((el) => { el.spellcheck = false; });
    await this.loadSettings();
    // 顶栏临时切换(主题/语言)经 localStorage 恢复:刷新/切页/重开浏览器均不丢,
    // 优先于 config.json 持久化设置(设置对话框保存后清除,见 saveSettings)
    // Restore the top-bar temporary toggles (theme/language) from localStorage:
    // survives refresh/page switches/browser restarts, takes precedence
    // over the persisted config.json settings (cleared on save from the settings dialog)
    const tmpUI = loadPageState("uitmp") || {};
    applyTheme(tmpUI.theme || this.settings.theme || "system");
    setLang(tmpUI.lang || this.settings.lang || "zh");
    this.renderNav();
    this.bindGlobal();
    this.connectSSE();
    await this.checkEnv();
    if (window.PAGE) {
      try { this.pageInit = await window.PAGE.init(); }
      catch (e) { console.error("页面初始化失败", e); }
    }
    this.renderStatus();
  },

  renderNav() {
    let cur = (location.pathname || "/").replace(/\/+$/, "") || "/";
    if (cur === "/index.html") cur = "/";
    document.querySelectorAll("nav.mainnav a").forEach((a) => {
      const key = a.dataset.nav;
      if (key) a.textContent = t("nav." + key);
      const href = (a.getAttribute("href") || "/").replace(/\/+$/, "") || "/";
      if (cur === href || (href !== "/" && cur.endsWith(href))) a.classList.add("active");
    });
    // 应用名
    // App name
    const brand = document.querySelector(".brand");
    if (brand && brand.dataset.i18n) brand.textContent = t("brand.name");
  },

  /* ================================ 设置 ================================ */
  // Settings

  async loadSettings() {
    try {
      this.settings = await api("/api/settings");
    } catch (e) {
      this.settings = { lang: "zh", theme: "system" };
    }
  },

  /* ================================ 环境 ================================ */
  // Environment

  async checkEnv() {
    try {
      this.env = await api("/api/env");
    } catch (e) {
      this.env = { missing: ["blastn"] };
    }
    const ok = !(this.env.missing && this.env.missing.length);
    const light = document.querySelector(".env-light");
    const label = document.querySelector(".env-label");
    if (light) light.classList.toggle("bad", !ok);
    if (label) label.textContent = t(ok ? "env.ok" : "env.bad");
    const banner = document.getElementById("env-banner");
    if (banner) banner.classList.toggle("show", !ok);
    document.body.classList.toggle("env-bad", !ok);
    document.body.classList.toggle("env-ok", ok);
    if (window.onEnvChanged) window.onEnvChanged(ok);
  },

  /* ================================ 任务 / SSE ================================ */
  // Tasks / SSE

  connectSSE() {
    const es = new EventSource("/api/tasks/events");
    this.sse = es;
    es.onmessage = (ev) => {
      let msg;
      try { msg = JSON.parse(ev.data); } catch (e) { return; }
      this.handleSSE(msg);
    };
    es.onerror = () => {
      // EventSource 自动重连;浏览器关闭连接时置灰状态
      // EventSource auto-reconnects; grey out the status dot while the browser keeps the connection closed
      document.querySelector(".status-dot")?.classList.remove("running");
    };
    es.onopen = () => { /* 重连后服务端会再发 snapshot */ };
    // The server re-sends a snapshot after reconnect
  },

  handleSSE(msg) {
    const type = msg.type || "";
    if (type === "snapshot") {
      this.tasks = new Map();
      this.logs = [];
      (msg.tasks || []).forEach((s) => {
        this.upsertTask(s);
        // 从快照回填历史日志,否则刷新后日志抽屉为空
        // Backfill historical logs from the snapshot, otherwise the log drawer would be empty after a refresh
        (s.logs || []).forEach((l) =>
          this.logs.push({ task_id: s.id, t: l.t, level: l.level || "info", msg: l.msg }));
      });
      // 快照中已成功的设计/比对任务 → 结果兜底持久化(任务完成时页面已
      // 切走,页面级 saveResult 从未执行;persistTaskResult 内部按 taskId
      // 去重,已保存的跳过)
      // Succeeded design/blast tasks in the snapshot -> fallback result
      // persistence (the page may have navigated away before its saveResult
      // ran; persistTaskResult dedupes by taskId and skips saved ones)
      (msg.tasks || []).forEach((s) => {
        if (s.status === "succeeded" && (s.kind === "design" || s.kind === "blast")) {
          this.persistTaskResult(s.id);
        }
      });
      this.renderLogs();
      this.renderStatus();
      return;
    }
    if (type === "ping") return;
    const id = msg.task_id;
    const payload = msg.payload || {};
    if (type === "task_started") { this.upsertTask(payload); }
    else if (type === "task_progress" && id) {
      const s = this.tasks.get(id);
      if (s) { s.progress = payload.progress; s.progress_label = payload.label || ""; }
    } else if (type === "task_log" && id) {
      this.pushLog(id, payload);
    } else if (type === "task_succeeded" && id) {
      // 任务完成即持久化结果:页面存活时其 saveResult 后写覆盖(带交互状态),
      // 页面已切走时此处兜底保证切回来结果仍在、k-mer 缓存已入 sessionStorage
      // Persist the result on completion: while the page lives its own
      // saveResult overwrites later (with interactive state); if the page is
      // gone this fallback keeps the result (and k-mer caches) available
      this.upsertTask(payload);
      this.persistTaskResult(id);
    } else if (type.startsWith("task_") && id) {
      this.upsertTask(payload);
    }
    this.renderStatus();
  },

  /* 任务结果兜底持久化:design/blast 任务成功时把结果写入 localStorage
     (页面级 saveResult 之外的保险 —— 任务在后台完成而用户已切页时,页面
     run() 的保存代码永远不会执行,结果与 k-mer 缓存就此丢失)。与页面级
     保存同结构,页面已保存(同 taskId)时跳过,避免覆盖交互状态。 */
  // Task-result fallback persistence: on design/blast success, write the
  // result to localStorage (a safety net for page-level saveResult — if the
  // user switched pages while the task ran in the background, the page's
  // save code never runs and the result plus k-mer caches are lost). Same
  // envelope as the page-level save; skipped when that taskId is already
  // saved so interactive state (resultIdx/selPair) is never clobbered.
  async persistTaskResult(taskId) {
    const s = this.tasks.get(taskId);
    if (!s || s.status !== "succeeded") return;
    const kind = s.kind;
    if (kind !== "design" && kind !== "blast") return;
    const key = kind === "design" ? "bp_design_result" : "bp_blast_result";
    try {
      const prev = await loadResult(key);
      if (prev && prev._v === 1 && prev.taskId === taskId) return;   // 页面已保存
      // Already saved by the page
    } catch (e) { /* 损坏数据:允许覆盖 */ }
    try {
      const { status, result } = await api(`/api/tasks/${taskId}/result`);
      if (status !== "succeeded" || !result) return;
      if (kind === "design") {
        // 瘦身后再持久化:大命中结果的 step1_hsps 可达数十 MB,不剔除则
        // 写不进 5 MB 配额,切页回来结果仍是空的(localStorage 超限还有
        // IndexedDB 兜底,见 persistResult)
        // Slim before persisting: step1_hsps on heavy hits reaches tens of MB;
        // without stripping it the 5 MB quota rejects the write and the
        // result stays empty after a page switch (persistResult additionally
        // falls back to IndexedDB when localStorage overflows)
        await persistResult(key, {
          _v: 1, taskId, result: slimDesignResult(result), resultIdx: 0, selPair: null,
        });
        // k-mer 分析一并兜底入 sessionStorage(原 storeKmerCaches 依赖页面
        // 存活;任务完成时页面已切走 → 缓存从未写入 → 再设计"不复用")
        // K-mer analyses also go to sessionStorage (the original
        // storeKmerCaches ran only while the page lived; a page switch
        // mid-task meant the caches were never written -> "not reused")
        (result.results || []).forEach((q) => {
          const c = q && q.kmer_cache;
          if (c && c.key) {
            try { sessionStorage.setItem("bp_kmer_cache:" + c.key, JSON.stringify(c)); }
            catch (e2) { /* quota: 放弃该条缓存 */ }
          }
        });
      } else {
        await persistResult(key, {
          _v: 1, taskId, result, rawShown: 0,
        });
      }
      // 通知页面:兜底写入完成(切页回来时页面 init 的同步 restoreResult 在
      // 写入前已执行 → 空结果;事件让页面此时再恢复渲染)
      // Notify the page that the fallback write finished (a page restored
      // after navigation ran its synchronous restoreResult before the async
      // write landed; the event lets it recover the result now)
      window.dispatchEvent(new CustomEvent("bp:result-persisted", {
        detail: { taskId, kind },
      }));
    } catch (e) { /* 拉取失败/超限额:静默,页面级保存兜底仍在 */ }
    // Fetch failure / quota exceeded: silent; the page-level save is still there
  },

  upsertTask(snapshot) {
    // 事件快照(with_logs=False)不带 logs 字段;整体替换会丢已累积的实时日志,需合并保留
    // Event snapshots (with_logs=False) carry no logs field; wholesale replacement would drop accumulated live logs, so merge and preserve
    const prev = this.tasks.get(snapshot.id);
    if (prev && prev.logs && snapshot.logs === undefined) {
      snapshot = { ...snapshot, logs: prev.logs };
    }
    this.tasks.set(snapshot.id, snapshot);
    if (snapshot.status === "running" || snapshot.status === "pending") {
      if (!this.activeTask) this.activeTask = snapshot.id;
    } else if (this.activeTask === snapshot.id) {
      this.activeTask = this.firstActive();
      const statusText = snapshot.status === "succeeded" ? t("status.done")
        : snapshot.status === "failed" ? t("status.failed")
        : snapshot.status === "cancelled" ? t("status.cancelled") : "";
      if (statusText) toast(`${this.taskLabel(snapshot)} — ${statusText}`);
    }
  },

  firstActive() {
    return [...this.tasks.values()].find(
      (t) => t.status === "running" || t.status === "pending")?.id || null;
  },

  taskLabel(s) {
    return s.title || s.kind || s.id;
  },

  pushLog(taskId, entry) {
    const line = { task_id: taskId, t: entry.t, level: entry.level || "info", msg: entry.msg };
    this.logs.push(line);
    if (this.logs.length > 3000) this.logs.splice(0, this.logs.length - 3000);
    // 同步写回任务快照,页面级轮询(如构建日志框)才能读到实时日志
    // Write back into the task snapshot in sync, so page-level polling (e.g. the build log box) can read live logs
    const s = this.tasks.get(taskId);
    if (s) {
      s.logs = s.logs || [];
      s.logs.push(entry);
    }
    // 日志批渲染:每帧只 append 一次 + 滚动一次,日志洪峰时不逐行触发布局
    // Batched log rendering: only one append + one scroll per frame, avoiding per-line layout triggers during log floods
    if (!this._logPending) this._logPending = [];
    this._logPending.push(line);
    if (!this._logFrame) {
      this._logFrame = requestAnimationFrame(() => {
        this._logFrame = null;
        const pending = this._logPending || [];
        this._logPending = null;
        const list = document.getElementById("log-list");
        if (!list) return;
        pending.forEach((l) => this.renderLogLine(l, false));
        list.scrollTop = list.scrollHeight;
      });
    }
  },

  /* ================================ 状态栏 ================================ */
  // Status bar

  renderStatus() {
    const bar = document.getElementById("status-bar");
    if (!bar) return;
    const dot = bar.querySelector(".status-dot");
    const text = bar.querySelector(".status-text");
    const prog = bar.querySelector(".progress-wrap");
    const pbar = bar.querySelector(".progress-bar");
    const cancelBtn = bar.querySelector("#btn-cancel-task");
    const active = this.activeTask ? this.tasks.get(this.activeTask) : null;
    if (active && (active.status === "running" || active.status === "pending")) {
      dot.className = "status-dot running";
      // 百分比直接拼进状态文字:进度条只有 6px 高、低进度时宽度很短,
      // 深色主题下轨道与底栏同色几乎隐形,文字里的百分比是最醒目的进度指示
      // Percentage goes straight into the status text: the bar is only 6px tall
      // and short at low progress, and its track blends into the dark footer,
      // so the text percentage is the most visible progress indicator
      const pct = active.progress > 0 ? ` ${Math.round(active.progress)}%` : "";
      text.textContent = `${this.taskLabel(active)}${pct}${active.progress_label ? " — " + active.progress_label : ""}`;
      prog.style.display = "";
      // 有真实百分比→定量进度条;否则(建库/BLAST 等无进度的长任务)→不确定进度动画
      // Real percentage -> determinate progress bar; otherwise (long tasks without progress such as makeblastdb/BLAST) -> indeterminate progress animation
      if (active.progress > 0) { pbar.classList.remove("indet"); pbar.style.width = active.progress + "%"; }
      else { pbar.classList.add("indet"); pbar.style.width = ""; }
      cancelBtn.style.display = "";
    } else {
      dot.className = "status-dot ok";
      text.textContent = t("status.idle");
      prog.style.display = "none";
      cancelBtn.style.display = "none";
    }
  },

  /* ================================ 日志抽屉 ================================ */
  // Log drawer

  async openLogs() {
    document.getElementById("log-drawer")?.classList.add("open");
    // 打开时拉取服务端最新任务快照(含完整日志),确保抽屉始终有历史可看
    // Pull the latest task snapshots from the server on open (including full logs) so the drawer always has history to show
    try {
      const data = await api("/api/tasks");
      if (!data.tasks) return;
      this.tasks = new Map(data.tasks.map((s) => [s.id, s]));
      this.activeTask = this.firstActive();
      this.logs = [];
      data.tasks.forEach((s) =>
        (s.logs || []).forEach((l) =>
          this.logs.push({ task_id: s.id, t: l.t, level: l.level || "info", msg: l.msg })));
      this.renderLogs();
      this.renderStatus();
    } catch (e) { /* 拉取失败时保留现有日志 */ }
    // Keep existing logs if the fetch fails
  },
  closeLogs() { document.getElementById("log-drawer")?.classList.remove("open"); },

  /* 下载全部任务日志为文本文件(含任务 id 与时间戳,按时间升序)。
     下载完整历史(不受抽屉的级别过滤/搜索影响)。 */
  // Download all task logs as a text file (task id + timestamp, ascending).
  // Downloads the full history (unaffected by the drawer's level filter/search).
  downloadLogs() {
    const lines = [...this.logs]
      .sort((a, b) => a.t - b.t)
      .map((l) => {
        const ts = new Date(l.t * 1000);
        const pad = (v) => String(v).padStart(2, "0");
        const time = `${ts.getFullYear()}-${pad(ts.getMonth() + 1)}-${pad(ts.getDate())} ` +
          `${pad(ts.getHours())}:${pad(ts.getMinutes())}:${pad(ts.getSeconds())}`;
        return `${time} [${(l.level || "info").toUpperCase()}] [${l.task_id || ""}] ${l.msg}`;
      });
    if (!lines.length) { toast(t("log.empty")); return; }
    const stamp = new Date();
    const pad = (v) => String(v).padStart(2, "0");
    const name = `blastprime_logs_${stamp.getFullYear()}${pad(stamp.getMonth() + 1)}${pad(stamp.getDate())}_` +
      `${pad(stamp.getHours())}${pad(stamp.getMinutes())}${pad(stamp.getSeconds())}.txt`;
    const blob = new Blob([lines.join("\n") + "\n"], { type: "text/plain;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url; a.download = name;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  },

  renderLogs() {
    const list = document.getElementById("log-list");
    if (!list) return;
    list.innerHTML = "";
    // 全量渲染按时间升序(旧→新):任务快照按 created_at 降序回填,直接渲染会让新任务日志排在上方
    // Full render in ascending time order (old -> new): task snapshots are backfilled in descending created_at order, so rendering directly would place new tasks' logs on top
    [...this.logs].sort((a, b) => a.t - b.t).forEach((l) => this.renderLogLine(l, false));
    list.scrollTop = list.scrollHeight;
  },

  logVisible(l) {
    const lv = document.getElementById("log-level-filter")?.value || "all";
    const q = (document.getElementById("log-search")?.value || "").trim().toLowerCase();
    if (lv !== "all" && l.level !== lv) return false;
    if (q && !l.msg.toLowerCase().includes(q)) return false;
    return true;
  },

  renderLogLine(l, append = true) {
    const list = document.getElementById("log-list");
    if (!list) return;
    if (!this.logVisible(l)) return;
    const div = document.createElement("div");
    div.className = "log-line " + (l.level || "info");
    const tm = new Date(l.t * 1000);
    const ts = `${String(tm.getHours()).padStart(2, "0")}:${String(tm.getMinutes()).padStart(2, "0")}:${String(tm.getSeconds()).padStart(2, "0")}`;
    div.innerHTML = `<span class="lt">${ts}</span>${escapeHtml(l.msg)}`;
    if (append) {
      list.appendChild(div);
      list.scrollTop = list.scrollHeight;
    } else {
      list.appendChild(div);
    }
  },

  /* ================================ 设置对话框 ================================ */
  // Settings dialog

  openSettings() {
    const s = this.settings || {};
    const m = document.getElementById("settings-modal");
    if (!m) return;
    m.querySelector("#set-lang").value = s.lang || "zh";
    m.querySelector("#set-theme").value = s.theme || "system";
    m.querySelector("#set-loglevel").value = s.loglevel || "INFO";
    m.querySelector("#set-logfile").value = s.logfile || "";
    m.querySelector("#set-blastdir").value = s.blast_bin_dir || "";
    m.classList.add("show");
  },

  closeSettings() {
    document.getElementById("settings-modal")?.classList.remove("show");
  },

  async saveSettings() {
    const m = document.getElementById("settings-modal");
    const body = {
      lang: m.querySelector("#set-lang").value,
      theme: m.querySelector("#set-theme").value,
      loglevel: m.querySelector("#set-loglevel").value,
      logfile: m.querySelector("#set-logfile").value,
      blast_bin_dir: m.querySelector("#set-blastdir").value,
    };
    try {
      this.settings = await api("/api/settings", { method: "POST", json: true, body });
      setLang(this.settings.lang);
      applyTheme(this.settings.theme);
      // 设置对话框保存后,顶栏临时切换按钮同步为持久化主题,并清除 localStorage 临时状态(回到跟随系统)
      // After saving from the settings dialog, sync the temporary top-bar toggle to the persisted theme and clear the localStorage temporary state (back to following the system)
      const tb = document.getElementById("theme-switch");
      if (tb) {
        renderThemeSwitch(this.settings.theme || "system");
        delete tb.dataset.manual;
      }
      try { localStorage.removeItem("bp_uitmp"); } catch (e) {}
      this.checkEnv();
      this.closeSettings();
      toast(t("settings.saved"));
    } catch (e) {
      toast(e.message || t("common.save_failed"));
    }
  },

  /* 清空 localStorage:页面参数/临时主题语言/比对与设计结果全部清除,回退到 config.json 持久化设置。
     先提示后重载(延时让 toast 可见),页面以初始状态重新启动 */
  // Clear localStorage: page params, temporary theme/language, blast and design results are all
  // wiped; falls back to the settings persisted in config.json. Toast first, then reload (with a
  // delay so the toast is visible); the page restarts in its initial state
  clearLocalStorage() {
    try { localStorage.clear(); } catch (e) {}
    applyTheme(this.settings?.theme || "system");
    setLang(this.settings?.lang || "zh");
    this.closeSettings();
    toast(t("settings.storage_cleared"));
    setTimeout(() => location.reload(), 800);
  },

  /* ================================ 全局绑定 ================================ */
  // Global bindings

  bindGlobal() {
    document.querySelectorAll("[data-i18n]").forEach((el) => {
      if (!el.dataset.i18nSet) {
        el.dataset.i18nSet = "1";
        // 只替换首个非空文本节点,保留子元素(如标题内嵌的状态徽章 span)
        // Replace only the first non-empty text node, keeping child elements (e.g. status badge spans embedded in titles)
        let done = false;
        for (const n of el.childNodes) {
          if (n.nodeType === Node.TEXT_NODE) { n.textContent = t(el.dataset.i18n); done = true; break; }
        }
        if (!done) el.textContent = t(el.dataset.i18n);
      }
    });
    document.querySelectorAll("[data-i18n-ph]").forEach((el) => {
      el.placeholder = t(el.dataset.i18nPh);
    });

    // 导航 & 顶层按钮
    // Navigation & top-level buttons
    document.querySelectorAll("[data-action]").forEach((el) => {
      el.addEventListener("click", () => {
        const act = el.dataset.action;
        if (act === "settings") this.openSettings();
        else if (act === "logs") this.openLogs();
        else if (act === "cancel-task") this.cancelActive();
        else if (act === "env-guide") this.showGuide();
        else if (act === "env-setdir") this.openSettings();
        else if (act === "save-project") this.saveProject();
        else if (act === "load-project") this.chooseProjectFile();
      });
    });

    // 语言切换菜单:选项从 LANGUAGES 单一来源填充(顶栏 + 设置对话框),新增语言只需加一行
    // Language switch menus: options filled from the LANGUAGES single source (top bar + settings dialog); adding a language only needs one line
    ["lang-switch", "set-lang"].forEach((id) => {
      const sel = document.getElementById(id);
      if (!sel) return;
      sel.innerHTML = "";
      LANGUAGES.forEach((l) => {
        const opt = document.createElement("option");
        opt.value = l.code;
        opt.textContent = l.label;
        sel.appendChild(opt);
      });
    });
    const langSel = document.getElementById("lang-switch");
    if (langSel) {
      // 临时语言优先于 config.json 持久化值(经 localStorage 恢复)
      // The temporary language takes precedence over the persisted config.json value (restored from localStorage)
      langSel.value = loadPageState("uitmp")?.lang || this.settings?.lang || "zh";
      langSel.title = t("top.lang_tmp_title");
      langSel.addEventListener("change", () => {
        this.settings.lang = langSel.value;
        setLang(langSel.value);
        // 临时切换写入 localStorage(刷新/切页/重开浏览器均保留)
        // Persist the temporary toggle to localStorage (survives refresh/page switches/browser restarts)
        savePageState("uitmp", { ...(loadPageState("uitmp") || {}), lang: langSel.value });
      });
      window.onLangChanged = () => {
        langSel.value = I18N_LANG;
        langSel.title = t("top.lang_tmp_title");
        document.querySelectorAll("[data-i18n]").forEach((el) => {
          if (el.dataset.i18n && !el.dataset.i18nStatic) el.textContent = t(el.dataset.i18n);
        });
        document.querySelectorAll("[data-i18n-ph]").forEach((el) => {
          el.placeholder = t(el.dataset.i18nPh);
        });
        // 顶栏主题按钮:文本随语言刷新,保持当前临时选择
        // Top-bar theme button: text refreshes with the language while keeping the current temporary choice
        const tb = document.getElementById("theme-switch");
        if (tb) renderThemeSwitch(tb.dataset.theme || this.settings?.theme || "system");
        // 动态渲染的内容(如数据库记录表)由各页面自行重绘
        // Dynamically rendered content (e.g. the database records table) is redrawn by each page itself
        window.PAGE?.onLangChanged?.();
      };
    }

    // 顶栏主题临时切换:点按在 浅色↔深色 间循环,仅本次会话生效,不写入 config.json(设置对话框保存才持久化)
    // Temporary top-bar theme toggle: clicking cycles light<->dark, effective for this session only, not written to config.json (persisted only when saved from the settings dialog)
    const themeBtn = document.getElementById("theme-switch");
    if (themeBtn) {
      // 临时主题优先于 config.json 持久化值(经 localStorage 恢复);恢复过临时主题即视为手动切换过
      // The temporary theme takes precedence over the persisted config.json value (restored from localStorage); restoring a temporary theme counts as a manual toggle
      const tmpUI = loadPageState("uitmp") || {};
      renderThemeSwitch(tmpUI.theme || this.settings?.theme || "system");
      if (tmpUI.theme) themeBtn.dataset.manual = "1";
      themeBtn.addEventListener("click", () => {
        const cur = themeBtn.dataset.theme || "light";
        const next = cur === "light" ? "dark" : "light";
        themeBtn.dataset.manual = "1";  // 手动切换过,系统主题变化不再自动跟随
        // Manually toggled, so system theme changes no longer auto-follow
        applyTheme(next);
        renderThemeSwitch(next);
        // 临时切换写入 localStorage(刷新/切页/重开浏览器均保留;设置对话框保存后清除,见 saveSettings)
        // Persist the temporary toggle to localStorage (survives refresh/page switches/browser restarts; cleared after the settings dialog saves, see saveSettings)
        savePageState("uitmp", { ...(loadPageState("uitmp") || {}), theme: next, manual: "1" });
      });
    }

    // 日志抽屉
    // Log drawer
    document.getElementById("log-level-filter")?.addEventListener("change", () => this.renderLogs());
    document.getElementById("log-search")?.addEventListener("input", () => this.renderLogs());
    document.querySelector("#log-drawer .btn-close")?.addEventListener("click", () => this.closeLogs());
    document.getElementById("btn-log-download")?.addEventListener("click", () => this.downloadLogs());

    // 设置对话框
    // Settings dialog
    document.getElementById("settings-modal")?.addEventListener("click", (e) => {
      if (e.target.id === "settings-modal") this.closeSettings();
    });
    document.querySelector("#settings-modal .btn-save")?.addEventListener("click", () => this.saveSettings());
    document.querySelector("#settings-modal .btn-cancel")?.addEventListener("click", () => this.closeSettings());
    // 清空 localStorage(三页同构,共享 handler)
    // Clear localStorage (same markup on all three pages, shared handler)
    document.querySelectorAll("#settings-modal .btn-clear-storage").forEach((btn) => {
      btn.addEventListener("click", () => this.clearLocalStorage());
    });
  },

  async cancelActive() {
    if (!this.activeTask) return;
    try { await api(`/api/tasks/${this.activeTask}/cancel`, { method: "POST" }); }
    catch (e) { toast(e.message || t("common.cancel_failed")); }
  },

  showGuide() {
    window.open("/static/guide.html", "_blank");
  },

  /* ================= 全局项目保存/加载(guide.md 9.4,保存 blast 与引物设计页全部输入输出) ================= */
  // Global project save/load (guide.md 9.4, saves all inputs/outputs of the blast and primer-design pages)

  saveProject() {
    if (window.PAGE?.saveProject) window.PAGE.saveProject();
    else toast(t("top.no_project"));
  },

  chooseProjectFile() {
    // 共享隐藏文件选择框:懒创建一次,复用读文件流程
    // Shared hidden file picker: lazily created once, reusing the read-file flow
    if (!this._projectInput) {
      this._projectInput = document.createElement("input");
      this._projectInput.type = "file";
      this._projectInput.accept = ".json,application/json";
      this._projectInput.style.display = "none";
      document.body.appendChild(this._projectInput);
      this._projectInput.addEventListener("change", () => {
        const f = this._projectInput.files && this._projectInput.files[0];
        this._projectInput.value = "";   // 允许重复选择同一文件
        // Allow re-selecting the same file
        if (!f) return;
        const rd = new FileReader();
        rd.onload = async () => {
          try {
            const data = await api("/api/project/load", { method: "POST", json: true, body: { content: String(rd.result) } });
            this.dispatchLoadedProject(data);
          } catch (e) {
            toast(e.message || t("common.load_failed"));
          }
        };
        rd.readAsText(f);
      });
    }
    this._projectInput.click();
  },

  /* 按 kind 分发加载到的项目:同页直接调页面 applyLoadedProject,跨页经 localStorage 通道跳转 */
  // Dispatch the loaded project by kind: same page calls the page's applyLoadedProject directly; cross-page goes through a localStorage channel before navigating
  dispatchLoadedProject(data) {
    const kind = data && data.kind;
    if (kind !== "blast" && kind !== "primer_design") {
      toast(t("top.unknown_project"));
      return;
    }
    const target = kind === "blast" ? "/blast.html" : "/design.html";
    const here = (location.pathname || "/").replace(/\/+$/, "") || "/";
    if (here === target) {
      window.PAGE?.applyLoadedProject?.(data);
      return;
    }
    try {
      if (kind === "blast") {
        // 表单状态 + 结果状态分开两条通道,blast.js 的 restoreState/restoreResult 各自读取
        // Form state and result state go through two separate channels; blast.js's restoreState/restoreResult each read their own
        const o = data.options || {};
        savePageState("blast", {
          seq: o.query_text || "",
          db: o.db || "",
          program: o.program || "",
          evalue: o.evalue != null ? String(o.evalue) : "",
          maxtargets: o.max_targets != null ? String(o.max_targets) : "",
          short: !!o.short_seq_mode,
          remote: !!o.remote,
          remoteDb: o.remote_db || "",
        });
        localStorage.setItem("bp_blast_result", JSON.stringify({
          _v: 1, taskId: null,
          result: { raw_output: data.raw_output || "", parsed: data.parsed || [], options: o },
          rawShown: 0,
        }));
      } else {
        localStorage.setItem("bp_design_load", JSON.stringify(data));
      }
    } catch (e) { /* 超 localStorage 限额:降级为仅跳页,页面按空态启动 */ }
    // Over localStorage quota: degrade to navigating only, page starts in empty state
    location.href = target;
  },
};

/* ================================ 工具函数 ================================ */
// Utility functions

async function api(url, opts = {}) {
  const init = {
    method: opts.method || "GET",
    headers: opts.json ? { "Content-Type": "application/json" } : undefined,
  };
  if (opts.body !== undefined) {
    if (opts.json) init.body = JSON.stringify(opts.body);
    else init.body = opts.body;
  }
  const res = await fetch(url, init);
  let data = null;
  const ct = res.headers.get("content-type") || "";
  if (ct.includes("json")) {
    try { data = await res.json(); } catch (e) { data = null; }
  } else {
    data = await res.text();
  }
  if (!res.ok) {
    const msg = (data && (data.detail || data.error || data.message)) || `HTTP ${res.status}`;
    throw new Error(typeof msg === "string" ? msg : JSON.stringify(msg));
  }
  return data;
}

/* 设计结果瘦身:剔除渲染不需要的大字段再入 localStorage/sessionStorage。
   step1_hsps 在大命中(高重复模板 × 全基因组库)时可含数十万 HSP,仅坐标
   版就 15~20 MB,加 qseq/sseq 全长序列实测 156 MB —— 超出 localStorage
   5 MB 配额,保存静默失败 → 切页/刷新后结果丢失。渲染只用 depth/stages/
   pairs/kmer_cache,step1_hsps 前端从不读取;kmer_cache 由调用方拆入
   sessionStorage。blast 页 raw_output 同理(见 blast.js saveResult)。 */
// Slim the design result before persisting: strip fields rendering never
// reads. step1_hsps on heavy hits (repeat-rich template × whole-genome db)
// can hold hundreds of thousands of HSPs — ~15-20 MB as coordinates alone,
// 156 MB measured with full qseq/sseq — beyond the 5 MB localStorage quota,
// so saving silently failed and results vanished on page switch/refresh.
// Rendering uses depth/stages/pairs/kmer_cache only; the caller splits
// kmer_cache into sessionStorage. blast.js handles raw_output the same way.
function slimDesignResult(result) {
  if (!result || !Array.isArray(result.results)) return result;
  return Object.assign({}, result, {
    results: result.results.map((q) => {
      if (!q || !("step1_hsps" in q)) return q;
      const { step1_hsps, ...rest } = q;
      return rest;
    }),
  });
}

/* ===== 结果持久化存储:localStorage 优先(同步、跨会话),超限回退 IndexedDB =====
   localStorage 5 MB 配额对设计结果(深度数组随模板线性增长)与 BLAST 结果
   (raw_output 全文可达数十 MB)都不够用,超限写入静默失败 → 结果丢失。
   IndexedDB 无配额限制,同一 get/set 接口对调用方透明:保存前先试
   localStorage(小结果走原路径,兼容旧会话),QuotaExceeded 时清掉旧值并
   落入 IndexedDB;读取先 localStorage 后 IndexedDB。 */
// Result persistence: localStorage first (synchronous, survives browser
// restarts), falling back to IndexedDB (no quota) when it overflows. The
// 5 MB localStorage quota cannot hold design results (depth arrays grow
// linearly with the template) nor BLAST results (raw_output runs to tens of
// MB), and a silent quota failure used to lose the result. IndexedDB has no
// quota; the get/set interface is transparent to callers: save tries
// localStorage first (small results keep the old path and stay compatible
// with older sessions), clears the stale value and falls into IndexedDB on
// QuotaExceededError; load checks localStorage then IndexedDB.
const IDB_NAME = "blastprime";
const IDB_STORE = "kv";
let _idbDb = null;
async function _idbOpen() {
  if (_idbDb) return _idbDb;
  return new Promise((resolve, reject) => {
    const req = indexedDB.open(IDB_NAME, 1);
    req.onupgradeneeded = () => {
      if (!req.result.objectStoreNames.contains(IDB_STORE)) {
        req.result.createObjectStore(IDB_STORE);
      }
    };
    req.onsuccess = () => { _idbDb = req.result; resolve(_idbDb); };
    req.onerror = () => reject(req.error);
  });
}
async function _idbSet(key, val) {
  const db = await _idbOpen();
  return new Promise((resolve, reject) => {
    const tx = db.transaction(IDB_STORE, "readwrite");
    tx.objectStore(IDB_STORE).put(val, key);
    tx.oncomplete = () => resolve();
    tx.onerror = () => reject(tx.error);
  });
}
async function _idbGet(key) {
  const db = await _idbOpen();
  return new Promise((resolve, reject) => {
    const tx = db.transaction(IDB_STORE, "readonly");
    const req = tx.objectStore(IDB_STORE).get(key);
    req.onsuccess = () => resolve(req.result === undefined ? null : req.result);
    req.onerror = () => reject(req.error);
  });
}
async function persistResult(key, envelope) {
  try {
    localStorage.setItem(key, JSON.stringify(envelope));
    return "local";
  } catch (e) {
    try {
      // 清掉 localStorage 里的旧值,避免 load 读到过期的小结果
      // Drop the stale localStorage value so load never reads an old small result
      localStorage.removeItem(key);
      await _idbSet(key, envelope);
      return "idb";
    } catch (e2) {
      return null;
    }
  }
}
async function loadResult(key) {
  try {
    const raw = localStorage.getItem(key);
    if (raw != null) return JSON.parse(raw);
  } catch (e) {
    localStorage.removeItem(key);
  }
  try { return await _idbGet(key); } catch (e) { return null; }
}
async function _idbDelete(key) {
  const db = await _idbOpen();
  return new Promise((resolve, reject) => {
    const tx = db.transaction(IDB_STORE, "readwrite");
    tx.objectStore(IDB_STORE).delete(key);
    tx.oncomplete = () => resolve();
    tx.onerror = () => reject(tx.error);
  });
}
/* 清结果:localStorage 与 IndexedDB 两处都要删(localStorage 超限时结果在
   IndexedDB 里,只删 localStorage 会让 load 读到陈旧结果)。 */
// Clear a result from both stores (on quota fallback the result lives in
// IndexedDB; deleting only localStorage would let load resurrect it)
async function clearResult(key) {
  try { localStorage.removeItem(key); } catch (e) {}
  try { await _idbDelete(key); } catch (e) {}
}

/* ===== 页面状态保存:localStorage 持久(切页/刷新/关浏览器均不丢) ===== */
// Persistent page-state save: localStorage (survives page switches, refreshes and browser restarts)
function savePageState(key, obj) {
  try { localStorage.setItem("bp_" + key, JSON.stringify(obj)); } catch (e) {}
}
function loadPageState(key) {
  try { return JSON.parse(localStorage.getItem("bp_" + key)) || null; } catch (e) { return null; }
}

/* 输入序列重名处理:仅提交/解析时改写(输入框保持原文)。完整标题行去重,
   取最小未用 " (N)" 后缀;含 `=` 的名称型行跳过(逐行独立解析);幂等。
   与后端 blast.dedupe_fasta_headers 同逻辑。 */
// Duplicate-header handling at submit/parse time only (the input box keeps
// the original text). Full-header dedupe with the smallest unused " (N)"
// suffix; name-type lines (containing "=") are skipped. Idempotent — mirrors
// the backend blast.dedupe_fasta_headers.
function dedupeFastaHeaders(text) {
  const seen = new Set();
  return (text || "").split("\n").map((line) => {
    if (!line.startsWith(">") || line.includes("=")) return line;
    const hdr = line.slice(1).trim();
    if (!seen.has(hdr)) { seen.add(hdr); return ">" + hdr; }
    let n = 1;
    while (seen.has(`${hdr} (${n})`)) n++;
    const cand = `${hdr} (${n})`;
    seen.add(cand);
    return ">" + cand;
  }).join("\n");
}

function toast(msg) {
  let el = document.getElementById("toast");
  if (!el) {
    el = document.createElement("div");
    el.id = "toast";
    document.body.appendChild(el);
  }
  el.textContent = msg;
  el.classList.add("show");
  clearTimeout(el._timer);
  el._timer = setTimeout(() => el.classList.remove("show"), 2600);
}

function escapeHtml(s) {
  return String(s ?? "").replace(/[&<>"']/g, (c) => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
  }[c]));
}

function applyTheme(theme) {
  // 缓存用户主题设置(而非解析后的值),供 <head> 预涂脚本消除首帧闪变
  // Cache the user's theme setting (not the resolved value) for the <head> pre-paint script to avoid first-frame flash
  try { localStorage.setItem("blastprime_theme", theme || "system"); } catch (e) {}
  if (theme === "system") {
    const dark = window.matchMedia("(prefers-color-scheme: dark)").matches;
    document.documentElement.dataset.theme = dark ? "dark" : "light";
  } else {
    document.documentElement.dataset.theme = theme;
  }
  // 通知当前页面重绘 canvas 等用 CSS 变量上色的绘制内容(如图形摘要的基因名底框/标尺),
  // 否则主题变量变化而位图不变(init 早期 PAGE 未就绪时可选链跳过)
  // Notify the current page to redraw canvas and other CSS-variable-colored drawings (e.g. the graphic summary's gene-name box/ruler),
  // otherwise the theme variables change while the bitmap stays stale (optional chaining skips early when PAGE is not ready during init)
  window.PAGE?.onThemeChanged?.();
}

/* 顶栏主题临时切换:点击在 浅色↔深色 间循环(按钮当前值存 dataset.theme,
   语言切换时据此重绘文本;不写入 config.json)。设置若为"跟随系统",
   按钮按当前系统偏好解析后的实际主题描述,点击后进入手动浅/深循环 */
// Temporary top-bar theme toggle: clicking cycles light<->dark (the button's current value lives in dataset.theme,
// used to redraw the text on language switch; not written to config.json). If the setting is "follow system",
// the button shows the actual theme resolved from the current system preference; clicking enters the manual light/dark cycle
const THEME_ICONS = { light: "☀", dark: "🌙" };

function resolveTheme(theme) {
  if (theme === "system") {
    return window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  }
  return theme;
}

function renderThemeSwitch(cur) {
  const btn = document.getElementById("theme-switch");
  if (!btn) return;
  const v = resolveTheme(cur);
  btn.dataset.theme = v;
  btn.textContent = `${THEME_ICONS[v] || "☀"} ${t("settings.theme_" + v)}`;
  btn.title = t("top.theme_tmp_title");
}

async function copyText(text) {
  try {
    await navigator.clipboard.writeText(text);
    toast(t("common.copied"));
  } catch (e) {
    const ta = document.createElement("textarea");
    ta.value = text;
    document.body.appendChild(ta);
    ta.select();
    document.execCommand("copy");
    ta.remove();
    toast(t("common.copied"));
  }
}

function downloadURL(url, fallbackName) {
  const a = document.createElement("a");
  a.href = url;
  a.download = fallbackName || "";
  document.body.appendChild(a);
  a.click();
  a.remove();
}

function modalShow(el) { el.classList.add("show"); }
function modalHide(el) { el.classList.remove("show"); }

/* 任务式请求:POST 启动任务 → 返回 task_id */
// Task-style request: POST to start a task -> returns task_id
async function startTask(url, formBody) {
  const data = await api(url, { method: "POST", body: formBody });
  return data.task_id;
}

/* 将选中的文件内容追加到文本框(不覆盖已有内容);支持多选,返回 {added, failed} */
// Append the selected files' contents to the textarea (without overwriting existing content); supports multi-select, returns {added, failed}
async function appendFilesToTextarea(fileInput, ta) {
  const files = [...(fileInput?.files || [])];
  let added = 0, failed = 0;
  for (const f of files) {
    try {
      const text = await f.text();
      if (!text.trim()) continue;
      if (ta.value && !ta.value.endsWith("\n")) ta.value += "\n";
      ta.value += text;
      if (!ta.value.endsWith("\n")) ta.value += "\n";
      added++;
    } catch (e) { failed++; }
  }
  if (fileInput) fileInput.value = "";  // 清空选择,允许重复选同一文件
  // Clear the selection so the same file can be picked again
  return { added, failed };
}

/* 将表单元素 + options JSON 组装为 multipart(与后端 Form/File 契约对齐) */
// Assemble form fields + options JSON into multipart (aligned with the backend Form/File contract)
function buildForm(fields, fileInput) {
  const fd = new FormData();
  Object.entries(fields).forEach(([k, v]) => {
    if (v !== undefined && v !== null && v !== "") fd.append(k, v);
  });
  if (fileInput && fileInput.files.length) {
    [...fileInput.files].forEach((f) => fd.append("files", f));
  }
  return fd;
}

/* 深色模式下的媒体查询变化(跟随系统时) */
// Media-query change in dark mode (when following the system)
window.matchMedia("(prefers-color-scheme: dark)").addEventListener?.("change", (e) => {
  if (App.settings?.theme === "system") {
    // 顶栏按钮手动切换过(含从 localStorage 恢复的临时主题)时,系统主题变化不覆盖临时选择
    // When the top-bar button was manually toggled (including a temporary theme restored from localStorage), system theme changes do not override the temporary choice
    const btn = document.getElementById("theme-switch");
    if (btn?.dataset.manual) return;
    document.documentElement.dataset.theme = e.matches ? "dark" : "light";
    // 顶栏临时按钮未手动切换过时,描述跟随系统更新
    // When the top-bar temporary button has not been manually toggled, keep its description following the system
    if (btn) renderThemeSwitch("system");
  }
});

document.addEventListener("DOMContentLoaded", () => App.init());
