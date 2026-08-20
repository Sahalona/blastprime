/* 数据库页:构建 / 记录 / 浏览 / 删除 */
// Database page: build / records / browse / delete
"use strict";

window.PAGE = {
  async init() {
    this.restoreState();
    this.bindBuild();
    this.bindModals();
    this.bindImport();
    this.bindScan();
    await this.loadRecords();
  },

  /* ---------------- 页面状态临时保存(切页复原,见 app.js savePageState) ---------------- */
  // Temporary page-state save (restored on page switch, see app.js savePageState)

  saveState() {
    savePageState("db", {
      paste: document.getElementById("db-paste").value,
      name: document.getElementById("db-name").value,
      type: document.getElementById("db-type").value,
    });
  },

  restoreState() {
    const s = loadPageState("db");
    if (!s) return;
    const set = (id, v) => { if (v != null) document.getElementById(id).value = v; };
    set("db-paste", s.paste);
    set("db-name", s.name);
    set("db-type", s.type);
    document.getElementById("db-paste").addEventListener("input", () => this.saveState());
    document.getElementById("db-name").addEventListener("input", () => this.saveState());
    document.getElementById("db-type").addEventListener("change", () => this.saveState());
  },

  /* ---------------- 扫描默认目录 ---------------- */
  // Scan default directory

  bindScan() {
    const btn = document.getElementById("btn-scan-dir");
    if (!btn) return;
    btn.addEventListener("click", async () => {
      btn.disabled = true;
      try {
        const { scanned } = await api("/api/db/scan", { method: "POST" });
        if (scanned > 0) toast(`${t("db.scan_found")}: ${scanned}`);
        else toast(t("db.scan_none"));
        this.loadRecords();
      } catch (e) {
        toast(e.message);
      } finally {
        btn.disabled = false;
      }
    });
  },

  /* ---------------- 导入已有数据库 ---------------- */
  // Import existing database

  bindImport() {
    const btn = document.getElementById("btn-import-db");
    const modal = document.getElementById("import-modal");
    if (!btn || !modal) return;
    btn.addEventListener("click", () => {
      document.getElementById("import-path").value = "";
      modalShow(modal);
    });
    modal.querySelector(".btn-cancel").addEventListener("click", () => modalHide(modal));
    modal.querySelector(".btn-do-import").addEventListener("click", async () => {
      const path = document.getElementById("import-path").value.trim();
      if (!path) { toast(t("db.path_empty")); return; }
      try {
        await api("/api/db/browse", { method: "POST", json: true, body: { index_file: path } });
        toast(t("db.import_ok"));
        modalHide(modal);
        this.loadRecords();
      } catch (e) { toast(t("db.import_fail") + ": " + e.message); }
    });
  },

  /* 浏览/删除/备注模态的关闭按钮 */
  // Close buttons for the browse/delete/note modals
  bindModals() {
    document.querySelector("#browse-modal .btn-close").addEventListener("click", () =>
      modalHide(document.getElementById("browse-modal")));
    document.querySelector("#del-modal .btn-cancel-del").addEventListener("click", () =>
      modalHide(document.getElementById("del-modal")));
    document.querySelector("#note-modal .btn-cancel-note").addEventListener("click", () =>
      modalHide(document.getElementById("note-modal")));
  },

  /* ---------------- 构建 ---------------- */
  // Build

  bindBuild() {
    const btn = document.getElementById("btn-build");
    btn.addEventListener("click", async () => {
      if (btn.disabled) return; // 构建中防重复提交
      // Prevent duplicate submission while building
      const name = document.getElementById("db-name").value.trim();
      if (!name) { toast(t("db.name_empty")); document.getElementById("db-name").classList.add("invalid"); return; }
      const type = document.getElementById("db-type").value;
      // 重名标题加后缀(仅提交/解析时改写,输入框保持原文)
      // Duplicate headers get " (N)" suffixes at submit/parse time only (the
      // input box keeps the original text)
      const paste = dedupeFastaHeaders(document.getElementById("db-paste").value).trim();
      const files = document.getElementById("db-files");
      if (!paste && !files.files.length) {
        toast(t("db.files_empty"));
        return;
      }
      const fd = buildForm({
        db_name: name, fasta_text: paste,
        // 类型:后端固定 auto,前端类型选择只影响展示提示
        // Type: the backend always uses auto; the frontend type selector only affects display hints
      }, files);
      // 点击立即进入"正在构建"(乐观状态):文件上传 + 后端接收可能耗时
      // 1-2s,不能留无反馈的空白;POST 返回真实任务 id 后由 watchBuild
      // 以真实任务接管(会先清掉乐观轮询,见 watchBuild 开头)
      // Show the running state immediately on click (optimistic): file
      // upload + backend acceptance can take 1-2s and must not leave a
      // feedback gap; once the POST returns the real task id, watchBuild
      // takes over with the real task (it clears the optimistic poll first)
      this.watchBuild(null, btn);
      try {
        const { task_id } = await api("/api/db/build", { method: "POST", body: fd });
        this.watchBuild(task_id, btn);
        toast(t("db.build_started") + ": " + name);
      } catch (e) {
        // 请求失败:恢复按钮与状态行
        // Request failed: restore the button and status line
        if (btn) btn.disabled = false;
        const status = document.getElementById("build-status");
        const hint = document.getElementById("build-status-hint");
        if (status) { status.className = "build-status err"; status.innerHTML = `✗ ${t("db.build_failed")}`; }
        if (hint) hint.textContent = "";
        toast(e.message);
      }
    });
  },

  /* 构建进行中:按钮禁用 + 状态行(spinner/秒表/结果),日志实时刷新。
     大库建库耗时较长,makeblastdb 无逐条进度,秒表让"仍在运行"可见可感。
     taskId 可为 null(乐观阶段,POST 尚未返回):秒表照走,拿到真实任务
     id 后再次调用本函数接管(开头会清掉前一个乐观轮询)。 */
  // While building: button disabled + status line (spinner/stopwatch/result), logs refresh live.
  // Building large databases takes a while and makeblastdb has no per-record progress, so the stopwatch makes "still running" visible and tangible.
  // taskId may be null (optimistic phase, the POST has not returned yet): the
  // stopwatch still runs; once the real id arrives, call this again to take
  // over (the previous optimistic poll is cleared at the top).
  watchBuild(taskId, btn) {
    const card = document.getElementById("build-log-card");
    const log = document.getElementById("build-log");
    const status = document.getElementById("build-status");
    const hint = document.getElementById("build-status-hint");
    card.style.display = "";
    log.textContent = "";
    if (btn) btn.disabled = true;
    if (this._buildTimer) clearInterval(this._buildTimer);
    this._buildTaskId = taskId;
    this._buildStarted = Date.now();
    this._buildElapsed = null;
    const started = Date.now();
    const setStatus = (cls, html) => {
      status.className = "build-status" + (cls ? " " + cls : "");
      status.innerHTML = html;
    };
    const fmt = (s) => (s < 60 ? s + "s" : Math.floor(s / 60) + "m" + (s % 60) + "s");
    const flush = () => {
      const s = App.tasks.get(taskId);
      if (!s) return;
      const tail = (s.logs || []).map((l) => l.msg).join("\n");
      if (tail) log.textContent = tail;
      log.scrollTop = log.scrollHeight;
    };
    const running = () => {
      const secs = Math.floor((Date.now() - started) / 1000);
      setStatus("", `<span class="spinner"></span>${t("db.build_running", { s: fmt(secs) })}`);
    };
    // 乐观状态:任务 id 未到时同样显示"正在构建"(上传窗口不留空白)
    // Optimistic state: show "building" even before the task id exists —
    // the upload window must not look blank
    running();
    hint.textContent = t("db.build_elsewhere");
    // 每 500ms 轮询任务快照(SSE 已在 App 中实时维护)+ 刷新秒表
    // Poll the task snapshot every 500ms (SSE is already maintained live in App) + refresh the stopwatch
    this._buildTimer = setInterval(() => {
      const s = taskId ? App.tasks.get(taskId) : null;
      if (taskId && !s) { clearInterval(this._buildTimer); return; }
      if (taskId) flush();
      if (!s || s.status === "running" || s.status === "pending") {
        running();
      } else {
        clearInterval(this._buildTimer);
        if (btn) btn.disabled = false;
        hint.textContent = "";
        const secs = Math.floor((Date.now() - started) / 1000);
        this._buildElapsed = secs;   // 固定耗时:语言切换重绘不再随时间增长
        if (s.status === "succeeded") {
          setStatus("done", `✓ ${t("db.build_done", { s: fmt(secs) })}`);
          this.loadRecords();
        } else if (s.status === "cancelled") {
          setStatus("err", `✗ ${t("db.build_canceled")}`);
        } else {
          setStatus("err", `✗ ${t("db.build_failed")}`);
        }
      }
    }, 500);
  },

  /* 语言切换后重绘构建状态行(动态 innerHTML 快照,setLang 不刷新) */
  // Re-render the build status line after a language switch (dynamic
  // innerHTML snapshots that setLang cannot refresh)
  _renderBuildStatus() {
    const s = this._buildTaskId ? App.tasks.get(this._buildTaskId) : null;
    if (!s) return;
    const status = document.getElementById("build-status");
    const hint = document.getElementById("build-status-hint");
    const fmt = (v) => (v < 60 ? v + "s" : Math.floor(v / 60) + "m" + (v % 60) + "s");
    if (!status) return;
    if (s.status === "running" || s.status === "pending") {
      if (hint) hint.textContent = t("db.build_elsewhere");
      const secs = Math.floor((Date.now() - (this._buildStarted || Date.now())) / 1000);
      status.className = "build-status";
      status.innerHTML = `<span class="spinner"></span>${t("db.build_running", { s: fmt(secs) })}`;
    } else {
      if (hint) hint.textContent = "";
      const secs = this._buildElapsed != null ? this._buildElapsed : 0;
      status.className = "build-status" + (s.status === "succeeded" ? " done" : " err");
      status.innerHTML = s.status === "succeeded" ? `✓ ${t("db.build_done", { s: fmt(secs) })}`
        : s.status === "cancelled" ? `✗ ${t("db.build_canceled")}` : `✗ ${t("db.build_failed")}`;
    }
  },

  /* ---------------- 记录列表 ---------------- */
  // Records list

  async loadRecords() {
    const wrap = document.getElementById("records-wrap");
    this.recordsWrap = wrap;
    try {
      const { records } = await api("/api/db/records");
      if (!records.length) {
        wrap.innerHTML = `<div class="empty">${t("db.records_empty")}</div>`;
        return;
      }
      const rows = await Promise.all(records.map(async (r) => {
        let info = { type: null, seq_count: null };
        try { info = await api(`/api/db/info?prefix=${encodeURIComponent(r.prefix)}`); }
        catch (e) { /* 索引缺失的记录由后端载入时剔除,这里容忍竞态 */ }
        // Records whose index files are missing are pruned by the backend on load; tolerate the race here
        return this.recordRow(r, info);
      }));
      wrap.innerHTML = `
        <table class="data db-table">
          <thead><tr>
            <th data-i18n="db.col_sort" title="${t("db.sort_hint")}">${t("db.col_sort")}</th>
            <th data-i18n="db.col_name">${t("db.col_name")}</th>
            <th data-i18n="db.col_type">${t("db.col_type")}</th>
            <th data-i18n="db.col_prefix">${t("db.col_prefix")}</th>
            <th data-i18n="db.col_created">${t("db.col_created")}</th>
            <th data-i18n="db.col_actions">${t("db.col_actions")}</th>
          </tr></thead>
          <tbody>${rows.join("")}</tbody>
        </table>`;
      rows.forEach((_, i) => this.bindRecordActions(records[i]));
      this.bindRowDrag();
    } catch (e) {
      wrap.innerHTML = `<div class="empty">${escapeHtml(e.message)}</div>`;
    }
  },

  /* 语言切换时重渲染:类型/来源等单元格是动态文本,不走 data-i18n 遍历 */
  // Re-render on language switch: type/source cells are dynamic text and don't go through the data-i18n pass
  onLangChanged() {
    this._renderBuildStatus();
    this.loadRecords();
  },

  recordRow(r, info) {
    const name = basename(r.prefix) || r.prefix;
    const type = info.type === "prot" ? t("db.type_prot") : t("db.type_nucl");
    const created = r.is_created ? `<span class="tag info">${t("db.created_by_app")}</span>` : "";
    const extra = info.seq_count ? `<span class="muted">${info.seq_count} seqs · ${info.total_len} bp</span>` : "";
    const note = r.note ? ` <span class="tag" title="${escapeHtml(r.note)}">${escapeHtml(r.note)}</span>` : "";
    return `
      <tr data-prefix="${escapeHtml(r.prefix)}" draggable="true">
        <td><span class="drag-handle" title="${t("db.sort_hint")}">⠿</span></td>
        <td><b>${escapeHtml(name)}</b>${note}<br><span class="muted">${extra}</span></td>
        <td>${escapeHtml(type)}</td>
        <td class="mono" style="font-size:11.5px; word-break:break-all">${escapeHtml(r.prefix)}</td>
        <td>${created}</td>
        <td>
          <button class="small act-note" data-i18n="db.note_btn">${t("db.note_btn")}</button>
          <button class="small act-browse" data-i18n="db.act_browse">${t("db.act_browse")}</button>
          <button class="small danger act-delete" data-i18n="db.act_delete">${t("db.act_delete")}</button>
        </td>
      </tr>`;
  },

  bindRecordActions(record) {
    const tr = document.querySelector(`tr[data-prefix="${CSS.escape(record.prefix)}"]`);
    if (!tr) return;
    tr.querySelector(".act-note").addEventListener("click", () => this.confirmNote(record.prefix, record.note));
    tr.querySelector(".act-browse").addEventListener("click", () => this.browse(record.prefix));
    tr.querySelector(".act-delete")?.addEventListener("click", () => this.confirmDelete(record.prefix));
  },

  /* ---------------- 手动排序(拖拽) ---------------- */
  // Manual sort (drag & drop)

  bindRowDrag() {
    const tbody = this.recordsWrap?.querySelector("tbody");
    if (!tbody) return;
    let dragged = null;   // 当前被拖动的行
    let prev = [];        // 拖动前的顺序(用于判断是否有变化)
    const order = () =>
      [...tbody.querySelectorAll("tr[data-prefix]")].map(tr => tr.dataset.prefix);
    const cleanup = () => {
      tbody.querySelectorAll(".dragging, .drop-hint").forEach(x =>
        x.classList.remove("dragging", "drop-hint"));
      dragged = null;
    };
    tbody.addEventListener("dragstart", (e) => {
      const tr = e.target.closest("tr[data-prefix]");
      if (!tr) return;
      dragged = tr;
      prev = order();
      tr.classList.add("dragging");
      e.dataTransfer.effectAllowed = "move";
      // Firefox 需要 setData 才允许 drop
      // Firefox requires setData for the drop to be allowed
      e.dataTransfer.setData("text/plain", tr.dataset.prefix);
    });
    tbody.addEventListener("dragover", (e) => {
      if (!dragged) return;
      e.preventDefault();
      e.dataTransfer.dropEffect = "move";
      const tr = e.target.closest("tr[data-prefix]");
      tbody.querySelectorAll(".drop-hint").forEach(x => x.classList.remove("drop-hint"));
      if (!tr || tr === dragged) return;
      // 光标落在目标行上半部插其前,下半部插其后;被拖行实时跟随
      // Insert before the target row when the cursor is in its upper half,
      // after it in the lower half; the dragged row follows live
      const r = tr.getBoundingClientRect();
      const before = e.clientY < r.top + r.height / 2;
      tr.classList.add("drop-hint");
      tbody.insertBefore(dragged, before ? tr : tr.nextSibling);
    });
    tbody.addEventListener("drop", async (e) => {
      if (!dragged) return;
      e.preventDefault();
      const next = order();
      cleanup();
      if (next.join("\n") === prev.join("\n")) return;   // 位置未变
      try {
        await api("/api/db/reorder", { method: "POST", json: true, body: { order: next } });
        toast(t("db.sort_saved"));
      } catch (err) {
        toast(err.message);
        this.loadRecords();   // 回滚到服务端顺序
      }
    });
    tbody.addEventListener("dragend", cleanup);
  },

  /* ---------------- 备注 ---------------- */
  // Note

  confirmNote(prefix, note) {
    const modal = document.getElementById("note-modal");
    const input = document.getElementById("note-input");
    const name = basename(prefix) || prefix;
    document.getElementById("note-prefix").innerHTML =
      `<b style="color:var(--fg)">${escapeHtml(name)}</b> <span class="mono">${escapeHtml(prefix)}</span>`;
    input.value = note || "";
    modal.querySelector(".btn-save-note").onclick = async () => {
      const val = input.value.trim();
      modalHide(modal);
      try {
        await api("/api/db/note", { method: "POST", json: true, body: { prefix, note: val } });
        toast(t("db.note_saved"));
        this.loadRecords();
      } catch (e) { toast(e.message); }
    };
    modalShow(modal);
    input.focus();
  },

  /* ---------------- 浏览 ---------------- */
  // Browse

  async browse(prefix) {
    const modal = document.getElementById("browse-modal");
    document.getElementById("browse-info").textContent = "…";
    document.getElementById("browse-entries").textContent = "";
    modalShow(modal);
    try {
      const [info, entries] = await Promise.all([
        api(`/api/db/info?prefix=${encodeURIComponent(prefix)}`).catch(() => null),
        api(`/api/db/entries?prefix=${encodeURIComponent(prefix)}`).catch(() => null),
      ]);
      const name = basename(prefix);
      const lines = [];
      if (info) {
        lines.push(`${name}  (${info.type === "prot" ? t("db.type_prot") : t("db.type_nucl")})`);
        if (info.seq_count != null) lines.push(`sequences: ${info.seq_count}`);
        if (info.total_len != null) lines.push(`total letters: ${info.total_len}`);
      }
      document.getElementById("browse-info").innerHTML = `<pre class="code">${escapeHtml(lines.join("\n"))}</pre>`;
      document.getElementById("browse-entries").textContent =
        (entries?.entries || []).slice(0, 500).join("\n");
    } catch (e) {
      document.getElementById("browse-info").textContent = e.message;
    }
  },

  /* ---------------- 删除(二选一:仅移除记录 / 删除数据库文件) ---------------- */
  // Delete (one of two: remove the record only / delete the database files)

  confirmDelete(prefix) {
    const modal = document.getElementById("del-modal");
    document.getElementById("del-prefix").textContent = prefix;
    modal.querySelector(".btn-remove-only").onclick = async () => {
      modalHide(modal);
      try {
        await api("/api/db/remove", { method: "POST", json: true, body: { prefix } });
        toast(t("db.removed"));
        this.loadRecords();
      } catch (e) { toast(e.message); }
    };
    modal.querySelector(".btn-delete-files").onclick = async () => {
      modalHide(modal);
      try {
        await api("/api/db/delete", { method: "POST", json: true, body: { prefix } });
        toast(t("db.deleted"));
        this.loadRecords();
      } catch (e) { toast(e.message); }
    };
    modalShow(modal);
  },
};
