/* 结果 JSON 后台线程解析:大结果(parsed 含每 HSP 全序列,可达数十 MB)
   JSON.parse 在主线程会阻塞数秒导致"页面没有响应",移到 Worker 里解析。 */
// Background-thread parsing of result JSON: large results (parsed contains full per-HSP sequences, up to tens of MB)
// would block the main thread for seconds during JSON.parse causing "page not responding"; parsing moved into the Worker.
"use strict";

self.onmessage = async (ev) => {
  const { url } = ev.data || {};
  try {
    const res = await fetch(url);
    if (!res.ok) { self.postMessage({ ok: false, error: "HTTP " + res.status }); return; }
    const data = await res.json();
    self.postMessage({ ok: true, data });
  } catch (e) {
    self.postMessage({ ok: false, error: String(e) });
  }
};
