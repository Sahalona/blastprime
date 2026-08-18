"""FastAPI 应用与路由:环境、配置、数据库、BLAST、引物分析/设计、任务 SSE、导出。

CLI 入口:uv run blastprime [--host --port --no-browser --loglevel --logfile --config]

FastAPI application and routes: environment, configuration, database, BLAST,
primer analysis/design, task SSE, and export.

CLI entry: uv run blastprime [--host --port --no-browser --loglevel --logfile --config]
"""

from __future__ import annotations

import asyncio
import csv
import io
import json
import logging
import os
import shutil
import sys
import tempfile
import threading
import time
import webbrowser
from pathlib import Path
from typing import Any

from fastapi import FastAPI, File, Form, HTTPException, Query, Request, UploadFile
from fastapi.responses import FileResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles

from Bio import SeqIO

from . import blast, database, seq_tools
from .config import DEFAULT_PRIMER_PARAMS, detect_system_lang, get_config
from .primer_metrics import analyze_short_sequence
from .primer_design import design_pipeline
from .tasks import STATUS_FAILED, manager as task_manager

log = logging.getLogger("blastprime")

app = FastAPI(title="BlastPrime Studio", version="1.0.0")


# ---------------------------------------------------------------- 静态托管 (Static hosting)

def _static_dir() -> Path:
    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        return Path(meipass) / "static"
    return Path(__file__).resolve().parent.parent / "static"


STATIC_DIR = _static_dir()
if STATIC_DIR.is_dir():
    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


def _fasta_queries(fa_path: str) -> list[dict]:
    """从查询 FASTA 提取 {name, seq}(与 outfmt0 的 Query= 名称一致,用完整 defline)。

    Extract {name, seq} from a query FASTA (names match the Query= of outfmt0,
    using the full defline).
    """
    try:
        return [{"name": rec.description.strip(), "seq": str(rec.seq).upper()}
                for rec in SeqIO.parse(fa_path, "fasta")]
    except Exception:
        return []


def _reconstruct_query(parsed: list[dict], query_name: str) -> str:
    """导入场景(-outfmt 0 文本)下从所选查询各 HSP 的 qseq 按查询坐标重建完整序列。
    仅当每个碱基都被至少一个 HSP 覆盖时返回,否则返回空串。

    Reconstruct the full sequence of the selected query from the qseq of its
    HSPs by query coordinates (import scenario, -outfmt 0 text). Returns it only
    when every base is covered by at least one HSP; otherwise returns empty.
    """
    for q in parsed:
        if query_name and q.get("name") != query_name:
            continue
        qlen = q.get("length") or 0
        if not qlen:
            continue
        buf = [""] * qlen
        covered: set[int] = set()
        for subj in q.get("subjects", []):
            for h in subj.get("hits", []):
                qs, qe = h.get("qstart"), h.get("qend")
                qseq = h.get("qseq") or ""
                if not qs or not qe or not qseq:
                    continue
                pos = int(qs)
                for ch in qseq:
                    if ch == "-":
                        continue
                    if 1 <= pos <= qlen:
                        buf[pos - 1] = ch
                        covered.add(pos)
                    pos += 1
        if len(covered) == qlen:
            return "".join(buf)
    return ""


@app.middleware("http")
async def _no_cache_static(request, call_next):
    """静态资源强制每次重新验证,避免浏览器启发式缓存旧版前端文件。

    Force static resources to be revalidated on every request, preventing the
    browser from heuristically caching stale frontend files.
    """
    response = await call_next(request)
    if request.url.path.startswith("/static"):
        response.headers["Cache-Control"] = "no-cache"
    return response


@app.get("/")
async def index() -> FileResponse:
    idx = STATIC_DIR / "index.html"
    if idx.exists():
        return FileResponse(str(idx))
    return Response("BlastPrime Studio: static/index.html 缺失", media_type="text/plain")


_PAGE_ROUTES = {
    "blast.html": "blast.html",
    "design.html": "design.html",
}


@app.get("/{page}")
async def page(page: str) -> FileResponse:
    """blast.html / design.html 直接可访问(导航链接指向的路径)。

    Make blast.html / design.html directly accessible (the paths the
    navigation links point to).
    """
    if page in _PAGE_ROUTES:
        f = STATIC_DIR / _PAGE_ROUTES[page]
        if f.exists():
            return FileResponse(str(f))
    raise HTTPException(404, "页面不存在")


# ---------------------------------------------------------------- 工具 (Utilities)

def _sse(data: dict) -> str:
    return f"data: {json.dumps(data, ensure_ascii=False)}\n\n"


def _save_uploads(uploads: list[UploadFile]) -> str:
    """把上传文件写入临时目录,返回目录路径(任务结束由 _cleanup_upload_dir 删除)。

    Write uploaded files to a temporary directory and return its path
    (removed by _cleanup_upload_dir when the task finishes).
    """
    td = Path(tempfile.mkdtemp(prefix="blastprime_up_"))
    for u in uploads:
        if u is None:
            continue
        dest = td / (u.filename or "upload.fa")
        with dest.open("wb") as f:
            shutil.copyfileobj(u.file, f)
    return str(td)


def _cleanup_upload_dir(path: str | None) -> None:
    if path:
        shutil.rmtree(path, ignore_errors=True)


def _parse_target(raw: Any) -> tuple[int, int] | None:
    """解析 target 参数(前端可传 [start, end] 或 null)。

    Parse the target parameter (the frontend may pass [start, end] or null).
    """
    if not raw:
        return None
    try:
        s, e = int(raw[0]), int(raw[1])
        if e >= s and s >= 1:
            return (s, e)
    except (TypeError, ValueError, IndexError):
        pass
    return None


# ---------------------------------------------------------------- 环境与设置 (Environment & settings)

@app.get("/api/env")
async def api_env() -> dict:
    env = blast.detect_env()
    cfg = get_config()
    env["blast_bin_dir"] = cfg.get("blast_bin_dir", "")
    return env


@app.get("/api/settings")
async def api_settings_get() -> dict:
    cfg = get_config()
    return {
        "lang": cfg.get("lang", "zh"), "theme": cfg.get("theme", "system"),
        "loglevel": cfg.get("loglevel", "INFO"), "logfile": cfg.get("logfile", ""),
        "blast_bin_dir": cfg.get("blast_bin_dir", ""),
        "db_records": cfg.db_records(), "primer_params": cfg.primer_params(),
    }


@app.post("/api/settings")
async def api_settings_post(body: dict) -> dict:
    cfg = get_config()
    for key in ("lang", "theme", "loglevel", "logfile"):
        if key in body:
            cfg.set(key, body[key])
    if "blast_bin_dir" in body:
        blast.set_manual_bin_dir(body["blast_bin_dir"] or None)
    if "primer_params" in body and isinstance(body["primer_params"], dict):
        cfg.set_primer_params(body["primer_params"])
    return await api_settings_get()


# ---------------------------------------------------------------- 数据库 (Database)

@app.post("/api/db/build")
async def api_db_build(
    db_name: str = Form(...),
    fasta_text: str = Form(""),
    files: list[UploadFile] = File(default=[]),
) -> dict:
    up = _save_uploads([f for f in files if f is not None])
    text_path = None
    if fasta_text.strip():
        text_path = os.path.join(up, "_pasted.fa")
        with open(text_path, "w", encoding="utf-8") as f:
            # 重名标题去重(提交/解析时改写,输入框保持原文)
            # Duplicate headers deduped at submit/parse time (input box keeps
            # the original text)
            f.write(blast.dedupe_fasta_headers(fasta_text))

    def worker(cancel, on_log, on_progress) -> dict:
        try:
            paths = [os.path.join(up, p) for p in os.listdir(up)]
            if not paths:
                raise blast.BlastError("请提供 FASTA 文件或粘贴序列")
            result = blast.build_database(paths, "auto", db_name, None,
                                          on_log=on_log, on_progress=on_progress, cancel=cancel)
            on_log(f"建库完成: {result['seq_count']} 条序列 → {result['prefix']}")
            return result
        finally:
            _cleanup_upload_dir(up)

    task = task_manager.start("build_db", f"构建数据库 {db_name}", worker)
    return {"task_id": task.id}


@app.get("/api/db/records")
async def api_db_records() -> dict:
    cfg = get_config()
    return {"records": cfg.db_records()}


@app.post("/api/db/reorder")
async def api_db_reorder(body: dict) -> dict:
    """按用户手动排序重排数据库记录(持久化到 config.json)。

    Reorder database records by the user's manual sort (persisted to
    config.json).
    """
    order = body.get("order")
    if not isinstance(order, list):
        raise HTTPException(400, "缺少排序列表")
    prefixes = [str(p) for p in order if isinstance(p, str) and p.strip()]
    database.reorder_records(prefixes)
    return {"ok": True}


@app.post("/api/db/scan")
async def api_db_scan() -> dict:
    """手动触发:扫描默认目录,自动识别复制进来的数据库(幂等)。

    Manual trigger: scan the default directory and auto-detect databases
    that were copied in (idempotent).
    """
    return {"scanned": database.scan_default_dir()}


@app.post("/api/db/browse")
async def api_db_browse(body: dict) -> dict:
    index_file = (body.get("index_file") or "").strip()
    if not index_file:
        raise HTTPException(400, "缺少索引文件路径")
    try:
        return database.browse_existing(index_file)
    except ValueError as e:
        raise HTTPException(400, str(e))


@app.post("/api/db/remove")
async def api_db_remove(body: dict) -> dict:
    prefix = (body.get("prefix") or "").strip()
    if not prefix:
        raise HTTPException(400, "缺少库前缀路径")
    cfg = get_config()
    removed = cfg.remove_db_record(prefix)
    return {"removed": removed}


@app.post("/api/db/note")
async def api_db_note(body: dict) -> dict:
    """给数据库记录设置备注(仅存配置,不落盘到库目录)。

    Set a note on a database record (stored in config only, never written
    into the database directory).
    """
    prefix = (body.get("prefix") or "").strip()
    if not prefix:
        raise HTTPException(400, "缺少库前缀路径")
    note = str(body.get("note") or "").strip()
    try:
        database.set_note(prefix, note)
    except KeyError as e:
        raise HTTPException(404, str(e))
    return {"ok": True}


@app.post("/api/db/delete")
async def api_db_delete(body: dict) -> dict:
    """物理删除(任意来源的库均可),前端已弹选择框由用户确认。

    Physically delete a database (any source works); the frontend already
    shows a confirmation dialog for the user.
    """
    prefix = (body.get("prefix") or "").strip()
    if not prefix:
        raise HTTPException(400, "缺少库前缀路径")
    deleted = database.physically_delete(prefix)
    cfg = get_config()
    cfg.remove_db_record(prefix)
    return {"deleted": deleted}


@app.get("/api/db/info")
async def api_db_info(prefix: str = Query(...)) -> dict:
    if not database.record_valid(prefix):
        raise HTTPException(404, "数据库索引文件不存在")
    info = blast.blastdbcmd_info(prefix)
    info["prefix"] = prefix
    info["type"] = database.db_type_of(prefix)
    return info


@app.get("/api/db/entries")
async def api_db_entries(prefix: str = Query(...)) -> dict:
    if not database.record_valid(prefix):
        raise HTTPException(404, "数据库索引文件不存在")
    return {"entries": blast.list_entries(prefix)}


# ---------------------------------------------------------------- 序列工具 (Sequence tools)

@app.post("/api/seq/revcomp")
async def api_seq_revcomp(body: dict) -> dict:
    return {"text": seq_tools.reverse_complement(body.get("text") or "")}


@app.post("/api/seq/clean")
async def api_seq_clean(body: dict) -> dict:
    return {"text": seq_tools.clean_sequence(body.get("text") or "")}


@app.get("/api/seq/tables")
async def api_seq_tables() -> dict:
    return {"tables": seq_tools.codon_tables()}


@app.post("/api/seq/translate6")
async def api_seq_translate6(body: dict) -> dict:
    return seq_tools.translate6(body.get("text") or "", int(body.get("table") or 1))


# ---------------------------------------------------------------- BLAST 比对 (BLAST alignment)

@app.post("/api/blast/run")
async def api_blast_run(options: str = Form(...),
                        files: list[UploadFile] = File(default=[])) -> dict:
    """启动比对任务。options 为 JSON 字符串,查询文件经 multipart 上传。

    Start an alignment task. options is a JSON string; query files are
    uploaded via multipart.
    """
    try:
        opts = json.loads(options) if options else {}
    except json.JSONDecodeError:
        raise HTTPException(400, "options 不是有效 JSON")
    if not isinstance(opts, dict):
        raise HTTPException(400, "options 必须是对象")
    up = _save_uploads([f for f in files if f is not None])
    qfa_own = None

    def worker(cancel, on_log, on_progress) -> dict:
        nonlocal qfa_own
        try:
            if opts.get("remote") and not (opts.get("remote_db") or "").strip():
                raise blast.BlastError("远程模式需要填写远程库名")
            # `#` 注释行自动忽略 + 重名标题去重(blast/design/建库统一行为,
            # 仅提交/解析时改写,输入框保持原文)
            # `#` comment lines auto-ignored + duplicate headers deduped
            # (shared by blast/design/db-build; submit/parse time only)
            query_text = blast.dedupe_fasta_headers(blast.strip_comment_lines(
                str(opts.get("query_text") or ""))).strip()
            qfa = None
            run_db = None
            if query_text:
                first_line = query_text.split("\n")[0].strip()
                # 名称型 query(除 gene 外全部可选:range 默认全长、database
                # 默认 targetbase、再默认当前查询库;targetbase 同时决定 blast
                # 运行库):
                # Name-style query (everything but the gene is optional: range
                # defaults to full length, database to targetbase, then to the
                # page db; targetbase also selects the blast run db):
                # 首行含 key=value 语法(逗号),或首行后无任何序列字符(裸 >gene
                # 即名称查询);普通 FASTA(>id + 序列)不受影响
                # The first line carries key=value syntax (comma-separated), or
                # has no sequence characters after it (a bare >gene is a name
                # query); plain FASTA (>id + sequence) is unaffected.
                qq = blast.parse_name_query(first_line) if first_line.startswith(">") else None
                rest = "".join(query_text.split("\n")[1:]) if "\n" in query_text else ""
                has_seq = any(c in "ACGTUNacgtun" for c in rest)
                if qq and ("," in first_line or not has_seq):
                    resolved = blast.resolve_name_query(query_text,
                                                        opts.get("db") or None)
                    # targetbase 同时作用于 blast 页:运行库 = targetbase 解析
                    # 结果(远程模式忽略 targetbase,只影响设计页特异性比对);
                    # targetbase=False 关闭特异性 → 无 targetbase 语义,回落
                    # 页面库(blast 比对页面无影响)
                    # targetbase also applies to the blast page: run db =
                    # resolved targetbase (remote mode ignores it);
                    # targetbase=False disables specificity -> same as no
                    # targetbase, falling back to the page db (the blast page
                    # is unaffected)
                    if qq["targetbase"] and not opts.get("remote") \
                            and not blast.targetbase_disabled(qq["targetbase"]):
                        run_db = blast.resolve_db_ref(qq["targetbase"])
                elif first_line.startswith(">"):
                    resolved = query_text
                else:
                    # 纯序列:自动补 >query 头
                    # Bare sequence: automatically prepend the >query header
                    seq = "".join(query_text.split())
                    if not seq:
                        raise blast.BlastError("查询序列为空")
                    resolved = f">query\n{seq}\n"
                if not resolved.strip():
                    raise blast.BlastError("名称解析失败")
                fd, qfa_own = tempfile.mkstemp(suffix=".fa", prefix="blastprime_q_")
                with os.fdopen(fd, "w", encoding="utf-8") as f:
                    f.write(resolved)
                qfa = qfa_own
            elif up and os.listdir(up):
                qfa = os.path.join(up, os.listdir(up)[0])
            else:
                raise blast.BlastError("请提供查询序列(文本或文件)")

            # 解析查询序列(名称 + 序列),供引物分析复用(不重新跑 BLAST)
            # Parse the query sequences (name + sequence) for reuse by primer
            # analysis (without re-running BLAST)
            queries = _fasta_queries(qfa)

            if run_db is None:
                run_db = opts.get("db") or None
            result = blast.run_blast(
                run_db, qfa,
                program=str(opts.get("program") or "auto"),
                query_type=str(opts.get("query_type") or "auto"),
                evalue=str(opts.get("evalue") or "10"),
                max_targets=int(opts.get("max_targets") or 500),
                matrix=str(opts.get("matrix") or ""),
                extra_args=str(opts.get("extra_args") or ""),
                short_seq_mode=bool(opts.get("short_seq_mode")),
                remote=bool(opts.get("remote")),
                remote_db=str(opts.get("remote_db") or "nr"),
                on_log=on_log, on_progress=on_progress, cancel=cancel, timeout=600,
            )
            result["options"] = {**opts, "db": run_db}
            result["queries"] = queries
            if on_log:
                on_log(f"比对完成: {len(result['parsed'])} 条查询, "
                       f"{sum(len(q.get('subjects', [])) for q in result['parsed'])} 条命中")
            return result
        finally:
            if qfa_own:
                try:
                    os.unlink(qfa_own)
                except OSError:
                    pass
            _cleanup_upload_dir(up)

    task = task_manager.start("blast", "BLAST 比对", worker)
    return {"task_id": task.id}


@app.post("/api/blast/import")
async def api_blast_import(body: dict) -> dict:
    content = (body.get("content") or "").strip()
    if not content:
        raise HTTPException(400, "缺少导入内容")
    if "Query=" not in content:
        raise HTTPException(400, "需要成对文本格式(-outfmt 0),其他格式请先转换")
    parsed = blast.parse_outfmt0(content)
    if not parsed:
        raise HTTPException(400, "解析失败:未找到有效的 Query 块,请检查文件")
    return {"raw_output": content, "parsed": parsed}


# ---------------------------------------------------------------- 模块三:引物分析 (Module 3: primer analysis)

@app.post("/api/analysis/primer")
async def api_analysis_primer(body: dict) -> dict:
    seq = (body.get("seq") or "").strip().upper()
    query_name = str(body.get("query_name") or "").strip()
    queries = body.get("queries") or []
    parsed = body.get("parsed")
    db_prefix = (body.get("db_prefix") or "").strip() or None

    # 序列来源:显式 seq → 本次比对结果中的查询序列 → 导入场景从 HSP 的 qseq 重建
    # Sequence sources: explicit seq → the query sequences from the current
    # alignment result → reconstructed from HSP qseq in the import scenario
    if not seq and query_name:
        for q in queries:
            if q.get("name") == query_name:
                seq = str(q.get("seq") or "").upper()
                break
    if not seq and parsed:
        seq = _reconstruct_query(parsed, query_name)
    seq = seq.strip()
    if not seq:
        if parsed:
            raise HTTPException(400, "无法从比对结果重建查询序列(无命中或覆盖不全)")
        raise HTTPException(400, "缺少序列")
    if any(c not in "ACGTUN" for c in seq):
        raise HTTPException(400, "序列含非法字符(引物分析仅支持核酸序列)")

    hsps = None
    if parsed:
        # 从 parsed(Query→Subject→HSP)扁平化所选查询的 HSP
        # Flatten the HSPs of the selected query from parsed
        # (Query → Subject → HSP)
        for q in parsed:
            if query_name and q.get("name") != query_name:
                continue
            flat: list[dict] = []
            for subj in q.get("subjects", []):
                for h in subj.get("hits", []):
                    h = dict(h)
                    h.setdefault("sseqid", subj.get("name", ""))
                    if "length" not in h:
                        qs, qe = h.get("qstart"), h.get("qend")
                        if qs is not None and qe is not None:
                            h["length"] = abs(qe - qs) + 1
                    if "mismatch" not in h and h.get("ident_frac") is not None:
                        h["mismatch"] = int(round(
                            (1.0 - h["ident_frac"]) * h.get("length", 0)))
                    flat.append(h)
            if flat:
                hsps = flat
            break
    elif db_prefix and database.record_valid(db_prefix):
        # 未带比对结果:用 blastn-short 现跑(短序列模式)
        # No alignment result provided: run blastn-short on the fly
        # (short-sequence mode)
        fd, qfa = tempfile.mkstemp(suffix=".fa", prefix="blastprime_aq_")
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                f.write(f">query\n{seq}\n")
            r = blast.run_blast_tabular(db_prefix, qfa, program="blastn",
                                        evalue=1000, max_targets=1000)
        finally:
            try:
                os.unlink(qfa)
            except OSError:
                pass
        if r:
            hsps = r
    report = analyze_short_sequence(seq, hsps)
    return report


# ---------------------------------------------------------------- 名称型模板提取(guide_sup1.md §5) (Name-style template extraction, guide_sup1.md §5)

def _extract_template(db_prefix: str, entry: str, smin: int, smax: int,
                      strand: str, flank: int) -> dict:
    """侧翼模板提取(§5):[smin−flank, smax+flank],严格同条目内截断,
    负链反向互补,返回模板 FASTA + 模板内目标区间 + 坐标换算信息。
    extract 端点与 /api/design/template 共用。

    Flanking template extraction (§5): [smin−flank, smax+flank], strictly
    truncated within the same entry; reverse-complemented on the minus
    strand; returns the template FASTA, the target range within the template,
    and coordinate conversion info. Shared with the /api/design/template
    extract endpoint.
    """
    if smin > smax:
        smin, smax = smax, smin
    minus = strand in ("minus", "-")

    entry_len = blast.entry_length(db_prefix, entry)
    req_s, req_e = smin - flank, smax + flank
    ext_s = max(1, req_s)
    ext_e = min(entry_len, req_e)
    if ext_e < ext_s:
        raise HTTPException(
            400, f"目标区间超出库条目 '{entry}' 范围(全长 {entry_len} bp)")
    fa = blast.fetch_region(db_prefix, entry, ext_s, ext_e)
    try:
        rec = SeqIO.read(io.StringIO(fa), "fasta")
    except Exception:
        raise HTTPException(400, f"无法从库条目 '{entry}' 提取区间 {ext_s}-{ext_e}")
    seq = str(rec.seq).upper()
    if minus:
        seq = blast.revcomp(seq)
    # 目标 HSP 在模板内的 1-based 区间(§5.3 换算规则的逆)
    # 1-based range of the target HSP within the template
    # (inverse of the §5.3 conversion rule)
    if minus:
        t_start = ext_e - smax + 1
        t_end = ext_e - smin + 1
    else:
        t_start = smin - ext_s + 1
        t_end = smax - ext_s + 1
    lines = [f">{entry}"]
    lines += [seq[i:i + 60] for i in range(0, len(seq), 60)]
    return {
        "entry": entry, "entry_len": entry_len,
        "strand": "minus" if minus else "plus",
        "range": [smin, smax],   # 目标区间(条目上的基因组坐标,1-based)
        # Target range (genomic coordinates on the entry, 1-based)
        "requested": {"start": req_s, "end": req_e},
        "extract": {"start": ext_s, "end": ext_e},
        "truncated": req_s < ext_s or req_e > ext_e,
        "flank_len": flank,
        "target": {"start": t_start, "end": t_end},
        "template_len": len(seq),
        "template_fasta": "\n".join(lines) + "\n",
    }


@app.post("/api/design/template")
async def api_design_template(body: dict) -> dict:
    """名称型模板解析(§5): 首行 `>[gene],range=...,database=...,name=...,
    targetbase=...` → 按 [range−flank, range+flank] 同条目截取返回模板;
    非名称型 → {named:false}。除 gene 外全部可选:range 缺省 = 全长;
    database 缺省 = targetbase,再缺省 = body.db_prefix;名称型恒 plus 链。
    响应新增 db_prefix(条目库)、spec_db(targetbase||database 解析结果,
    null = 页面库)、spec_skip(targetbase=False → 设计跳过特异性比对)、
    display_name(name||gene,模板 FASTA 头)。

    Name-style template resolution (§5): `>[gene],range=...,database=...,
    name=...,targetbase=...` → truncate the template from the same entry
    over [range−flank, range+flank]; non-name-style → {named:false}. Everything
    but the gene is optional: range defaults to full length; database defaults
    to targetbase, then to body.db_prefix. Always the plus strand.
    targetbase=False disables specificity for the design (spec_skip:true).
    """
    text = (body.get("template_text") or "").strip()
    if not text:
        raise HTTPException(400, "缺少模板序列")
    first_line = text.splitlines()[0].strip()
    q = blast.parse_name_query(first_line) if first_line.startswith(">") else None
    # 名称型判定与 blast 页同 gate:首行含逗号(key=value 语法),或首行后
    # 无任何序列字符(裸 >gene 即名称查询);普通 FASTA(>id + 序列)不受影响
    # Name-type detection shares the blast page's gate: the first line has a
    # comma (key=value syntax), or no sequence characters follow it (a bare
    # >gene is a name query); plain FASTA (>id + sequence) is unaffected.
    rest = "".join(text.split("\n")[1:]) if "\n" in text else ""
    has_seq = any(c in "ACGTUNacgtun" for c in rest)
    if not q or not ("," in first_line or not has_seq):
        return {"named": False}
    rng = q["range"]
    db_prefix = (body.get("db_prefix") or "").strip()
    try:
        flank = max(0, int(body.get("flank_len") or 0))
    except (TypeError, ValueError):
        raise HTTPException(400, "侧翼长度参数非法")
    try:
        # targetbase=False → 关闭特异性比对(spec_skip),不解析为库:
        # db_prefix 回落 database= / 页面库;spec_db 保持 null(深度统计照常
        # 走页面库,仅跳过 blastn-short 逆向验证 —— 与 skip_spec_eval 同语义)
        # targetbase=False -> disable specificity (spec_skip), never resolve it
        # as a db: db_prefix falls back to database= / the page db; spec_db
        # stays null (depth counting still runs on the page db; only the
        # blastn-short reverse validation is skipped — same semantics as
        # skip_spec_eval)
        spec_skip = blast.targetbase_disabled(q["targetbase"])
        # database= 存在 → 解析;否则 targetbase → 解析;再否则 body.db_prefix
        # database= if present, else targetbase, else body.db_prefix
        if q["database"]:
            db_prefix = blast.resolve_db_ref(q["database"])
        elif q["targetbase"] and not spec_skip:
            db_prefix = blast.resolve_db_ref(q["targetbase"])
        # spec_db:特异性比对库 = targetbase||database 的解析结果,null = 页面库
        # spec_db: the specificity-comparison db = resolved targetbase||database,
        # null means the page db
        spec_db = None
        if q["targetbase"] and not spec_skip:
            spec_db = blast.resolve_db_ref(q["targetbase"])
        elif q["database"]:
            spec_db = blast.resolve_db_ref(q["database"])
    except blast.BlastError as e:
        raise HTTPException(400, str(e))
    if not db_prefix or not database.record_valid(db_prefix):
        raise HTTPException(400, "请选择有效的核酸基因组数据库")
    try:
        # 先解析出真实条目名(resolve_entry 支持模糊包含匹配,entry_length 不支持模糊名)
        # Resolve the real entry name first (resolve_entry supports fuzzy
        # substring matching; entry_length does not accept fuzzy names)
        fa = blast.resolve_entry(db_prefix, q["gene"])
        rec = SeqIO.read(io.StringIO(fa), "fasta")
        entry = rec.id
        entry_len = len(rec.seq)
        # range 解析(缺省 = 全长;支持 * 通配,1-based 闭区间)
        # Range parsing (default = full length; supports * wildcards, 1-based
        # closed interval)
        start, end = 1, entry_len
        if rng:
            a, b = rng.split("-")
            a, b = a.strip(), b.strip()
            try:
                start = 1 if a == "*" else int(a)
                end = entry_len if b == "*" else int(b)
            except ValueError:
                raise HTTPException(400, f"范围参数非法: range={rng}")
            if start < 1 or end > entry_len or start > end:
                raise HTTPException(
                    400, f"范围越界: range={rng} (条目 '{entry}' 全长 {entry_len} bp)")
        out = {"named": True, "db_prefix": db_prefix, "spec_db": spec_db,
               "spec_skip": spec_skip,   # targetbase=False:设计跳过特异性比对
               "display_name": q["name"] or q["gene"],
               **_extract_template(db_prefix, entry, start, end, "plus", flank)}
    except blast.BlastError as e:
        # 基因不存在/取序列失败 → 4xx 带后端中文消息,前端 toast 展示
        # Missing gene / sequence fetch failure -> 4xx with the backend message
        raise HTTPException(400, str(e))
    # 模板 FASTA 头 = 显示名(name||gene)
    # Template FASTA header = the display name (name||gene)
    out["template_fasta"] = ">" + out["display_name"] + "\n" + \
        out["template_fasta"].split("\n", 1)[1]
    return out


@app.post("/api/design/extract")
async def api_design_extract(body: dict) -> dict:
    """侧翼模板提取(§5),见 _extract_template;保留兼容旧调用方。

    Flanking template extraction (§5); see _extract_template. Kept for
    backward compatibility with older callers.
    """
    db_prefix = (body.get("db_prefix") or "").strip()
    entry = (body.get("entry") or "").strip()
    if not db_prefix or not database.record_valid(db_prefix):
        raise HTTPException(400, "请选择有效的核酸基因组数据库")
    if not entry:
        raise HTTPException(400, "缺少目标库条目名")
    try:
        smin, smax = int(body.get("sbjct_min")), int(body.get("sbjct_max"))
        flank = max(0, int(body.get("flank_len") or 0))
    except (TypeError, ValueError):
        raise HTTPException(400, "坐标或侧翼长度参数非法")
    return _extract_template(
        db_prefix, entry, smin, smax,
        str(body.get("strand") or "plus"), flank)


# ---------------------------------------------------------------- 引物设计 (Primer design)

@app.post("/api/design/run")
async def api_design_run(body: dict) -> dict:
    # 定位模式(locate_jobs)自带模板,无需顶层 template(§5);缺省则走主流程
    # Locate mode (locate_jobs) carries its own templates, so no top-level
    # template is needed (§5); the default falls through to the main flow
    locate_jobs = body.get("locate_jobs")
    template = (body.get("template") or "").strip()
    if not template and not locate_jobs:
        raise HTTPException(400, "缺少模板序列")
    db_prefix = (body.get("db_prefix") or "").strip()
    if not db_prefix or not database.record_valid(db_prefix):
        raise HTTPException(400, "请选择有效的核酸基因组数据库")
    params = {**DEFAULT_PRIMER_PARAMS, **((body.get("params") or {}))}
    target = _parse_target(body.get("target"))
    # 特异性比对库(可选):name=/targetbase= 解析结果;缺省 = db_prefix
    # Specificity-comparison db (optional): the name-style targetbase= result;
    # defaults to db_prefix
    spec_db = (body.get("spec_db") or "").strip() or None
    if spec_db and not database.record_valid(spec_db):
        raise HTTPException(400, "spec_db 不是有效的库记录")
    # k-mer 结果缓存(R17):前端 sessionStorage/项目文件带回的既往 k-mer 分析
    # (键含序列/库/尺度/范围,后端校验,命中即跳过第一步 blastn 与计数)
    # K-mer result cache (R17): prior k-mer analyses carried back from the
    # frontend sessionStorage/project file (keyed by seq/db/scales/ranges;
    # the backend validates and skips step-1 blastn + counting on a hit)
    kmer_caches = body.get("kmer_caches") or None

    def worker(cancel, on_log, on_progress) -> dict:
        if locate_jobs:
            results = []
            succeeded = 0
            total = len(locate_jobs)
            for i, job in enumerate(locate_jobs):
                if cancel is not None and cancel.cancelled:
                    break
                job = job or {}
                tmpl = str(job.get("template") or "").strip()
                ctx = job.get("ctx") or {}
                tgt = _parse_target(job.get("target"))
                # 每查询独立参数(§6 F/R 设计范围、相对产物换算等)
                # Per-query independent params (§6: F/R design ranges,
                # relative product conversion, etc.)
                job_params = {**params, **((job.get("params") or {}))}
                if on_log:
                    on_log(f"[{i + 1}/{total}] 定位模板设计: {ctx.get('entry', '?')}")
                # 进度分段:索引 0-10 在 design_pipeline 内联;查询级 10-100
                # 按 job 序号均摊(locate 模式逐 job 线性推进)。design_pipeline
                # 的 on_progress 输出已是 0~100(f ∈ [0,100]),须除以 100 还原成
                # 0~1 再映射 —— 直接乘会把 f=2.5 变成 235,任务进度恒钳满 100,
                # 状态栏进度条看起来"完全没用上"
                # Progress mapping: the 0-10 index phase is inline inside
                # design_pipeline; the query-level 10-100 is spread across
                # jobs by index (locate mode advances linearly per job).
                # design_pipeline's on_progress already emits 0~100 (f ∈ [0,100]);
                # divide by 100 before mapping — feeding 2.5 straight in yields
                # 235, clamping the task progress at a permanent 100 and making
                # the status-bar progress bar look entirely unused
                job_prog = (lambda f: on_progress(10.0 + 90.0 * (i + f / 100.0) / total)
                            if on_progress else None)
                try:
                    # design_pipeline 本身按查询分组返回 {"results": [...]}(支持多 FASTA),
                    # design_pipeline returns {"results": [...]} grouped by
                    # query (supports multi-FASTA);
                    # 此处每 job 单条模板,取回第一条(即本 job 的设计结果)
                    # each job here has a single template, so take the first
                    # entry (that job's design result)
                    res = design_pipeline(tmpl, db_prefix, tgt, job_params,
                                          on_log=on_log, cancel=cancel,
                                          locate_ctx=ctx, spec_db=spec_db,
                                          on_progress=job_prog,
                                          kmer_caches=kmer_caches)
                    r0 = res["results"][0]
                    results.append(r0)
                    # 仅设计成功(阶段一~四产出引物)才计入 succeeded(此前无条件 +1,
                    # 全失败也误报"N/N 成功")
                    # Only designs that actually produced primers (stages 1-4)
                    # count toward succeeded; previously the counter incremented
                    # unconditionally, falsely reporting "N/N succeeded" even
                    # when all designs failed
                    succeeded += 1 if r0.get("success") else 0
                except Exception as e:
                    if cancel is not None and cancel.cancelled:
                        raise
                    results.append({"query": ctx.get("entry", "?"), "success": False,
                                    "error": str(e), "stage_reached": 0,
                                    "template_len": 0})
            result = {"results": results, "total": total, "succeeded": succeeded,
                      "mode": params.get("mode", "standard"), "params": params,
                      "locate_mode": True}
        else:
            result = design_pipeline(template, db_prefix, target, params,
                                     on_log=on_log, cancel=cancel,
                                     spec_db=spec_db, on_progress=on_progress,
                                     kmer_caches=kmer_caches)
        on_log(f"设计结束: {result['succeeded']}/{result['total']} 条查询成功")
        return result

    task = task_manager.start("design", f"引物设计({params.get('mode', 'standard')})", worker)
    return {"task_id": task.id}


# ---------------------------------------------------------------- 任务与 SSE (Tasks & SSE)

@app.get("/api/tasks")
async def api_tasks() -> dict:
    return {"tasks": [t.to_snapshot(with_logs=True) for t in task_manager.list()]}


@app.get("/api/tasks/events")
async def api_task_events(request: Request):
    queue: asyncio.Queue = asyncio.Queue()

    async def gen():
        try:
            yield _sse({"type": "snapshot",
                        "tasks": [t.to_snapshot(with_logs=True) for t in task_manager.list()]})
            while True:
                if await request.is_disconnected():
                    break
                try:
                    ev = await asyncio.wait_for(queue.get(), timeout=20)
                    yield _sse(ev)
                except asyncio.TimeoutError:
                    yield _sse({"type": "ping"})
        finally:
            task_manager.unsubscribe(queue)

    task_manager.subscribe(queue, asyncio.get_running_loop())
    return StreamingResponse(gen(), media_type="text/event-stream",
                             headers={"Cache-Control": "no-cache",
                                      "X-Accel-Buffering": "no"})


@app.get("/api/tasks/{task_id}/result")
async def api_task_result(task_id: str) -> dict:
    task = task_manager.get(task_id)
    if task is None:
        raise HTTPException(404, "任务不存在")
    if task.status == STATUS_FAILED:
        return {"status": "failed", "error": task.error or "任务失败"}
    if task.status != "succeeded":
        return {"status": task.status, "progress": task.progress,
                "progress_label": task.progress_label}
    return {"status": "succeeded", "result": task.result}


@app.post("/api/tasks/{task_id}/cancel")
async def api_task_cancel(task_id: str) -> dict:
    ok = task_manager.cancel(task_id)
    return {"cancelled": ok}


# ---------------------------------------------------------------- 导出 (Export)

def _export_response(content: str, filename: str, media_type: str) -> Response:
    return Response(content, media_type=media_type, headers={
        "Content-Disposition": f'attachment; filename="{filename}"'})


def _csv_response(rows: list[list[Any]]) -> Response:
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerows(rows)
    return _export_response("﻿" + buf.getvalue(), "stats.csv",
                            "text/csv; charset=utf-8")


@app.get("/api/export/blast/{task_id}/raw")
async def api_export_raw(task_id: str) -> Response:
    r = _blast_result(task_id)
    return _export_response(r.get("raw_output") or "", "blast_raw.txt", "text/plain")


@app.get("/api/export/blast/{task_id}/aln")
async def api_export_aln(task_id: str) -> Response:
    r = _blast_result(task_id)
    lines = ["BlastPrime Studio 比对结果", "=" * 40, ""]
    for q in r.get("parsed", []):
        lines.append(f"> {q['name']}  (length={q.get('length')})")
        for subj in q.get("subjects", []):
            lines.append(f"  {subj['name']}  (length={subj.get('length')})")
            for h in subj.get("hits", []):
                lines.append(f"    score={h.get('score')}  evalue={h.get('evalue')}  "
                             f"identity={h.get('identity')}  "
                             f"q:{h.get('qstart')}-{h.get('qend')}  "
                             f"s:{h.get('sstart')}-{h.get('send')}")
        lines.append("")
    return _export_response("\n".join(lines), "blast_alignment.txt", "text/plain")


@app.get("/api/export/blast/{task_id}/stats.csv")
async def api_export_stats(task_id: str) -> Response:
    r = _blast_result(task_id)
    rows = [["Query_ID", "Subject_ID", "Identity", "E-value", "Bit_Score",
             "Subject_Length", "Query_Start", "Query_End"]]
    for q in r.get("parsed", []):
        for subj in q.get("subjects", []):
            for h in subj.get("hits", []):
                rows.append([q["name"], subj.get("name", ""), h.get("identity", ""),
                             h.get("evalue", ""), h.get("bitscore", h.get("score", "")),
                             subj.get("length", ""), h.get("qstart", ""), h.get("qend", "")])
    return _csv_response(rows)


@app.get("/api/export/blast/{task_id}/hits.fasta")
async def api_export_hits(task_id: str) -> Response:
    r = _blast_result(task_id)
    parts: list[str] = []
    for q in r.get("parsed", []):
        hsp_no = 0
        for subj in q.get("subjects", []):
            for h in subj.get("hits", []):
                hsp_no += 1
                seq = h.get("sseq") or ""
                parts.append(f">subject_id_{subj.get('name', '')}_hsp{hsp_no} "
                             f"[from_query:{q['name']} score:{h.get('score', '')}]")
                for i in range(0, len(seq), 60):
                    parts.append(seq[i:i + 60])
    return _export_response("\n".join(parts) + "\n", "blast_hits.fasta",
                            "text/plain")


def _locate_segment(res: dict) -> str:
    """定位模式结果段的 # 注释行(§10:目标条目/坐标/链/提取范围/模板长度)。

    The # comment line for a locate-mode result segment (§10: target entry /
    coordinates / strand / extract range / template length).
    """
    loc = res.get("locate")
    if not loc:
        return ""
    g1, g2 = loc.get("target_genomic", ["", ""])
    return (f"# target={loc.get('entry', '')} genomic={g1}-{g2} "
            f"strand={loc.get('strand', '')} extract={loc.get('extract_start')}-"
            f"{loc.get('extract_end')} template_len={res.get('template_len', '')}")


@app.get("/api/export/design/{task_id}/pairs.fasta")
async def api_export_design_fasta(task_id: str) -> Response:
    r = _design_result(task_id)
    mode = r.get("mode", "standard")
    parts: list[str] = []
    for res in r.get("results", []):
        seg = _locate_segment(res)
        if seg:
            parts.append(seg)
        qname = res.get("query", "query")
        for i, p in enumerate(res.get("pairs", []), 1):
            left = p.get("left", {})
            seq = left.get("seq", "")
            covers = "" if p.get("covers_target") is None else \
                "covers=yes" if p.get("covers_target") else "covers=no"
            if mode == "standard":
                right = p.get("right", {})
                parts.append(f">Q_{qname}_Pair_{i}_F len={len(seq)} tm={left.get('tm', '')} "
                             f"gc={left.get('gc', '')}% abs={left.get('abs_start', '')}-"
                             f"{left.get('abs_end', '')} {covers}".rstrip())
                parts += [seq[j:j + 60] for j in range(0, len(seq), 60)]
                rseq = right.get("seq", "")
                parts.append(f">Q_{qname}_Pair_{i}_R len={len(rseq)} tm={right.get('tm', '')} "
                             f"gc={right.get('gc', '')}% abs={right.get('abs_start', '')}-"
                             f"{right.get('abs_end', '')}".rstrip())
                parts += [rseq[j:j + 60] for j in range(0, len(rseq), 60)]
            elif mode == "sgrna":
                parts.append(f">Q_{qname}_Guide_{i} len={len(seq)} gc={left.get('gc', '')}% "
                             f"PAM={p.get('pam', '')} strand={p.get('strand', '')}")
                parts += [seq[j:j + 60] for j in range(0, len(seq), 60)]
            else:  # single
                parts.append(f">Q_{qname}_Primer_{i} len={len(seq)} tm={left.get('tm', '')} "
                             f"gc={left.get('gc', '')}% abs={left.get('abs_start', '')}-"
                             f"{left.get('abs_end', '')}".rstrip())
                parts += [seq[j:j + 60] for j in range(0, len(seq), 60)]
    return _export_response("\n".join(parts) + "\n", "primers.fasta", "text/plain")


_LOCATE_CSV_COLS = ["Target_Entry", "Target_Genomic_Coords", "Strand",
                    "Extract_Range", "Template_Len", "Covers_Target"]


def _locate_csv_row(res: dict, p: dict) -> list[Any]:
    """定位模式 CSV 追加列(§10):目标条目/基因组坐标/链/提取范围/模板长度/覆盖。

    Extra CSV columns for locate mode (§10): target entry / genomic
    coordinates / strand / extract range / template length / coverage.
    """
    loc = res.get("locate")
    if not loc:
        return []
    g1, g2 = loc.get("target_genomic", ["", ""])
    cov = "" if p.get("covers_target") is None else \
        ("yes" if p.get("covers_target") else "no")
    return [loc.get("entry", ""), f"{g1}-{g2}", loc.get("strand", ""),
            f"{loc.get('extract_start')}-{loc.get('extract_end')}",
            res.get("template_len", ""), cov]


@app.get("/api/export/design/{task_id}/stats.csv")
async def api_export_design_csv(task_id: str) -> Response:
    r = _design_result(task_id)
    mode = r.get("mode", "standard")
    located = any(res.get("locate") for res in r.get("results", []))
    rows: list[list[Any]] = []
    if mode == "sgrna":
        rows.append(["Query_ID", "Guide_Index", "sgRNA_Seq", "PAM", "Strand",
                     "GC", "Tm", "Off_Target_Sites"])
    elif mode == "single":
        rows.append(["Query_ID", "Primer_Index", "Direction", "Primer_Seq",
                     "Primer_Len", "Tm", "GC", "Off_Target_Sites",
                     "Specificity_Level", "Composite_Score"])
    else:
        rows.append(["Query_ID", "Pair_Index", "Forward_Seq", "Forward_Len",
                     "Forward_Tm", "Forward_GC", "Forward_Start", "Forward_End",
                     "Reverse_Seq", "Reverse_Len", "Reverse_Tm", "Reverse_GC",
                     "Reverse_Start", "Reverse_End", "Product_Len",
                     "Max_Dimer_Consec", "Max_Dimer_Total",
                     "Off_Target_Sites", "Specificity_Level", "Composite_Score"])
    if located:
        rows[0] += _LOCATE_CSV_COLS
    for res in r.get("results", []):
        qname = res.get("query", "query")
        for i, p in enumerate(res.get("pairs", []), 1):
            spec = p.get("specificity", {})
            extra = _locate_csv_row(res, p)
            if mode == "sgrna":
                l = p.get("left", {})
                rows.append([qname, i, l.get("seq", ""), p.get("pam", ""),
                             p.get("strand", ""), l.get("gc", ""), l.get("tm", ""),
                             spec.get("off_target_sites", 0)] + extra)
                continue
            if mode == "single":
                l = p.get("left", {})
                rows.append([qname, i, "F", l.get("seq", ""), l.get("len", ""),
                             l.get("tm", ""), l.get("gc", ""),
                             spec.get("off_target_sites", 0),
                             spec.get("label", ""), p.get("composite_score", "")] + extra)
                continue
            l, rt = p.get("left", {}), p.get("right", {})
            dimer = p.get("dimer", {})
            rows.append([qname, i, l.get("seq", ""), l.get("len", ""), l.get("tm", ""),
                         l.get("gc", ""), l.get("start", ""),
                         (l.get("start") or 0) + (l.get("len") or 0) - 1,
                         rt.get("seq", ""), rt.get("len", ""), rt.get("tm", ""),
                         rt.get("gc", ""), rt.get("start", ""),
                         (rt.get("start") or 0) + (rt.get("len") or 0) - 1,
                         p.get("product_len", ""), dimer.get("max_consec", ""),
                         dimer.get("max_total", ""), spec.get("off_target_sites", 0),
                         spec.get("label", ""), p.get("composite_score", "")] + extra)
    return _csv_response(rows)


def _blast_result(task_id: str) -> dict:
    task = task_manager.get(task_id)
    if task is None or task.status != "succeeded" or task.kind != "blast":
        raise HTTPException(404, "结果不存在或任务未完成")
    return task.result or {}


def _design_result(task_id: str) -> dict:
    task = task_manager.get(task_id)
    if task is None or task.status != "succeeded" or task.kind != "design":
        raise HTTPException(404, "结果不存在或任务未完成")
    return task.result or {}


# ---------------------------------------------------------------- 项目保存/加载 (Project save/load)

@app.post("/api/project/save")
async def api_project_save(body: dict) -> dict:
    kind = body.get("kind")
    if kind not in ("blast", "primer_design"):
        raise HTTPException(400, "kind 必须为 blast 或 primer_design")
    payload = {
        "app": "BlastPrimeStudio", "version": 1,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "kind": kind,
        "raw_output": body.get("raw_output", ""),
        "parsed": body.get("parsed", []),
        "db_prefix": body.get("db_prefix", ""),
        "db_type": body.get("db_type", ""),
        "options": body.get("options", {}),
    }
    if kind == "primer_design":
        payload["primer_design"] = body.get("primer_design", {})
    return payload


@app.post("/api/project/load")
async def api_project_load(body: dict) -> dict:
    content = body.get("content")
    if isinstance(content, str):
        try:
            data = json.loads(content)
        except json.JSONDecodeError:
            raise HTTPException(400, "项目文件损坏:不是有效 JSON")
    elif isinstance(content, dict):
        data = content
    else:
        raise HTTPException(400, "缺少项目内容")
    if data.get("app") != "BlastPrimeStudio":
        raise HTTPException(400, "不是 BlastPrimeStudio 项目文件")
    if data.get("version") != 1:
        raise HTTPException(400, f"不支持的项目版本: {data.get('version')}")
    return data


# ---------------------------------------------------------------- 序列工具 (Sequence tools)

@app.post("/api/seq/revcomp")
async def api_seq_revcomp(body: dict) -> dict:
    seq = (body.get("seq") or "").strip().upper()
    return {"revcomp": blast.revcomp(seq)}


@app.post("/api/seq/translate")
async def api_seq_translate(body: dict) -> dict:
    seq = (body.get("seq") or "").strip().upper()
    code = (body.get("code") or "Standard").strip()
    frames = blast.six_frame_translate(seq, code)
    return {"frames": frames}


# ---------------------------------------------------------------- CLI 入口 (CLI entry)

# CLI 帮助文本(按系统语言选择 zh/en)
# CLI help texts (zh/en chosen by system language)
_CLI_LANG_TEXT = {
    "zh": {
        "desc": "BlastPrime Studio 本地生物信息学工作站",
        "host": "监听地址(默认 127.0.0.1)",
        "port": "监听端口(默认 8686)",
        "no_browser": "不自动打开浏览器",
        "loglevel": "日志级别 DEBUG/INFO/WARNING/ERROR",
        "logfile": "日志文件路径",
        "config": "配置文件路径(默认 数据目录/config.json)",
        "started": "BlastPrime Studio 启动: %s  (Ctrl+C 退出)",
        "open_fail": "自动打开浏览器失败(可手动访问 %s)",
    },
    "en": {
        "desc": "BlastPrime Studio — local bioinformatics workbench",
        "host": "listen address (default 127.0.0.1)",
        "port": "listen port (default 8686)",
        "no_browser": "do not auto-open a browser",
        "loglevel": "log level DEBUG/INFO/WARNING/ERROR",
        "logfile": "log file path",
        "config": "config file path (default: data dir/config.json)",
        "started": "BlastPrime Studio started: %s  (Ctrl+C to quit)",
        "open_fail": "failed to open browser (visit %s manually)",
    },
}


def _pick_port(host: str, port: int, max_tries: int = 20) -> int:
    """尝试绑定端口;失败(占用 / Windows 端口排除范围等权限错误)则
    顺延到下一端口,返回最终可用端口。全部失败返回原端口(交给
    uvicorn 报原始错误,便于诊断)。

    Probes the port and falls forward on bind failure (in use, or Windows
    excluded-port-range permission errors), returning the first free port;
    falls back to the requested one after max_tries so uvicorn surfaces the
    original error for diagnosis.
    """
    import socket
    for p in range(port, port + max_tries):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
                s.bind((host, p))
            return p
        except OSError:
            continue
    return port


def main() -> None:
    import argparse

    lang = detect_system_lang()
    H = _CLI_LANG_TEXT[lang]
    parser = argparse.ArgumentParser(prog="blastprime", description=H["desc"])
    parser.add_argument("--host", default="127.0.0.1", help=H["host"])
    parser.add_argument("--port", type=int, default=8686, help=H["port"])
    parser.add_argument("--no-browser", action="store_true", help=H["no_browser"])
    parser.add_argument("--loglevel", default=None, help=H["loglevel"])
    parser.add_argument("--logfile", default=None, help=H["logfile"])
    parser.add_argument("--config", default=None, help=H["config"])
    args = parser.parse_args()

    cfg = get_config()
    if args.config:
        from .config import Config
        cfg = Config(args.config)
    if args.loglevel:
        cfg.set("loglevel", args.loglevel)
    if args.logfile:
        cfg.set("logfile", args.logfile)

    logging.basicConfig(
        level=getattr(logging, str(cfg.get("loglevel", "INFO")).upper(), logging.INFO),
        filename=cfg.get("logfile") or None,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    # 端口顺延(帮助文本已声明"若被占用则自动顺延"):先探测再启动,
    # 避免 Windows 端口排除范围/占用直接抛 WinError 10013/10048 退出。
    # Resolve the port first (probe + fall forward), so excluded ranges /
    # in-use ports on Windows do not abort startup with 10013/10048.
    port = _pick_port(args.host, args.port)
    url = f"http://{args.host}:{port}/"
    if not args.no_browser:
        def _open() -> None:
            try:
                webbrowser.open(url)
            except Exception:  # noqa: BLE001
                log.warning(H["open_fail"], url)
        threading.Timer(1.2, _open).start()

    import uvicorn
    print(H["started"] % url)
    uvicorn.run(app, host=args.host, port=port, log_level="info")


if __name__ == "__main__":
    main()
