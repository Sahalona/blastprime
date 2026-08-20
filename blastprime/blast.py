"""BLAST 程序定位、建库、检索与结果解析。

- 工具定位顺序:打包模式 sys._MEIPASS/bin/ → exe 同级 bin/;源码模式:项目根目录 → 根目录 bin/ → 系统 PATH;再叠加手动指定的 blast_bin_dir。
- 结果解析:支持 -outfmt 0 成对文本与 tabular(含 qseq/sseq/sstrand)。

BLAST tool location, database construction, search, and result parsing.

- Tool lookup order: packaged mode sys._MEIPASS/bin/ → bin/ next to the exe; source mode: project root → root bin/ → system PATH; the manually specified blast_bin_dir is prepended.
- Result parsing: supports -outfmt 0 pairwise text and tabular (including qseq/sseq/sstrand).
"""

from __future__ import annotations

import logging
import os
import queue
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from Bio import SeqIO
from Bio.Seq import Seq

from .config import get_config

log = logging.getLogger("blastprime")

TOOLS = ["makeblastdb", "blastn", "blastp", "blastx", "tblastn", "tblastx", "blastdbcmd"]

_tool_cache: dict[str, str | None] = {}
_cache_stamp: int = 0


class BlastError(Exception):
    """BLAST 相关错误,message 面向用户展示。

    BlastError for BLAST-related errors; message is shown to the user.
    """


# ---------------------------------------------------------------- 程序定位
# ---------------------------------------------------------------- Tool location

def _candidates() -> list[Path]:
    """返回候选 bin 目录(优先级从高到低)。

    Return candidate bin directories (highest priority first).
    """
    dirs: list[Path] = []
    if getattr(sys, "frozen", False):  # PyInstaller
        meipass = getattr(sys, "_MEIPASS", None)
        if meipass:
            dirs.append(Path(meipass) / "bin")
        dirs.append(Path(sys.executable).parent / "bin")
    # 源码模式
    # Source mode
    root = Path(__file__).resolve().parent.parent
    dirs.append(root)
    dirs.append(root / "bin")
    manual = get_config().get("blast_bin_dir", "")
    if manual:
        dirs.insert(0, Path(manual))
    return dirs


def find_tool(name: str, refresh: bool = False) -> str | None:
    """按固定顺序查找 BLAST 可执行文件,结果内存缓存。

    Locate a BLAST executable in the fixed order, caching the result in memory.
    """
    global _tool_cache, _cache_stamp
    if name not in _tool_cache or refresh:
        path: str | None = None
        for d in _candidates():
            exe = d / (name + (".exe" if os.name == "nt" else ""))
            if exe.is_file() and os.access(exe, os.X_OK):
                path = str(exe)
                break
        if path is None:
            path = shutil.which(name)
        _tool_cache[name] = path
    return _tool_cache.get(name)


def detect_env() -> dict:
    """环境检测:返回 {tool: 路径或 None}。

    Environment detection: return {tool: path or None}.
    """
    result = {}
    for t in TOOLS:
        result[t] = find_tool(t, refresh=True)
    result["missing"] = [t for t in TOOLS if result[t] is None]
    return result


def set_manual_bin_dir(path: str | None) -> dict:
    """手动指定 BLAST 目录,立即重新检测并持久化。

    Set the BLAST directory manually, re-detect immediately, and persist the setting.
    """
    cfg = get_config()
    cfg.set("blast_bin_dir", (path or "").strip())
    _tool_cache.clear()
    env = detect_env()
    if env["missing"]:
        raise BlastError("指定目录中未找到完整的 BLAST 程序: " + ", ".join(env["missing"]))
    return env


def require_tool(name: str) -> str:
    path = find_tool(name)
    if not path:
        raise BlastError(f"未找到 BLAST 程序 {name},请先安装或手动指定 BLAST 目录")
    return path


# ---------------------------------------------------------------- 子进程运行
# ---------------------------------------------------------------- Subprocess execution

@dataclass
class CancelFlag:
    cancelled: bool = False


def run_proc(
    cmd: list[str],
    on_log: Callable[[str, str], None] | None = None,
    cancel: CancelFlag | None = None,
    timeout: int = 300,
    cwd: str | None = None,
    log_stdout: bool = True,
) -> subprocess.CompletedProcess:
    """运行子进程,流式转发 stdout 到 on_log(msg, level)(log_stdout=False 时仅收集),
    支持取消与超时。分析型 BLAST(tabular 含全序列列)应关掉 stdout 日志。
    启动命令与失败详情(exit code + 完整输出)同时写入日志文件,失败时还把
    输出尾部以 error 级推入任务日志 —— 建库/比对失败不再无声无息。

    Run a subprocess, streaming stdout to on_log(msg, level) (only collected
    when log_stdout=False); supports cancellation and timeout. Analytical
    BLAST (tabular with full-sequence columns) should disable stdout logging.
    The command line and failure detail (exit code + full output) are written
    to the log file; on failure the output tail is pushed to the task log at
    error level — so database-build/alignment failures are never silent.
    """
    # 命令行走双链路:日志文件(log.info)与网页任务日志(on_log)各一份。
    # 注意两者是独立的通道 —— 只写 logger 的话日志抽屉看不到命令。
    # The command line goes through both channels: the log file (log.info)
    # and the web task log (on_log). They are independent — logger-only
    # would leave the task drawer without the command.
    log.info("运行命令: %s", " ".join(cmd))
    if on_log:
        on_log(f"运行: {' '.join(cmd)}")
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        stdin=subprocess.DEVNULL,
        text=True,
        cwd=cwd,
        bufsize=1,
    )
    lines: list[str] = []
    deadline = timeout if timeout > 0 else 0
    # 取消/超时轮询:BLAST 长查询的搜索阶段可能长时间无 stdout 输出,
    # 若在主线程按行阻塞读取,取消与超时都无法生效。stdout 由读线程
    # 转发到队列,主循环限时取队列 —— 跨平台一致(Windows 的 select 只
    # 支持套接字,对管道调用会抛 WinError 10038,exe 上建库/比对全崩)。
    # Cancel/timeout polling: a long BLAST search can produce no stdout for
    # a while, so a blocking line loop in the main thread would make
    # cancellation/timeout ineffective. A reader thread forwards stdout into
    # a queue and the main loop polls it with a timeout — portable across
    # platforms (Windows select() supports sockets only; calling it on a
    # pipe raises WinError 10038, crashing every build/alignment on the exe).
    q: queue.Queue = queue.Queue()

    def _reader() -> None:
        try:
            assert proc.stdout is not None
            for ln in proc.stdout:
                q.put(ln.rstrip("\n"))
        finally:
            q.put(None)   # EOF 哨兵:正常结束与异常都必须送,否则主循环空转 (EOF sentinel; sent on both normal exit and errors, otherwise the main loop spins)

    threading.Thread(target=_reader, daemon=True,
                     name="proc-stdout-reader").start()
    try:
        start = time.monotonic()
        while True:
            if cancel is not None and cancel.cancelled:
                proc.kill()
                proc.wait()
                log.warning("命令已取消: %s ...", " ".join(cmd[:4]))
                raise BlastError("任务已取消")
            if deadline and time.monotonic() - start > deadline:
                proc.kill()
                proc.wait()
                log.warning("命令超时(%ss),已终止: %s ...", timeout, " ".join(cmd[:4]))
                raise BlastError(f"命令超时(>{timeout}s),已终止: {' '.join(cmd[:4])} ...")
            try:
                line = q.get(timeout=0.2)
            except queue.Empty:
                continue
            if line is None:
                break
            lines.append(line)
            if on_log and log_stdout:
                on_log(line)
        proc.wait(timeout=deadline if deadline else None)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()
        log.warning("命令超时(%ss),已终止: %s ...", timeout, " ".join(cmd[:4]))
        raise BlastError(f"命令超时(>{timeout}s),已终止: {' '.join(cmd[:4])} ...")
    if proc.returncode != 0:
        tail = "\n".join(lines[-40:])
        summary = f"命令失败 (exit {proc.returncode}): {' '.join(cmd[:4])} ..."
        # 完整输出(可能数万行)写日志文件,尾部 40 行同时进任务日志(ERROR 级,
        # 前端日志抽屉可见)—— stderr 已被合并进 stdout(subprocess.STDOUT)
        # Full output (potentially tens of thousands of lines) goes to the log
        # file; the last 40 lines also go to the task log at ERROR level (visible
        # in the frontend log drawer) — stderr is already merged into stdout
        log.error("%s\n--- 完整输出(共 %d 行) ---\n%s",
                  summary, len(lines), "\n".join(lines))
        if on_log:
            on_log(summary, "error")
            for ln in lines[-40:]:
                on_log(ln, "error")
        raise BlastError(f"{summary}\n{tail}")
    return subprocess.CompletedProcess(cmd, 0, "\n".join(lines), "")


# ---------------------------------------------------------------- 建库
# ---------------------------------------------------------------- Database construction

def merge_fasta_files(paths: list[str], out_fa: Path) -> int:
    """合并多个 FASTA 到临时文件(不修改原文件),返回记录数。

    Merge multiple FASTA files into a temporary file (originals untouched), returning the record count.
    """
    count = 0
    used: set[str] = set()
    with open(out_fa, "w", encoding="utf-8") as fout:
        for p in paths:
            try:
                recs = list(SeqIO.parse(p, "fasta"))
            except Exception as e:
                raise BlastError(f"无法解析 FASTA 文件 {p}: {e}")
            if not recs:
                raise BlastError(f"文件 {p} 中没有序列记录")
            for r in recs:
                fout.write(f">{_dedupe_id(used, r.description)}\n")
                seq = str(r.seq).replace(" ", "").replace("\t", "")
                for i in range(0, len(seq), 80):
                    fout.write(seq[i:i + 80] + "\n")
                count += 1
    return count


def detect_seq_type(paths: list[str]) -> str:
    """自动检测:取前 2000 碱基样本,ACGTUN 占比 ≥90% 判核酸,否则蛋白。

    Auto-detection: sample the first 2000 bases; if ACGTUN makes up ≥90%, classify as nucleotide, otherwise protein.
    """
    total = 0
    acgtun = 0
    for p in paths:
        for rec in SeqIO.parse(p, "fasta"):
            s = str(rec.seq)[:2000]
            total += len(s)
            acgtun += sum(1 for c in s.upper() if c in "ACGTUN")
            if total >= 2000:
                break
        if total >= 2000:
            break
    if total == 0:
        raise BlastError("无法读取序列内容")
    if acgtun / total >= 0.90:
        return "nucl"
    return "prot"


def build_database(
    fasta_paths: list[str],
    dbtype: str,          # auto | nucl | prot
    name: str,
    out_dir: str | None,
    on_log: Callable[[str, str], None] | None = None,
    on_progress: Callable[[float, str], None] | None = None,
    cancel: CancelFlag | None = None,
) -> dict:
    """构建 BLAST 数据库,返回 {prefix, type, seq_count, dir}。

    Build a BLAST database, returning {prefix, type, seq_count, dir}.
    """
    from .database import validate_db_name, default_db_dir, add_record

    name = validate_db_name(name)
    if dbtype == "auto":
        dbtype = detect_seq_type(fasta_paths)
        if on_log:
            on_log(f"[auto-detect] 序列类型: {'核酸' if dbtype == 'nucl' else '蛋白'}")

    out_dir_path = Path(out_dir) if out_dir and out_dir.strip() else default_db_dir()
    out_dir_path.mkdir(parents=True, exist_ok=True)
    prefix = out_dir_path / name

    with tempfile.TemporaryDirectory(prefix="blastprime_build_") as td:
        merged = Path(td) / "merged.fa"
        if on_progress:
            on_progress(0, "合并 FASTA 文件")
        n = merge_fasta_files(fasta_paths, merged)
        if on_log:
            on_log(f"已合并 {len(fasta_paths)} 个文件,共 {n} 条记录")

        mkb = require_tool("makeblastdb")
        cmd = [mkb, "-in", str(merged), "-dbtype", dbtype, "-out", str(prefix),
               "-parse_seqids", "-title", name]
        # makeblastdb 无逐条进度,阶段 label 让状态栏感知仍在推进
        # makeblastdb reports no per-record progress; the stage label keeps the status bar aware that work is still in progress
        if on_progress:
            on_progress(5, "makeblastdb 构建索引中")
        run_proc(cmd, on_log=on_log, cancel=cancel, timeout=300)
        if on_log:
            on_log(f"数据库构建完成: {prefix}")

    add_record(str(prefix), is_created=True)
    if on_progress:
        on_progress(100, "数据库构建完成")
    return {"prefix": str(prefix), "type": dbtype, "seq_count": n, "dir": str(out_dir_path)}


# ---------------------------------------------------------------- blastdbcmd

def blastdbcmd_info(prefix: str) -> dict:
    """库详情:blastdbcmd -info 输出解析。

    Database details: parse the output of blastdbcmd -info.
    """
    tool = require_tool("blastdbcmd")
    p = Path(prefix)
    if p.is_file() and p.suffix.lstrip(".") in ("nin", "pin", "nsq", "psq"):
        p = p.with_suffix("")
    if not (Path(str(p) + ".nin").exists() or Path(str(p) + ".pin").exists()):
        raise BlastError(f"数据库不存在或索引缺失: {p}")
    try:
        r = run_proc([tool, "-db", str(p), "-info"], timeout=120)
        text = r.stdout
        info: dict = {"raw": text}
        # 兼容两种输出:旧版行首无 #,新版(BLAST 2.14+)带 #;
        # 数字可能带千位逗号(如 "460,334,017 total bases"),取值后剥离再转 int
        # Compatible with both output styles: older versions have no leading #, newer ones (BLAST 2.14+) do;
        # numbers may contain thousands separators (e.g. "460,334,017 total bases"); strip them before converting to int
        m = re.search(r"^[#\s]*([\d,]+)\s+sequences?\s*;", text, re.MULTILINE)
        if m:
            info["seq_count"] = int(m.group(1).replace(",", ""))
        m = re.search(r"([\d,]+)\s+total (?:bases|letters)", text)
        if m:
            info["total_len"] = int(m.group(1).replace(",", ""))
    except BlastError:
        # 索引文件不完整(如手动复制缺 .ndb 等)时降级:仅报告类型,序列数留空
        # Degrade gracefully when index files are incomplete (e.g. a manually copied db missing .ndb): report type only, leave sequence count blank
        info = {"raw": "", "seq_count": None, "total_len": None}
    info["type"] = "nucl" if Path(str(p) + ".nin").exists() else "prot"
    return info


def list_entries(prefix: str) -> list[str]:
    """库条目名列表(blastdbcmd -entry all)。

    List database entry names (blastdbcmd -entry all).
    """
    tool = require_tool("blastdbcmd")
    r = run_proc([tool, "-db", str(prefix), "-entry", "all", "-outfmt", "%a"], timeout=300)
    return [ln.strip() for ln in r.stdout.splitlines() if ln.strip()]


def fetch_entry(prefix: str, entry: str) -> str:
    """取库条目序列(FASTA 文本)。

    Fetch an entry's sequence as FASTA text.
    """
    tool = require_tool("blastdbcmd")
    r = run_proc([tool, "-db", str(prefix), "-entry", entry, "-outfmt", "%f"], timeout=300)
    return r.stdout


def entry_length(prefix: str, entry: str) -> int:
    """库条目全长(blastdbcmd -outfmt %l)。侧翼提取边界判定用。

    Full length of a database entry (blastdbcmd -outfmt %l). Used for flanking-extraction boundary checks.
    """
    tool = require_tool("blastdbcmd")
    r = run_proc([tool, "-db", str(prefix), "-entry", entry, "-outfmt", "%l"], timeout=300)
    text = r.stdout.strip()
    try:
        return int(text)
    except ValueError:
        raise BlastError(
            f"无法获取库条目 '{entry}' 的长度(blastdbcmd 输出: {text[:80]})")


def fetch_region(prefix: str, entry: str, start: int, end: int) -> str:
    """取库条目指定区间(1-based 闭区间)的序列,FASTA 文本。

    严格截取于同一库条目内部(§5.1),调用方负责把坐标钳制在 [1, entry_length] 内。

    Fetch the sequence of a specified region (1-based closed interval) of a database entry as FASTA text.

    The region is strictly confined within a single database entry (§5.1); the caller is
    responsible for clamping coordinates to [1, entry_length].
    """
    tool = require_tool("blastdbcmd")
    r = run_proc([tool, "-db", str(prefix), "-entry", entry,
                  "-range", f"{start}-{end}", "-outfmt", "%f"], timeout=300)
    return r.stdout


# ---------------------------------------------------------------- 序列文本预处理
# ---------------------------------------------------------------- Sequence text preprocessing


def strip_comment_lines(text: str) -> str:
    """去掉 `#` 开头的注释行(序列输入框自动忽略,blast/design 共用)。

    逐行过滤:空行保留(SeqIO 容错),仅 lstrip 后以 `#` 起始的行删除。

    Strip lines starting with `#` (auto-ignored in sequence input boxes).
    """
    return "\n".join(l for l in text.splitlines() if not l.lstrip().startswith("#"))


def _dedupe_id(used: set[str], ident: str) -> str:
    """seen 集合中取最小未用后缀:"ident (1)"、"ident (2)"…

    Smallest unused " (N)" suffix for a duplicate identifier.
    """
    if ident not in used:
        used.add(ident)
        return ident
    n = 1
    while f"{ident} ({n})" in used:
        n += 1
    cand = f"{ident} ({n})"
    used.add(cand)
    return cand


def dedupe_fasta_headers(text: str) -> str:
    """多 FASTA 重名标题追加 " (1)"/" (2)" 后缀(仅提交/解析时改写)。

    按完整标题行去重,取最小未用后缀;含 `=` 的名称型行跳过(名称型逐行
    独立解析,不参与去重)。幂等:再跑一遍不再追加。

    Deduplicate FASTA headers by appending " (N)" suffixes (submit/parse
    time only, the input box keeps the original text). Name-type lines
    (headers containing `=`) are skipped — each resolves independently.
    Idempotent.
    """
    used: set[str] = set()
    out = []
    for line in text.splitlines():
        if line.startswith(">") and "=" not in line:
            out.append(">" + _dedupe_id(used, line[1:].strip()))
        else:
            out.append(line)
    return "\n".join(out)


# ---------------------------------------------------------------- 名称型 query
# ---------------------------------------------------------------- Named query

# 名称型参数:key 大小写不敏感,顺序自由;未知 key 忽略;除 gene 外全部可选
# Named-query params: case-insensitive keys, free order, unknown keys ignored;
# everything except the gene is optional.
_NAME_KEYS = {"range", "database", "name", "targetbase"}


def targetbase_disabled(v) -> bool:
    """targetbase=False(大小写不敏感)关闭特异性比对:设计页不做 blastn-short
    逆向验证;blast 页运行库回落页面库(与不写 targetbase 完全一致,无影响)。

    targetbase=False (case-insensitive) disables specificity evaluation: the
    design page skips the blastn-short reverse validation; the blast page
    falls back to the page db (identical to omitting targetbase — no effect).
    """
    return v is not None and str(v).strip().lower() == "false"


def parse_name_query(text: str) -> dict | None:
    """解析名称型 query,返回参数字典或 None。

    语法: `>gene[,range=100-200|100-*|*-200][,database=库][,name=显示名]
        [,targetbase=库]` —— 首 token(空格/逗号分隔)为 gene;其余逗号分
    token 逐个 `key=value`(key 大小写不敏感,未知 key 忽略,顺序自由)。
    除 gene 外均可省略。返回 {gene, range, database, name, targetbase}
    (缺省 None);首 token 为空或含 `=` 时返回 None(普通 FASTA)。

    Parse a named query into a params dict, or None (plain FASTA). The first
    token is the gene; remaining comma tokens are `key=value` pairs with
    case-insensitive, free-order, unknown-ignored keys. Returns
    {gene, range, database, name, targetbase} (missing → None).
    """
    stripped = text.strip().lstrip(">").strip()
    tokens = [t.strip() for t in stripped.split(",")]
    gene = (tokens[0].split() or [""])[0] if tokens else ""
    if not gene or "=" in gene:
        return None
    out = {"gene": gene, "range": None, "database": None,
           "name": None, "targetbase": None}
    for t in tokens[1:]:
        if "=" in t:
            k, _, v = t.partition("=")
            k, v = k.strip().lower(), v.strip()
            if k in _NAME_KEYS and v:
                out[k] = v
    return out


def _slice_range(seq: str, rng: str) -> tuple[str, int, int]:
    """range=起-止,支持 * 通配;1-based 闭区间。返回 (切片, 起, 止)。

    range=start-end, supporting * wildcards; 1-based closed interval. Returns (slice, start, end).
    """
    try:
        a, b = rng.split("-")
        a, b = a.strip(), b.strip()
        start = 1 if a == "*" else int(a)
        end = len(seq) if b == "*" else int(b)
    except ValueError:
        raise BlastError(f"范围格式无效: range={rng} (应为 起-止 或 *-* )")
    if start < 1 or end > len(seq) or start > end:
        raise BlastError(f"范围越界: range={rng} (序列长度 {len(seq)})")
    return seq[start - 1:end], start, end


def resolve_db_ref(ref: str) -> str:
    """按引用名解析历史库,返回库前缀。

    先对 prefix 的 basename 子串匹配(大小写不敏感),再对 note(备注名)
    子串匹配;取首个命中,找不到 → BlastError。

    Resolve a database reference (basename or note) against history records,
    returning the db prefix. Basename substring match first, then note
    substring match; first hit wins; BlastError when not found.
    """
    from .database import list_records
    ref_l = ref.lower()
    recs = list_records()
    for rec in recs:
        if ref_l in Path(rec["prefix"]).name.lower():
            return rec["prefix"]
    for rec in recs:
        note = (rec.get("note") or "").lower()
        if note and ref_l in note:
            return rec["prefix"]
    raise BlastError(f"历史记录中未找到数据库 '{ref}'")


def resolve_entry(db_prefix: str, name: str) -> str:
    """名称解析:精确 → 去 lcl| 前缀 → 模糊包含匹配(blastdbcmd -list)。
    返回该条目完整 FASTA(seq 全长,未切片)。

    Name resolution: exact → strip lcl| prefix → fuzzy substring match (blastdbcmd -list).
    Returns the entry's full FASTA (full-length seq, unsliced).
    """
    fa = None
    err = None
    for attempt in (name, name[4:] if name.lower().startswith("lcl|") else name):
        try:
            fa = fetch_entry(db_prefix, attempt)
            break
        except BlastError as e:
            err = e
    if fa is None:
        # 模糊包含匹配
        # Fuzzy substring match
        try:
            entries = list_entries(db_prefix)
        except BlastError as e:
            raise BlastError(f"无法列出库条目: {e}")
        hits = [e for e in entries if name.lower() in e.lower()]
        if not hits:
            raise BlastError(
                f"库中找不到该基因: '{name}' {('(索引不存在)' if err else '')}"
            )
        fa = fetch_entry(db_prefix, hits[0])
    return fa


def resolve_name_query(name_query: str, db_prefix: str | None = None) -> str:
    """解析名称型 query,返回回填后的 FASTA 文本。

    语法: >gene[,range=100-200|100-*|*-200][,database=库][,name=显示名]
        [,targetbase=库]
    - range 缺省 = 全长;database 缺省 = targetbase,再缺省 = 页面库
      (db_prefix);
    - 条目库(database 或 targetbase)与 db_prefix 独立:特异性比对/回填
      在条目库进行,blast 运行库由调用方另行决定(targetbase 同时作用
      于 blast 页);
    - 标题 = name || gene(显示名,避免把 ,range= 等语法留在查询名里)。

    Resolve a named query, returning the backfilled FASTA text. range
    defaults to the full length; database defaults to targetbase, then to
    the page db (db_prefix). The header carries the display name
    (name || gene).
    """
    q = parse_name_query(name_query)
    if not q:
        raise BlastError(f"无法解析名称型 query: {name_query}")
    entry_db = None
    if q["database"]:
        entry_db = resolve_db_ref(q["database"])
    elif q["targetbase"] and not targetbase_disabled(q["targetbase"]):
        # targetbase=False → 不解析为库,回落页面库(blast 页无影响)
        # targetbase=False -> never resolve it as a db; fall back to the page db
        entry_db = resolve_db_ref(q["targetbase"])
    if not entry_db:
        entry_db = db_prefix
    if not entry_db:
        raise BlastError("名称型 query 需要先选择本地数据库")

    fa = resolve_entry(entry_db, q["gene"])

    import io
    rec = SeqIO.read(io.StringIO(fa), "fasta")
    seq = str(rec.seq)
    if q["range"]:
        seq, start, end = _slice_range(seq, q["range"])
    display = q["name"] or q["gene"]
    rng = q["range"]
    header = f">{display}{f'[{rng}]' if rng else ''}"
    out = []
    for i in range(0, len(seq), 60):
        out.append(seq[i:i + 60])
    return header + "\n" + "\n".join(out) + "\n"


# ---------------------------------------------------------------- 比对执行
# ---------------------------------------------------------------- BLAST execution

def auto_program(db_type: str, query_type: str) -> str:
    """程序自动选择:库类型 × 查询类型。

    Automatic program selection: database type × query type.
    """
    if db_type == "nucl":
        return "blastn" if query_type != "prot" else "tblastn"
    return "blastx" if query_type != "prot" else "blastp"


def run_blast(
    db_prefix: str | None,
    query_fa: str,
    program: str,          # auto | blastn | ...
    query_type: str,       # auto | nucl | prot
    evalue: str,
    max_targets: int,
    matrix: str,
    extra_args: str,
    short_seq_mode: bool,
    remote: bool = False,
    remote_db: str = "nr",
    on_log: Callable[[str, str], None] | None = None,
    on_progress: Callable[[float, str], None] | None = None,
    cancel: CancelFlag | None = None,
    timeout: int = 600,
) -> dict:
    """执行 BLAST,返回 {raw_output, parsed, program, db, short_mode}。

    raw_output 为 -outfmt 0 成对文本;parsed 由 parse_outfmt0 解析。

    Run BLAST, returning {raw_output, parsed, program, db, short_mode}.

    raw_output is -outfmt 0 pairwise text; parsed is produced by parse_outfmt0.
    """
    if query_type == "auto":
        query_type = detect_seq_type([query_fa])
    if program == "auto":
        db_type = None
        if remote:
            db_type = "nucl" if query_type != "prot" else "prot"
        else:
            from .database import db_type_of
            db_type = db_type_of(db_prefix) if db_prefix else None
            if db_type is None:
                raise BlastError("无法确定数据库类型,请检查数据库文件")
        program = auto_program(db_type, query_type)

    tool = require_tool(program)
    cmd = [tool, "-query", query_fa]
    if remote:
        cmd += ["-db", remote_db, "-remote"]
    else:
        if not db_prefix:
            raise BlastError("请选择本地数据库")
        cmd += ["-db", db_prefix]
    cmd += ["-outfmt", "0", "-evalue", str(evalue), "-max_target_seqs", str(max_targets)]
    if short_seq_mode:
        if program in ("blastn", "blastx", "tblastn"):
            cmd += ["-task", "blastn-short", "-word_size", "7"]
        elif program == "blastp":
            cmd += ["-task", "blastp-short", "-word_size", "2"]
    if program in ("blastp", "blastx", "tblastn", "tblastx") and matrix:
        cmd += ["-matrix", matrix]
    if extra_args and extra_args.strip():
        cmd += extra_args.split()
    cmd += ["-out", "-"]  # stdout

    if on_progress:
        # BLAST 无逐条进度,阶段化推进(启动→子进程结束→解析→完成)让状态栏
        # 进度条定量移动而非恒 0;远程比对可能数小时,阶段标签同样可见
        # BLAST reports no per-record progress; staged advances (start → child
        # done → parse → done) keep the status-bar bar moving quantitatively
        # instead of stuck at 0; remote alignments may take hours, and the
        # stage label stays visible there too
        on_progress(0, f"{program} 比对运行中")
    if on_log:
        # 成对文本输出可达数万行,逐行进 SSE 会刷爆前端主线程导致页面无响应;
        # 只收集不转发(失败时仍随错误消息附最后 40 行),原始输出完成后经结果接口提供
        # Pairwise text output can reach tens of thousands of lines; forwarding line-by-line over SSE
        # would flood the frontend main thread and freeze the page. Collect without forwarding
        # (on failure, the last 40 lines are still attached to the error message); the raw output
        # is delivered via the results endpoint once complete.
        # (完整命令已由 run_proc 写入日志文件,并推入任务日志)
        # (run_proc already wrote the full command to the log file and pushed
        # it into the task log)
        on_log("比对进行中,完整原始输出将在完成后提供…")
    r = run_proc(cmd, on_log=on_log, cancel=cancel, timeout=timeout, log_stdout=False)
    raw = r.stdout
    if on_progress:
        on_progress(90, "解析比对结果")
    parsed = parse_outfmt0(raw)
    if on_progress:
        on_progress(100, "比对完成")
    return {
        "raw_output": raw,
        "parsed": parsed,
        "program": program,
        "db": remote_db if remote else db_prefix,
        "remote": remote,
        "short_mode": short_seq_mode,
    }


def run_blast_tabular(
    db_prefix: str,
    query_fa: str,
    program: str = "blastn",
    evalue: float = 10.0,
    max_targets: int = 5000,
    extra_args: list[str] | None = None,
    on_log: Callable[[str, str], None] | None = None,
    cancel: CancelFlag | None = None,
    timeout: int = 600,
) -> list[dict]:
    """tabular 输出(含 qseq/sseq/sstrand),供匹配深度与特异性分析。

    Tabular output (including qseq/sseq/sstrand), for match-depth and specificity analysis.
    """
    tool = require_tool(program)
    cmd = [tool, "-query", query_fa, "-db", db_prefix,
           "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore qseq sseq",
           "-evalue", str(evalue), "-max_target_seqs", str(max_targets)]
    if extra_args:
        cmd += extra_args
    # tabular 输出含完整 qseq/sseq(最长可达数千碱基),不进任务日志
    # Tabular output contains full qseq/sseq (up to thousands of bases), so keep it out of the task log
    r = run_proc(cmd, on_log=on_log, cancel=cancel, timeout=timeout, log_stdout=False)
    return parse_tabular(r.stdout)


# ---------------------------------------------------------------- 解析:-outfmt 0
# ---------------------------------------------------------------- Parsing: -outfmt 0

_SUBJECT_LINE_RE = re.compile(r"^>+\s*(?P<name>.+?)\s*$")
_SCORE_RE = re.compile(r"Score\s*=\s*(?P<score>[\d.]+)\s*bits?\s*\((?P<bits>\d+)\),\s*Expect\s*=\s*(?P<evalue>[\d.eE+-]+)")
_IDENT_RE = re.compile(r"Identities\s*=\s*(?P<n>\d+)/(?P<den>\d+)\s*\((?P<pct>[\d.]+)%\)")
_GAPS_RE = re.compile(r"Gaps\s*=\s*(?P<n>\d+)/(?P<den>\d+)\s*\((?P<pct>[\d.]+)%\)")
_STRAND_RE = re.compile(r"Strand=(?P<strand>\S+)")
_FRAME_RE = re.compile(r"Frame=(?P<frame>\S+)")
_QUERY_LINE_RE = re.compile(r"^Query\s+(?P<start>\d+)\s+(?P<seq>\S+)\s+(?P<end>\d+)")
_SBJCT_LINE_RE = re.compile(r"^Sbjct\s+(?P<start>\d+)\s+(?P<seq>\S+)\s+(?P<end>\d+)")


def parse_outfmt0(text: str) -> list[dict]:
    """解析 -outfmt 0 成对输出为 Query → Subject → HSP 结构。

    Parse -outfmt 0 pairwise output into a Query → Subject → HSP structure.
    """
    # 按 "Query= " 切分
    # Split by "Query= "
    sections = re.split(r"(?m)^Query=\s*", text)
    queries: list[dict] = []
    for sec in sections[1:]:
        lines = sec.splitlines()
        name = lines[0].strip() if lines else ""
        length = 0
        m = re.search(r"Length=(\d+)", sec)
        if m:
            length = int(m.group(1))
        if "No hits found" in sec:
            queries.append({"name": name, "length": length, "subjects": []})
            continue
        subjects: list[dict] = []
        # 逐块解析 "> subject" 块
        # Parse each "> subject" block
        blocks = re.split(r"(?m)^>", sec)[1:]
        for blk in blocks:
            blk_lines = blk.splitlines()
            subj_name = _SUBJECT_LINE_RE.match(blk_lines[0]).group("name") if blk_lines and _SUBJECT_LINE_RE.match(blk_lines[0]) else blk_lines[0].strip()
            subj_len = 0
            m = re.search(r"Length=(\d+)", blk)
            if m:
                subj_len = int(m.group(1))
            hits: list[dict] = []
            # HSP 单元:Score 行开始,到空行结束
            # An HSP unit starts at a Score line and ends at a blank line
            current: dict | None = None
            aln_rows: list[str] = []
            in_aln = False
            for ln in blk_lines[1:]:
                ms = _SCORE_RE.search(ln)
                if ms:
                    if current:
                        _finish_hsp(current, aln_rows)
                        hits.append(current)
                    current = {
                        "score": float(ms.group("score")),
                        "bits": int(ms.group("bits")),
                        "bitscore": int(ms.group("bits")),   # 与 bits 同值,前端契约名 (same value as bits; frontend contract name)
                        "evalue": float(ms.group("evalue")),
                        "qstart": None, "qend": None, "sstart": None, "send": None,
                        "strand": None, "frame": None, "identity": None,
                        "ident_frac": None, "gaps": None,
                        "qseq": "", "sseq": "", "mid": "",
                        "alignment": "",
                    }
                    aln_rows = []
                    in_aln = False
                    continue
                if current is None:
                    continue
                mi = _IDENT_RE.search(ln)
                if mi:
                    current["identity"] = f"{mi.group('n')}/{mi.group('den')} ({mi.group('pct')}%)"
                    current["ident_frac"] = float(mi.group("n")) / float(mi.group("den"))
                    continue
                mg = _GAPS_RE.search(ln)
                if mg:
                    current["gaps"] = f"{mg.group('n')}/{mg.group('den')} ({mg.group('pct')}%)"
                    continue
                mst = _STRAND_RE.search(ln)
                if mst:
                    current["strand"] = mst.group("strand")
                    continue
                mf = _FRAME_RE.search(ln)
                if mf:
                    current["frame"] = mf.group("frame")
                    continue
                if _QUERY_LINE_RE.match(ln) or _SBJCT_LINE_RE.match(ln) or (in_aln and (ln.startswith(" ") or ln.startswith("|") or ln.startswith("-"))):
                    in_aln = True
                    aln_rows.append(ln)
                    continue
                if not ln.strip():
                    in_aln = False
            if current:
                _finish_hsp(current, aln_rows)
                hits.append(current)
            if hits:
                for h in hits:
                    h.setdefault("sseqid", subj_name)  # 前端 hover 提示需要 (needed for the frontend hover tooltip)
                subjects.append({"name": subj_name, "length": subj_len, "hits": hits})
        queries.append({"name": name, "length": length, "subjects": subjects})
    return queries


def _finish_hsp(hsp: dict, aln_rows: list[str]) -> None:
    """从三行式比对行中提取坐标、序列与纯文本比对。

    注意:反向链比对中 Sbjct 坐标递减(首行坐标 > 末行坐标),
    基因组坐标 = (末行 end, 首行 start),且 sseq 需反向互补才是真实序列。

    Extract coordinates, sequences, and the plain-text alignment from three-line alignment rows.

    Note: in minus-strand alignments, Sbjct coordinates decrease (first-row coordinate >
    last-row coordinate); genomic coordinates = (last-row end, first-row start), and sseq
    must be reverse-complemented to obtain the true sequence.
    """
    qseq = ""
    sseq = ""
    mid = ""
    q_first = q_last = s_first = s_last = None
    q_row_len = s_row_len = 60  # blast 三行式固定列宽 (fixed column width of BLAST's three-line alignment rows)
    for ln in aln_rows:
        mq = _QUERY_LINE_RE.match(ln)
        ms = _SBJCT_LINE_RE.match(ln)
        if mq:
            if q_first is None:
                q_first = int(mq.group("start"))
            q_last = int(mq.group("end"))
            qseq += mq.group("seq")
            q_row_len = len(mq.group("seq"))
        elif ms:
            if s_first is None:
                s_first = int(ms.group("start"))
            s_last = int(ms.group("end"))
            sseq += ms.group("seq")
            s_row_len = len(ms.group("seq"))
        else:
            stripped = ln.strip()  # 中间行(匹配/错配/缺口符号) (middle line: match/mismatch/gap symbols)
            # 整行空格(完全错配区)strip 后为空,必须按列宽补空格,否则与 qseq/sseq 不等长导致前端按行切片错位
            # A fully mismatched row strips to empty; pad it to the row width, otherwise qseq/sseq
            # would differ in length and the frontend's per-row slicing would misalign
            mid += stripped if stripped else " " * max(q_row_len, s_row_len)
    hsp["qseq"] = qseq
    hsp["sseq"] = sseq
    hsp["mid"] = mid
    hsp["alignment"] = "\n".join(aln_rows)
    if q_first is not None:
        hsp["qstart"], hsp["qend"] = q_first, q_last
    if s_first is not None:
        minus = bool(hsp.get("strand")) and "Minus" in hsp["strand"]
        if minus:
            hsp["sstart"], hsp["send"] = s_last, s_first
        else:
            hsp["sstart"], hsp["send"] = s_first, s_last
        hsp["sbjct_genomic"] = str(Seq(sseq).reverse_complement()) if minus else sseq


# ---------------------------------------------------------------- 解析:tabular
# ---------------------------------------------------------------- Parsing: tabular

def parse_tabular(text: str) -> list[dict]:
    """解析自定义 tabular 输出(含 qseq/sseq/sstrand)。

    Parse custom tabular output (including qseq/sseq/sstrand).
    """
    hsps: list[dict] = []
    for ln in text.splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) < 15:
            f = ln.split()
        if len(f) < 15:
            continue
        try:
            hsps.append({
                "qseqid": f[0], "sseqid": f[1],
                "pident": float(f[2]), "length": int(f[3]),
                "mismatch": int(f[4]), "gapopen": int(f[5]),
                "qstart": int(f[6]), "qend": int(f[7]),
                "sstart": int(f[8]), "send": int(f[9]),
                "sstrand": f[10],
                "evalue": float(f[11]), "bitscore": float(f[12]),
                "qseq": f[13].upper(), "sseq": f[14].upper(),
            })
        except (ValueError, IndexError):
            continue
    return hsps


# ---------------------------------------------------------------- 序列工具
# ---------------------------------------------------------------- Sequence utilities

def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def clean_sequence(seq: str) -> str:
    """序列清洗:保留标题行,序列行仅保留字母、*、-。

    Sequence cleaning: keep header lines; in sequence lines keep only letters, * and -.
    """
    out = []
    for ln in seq.splitlines():
        if ln.startswith(">"):
            out.append(ln)
        else:
            out.append(re.sub(r"[^A-Za-z*\-]", "", ln))
    return "\n".join(out)


GENETIC_CODES = {
    "Standard": "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "Vertebrate Mitochondrial": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    "Yeast Mitochondrial": "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "Invertebrate Mitochondrial": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
}


def six_frame_translate(seq: str, code_name: str = "Standard") -> list[dict]:
    """6 框翻译:正链 3 框 + 负链 3 框。返回 [{frame, strand, seq}]。

    Six-frame translation: 3 frames on the forward strand + 3 on the reverse strand.
    Returns [{frame, strand, seq}].
    """
    from Bio.Seq import translate
    table = GENETIC_CODES.get(code_name, GENETIC_CODES["Standard"])
    frames = []
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        for offset in range(3):
            seg = s[offset:]
            prot = str(translate(seg, table=table, stop_symbol="*"))
            frames.append({"strand": strand, "frame": offset + 1, "seq": prot})
    return frames
