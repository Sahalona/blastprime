"""后台任务调度、取消与 SSE 进度推送。

所有长任务(BLAST 子进程、建库、引物设计)在后台线程执行,不阻塞 HTTP;
进度/日志经 SSE 实时推送到前端,支持断线重连(新连接先收到活动任务快照)。
取消任务 = 设置 CancelFlag → run_proc 杀子进程 → 上下文管理器清理临时文件。

Background task scheduling, cancellation, and SSE progress push.

All long-running tasks (BLAST subprocesses, database building, primer design)
run on background threads without blocking HTTP; progress/logs are pushed to the
frontend in real time over SSE, with reconnection support (new connections first
receive a snapshot of active tasks). Cancelling a task = setting CancelFlag →
run_proc kills the subprocess → the context manager cleans up temporary files.
"""

from __future__ import annotations

import asyncio
import logging
import threading
import time
import uuid
from dataclasses import dataclass, field
from typing import Any, Callable, Coroutine

from .blast import BlastError, CancelFlag

log = logging.getLogger("blastprime")

STATUS_PENDING = "pending"
STATUS_RUNNING = "running"
STATUS_SUCCEEDED = "succeeded"
STATUS_FAILED = "failed"
STATUS_CANCELLED = "cancelled"

MAX_TASK_HISTORY = 50


@dataclass
class Task:
    id: str
    kind: str                      # build_db / blast / analysis / design
    title: str                     # 显示标题 (Display title)
    status: str = STATUS_PENDING
    progress: float = 0.0          # 0~100
    progress_label: str = ""
    logs: list[dict] = field(default_factory=list)   # [{t, level, msg}]
    result: dict | None = None
    error: str | None = None
    created_at: float = field(default_factory=time.time)
    finished_at: float | None = None
    _flag: CancelFlag = field(default_factory=CancelFlag)
    _thread: threading.Thread | None = None

    def to_snapshot(self, with_logs: bool = True) -> dict:
        d = {
            "id": self.id, "kind": self.kind, "title": self.title,
            "status": self.status, "progress": round(self.progress, 1),
            "progress_label": self.progress_label,
            "created_at": self.created_at, "finished_at": self.finished_at,
            "error": self.error,
        }
        if with_logs:
            d["logs"] = self.logs
        return d


WorkerFn = Callable[[CancelFlag, Callable[[str, str], None],
                     Callable[[float, str], None]], Any]


class TaskManager:
    """线程安全的任务注册表 + SSE 事件发布。

    Thread-safe task registry + SSE event publishing.
    """

    def __init__(self) -> None:
        self._tasks: dict[str, Task] = {}
        self._lock = threading.Lock()
        self._queues: set[asyncio.Queue] = set()
        self._loop: asyncio.AbstractEventLoop | None = None

    # ------------------------------------------------------------ 任务控制
    # ------------------------------------------------------------ Task control

    def start(self, kind: str, title: str, worker: WorkerFn) -> Task:
        task = Task(id=uuid.uuid4().hex[:12], kind=kind, title=title)
        task.status = STATUS_RUNNING
        with self._lock:
            self._tasks[task.id] = task
            if len(self._tasks) > MAX_TASK_HISTORY:
                # 清理最旧的已完成任务
                # Prune the oldest finished tasks.
                finished = sorted(
                    (t for t in self._tasks.values() if t.finished_at),
                    key=lambda t: t.finished_at or 0)
                for t in finished[:len(self._tasks) - MAX_TASK_HISTORY]:
                    self._tasks.pop(t.id, None)

        def _run() -> None:
            log.info("任务开始 %s(%s): %s", task.id, kind, title)
            try:
                result = worker(task._flag,
                                lambda msg, level="info": self.log(task.id, msg, level),
                                lambda v, label="": self.progress(task.id, v, label))
                if task._flag.cancelled:
                    self._finish(task.id, status=STATUS_CANCELLED, error="任务已取消")
                else:
                    self._finish(task.id, status=STATUS_SUCCEEDED, result=result)
            except BlastError as e:
                status = STATUS_CANCELLED if task._flag.cancelled else STATUS_FAILED
                msg = str(e)
                if status == STATUS_FAILED:
                    # 失败详情(含 BLAST 子进程输出尾部)进任务日志(ERROR 级)与日志文件,
                    # 前端日志抽屉直接可见,不再只有任务状态栏里一行错误
                    # Failure detail (incl. the BLAST child-process output tail)
                    # goes to the task log (ERROR level) and the log file — the
                    # frontend log drawer shows it, not just one error line
                    self.log(task.id, f"任务失败: {msg}", "error")
                    log.error("任务失败 %s(%s): %s", task.id, kind, msg)
                self._finish(task.id, status=status, error=msg)
            except Exception as e:  # noqa: BLE001 — 兜底记录一切失败 (catch-all that records any failure)
                # traceback 写日志文件(未预期异常的唯一排查途径)
                # the traceback goes to the log file — the only way to debug
                # unexpected exceptions
                log.exception("任务异常 %s(%s)", task.id, kind)
                msg = str(e) or "未知异常(详见日志文件)"
                self.log(task.id, f"任务异常: {msg}", "error")
                self._finish(task.id, status=STATUS_FAILED, error=msg)

        task._thread = threading.Thread(target=_run, daemon=True, name=f"task-{task.id}")
        task._thread.start()
        self._publish({"type": "task_started", "task_id": task.id,
                       "payload": task.to_snapshot(with_logs=False)})
        return task

    def get(self, task_id: str) -> Task | None:
        with self._lock:
            return self._tasks.get(task_id)

    def list(self) -> list[Task]:
        with self._lock:
            return sorted(self._tasks.values(), key=lambda t: t.created_at, reverse=True)

    def active(self) -> list[Task]:
        return [t for t in self.list() if t.status in (STATUS_PENDING, STATUS_RUNNING)]

    def cancel(self, task_id: str) -> bool:
        task = self.get(task_id)
        if task is None:
            return False
        task._flag.cancelled = True
        return True

    # ------------------------------------------------------------ 状态更新
    # ------------------------------------------------------------ Status updates

    def log(self, task_id: str, msg: str, level: str = "info") -> None:
        task = self.get(task_id)
        if task is None:
            return
        entry = {"t": time.time(), "level": level, "msg": msg}
        with self._lock:
            task.logs.append(entry)
        self._publish({"type": "task_log", "task_id": task_id, "payload": entry})

    def progress(self, task_id: str, value: float, label: str = "") -> None:
        task = self.get(task_id)
        if task is None:
            return
        task.progress = max(0.0, min(100.0, value))
        if label:
            task.progress_label = label
        self._publish({"type": "task_progress", "task_id": task_id,
                       "payload": {"progress": task.progress, "label": task.progress_label}})

    def _finish(self, task_id: str, status: str, result: dict | None = None,
                error: str | None = None) -> None:
        task = self.get(task_id)
        if task is None:
            return
        task.status = status
        task.finished_at = time.time()
        task.progress = 100.0 if status == STATUS_SUCCEEDED else task.progress
        if result is not None:
            task.result = result
        if error:
            task.error = error
        self._publish({"type": f"task_{status}", "task_id": task_id,
                       "payload": task.to_snapshot(with_logs=False)})
        self._prune_results()

    def _prune_results(self, keep: int = 3) -> None:
        """只保留最近 keep 个已完成任务的结果。比对结果 parsed 含每 HSP 全序列
        (qseq/sseq,数千 HSP 可达数十 MB),多次比对不清理会堆积数 GB 内存,
        拖慢后续所有操作(任务历史上限 MAX_TASK_HISTORY 只限条数,不限内存)。

        Keep only the results of the most recent keep finished tasks. The parsed
        BLAST results contain full sequences per HSP (qseq/sseq; thousands of
        HSPs can reach tens of MB), so repeated alignments without cleanup would
        accumulate gigabytes of memory and slow down everything else
        (MAX_TASK_HISTORY only caps the entry count, not memory).
        """
        with self._lock:
            finished = sorted(
                (t for t in self._tasks.values()
                 if t.finished_at is not None and t.result is not None),
                key=lambda t: t.finished_at or 0)
            for t in finished[:-keep]:
                t.result = None

    # ------------------------------------------------------------ SSE

    def subscribe(self, queue: asyncio.Queue, loop: asyncio.AbstractEventLoop) -> None:
        if self._loop is None:
            self._loop = loop
        with self._lock:
            self._queues.add(queue)

    def unsubscribe(self, queue: asyncio.Queue) -> None:
        with self._lock:
            self._queues.discard(queue)

    def _publish(self, event: dict) -> None:
        """把事件投递到所有 SSE 订阅队列(线程 → asyncio 桥)。

        Deliver events to all SSE subscription queues (thread → asyncio bridge).
        """
        with self._lock:
            queues = list(self._queues)
            loop = self._loop
        if not queues or loop is None:
            return

        async def _put(q: asyncio.Queue) -> None:
            await q.put(event)

        try:
            for q in queues:
                asyncio.run_coroutine_threadsafe(_put(q), loop)
        except (RuntimeError, ValueError):
            pass  # loop 已关闭(服务器退出),忽略 (loop already closed, server exiting; ignore)


manager = TaskManager()
