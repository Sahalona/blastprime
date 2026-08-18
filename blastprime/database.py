"""数据库历史记录、浏览已有库、删除等管理能力。

Database history records, browsing existing databases, deletion, and other
management capabilities.
"""

from __future__ import annotations

import random
import re
import shutil
import string
from pathlib import Path

from .config import DATA_DIR, get_config

INDEX_EXTS = [
    "nin", "nsq", "nhr", "pin", "psq", "phr",
    "nsd", "nsi", "nnd", "nni", "njs", "nog", "pog", "psd", "psi", "pnd", "pni", "pjs",
]

DB_NAME_RE = re.compile(r"^[A-Za-z0-9_]+$")


def default_db_dir() -> Path:
    """未指定输出目录时的默认存储:local_blast_dbs/blastdb_<随机后缀>/

    Default storage when no output directory is specified:
    local_blast_dbs/blastdb_<random suffix>/
    """
    suffix = "".join(random.choices(string.ascii_lowercase + string.digits, k=8))
    return DATA_DIR / f"blastdb_{suffix}"


def validate_db_name(name: str) -> str:
    name = name.strip()
    if not name:
        raise ValueError("数据库名称不能为空")
    if not DB_NAME_RE.match(name):
        raise ValueError("数据库名称仅允许字母、数字、下划线")
    return name


def record_valid(prefix: str) -> bool:
    from .config import Config
    return Config._record_valid({"prefix": prefix})


def list_records() -> list[dict]:
    return get_config().db_records()


def add_record(prefix: str, is_created: bool, note: str = "") -> None:
    get_config().add_db_record(prefix, is_created, note)


def set_note(prefix: str, note: str) -> None:
    get_config().set_db_note(prefix, note)


def reorder_records(prefixes: list[str]) -> None:
    """按用户手动排序重排历史记录,持久化到 config.json。

    Reorder history records by the user's manual sort, persisted to config.json.
    """
    get_config().reorder_db_records([p for p in prefixes if p])


def remove_record(prefix: str) -> None:
    """仅移除历史记录,保留磁盘文件。

    Remove only the history record, keeping the files on disk.
    """
    get_config().remove_db_record(prefix)


def physically_delete(prefix: str) -> None:
    """物理删除:删除全部索引文件;默认目录下的库连空父目录一并删除。

    Physical deletion: remove all index files; databases under the default
    directory also get their empty parent directory removed.
    """
    p = Path(prefix)
    deleted = 0
    if p.is_file():
        # 用户给的可能是文件名(如 xxx.nin);若指向索引文件则删除该库的成套索引
        # The user may provide a file name (e.g. xxx.nin); if it points to an
        # index file, delete that database's full set of indexes.
        ext = p.suffix.lstrip(".")
        if ext in INDEX_EXTS:
            base = p.with_suffix("")
            for e in INDEX_EXTS:
                f = Path(str(base) + "." + e)
                if f.exists():
                    f.unlink()
                    deleted += 1
            parent = p.parent
            if parent.name.startswith("blastdb_") and not any(parent.iterdir()):
                shutil.rmtree(parent, ignore_errors=True)
        else:
            p.unlink()
            deleted += 1
    elif p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
        deleted = 1
    else:
        parent = p.parent
        found = False
        # 1) prefix 的父目录正是默认数据目录下的 blastdb_* 目录(避免同名库删错)
        # 1) prefix's parent is a blastdb_* directory under the default data
        #    directory (avoid deleting a same-named database by mistake)
        if parent.name.startswith("blastdb_") and parent.parent == DATA_DIR and parent.exists():
            shutil.rmtree(parent, ignore_errors=True)
            found = True
            deleted = 1
        # 2) 外部目录中的库(如手动导入):删除该库成套索引,不删父目录
        # 2) Databases in external directories (e.g. manually imported): delete
        #    the full index set but keep the parent directory.
        elif parent.is_dir() and any(Path(str(p) + "." + e).exists() for e in INDEX_EXTS):
            for e in INDEX_EXTS:
                f = Path(str(p) + "." + e)
                if f.exists():
                    f.unlink()
                    deleted += 1
            found = True
        # 3) 兼容:前缀为 basename 形式(库名),尝试在默认数据目录下找 blastdb_* 匹配
        # 3) Compatibility: prefix is a basename (database name); try to find a
        #    matching blastdb_* under the default data directory.
        if not found and DATA_DIR.is_dir():
            for d in DATA_DIR.iterdir():
                if d.name.startswith("blastdb_") and (d / (p.name + ".nin")).exists():
                    shutil.rmtree(d, ignore_errors=True)
                    found = True
                    deleted = 1
                    break
        if not found:
            raise ValueError("未找到该数据库的索引文件")
    return {"deleted_files": deleted}


def scan_default_dir() -> int:
    """扫描默认数据目录,自动识别手动复制进来的 BLAST 索引库(幂等)。

    遍历 DATA_DIR/blastdb_*/ 下所有 *.nin/*.pin,配套文件齐全且未在历史记录中
    的库自动加入记录(is_created=False),返回本次新发现的数量。

    Scan the default data directory and automatically recognize BLAST index
    databases copied in manually (idempotent).

    Iterate all *.nin/*.pin under DATA_DIR/blastdb_*/; databases with complete
    companion files that are not yet in the history are added automatically
    (is_created=False). Returns the number newly discovered.
    """
    if not DATA_DIR.is_dir():
        return 0
    cfg = get_config()
    known = {r["prefix"] for r in cfg.db_records()}
    required = {"nin": "nsq", "pin": "psq"}
    found = 0
    for d in sorted(DATA_DIR.iterdir()):
        if not d.is_dir() or not d.name.startswith("blastdb_"):
            continue
        for idx in sorted(d.iterdir()):
            ext = idx.suffix.lstrip(".")
            if ext not in required:
                continue
            base = idx.with_suffix("")
            if not Path(str(base) + "." + required[ext]).exists():
                continue
            prefix = str(base)
            if prefix in known:
                continue
            add_record(prefix, is_created=False)
            known.add(prefix)
            found += 1
    return found


def browse_existing(index_file: str) -> dict:
    """手动选择磁盘上的索引文件(*.nin/*.pin)加入历史。

    Manually select an index file (*.nin/*.pin) on disk to add to the history.
    """
    p = Path(index_file)
    ext = p.suffix.lstrip(".")
    if ext not in ("nin", "pin"):
        raise ValueError("请选择 BLAST 索引文件(*.nin / *.pin)")
    if not p.is_file():
        raise ValueError(f"索引文件不存在: {index_file}")
    base = p.with_suffix("")
    # 校验配套文件存在
    # Verify the companion files exist.
    required = {"nin": "nsq", "pin": "psq"}
    if not Path(str(base) + "." + required[ext]).exists():
        raise ValueError(f"索引文件不完整:缺少 {base}.{required[ext]}")
    prefix = str(base)
    add_record(prefix, is_created=False)
    return {"prefix": prefix, "is_created": False}


def db_type_of(prefix: str) -> str | None:
    """根据存在的索引文件判断库类型:nucl / prot / None。

    Determine the database type from the index files present: nucl / prot / None.
    """
    p = Path(prefix)
    if p.is_file() and p.suffix.lstrip(".") in INDEX_EXTS:
        p = p.with_suffix("")
    if (Path(str(p) + ".nin")).exists():
        return "nucl"
    if (Path(str(p) + ".pin")).exists():
        return "prot"
    return None
