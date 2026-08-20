"""配置与状态持久化:config.json 读写、损坏备份、失效记录清理。

Configuration and state persistence: read/write of config.json, corruption
backup, and cleanup of stale records.
"""

from __future__ import annotations

import copy
import json
import logging
import os
import shutil
import sys
import threading
from pathlib import Path
from typing import Any

log = logging.getLogger("blastprime")


def default_data_dir() -> Path:
    """默认程序数据目录(存放 config.json 与 blastdb_* 数据库目录):
    - 源码模式:项目根目录下 local_blast_dbs/;
    - 打包模式:exe 同级 local_blast_dbs/(可写时),否则
      %APPDATA%/BlastPrimeStudio/local_blast_dbs/。
    注意:onefile 的 sys._MEIPASS 是每次启动重建的临时解包目录,把数据写进去
    会在进程退出时丢失 —— 绝不能作为数据目录。

    Default program data dir (holds config.json and blastdb_* database dirs):
    - source mode: local_blast_dbs/ under the project root;
    - packaged mode: local_blast_dbs/ next to the exe (when writable), else
      %APPDATA%/BlastPrimeStudio/local_blast_dbs/.
    sys._MEIPASS (onefile) is a temp extraction dir rebuilt on every start;
    data written there is lost on exit — never used as the data dir.
    """
    if getattr(sys, "frozen", False):
        exe_dir = Path(sys.executable).resolve().parent
        cand = exe_dir / "local_blast_dbs"
        try:
            cand.mkdir(parents=True, exist_ok=True)
            probe = cand / ".write_probe"
            probe.write_text("ok", encoding="utf-8")
            probe.unlink()
            return cand
        except OSError:
            pass
        base = (os.environ.get("APPDATA")
                or os.environ.get("LOCALAPPDATA")
                or str(Path.home()))
        return Path(base) / "BlastPrimeStudio" / "local_blast_dbs"
    return Path(__file__).resolve().parent.parent / "local_blast_dbs"

DEFAULT_PRIMER_PARAMS: dict[str, Any] = {
    "mode": "standard",            # standard | sgrna | single
    "tm_min": 55.0, "tm_opt": 60.0, "tm_max": 65.0,
    "gc_min": 30.0, "gc_max": 70.0,
    "primer_len_min": 18, "primer_len_max": 25,
    "product_len_mode": "absolute", # absolute | relative | unlimited
    "product_len_min": 150, "product_len_max": 300,
    "product_offset1": 0.0,         # 相对产物偏移下限(基准=目标长度/模板长度) (Lower bound of relative product offset; base = target length / template length)
    "product_offset2": 300.0,       # 相对产物偏移上限 (Upper bound of relative product offset)
    "flank_extension": 150,         # 侧翼延伸长度 bp (Flank extension length in bp)
    "max_dimer": 5,                 # 最大引物二聚体 bp (Maximum primer dimer length in bp)
    "max_tm_diff": 2.0,             # 上下游最大 Tm 差 (Maximum Tm difference between forward and reverse primers)
    "max_gc_clamp_3p": 3,           # 3' 端 5 bp 最大 GC 数 (Maximum GC count in the terminal 5 bp at the 3' end)
    "high_risk_threshold": 8,       # 高危 3' 端阈值 (High-risk 3' end threshold)
    "level2_global_th": 0.6934,     # 阶段二放行阈值(global/3',count 2-3 边界,3^(-1/3))(Stage-2 release threshold)
    "level3_global_th": 0.5503,     # 阶段三放行阈值(global/3',count 4-6 边界,6^(-1/3))(Stage-3 release threshold)
    "buffer_len": 8,                # 可用区两侧缓冲:连续可用区上下游各外扩 bp (Both-sides available-region buffer, bp)
    "stage4_pool_min": 20, "stage4_pool_max": 50,   # 旧算法保留 (legacy, unused by the new pipeline)
    "offtarget_product_min": 50, "offtarget_product_max": 4000,
    "target_buffer": 50,            # 目标区内命中判定缓冲 ±50 bp (Hit-determination buffer of ±50 bp within the target region)
    # ---- 新特异性引擎(guide_sup2.md) ----
    "specificity_kmers": [8, 10, 12, 15],      # 多尺度 k-mer(§9)(Multi-scale k-mers)
    "three_prime_windows": [8, 10, 12, 15],    # 3' 端窗口尺度(§18)(3'-end window scales)
    "prefilter_kmers": [10, 12, 15],           # seed 预筛尺度(§26):不含 8-mer,8-mer 在真实基因组必然高频(2N/4^8≈396@26Mb),预筛只拦"明显不可能特异" (Prefilter scales; no 8-mer: its expected count 2N/4^8≈396@26Mb exceeds max_hits, the prefilter only rejects "clearly impossible")
    "seed_k": 12,                  # binding-site 种子长度(§26)(Binding-site seed length)
    "prefilter_max_hits": 200,     # seed 预筛最大命中数:超出即淘汰/截断(§26)
    "binding_min_identity": 80,    # 结合位点最小身份比 %(§27,低于此 = 低 Tm 无效短命中,非真实结合位点;R24:50% 放行 11 bp 窗口命中导致大量假脱靶位点 → truncated 误淘汰)(Minimum binding-site identity %; below it = low-Tm invalid short hit, not a real binding site)
    "candidate_count": 50,         # 每级 primer3 候选池大小(§25)
    "score_physical_weight": 0.5,  # 综合分物理权重(§39)(Physical weight of the composite score)
    "score_specificity_weight": 0.5,  # 综合分特异性权重(§39)
    "blast_evalue": 10.0,
    "blast_max_targets": 5000,
    "timeout_sec": 600,
    "sgrna_len": 20, "sgrna_pam": "NGG",
    "sgrna_target_only": True,   # sgRNA 引导序列必须位于目标区段内(§7,默认开启) (sgRNA guide sequence must lie within the target region (§7, enabled by default))
    "skip_spec_eval": False,     # 是否跳过 blastn-short 逆向特异性验证(特异性分按 k-mer 深度剖面评分)
    "skip_kmer_scoring": False,  # 是否跳过 k-mer 评分过程(不建深度/剖面,四段式按全放行设计)
}

def detect_system_lang() -> str:
    """按系统语言返回 zh/en:优先环境变量(LC_ALL/LC_MESSAGES/LANG),
    Windows 无这些变量时读系统 UI 语言;都无法判定则默认英文。

    Return "zh"/"en" by system language: prefer environment variables
    (LC_ALL/LC_MESSAGES/LANG); on Windows without these, read the system UI
    language; default to English if undeterminable.
    """
    for var in ("LC_ALL", "LC_MESSAGES", "LANG"):
        v = (os.environ.get(var) or "").lower()
        if v.startswith("zh"):
            return "zh"
        if v.startswith("en"):
            return "en"
    try:  # Windows 控制台通常无 LANG 环境变量,读取系统 UI 语言(主语言 0x04 = 中文) (Windows consoles usually lack LANG; read the system UI language, primary language 0x04 = Chinese)
        import ctypes
        lang_id = ctypes.windll.kernel32.GetUserDefaultUILanguage() & 0x3FF
        return "zh" if lang_id == 0x04 else "en"
    except Exception:  # noqa: BLE001 — Linux 等无 windll 的环境 (environments without windll, e.g. Linux)
        return "en"


DEFAULT_CONFIG: dict[str, Any] = {
    "app": "BlastPrimeStudio",
    "version": 1,
    "lang": detect_system_lang(),  # zh | en:创建时按系统语言,无法判定默认英文 (zh | en: set by system language at creation, defaults to English if undeterminable)
    "theme": "system",             # system | light | dark
    "loglevel": "INFO",
    "logfile": "",
    "blast_bin_dir": "",           # 手动指定的 BLAST 目录 (Manually specified BLAST binary directory)
    "data_dir": "",                # 默认数据库存储路径(留空 = default_data_dir())(Default database storage dir; empty = default_data_dir())
    "db_records": [],              # [{"prefix": str, "is_created": bool}]
    "primer_params": copy.deepcopy(DEFAULT_PRIMER_PARAMS),
}


def get_data_dir() -> Path:
    """当前生效的数据目录:config 中 data_dir 设置优先,未设置 = 默认。
    设置页可修改;修改后新建数据库立即落到新目录。

    Effective data dir: the config's data_dir setting wins, else the default.
    Editable from the settings dialog; new databases land there immediately.
    """
    v = (get_config().get("data_dir") or "").strip()
    return Path(v).expanduser() if v else default_data_dir()


class Config:
    """线程安全的配置对象,带损坏备份与自动清理。

    Thread-safe configuration object with corruption backup and automatic cleanup.
    """

    def __init__(self, path: str | Path | None = None):
        self.path = Path(path) if path else default_data_dir() / "config.json"
        self._lock = threading.Lock()
        self.data: dict[str, Any] = copy.deepcopy(DEFAULT_CONFIG)
        self.loaded_ok = True
        self.backup_made = False
        self.last_save_error: str | None = None
        self._load()

    def _load(self) -> None:
        if not self.path.exists():
            # 首次启动:创建失败(如目录不可写)不中断启动,内存配置继续可用,
            # 失败细节记录到日志;后续 set() 写入失败会抛 OSError 通知用户
            # First start: a failed create (e.g. unwritable dir) must not abort
            # startup — the in-memory config keeps working and the failure is
            # logged; later set() writes raise OSError to surface the problem
            try:
                self._save()
            except OSError:
                pass
            return
        try:
            raw = json.loads(self.path.read_text(encoding="utf-8"))
            if not isinstance(raw, dict):
                raise ValueError("config 根节点不是对象")
            self.data = self._merge(copy.deepcopy(DEFAULT_CONFIG), raw)
            # 剔除废弃键(旧版设置页"配置文件路径"的引导链记录;配置迁移改为
            # 设置页导入/导出 + --config 参数,不再有 config_path 引导)
            # Drop the deprecated key (the old settings-page "config file path"
            # bootstrap record; config relocation is now import/export + --config)
            self.data.pop("config_path", None)
            # 剔除 primer_params 遗留键(全库 k-mer 索引 pickle 缓存已随
            # 设计流程改版移除,index_version/index_cache_enabled 不再使用)
            # Drop legacy primer_params keys (the full-db k-mer index pickle
            # cache was removed with the pipeline rework; index_version /
            # index_cache_enabled are no longer used)
            pp = self.data.get("primer_params")
            if isinstance(pp, dict):
                pp.pop("index_version", None)
                pp.pop("index_cache_enabled", None)
        except Exception:
            # 配置损坏:备份原文件,恢复默认
            # Config corrupted: back up the original file and restore defaults.
            self.backup_made = True
            try:
                bak = self.path.with_suffix(self.path.suffix + ".bak")
                shutil.copy2(self.path, bak)
            except OSError:
                pass
            log.warning("配置文件损坏,已备份为 .bak 并恢复默认: %s", self.path)
            self.data = copy.deepcopy(DEFAULT_CONFIG)
            try:
                self._save()
            except OSError:
                pass
        self._cleanup_stale_records()

    @staticmethod
    def _merge(base: dict, override: dict) -> dict:
        out = copy.deepcopy(base)
        for k, v in override.items():
            if k in out and isinstance(out[k], dict) and isinstance(v, dict):
                out[k] = Config._merge(out[k], v)
            else:
                out[k] = v
        return out

    def _cleanup_stale_records(self) -> None:
        """载入配置时自动剔除索引文件已不存在的历史条目,并即时重写。

        On config load, drop history entries whose index files no longer exist
        and rewrite immediately.
        """
        keep = [r for r in self.data.get("db_records", []) if self._record_valid(r)]
        if len(keep) != len(self.data.get("db_records", [])):
            self.data["db_records"] = keep
            self._save()

    @staticmethod
    def _record_valid(rec: dict) -> bool:
        prefix = rec.get("prefix", "")
        if not prefix:
            return False
        p = Path(prefix)
        for ext in ("nin", "pin", "nsq", "psq"):
            if (p.with_suffix("." + ext)).exists():
                return True
        return False

    def get(self, key: str, default: Any = None) -> Any:
        with self._lock:
            return self.data.get(key, default)

    def set(self, key: str, value: Any) -> None:
        with self._lock:
            self.data[key] = value
            self._save()

    def primer_params(self) -> dict[str, Any]:
        with self._lock:
            return copy.deepcopy(self.data.get("primer_params", DEFAULT_PRIMER_PARAMS))

    def set_primer_params(self, params: dict[str, Any]) -> None:
        with self._lock:
            merged = self._merge(copy.deepcopy(DEFAULT_PRIMER_PARAMS), params)
            self.data["primer_params"] = merged
            self._save()

    def snapshot(self) -> dict[str, Any]:
        """当前配置的完整深拷贝(下载配置用,剔除废弃键)。

        Deep copy of the whole config (for the download-config feature;
        deprecated keys stripped).
        """
        with self._lock:
            data = copy.deepcopy(self.data)
        data.pop("config_path", None)
        return data

    def replace_data(self, raw: dict[str, Any]) -> None:
        """整体替换配置数据(设置页"读取配置"导入):与默认值合并,剔除废弃键,
        校验类型,清理失效记录,保存到当前配置文件位置。导入内容里的任何
        路径字段都不改变程序对当前配置文件的定位。

        Replace the whole config data (settings-dialog "import config"):
        merge with defaults, strip deprecated keys, validate types, clean up
        stale records, and persist to the current config file location. Path
        fields inside the imported content never change where the program
        looks for its config file.
        """
        merged = self._merge(copy.deepcopy(DEFAULT_CONFIG), raw)
        merged.pop("config_path", None)
        if not isinstance(merged.get("primer_params"), dict):
            merged["primer_params"] = copy.deepcopy(DEFAULT_PRIMER_PARAMS)
        if not isinstance(merged.get("db_records"), list):
            merged["db_records"] = []
        with self._lock:
            self.data = merged
            self._save()
        self._cleanup_stale_records()

    def db_records(self) -> list[dict]:
        with self._lock:
            return copy.deepcopy(self.data.get("db_records", []))

    def add_db_record(self, prefix: str, is_created: bool, note: str = "") -> None:
        with self._lock:
            recs = [r for r in self.data.get("db_records", []) if r.get("prefix") != prefix]
            recs.insert(0, {"prefix": prefix, "is_created": bool(is_created), "note": note or ""})
            self.data["db_records"] = recs
            self._save()

    def set_db_note(self, prefix: str, note: str) -> None:
        with self._lock:
            for r in self.data.get("db_records", []):
                if r.get("prefix") == prefix:
                    r["note"] = note or ""
                    self._save()
                    return
            raise KeyError(f"记录不存在: {prefix}")

    def remove_db_record(self, prefix: str) -> None:
        with self._lock:
            self.data["db_records"] = [
                r for r in self.data.get("db_records", [])
                if r.get("prefix") != prefix
            ]
            self._save()

    def reorder_db_records(self, prefixes: list[str]) -> None:
        """按用户手动排序重排记录(持久化);未列出的记录保持原相对顺序收尾。

        Reorder records to follow the user's manual order (persisted);
        records not listed keep their original relative order at the end.
        """
        with self._lock:
            recs = self.data.get("db_records", [])
            by_prefix = {r.get("prefix"): r for r in recs}
            ordered = []
            for p in prefixes:
                r = by_prefix.pop(p, None)
                if r is not None:
                    ordered.append(r)
            ordered.extend(by_prefix.values())
            self.data["db_records"] = ordered
            self._save()

    def is_created_record(self, prefix: str) -> bool:
        for r in self.data.get("db_records", []):
            if r.get("prefix") == prefix:
                return bool(r.get("is_created", False))
        return False

    def _save(self) -> None:
        """写回配置文件;失败记录到日志并抛 OSError(调用方可见,不再静默)。

        此前 OSError 被吞掉,exe 打包模式下数据目录解析到一次性临时目录时,
        保存设置界面仍提示"已保存",实际全部丢失 —— 静默失败正是排查难的根源。

        Persist the config file; on failure log the error and raise OSError
        instead of silently swallowing it (the old `except OSError: pass`
        made the UI claim "saved" while data vanished into a temp dir in
        packaged mode — silent failure is exactly what made debugging hard).
        """
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.path.write_text(
                json.dumps(self.data, ensure_ascii=False, indent=2), encoding="utf-8"
            )
            self.last_save_error = None
        except OSError as e:
            self.last_save_error = f"{self.path}: {e}"
            log.error("写入配置文件失败: %s", self.last_save_error)
            raise


# 全局单例
# Global singleton
config: Config | None = None


def get_config() -> Config:
    global config
    if config is None:
        config = Config()
    return config


def set_config_path(path: str) -> Config:
    """用指定路径重新加载配置并切换全局单例(`--config` 运行参数)。
    显式指定优先:只影响本次运行的读写位置,不改写文件内容、不持久化;
    下次不带参数启动回到数据目录的默认配置。配置的迁移/备份走设置页的
    读取配置/下载配置(导入导出),不设引导链。

    Load the config from an explicit path and swap the global singleton
    (the `--config` CLI flag). An explicit path wins: it only changes where
    this run reads/writes, never rewrites the file or persists anything; the
    next run without the flag returns to the default config in the data dir.
    Config migration/backup goes through the settings dialog's import/export
    (read/download config); there is no bootstrap chain.
    """
    global config
    cfg = Config(path)
    config = cfg
    return cfg
