"""多尺度 k-mer 特异性索引与逐碱基特异性剖面(guide_sup2.md §8-§19)。

设计要点
--------
- 计数索引:8/10/12-mer 用 2-bit 编码定长数组 ``array('I')``(4^k 定界),
  15-mer 用 dict 只存**重复 ≥2 次**的 canonical k-mer(数学依据:15-mer
  重复 ⟹ 其任一子串同位置必重复,因此只需在较短尺度重复时统计 15-mer)。
  唯一 k-mer 缺席即得 1.0 分(缺席 == 全库仅此一处),故计数索引无需存储
  任何唯一 k-mer —— 随机基因组下索引近乎为空,内存与构建近零开销。
- 位置索引:seed_k(默认 12)mer 的 canonical → [(seq_idx, pos, strand)],
  同样只存重复 ≥2 次的 k-mer;供 3' 端 seed 预筛后的全长局部验证使用。
- canonical 化:``canonical(kmer) = min(code, rc_code(code))``,occurrence
  保留实际方向(查询端按 code == canonical 判定正链,否则负链)。
- N 处理:任何含非 ACGT 碱基的 k-mer 一律跳过(不计数、不进位置索引)。
- 构建性能:序列预编码为 2-bit 字节表,滚动窗口一次产出 k_eff-mer 及其
  全部前缀 k-mer(同一开始位置的各尺度 k-mer 互为前缀),反向互补用
  预计算查表(5/8/10-mer)或 O(k) 函数(12/15-mer)。
- pickle 缓存:DATA_DIR/.primer_index/<sha1>.pkl,签名 = index_version +
  kmer_set + seed_k + 各库索引文件 (mtime_ns, size);库变化自动失效。

复杂度(§67):index build O(KD)(两趟扫描);target 扫描 O(KN);
per-base profile O(KN)(差分数组);局部验证 O(HL)。内存 O(KD),唯一区近零。

Multi-scale k-mer specificity index and per-base specificity profile
(guide_sup2.md §8-§19).
"""

from __future__ import annotations

import hashlib
import io
import math
import os
import pickle
import tempfile
from array import array
from pathlib import Path
from typing import Callable

from Bio import SeqIO

from . import blast
from .blast import CancelFlag
from .config import DATA_DIR

INDEX_VERSION = 1

# 2-bit 编码:00=A 01=C 10=G 11=T;255 = 非 ACGT(N 等,不参与计数)
# 2-bit encoding: 00=A 01=C 10=G 11=T; 255 = non-ACGT (N etc., not counted)
_BASE_TABLE = bytearray([255] * 256)
for _ch, _code in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
    _BASE_TABLE[ord(_ch)] = _code
_BASE_TABLE = bytes(_BASE_TABLE)

_RC_TABLES: dict[int, array] = {}


def _rc_table(k: int) -> array:
    t = _RC_TABLES.get(k)
    if t is None:
        t = array("I", [0]) * (1 << (2 * k))
        for code in range(1 << (2 * k)):
            t[code] = rc_code(code, k)
        _RC_TABLES[k] = t
    return t


def _base_codes(s: str) -> bytes:
    return s.encode("ascii").translate(_BASE_TABLE)


def rc_code(code: int, k: int) -> int:
    """2-bit 整数的反向互补(补码 3-b + 逆序)。

    Reverse complement of a 2-bit integer (complement 3-b + reversal).
    """
    out = 0
    for _ in range(k):
        out = (out << 2) | (3 - (code & 3))
        code >>= 2
    return out


def canonical_code(code: int, k: int) -> int:
    """canonical(kmer) = min(kmer, revcomp(kmer))(§9);短 k 用查表。"""
    if k in (5, 8, 10):
        return min(code, _rc_table(k)[code])
    return min(code, rc_code(code, k))


def encode_string(kmer: str) -> tuple[int, bool]:
    """kmer 字符串 → 2-bit 整数;含非 ACGT 返回 ok=False。"""
    code = 0
    for ch in kmer:
        b = _BASE_TO_CODE.get(ch)
        if b is None:
            return 0, False
        code = (code << 2) | b
    return code, True


_BASE_TO_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def kmer_occurrence_score(count: int) -> float:
    """k-mer occurrence score(§14,独立策略函数,方案 B-k3)。

    幂函数 score = count^(-1/3)(count≥2),count≤1(缺席/唯一)→ 1.0。
    旧式 1/(1+log2(count)) 在整数 count 下分数只有 {1.0, 0.5, 0.387, ...}
    两档(0.5~1.0 之间无整数 count 对应),global 剖面两极分化、四段式
    L1/L2/L3 阈值在 global 维度全部等价(都要求 count=1);幂函数使
    count=2 → 0.794、count=3 → 0.693、count=8 → 0.5,产生平滑过渡区,
    且 0.80 阈值恰好对应 count=2 —— L3 首次获得真实区分度。
    唯一 → 1.0,单调递减,可解释。

    k-mer occurrence score (§14, independent policy function, option B-k3).

    Power function score = count^(-1/3) for count ≥ 2; count ≤ 1 (absent /
    unique) → 1.0. The old 1/(1+log2(count)) yields only {1.0, 0.5, 0.387,
    ...} for integer counts (no count maps into 0.5~1.0), polarising the
    global profile and making the L1/L2/L3 thresholds equivalent in the
    global dimension (all require count=1). The power function gives
    count=2 → 0.794, count=3 → 0.693, count=8 → 0.5 — a smooth transition
    band, and 0.80 maps exactly to count=2, giving L3 real discrimination.
    Unique → 1.0, monotone decreasing, explainable.
    """
    if count <= 1:
        return 1.0
    return count ** (-1.0 / 3.0)


class KmerIndex:
    """数据库多尺度 k-mer 索引。

    ``db_prefix`` 路径从真实 BLAST 库构建(带 pickle 缓存);
    ``sequences`` 路径直接从序列列表构建(单元测试/无 blastdbcmd 环境用)。
    """

    def __init__(
        self,
        db_prefix: str | None = None,
        sequences: list[str] | None = None,
        seq_ids: list[str] | None = None,
        kmer_set: tuple[int, ...] = (8, 10, 12, 15),
        seed_k: int = 12,
        three_prime_windows: tuple[int, ...] = (8, 10, 12, 15),
        use_cache: bool = True,
        on_log: Callable[[str], None] | None = None,
        on_progress: Callable[[float], None] | None = None,
        cancel: CancelFlag | None = None,
    ):
        self.db_prefix = db_prefix
        self.kmer_set = tuple(sorted(set(kmer_set)))
        self.seed_k = seed_k
        self.three_prime_windows = tuple(sorted(set(three_prime_windows)))
        self.sequences: list[str] = []
        self.seq_ids: list[str] = []
        # 计数尺度 = 配置尺度 ∪ 3' 端窗口尺度(5/8/10 供 §18 剖面使用)
        # Count scales = configured scales ∪ 3'-end window scales (5/8/10 serve §18)
        self.count_ks = tuple(sorted(set(self.kmer_set) | set(self.three_prime_windows)))
        self.counts: dict[int, object] = {}   # k -> array('I') 或 dict(canonical->count)
        self.positions: dict[int, list[tuple[int, int, int]]] = {}
        self.cached = False

        if sequences is not None:
            self._build_from_sequences(sequences, seq_ids or [],
                                       on_log=on_log, on_progress=on_progress,
                                       cancel=cancel)
            return
        if not db_prefix:
            raise ValueError("KmerIndex 需要 db_prefix 或 sequences 之一")

        sig = self._signature(db_prefix)
        cache_path = self._cache_path(db_prefix, sig)
        if use_cache and cache_path.exists():
            try:
                obj = pickle.loads(cache_path.read_bytes())
                if (obj.get("index_version") == INDEX_VERSION
                        and obj.get("signature") == sig
                        and tuple(obj.get("kmer_set", ())) == self.kmer_set
                        and obj.get("seed_k") == seed_k):
                    self.sequences = obj["sequences"]
                    self.seq_ids = obj["seq_ids"]
                    self.counts = obj["counts"]
                    self.positions = obj["positions"]
                    self.cached = True
                    if on_log:
                        on_log(f"命中 k-mer 索引缓存: {cache_path.name}")
                    return
            except Exception:
                pass  # 缓存损坏 → 重建 (corrupt cache → rebuild)

        self._build_from_db(db_prefix, on_log=on_log, on_progress=on_progress,
                            cancel=cancel)
        if use_cache:
            try:
                cache_path.parent.mkdir(parents=True, exist_ok=True)
                payload = {
                    "index_version": INDEX_VERSION,
                    "signature": sig, "kmer_set": list(self.kmer_set),
                    "seed_k": seed_k, "seq_ids": self.seq_ids,
                    "sequences": self.sequences, "counts": self.counts,
                    "positions": self.positions,
                }
                fd, tmp = tempfile.mkstemp(dir=str(cache_path.parent),
                                           prefix=".idx_", suffix=".pkl")
                with os.fdopen(fd, "wb") as f:
                    pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
                os.replace(tmp, cache_path)   # 原子替换,并发双建无害 (atomic; racing double-build harmless)
            except OSError:
                pass

    # ------------------------------------------------------------ 签名与缓存

    @staticmethod
    def _signature(db_prefix: str) -> str:
        """库文件 (mtime_ns, size) 签名;任一索引文件变化即失效。

        Database-file (mtime_ns, size) signature; any index-file change
        invalidates the cache.
        """
        parts = [f"v{INDEX_VERSION}"]
        p = Path(db_prefix)
        for ext in ("nin", "pin", "nsq", "psq", "nhr", "phr"):
            f = p.with_suffix("." + ext)
            try:
                st = f.stat()
                parts.append(f"{ext}:{st.st_mtime_ns}:{st.st_size}")
            except OSError:
                pass
        return hashlib.sha1("|".join(parts).encode()).hexdigest()

    @staticmethod
    def _cache_path(db_prefix: str, sig: str) -> Path:
        name = Path(db_prefix).name
        return DATA_DIR / ".primer_index" / f"{name}_{sig}.pkl"

    # ------------------------------------------------------------ 构建 (Build)

    def _build_from_sequences(self, sequences, seq_ids, on_log, on_progress, cancel):
        self.sequences = [s.upper() for s in sequences]
        self.seq_ids = list(seq_ids) or [f"seq{i + 1}" for i in range(len(sequences))]
        if on_log:
            on_log(f"构建 k-mer 索引: {len(self.sequences)} 条序列, "
                   f"{sum(map(len, self.sequences))} bp, k={self.count_ks}")
        self._scan(0.0, 1.0, on_progress, cancel)

    def _build_from_db(self, db_prefix, on_log, on_progress, cancel):
        if on_log:
            on_log(f"构建 k-mer 索引: {db_prefix} ...")
        # 取全库序列(一次 blastdbcmd,§37 序列访问)
        # Fetch all sequences in one blastdbcmd call (§37 sequence access)
        fa = blast.fetch_entry(db_prefix, "all")
        seqs: list[str] = []
        ids: list[str] = []
        for rec in SeqIO.parse(io.StringIO(fa), "fasta"):
            if cancel is not None and cancel.cancelled:
                raise blast.BlastError("任务已取消")
            seqs.append(str(rec.seq).upper())
            ids.append(str(rec.id))
        if on_log:
            on_log(f"库序列 {len(seqs)} 条 / {sum(map(len, seqs))} bp")
        self.sequences = seqs
        self.seq_ids = ids
        self._scan(0.0, 1.0, on_progress, cancel)

    def _scan(self, t0: float, t1: float, on_progress, cancel) -> None:
        """两趟扫描:

        趟 1:各尺度定长数组计数(滚动窗口 + 前缀提取,§12 2-bit);
        趟 2:重复 seed_k-mer 的位置索引 + 重复 15-mer 计数
        (15-mer 重复 ⟹ 其子串同位置重复,故只在较短尺度重复时统计)。
        """
        count_ks = self.count_ks
        arrays: dict[int, array] = {}
        for k in count_ks:
            if k <= 12:
                arrays[k] = array("I", [0]) * (1 << (2 * k))
        use15 = 15 in count_ks
        kmax = 15 if use15 else max(count_ks)
        total = sum(len(s) for s in self.sequences) or 1
        done = 0
        for si, s in enumerate(self.sequences):
            if cancel is not None and cancel.cancelled:
                raise blast.BlastError("任务已取消")
            n = len(s)
            if n >= min(count_ks):
                self._scan_counts(s, arrays, count_ks, kmax)
            done += n
            if on_progress and (si & 0x3F) == 0:
                on_progress(t0 + 0.7 * done / total)
        # 趟 2:重复 seed_k 位置 + 重复 15-mer 计数
        seed_arr = arrays.get(self.seed_k)
        positions: dict[int, list] = {}
        counts15: dict[int, int] = {}
        done = 0
        for si, s in enumerate(self.sequences):
            if cancel is not None and cancel.cancelled:
                raise blast.BlastError("任务已取消")
            n = len(s)
            if seed_arr is not None and n >= self.seed_k:
                self._scan_positions(s, si, seed_arr, positions)
            if use15 and n >= 15:
                self._scan_15mers(s, arrays, counts15)
            done += n
            if on_progress and (si & 0x3F) == 0:
                on_progress(t0 + 0.7 + 0.3 * done / total)
        self.counts = arrays
        if counts15:
            # 只保留真正重复的 15-mer(§14:缺席=唯一)
            self.counts[15] = {c: v for c, v in counts15.items() if v >= 2}
        self.positions = positions

    @staticmethod
    def _scan_counts(s: str, arrays: dict, count_ks: tuple, kmax: int) -> None:
        """趟 1 单序列:滚动 k_eff-mer,各尺度 k-mer 取其高位前缀(同一起点)。

        N 门控按前缀判定:窗口含 N 时,其 N 之前的干净前缀仍计数(整窗跳过
        会漏掉 N 上游的合法 k-mer);快路径(整窗干净)一次加全部前缀。
        """
        n = len(s)
        k_eff = min(kmax, n)
        bs = _base_codes(s)
        mask = (1 << (2 * k_eff)) - 1
        code = 0
        bad = 0
        for j in range(k_eff):
            b = bs[j]
            code = (code << 2) | (b & 3)
            if b == 255:
                bad += 1
        if bad == 0:
            _add_prefix_counts(code, k_eff, arrays, count_ks)
        else:
            _add_clean_prefixes(code, k_eff, bs, 0, arrays, count_ks)
        for i in range(1, n - k_eff + 1):
            b_in = bs[i + k_eff - 1]
            b_out = bs[i - 1]
            code = ((code << 2) | (b_in & 3)) & mask
            if b_out == 255:
                bad -= 1
            if b_in == 255:
                bad += 1
            if bad == 0:
                _add_prefix_counts(code, k_eff, arrays, count_ks)
            else:
                _add_clean_prefixes(code, k_eff, bs, i, arrays, count_ks)
        # 尾部短 k-mer:起点在最后一个 k_eff 窗口之后(n-k_eff+1 .. n-k)的
        # k-mer 是窗口内部片段,主循环不计数 —— 漏计会让序列末尾的重复
        # k-mer 被误判为唯一(score 1.0)。规格(§10)要求 every valid position。
        # Tail short k-mers: starts beyond the last k_eff-window (n-k_eff+1
        # .. n-k) are interior fragments, never counted by the main loop —
        # repeats near the sequence end would be mis-scored as unique. The
        # spec (§10) requires every valid position.
        for k in count_ks:
            if k not in arrays or k >= k_eff:
                continue
            for i in range(n - k_eff + 1, n - k + 1):
                if 255 in bs[i:i + k]:
                    continue
                code = 0
                for b in bs[i:i + k]:
                    code = (code << 2) | (b & 3)
                c = canonical_code(code, k)
                arrays[k][c] += 1

    def _scan_positions(self, s: str, si: int, seed_arr: array,
                        positions: dict) -> None:
        """趟 2 单序列:重复 seed_k-mer 的位置(含方向)。"""
        n = len(s)
        k = self.seed_k
        bs = _base_codes(s)
        mask = (1 << (2 * k)) - 1
        code = 0
        bad = 0
        for j in range(k):
            b = bs[j]
            code = (code << 2) | (b & 3)
            if b == 255:
                bad += 1
        if bad == 0:
            _maybe_record(code, k, si, 0, seed_arr, positions)
        for i in range(1, n - k + 1):
            b_in = bs[i + k - 1]
            b_out = bs[i - 1]
            code = ((code << 2) | (b_in & 3)) & mask
            if b_out == 255:
                bad -= 1
            if b_in == 255:
                bad += 1
            if bad == 0:
                _maybe_record(code, k, si, i, seed_arr, positions)

    def _scan_15mers(self, s: str, arrays: dict, counts15: dict) -> None:
        """趟 2 单序列:15-mer 计数(仅当较短尺度子串重复时,必要条件门控)。"""
        n = len(s)
        bs = _base_codes(s)
        # 用最长的 <15 计数尺度做门控(必要条件:任一子串重复 ⟹ 15-mer 可能重复)
        gates = [k for k in sorted(self.count_ks) if 15 > k]
        gate_k = gates[-1] if gates else None
        mask = (1 << 30) - 1
        code = 0
        bad = 0
        for j in range(15):
            b = bs[j]
            code = (code << 2) | (b & 3)
            if b == 255:
                bad += 1
        if bad == 0:
            _maybe_count15(code, gate_k, arrays, counts15)
        for i in range(1, n - 14):
            b_in = bs[i + 14]
            b_out = bs[i - 1]
            code = ((code << 2) | (b_in & 3)) & mask
            if b_out == 255:
                bad -= 1
            if b_in == 255:
                bad += 1
            if bad == 0:
                _maybe_count15(code, gate_k, arrays, counts15)

    # ------------------------------------------------------------ 查询 (Query)

    def count_of(self, kmer: str) -> int:
        """kmer(两条链合计)在库中的 occurrence 数;缺席(唯一)返回 0。

        Occurrence count of a k-mer (both strands combined); 0 = absent
        (unique, score 1.0).
        """
        code, ok = encode_string(kmer)
        if not ok:
            return 0
        k = len(kmer)
        c = canonical_code(code, k)
        cnt = self.counts.get(k)
        if cnt is None:
            return 0
        if isinstance(cnt, dict):
            return cnt.get(c, 0)
        return cnt[c]

    def positions_of(self, kmer: str) -> list[tuple[int, int, int]]:
        """kmer 的全部重复位置 (seq_idx, pos0, strand),方向以查询为准;
        仅覆盖重复 ≥2 次的 k-mer(缺席返回空)。

        canonical 化把正/反向互补出现合并为一条记录:查询即 canonical 时
        直接返回全部记录;查询为其反向互补时,两条链的出现互换——
        canonical 正链出现 = 查询负链出现,canonical 负链出现 = 查询正链
        出现(方向标签翻转,位置不变)。两种形式的出现全部返回。

        All repeated positions of a k-mer (seq_idx, pos0, strand in the
        query's orientation); only k-mers occurring ≥ 2 times are indexed
        (absent → empty). Canonicalization merges forward/reverse-complement
        occurrences into one entry: when the query equals the canonical form
        all entries are returned as-is; when the query is its reverse
        complement the strand labels flip (canonical-forward == query-reverse
        and vice versa), and both occurrences are still returned.
        """
        if len(kmer) != self.seed_k:
            raise ValueError(f"positions 索引仅支持 seed_k={self.seed_k} 长度的 k-mer")
        code, ok = encode_string(kmer)
        if not ok:
            return []
        c = canonical_code(code, self.seed_k)
        q_is_canonical = code == c
        occ = self.positions.get(c)
        if not occ:
            return []
        if q_is_canonical:
            return list(occ)
        return [(o[0], o[1], 1 - o[2]) for o in occ]


def _add_prefix_counts(code: int, k_eff: int, arrays: dict, count_ks: tuple) -> None:
    for k in count_ks:
        if k > k_eff or k not in arrays:
            continue  # 15-mer 等 dict 尺度由趟 2 单独计数
        c = canonical_code(code >> (2 * (k_eff - k)), k)
        arrays[k][c] += 1


def _add_clean_prefixes(code: int, k_eff: int, bs: bytes, i: int,
                        arrays: dict, count_ks: tuple) -> None:
    """窗口含 N 时逐前缀判定:N 之前的干净前缀仍计数。

    Count each prefix separately when the window contains N: clean prefixes
    (no N in their own span) still count.
    """
    for k in count_ks:
        if k > k_eff or k not in arrays:
            continue
        if 255 in bs[i:i + k]:
            continue
        c = canonical_code(code >> (2 * (k_eff - k)), k)
        arrays[k][c] += 1


def _maybe_record(code: int, k: int, si: int, pos: int, seed_arr: array,
                  positions: dict) -> None:
    c = canonical_code(code, k)
    if seed_arr[c] >= 2:
        positions.setdefault(c, []).append((si, pos, 1 if code == c else 0))


def _maybe_count15(code: int, gate_k: int | None, arrays: dict,
                   counts15: dict) -> None:
    if gate_k is not None:
        arr = arrays.get(gate_k)
        if arr is not None:
            # 15-mer 的末尾 gate_k 碱基(同起点前缀提取:code >> (2*(15-gate_k)))
            cg = canonical_code(code >> (2 * (15 - gate_k)), gate_k)
            if arr[cg] < 2:
                return
    c15 = canonical_code(code, 15)
    counts15[c15] = counts15.get(c15, 0) + 1


# ---------------------------------------------------------------- 剖面
# ---------------------------------------------------------------- Profiles

def combine_kmer_scores(values: list[float], weights: list[float] | None = None) -> float:
    """多尺度 combine(§17 独立函数)。

    取 max —— "该位置最强特异性证据"。与文档的加权和不同,原因:
    同一开始位置的 k-mer 分数随 k 单调不减(更长 k-mer 唯一 ⟹ 其同位置
    子串必唯一),加权和会被短 k-mer 的高重复计数压低到阈值之下(数学
    必然);max 在保持"长 k-mer 占优"的同时保证 0.95 级阈值可达。
    weights 参数保留用于未来调参(若改加权和需同步调低阈值)。

    Multi-scale combine (§17, independent function). Uses max — "the
    strongest specificity evidence at the position". Differs from the spec's
    weighted sum: k-mer scores at the same start are monotone non-decreasing
    in k (a unique longer k-mer implies unique substrings at the same
    position), so a weighted sum is mathematically forced below the 0.95
    thresholds by short-k repeat counts; max keeps long-k dominance while
    keeping the thresholds reachable. The weights parameter is reserved for
    future tuning (switching to a weighted sum requires lowering the
    thresholds).
    """
    if not values:
        return 0.0
    return max(values)


def _covering_count(j: int, k: int, n: int) -> int:
    """覆盖 0-based 位置 j 的 k-mer 起点个数(差分数组除法用)。"""
    return min(j + 1, k, n - j)


def compute_profiles(
    index: KmerIndex,
    template: str,
    kmer_set: tuple[int, ...] | None = None,
    three_prime_windows: tuple[int, ...] = (8, 10, 12, 15),
    on_progress: Callable[[float], None] | None = None,
    cancel: CancelFlag | None = None,
) -> dict:
    """逐碱基特异性剖面(§16-§19)。

    返回:
      per_k:      {k: [float]*n}  各尺度剖面(差分数组 O(N) 构建;j 处 =
                 覆盖 j 的全部 k-mer 分数均值)
      global:     [float]*n       combine_kmer_scores(max) 合并
      three_prime:[float]*n       3' 端剖面:窗口 {8,10,12,15} 以 j 结尾的
                 k-mer 分数取 max(模板边缘窗口不足时为 0,保守)

    3' 窗口含 15-mer 的理由(实测):真实基因组(≥2 Mb)中 5/8/10-mer 的
    期望命中数 = 2N/4^k,10-mer 在 26 Mb 库中即 ~50 次 → 分数塌到 ~0.15,
    3' 阈值(Level 1: 0.95)永远不可达,四等级设计退化且预筛误杀。加入
    15-mer(26 Mb 库期望 0.008,唯一概率 >99%)后 3' 分数恢复可达性;
    短窗口保留,使重复区 3' 分仍低(区分度来自最长的唯一窗口)。

    Per-base specificity profiles (§16-§19).

    Returns:
      per_k:      {k: [float]*n}  per-scale profile (built O(N) via
                 difference arrays; at j = mean score of all k-mers
                 covering j)
      global:     [float]*n       combined via combine_kmer_scores (max)
      three_prime:[float]*n       3'-end profile: max over windows
                 {8,10,12,15} ending at j (0 at template edges where the
                 window is short). The 15-mer scale is included because on
                 real genomes (>=2 Mb) shorter seeds hit dozens of times
                 and the 3' thresholds would be unreachable.
    """
    n = len(template)
    ks = tuple(k for k in (kmer_set or index.kmer_set) if k <= n)
    per_k: dict[int, list[float]] = {}
    for ki, k in enumerate(ks):
        if cancel is not None and cancel.cancelled:
            raise blast.BlastError("任务已取消")
        diff = [0.0] * (n + 1)
        cnt = index.counts.get(k)
        for i in range(n - k + 1):
            code, ok = encode_string(template[i:i + k])
            if not ok:
                continue  # N:不贡献证据(score 0)
            c = canonical_code(code, k)
            count = cnt.get(c, 0) if isinstance(cnt, dict) else (cnt[c] if cnt is not None else 0)
            s = kmer_occurrence_score(count)
            diff[i] += s
            diff[i + k] -= s
        prof = [0.0] * n
        acc = 0.0
        for j in range(n):
            acc += diff[j]
            c = _covering_count(j, k, n)
            prof[j] = round(acc / c, 4) if c else 0.0
        per_k[k] = prof
        if on_progress:
            on_progress(0.35 * (ki + 1) / len(ks))

    global_profile = [
        combine_kmer_scores([per_k[k][j] for k in per_k]) for j in range(n)
    ]

    three_prime = [0.0] * n
    for j in range(n):
        if cancel is not None and cancel.cancelled and (j & 0xFF) == 0:
            raise blast.BlastError("任务已取消")
        vals = []
        for w in three_prime_windows:
            if j + 1 < w:
                continue
            cnt = index.counts.get(w)
            if cnt is None:
                continue  # 窗口无计数(配置未含该尺度)→ 不贡献证据,非"唯一"
            code, ok = encode_string(template[j - w + 1:j + 1])
            if not ok:
                continue
            c = canonical_code(code, w)
            count = cnt.get(c, 0) if isinstance(cnt, dict) else cnt[c]
            vals.append(kmer_occurrence_score(count))
        three_prime[j] = round(combine_kmer_scores(vals), 4)
        if on_progress and (j & 0xFF) == 0:
            on_progress(0.6 + 0.4 * (j + 1) / n)
    return {"per_k": per_k, "global": global_profile,
            "three_prime": three_prime}
