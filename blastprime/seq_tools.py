"""序列工具:反向互补、序列清洗、6 框翻译(按 guide.md 5.3)。

Sequence utilities: reverse complement, sequence cleaning, and 6-frame
translation (per guide.md 5.3).
"""

from __future__ import annotations

from Bio.Data import CodonTable
from Bio.Seq import Seq

# 互补映射:A↔T、C↔G(U↔A),大小写各保留
# Complement mapping: A↔T, C↔G (U↔A), preserving case.
_COMPL = str.maketrans("ATGCUNatgcun", "TACGANtacgna")

# 常用遗传密码表(标准表 + 线粒体各表等,来自 Biopython 的 NCBI 表)
# Common genetic code tables (standard table plus mitochondrial tables, etc.,
# from Biopython's NCBI tables).
CODON_TABLE_IDS = [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33]


def codon_tables() -> list[dict]:
    """遗传密码表清单 [{id, name}],供前端下拉。

    List of genetic code tables [{id, name}] for the frontend dropdown.
    """
    return [
        {"id": i, "name": CodonTable.unambiguous_dna_by_id[i].names[0]}
        for i in CODON_TABLE_IDS
    ]


def reverse_complement(text: str) -> str:
    """对序列行做 A↔T、C↔G(U↔A) 的反向互补;FASTA 标题行追加 (rev_comp)。

    Reverse-complement sequence lines via A↔T, C↔G (U↔A); FASTA header lines
    get " (rev_comp)" appended.
    """
    out = []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            out.append(line)
        elif s.startswith(">"):
            out.append(s + " (rev_comp)")
        else:
            out.append(s.translate(_COMPL)[::-1])
    return "\n".join(out)


def clean_sequence(text: str) -> str:
    """序列清洗:仅保留字母、*、-,剔除数字/空格等;FASTA 标题行保留。

    Clean sequences: keep only letters, "*" and "-"; drop digits, spaces, etc.
    FASTA header lines are kept.
    """
    out = []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith(">"):
            out.append(line)
        else:
            out.append("".join(ch for ch in s if ch.isalpha() or ch in "*-"))
    return "\n".join(out)


def _parse_fasta(text: str) -> list[tuple[str, str]]:
    """解析 FASTA 为 [(name, seq)];纯序列(无标题)视为单条。

    Parse FASTA into [(name, seq)]; bare sequence (no header) is treated as a
    single record.
    """
    seqs: list[tuple[str, str]] = []
    name = ""
    buf: list[str] = []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith(">"):
            if name or buf:
                seqs.append((name or "seq1", "".join(buf)))
            name = s[1:].strip()
            buf = []
        else:
            buf.append(s)
    if name or buf:
        seqs.append((name or "query", "".join(buf)))
    return seqs


def translate6(text: str, table: int) -> dict:
    """对每条 DNA 序列做正链 3 框 + 负链 3 框翻译。

    返回 {seqs: [{name, length, frames: {f1..f3, r1..r3}}]},
    负链先反向互补再按 3 个读框翻译。

    Translate each DNA sequence in 3 forward frames + 3 reverse frames.

    Returns {seqs: [{name, length, frames: {f1..f3, r1..r3}}]}; the reverse
    strand is reverse-complemented first, then translated in 3 reading frames.
    """
    seqs = _parse_fasta(text)
    if not seqs or all(not s for _, s in seqs):
        return {"seqs": []}
    out = []
    for name, dna_raw in seqs:
        # 仅保留碱基字母并大写,U 视作 T
        # Keep only base letters, uppercase; treat U as T.
        dna = "".join(ch for ch in dna_raw.upper() if ch in "ACGTUN").replace("U", "T")
        if not dna:
            continue
        rev = dna.translate(_COMPL)[::-1]
        frames: dict[str, str] = {}
        for fname, strand in (("f1", dna), ("f2", dna[1:]), ("f3", dna[2:]),
                              ("r1", rev), ("r2", rev[1:]), ("r3", rev[2:])):
            # 截齐到 3 的倍数,避免 Biopython 对末尾不完整密码子告警
            # Truncate to a multiple of 3 to avoid Biopython warnings about a
            # trailing incomplete codon.
            strand = strand[: len(strand) - len(strand) % 3]
            frames[fname] = str(Seq(strand).translate(table=table, to_stop=False))
        out.append({"name": name, "length": len(dna), "frames": frames})
    return {"seqs": out}
