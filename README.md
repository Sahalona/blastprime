# BlastPrime Studio

[English](README.md) | [简体中文](README_zh.md)

> A desktop-grade application for local bioinformatics search and experimental design — local BLAST database building, alignment & visualization, primer/short-sequence analysis, and four-level primer design

BlastPrime Studio is a **fully local** BLAST workbench: a FastAPI backend with a pure static frontend that opens automatically in your browser on startup. Your sequence data never leaves your machine (except remote BLAST), the UI is bilingual (中文/English) with instant switching, light/dark/system themes, and long-running tasks execute in the background — cancellable, with live progress and logs.

---

## ✨ Feature Overview

| Module | Capability |
|---|---|
| **① Database Management** | Multi-FASTA merge & build, auto sequence-type detection, history, physical delete |
| **② BLAST Alignment & Visualization** (independent top-level feature) | Automatic program selection, name-style queries, graphic summary, linked result tree/detail, lazy-loaded raw output, four export formats |
| **③ Primer / Short-Sequence Analysis** (≤100 bp) | GC / Tm / 3′-end GC clamp / self-dimer / BLAST specificity — six metrics with red-yellow-green verdicts |
| **④ Primer Design** | Whole-template genome alignment → per-base match depth → progressive masking + primer3 → blastn-short reverse specificity evaluation, four stages to primers |

---

## 📦 Installation & Startup

### Dependencies

- **Python 3.10+**
- **NCBI BLAST+ 2.12+** (optional — see “BLAST+ Discovery Order” below)

### Run from source (Linux / macOS / Windows)

**Option A: uv (recommended, not required)** — [https://docs.astral.sh/uv/](https://docs.astral.sh/uv/)

```bash
# 1. Install dependencies (first time)
uv sync

# 2. Start
uv run blastprime
```

**Option B: standard Python virtualenv + pip**

```bash
# 1. Create and activate a virtualenv
python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# 2. Install dependencies (first time)
pip install -r requirements.txt

# 3. Start
python -m blastprime.app         # or: blastprime if the entry script is available
```

Your browser opens automatically (default `http://127.0.0.1:8686`); if it doesn't, visit the URL printed in the terminal.

> On Linux without a graphical session (no `DISPLAY`/`WAYLAND_DISPLAY`), only the URL is printed and the browser is not launched; use `--no-browser` for headless/testing scenarios.

### Command-line options

| Option | Description | Default |
|---|---|---|
| `--host` | Bind address | `127.0.0.1` |
| `--port` | Listen port; auto-increments if occupied (logged) | `8000` |
| `--no-browser` | Do not auto-open a browser | auto-open |
| `--loglevel` | `DEBUG` / `INFO` / `WARNING` / `ERROR` | config value, else `INFO` |
| `--logfile <path>` | Also write logs to a file | config value, else none |
| `--config <path>` | Custom config file path | `local_blast_dbs/config.json` |

### NCBI BLAST+ discovery order

The app locates the BLAST+ executables in a fixed order:

1. **Packaged (Windows .exe)**: `sys._MEIPASS/bin/` (PyInstaller bundle dir) → `bin/` next to the exe
2. **From source**: project root → `bin/` under the project root → system `PATH`

If not found: a banner warning + “Installation guide” (official NCBI downloads, conda/apt/brew commands) + “Set BLAST dir manually” (persisted); build/align/design buttons stay disabled until resolved.

```bash
# Installing BLAST+ (examples)
sudo apt install ncbi-blast+          # Debian/Ubuntu
brew install blast                    # macOS
conda install -c bioconda blast       # conda users
```

---

## 🗄️ Module 1: Local BLAST Database Build & Management

- **Multi-file merge & build**: pick multiple FASTA files (.fa/.fasta/.fna/.faa/.txt) or paste text; files are merged automatically (originals never modified)
- **Auto type detection**: samples the first 2000 bases — if ACGTUN ≥ 90% it's nucleotide, otherwise protein; you can also override manually
- **`-parse_seqids` added automatically** for nucleotide DBs (enables name-style queries later)
- **Default storage**: `local_blast_dbs/blastdb_<random-suffix>/`; rebuilding never overwrites an existing DB
- **History**: lists every built/imported DB; browse DB info and entries (first 500)
- **Delete**: DBs created by the app can be **physically deleted** (all index files — .nin/.nsq/.nhr/.pin/.psq/.phr etc., including empty parent dirs) after a confirmation dialog

---

## 🔍 Module 2: BLAST Alignment & Visualization

### Parameters

| Parameter | Description |
|---|---|
| Target DB | Dropdown of local DBs, or enable remote and enter a remote DB name (e.g. nr/nt) |
| Program | Auto-selected: nucleotide×nucleotide→blastn, nucleotide×protein→tblastn, protein×nucleotide→blastx, protein×protein→blastp |
| E-value | Default 10 |
| Max targets | Default 500 |
| Short-sequence mode | Recommended for queries ≤30 bp: adds `-task blastn-short -word_size 7` (blastp: `blastp-short -word_size 2`), E-value relaxed to 1000 |

### Name-style queries

Besides pasting FASTA, you can query by name — the sequence is fetched and sliced from the DB automatically:

```
>gene                        # exact match (lcl| prefix stripped), fuzzy-contains fallback
>gene,range=100-200          # slice 100–200 (1-based; 100-* and *-200 supported)
>gene,range=1-100,database=mydb   # specify the DB
```

Name-style queries are not supported in remote mode. Queries ≤30 bp without short-sequence mode trigger a warning before submission.

### Results

- **Graphic summary**: bar-chart strip along the query axis; color by Score — ≥200 red / ≥80 magenta / ≥50 green / ≥40 blue / else black; **Ctrl+wheel to zoom horizontally**; hover for details; click to link the result tree and detail panel
- **Result tree**: query → subject → HSP three-level expansion, sortable/filterable
- **Alignment detail**: full HSP view (Score/E-value, Identity/Gaps, coordinates, Strand/Frame, qseq/sseq two-line alignment)
- **Raw output**: `-outfmt 0` pairwise format **lazy-loaded** in batches of 1000 lines; you can also import external raw output to rebuild the view
- **Exports** (see [Export Formats](#export-formats)): raw `.txt`, alignment `.txt`, stats `.csv`, hit sequences `.fasta`

### Remote BLAST warning

Remote BLAST sends your query to NCBI servers and can take hours — a confirmation dialog is required before submission.

---

## 🧪 Module 3: Primer / Short-Sequence Analysis (≤100 bp)

Enter a primer, probe, or short sequence ≤100 bp (optionally combined with the alignment above for specificity). Six metrics:

| Metric | Green (OK) | Yellow (Risk) | Red (Bad) |
|---|---|---|---|
| GC content | 40–60% | 30–40% or 60–70% | otherwise |
| Tm | 50–65 ℃ | 45–50 or 65–70 ℃ | otherwise |
| 3′-end GC clamp (last 5 bp) | 1–3 | 0 or 4 | ≥5 |
| Self-complementary dimer (max consecutive pairs) | ≤3 | 4 | ≥5 |
| BLAST specificity | exactly 1 target | potential off-targets | high-risk off-targets |

- **Specificity rule**: beyond the first HSP (sorted by E-value), every hit with E-value ≤ 10 is scored — `L ≥ 8 + 2k + 3×errors` → **high-risk**, `L ≥ max(8,6) + 2k + 3×errors` → **potential** (k = base offset from the query's 3′ end, errors = mismatch count); multiple binding sites within the same DB entry also count as potential risk
- Tm: `<14 bp` uses `2×(A+T)+4×(G+C)`; longer sequences use salt-corrected or primer3 physics models (source noted in the UI)
- **Overall verdict**: any red → unusable; no red but yellow → risky/pending; all green → usable

---

## 🧬 Module 4: Primer Design
### Input

- **Template sequence**: paste FASTA, or name-style query `>gene,range=1-500,database=mydb,name=display`
- **Nucleotide genome DB (required)**: the reference for specificity evaluation

### The four levels

| Level | Strategy |
|---|---|
| **Step 1** | blastn of the whole template (E-value 10, max_targets 5000, **DUST disabled** so repeats stay visible) → per-base k-mer counting (8/10/12/15-mer windows) → depth d(i) → specificity profile (score = count^(−1/3)) |
| **Level 1** | Only count=1 (unique) positions allowed; success if any pairs are produced |
| **Level 2** | count 2–3 released, accumulating with level 1 |
| **Level 3** | count 4–6 released, same accumulation |
| **Level 4** | count≥7 all released → per-level candidates (candidate_count, default 50) → every primer checked with blastn-short (word_size 7, E-value 1000) → binding-site parsing (full-length identity ≥80%, truncated tier) → 3′-end specificity tiering: ① unique genome match ② 3′-terminal 1–3 bp mismatch exemption ③ cannot form an amplifying pair outside the target (PCR orientation + product 50–4000 bp); primers with ≥3 bp perfect 3′-end pairing that can co-amplify in the same region are eliminated; if all are eliminated, a failure diagnosis with suggestions is shown |

**Both-sides available-region buffer** (throughout, default 8 bp): every contiguous available run is extended by buffer_len on each side to form the designable region — the forward primer's 5′ end extends on the left, the reverse primer's 5′ end on the right; both 3′ ends must stay inside the original available runs (specific territory), while the 5′ end may sit in the non-specific buffer band; the extension is clipped to the F/R design ranges.

### Design parameters (defaults)

| Parameter | Default |
|---|---|
| Tm range (min/opt/max) | 55 / 60 / 65 ℃ |
| GC range | 30% – 70% |
| Primer length | 18 – 25 bp |
| Product length | 150 – 300 bp (absolute / relative to sequence length / unlimited) |
| Flank extension | 150 bp (never crosses DB entries) |
| Max primer dimer | 5 bp |
| Max Tm difference (forward/reverse) | 2 ℃ |
| Max GC in last 5 bp of 3′ end | 3 |
| Available-region side buffer | 8 bp |
| Candidates per stage | 50 |
| Stage-2/3 release threshold (global) | 0.6934 / 0.5503 |
| Off-target amplicon size range | 50–4000 bp |
| Target-region hit buffer | ±50 bp |

All parameters persist; “Restore defaults” is provided.

### Composite score (0–100)

- **Physics 60%**: normalized primer3 penalty + Tm/GC bonuses + hairpin/dimer deductions
- **Specificity 40%**: unique match = 100, 3′-end exemption = 80, unpaired = 60; each off-target site deducts points
- Sorted descending; when stages 1–3 succeed, the specificity score is scaled by the single-copy fraction of the target region

### Additional modes

- **sgRNA**: 20 bp + PAM (default NGG), scans both strands, GC 40%–60%
- **Single primer**: product length “unlimited”, no pair-amplification check, specificity evaluated by the single-primer 3′-end rule

### Results & exports

- **Match-depth visualization**: per-base depth chart with a legend distinguishing usable/buffer/masked/target/forward-reverse primer regions; two-way linking with the results table and per-primer details
- Primer details: specificity level, off-target sites, score formula, product length, forward/reverse sequences (Tm/GC/hairpin/start positions)
- Exports: details `.txt`, primer sequences `.fasta`, stats `.csv`

---

## 💾 Data & Configuration

```
local_blast_dbs/              # program data directory
├── config.json               # global config (lang/theme/loglevel/logfile/blast_bin_dir/
│                             #   db_records/primer_params)
└── blastdb_<random-suffix>/  # one directory per DB (contains .gitignore with `*`)
```

- **config.json**: language, theme, log level/file, BLAST dir, DB history, design parameters; if corrupted it is backed up as `config.json.bak` and the app starts with defaults
- **Project files `.json`**: BLAST results and primer designs can be saved as project files (`app: BlastPrimeStudio`, `version: 1`, `kind: blast|primer_design`) and reloaded; corrupted content is rejected

### Export formats

**BLAST stats .csv** (UTF-8 with BOM):

```
Query_ID, Subject_ID, Identity, E-value, Bit_Score, Subject_Length, Query_Start, Query_End
```

**Primer stats .csv** (UTF-8 with BOM):

```
Query_ID, Pair_Index, Forward_Seq, Forward_Len, Forward_Tm, Forward_GC, Forward_Start, Forward_End,
Reverse_Seq, Reverse_Len, Reverse_Tm, Reverse_GC, Reverse_Start, Reverse_End,
Product_Len, Max_Dimer_Consec, Max_Dimer_Total, Off_Target_Sites, Specificity_Level, Composite_Score
```

(Single-primer mode replaces the F/R columns with `Direction`; sgRNA mode adds `sgRNA_Seq, PAM, Strand, GC, Tm, Off_Target_Sites`.)

**Hit sequences .fasta**: headers `>subject_id_hspN [from_query:xxx score:xxx]`, wrapped at 60 characters.
**Primer sequences .fasta**: headers `>Q_<query>_Pair_<n>_F len=xx tm=xx.x gc=xx.x%` (and `_R`).

---

## 📦 Packaging for Distribution (Windows)

Build a double-clickable Windows exe with PyInstaller, **bundling the full NCBI BLAST+ toolset and frontend assets**:

```bash
pyinstaller --onefile \
  --add-data "static:static" \
  --add-data "bin:bin" \
  --name "BlastPrimeStudio-<version>-win64" \
  blastprime/__main__.py
```

- Artifact naming: `BlastPrimeStudio-<version>-win64(.exe)`, semantic versioning
- In packaged mode the data directory lives in a writable location next to the exe or in `%APPDATA%/BlastPrimeStudio/` (path shown on the Settings page)
- Packaging verification cases: launch on a clean machine (no Python/BLAST) → build DB → BLAST → visualize; full four-level primer design; run from a path containing Chinese characters/spaces

---

## 🧪 Testing

```bash
# API end-to-end (build DB → BLAST → analyze → design 3 modes → project → delete, 29 checks)
python3 /tmp/bptest/e2e.py          # server must be running (--no-browser --port 8899)

# Frontend logic (jsdom, 18 checks)
node /tmp/bptest/smoke.js
```

---

## ❓ FAQ

**Q: The banner says “NCBI BLAST+ tools not found”**
A: Install BLAST+ (see above), refresh, or click “Set BLAST dir manually” and point it at the bin directory.

**Q: Build/alignment is slow or stuck?**
A: Long tasks run in the background — click “Cancel” in the status bar; the task-log drawer (top-right “Logs”) shows live progress.

**Q: Primer design reports no usable regions / failure?**
A: Read the failure diagnosis: with a high repeat fraction, switch to a longer, less repetitive template region, or relax E-value / stage-2–3 release thresholds / product length / Tm & GC ranges.

**Q: Port already in use?**
A: The app auto-increments the port and logs the actual one; or pick one explicitly with `--port`.

**Q: What should I know about remote BLAST?**
A: Your query data is sent to NCBI servers and may take hours — confirm the dialog before submitting, and prefer local DBs when possible.

**Q: Where is my data stored?**
A: Databases live under `local_blast_dbs/` (packaged mode: next to the exe or `%APPDATA%/BlastPrimeStudio/`); settings in `config.json`.

---

## Acknowledgements

- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
- [primer3-py](https://github.com/libnano/primer3-py)
- [Biopython](https://biopython.org/)
- [FastAPI](https://fastapi.tiangolo.com/)
- [DeepSeek](https://www.deepseek.com/)
