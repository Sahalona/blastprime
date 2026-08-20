# BlastPrime Studio

[English](README.md) | [简体中文](README_zh.md)

> 面向生物信息学本地检索与实验设计的桌面级应用 —— 本地 BLAST 数据库构建、比对与可视化、引物/短序列分析、四段式引物设计

BlastPrime Studio 是一个**完全本地运行**的 BLAST 工作台：FastAPI 后端 + 纯静态前端，启动后自动打开浏览器。序列数据不出本机（远程比对除外），界面中英双语、浅色/深色主题即时切换，长任务后台执行、可取消、进度与日志实时推送。

---

## ✨ 功能总览

| 模块 | 能力 |
|---|---|
| **① 数据库管理** | 多 FASTA 合并建库、自动序列类型检测、历史记录、物理删除 |
| **② BLAST 比对与可视化**（一级独立功能） | 程序自动选择、名称型查询、图形摘要、结果树/详情联动、原始输出懒加载、四种导出 |
| **③ 引物/短序列分析**（≤100 bp） | GC / Tm / 3′ 端 GC 夹子 / 自互补二聚体 / BLAST 特异性 六项指标与红黄绿判定 |
| **④ 引物设计** | 全序列基因组比对 → k-mer 计数逐碱基深度 → 特异性剖面四段式分级 + primer3 → blastn-short 逆向特异性评估 |

---

## 📦 安装与启动

### 依赖

- **Python 3.10+**
- **NCBI BLAST+ 2.12+**（可选，见下方"BLAST 定位规则"）

### 源码运行（Linux / macOS / Windows）

**方式一：uv（推荐，但不是必须）**——[https://docs.astral.sh/uv/](https://docs.astral.sh/uv/)

```bash
# 1. 安装依赖（首次）
uv sync

# 2. 启动
uv run blastprime
```

**方式二：标准 Python 虚拟环境 + pip**

```bash
# 1. 创建并激活虚拟环境
python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# 2. 安装依赖（首次）
pip install -r requirements.txt

# 3. 启动（任选其一）
python -m blastprime.app         # 模块入口
python run.py                    # 顶层入口（exe 打包用的同一入口）
blastprime                       # 仅当 pip 安装过入口脚本时可用
```

启动后自动打开浏览器（默认 `http://127.0.0.1:8686`）；未自动打开时访问终端打印的地址即可。

> Linux 无图形环境（无 `DISPLAY`/`WAYLAND_DISPLAY`）时只打印 URL 不打开浏览器；测试/无头场景可用 `--no-browser`。

### 命令行参数

| 参数 | 说明 | 默认 |
|---|---|---|
| `--host` | 监听地址 | `127.0.0.1` |
| `--port` | 监听端口；被占用则自动顺延并在日志告知 | `8686` |
| `--no-browser` | 不自动打开浏览器 | 自动打开 |
| `--loglevel` | `DEBUG` / `INFO` / `WARNING` / `ERROR` | 配置值，缺省 `INFO` |
| `--logfile <path>` | 同时写入日志文件 | 配置值，缺省不写 |
| `--config <path>` | 本次运行读取指定配置文件（不持久化） | `local_blast_dbs/config.json` |

### NCBI BLAST+ 定位规则

程序按固定顺序查找 BLAST+ 可执行文件：

1. **打包模式（Windows .exe）**：`sys._MEIPASS/bin/`（PyInstaller 解包目录）→ exe 同级 `bin/`
2. **源码模式**：项目根目录 → 项目根目录下 `bin/` → 系统 `PATH`

找不到时：顶部横幅警告 + "打开安装说明"（内置弹窗）+ "手动指定 BLAST 目录"（持久化）；警告消除前建库/比对/设计按钮置灰。

### 安装 NCBI BLAST+（未随程序携带时）

全部安装包见 https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ （当前最新 2.17.0）。

**Windows**
- 运行 `ncbi-blast-2.17.0+-win64.exe` 安装程序——安装器会自动配置 PATH。
- 或解压 `ncbi-blast-2.17.0+-x64-win64.tar.gz`，在**设置 → 手动指定 BLAST 目录**指向其中的 `bin` 文件夹。
- 手动添加 PATH（示例 `C:\Program Files\NCBI\blast-2.17.0+\bin`）：设置 → 系统 → 关于 → 高级系统设置 → 环境变量 → 双击 `Path` → 新建 → 粘贴 bin 路径 → 确定。

**Linux**
```bash
sudo apt install ncbi-blast+      # Debian / Ubuntu
sudo dnf install ncbi-blast+      # Fedora / RHEL
sudo pacman -S blast              # Arch
```
或解压 `ncbi-blast-2.17.0+-x64-linux.tar.gz`，在设置页手动指定 bin 目录。

**macOS**
```bash
brew install blast                # Homebrew
```
或解压 `ncbi-blast-2.17.0+-universal-macosx.tar.gz`，在设置页手动指定 bin 目录。

在**设置页手动指定目录**适用于所有平台，无需修改系统 PATH。

---

## 🗄️ 模块一：本地 BLAST 数据库构建与管理

- **多文件合并建库**：选择多个 FASTA（.fa/.fasta/.fna/.faa/.txt）或直接粘贴文本，自动合并（不修改原文件）
- **自动类型检测**：取前 2000 碱基样本，ACGTUN 占比 ≥90% 判核酸，否则蛋白；也可手动指定
- **核酸库自动 `-parse_seqids`**（保证后续名称查询可用）
- **默认存储目录**：`local_blast_dbs/blastdb_<随机后缀>/`，重建不会覆盖旧库
- **历史记录**：列出全部已建/已导入库；浏览库信息与条目（前 500 条）
- **删除**：本程序创建的库可**物理删除**全部索引文件（.nin/.nsq/.nhr/.pin/.psq/.phr 等，含空父目录），删除前必须确认

---

## 🔍 模块二：BLAST 比对与可视化

### 比对参数

| 参数 | 说明 |
|---|---|
| 目标数据库 | 本地库下拉选择，或勾选远程比对填写远程库名（如 nr/nt） |
| 比对程序 | 自动选择：核酸×核酸→blastn、核酸×蛋白→tblastn、蛋白×核酸→blastx、蛋白×蛋白→blastp |
| E-value | 默认 10 |
| 最多命中数 | 默认 500 |
| 短序列模式 | 查询 ≤30 bp 建议勾选：blastn 加 `-task blastn-short -word_size 7`（blastp 对应 `blastp-short -word_size 2`），E-value 放宽到 1000 |

### 名称型查询

粘贴普通 FASTA 之外，支持名称型查询，自动从库中取序列切片：

```
>gene                        # 精确匹配（自动去 lcl| 前缀），模糊包含三级解析
>gene,range=100-200          # 取 100~200 位（1-based，支持 100-* / *-200）
>gene,range=1-100,database=mydb   # 指定库
```

远程模式不支持名称型查询。查询 ≤30 bp 且未勾选短序列模式时，提交前给出警告。

### 结果展示

- **图形摘要**：查询坐标轴条带图，颜色按 Score 分级 —— ≥200 红 / ≥80 品红 / ≥50 绿 / ≥40 蓝 / 其余黑；**Ctrl+滚轮横向缩放**；悬停查看详情；点击联动结果树与详情面板
- **结果树**：查询 → 命中序列 → HSP 三级展开，排序/过滤
- **比对详情**：HSP 完整比对（Score/E-value、Identity/Gaps、坐标、Strand/Frame、qseq/sseq 双行比对）
- **原始输出**：`-outfmt 0` 成对格式**懒加载**（每批 1000 行），也可导入外部原始输出重建视图
- **导出**（见 [导出格式](#导出格式)）：原始 `.txt`、比对 `.txt`、统计 `.csv`、命中序列 `.fasta`

### 远程比对警告

远程比对会把查询数据发送到 NCBI 服务器，可能耗时数小时——提交前必须确认弹窗。

---

## 🧪 模块三：引物/短序列分析（≤100 bp）

输入 ≤100 bp 的引物、探针或短序列，可结合上方比对结果做特异性评估。六项指标与判定：

| 指标 | 绿色（合格） | 黄色（风险） | 红色（不合格） |
|---|---|---|---|
| GC 含量 | 40–60% | 30–40% 或 60–70% | 其余 |
| Tm 值 | 50–65 ℃ | 45–50 或 65–70 ℃ | 其余 |
| 3′ 端 GC 夹子（末 5 bp） | 1–3 | 0 或 4 | ≥5 |
| 自互补二聚体（最大连续配对） | ≤3 | 4 | ≥5 |
| BLAST 特异性 | 仅 1 个靶标 | 潜在脱靶 | 高危脱靶 |

- **特异性判定规则**：首个 HSP（按 E-value 排序）之外的其余命中，凡 E-value ≤ 10 者按 `L ≥ 8 + 2k + 3×errors` 判**高危**、`L ≥ max(8,6) + 2k + 3×errors` 判**潜在**（k = 距查询 3′ 端的碱基偏移，errors = 错配数）；同一库条目内多处结合区亦判潜在风险
- Tm：<14 bp 用 `2×(A+T)+4×(G+C)`；较长序列用盐校正或 primer3 物理模型（界面注明来源）
- **总体判定**：任一红 → 不可用；无红有黄 → 风险/待定；全绿 → 可用

---

## 🧬 模块四：引物设计

### 输入

- **模板序列**：粘贴 FASTA，或名称型查询 `>gene,range=1-500,database=mydb,name=显示名`
- **核酸基因组库（必选）**：特异性评估的比对对象

### 四段式流程

| 阶段 | 策略 |
|---|---|
| **第一步** | 整段模板 blastn（E-value 10、max_targets 5000）→ 窗口 k-mer（8/10/12/15-mer 逐碱基）计数 → 逐碱基深度 d(i) → 特异性剖面（得分 = count^(-1/3)） |
| **Level 1** | count=1（唯一）才放行；有产出即成功 |
| **Level 2** | count 2-3 放行，与 Level 1 累积 |
| **Level 3** | count 4-6 放行，同样累积 |
| **Level 4** | count≥7 全放行 → 每级候选（candidate_count 默认 50）→ 全部引物 blastn-short（word_size 7、E-value 1000）→ 结合位点解析（全长身份 ≥80%、truncated 独立档）→ 3′ 端特异性三级排序：① 基因组唯一匹配 ② 3′ 端末 1~3 bp 错配豁免 ③ 非目标区无法成对扩增（PCR 取向 + 产物 50~4000 bp）；3′ 端 ≥3 bp 完全配对且同区可成对扩增者淘汰；全淘汰则给出失败诊断与建议 |

**可用区两侧缓冲**（贯穿各阶段，默认 8 bp）：每个连续可用区上下游各外扩为设计允许区——正向引物 5′ 端在左（左扩）、反向引物 5′ 端在右（右扩）；两条引物 3′ 端必须落在原始可用区内（特异区），5′ 端可落入非特异缓冲带；外扩裁剪到 F/R 设计范围内。

### 设计参数（默认值）

| 参数 | 默认 |
|---|---|
| Tm 范围（最小/最优/最大） | 55 / 60 / 65 ℃ |
| GC 范围 | 30% ~ 70% |
| 引物长度 | 18 ~ 25 bp |
| 产物长度 | 150 ~ 300 bp（绝对值 / 相对序列长度 / 不限制） |
| 侧翼延伸 | 150 bp（不跨库条目） |
| 最大引物二聚体 | 5 bp |
| 最大 Tm 差（上下游） | 2 ℃ |
| 3′ 端 5 bp 最大 GC 数 | 3 |
| 可用区两侧缓冲 | 8 bp |
| 每阶段候选对数 | 50 |
| 阶段二/三放行阈值(global) | 0.6934 / 0.5503 |
| 非目标扩增产物判定 | 50~4000 bp |
| 目标区命中判定缓冲 | ±50 bp |

所有参数持久化，提供"恢复默认"。

### 综合评分（0~100）

- **物理分 50%**：primer3 惩罚分归一 + Tm/GC 达标加减分 + 发夹/二聚体扣分
- **特异性分 50%**：唯一匹配 = 100、3′ 端豁免 = 80、不可成对 = 60；每个脱靶位点扣分（权重可配置）
- 降序排列输出；阶段一~三成功时特异性分按目标区单拷贝程度折算

### 附加模式

- **sgRNA**：20 bp + PAM（默认 NGG），两条链扫描，GC 40%~60%
- **单引物**：产物长度"不限制"，不进行成对扩增判定，特异性按单引物 3′ 端法则评估

### 结果展示与导出

- **匹配深度图**：逐碱基深度可视化，图例区分可用/缓冲/屏蔽/目标区/正反向引物；与结果表、引物详情双向联动
- 引物详情：特异性等级、脱靶位点、综合分公式、产物长度、正/反向序列（Tm/GC/发夹/起始位置）
- 导出：详情 `.txt`、引物序列 `.fasta`、统计 `.csv`

---

## 💾 数据与配置

```
local_blast_dbs/              # 程序数据目录
├── config.json               # 全局配置（lang/theme/loglevel/logfile/blast_bin_dir/
│                             #   data_dir/db_records/primer_params）
└── blastdb_<随机后缀>/       # 每个库一个目录（含 .gitignore 写 `*`）
```

- **config.json**：语言、主题、日志级别/文件、BLAST 目录、默认数据库存储路径（`data_dir`）、数据库历史、设计参数；损坏时自动备份为 `config.json.bak` 并以默认配置启动
- **设置页**：可改 `data_dir`（立即生效）；"导入配置/下载配置"整体导入导出（迁移/备份，无引导链）；`--config` 仅本次运行读取指定文件
- **项目文件 `.json`**：BLAST 结果或引物设计可保存为项目文件（`app: BlastPrimeStudio`、`version: 1`、`kind: blast|primer_design`），可重新加载；损坏内容拒绝加载

### 导出格式

**BLAST 统计 .csv**（UTF-8 with BOM）：

```
Query_ID, Subject_ID, Identity, E-value, Bit_Score, Subject_Length, Query_Start, Query_End
```

**引物统计 .csv**（UTF-8 with BOM）：

```
Query_ID, Pair_Index, Forward_Seq, Forward_Len, Forward_Tm, Forward_GC, Forward_Start, Forward_End,
Reverse_Seq, Reverse_Len, Reverse_Tm, Reverse_GC, Reverse_Start, Reverse_End,
Product_Len, Max_Dimer_Consec, Max_Dimer_Total, Off_Target_Sites, Specificity_Level, Composite_Score
```

（单引物模式以 `Direction` 列替代 F/R；sgRNA 模式含 `sgRNA_Seq, PAM, Strand, GC, Tm, Off_Target_Sites`。）

**命中序列 .fasta**：标题 `>subject_id_hspN [from_query:xxx score:xxx]`，每 60 字符换行。
**引物序列 .fasta**：标题 `>Q_<query>_Pair_<n>_F len=xx tm=xx.x gc=xx.x%`（R 同理）。

---

## 📦 打包分发（Windows）

使用 PyInstaller 打包为双击即用的 Windows exe，**内置 NCBI BLAST+ 全套与前端资源**：

```bat
rem 按你的环境修改 P3DIR（primer3 所在目录，可运行 python -c "import primer3,os;print(os.path.dirname(primer3.__file__))" 查询）
set P3DIR=%LOCALAPPDATA%\Programs\Python\Python314\Lib\site-packages\primer3
pyinstaller --noconfirm --clean --onefile --name BlastPrimeStudio --icon static\favicon.ico --add-data "static;static" --add-data "bin;bin" --add-data "%P3DIR%\p3helpers.cp314-win_amd64.pyd;primer3" --add-data "%P3DIR%\src\libprimer3\primer3_config;primer3\src\libprimer3\primer3_config" run.py
```

- 产物命名 `BlastPrimeStudio-<版本>-win64(.exe)`，版本号遵循语义化版本
- 打包模式下数据目录建在 exe 同级的可写位置或 `%APPDATA%/BlastPrimeStudio/`（路径显示在设置页）
- 打包验证用例：全新机器（无 Python/BLAST）双击启动 → 建库 → BLAST → 可视化；四段式引物设计全流程；中文/含空格路径运行

---

## ❓ 常见问题

**Q: 顶部横幅提示"未找到 NCBI BLAST+ 工具"**
A: 安装 BLAST+（见上文）后刷新，或点"手动指定 BLAST 目录"选择 bin 目录。

**Q: 建库/比对很慢或卡住？**
A: 长任务后台执行，可在状态栏点击"取消"；任务日志抽屉（右上"日志"）可查看实时进度。

**Q: 引物设计提示无可用区/失败？**
A: 查看失败诊断：重复区占比过高时建议换用更长、重复度更低的模板区域，或放宽 E-value / 阶段二、三放行阈值 / 产物长度 / Tm 与 GC 范围。

**Q: 端口被占用？**
A: 程序自动顺延端口并在日志中告知实际端口；也可用 `--port` 指定。

**Q: 远程比对要注意什么？**
A: 查询数据会发送到 NCBI 服务器且可能耗时数小时；提交前务必确认弹窗，建议优先使用本地库。

**Q: 数据存在哪里？**
A: 库文件在 `local_blast_dbs/`（打包模式在 exe 同级或 `%APPDATA%/BlastPrimeStudio/`），配置在 `config.json`。

---

## 致谢

- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
- [primer3-py](https://github.com/libnano/primer3-py)
- [Biopython](https://biopython.org/)
- [FastAPI](https://fastapi.tiangolo.com/)
- [DeepSeek](https://www.deepseek.com/)


