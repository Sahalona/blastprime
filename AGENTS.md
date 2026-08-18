# CLAUDE.md — BlastPrime Studio 项目工作规范

> 本文件是持久化项目规范。**完整且唯一的规格依据是仓库根目录的 `GUIDE.md`（BlastPrime Studio 重构规范文档）**，本文件只是可执行要点摘要。一切行为以 `GUIDE.md` 为准；本摘要与 GUIDE.md 冲突时以 GUIDE.md 为准。

## 项目定位

BlastPrime Studio：面向生物信息学本地检索与实验设计的桌面级应用（FastAPI 后端 + 本地端口 + 自动开浏览器 + 纯静态前端）。三大能力：

1. **本地 BLAST 数据库构建与管理**（模块一）
2. **BLAST 比对与可视化**（模块二，一级独立功能，不引导引物设计）
3. **引物设计**（模块四：全序列基因组比对 → k-mer 计数逐碱基深度 → 特异性剖面四段式阈值分级 + primer3-py → blastn-short 逆向特异性评估）
4. 模块三：短序列/引物分析（≤100 bp 查询的物理指标 + 特异性评估）

## 技术栈与架构

- Python 3.10+（**无 tkinter、无原生扩展**）；FastAPI + uvicorn；Biopython；primer3-py
- 环境管理：`uv`（pyproject.toml / uv.lock，推荐但不必须）；亦支持标准虚拟环境 + `pip install -r requirements.txt`；入口 `uv run blastprime`
- 前端：static/ 下原生 HTML/CSS/JS，**全部资源本地打包，禁止外部 CDN**
- 双语（中/英，i18n 表即时切换）、浅/深主题（默认跟随系统）、设置持久化到 config.json
- 所有长任务后台执行、可取消、进度/日志经 SSE/WebSocket 实时推送（支持断线重连）
- 进程模型：BLAST/primer3 子进程在后台线程跑，不阻塞 HTTP；停止任务必须杀子进程并清理临时文件

## 目录结构（约定）

```
BlastPrimeStudio/
├── pyproject.toml / uv.lock
├── blastprime/          # 后端包：app.py(路由) blast.py(BLAST封装) database.py(历史)
│   ├── config.py        # 配置管理
│   ├── primer_metrics.py# primer3-py 封装
│   ├── primer_design.py # 四段式设计流程
│   ├── masking.py       # 匹配深度与屏蔽掩码
│   ├── spec_eval.py     # blastn-short 逆向特异性评估
│   └── tasks.py         # 后台任务/取消/SSE
├── static/              # index.html(数据库页) blast.html design.html css/ js/
├── bin/                 # 打包携带的 NCBI BLAST+ 可执行文件
├── local_blast_dbs/     # 程序数据目录：config.json + blastdb_<随机后缀>/（含 .gitignore 写 `*`）
└── README.md / README_zh.md
```

## BLAST 程序定位（关键规则）

固定顺序：① 用户手动指定的 `blast_bin_dir`（持久化）优先于一切候选；② 打包模式 `sys._MEIPASS/bin/` → exe 同级 `bin/`；③ 源码模式：项目根目录 → 根目录 `bin/` → 系统 PATH。
找不到时：顶部横幅警告 + "打开安装说明" + "手动指定 BLAST 目录"，相关按钮置灰。

## 数据库模块要点

- `makeblastdb` 必须带 `-parse_seqids`；超时 300s；失败弹窗展示 stderr 全文
- 多 FASTA 自动合并（不修改原文件）；序列类型自动检测：前 2000 碱基 ACGTUN 占比 ≥90% 判核酸
- 默认存储 `local_blast_dbs/blastdb_<随机后缀>/`，重建不覆盖旧库
- 历史记录字段：`prefix`（库前缀路径）+ `is_created`；载入配置时自动剔除索引文件已不存在的条目
- 物理删除仅限 `is_created=true`，删除全部索引文件（.nin/.nsq/.nhr/.pin/.psq/.phr 及 .nsd .nsi .nnd .nni .njs .nog .pog .psd .psi .pnd .pni .pjs 等），空父目录一并删除，必须先弹确认框

## BLAST 比对模块要点

- 程序自动选择：核酸库×核酸→blastn、核酸库×蛋白→tblastn、蛋白库×核酸→blastx、蛋白库×蛋白→blastp
- 短序列模式：blastn 加 `-task blastn-short -word_size 7`，blastp `-task blastp-short -word_size 2`，E-value 放宽到 1000；查询 ≤30 bp 未勾选时提交前警告
- 名称型 query：`>gene[,range=100-200|100-*|*-200][,database=库basename]`，精确→去 `lcl|`→模糊包含三级解析，取序列切片回填为 FASTA；远程模式不支持
- 图形摘要：查询坐标轴条带图，颜色按 Score 分级（≥200 红、≥80 品红、≥50 绿、≥40 蓝、其余黑）；Ctrl+滚轮横向缩放；悬停提示；点击联动结果树/详情
- 原始输出 `-outfmt 0` 懒加载（每批 1000 行）；导入校验含 `Query=` 且为成对格式，其他格式拒绝
- 导出：原始 .txt、比对 .txt、统计 .csv（UTF-8 with BOM，表头见 GUIDE.md 9.5）、命中序列 .fasta（标题 `>subject_id_hspN [from_query:xxx score:xxx]`，60 字符换行）
- 远程比对：需填写远程库名，提交前必须弹窗警告"数据离开本机、可能耗时数小时"

## 引物设计四段式流程（核心）

输入：模板序列（粘贴 FASTA 或 名称+range）+ 必选核酸基因组库 + 可选目标区域。

- **第一步**：整段模板 blastn（E-value 10、max_target_seqs 5000、tabular 带坐标和 qseq/sseq/sstrand）→ 窗口 k-mer（8/10/12/15-mer 逐碱基）由 kmer_count 纯 Python 计数 → 逐碱基深度 d(i)：≥2 重复区、=1 单拷贝、=0 无匹配
- **特异性剖面**：得分 S = count^(-1/3)（count≤1 → 1.0，count=2 → 0.794，count=8 → 0.5）；global=四尺度 max、three_prime=以 j 结尾 k-mer 分数 max；与自重复分量逐碱基 min 合并
- **四段式分级（count 分档）**：L1（1.0/1.0，count=1）→ L2（0.6934/0.6934，count 2-3）→ L3（0.5503/0.5503，count 4-6）→ L4 全放行（count≥7）；L2/L3 阈值可由高级参数"阶段二/三放行阈值"覆盖（level2/3_global_th，可视化/图例跟随）；每级可用区（global≥g 且 3'≥t）交 primer3 出候选（candidate_count 50/级，高级参数"每阶段候选对数"）→ 级内后置过滤（3' 端回约束、GC 夹子、二聚体）→ 3' 12-mer 深度预筛（>200 淘汰）→ 有产出即成功
- **可用区两侧缓冲（默认 8 bp）**：每个连续可用区上下游各外扩为设计允许区；方向性——正向引物 5' 端在左（左扩）、反向引物 5' 端在右（右扩）；两条引物 3' 端由级内后置过滤约束在原始可用区内；外扩裁剪到 F/R 设计范围（避免引物超范围）
- **逆向特异性评估**：全部引物 blastn-short（word_size 7、E-value 1000、DUST off）→ 结合位点解析（全长身份 ≥ binding_min_identity 80%、修剪端判定、有效位点 >200 截断）→ 3' 端特异性三级排序：
  1. 基因组唯一匹配；2. 3' 端最后 1~3 bp 存在错配豁免；3. 非目标区无法成对扩增（同条目正/反向位点需 PCR 取向、产物 50~4000 bp）
  - 淘汰：3' 端 ≥3 bp 完全配对 + 同区可成对扩增；truncated（命中数超上限）独立档位淘汰；全淘汰则设计失败并给失败诊断
- **综合评分 0~100**：物理分 + 特异性分（权重可配置，默认各 50%；唯一=100、3'豁免=80、不可成对=60、k-mer 深度档按剖面折算、truncated=0；脱靶位点数每个扣分），降序排列
- **跳过开关**：跳过特异性查询（skip_spec_eval，不做 blastn-short，按 k-mer 剖面评分）；跳过 k-mer 评分过程（skip_kmer_scoring，不建深度/剖面，四段式全放行，不查/不写 k-mer 缓存，特异性分 0 注明）
- 目标区内命中判定缓冲 ±50 bp（仅评估判定，不参与几何设计）
- 附加模式：**sgRNA**（20 bp + PAM 默认 NGG，两条链扫描，GC 40%~60%；同样参与 k-mer 评分与四段式分级，skip_spec 时按 3' 剖面评分）、**单引物**（产物长度"不限制"，不进行成对扩增判定）
- 设计参数默认值见 GUIDE.md 7.6 表，全部持久化，提供"恢复默认"

## 引物分析（模块三）指标

GC（40-60 绿/30-40 或 60-70 黄/其余红）、Tm（50-65 绿/45-50 或 65-70 黄/其余红）、3' 端 GC 夹子（末 5 bp，1-3 绿/0 或 4 黄/≥5 红）、自互补二聚体最大连续配对（≤3 绿/4 黄/≥5 红）。Tm：<14 bp 用 2×(A+T)+4×(G+C)，较长用盐校正或 primer3 物理模型（注明来源）。特异性：首个 HSP 外的其余 HSP 按 `L ≥ 8 + 2k + 3×errors`（E-value ≤10）判高危红、`L ≥ max(8,6) + 2k + 3×errors` 判潜在黄。总体：任一红→不可用；无红有黄→风险/待定；全绿→可用。

## 文件与配置格式

- 项目文件 .json（GUIDE.md 9.4）：`app: "BlastPrimeStudio"`、`version: 1`、`kind: blast|primer_design`、raw_output/parsed/options/primer_design 段；加载校验 app 与 version，损坏不崩溃
- config.json：lang/theme/loglevel/logfile/blast_bin_dir/db_records/primer_params；损坏时备份为 `config.json.bak` 并默认配置启动
- CSV 均 UTF-8 with BOM；引物统计表表头见 GUIDE.md 9.5

## 错误处理要点

BLAST 缺失→横幅+置灰；库文件缺失→错误弹窗含完整路径；FASTA 非法→红框+定位非法行；建库失败→stderr 全文；设计失败→失败阶段+重复度统计+建议放宽参数；任务取消→"已取消"+清理进程；配置损坏→备份+恢复默认+状态栏提示。

## 打包

PyInstaller，Windows exe 内置 bin/（BLAST+ 全套）与 static/；数据目录放 exe 同级可写位置或 %APPDATA%/BlastPrimeStudio/；产物命名 `BlastPrimeStudio-<版本>-win64(.exe)`。打包验证用例见 GUIDE.md 11.2。

## 验收清单

GUIDE.md 第 13 章 A~E 五组勾选框为最终验收标准：通用（uv 启动/打包/定位/双语主题/任务管理）、数据库、BLAST 可视化、引物分析、引物设计（四段式、评分、联动、附加模式、导出）。

## 开发环境（当前机器）

- Linux WSL2 (Ubuntu 24.04)，uv 0.11.22，Python 3.12.3
- NCBI BLAST+ 2.12.0 已系统安装于 /usr/bin（blastn/makeblastdb/blastdbcmd 等），源码模式直接可用
- DISPLAY=:0 存在，可自动开浏览器；测试用 `--no-browser` + curl 即可
