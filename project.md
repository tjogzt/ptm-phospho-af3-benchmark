# project.md — PTM 磷酸化结构效应基准（唯一权威记录）

> 本文件是项目的**唯一权威记录**（single source of truth）。所有成员必须遵循本文件的
> 版本控制、更新维护与质量保障规范。数值一律取自 `results/` 下的**冻结产物**，
> 禁止手工誊写或口头传递；引用产物时须给出路径。

- **课题代号 / Project ID**：PTM-phospho-AF3-benchmark
- **负责人 / Owner**：Tao Zhu（通讯，zhutao@tjh.tjmu.edu.cn）
- **在研稿件 / Manuscript**：*Data purity, not model capability, can dominate apparent
  AlphaFold3 performance on phosphorylation benchmarks*（目标期刊 **Bioinformatics**）
- **代码仓库 / Repo**：<https://github.com/tjogzt/ptm-phospho-af3-benchmark>
- **归档 / Archive**：Zenodo DOI [10.5281/zenodo.24087962](https://doi.org/10.5281/zenodo.24087962) · MIT
- **当前版本 / Version**：v0.3-review3（2026-11-13）

---

## 目录索引 / Index

| 章节 | 关键词 |
|---|---|
| [1. 课题目的](#1-课题目的) | 研究背景 / 目标 / 预期成果 / negative result |
| [2. 详细实验路线](#2-详细实验路线) | P1–P5 门禁 / AF3 / ESMFold / 度量 / 噪声底 |
| [3. 实验结果](#3-实验结果) | 漏斗 24→9 / headline r / 分层 / 残基谱 / 评审轮次 |
| [4. 图表文件](#4-图表文件) | Fig1–6 / FigS1 / Table1·2·S1（含生成日期） |
| [5. 研究计划完成情况](#5-研究计划完成情况) | 进度表 / 未完成项 |
| [6. 版本控制规范](#6-版本控制规范) | ▲▼■⚠ 标记 / Git / 变更说明 |
| [7. 更新维护要求](#7-更新维护要求) | 24h 更新 / 变更说明栏 / 索引 |
| [8. 质量保障措施](#8-质量保障措施) | 周核查 / 会审 / 数据-代码关联 |

---

## 1. 课题目的

### 1.1 研究背景

AlphaFold3（AF3）接受显式的翻译后修饰（PTM）输入，这常被理解为它能够建模磷酸化诱导的
构象后果。然而已有基准很少审核所配对的「磷酸化 / 非磷酸化」实验结构是否**真正可比**，
因此复合物组成等协变量（oligomeric state、bound ligand、partner chain）诱导的构象差异
可能被误记为磷酸化效应。现有评估（`nam2026alphafold3`、`ramasamy2026phosphoaf3`）多聚焦
修饰残基的**局部几何**，对**整体构象响应**是否被捕捉仍不明确。

### 1.2 研究目标

1. 建立一套**生物纯度门禁（P1–P5）**，系统审计磷酸化结构对的可比性；
2. 在通过门禁的「干净」队列上检验 AF3 预测的构象响应是否追踪实验观测响应；
3. 量化 AF3 内部置信度信号（pLDDT / PAE）是否可作正确性代理；
4. 以 ESMFold（sequence-only）作为独立对照线。

### 1.3 预期成果

- 一个可复现的磷酸化结构效应基准（冻结 manifest + 分析脚本 + 中间产物）；
- 一个**负面结果（negative result）**及其方法学含义：在配对结构变化基准上，
  **配对纯度与同态噪声标定**（而非模型能力）可以主导表观性能；
- 面向 Bioinformatics 的可投稿稿件 `paper/main.tex`（review 版草稿）+ `paper/main_oup.tex`
  （OUP 官方模板投稿版，两版共享正文 `paper/sections.tex`）（+ 归档仓库与 Zenodo DOI）。

---

## 2. 详细实验路线

> 复现命令总表见 [README.md §复现步骤](./README.md#复现步骤--reproduction)；
> 脚本→产物映射见 [README.md §脚本—产物对应](./README.md#脚本产物对应--scripts--artifacts)。

### 步骤 1 — 基准构建与冻结（B1）
- 唯一蛋白定义源：`data/benchmark_dataset.py`（paired / unphos / phos / 位点 / 链）。
- `scripts/build_manifest.py` → `results/manifest_frozen.json`（**唯一权威**）+
  `audit/MANIFEST_FROZEN.md`。
- PDB 结构存于 `data/pdb/`；UniProt canonical 序列与 PDB↔UniProt 映射存于 `data/curated/`。

### 步骤 2 — 五项生物纯度门禁 P1–P5
| 门禁 | 判据 |
|---|---|
| P1 identity | 两态映射到同一 UniProt 条目且序列一致性高 |
| P2 genuine MODRES | phos 态含至少一个 MODRES 记录的磷酸化残基（SEP/TPO/PTR） |
| P3 two-state coverage | 两态各自覆盖 UniProt 序列 ≥ 50% |
| P4 overlap | 两态共享 ≥ 20 个残基且共享区覆盖目标位点 |
| P5 composition | 两态链组成一致（链数/伙伴） |

P1/P2 失败者**直接剔除**；P3/P4/P5 记为质量 flag，主队列 = 通过 P1–P2 且两态 AF3 模型
齐全者（15 对），P3 失败者另报表 S1，P5 失败者另作分层。

### 步骤 3 — 预测
- **AF3**：AlphaFold Server，全长 UniProt 序列 + CCD 指定修饰残基；每态 **15 模型**
  （5 diffusion samples × 3 random seeds）。
- **ESMFold**：本地推理（macOS / Apple Silicon / MPS）；因不解析 CCD 修饰残基，
  phos 态以 phosphomimetic 替换表示（pSer/pThr→Asp, pTyr→Glu），故该线为**定性**对照。

### 步骤 4 — 度量与主分析
- 统一坐标协议：以 **unphos 实验结构为参考系**；两态**各自独立** Kabsch 叠合（`indep`）。
- `scripts/analysis/compare_af3_vs_pdb.py` → `af3_vs_pdb_summary.json`、
  `af3_vs_pdb_per_residue.csv`、`projection_sensitivity.json`；
  `compare_esmfold_vs_pdb.py` → `esmfold_vs_pdb_summary.json`。
- 指标：per-protein mean Cα displacement、amplitude ratio、direction cosine、magnitude correlation。

### 步骤 5 — 共变量审计与稳健性
- `scripts/analysis/audit_confounds.py` → `pair_confound_audit.json`：从 PDB 头
  （TITLE/HET/HETATM/MODRES/SEQADV/DBREF）抽取配体、点突变、完整 MODRES 位点集；
  位点按**序列比对**由 PDB 编号映射到 **UniProt 编号**，且仅计入与目标蛋白同源的链。
- 分层不确定性：`purity_strata` 95% bootstrap CI + 差异置换检验（20k reshuffles, seed 7）。

### 步骤 6 — 机制分析
- `geometry_decomposition.py`（几何/铰链分解）、`residue_spectrum.py`（残基谱，5,148 残基）、
  `noise_floor.py`（同态噪声底）、`explain_performance.py`（pLDDT）、`pae_analysis.py`（ΔPAE）。

### 步骤 7 — 出图 / 制表 / 编译
- `plot_purity_funnel.py`、`plot_af3_vs_pdb.py` → `results/figures/`；
  `make_tables.py` → `results/tables/`；
- review 版草稿：`latexmk -pdf -halt-on-error main.tex`（TinyTeX，14 页）；
- OUP 投稿版（cls/bst 在 `paper/oup/`）：
  `TEXINPUTS=".:./oup:" BSTINPUTS=".:./oup:" latexmk -pdf -interaction=nonstopmode main_oup.tex`（10 页）。
- 两版正文共享 `paper/sections.tex`；改正文只需一处。

---

## 3. 实验结果

> 下表按时间顺序记录。所有数值取自冻结产物，未手工誊写。

### 3.1 基准漏斗（冻结）
`results/manifest_frozen.json`：

| 阶段 | 数量 |
|---|---|
| 筛选候选 screened | 24 |
| 因无 MODRES 磷酸残基剔除 removed | 6 |
| 冻结条目 curated | 18（paired 16 + single 2） |
| 成对结构 paired | 16 |
| 两态 AF3 模型齐全（主队列 primary） | 15 |
| 通过 P1–P5（pure 队列） | 9 |
| 成对且 P1–P5 全过 / WARN / FAIL | 10 / 4 / 2 |

### 3.2 主结果（零相关，2026-11-12 冻结）
`results/analysis/af3_vs_pdb_summary.json`：

| 量 | 值 |
|---|---|
| Pearson $r$（AF3） | **+0.1132**（95% CI [−0.3433, +0.5744]，$n=15$） |
| Spearman $\rho$ | +0.1536（$p = 0.58$） |
| slope | +0.053 Å/Å |
| 预测位移 mean ± SD | 0.521 ± 0.359 Å |
| 实验位移范围 | 0.25 – 2.51 Å（median 1.50） |
| 剔除最大杠杆点（Src）后 $r$ | +0.2202 |
| 剔除 2 个 P3 失败对（$n=13$）后 $r$ | +0.2134 |
| direct 叠合约定 $r$ | +0.1214 |
| ESMFold $r$ | −0.3622（CI [−0.7175, +0.3079]） |

**观察**：预测位移近似恒定（≈0.5 Å），不随实验位移变化；该近恒定性为**队列均值**语义，
逐蛋白不稳定（within-protein CV median 0.332；8/15 > 0.3）。

### 3.3 纯度分层（符号翻转）
- pure $r = +0.2213$（$n=9$，95% CI [−0.3895, +0.9109]）；
- impure $r = −0.2385$（$n=6$，95% CI [−1.0, +0.8369]）；
- 差异置换检验 $p = 0.4007$ → 两层各自与 0 不可区分、差异不显著；
  **混合池不稳定**，非已证的两总体差异。

### 3.4 共变量审计（`results/analysis/pair_confound_audit.json`）
- 9 对 pure 中 **7 对**配体不同；**2 对**含点突变（AURKA T287D；PAK1 K299R）；
- 严格 composition-matched（配体相同且无替换）仅 **2/16**（ERK2、p38α）；
- 额外磷酸位点 **2 对**（CDK7 pS164、IKBKB pS257）+ **1 对**异源链（AMPKα2 pS108 位于 β 亚基）。

### 3.5 残基谱与机制（`residue_spectrum.json` 等）
- 残基级：**5,148 残基 / 15 蛋白**；幅度分箱 median ratio 严格单调
  0.498 / 0.308 / 0.273 / 0.224 / 0.168 / 0.067；
- 距离分解 0.15–0.31 非单调，**不**支持铰链几何解释；
- 噪声底：7/7 有 replicate 的蛋白超越自身 floor；预测器自洽 floor max 3.033 Å（AKT1）。

### 3.6 自主评审循环（auto-review-loop，Run ID `run_20261112_882bced5`）

| 轮次 | 日期 | reviewer | Score | verdict | 处置 |
|---|---|---|---|---|---|
| Round 1 | 2026-11-12 | GLM-5.3 (zhipu) | 5.0 / 10 | not ready | W1–W10 修复 |
| Round 2 | 2026-11-13 | GLM-5.3 | 6.5 / 10 | not ready | R1–R7 修复 |
| Round 3 | 2026-11-13 | GLM-5.3 | **7.5 / 10** | **almost** | W1–W7 残余修复；停止门 `stop` |

- 后端 `llm-chat`（SiliconFlow）；executor `deepseek-v4-pro`（deepseek）
  vs reviewer `zai-org/GLM-5.3`（zhipu）→ 跨家族独立（`independence_verified: true`）。
- 停止门 `review_gate.py` → `decision: stop`（score ≥ 6 且 verdict ∈ {ready, almost}）；
  循环在 MAX_ROUNDS=4 的**第 3 轮正规终止**。
- 全过程留档：`review-stage/AUTO_REVIEW.md`（累积日志）、`REVIEWER_MEMORY.md`（跨轮 memory）、
  `ROUND_*.response.md`、`REVIEW_STATE.json`。

---

## 4. 图表文件

> 全部为矢量 PDF（首选 Cairo），预览 PNG 同目录；图内元素英文标注。
> 日期为文件生成日期；版本随 `results/` 冻结（见 §6）。

| 文件 | 内容 | 生成日期 |
|---|---|---|
| `results/figures/Fig1_purity_funnel.pdf` | 配对纯度漏斗 + 分层对比（hero） | 2026-11-12 |
| `results/figures/Fig2_af3_vs_pdb.pdf` | AF3 预测位移 vs 实验位移（零相关主结果） | 2026-11-12 |
| `results/figures/Fig3_residue_spectrum.pdf` | 残基谱：幅度比随距离/幅度衰减 | 2026-11-13 |
| `results/figures/Fig4_plddt_vs_amplitude.pdf` | 位点区 pLDDT vs 幅度比 | 2026-11-12 |
| `results/figures/Fig5_pae_leverage.pdf` | PAE 信号与杠杆点分析 | 2026-11-12 |
| `results/figures/Fig6_geometry_decomposition.pdf` | 几何分解（全长三方对齐） | 2026-11-12 |
| `results/figures/FigS1_noise_floor.pdf` | 预测器 / 实验自洽噪声底 | 2026-11-12 |
| `results/tables/Table1_benchmark.{csv,tex}` | 逐蛋白门禁判定与队列归属 | 2026-11-13 |
| `results/tables/Table2_metrics.{csv,tex}` | 逐蛋白幅度比 / 方向余弦 / 幅度相关 | 2026-11-13 |
| `results/tables/TableS1_stats.{csv,tex}` | 队列统计与稳健性 | 2026-11-13 |

被测数据产物（图/表的输入）：`results/analysis/*.json`（11 个，见 §2 步骤 4–6）。
历史版本：`results/figures/archive/`、`results/analysis/archive/`。

---

## 5. 研究计划完成情况

| # | 阶段 | 产物 | 状态 |
|---|---|---|---|
| 1 | 基准构建与冻结（B1） | `manifest_frozen.json` | ✅ 完成 |
| 2 | P1–P5 门禁审计 | `Table1_benchmark.*` | ✅ 完成 |
| 3 | AF3 预测（15 模型/态） | `results/af3/` | ✅ 完成 |
| 4 | ESMFold 对照线 | `esmfold_vs_pdb_summary.json` | ✅ 完成 |
| 5 | 主分析（零相关） | `af3_vs_pdb_summary.json` | ✅ 完成 |
| 6 | 共变量审计 | `pair_confound_audit.json` | ✅ 完成 |
| 7 | 机制分析（几何/残基谱/噪声底/置信度） | `geometry_metrics.json` etc. | ✅ 完成 |
| 8 | 出图与制表 | §4 清单 | ✅ 完成 |
| 9 | 论文初稿（自包含草稿） | `paper/main.tex`（14 页） | ✅ 完成 |
| 10 | 去 AI 化润色 | git `a183f0e`, `97c4bc1` | ✅ 完成 |
| 11 | 自主评审循环 R1–R3 | `review-stage/` | ✅ 完成（stop） |
| 12 | 投稿格式转换（OUP 官方模板） | `paper/main_oup.tex`（10 页） | ✅ 完成 |
| 13 | BioRxiv 预印本占位 | `PAPER_PLAN.md` | ⬜ 未完成 |

> 说明：用户所述 `bioinfo.cls` 为旧版；OUP 现行官方投稿类为 CTAN 的
> `oup-authoring-template.cls`（v1.5）。Bioinformatics 的官方 class 声明为
> `\documentclass[numsec,webpdf,modern,large]` + author-date 引用制（见 OUP
> `OUP_LaTeX_supported_jnls.xlsx` 第 366 行）。
>
> 已完成：`paper/main.tex`（review 版草稿）与 `paper/main_oup.tex`（OUP 投稿版）
> **共享正文** `paper/sections.tex`，两版均编译通过（14 页 / 10 页）。正文改动只需
> 维护一处，杜绝「修复不传播」。
>
> 下一步：据 `PAPER_PLAN.md` 规划 BioRxiv 预印本。

---

## 6. 版本控制规范

### 6.1 标记格式（标准）

| 标记 | 含义 |
|---|---|
| ▲ | 新增（added） |
| ▼ | 删除（removed） |
| ■ | 修改（modified） |
| ⚠ | 注意（caution / breaking） |

### 6.2 Git 规范
- 使用 Git 管理全部文档、代码、数据与产物；**每次修改必须提交并附变更说明**。
- 提交信息聚焦「为什么」；构建产物与中间物按 `.gitignore` 规则跟踪（仅保留 `main.pdf`）。
- 保留历史版本：产物级历史放 `*/archive/`；正文级历史由 Git 承担，可追溯、可回滚。

### 6.3 变更说明栏（Change log）

| 日期 | 修改人 | 标记 | 变更内容 | 原因 | 关联提交/产物 |
|---|---|---|---|---|---|
| 2026-11-13 | T. Zhu | ▲ | 创建 `project.md`（本文件） | 规则 13：建立唯一权威记录 | — |
| 2026-11-13 | T. Zhu | ■ | 评审循环收尾：`AUTO_REVIEW.md` 增 Round 3 段与运行总结；`REVIEW_STATE.json` → `completed/round=3/score=7.5/almost` | 停止门 `decision: stop`，循环终止 | `review-stage/` |
| 2026-11-13 | T. Zhu | ■ | Round 3 残余修复 W1–W7（IKBKB/PAK1 事实更正、§2.2 论证方向、Conclusion 动词、§4.3 标题、无引用断言、P2 交叉引用、model-0 例外声明） | reviewer 要求「honestly fixed」 | `paper/main.tex` |
| 2026-11-13 | T. Zhu | ▼ | `model_stability.json` 由 `results/analysis/` 移出至 `archive/` | 孤儿产物且与 `af3_vs_pdb_summary.json` 冲突 | `results/analysis/archive/` |
| 2026-11-13 | T. Zhu | ■ | `audit_confounds.py` 按 dbAccession 过滤 SEQADV | CDK1 的 partner chain（cyclin B1）被误计为突变 | `pair_confound_audit.json` |
| 2026-11-13 | T. Zhu | ■ | `make_tables.py` 增 Table 1 编号约定脚注后重跑 | UniProt vs PDB MODRES 编号易混淆 | `Table1_benchmark.*` |
| 2026-11-13 | T. Zhu | ▲ | 新建 `paper/main_oup.tex`（OUP 官方模板投稿版，10 页）与 `paper/sections.tex`（两版共享正文） | 规则 13/任务 (b)：OUP 投稿格式转换；DRY 重构杜绝「修复不传播」 | `paper/main_oup.tex`, `paper/sections.tex`, `paper/oup/` |
| 2026-11-13 | T. Zhu | ■ | `paper/main.tex` 正文改为 `\input{sections}`（与 OUP 版共享） | 同上 | `paper/main.tex`（14 页） |
| 2026-11-13 | T. Zhu | ■ | `make_tables.py` 表缩放 `\textwidth`→`\linewidth` 并重跑三表 | OUP 双栏下 `table` 浮动体宽度为 `\columnwidth`，`\textwidth` 致宽表溢出约一栏宽 | `results/tables/*.tex` |

---

## 7. 更新维护要求

1. **时效性**：实验内容变更后 **24 小时内**必须更新本文件。
2. **重大调整**：须在 §6.3 变更说明栏登记（日期 / 修改人 / 原因）。
3. **索引与关键词**：新增章节须同步 §目录索引；关键术语在正文首现处以英文原名标注。
4. **溯源**：文档变更须关联对应的**实验数据 / 代码提交记录**（在变更说明栏填「关联提交/产物」列）。
5. **禁止手工誊写数值**：一律引用冻结产物路径；若产物重跑，须同步复核本文件数值。

---

## 8. 质量保障措施

1. **周核查**：每周由项目负责人核查本文件完整性（章节齐备、变更说明栏无缺口）。
2. **关键节点会审**：基准冻结、投稿前、评审残余修复后，须团队会审确认文档与产物一致。
3. **数据-代码关联**：任何文档数值变更，须能对应到 `results/` 产物或 Git 提交。
4. **fail-closed 评审**：`auto-review-loop` 的每轮判定与停止门均留痕（`review-stage/`），
   评审不可用时不静默通过（见 `review-stage/REVIEW_UNAVAILABLE.md` 历史记录）。
5. **可复现性**：分析对冻结 manifest 确定性；随机种子固定；报告与图表可由 §2 命令一键重跑。
