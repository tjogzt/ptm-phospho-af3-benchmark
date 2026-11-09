# Paper Plan — 投稿 Bioinformatics 组稿方案

**工作标题（备选）**
1. *Data purity, not model capability, dominates apparent AlphaFold3 performance on phosphorylation benchmarks*
2. *AlphaFold3 does not capture phosphorylation-induced conformational change: a benchmark, a data-purity warning, and two falsified confidence hypotheses*
3. *When confidence misleads: benchmark pitfalls and a cautionary negative result for AlphaFold3 phospho-structure prediction*

**一句话贡献**：在一个经配对纯度（P1–P5）严格清洗的磷酸化结构基准上，AlphaFold3 预测的构象位移与实验位移无相关（r=+0.11, 95% CI [−0.34,+0.57], n=15），且其两类内部置信度指标（pLDDT、PAE）与预测差异的关联均由单个杠杆点（IRF3）制造——由此给出「配对纯度决定结论、内部置信度≠预测正确性」的双重方法学警示。

**目标期刊**：Bioinformatics (Oxford University Press)
**文章类型**：Original Paper（上限 7 个印刷页，约 5,000 词，不含图表；CreateSpace 双栏）
**摘要格式**：结构化 → `Motivation / Results / Availability and Implementation / Contact / Supplementary information`（官方 author-guidelines 页面**未给出摘要词数上限**；当前草稿 215 词保留）
**参考文献**：Vancouver 数字制
**数据与代码**：须存放在稳定 URL 并归档（Bioinformatics 硬性要求）
**日期**：2026-11-08（初版）｜ 2026-11-09（数字回写 + 进度更新，见 §7）
**正文预算**：≈5,000 词（Introduction→Conclusion），图表不计入
**章节数**：6（+ Availability + Supplementary）

---

## 0. 投稿前置：硬阻塞清单（必须先做，否则结论不可引用）

| # | 阻塞项 | 依据 | 影响章节 |
|---|---|---|---|
| B1 | 冻结唯一 manifest（蛋白、UniProt、位点、phos/unphos PDB+链、P1–P5 判定、效应标签） | 仓库内纯度计数存在 n=2/6/10 三套互斥口径 | 全部 |
| B2 | 修复 Cα 提取，纳入 SEP/TPO/PTR；位点窗统一 ±10 | `res.id[0]==" "` 过滤掉 HETATM 磷酸化残基 | §3.4, §4 |
| B3 | 统一 geometry 与 compare 的 PDB/链选择逻辑，全蛋白审计 | CDK7 两管线数值矛盾（1.006/32.726 vs 0.53/2.05） | §4全部、Fig 6 |
| B4 | 重写 `main()` 输出核心 r 及其 bootstrap CI；澄清每蛋白模型数（5 vs 15） | 活代码 `n_models=5`，文档称 15 | §4, Fig 2 |
| B5 | 补「同态噪声地板」对照 | 实验位移 0.25–2.51 Å 是否超出噪声未裁决 | §4, §5 |
| B6 | 重算所有终点并覆盖旧数字 | 旧结论（69% 低估、0.75 纯净配对）已撤回 | 全部 |
| B7 | 实现**磷酸化感知的 ESMFold 推理**并纳入对比线（决策已定）**[已完成]** | 现有 `scripts/esm/run_esmfold.py` 直喂原始序列、未标注磷酸化位点；ESMFold 不解析 CCD_ 修饰，须设计 phosphomimetic 或位点掩蔽方案 | §3.3, §4.2, Fig 2 |

> 裁决依据见 `/Volumes/thinkplus/PTM/PTM_BioHex_Review_Report.pdf` 第二、四节。

---

## 1. Claims–Evidence 矩阵（论文骨架）

| # | Claim | Evidence（证据文件） | 状态 | 章节 |
|---|---|---|---|---|
| C1 | 公开 PDB 磷酸化"配对"普遍不纯：同时改变复合物组成/配体/构建体范围 | `audit/PAIR_AUDIT_REPORT.md`、`audit/CONFOUNDING_ANALYSIS.md`、`pair_purity_analysis.json` | Supported（须以 B1 manifest 复算计数） | §3.2, §4.1 |
| C2 | 全队列"幅度比/低估"结论随配对纯度分层而改变，须以冻结 manifest 复算后表述 | 冻结复算（`af3_vs_pdb_summary.json` + `manifest_frozen.json`）：纯净 n=9 中位幅度比 0.46、magnitude r=+0.22 vs 混杂 n=6 中位 0.30、r=−0.24 | Supported（作为"反面教材"论点） | §4.1, §5 |
| C3 | 在严格清洗队列上，AF3 预测位移与实验位移无相关 | `af3_vs_pdb_summary.json`（r=+0.1132, 95% CI [−0.3433,+0.5744], n=15；Spearman +0.1536, p=0.5848；slope +0.0530；最大杠杆点 Src 剔除后 r=+0.2202） | Supported | §4.2, Fig 2 |
| C4 | AF3 预测位移近似常数（0.52±0.36 Å），与实验尺度（0.25–2.51 Å）脱耦 | `af3_vs_pdb_summary.json`、`noise_floor.json` | Supported | §4.2, Fig 2 |
| C5 | 位点区 pLDDT 与预测差异的关联由单个杠杆点制造（r=+0.478 → 剔除 IRF3 后 +0.214） | `performance_explained.json`、`Fig4_plddt_vs_amplitude.pdf` | Supported（作"证伪"论点） | §4.3, Fig 4 |
| C6 | ΔPAE 与预测差异的关联由单个杠杆点制造（r=+0.601 → 剔除 IRF3 后 +0.142） | `pae_metrics.json`、`Fig5_pae_leverage.pdf` | Supported | §4.3, Fig 5 |
| C7 | 配对纯度是"必要非充分"条件：纯净配对扩样后亦出现失败案例（AURKA 0.12、Src 0.10，冻结幅度比） | `af3_vs_pdb_summary.json`、`manifest_frozen.json` | Supported | §4.1, §4.4 |
| C8 | 铰链型运动是否更难预测——**未定论** | `geometry_metrics.json`（幅度比 vs hinge r=−0.398, p=0.142，不显著；对照：vs 整体形变 RMSD r=−0.602, p=0.017；vs 平均位移 r=−0.616, p=0.014） | **Unsupported → 只作探索性描述** | §4.4（明确标注不显著） |
| C9 | 第二条、仅依赖序列的 ESMFold 线（phosphomimetic）同样不追踪实验位移，且其预测位移更小（0.25±0.41 Å）→ 两条预测器、两种输入模态、两种磷酸表示法共享"非响应" | `esmfold_vs_pdb_summary.json`（r=−0.3622, 95% CI [−0.7175,+0.3079], n=15；Spearman −0.2359, p=0.3973；slope −0.1930；预测位移 0.253±0.409 Å；最大杠杆点 PAK1 剔除后 r=−0.4705）；ESMFold 与 AF3 线的共同残基/实验位移一致（range [0.25,2.51]） | Supported | §4.2, Fig 2 |

**禁止写入的表述**：任何形式的"AF3 系统性低估 X%"、"纯净配对预测良好"、"AF3 无法预测磷酸化效应"（超出证据的普遍化）。

---

## 2. 章节结构（逐节规划）

### Abstract（结构化；官方未限词数，当前 215 词）
- **Motivation**：AF3 支持显式磷酸化输入，但其对磷酸化诱导构象变化的真实能力边界未明；现有基准多受配对纯度污染。
- **Results**：清洗后的 15 组基准上 r=+0.11（95% CI [−0.34,+0.57]）；预测位移近似常数（0.52±0.36 Å）对实验跨度 0.25–2.51 Å；第二条、仅依赖序列的 ESMFold 线（phosphomimetic）同样无信息（r=−0.36）；两类置信度关联（pLDDT +0.48、ΔPAE +0.60）均由单点（IRF3）制造（剔除后 +0.21 / +0.14）；全队列"低估"随配对纯度分层而改变。
- **Availability and Implementation**：Python 脚本 + 清洗后基准数据，GitHub/Zenodo（稳定 URL + DOI），语言/OS。
- **Contact**：通讯作者邮箱（占位）。
- **Supplementary information**：补充表/图 URL。

### §1 Introduction（≈900 词，1.2 页）
- Hook：磷酸化是最普遍的调控型 PTM；AF3 引入显式修饰残基输入被视为可用其做"机制性"预测。
- Gap：① 鲜有对"配对是否纯净"的审计；② 现有评估常把复合物效应误记为磷酸化效应；③ 内部置信度被当作正确性的代理。
- One-sentence contribution（见上）。
- 贡献列表（3 条，对应 C1/C3+C4/C5+C6）：方法学审计（P1–P5 判据）；清洗队列上的零相关阴性结果；两次置信度假设证伪。
- 结果预览 + Hero figure（Fig 1）。
- 关键引文槽位（[VERIFY]）：AF3 原文、AF3 对 PTM/修饰的评估类工作、蛋白构象变化基准、Kabsch、hinge 分解。

### §2 Background and Related Work（≈450 词，可并为 §1 子节）
- 2.1 AF3 与显式修饰残基；2.2 结构变化基准与度量；2.3 内部置信度作为代理的文献现状与混淆。
- 组织方式：按"方法族/假设"而非逐篇罗列。

### §3 Materials and Methods（≈1,400 词，2 页）
- 3.1 基准构建：数据源（PDB/UniProt）、效应类型标签（local_rearrangement / allosteric / disorder_order）。
- 3.2 配对纯度门禁 P1–P5（**方法学核心**）：P1 同 UniProt、P2 真实磷酸化残基（MODRES）、P3 两态覆盖 ≥50%、P4 共有残基 ≥20 且含目标位点、P5 链组成/配体/构建体范围一致（失败 P1–P4 排除；仅失败 P5 保留并标记 `impure`）。
- 3.3 预测协议（两条对比线）：
  - **AF3**：alphafoldserver.com、UniProt 全长序列 + 核准位点、CCD_ 修饰、多种子多模型。
  - **ESMFold**：因 ESMFold 不解析 CCD_ 修饰，采用 **phosphomimetic 替换**（pS/pT→Asp、pY→Glu）作为磷酸化感知近似；实现为 HuggingFace `transformers` 的 `EsmForProteinFolding`（`facebook/esmfold_v1`，纯 PyTorch，单模型、fp32；无 MSA 输入、无多种子平均）；两条线共用同一坐标框架与共同残基集。
- 3.4 度量：Kabsch 叠合（以 PDB unphos 为共同框架、两态独立叠合）、逐残基 Cα 位移、幅度比、方向余弦、hinge 指数（分块叠合）、严格定义位点窗（±10）。
- 3.5 统计：以**蛋白**为独立单位（明确种子非生物学重复）、bootstrap CI、杠杆点诊断（Cook's D / jackknife）、多重比较处理。
- 3.6 同态噪声地板（B5）。
- 3.7 复现：脚本、随机种子、环境、归档 DOI。

### §4 Results（≈1,500 词，2 页）
- 4.1 基准产出与配对纯度：18 条入选 → 16 两态配对 → 15 有 AF3 双态模型 → 9 通过全门禁 P1–P5（6 未通过：4×P5，2×P3/P5）；按纯度分层后结论如何改变（C1、C2、C7）。
- 4.2 主结果：零相关（C3）+ 预测位移近似常数（C4）+ **AF3 vs ESMFold 双线对比**（ESMFold 同样无响应：r=−0.36、位移 0.25±0.41 Å，C9）。Fig 2。
- 4.3 两次置信度证伪：pLDDT（C5，Fig 4）与 ΔPAE + 杠杆点（C6，Fig 5）。
- 4.4 探索性：效应类型分层与 hinge 差异（C8，明确"不显著、待扩大样本"）。

### §5 Discussion（≈600 词，0.8 页）
- 解释：为何零相关（配对纯度、度量敏感性、AF3 对修饰的响应机制）——**区分"未观测到"与"不存在"**。
- 局限（诚实清单）：n=15、激酶为主、单一叠合度量；ESMFold 对比线使用 phosphomimetic 近似（pS/pT→Asp、pY→Glu）而非真实磷酸基团，且无 MSA 输入，故该线为定性对照（尽管其同样显示不追踪：r=−0.36）；未做 MD。
- 对 AF3 用户的实践建议：先做纯度审计；不要用 pLDDT/PAE 代理 PTM 正确性。

### §6 Conclusion（≈150 词）
- 重申三条贡献与边界；不给超证据承诺。

### Availability（Bioinformatics 必填）
- 代码、清洗后基准、中间结果、归档 DOI、许可、语言/OS。

---

## 3. Figure / Table 计划

| ID | 类型 | 内容 | 数据来源 | 优先级 |
|---|---|---|---|---|
| Fig 1 | Hero（流程 + 纯度漏斗） | 候选池 → P1–P5 门禁 → 可用配对；对比"未清洗 vs 清洗"的结论差异 | `audit_pairs.py` 输出 + manifest | HIGH |
| Fig 2 | 1:1 散点（双层：|向量| 与方向）+ 双线对比 | AF3 与 ESMFold 预测 vs 实验位移；叠加 r 与 bootstrap CI；常数带 | `af3_vs_pdb_summary.json` + `esmfold_vs_pdb_summary.json` | HIGH（已含 ESMFold：Panel a 三角、Panel b 菱形） |
| Fig 3 | 逐残基谱 | 按蛋白/效应类型的位移谱与位点位置 | `af3_vs_pdb_per_residue.csv` | MEDIUM |
| Fig 4 | 证伪图（v2, n=11） | 位点 pLDDT vs 幅度比；标注 n=6→11 漂移 | `FIG4_LEGEND.md` / `performance_explained.json` | HIGH |
| Fig 5 | 杠杆点图 | ΔPAE vs 幅度比，IRF3 高亮 + 剔除后 r | `pae_metrics.json` | HIGH |
| Fig 6 | 核心结果汇总 | 幅度比/方向余弦/hinge 分层（B3 已修复） | `geometry_metrics.json` | HIGH |
| Table 1 | 基准表 | 蛋白 / UniProt / 位点 / PDB / P1–P5 | manifest | HIGH |
| Table 2 | 逐蛋白指标 | 共有残基、实验/预测位移、幅度比、方向余弦、幅度 r（含 ESM 列） | `af3_vs_pdb_summary.json` + `esmfold_vs_pdb_summary.json` | HIGH |
| Table S1 | 完整统计 | 相关、CI、p、杠杆点诊断（含 ESMFold 行） | 复算产物 | 补充 |

> **图表规范**：单栏宽 86 mm；矢量 PDF（Cairo）；图内元素英文标注；图注自洽；Fig 6 标题禁用 "uncorrelated with reality"/"harder" 等超证据表述。

---

## 4. Citation 计划（槽位，全部 [VERIFY]，禁止凭记忆生成 BibTeX）

- §1：AF3 主文；AF3 应用/局限述评；磷酸化构象变化机制综述；结构变化基准（RMSD/Kabsch）。
- §2.1：AF3 显式修饰输入相关；AF3 非蛋白/修饰组分评估。
- §2.3：以 pLDDT/PAE 为置信度代理的文献；置信度被误作正确性的讨论。
- §3.2：PDB/UniProt 注释规范；MODRES；配对可比性。
- §3.4：Kabsch（1976）；hinge/分块叠合方法；构象变化度量。
- §3.5：bootstrap、Cook's D、多重比较（BH）。
- §4/§5：与本文结论一致的阴性/警示类工作（用于定位，非背书）。

---

## 5. 范围决策（已于 2026-11-08 确认）

| 决策 | 选定 | 说明 |
|---|---|---|
| A. 文章类型 | **Original Paper** | ≤7 印刷页 / ≈5,000 词 |
| B. 论文主线 | **配对纯度警示** | 标题以「数据纯度主导表观性能」为主语；C5 判据为方法学卖点 |
| C. 对比线范围 | **AF3 + ESMFold 双线** | 需完成 B7（磷酸化感知 ESMFold 推理） |
| D. 写作顺序 | **先修复再写** | 完成 B1–B7 并重算后，再据实数字起草 |

**拟定标题**：*Data purity, not model capability, dominates apparent AlphaFold3 performance on phosphorylation benchmarks*（副标题可视 ESMFold 结果调整，如加入 "*an AF3–ESMFold comparison*"）。

---

## 6. Next Steps（按「先修复再写」排序）

- [x] **B1** 构建并冻结唯一 manifest（`scripts/build_manifest.py` → `manifest_frozen.json`）
- [x] **B2** 修复 Cα 提取（纳入 SEP/TPO/PTR）+ 位点窗统一 ±10
- [x] **B3** 统一 geometry/compare 的 PDB/链选择，全蛋白一致性审计
- [x] **B4** 重写 `main()`：输出核心 r + bootstrap CI；澄清每蛋白模型数
- [x] **B5** 补同态噪声地板对照（`noise_floor.json`）
- [x] **B6** 全量重算，产出冻结结果集（Table 1/2/S1、Fig 1–6）
- [x] **B7** 实现磷酸化感知 ESMFold 推理并纳入对比线（`run_esmfold.py`：HF `EsmForProteinFolding` + phosphomimetic、MPS/fp32；全队列 30/30 PDB 无错误；`compare_esmfold_vs_pdb.py` 得 r=−0.3622、CI [−0.7175,+0.3079]、位移 0.25±0.41 Å、max leverage PAK1→−0.4705；已并入 Fig 2、Table 2、TableS1 并复编）
- [x] 创建 `paper/refs.bib`（8 条，经 Crossref 元数据核验）并完成 `main.tex` 全部 `REF_*` 引用替换
- [x] `/paper-compile` 出 PDF（`paper/main.pdf`，13 页，无 multiply-defined / undefined 标签）
- [x] 补全 `main.tex` 的 `\verify{}` 占位：**全部完成**——作者列表（7 人）、单位（1+2 双标注，全部作者）、通讯邮箱（朱涛 zhutao@tjh.tjmu.edu.cn）、ORCID、资助（NSFC 82403759/82403616）、OS/硬件（macOS/Apple Silicon MPS/transformers）、仓库 URL（https://github.com/tjogzt/ptm-phospho-af3-benchmark）、MIT 许可、**Zenodo DOI（10.5281/zenodo.24087962）**；`main.tex` 中已无 `\verify{}` 残留
- [x] ESMFold 结果并入 Fig 2、Table 2（含 ESM 列）、TableS1 后复编（`main.pdf` 13 页）
- [x] `/citation-audit` 核验引文（8/8 条经 Crossref REST API 逐条核验：DOI、期刊、卷/期/页、年份全部与 `refs.bib` 一致；正文 14 处引用上下文均被原文献支持）
- [ ] 投稿前：BioRxiv 预印本占位

---

## 7. 变更记录 / Changelog

| 日期 | 变更 | 说明 |
|---|---|---|
| 2026-11-08 | 初版 | 组稿方案（B1–B7 阻塞清单、C1–C8、6 章结构） |
| 2026-11-09 | 数字回写 | 将旧/作废数字替换为冻结产物数字：C3 r=0.0016→+0.1132（CI [−0.3433,+0.5744]）；C4 0.49±0.33→0.52±0.36 Å、0.25–2.37→0.25–2.51 Å；C5 改为 n=15 单点证伪（+0.478→+0.214）；C6 改为 +0.601→+0.142；C2 改为按冻结纯度分层（纯净 n=9 中位 0.46 vs 混杂 n=6 0.30）；C7/C8 证据文件与数字更新为 `af3_vs_pdb_summary.json` / `geometry_metrics.json`；判据命名统一为 P1–P5 |
| 2026-11-09 | 进度更新 | B1–B6 完成；B7 推理跑通（MPS + fp32，解决 Metal 缓存 root 权限与 fp16 compute_tm 崩溃）；论文引用落地并成功编译 |
| 2026-11-09 | B7 收尾 + 双线合稿 | ESMFold 全队列 15 蛋白 × 2 态推理完成（30/30 PDB 无错误）；修复 `compare_esmfold_vs_pdb.py` 的 `cohort_stats`（键名与杠杆点逻辑对齐 AF3 线）与映射方式（phosphomimetic 替换位点改按位置 1:1 映射），使两线共同残基/实验位移一致（range [0.25,2.51]）；ESM 线结果 r=−0.3622、CI [−0.7175,+0.3079]、Spearman −0.2359(p=0.3973)、位移 0.253±0.409 Å、max leverage PAK1→−0.4705；`plot_af3_vs_pdb.py` Panel b 加 ESMFold 菱形叠加；重算 Table 2（含 ESM 列）/TableS1（含 ESMFold 行）；`main.tex` 加摘要/贡献/§4 新子节/Fig2 caption/Limitations 后干净复编 13 页（无 multiply-defined/undefined，仅 0.29 pt 溢出）；核实 Bioinformatics 官方未限摘要词数（删去"≤150 词"） |
| 2026-11-09 | 引文核验 | `/citation-audit`：`refs.bib` 8 条经 Crossref REST API 逐条比对，DOI/期刊/卷期页/年份全部一致；正文 14 处 `\citep{}` 上下文均匹配原文献主题。遗留：`nam2025alphafold3` 的引用键含 "2025" 而 `year={2026}`（Crossref published-print 2026-01），仅为键名，年份字段正确 |
| 2026-11-09 | 作者/单位/基金落稿 | 按实验室要求填入 7 位作者（墨青青一作、朱涛通讯）及**双单位完整标注**（同济医院妇产科/国家临床医学研究中心；肿瘤侵袭转移教育部重点实验室/湖北省重点实验室，全部作者 1+2）；通讯邮箱 zhutao@tjh.tjmu.edu.cn、ORCID 全列（注释块存档）；资助 NSFC 82403759(C.X.)/82403616(Y.W.)；OS/硬件填 macOS/Apple Silicon MPS；修复 `article` 类 `\author` 单行 tabular 导致的长单位名溢出（改用 `\parbox`）；复编 13 页无 Overfull/Underfull、无 undefined/multiply-defined。仅剩 3 处 `\verify{}`（仓库 URL、归档 DOI、许可） |
| 2026-11-09 | 建仓发布 | 新建公开仓库 `tjogzt/ptm-phospho-af3-benchmark`（main 分支）；写 `.gitignore` 仅纳入可复现最小集（curated/fasta/metadata、frozen manifest、scripts、figures、tables、audit），排除可重建原始大数据（`data/pdb/`、`results/af3/downloads|models/`、`results/md/`、`results/esm/pdb/`）与工具态（`.aider-*`、`.codegraph/`）及内部审查报告；`git rm --cached` 解除已跟踪的 PDB/aider 文件；新增 `LICENSE`(MIT) 与 `CITATION.cff`（7 作者 + ORCID + 双单位）；README 更新 ESMFold 状态为已完成并加仓库链接；暂存 244 文件/40 MB，推送成功；补 topics。`main.tex` 仓库 URL 与 MIT 许可占位回填后复编 13 页（0 Overfull/0 undefined），仅剩 Zenodo DOI 占位 |
| 2026-11-09 | Zenodo 归档 | 经 REST API 创建 deposition 并发布：**DOI 10.5281/zenodo.24087962**（concept DOI 10.5281/zenodo.24087961，record https://zenodo.org/records/24087962）；快照 `ptm-phospho-af3-benchmark-1.0.0.tar.gz`（`git archive` HEAD，246 文件，2.38 MB，md5 8be77056…）；upload_type=software、license=MIT、7 creators（含 ORCID 与双单位）、related_identifier→GitHub 仓库。回填 `main.tex`（2 处，`\verify{}` 全部清零）、`README.md`、`CITATION.cff`（`doi:` 字段）后复编 13 页干净并推送 |
