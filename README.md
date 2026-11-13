# PTM 结构效应基准（Phosphorylation Structural Effect Benchmark）

> **论文 / Paper**: *Data purity, not model capability, can dominate apparent
> AlphaFold3 performance on phosphorylation benchmarks*（投稿 Bioinformatics）
> **代码 / Code**: <https://github.com/tjogzt/ptm-phospho-af3-benchmark>
> **归档 / Archive**: Zenodo DOI [10.5281/zenodo.24087962](https://doi.org/10.5281/zenodo.24087962) · **许可 / Licence**: MIT

评估 AlphaFold3 / ESMFold 对**磷酸化诱导构象变化**的建模拟合能力。
基准以 [`data/benchmark_dataset.py`](./data/benchmark_dataset.py) 为唯一蛋白定义源，
经五项生物纯度门禁（P1–P5）审核后冻结为
[`results/manifest_frozen.json`](./results/manifest_frozen.json)（**唯一权威记录**）。

## 核心结果 / Headline result

在主分析队列（15 对：paired 且两态 AF3 模型齐全）上，AF3 预测的每蛋白平均
Cα 位移与实验观测位移**无相关**：

| 量 / Quantity | 值 |
|---|---|
| Pearson $r$ | **+0.1132**（95% CI [−0.3433, +0.5744]） |
| Spearman $\rho$ | +0.1536（$p = 0.58$） |
| 预测位移 mean ± SD | 0.521 ± 0.359 Å |
| 实验位移范围 | 0.25 – 2.51 Å |
| 剔除最大杠杆点（Src）后 $r$ | +0.2202 |
| 剔除 2 个 P3 失败对（$n=13$）后 $r$ | +0.2134 |
| 改用 direct 叠合约定后 $r$ | +0.1214 |
| ESMFold（phosphomimetic 替换线）$r$ | −0.3622 |

结论对队列定义稳健：AF3 的"近恒定约 0.5 Å"预测位移既不随实验位移变化，也仅略高于
预测器自身的模型采样噪声底 —— 即该基准上的表观性能主要由**配对纯度**而非模型能力决定。
全部数值取自冻结产物 `results/analysis/` 与 `results/tables/`，见 [Table S1](./results/tables/TableS1_stats.csv)。

## 基准规模 / Benchmark size（冻结 manifest）

| 项 | 数量 |
|---|---|
| 基准条目 benchmark_total | 18 |
| 成对结构 paired | 16 |
| 单结构 single | 2 |
| 主分析队列 primary（paired + AF3 齐全） | 15 |
| 纯净队列 pure（P1–P5 全过 + AF3 齐全） | 9 |
| 成对且 P1–P5 全过 / WARN / FAIL | 10 / 4 / 2 |
| 已移除 removed | 6 |

**五项生物纯度门禁（P1–P5）**：P1 身份核验 HIT；P2 phos 态含真实 MODRES 磷酸化
（SEP/TPO/PTR）；P3 两态 UniProt 覆盖度均 ≥ 0.5；P4 两态残基范围重叠 ≥ 20 且覆盖目标
位点；P5 两态链组成完全一致。逐蛋白判定见 [MANIFEST_FROZEN.md](./audit/MANIFEST_FROZEN.md)。

> **口径说明**：P5 只核验链组成，不含配体 / 点突变 / 修饰集差异。这些共变量已单独
> 审计并记录在 [`results/analysis/pair_confound_audit.json`](./results/analysis/pair_confound_audit.json)。
> 因此本基准的 "pure" 必须理解为 **chain-composition-matched**（P5），而**不是**
> "仅磷酸化不同"：审计显示 9 对 pure 中 7 对两态配体不同，2 对含点突变（AURKA 的
> unphos 态携带 T287D；PAK1 两态均携带 K299R）；非 pure 的 IKBKB 其 unphos 参考
> 4E3C 为 S177E/S181E 拟磷酸双突变；另有 3 个 phos 态含 AF3 未建模的额外磷酸位点；
> **严格意义（配体相同且无替换）下仅 16 对中的 2 对（ERK2、p38α）成立**。
> 历史文档（`CONFOUNDING_ANALYSIS.md` / `PURE_PAIR_REPORT.md` / `EXPANSION_REPORT.md`）中的
> n=2/6/10/11 口径均以 `manifest_frozen.json` 为准作废。

## 目录结构 / Layout

```
data/
  benchmark_dataset.py        18 蛋白唯一权威定义（paired / unphos / phos / 位点 / 链）
  pdb/                        PDB 结构
  curated/                    UniProt canonical 序列 + 位点校验 + PDB↔UniProt 映射
  metadata/  fasta/           元数据；旧版 FASTA（已知缺陷，仅追溯）
scripts/
  build_manifest.py           冻结唯一权威 manifest（B1）
  audit_pairs.py / find_*.py  配对检索与 C1–C5 门禁
  af3/  esm/  md/             AF3 作业生成与结果解析、ESMFold 推理、MD 体系准备
  analysis/                   论文全部分析与出图（见下表）
results/
  manifest_frozen.json        冻结 manifest（唯一权威）
  af3/                        作业汇总 af3_summary.csv + 解出的 CIF + 生命周期台账
  esm/                        ESMFold 汇总与结构
  analysis/                   冻结分析产物（*.json / *_per_residue.csv）
  tables/                     Table1 / Table2 / TableS1（csv + tex）
  figures/                    论文级图（PDF 矢量 + PNG 预览）
audit/                        审计报告（MANIFEST_FROZEN.md / PAIR_AUDIT_REPORT.md / …）
paper/                        main.tex + main.pdf（Bioinformatics 投稿稿）
review-stage/                 auto-review-loop 评审留档（AUTO_REVIEW.md / REVIEWER_MEMORY.md）
```

## 复现步骤 / Reproduction

```bash
# 1. 冻结基准 manifest（唯一权威，下游全部分析以此为源）
python3 scripts/build_manifest.py

# 2. 主分析：AF3 / ESMFold 预测位移 vs 实验位移 + 配对共变量审计
python3 scripts/analysis/compare_af3_vs_pdb.py
python3 scripts/analysis/compare_esmfold_vs_pdb.py
python3 scripts/analysis/audit_confounds.py

# 3. 论文表格
python3 scripts/analysis/make_tables.py

# 4. 机制分析：几何分解 / 残基谱 / 噪声底 / 置信度 / PAE
python3 scripts/analysis/geometry_decomposition.py
python3 scripts/analysis/residue_spectrum.py
python3 scripts/analysis/noise_floor.py
python3 scripts/analysis/explain_performance.py
python3 scripts/analysis/pae_analysis.py

# 5. 出图
python3 scripts/analysis/plot_purity_funnel.py
python3 scripts/analysis/plot_af3_vs_pdb.py

# 6. 编译论文（需 TinyTeX / TeX Live 的 latexmk）
cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

依赖：Python 3.12 · numpy · scipy · matplotlib · biopython（AF3 模型已随仓库归档，无需重跑
预测即可复现全部分析与出图）。

## 脚本—产物对应 / Scripts → artifacts

| 脚本 | 产物 |
|---|---|
| `scripts/build_manifest.py` | `results/manifest_frozen.json`、`audit/MANIFEST_FROZEN.md` |
| `scripts/analysis/compare_af3_vs_pdb.py` | `af3_vs_pdb_summary.json`、`af3_vs_pdb_per_residue.csv`、`projection_sensitivity.json` |
| `scripts/analysis/compare_esmfold_vs_pdb.py` | `esmfold_vs_pdb_summary.json`、`esmfold_vs_pdb_per_residue.csv` |
| `scripts/analysis/audit_confounds.py` | `pair_confound_audit.json` |
| `scripts/analysis/make_tables.py` | `results/tables/Table1_benchmark.*`、`Table2_metrics.*`、`TableS1_stats.*` |
| `scripts/analysis/geometry_decomposition.py` | `geometry_metrics.json`、Fig 6 |
| `scripts/analysis/residue_spectrum.py` | `residue_spectrum.json`、Fig 3 |
| `scripts/analysis/noise_floor.py` | `noise_floor.json`、Fig S1 |
| `scripts/analysis/explain_performance.py` | `performance_explained.json`、Fig 4 |
| `scripts/analysis/pae_analysis.py` | `pae_metrics.json`、Fig 5 |
| `scripts/analysis/plot_purity_funnel.py` | `Fig1_purity_funnel.*` |
| `scripts/analysis/plot_af3_vs_pdb.py` | `Fig2_af3_vs_pdb.*` |

## 图表索引 / Figures & tables

| 文件 | 内容 |
|---|---|
| `Fig1_purity_funnel.pdf` | 配对纯度漏斗 + 纯度分层对比（hero 图） |
| `Fig2_af3_vs_pdb.pdf` | AF3 预测位移 vs 实验位移（零相关主结果） |
| `Fig3_residue_spectrum.pdf` | 残基谱：位移幅度比随位点距离 / 幅度衰减 |
| `Fig4_plddt_vs_amplitude.pdf` | 位点区 pLDDT vs 幅度比 / 全局 pLDDT |
| `Fig5_pae_leverage.pdf` | PAE 信号与杠杆点分析 |
| `Fig6_geometry_decomposition.pdf` | 几何分解（全长三方对齐） |
| `FigS1_noise_floor.pdf` | 预测器 / 实验自洽噪声底 |
| `results/tables/Table1_benchmark.*` | 逐蛋白门禁判定与队列归属 |
| `results/tables/Table2_metrics.*` | 逐蛋白幅度比 / 方向余弦 / 幅度相关 |
| `results/tables/TableS1_stats.*` | 队列统计与稳健性（含剔除 / 投影 / 杠杆点敏感性） |

## 引用 / Citation

见 [`paper/main.tex`](./paper/main.tex) 与 Zenodo 归档页。
