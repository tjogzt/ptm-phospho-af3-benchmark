# PTM 结构效应研究（Phosphorylation Structural Effect Benchmark）

> **论文 / Paper**: *Data purity, not model capability, dominates apparent AlphaFold3
> performance on phosphorylation benchmarks*（投稿 Bioinformatics）
> **代码 / Code**: <https://github.com/tjogzt/ptm-phospho-af3-benchmark>
> **归档 / Archive**: Zenodo DOI（待发布） · **许可 / Licence**: MIT

评估 AlphaFold3 / ESMFold / 分子动力学对**磷酸化修饰诱导结构变化**的建模拟合能力。
基准含 18 个蛋白：9 组成对结构（phos/unphos 均有 PDB）+ 7 组单结构 + 2 组阴性对照，
按局部重排（local rearrangement）/ 别构（allosteric）/ 无序→有序（disorder→order）三类效应标签组织。

## 当前状态（2026-10-22）

| 技术线 | 状态 | 产物 |
|---|---|---|
| AF3 结构预测 | 88/114 作业已返回（13/18 蛋白） | `results/af3/af3_summary.csv`（88 行） |
| AF3 作业 v2 | 60 个已生成；**42 个已返回**（batch 1–2） | `af3_jobs/af3_jobs_v2.json` |
| AF3-vs-实验验证 | **15 组完成**（10 纯净 + 5 新增） | `Fig2`、`Fig3` |
| 成因机制分析 | **Fig4 v2**：n=11，原"位点置信度"假设**已证伪**（r=0.31） | `Fig4_why_variable.pdf`、`audit/FIG4_HYPOTHESIS_REFUTED.md` |
| ESMFold | **已完成**（phosphomimetic 替换线，n=15） | `results/analysis/esmfold_vs_pdb_summary.json` |
| GROMACS MD | **阻断**（磷酸化拓扑缺失，见下） | — |
| 统计分析 | 第一版已完成 | `results/figures/Fig1_af3_ptm_effect.pdf` |
| **AF3 vs PDB 验证** | **已完成（4 组有效配对）** | `results/figures/Fig2_af3_vs_pdb.pdf` |

## 目录结构

```
data/
  benchmark_dataset.py     18 蛋白定义（原位点编号含成熟/结构域编号混用，见 audit/）
  pdb/                     PDB 结构（27 个）
  fasta/                   旧版 FASTA（自观测坐标提取，已知缺陷，勿用）
  curated/                 【新】UniProt canonical 全长序列 + 位点校验 + PDB↔UniProt 映射
  metadata/                benchmark 元数据
scripts/
  rebuild_sequences.py             P0-1 序列重建（UniProt 权威源）
  audit_pairs.py                   配对质量门禁（C1–C4 四项硬核验）
  find_phos_structures.py          按真实磷酸化修饰检索候选 PDB
  af3/generate_af3_jobs.py         AF3 作业生成 v2（UniProt 序列 + CCD_ 格式）
af3_jobs/
  af3_jobs_v2.json                 60 个待提交作业（可直接上传 alphafoldserver.com）
  manifest_v2.json                 作业清单与位点对照
  download_pdbs.py                 PDB 批量下载
  af3/generate_af3_jobs.py         作业生成（注意：见 audit/FINDINGS_P0.md 发现1）
  af3/download_af3_results.py      【重写】幂等解析 + 生命周期台账
  analysis/analyze_af3_ptm_effect.py  【新】第一版统计与主图
  analysis/compare_af3_vs_pdb.py      【新】AF3 预测位移 vs 实验位移（含多模型平均）
  analysis/plot_af3_vs_pdb.py         【新】Fig2 出图
  esm/run_esmfold.py               ESMFold 推理（含已知缺陷，待重写）
  md/prepare_md_systems.py         MD 体系准备（磷酸化拓扑阻断，待修复）
results/
  af3/downloads/           原始结果 zip（11 个，6.1 GB，不入 git）
  af3/models/              解出的 CIF（792 个）
  af3/af3_summary.csv      作业级置信度汇总（幂等重建，含 data_version 列）
  af3/models/v1/           v1 作业解出的 CIF（位点部分有误，仅供追溯）
  af3/models/v2/           v2 作业解出的 CIF（位点已核准，分析只用此目录）
  af3/af3_ledger.json      生命周期台账
  figures/                 论文级图（PDF 矢量 + PNG 预览）
audit/
  FINDINGS_P0.md           P0 审计发现（5 项关键问题）
  site_audit_final.json    逐作业磷酸位点裁决证据链
  af3_tosubmit.json        待提交作业（位点已校正）
  backup_20261022/         原始文件备份
```

## 复现步骤

```bash
# 1. 重建权威序列（需网络访问 UniProt）
python3 scripts/rebuild_sequences.py

# 2. 重建结果汇总（幂等，可反复运行）
python3 scripts/af3/download_af3_results.py --extract

# 3. 第一版分析与主图
python3 scripts/analysis/analyze_af3_ptm_effect.py

# 4. AF3 预测 vs 实验位移（多模型平均）
python3 scripts/analysis/compare_af3_vs_pdb.py --all-models
python3 scripts/analysis/plot_af3_vs_pdb.py
```

## 基准修复（进行中）

`scripts/find_phos_structures.py` 按"真实磷酸化修饰"重新检索了全部候选结构，
发现 6 个蛋白可替换、2 个蛋白无合适实验结构。完整方案见 **`audit/REPAIR_PLAN.md`**。

| 蛋白 | 原 phos PDB | 问题 | 替代 |
|---|---|---|---|
| IRF3 | 1T2K | DNA 结合域，非磷酸化 | **3A77**（1.8 Å, pS386） |
| P38α | 3D83 | 抑制剂复合物 | **3NNX**（2.28 Å, pS182） |
| c-Src | 1Y57 | 标题明写 UNPHOSPHORYLATED | **1YI6**（2.0 Å, pY416） |
| AMPKα2 | 2H6D | 无磷酸化 | **6B1U**（2.77 Å, pT172） |
| CDK7 | 7ML6 | 细菌蛋白 CALU17 | **8PYR**（2.15 Å, pT170） |
| PKA-Cα | 1ATP | 链应为 E | **6QJ7**（1.69 Å, pT197） |
| Hsp27 | 4MHR | 细菌 EctD + 无磷酸化 | **无候选**（0/7） |
| PhKγ | 1PHK | 无磷酸化 | **无候选**（仅 1 entry） |
| CFTR | 5UAK | 标题明写 DEPHOSPHORYLATED | 仅有 R 域肽段-14-3-3 复合物 |

## 基准 v2（已落地）

`data/benchmark_dataset.py` 已更新为 v2（v1 备份于 `audit/backup_20261022/`），
说明见 `audit/BENCHMARK_V2.md`。

- **成对结构 10 组（数据齐全）**：CDK2, ERK2, Aurora-A, IRF3, p38α, CDK7, c-Src, AMPKα2, PYGM, CDK1
- **单结构组 2 组**：PKA-Cα(6QJ7), CTD(4JXT) —— 无可用非磷酸化对照
- **移出 6 个**：Hsp27, PhKγ, CFTR, CheY, Lysozyme, Ubiquitin（前 4 个无合适结构；后 2 个为对照，保留记录）
- **AF3 作业 v2**：72 个（12 蛋白 × 2 状态 × 3 种子），已提交 42 个并返回，**32 个待补**
  - ✅ `required_paired.json`（22 个）已提交并返回
  - `optional_single.json`（10 个）仍未提交 — 单结构组，无实验对照，可选
- **数据版本分离**：`results/af3/models/v1/` 与 `models/v2/`，分析只用 v2（位点已核准）

判定五条：phos 须真含 SEP/TPO/PTR；unphos 须无磷酸化；链 ID 逐一指定；
位点在覆盖范围内；与 UniProt 覆盖度 >40%（短肽豁免）。

## 配对质量门禁（必读）

`python3 scripts/audit_pairs.py` 对 18 个基准蛋白做四项硬核验（身份/真实磷酸化/序列/结构域重叠）。
**最新结果：18 个中仅 6 个的"磷酸化结构"真实含 SEP/TPO/PTR。**

- 完整报告：`audit/PAIR_AUDIT_REPORT.md`
- 机器可读：`audit/pair_audit.json`
- **可用于 AF3-vs-实验验证的配对仅 3 组**：CDK2(1HCL/1JST)、ERK2(1ERK/2ERK)、Aurora-A(1MQ4/1OL5)

任何新增配对本必须通过该门禁才能入库。

**判据 C5（配对纯度）**：两态链组成须一致。当前 10 组中仅 ERK2、p38α 为纯净配对，
其余 8 组混入复合物/配体差异——这显著影响结论（见 `audit/CONFOUNDING_ANALYSIS.md`）。

## 已知关键问题（必须知晓）

1. **存在两套互斥的作业 JSON**：`af3_batch_jobs.json`（88 job，**实际提交**）vs
   `af3_batch_jobs_part1+part2.json`（114 job，位点错误，从未提交）。
   证据见 `audit/FINDINGS_P0.md` 发现 1。
2. **位点编号体系混用**：benchmark 的 `phos_sites` 混用了 UniProt 全长编号、成熟蛋白编号、
   结构域相对编号三套体系。已改正 5 个（PKA/EGFR/LYZ/POLR2A/PHKG1）。
3. **MD 线阻断**：`amber99sb-ildn` 不含 SEP/TPO/PTR 参数，`pdb2gmx` 必然失败；
   需改用 CHARMM36m（CHARMM-GUI）或 ACP4MD 生成磷酸化拓扑。
4. **UBB 构建体是 685 aa 多聚泛素前体**（非 76 aa 成熟泛素）；
   **CFTR 是全长 1480 aa**（非 metadata 声明的 200 aa R 域片段）。论文表述需据此修正。
5. **多数"磷酸化结构"实为非磷酸化**：p38α(3D83) 是抑制剂复合物、c-Src(1Y57) 标题明写
   UNPHOSPHORYLATED、CFTR(5UAK) 标题明写 DEPHOSPHORYLATED、Hsp27(4MHR) 实为细菌 EctD。
   详见 `audit/PAIR_AUDIT_REPORT.md`。
6. **盲测 PDB 7ML6 不是 CDK7**，实为细菌蛋白 CALU17 —— 唯一盲测样本失效。
7. **IRF3 的 1QWT/1T2K 分属 C 端 TAD 与 N 端 DBD**，不可配对，allosteric 标签不成立。
   以上见 `audit/FINDINGS_P0.md` 发现 6/7。

## 待办（P1）

- [ ] 提交 `audit/af3_tosubmit.json` 中 8 个作业至 alphafoldserver.com
- [ ] 修复 MD 拓扑管线，跑 CDK2 / UBB / CFTR-R 三体系 × 2 状态 × 3 replicas
- [ ] 重写 `run_esmfold.py`（模型复用、序列去重、仅 unphos 参照）
- [ ] **重建基准配对集**：为 9 个结构不成立的蛋白重新检索真实磷酸化 PDB
      （以磷酸化注释为检索条件，非按蛋白名）
- [ ] 修正 PKA-Cα 链 ID（1ATP/1CTP 应为链 E）
- [ ] **提交 PKA 的 AF3 作业**（当前缺失，是 allosteric 类唯一可用配对）
- [ ] 扩充盲测子集（现有唯一样本已失效）
