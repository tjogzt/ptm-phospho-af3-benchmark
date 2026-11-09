# 产物清单与校验和（P0 执行后，2026-10-22）

## 新增/修改的核心文件

| 文件 | 说明 | 行数/大小 |
|---|---|---|
| `scripts/rebuild_sequences.py` | P0-1 序列重建（UniProt 权威源 + 位点校验 + PDB 映射） | 新建 |
| `scripts/af3/download_af3_results.py` | P0-2b 重写：幂等 + 生命周期台账 | 覆盖原版 |
| `scripts/analysis/analyze_af3_ptm_effect.py` | P0-4 第一版统计与主图 | 新建 |
| `data/curated/uniprot_sequences.json` | 18 蛋白 UniProt canonical 全长序列 | — |
| `data/curated/sequence_audit.json` | 逐蛋白位点校验 + 长度核对 | — |
| `data/curated/pdb_uniprot_mapping.json` | PDB→UniProt 编号偏移 | — |
| `results/af3/af3_summary.csv` | 88 作业完整汇总（原 30 行 → 88 行） | 88 行 |
| `results/af3/af3_ledger.json` | 生命周期台账（含逐蛋白 Δ 指标） | — |
| `results/figures/Fig1_af3_ptm_effect.pdf` | 主图（双面板，论文级矢量） | 25.9 KB |
| `results/figures/Fig1_af3_ptm_effect_stats.csv` | 逐蛋白统计表 | 13 行 |
| `audit/FINDINGS_P0.md` | 5 项关键审计发现 | — |
| `audit/site_audit_final.json` | 位点裁决证据链（构建体→UniProt） | 13 蛋白 |
| `audit/af3_tosubmit.json` | 8 个待提交作业（位点已校正） | — |
| `README.md` | 项目总览与复现步骤 | 77 行 |
| `requirements.txt` | 依赖声明 | 15 行 |
| `.gitignore` | 排除 AppleDouble / 大文件 / 工具状态 | — |

## 清理动作

- 删除 126 个 AppleDouble（`._*`）垃圾文件（分布于 `data/`、`results/`、`.git` 外各层）
- 原始文件备份于 `audit/backup_20261022/`

## 验证记录

| 检查项 | 结果 |
|---|---|
| `download_af3_results.py` 幂等性（连跑两次） | ✓ 字节级一致（md5 `c7d468c8...`） |
| 待提交作业位点残基自检（8 作业） | ✓ 0 处不匹配 |
| 已返回作业位点裁决（13 蛋白） | 12 OK / 1 ERROR（`phkg1`） |
| CSV 覆盖率 | 88/88 作业，13 蛋白（原 30/88，5 蛋白） |
| 图表渲染 | ✓ PDF 矢量 + PNG 300 dpi，标签无重叠 |

## 核心分析结果（第一版）

| 效应类型 | n 蛋白 | ΔpLDDT (phos − unphos) |
|---|---|---|
| Local rearrangement | 6 | **+1.67 ± 2.42** |
| Allosteric | 3 | −0.86 ± 3.85 |
| Disorder→order | 3 | −0.82 ± 1.55 |
| Negative control | 1 | −4.43 |

ΔpLDDT 与 ΔpTM 相关：**r = +0.655**（n=13）

**解读（需谨慎）**：局部重排类（激酶激活环磷酸化）AF3 给出置信度**提升**，
而别构类与无序→有序类表现为**下降或持平**——与"磷酸化诱导有序化"的直觉预期相反。
该模式若在扩容队列上复现，是本项目最有价值的科学发现。

---

# P1 追加产物（ground-truth 验证与 benchmark 纠错）

| 文件 | 说明 |
|---|---|
| `scripts/analysis/compare_af3_vs_pdb.py` | AF3 预测位移 vs PDB 实验位移（序列对齐 + 多模型平均） |
| `scripts/analysis/plot_af3_vs_pdb.py` | Fig2 出图 |
| `results/figures/Fig2_af3_vs_pdb.pdf` | 双面板：位移幅度 1:1 图 + 方向/幅度一致性 |
| `results/analysis/af3_vs_pdb_summary.json` | 逐蛋白位移指标 |
| `results/analysis/af3_vs_pdb_per_residue.csv` | 逐残基位移明细 |

## 关键纠错（benchmark 层面）

| 条目 | 问题 | 处置 |
|---|---|---|
| CDK7 / 7ML6 | 7ML6 实为细菌蛋白 CALU17，非 CDK7 | 从比较中剔除，需替换 PDB |
| IRF3 / 1QWT-1T2K | 分属 C 端 TAD 与 N 端 DBD，且 1T2K 非磷酸化 | 剔除，需重找配对 |
| PYGM / 1GPA | unphos 与 phos 为同一结构 | 仅用于 AF3 内部一致性 |

## 有效配对结果（4 组，15 模型平均）

| 蛋白 | 实验位移 | AF3 位移 | 方向余弦 | 幅度 r |
|---|---|---|---|---|
| CDK2 | 2.70 Å | 0.63 Å | +0.440 | +0.382 |
| ERK2 | 1.61 Å | 1.00 Å | +0.661 | +0.580 |
| Aurora-A | 0.80 Å | 0.27 Å | +0.471 | +0.805 |
| p38α | 0.66 Å | 1.19 Å | +0.309 | +0.052 |

结论：AF3 能捕捉磷酸化引起的位移**方向**（余弦 0.31–0.66），但系统性**低估幅度**（4 组中 3 组落在 1:1 线下方）。

## 已修正的方法学 bug

1. Kabsch 旋转矩阵写反（`Vt.T@D@U.T` → `U@D@Vt`），曾产生 28 Å 虚假 RMSD
2. AF3 两态共用同一变换导致刚体位移混入，曾产生 25–39 Å 假信号
   修正后单元测试 RMSD = 5.5e-15；AF3 unphos 与晶体结构 RMSD 落在 1.0–4.3 Å 合理区间。

---

# v2 数据回收与验证（2026-10-24）

## 数据回收

| 文件 | 大小 | 作业数 | 内容 |
|---|---|---|---|
| `folds_2026_10_23_18_06.zip` | 173 MB | 22 | 6 个成对蛋白的 phos（3 种子）+ 2 个 unphos |
| `folds_2026_10_24_07_41.zip` | 182 MB | 20 | 6 个成对蛋白的 unphos（2 种子补齐）+ 4 个单结构组 phos |

合计新增 42 个 v2 作业，全部为 UniProt 全长序列 + 核准位点。

## 解析器修复（本轮）

1. **重复提交兼容**：AF3 server 对重名作业追加 `_2`/`_3` 后缀，现归一化映射回原作业名
2. **版本标记**：新增 `data_version` 列（v1/v2），合并时 **v2 优先**
3. **已移除蛋白标记**：`REMOVED:<key>` 而非 `UNKNOWN`，避免与真未知混淆
4. **CIF 分版本目录**：`models/v1/` 与 `models/v2/`，防止 v1 覆盖 v2（曾导致分析读到旧位点模型）
5. **补齐 PYGM**：v2 配置遗漏，已按 1GPA 的 SEP@14 → UniProt Ser15 补入

## AF3 预测位移 vs 实验位移（6 组成对蛋白，15 模型平均）

| 蛋白 | 共有残基 | 实验位移 | AF3 预测 | 方向余弦 | 幅度 r | top20 重合 |
|---|---|---|---|---|---|---|
| CDK2 | 293 | 2.70 Å | 0.61 Å | +0.471 | +0.393 | 3/20 |
| ERK2 | 350 | 1.61 Å | 0.99 Å | +0.658 | +0.579 | 10/20 |
| Aurora-A | 259 | 0.80 Å | 0.28 Å | +0.457 | +0.800 | 7/20 |
| IRF3 | 230 | 0.32 Å | 0.86 Å | +0.087 | +0.386 | 12/20 |
| **p38α** | 336 | 1.43 Å | 1.27 Å | **+0.690** | **+0.959** | 10/20 |
| CDK7 | 259 | 0.51 Å | 0.17 Å | +0.095 | +0.406 | 6/20 |

**主要结论**：
1. AF3 能捕捉位移**方向**（余弦 0.09–0.69，中位 0.46）
2. 系统性**低估幅度**：6 组中 4 组落在 1:1 线下方（CDK2 最明显，2.70 → 0.61 Å）
3. **p38α 表现最佳**（r=0.959, cos=0.690），是唯一幅度接近实验的体系
4. IRF3 与 CDK7 方向一致性弱（cos ≈ 0.09），可能与两者实验位移本身很小（0.32/0.51 Å）有关，处于噪声量级

## 置信度分析（v2，6 蛋白）

| 效应类型 | n | ΔpLDDT |
|---|---|---|
| Local rearrangement | 5 | **+1.90 ± 2.60** |
| Allosteric | 1 | −2.88 |

ΔpLDDT 与 ΔpTM 相关：**r = +0.814**（n=6）

## 待办

- [ ] 提交 batch_3（剩余 18 作业，主要是单结构组种子）
- [ ] 单结构组（c-Src/AMPK/PKA/CTD）缺非磷酸化实验对照，暂无法做位移验证

---

# 磁盘清理记录（2026-10-24）

## 已删除

| 对象 | 数量 | 释放 | 说明 |
|---|---|---|---|
| `results/af3/models/v1/**/*.cif` | 450 个 | 201 MB | v1 作业位点有误，分析只使用 v2 |

**可重建性**：v1 CIF 可从 `results/af3/downloads/*.zip` 重新解出：

```bash
python3 scripts/af3/download_af3_results.py --extract
```

该命令会遍历全部 zip，按 `models/<version>/<gene>/<state>/` 结构解压。

## 保留

| 对象 | 数量 | 占用 | 理由 |
|---|---|---|---|
| `results/af3/downloads/*.zip` | 13 个 | 1.1 GB | **唯一原始数据**，重新解析的唯一来源 |
| `results/af3/models/v2/**/*.cif` | 360 个 | 114 MB | 当前分析实际使用 |

## .gitignore 修正

原规则 `results/af3/*/phos/*.cif` 在目录结构改为 `models/<version>/<gene>/<state>/`（四层）后**已失效**，
导致 CIF 可能被误纳入 git。已改为递归匹配 `results/af3/models/**/*.cif`。

## 删除后验证

| 检查项 | 结果 |
|---|---|
| `compare_af3_vs_pdb.py --all-models` | ✓ 正常运行，6 组配对出结果 |
| `analyze_af3_ptm_effect.py` | ✓ 正常运行，6 蛋白 / 90 作业 |
| v2 CIF 完整性 | ✓ 360 个全部保留 |
