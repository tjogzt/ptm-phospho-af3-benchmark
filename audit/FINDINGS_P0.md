# P0 数据审计发现（2026-10-22 执行）

## 发现 1：存在两个互斥的作业 JSON 版本 —— 实际提交的是 88-job 版

| 文件 | job 数 | src 位点 | mapk1 位点 | phkg1 位点 | 状态 |
|---|---|---|---|---|---|
| `af3_batch_jobs.json` | 88 | build 335 (→UniProt 419) ✓ | 185,187 ✓ | 2 ✗ | **实际提交（服务端存档证实）** |
| `af3_batch_jobs_part1+part2.json` | 114 | 416 ✗ | 183,185 ✗ | 10 ✗ | 未提交的错误版本 |

证据：从 11 个结果 zip 内的 `fold_*_job_request.json`（AF3 服务端返回的提交存档）逐条提取，
与 `af3_batch_jobs.json` 一致，与 `part1/part2` 不一致。

**结论**：审查报告 2.2 节"代码基线无法复现已提交作业"的判断成立且更严重——
`part1/part2` 不仅格式不同，**位点定义本身是错误的**（用成熟蛋白编号直接当 UniProt 编号）。
若按 `part1/part2` 提交，7 个蛋白的磷酸基团会挂错残基。

## 发现 2：位点正确性终裁（基于 AF3 服务端存档 + UniProt 官方注释）

以 UniProt REST `features[type=Modified residue]` 为权威真值，用序列比对确定构建体→UniProt 偏移：

- **OK (12/13)**：aurka 288、cdk2 160、cdk7 170、cftr 660/700/737/813、hspb1 15/78/82、
  irf3 386/396、mapk1 185/187、mapk14 180/182、prkaa2 172、pygm 15、src 419（构建体 335+84）、ubb 65
- **ERROR (1/13)**：`phkg1` 提交在 UniProt[2]（残基 T，但 UniProt 未注释该位点为磷酸化；
  该蛋白唯一注释位点是 S345）。`benchmark_dataset.py` 写的是 `Thr10`，既非 2 也非 345。

## 发现 3：原 benchmark 位点编号混用三套坐标系统

`benchmark_dataset.py` 中的 `phos_sites` 存在三种来源混用：
1. UniProt 全长编号（多数，正确）
2. 成熟蛋白编号（EGFR Y845 实际为 UniProt Y869，+24 信号肽；Src Y416 实际为 Y419，+3）
3. heptad/结构域内相对编号（POLR2A CTD 的 Ser2/Ser5 是 YSPTSPS 重复内序号，
   对应 UniProt 需映射到 CTD 区 1840+ 一带）

**影响**：这些错误**未污染已提交结果**（实际提交版本恰好正确），但若按 `benchmark_dataset.py`
字面编号重新生成作业，将导致 EGFR/Src/POLR2A/LYZ/PHKG1 五个蛋白全部挂错位点。

## 发现 4：UBB 构建体确认为 685 aa（多聚泛素前体），非 76 aa 成熟泛素

UniProt P0CG48 全长 685 aa（9×76 泛素重复 + 尾部）。AF3 提交的 685 aa 与 UniProt 完全一致
（identity 100%），位点 S65 正确（第一个重复单元的 S65）。
但 `benchmark_dataset.py` 声明 length=76，且 `data/fasta/UBB_phos_5OXL.fasta` 是 448 aa 的
PINK1-泛素复合物链 —— 三方不一致。**结论**：AF3 实际预测的是多聚泛素前体，
"最小结构效应负对照"的科学表述需重写。

## 发现 5：CFTR 确认为全长 1480 aa

AF3 提交 1480 aa = UniProt P13569 全长，与 metadata 声明"R-domain 片段 200 aa"不符。
位点 660/700/737/813 全部正确（均在 R 域内）。**结论**：实验做的是全长 CFTR，
论文中的"disorder_order 类 R 域片段"表述必须改为"全长 CFTR 含 R 域"。

## 数据修复动作（本轮已执行）

1. 清理 126 个 AppleDouble 垃圾文件，建立 `.gitignore`
2. 建立 `data/curated/`：18 个蛋白的 UniProt canonical 全长序列 + 位点校验报告 +
   PDB→UniProt 偏移映射表
3. 写入 `audit/site_audit_final.json`：逐作业位点裁决证据链

---

# P1 追加发现（ground-truth 比较阶段）

## 发现 6：盲测子集 PDB ID 错误 —— 7ML6 不是 CDK7

`benchmark_dataset.py` 将 CDK7 的磷酸化激活态设为 **7ML6** 并标记为唯一盲测样本
（声称"2022 年发布，晚于 AF3 训练截断"）。实际核验 PDB 文件头：

```
HEADER    UNKNOWN FUNCTION              27-APR-21   7ML6
TITLE     STRUCTURE OF CALU17 FROM THE CALICHEAMICIN BIOSYNTHESIS PATHWAY OF
          MICROMONOSPORA ECHINOSPORA
```

7ML6 是**细菌蛋白 CALU17**，与 CDK7（人源细胞周期依赖性激酶 7）毫无关系。
序列比对证实：7ML6 chain A 与 CDK7 全长仅 49/292 偶然匹配。

**影响**：项目唯一的盲测证据失效；CDK7 的任何 phos/unphos 对比都不能用 7ML6。
需替换为真实的 CDK7 磷酸化结构（或改用其他 2022 年后发布的磷酸化 PDB 作为盲测）。

## 发现 7：IRF3 成对结构分属不同结构域，不可配对

| PDB | 实际内容 | 残基范围 | 结构域 |
|---|---|---|---|
| 1QWT | IRF3 反式激活域（自抑制态） | 189–427 | C 端 TAD |
| 1T2K | IRF3/ATF-2/JUN **DNA 结合域复合物** | 3–112 | N 端 DBD |

两者覆盖 IRF3 的**完全不重叠区域**，且 1T2K 并非磷酸化结构（标题明确为 DNA 结合域复合物）。
序列比对：AF3 IRF3（427 aa）与 1QWT 匹配 239 残基，与 1T2K 匹配 0 残基。

**影响**：IRF3 的"allosteric（别构）"标签不成立——实验参照不覆盖磷酸化位点 Ser386/Ser396
所在的 C 端区域与 N 端区域的连接关系。该蛋白必须从成对比较中剔除，或另找同结构域的
磷酸化/非磷酸化配对。

## 发现 8：AF3 对磷酸化诱导的构象变化存在系统性低估（初步）

在 4 个校验通过的成对蛋白上（CDK2 / ERK2 / Aurora-A / p38α）：

| 蛋白 | 实验位移均值 | AF3 预测位移均值 | 位移方向余弦 |
|---|---|---|---|
| CDK2 | 2.70 Å | 0.42 Å | +0.594 |
| ERK2 | 1.61 Å | 0.97 Å | +0.764 |
| Aurora-A | 0.80 Å | 0.23 Å | +0.264 |
| p38α | 0.66 Å | 1.00 Å | +0.354 |

AF3 预测的位移幅度普遍小于实验值（除 p38α），但方向一致性在 ERK2 达 0.764、CDK2 达 0.594。
提示：**AF3 能捕捉磷酸化效应的方向，但显著低估幅度**。若在扩容队列复现，是有价值的核心结论。

## 方法学教训（已修正的 bug）

初版比较脚本存在两处严重错误，均由交叉验证发现并修正：
1. **Kabsch 旋转矩阵写反**（`Vt.T@D@U.T` 应为 `U@D@Vt`）→ 曾导致 28 Å 虚假 RMSD
2. **AF3 两态共用同一变换** → AF3 每次预测坐标系任意，差值混入刚体位移，曾产生 25–39 Å 假信号

修正后单元测试 RMSD = 5.5e-15，AF3 unphos 与晶体结构 RMSD 落在 1.0–4.3 Å 的合理区间。
