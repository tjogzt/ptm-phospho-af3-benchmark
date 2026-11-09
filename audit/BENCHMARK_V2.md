# 基准 v2（修复版）— 落地记录

执行日期：2026-10-23　依据：`audit/PAIR_AUDIT_REPORT.md` + `audit/PHOS_CANDIDATES.md`
机器可读配置：**`audit/benchmark_v2.json`**

## 判定标准（五条硬性）

1. phos 结构须经 MODRES 核验**真含 SEP/TPO/PTR**
2. unphos 结构须**不含任何磷酸化修饰**
3. 每个结构**单独指定链 ID**（同一蛋白不同 PDB 链号可能不同）
4. 声明磷酸化位点须落在**结构覆盖范围**内
5. 与 UniProt 序列覆盖度 > 40%（短肽 <30 aa 豁免）

## 一、成对结构（6 组）— 从 3 组提升而来

| 蛋白 | unphos | phos | 覆盖度 | phos 位点 | 变更 |
|---|---|---|---|---|---|
| CDK2 | 1HCL:A | 1JST:A | 99% / 100% | p160 | 保留 |
| ERK2 | 1ERK:A | 2ERK:A | 99% / 97% | p183,185 | 保留 |
| Aurora-A | 1MQ4:A | 1OL5:A | 65% / 65% | p287,288 | 保留 |
| **IRF3** | 1QWT:A | **3A77:A** | 56% / 54% | p386,402 | 1T2K → 3A77 |
| **p38α** | 1WFC:A | **3PY3:A** | 94% / 97% | p180,182 | 3D83 → 3PY3 |
| **CDK7** | **7B5O:J** | **8PYR:A** | 86% / 76% | p164,170 | 7ML6+1UA2 → 8PYR+7B5O |

**关键改进**：IRF3 与 3A77 均为 189–427 同结构域；3PY3 含完整双磷酸化（pT180+pY182）；
7B5O 为 2.5 Å cryo-EM 且无磷酸化，可作 CDK7 的非磷酸化参照（原 1UA2 竟已含 pT170）。

## 二、单结构组（4 组，仅 phos 有实验结构）

| 蛋白 | phos | 覆盖度 | 位点 | 变更 |
|---|---|---|---|---|
| **c-Src** | **1YI6:A** | 51% | p335,416 | 1Y57 → 1YI6（1Y57 标题明写 UNPHOSPHORYLATED） |
| **AMPKα2** | **6B1U:A** | 80% | p108,172 | 2H6D → 6B1U（全长 AMPK） |
| **PKA-Cα** | **6QJ7:A** | 96% | p197,338 | 1ATP → 6QJ7（1.69 Å） |
| CTD | 4JXT:B | 短肽 | p1621,1628 | 保留（10 残基 CTD 肽 YSPTSPS） |

**PKA 说明**：T197 磷酸化是激酶活性所必需，1CTP 亦含 pT197，实验上不存在非磷酸化 PKA
晶体结构，故降级为单结构组。

## 三、移出基准（4 个）

| 蛋白 | 原因 |
|---|---|
| Hsp27 (HSPB1) | 7 个 entry 中 0 个含磷酸化修饰 |
| PhKγ (PHKG1) | UniProt 仅关联 1 个 entry，无磷酸化结构 |
| CFTR (P13569) | 8/57 含磷酸化，但全为 R 域肽段-14-3-3 复合物，无全长结构 |
| CheY | Asp 磷酸化 AF3 不支持；1FQW 为 BeF3 激活态模拟 |

## 四、已下载并核验的结构文件

新增 8 个：`3A77` `3PY3` `1YI6` `6B1U` `8PYR` `6QJ7` `7B5O` `3NNX`(备用)
均位于 `data/pdb/`，已通过 `scripts/audit_pairs.py` 的五条判据。

## 五、下一步

1. 用 `audit/benchmark_v2.json` 替换 `data/benchmark_dataset.py`（需你确认后执行）
2. 为 6 组配对重新生成 AF3 作业（uniprot 全序列 + 已核验位点）
3. 重跑 Fig2 位移比较（预期从 3 组扩到 6 组，统计功效显著提升）
4. 单结构组无法做 AF3-vs-实验位移验证，仅可用于 AF3 内部 phos/unphos 差异分析
