# 扩样报告：纯净配对 6 → 10 组

日期：2026-10-25　脚本：`scripts/batch_search_pairs.py`

## 一、搜索策略

对 20 个磷酸化调控蛋白（主要是激酶激活环）批量检索，判据：
- 组成键 = **UniProt ID 集合**（比分子名可靠）
- 同组成组内须同时存在 phos 与 unphos 成员
- phos 的磷酸化位点须与目标激活环位点相符
- 两态均覆盖位点、序列覆盖度 > 40%

## 二、命中结果（6 个新增）

| 蛋白 | UniProt | phos | unphos | 位点 | 覆盖度 |
|---|---|---|---|---|---|
| ERK1 | P27361 | 2ZOQ | 4QTB | T202/Y204 | 93% |
| GSK3β | P49841 | 6Y9S | 1J1B | Y216 | 84% |
| PAK1 | Q13153 | 3FXZ | 1YHW | T423 | 54% |
| AKT1 | P31749 | 4GV1 | 3O96 | T308 | 69% |
| TBK1 | Q9UHD2 | 4EUU | 4IM0 | S172 | 42% |
| IKKβ | O14920 | 4KIK | 4E3C | S177/S181 | 86% |

## 三、未命中的 14 个

Chk1、Aurora B、PLK1、ABL1、CDK5、AMPKα1、p90RSK、Chk2、Wee1、PDK1、LKB1、JNK1、JNK2、p38β
—— 均无可用的同组成配对（多为单结构或复合物组成不一致）。

## 四、当前基准构成（18 蛋白 / 16 组成对）

| 类别 | 数量 | 蛋白 |
|---|---|---|
| **纯净配对**（C5 PASS） | **10** | CDK2, ERK2, Aurora-A, p38α, c-Src, PYGM, ERK1, GSK3β, PAK1, AKT1 |
| 组分相同但寡聚不同（WARN） | 4 | IRF3, CDK7, CDK1, IKKβ |
| 序列覆盖不足（FAIL） | 2 | AMPKα2, TBK1 |
| 单结构组 | 2 | PKA-Cα, CTD |

## 五、待提交作业

`af3_jobs/batches/required_paired.json` — **36 个**：
- AKT1, GSK3β, IKKβ, ERK1, PAK1, TBK1 各 6 个（phos×3 + unphos×3）

提交后 10 组纯净配对全部可算，Fig4 相关性样本量从 6 → 10。

## 六、预期

原 Fig4：位点 pLDDT vs 幅度比，r = +0.649（n=6）
扩样后：n=10，若趋势保持，相关系数将具备统计意义（p < 0.05 需 r > 0.63）
