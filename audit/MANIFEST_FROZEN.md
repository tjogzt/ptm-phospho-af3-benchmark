# 基准 Manifest（冻结版）

冻结时间：2026-11-12T05:44:30+08:00　来源：`data/benchmark_dataset.py` (sha256 `947f20a070d0ecd7`)、`audit/pair_audit.json`

> **本文件为基准唯一权威记录。** 论文 Table 1 与全部下游分析必须引用此 manifest；
> 历史文档（CONFOUNDING_ANALYSIS / PURE_PAIR_REPORT / EXPANSION_REPORT）中的 n=2/6/10/11 
> 口径均以本表为准作废。

## 计数 / Counts

| 项 | 数量 |
|---|---|
| 基准总数 benchmark_total | 18 |
| 成对 paired | 16 |
| 单结构 single | 2 |
| 主分析队列 primary（paired + AF3 齐全） | 15 |
| 纯净队列 pure（成对 + P1–P5 全过 + AF3 齐全） | 9 |
| 成对且 P1–P5 全过 / WARN / FAIL | 10 / 4 / 2 |
| 已移除 removed | 6 |

## 判据 / Criteria

- P3 覆盖度阈值 = 0.5
- P4 两态重叠阈值 = 20 残基
- P1 身份核验；P2 phos 态含 SEP/TPO/PTR；P5 两态链组成完全一致

## 逐蛋白 / Per-protein

| 蛋白 | 基因 | UniProt | 位点 | phos PDB | unphos PDB | 判定 | P1 | P2 | P3 | P4 | P5 | 主队列 | 纯净 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CDK2 | CDK2 | P24941 | 160 | 1H1R | 2C5N | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| ERK2 | MAPK1 | P28482 | 185,187 | 2ERK | 1ERK | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Aurora_A | AURKA | O14965 | 288 | 5DNR | 4J8M | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| IRF3 | IRF3 | Q14653 | 386,402 | 3A77 | 1QWT | WARN | ✓ | ✓ | ✓ | ✓ | · | ✓ | · |
| P38alpha | MAPK14 | Q16539 | 180,182 | 3PY3 | 1WFC | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| CDK7 | CDK7 | P50613 | 170 | 8PYR | 8P79 | WARN | ✓ | ✓ | ✓ | ✓ | · | ✓ | · |
| c-Src_kinase | SRC | P12931 | 419 | 9NS0 | 7NG7 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| AMPK_alpha2 | PRKAA2 | P54646 | 172 | 6B1U | 2H6D | FAIL | ✓ | ✓ | · | ✓ | · | ✓ | · |
| PKA_Calpha | PRKACA | P17612 | 198 | 6QJ7 | — | PASS | ✓ | ✓ | ✓ | · | · | · | · |
| CTD_kinase | POLR2A | P24928 | 1614,1621 | 4JXT | — | PASS | ✓ | ✓ | · | · | · | · | · |
| Glycogen_Phosphorylase | PYGM | P00489 | 15 | 1LWO | 2GJ4 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| CDK1 | CDK1 | P06493 | 161 | 9SKQ | 4YC3 | WARN | ✓ | ✓ | ✓ | ✓ | · | ✓ | · |
| ERK1 | MAPK3 | P27361 | 202,204 | 2ZOQ | 4QTB | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | · | · |
| GSK3B | GSK3B | P49841 | 216 | 6Y9S | 1J1B | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| PAK1 | PAK1 | Q13153 | 423 | 3FXZ | 1YHW | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| AKT1 | AKT1 | P31749 | 308 | 4GV1 | 3O96 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| TBK1 | TBK1 | Q9UHD2 | 172 | 4EUU | 4IM0 | FAIL | ✓ | ✓ | · | ✓ | · | ✓ | · |
| IKBKB | IKBKB | O14920 | 177,181 | 4KIK | 4E3C | WARN | ✓ | ✓ | ✓ | ✓ | · | ✓ | · |

## 队列成员 / Cohort membership

- **primary (n=15)**：AKT1, AMPK_alpha2, Aurora_A, CDK1, CDK2, CDK7, ERK2, GSK3B, Glycogen_Phosphorylase, IKBKB, IRF3, P38alpha, PAK1, TBK1, c-Src_kinase
- **pure (n=9)**：AKT1, Aurora_A, CDK2, ERK2, GSK3B, Glycogen_Phosphorylase, P38alpha, PAK1, c-Src_kinase
- **single (n=2)**：CTD_kinase, PKA_Calpha
- **fail (n=2)**：AMPK_alpha2, TBK1
- **removed (n=6)**：CFTR_R_domain, CheY, Hsp27, Lysozyme, PhK_gamma, Ubiquitin
