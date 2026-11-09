# 基准配对质量筛查报告

筛选阈值：与 UniProt 覆盖度 ≥ 50%；两态重叠 ≥ 20 残基。

## 汇总

| 蛋白 | 判定 | phos PDB | 真实磷酸化 | unphos PDB | 重叠残基 | 问题 |
|---|---|---|---|---|---|---|
| CDK2 | PASS | 1H1R | 是 | 2C5N | 296 | — |
| ERK2 | PASS | 2ERK | 是 | 1ERK | 353 | — |
| Aurora_A | PASS | 5DNR | 是 | 4J8M | 265 | — |
| IRF3 | WARN | 3A77 | 是 | 1QWT | 239 | C5_COMPOSITION |
| P38alpha | PASS | 3PY3 | 是 | 1WFC | 347 | — |
| CDK7 | WARN | 8PYR | 是 | 8P79 | 265 | C5_COMPOSITION |
| c-Src_kinase | PASS | 9NS0 | 是 | 7NG7 | 277 | — |
| AMPK_alpha2 | FAIL | 6B1U | 是 | 2H6D | 269 | C3_SEQ_MISMATCH; C5_COMPOSITION |
| PKA_Calpha | PASS | 6QJ7 | 是 | — | — | — |
| CTD_kinase | PASS | 4JXT | 是 | — | — | — |
| Glycogen_Phosphorylase | PASS | 1LWO | 是 | 2GJ4 | 824 | — |
| CDK1 | WARN | 9SKQ | 是 | 4YC3 | 286 | C5_COMPOSITION |
| ERK1 | PASS | 2ZOQ | 是 | 4QTB | 351 | — |
| GSK3B | PASS | 6Y9S | 是 | 1J1B | 349 | — |
| PAK1 | PASS | 3FXZ | 是 | 1YHW | 293 | — |
| AKT1 | PASS | 4GV1 | 是 | 3O96 | 287 | — |
| TBK1 | FAIL | 4EUU | 是 | 4IM0 | 308 | C3_SEQ_MISMATCH; C5_COMPOSITION |
| IKBKB | WARN | 4KIK | 是 | 4E3C | 648 | C5_COMPOSITION |

## 逐条问题明细

### IRF3（WARN）
- **[WARN] C5_COMPOSITION** — 两态链组成不同 unphos=['A', 'B'] vs phos=['A', 'B', 'C', 'D']（复合物/构建体差异可能混入磷酸化效应）

### CDK7（WARN）
- **[WARN] C5_COMPOSITION** — 两态链组成不同 unphos=['H', 'I', 'J'] vs phos=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']（复合物/构建体差异可能混入磷酸化效应）

### AMPK_alpha2（FAIL）
- **[FAIL] C3_SEQ_MISMATCH** — 2H6D 链 A 与 UniProt 覆盖仅 46%（阈值 50%）
- **[WARN] C5_COMPOSITION** — 两态链组成不同 unphos=['A'] vs phos=['A', 'B', 'C', 'D', 'E', 'F']（复合物/构建体差异可能混入磷酸化效应）

### CDK1（WARN）
- **[WARN] C5_COMPOSITION** — 两态链组成不同 unphos=['A', 'B', 'C'] vs phos=['A', 'B', 'H', 'I', 'J']（复合物/构建体差异可能混入磷酸化效应）

### TBK1（FAIL）
- **[FAIL] C3_SEQ_MISMATCH** — 4EUU 链 A 与 UniProt 覆盖仅 42%（阈值 50%）
- **[WARN] C5_COMPOSITION** — 两态链组成不同 unphos=['A'] vs phos=['A', 'B']（复合物/构建体差异可能混入磷酸化效应）

### IKBKB（WARN）
- **[WARN] C5_COMPOSITION** — 两态链组成不同 unphos=['A', 'B', 'C', 'D', 'E', 'F'] vs phos=['A', 'B']（复合物/构建体差异可能混入磷酸化效应）
