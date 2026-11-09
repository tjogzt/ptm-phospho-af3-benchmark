# AF3 待提交作业（v2）

已提交 42 个（`folds_2026_10_23`、`folds_2026_10_24`），以下为**尚未提交**的 32 个。

## 必需：required_paired.json（22 作业）

配对蛋白的完整作业，跑完后位移验证可从 6 组扩到 10 组。

| 蛋白 | 待补作业 | 已提交 |
|---|---|---|
| CDK1 | 6（phos×3 + unphos×3） | 无 |
| PYGM | 6（phos×3 + unphos×3） | 无 |
| c-Src (src) | 5（phos_seed2,3 + unphos×3） | phos_seed1 |
| AMPKα2 (prkaa2) | 5（phos_seed2,3 + unphos×3） | phos_seed1 |

## 可选：optional_single.json（10 作业）

单结构组（PKA-Cα、CTD），**找不到非磷酸化实验对照**，只能看 AF3 内部两态差异，
不参与"AF3 vs 实验"位移验证。配额紧张时可跳过。

## 提交方式

同前：alphafoldserver.com → Upload JSON → 选择文件 → 核对 modifications → Submit。

上传后把 `folds_*.zip` 放入 `results/af3/downloads/`，然后执行：

```bash
python3 scripts/af3/download_af3_results.py --extract
python3 scripts/analysis/compare_af3_vs_pdb.py --all-models
python3 scripts/analysis/plot_af3_vs_pdb.py
```
