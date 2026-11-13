# Fig4 图注（v2，2026-10-27，n=11）

## 文件

- `Fig4_why_variable.pdf` — 矢量版（投稿用）
- `Fig4_why_variable.png` — 300 dpi 位图
- 归档：`archive/Fig4_20261025_2226.*`（n=6 旧版）、`archive/Fig4_*_n11.pdf`（当前版）

## 重要变更（v1 → v2）

v1（n=6）曾观察到位点 pLDDT 与幅度比相关 r=+0.649，并写入图注。
扩样至 n=11 后 **r 降至 +0.310（p≈0.35），原假设未获支持**。
本版已撤回原表述，标题改为中性。详见 `audit/FIG4_HYPOTHESIS_REFUTED.md`。

## 英文图注

**Figure 4. Site-region confidence alone does not explain AlphaFold3's variable
accuracy on PTM-induced conformational shifts.**

**(a)** AF3 confidence at the phosphosite region (mean pLDDT over ±10 residues of
the modified site, unphosphorylated state) versus the amplitude ratio of predicted
to experimentally observed displacement, for eleven composition-matched
phospho/unphospho pairs. Dashed line: linear least-squares fit
(Pearson r = +0.31, n = 11, p ≈ 0.35). Red open circles mark the two proteins that
most clearly contradict the expected trend: GSK3β (site pLDDT 94.0, amplitude
ratio 0.25) and TBK1 (45.0, 0.27). Dotted horizontal line marks 1:1.

**(b)** AF3 pLDDT averaged over the phosphosite region (red) versus the whole
protein (grey). Phosphosite regions are lower than the global average in 9 of 11
proteins (mean difference −21 pLDDT), with GSK3β (+7.8) and PAK1 (+13.7) as
notable exceptions.

## 结论表述（写作时使用）

✅ 可用：
- "Site-region confidence explains little of the variance in AF3 accuracy
   across proteins (r = 0.31, n = 11)."
- "Phosphosite regions are generally modeled with lower confidence than the
   protein average, but this does not predict how well AF3 captures the PTM-induced shift."

❌ 不可用：
- "AF3 accuracy is significantly correlated with site-region confidence."
- 任何暗示因果关系的表述

## 仍然成立的三条观察

1. 幅度普遍被低估：11 组中 10 组幅度比 < 1
2. 位点区通常低置信：9/11 组位点 pLDDT 低于全局
3. 方向与幅度脱耦：TBK1 方向 0.559 但幅度 0.27；AKT1 幅度 0.30 但方向 −0.145

## 下一步可检验的替代假设

- 位移的几何类型（loop 重排 / 侧链转动 / 域铰链）
- 磷酸化位点溶剂可及性与静电环境
- 实验结构的分辨率与构象纯度
- AF3 对该蛋白家族的建模水平
