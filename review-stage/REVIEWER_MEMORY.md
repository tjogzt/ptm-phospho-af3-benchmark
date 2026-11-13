# Reviewer Memory — append-only

> 本文件为评审器的跨轮记忆（append-only）。技能 `auto-review-loop`
> 在 `REVIEWER_DIFFICULTY=medium` 下仍要求保存每轮评审原文，供后续轮次与
> 终止文档引用。**禁止删改历史条目**；新轮次只能追加。
>
> This file is the reviewer's round-to-round memory (append-only). Never edit or
> delete prior entries; only append new rounds.

- **Run ID**: `run_20261112_882bced5`
- **Reviewer backend**: `llm-chat`
- **Reviewer model**: `zai-org/GLM-5.3` (family: zhipu, backend-reported)
- **Executor model**: `deepseek-v4-pro` (family: deepseek, caller-declared)

---

## Round 1 — 2026-11-12

- Transport: llm-chat (local HTTP client)
- Endpoint: https://api.siliconflow.cn/v1
- Reviewer model: `zai-org/GLM-5.3` (backend-reported) → family `zhipu`
- Executor model: `deepseek-v4-pro` (caller-declared) → family `deepseek`
- Independence verified: true (different family)
- Files supplied: 33 primary artifacts (verbatim, ≈350 KB user message)
- Usage: prompt 123,524 tok / completion 71,537 tok (reasoning 68,167)
- Elapsed: 776.91 s
- Gate (`review_gate.py --round-backend llm-chat`): `{"decision":"continue","identity_assurance":"unverified","next_backend":"llm-chat","reason":"positive threshold not met","requires_external_acquittal":false}`
- Reviewer verdict: **Score 5/10 — not ready** (major revision)

### Reviewer response (verbatim)

# Review: "Data purity, not model capability, dominates apparent AlphaFold3 performance on phosphorylation benchmarks"

## 1. Score: **5 / 10** (major revision; not acceptable as-is, but salvageable without new experiments)

I recomputed every checkable number in the manuscript from the shipped JSON/CSV. The headline statistics are faithful to the artifacts — which is rare and to the authors' credit — but I found (a) three sentences the shipped result files directly contradict, (b) a Methods section that describes a *different protocol than the code implements*, (c) a stated inclusion rule that was not applied, and (d) a purity gate that fails on confounds visible in the authors' own manifest. None of these requires new predictions or wet-lab work to fix.

## What I verified and found correct (so the authors know the audit was real)

- **Cohort statistics** (af3_vs_pdb_summary.json): r=+0.1132, CI [−0.3433, +0.5744], slope +0.053, Spearman +0.1536/p=0.5848, mean pred 0.521±0.359, exp range [0.25, 2.51], LOO max leverage Src→+0.2202 — all match the text and Table S1 exactly.
- **Confidence hypotheses**: I recomputed site-pLDDT vs amplitude ratio from performance_explained.json: r=+0.478 → +0.215 without IRF3 (paper: +0.48 → +0.21). ΔPAE(site→rest) from pae_metrics.json: r=+0.601 → +0.142 (paper: +0.60 → +0.14); global ΔPAE +0.494 → +0.063; site-minus-global +0.429. All correct.
- **Purity stratification**: pure median ratio 0.460 (n=9), impure 0.303 (n=6); stratum correlations +0.221 / −0.239 — match Fig 1b claims.
- **ESMFold cohort**: r=−0.3622, CI [−0.7175, +0.3079], mean 0.253±0.409, PAK1 leverage → −0.4705 — all match.
- **Noise floor summary values** (0.335 / 1.309 / 7 proteins / median 1.50) and **geometry exploratory correlations** (r=−0.602 p≈0.017 vs r_full; −0.616 p≈0.015 vs mean disp; hinge −0.398 p≈0.14) reproduce from geometry_metrics.json.
- **Table 2**: I spot-checked all 15 rows against both summary JSONs; every cell matches.
- The B2 fix (SEP/TPO/PTR retained via `struct_utils`) and chain-audit (no fallbacks in either log) are real and verified in code.
- The negative result is **robust**: pure-only (n=9, r≈+0.22), and my own recomputation excluding the two P3-failing pairs (n=13) gives r≈+0.21 — still null. This robustness is currently undocumented and would *strengthen* the paper if reported.

## 2. Critical weaknesses, ranked

**W1. The Methods describe a protocol the code does not implement (reproducibility-critical).**
§2.4 states predicted fields are obtained "by superposing each predictor's non-phospho model onto the same reference **and applying that transform to its phospho model**… **without ever superposing predicted onto experimental coordinates**." The shipped code (`compare_af3_vs_pdb.py`, `compare_esmfold_vs_pdb.py`) does the opposite on both counts: it Kabsch-fits **each state independently** onto the experimental unphos reference (`R_u,t_u = kabsch(P_au,Q)` *and* `R_p,t_p = kabsch(P_ap,Q)`), and the JSON ships `rmsd_af3unphos_vs_pdb` for every protein — direct evidence of predicted-to-experimental superposition. The ESM script's own docstring says "用同一变换" while the code computes an independent transform. Since the central quantitative claim (near-constant 0.52 Å predicted displacement) depends on the fitting convention, this must be reconciled.
*Minimum fix:* rewrite §2.4 to describe the implemented (and defensible) protocol — independent per-state fits remove each prediction's arbitrary frame — and add a sensitivity analysis using (i) direct phos→unphos superposition of the two predictions and (ii) the shared-transform variant currently described. Re-analysis of shipped models only.

**W2. Manuscript claims contradicted by the shipped result files.**
- §4.6: the amplitude-binned ratio "falls **monotonically** from 0.54 … to 0.05." residue_spectrum.json gives 0.535, 0.282, 0.200, 0.129, **0.147**, 0.049 — non-monotonic at the 3–5 Å bin.
- §4.6: "residues far from the site (>80 residues away) show the **weakest and least correlated** response." The >80 bin has ratio 0.315 (second *highest*) and corr 0.166 (mid-pack). The weakest ratio is the 0–5 bin (0.121); the weakest correlation is 20–40 (−0.079). The sentence is backwards relative to the shipped data.
- §4.3: "ESMFold's … per-residue magnitude correlations are **uninformative**." Table 2 / esmfold JSON: ESM magnitude r averages ≈**0.50** (PYGM 0.99, AMPK 0.77, Src 0.72) — slightly *higher* than AF3's ≈0.46. And ESM direction cosines average ≈0.19 vs AF3's ≈0.22 — not "centred near zero" in any way that distinguishes it from AF3.
*Minimum fix:* correct all three sentences against the existing JSON. The defensible ESM statement is "shape correlation persists at ~25× amplitude under-prediction," not "uninformative."

**W3. The stated gate rule was not applied.**
§2.2: "A pair that fails any of P1–P4 **is excluded**." The manifest shows AMPKα2 and TBK1 fail P3 (coverage 0.464 / 0.418 < 0.5) yet both have `in_primary_cohort: true` and contribute to every headline statistic (n=15). The Results text even calls them "excluded on P3/P5" while Table 1 shows "In cohort: Y."
*Minimum fix:* correct §2.2 to describe the actual rule and report the n=13 sensitivity (r≈+0.21 by my recomputation — the conclusion survives, so this fix costs nothing but honesty), or actually exclude them and regenerate.

**W4. The purity gate does not control the confounds its own manifest advertises.**
P5 checks chain IDs only. From the frozen manifest's own titles: the "pure" PAK1 unphos structure 1YHW is a **K299R point mutant**; AURKA pairs CD532 (inhibitor) with ATP; GSK3B pairs AMPPNP with an inhibitor; CDK2 pairs two different inhibitors; PYGM pairs two different ligands; Src pairs an inhibitor-trapped *inactive* conformation with the ground state; AKT1 pairs two different inhibitors. CDK7's experimental phospho state is **dual-phosphorylated (pS164+pT170)** while AF3 modeled only pT170. Moreover, AF3 was run on apo full-length sequences, so ligand-induced experimental differences are invisible to the predictor by construction. "Pure" therefore means "composition-matched," not "phosphorylation-isolated" — which undercuts the title's causal claim.
*Minimum fix (no new predictions):* add ligand / mutation / modification-set columns to the audit (all extractable from the shipped PDB metadata), rename "pure" → "composition-matched" or add a P6 criterion and re-stratify Fig 1b, and state explicitly that apo predictions are being compared to ligand-bearing crystal pairs. This reframing actually *supports* the paper's thesis that the standard protocol cannot demonstrate the capability.

**W5. Title/abstract overclaim relative to the evidence.**
"Data purity, not model capability, **dominates**" / "pairing purity … **governs** the apparent result" rests on a 9-vs-6 stratification with no uncertainty: impure r=−0.24 (n=6) has a Fisher CI of roughly [−0.88, +0.71]; pure r=+0.22 (n=9) roughly [−0.52, +0.77]. The sign flip is noise-level, and the direction-cosine medians are essentially identical across strata (≈0.22 vs ≈0.20 by my computation), so the "purity effect" is metric-dependent. The Results section itself pledges to "avoid any claim whose sign depends on the mixture" — the title violates that pledge.
*Minimum fix:* soften to "confounds / can invert," and add CIs or a permutation test for the stratified statistics.

**W6. The noise floor is the wrong kind of floor, and the claim overgeneralizes.**
The experimental floor uses NCS copies *within one crystal* — but the comparison is *cross-crystal* (different space groups, packing, conditions), where same-state variability is systematically larger. Only 7/15 proteins have any floor (8 have `n_replicates: 0`), and PAK1's 0.25 Å "signal" — a pure-cohort leverage point — is *below* the cohort's within-crystal floor median (0.335). "The experimentally observed change therefore exceeds the upper bound of same-state noise" is supported only cohort-median-wise, and for IRF3 only marginally (0.33 vs 0.288).
*Minimum fix:* state the per-protein comparison for the 7 proteins with replicates, flag the 8 without, soften the sentence, and — cheap and high-value — add a cross-crystal same-state floor from additional same-state PDB entries (downloads only, no predictions).

**W7. Two shipped pipelines disagree on the same quantity.**
geometry_metrics.json vs af3_vs_pdb_summary.json for the *experimental* displacement: CDK7 mean 1.001 Å / max 32.7 vs 0.53 / 2.05; AKT1 2.214 / 41.4 vs 2.05 / 34.17. The geometry metrics feeding §4.6's correlations include 30–40 Å outlier residues excluded by the main pipeline — precisely the "amplitude summaries dominated by a few mobile residues" artifact the Introduction warns about.
*Minimum fix:* recompute geometry on the identical common-residue set used by `compare`, reconcile CDK7/AKT1, and re-run the §4.6 correlations on consistent inputs.

**W8. The reproducibility entry point (README) describes a different, older project.**
It states only 3 valid pairs (with wrong PDB IDs: 1HCL/1JST, 1MQ4/1OL5), "only ERK2/p38 pure," 10 paired, "10 pure," and reproduction steps covering 2 of ~9 analysis scripts — none of which produce the paper's tables (make_tables.py is absent), and pointing to a figure (Fig1_af3_ptm_effect) not in the paper. A stale `pair_purity_analysis.json` (n=6/4/10) ships alongside the current artifacts.
*Minimum fix:* regenerate the README from manifest_frozen.json; list the full script chain; delete or date-stamp stale artifacts.

**W9. Table 1 contains internal contradictions.**
PRKACA and POLR2A show verdict **PASS** with N in P3/P4/P5; MAPK3 (ERK1) passes all five criteria yet shows **Pure = N**. The manifest's `pure_pass: 12` also conflicts with the paper's "only 9 pass."
*Minimum fix:* mark single-state rows "n/a," set Pure=Y for ERK1 with In-cohort=N, fix the counts.

**W10. Smaller but real.**
(i) Abstract calls r=+0.11 the correlation of "per-residue Cα displacements" — it is per-protein *means* (n=15 proteins, not 3,355 residues). (ii) Fig 4's caption promises an IRF3-removed comparison; the generating code (`explain_performance.plot`) draws no such panel. (iii) The residue spectrum uses model 0 only — unstated. (iv) Within-protein model instability is unreported despite being shipped: 6/15 proteins have CV>0.3 (CDK2 spans 0.27–1.95 Å across its 15 models; AKT1 CV 0.95) — this qualifies the "near-constant" narrative. (v) Residual site-numbering errors: PKA listed T198 vs 6QJ7 MODRES TPO197 (+SEP338); CTD 1614/1621 vs MODRES 1621/1628 (no effect on the primary cohort, but the audit trail still has off-by-ones). (vi) "Existing benchmarks rarely audit…" is asserted with zero citations; the related-work section (7 refs) does not engage the literature being critiqued. (vii) `nam2025alphafold3` key vs year 2026. (viii) Funnel's first bar (24, "Curated candidate pool") vs text "18 curated entries." (ix) "Effect-type strata… reported only descriptively" — no strata are reported anywhere.

## 3. Minimum-fix summary
Every fix above is text correction, re-analysis of shipped artifacts, or PDB-metadata auditing. The only quasi-new data I would ask for are cross-crystal same-state floors (W6) — public downloads, no modeling. If the honest fix for W4 is "narrow the claim from *pure* to *composition-matched*," say so; that is in scope and, in my view, makes the paper *stronger*.

## 4. Ready for submission? **No.**
The empirical core — a purity-audited 15-pair benchmark on which AF3's predicted displacement is uncorrelated with experiment (r=+0.11, robust across cohort definitions), near-constant, and barely above the predictor's own sampling floor, with both confidence signals reducible to one leverage point — is real, verifiable, and useful. But a Bioinformatics reviewer who checks the shipped JSONs (as I did) will find three contradicted sentences, a Methods section describing a different protocol than the code, an unapplied exclusion rule, and a purity gate blind to confounds its own manifest displays. That is one honest revision cycle — perhaps 2–4 weeks, no new predictions — away from a credible submission (I would then rate it 7–8). As shipped: **5/10, major revision, not ready.**

---

## Round 2 — 2026-11-13

- Transport: llm-chat (local HTTP client)
- Endpoint: https://api.siliconflow.cn/v1
- Reviewer model: `zai-org/GLM-5.3` (backend-reported) -> family `zhipu`
- Executor model: `deepseek-v4-pro` (caller-declared) -> family `deepseek`
- Independence verified: true (different family)
- Files supplied: 36 primary artifacts (verbatim, incl. full Round 1 memory)
- Usage: prompt 142,257 tok / completion 42,123 tok (reasoning 37,393)
- Elapsed: 552.98 s
- Gate (`review_gate.py --round-backend llm-chat --score 6.5 --verdict "not ready"`): `{"decision":"continue","identity_assurance":"unverified","next_backend":"llm-chat","reason":"positive threshold not met","requires_external_acquittal":false}`
- Reviewer verdict: **Score 6.5/10 — not ready** (major-minor revision; "No — but close")

### Reviewer response (verbatim)

# Round 2 Review: "Data purity, not model capability, can dominate apparent AlphaFold3 performance on phosphorylation benchmarks"

## 1. Score: **6.5 / 10** (major–minor revision; one more focused pass, no new predictions needed)

This is a genuinely improved revision. The hardest reproducibility problems from Round 1 (W1, W3, W7, W8, W9) are fixed and verifiably so; the headline statistics still reproduce exactly from the shipped JSON; and the new artifacts (`projection_sensitivity.json`, `pair_confound_audit.json`, `excl_p3_failures`) are exactly the right kind of evidence. But the revision has a pattern: **artifact-level fixes that were never propagated into the manuscript text.** The confound audit exists but the paper still sells "pure" without it; the model-instability numbers are now computed and shipped but unmentioned; and one sentence I explicitly flagged in Round 1 as contradicted by the shipped data (the ESMFold per-protein claim) survives essentially verbatim. I also found one new manuscript-vs-artifact discrepancy this round that I missed in Round 1 (the residue-level analysis silently covers only 10 of 15 proteins).

## What I verified holds (so the authors know the audit was real)

- **All headline numbers reproduce**: r=+0.1132, CI [−0.3433, +0.5744], slope +0.053, Spearman +0.1536/p=0.5848, mean 0.521±0.359, range [0.25, 2.51], Src leverage → +0.2202; ESM r=−0.3622, CI [−0.7175, +0.3079], mean 0.253±0.409, PAK1 leverage → −0.4705. Table 2 spot-checks (all 15 rows, both predictors) match the JSONs cell-for-cell.
- **W7 reconciliation is real and complete.** geometry_metrics.json now matches compare for all 15 proteins (CDK7 0.531/2.052 vs 0.53/2.05; AKT1 2.052/34.17 vs 2.05/34.17; I checked all 15). I recomputed the §4.6 correlations from the *new* JSON: amplitude ratio vs r_full **r=−0.536** (paper −0.54, p=0.040 ✓), vs mean displacement **−0.581** (−0.58, p=0.023 ✓), vs r_block **−0.532** (−0.53, p=0.041 ✓), vs hinge **−0.384** (−0.38, p=0.16 ✓).
- **Purity stratification reproduces**: pure median ratio 0.460 (n=9), impure 0.303 (n=6); I recomputed pure r=**+0.221** and impure r=**−0.239** from the per-protein JSONs — matches the text (+0.22/−0.24).
- **P3-exclusion sensitivity** (n=13, r=+0.2134, CI [−0.2812, +0.7382]) matches my Round 1 independent recomputation (≈+0.21).
- **Projection sensitivity**: indep r=+0.1132 vs direct r=+0.1214, mean 0.521 vs 0.525 — the "near-constant, uncorrelated" conclusion is convention-robust, and Table S1 reports it.
- **Confidence analyses** (pLDDT +0.48→+0.21; ΔPAE site→rest +0.60→+0.14; global +0.49→+0.06) still reproduce from the unchanged JSONs.
- **Chain resolution**: both chain logs show `fallback: false` for every file, including CDK7 8P79 chain J and IKBKB chains B.
- **Noise floor values** (0.335/1.309/7 proteins; predictor 0.389) match; the new §2.5 statement that the predictor floor is restricted to experimentally resolved residues matches `noise_floor.py`.

## 2. Disposition of Round 1 weaknesses

**W1 (Methods described a protocol the code does not implement) — ADDRESSED (90%).** §2.4 now describes independent per-state superposition, which is what the code does; the direct-convention sensitivity is shipped and reported; the ESM docstring contradiction is fixed. *Residual:* the new text says each predicted state is superposed "onto its **experimental counterpart**." The code fits **both** predicted states onto the experimental **unphosphorylated** reference (`Q = r_unph[0][ip_up]`; `kabsch(P_ap, Q)`). A literal implementation of the sentence (AF3-phos→PDB-phos, AF3-unphos→PDB-unphos) would put the two aligned predictions in *different* frames and make Δ_pred meaningless. One clause fixes it: "both predicted states independently superposed onto the common experimental non-phospho reference frame."

**W2 (claims contradicted by shipped files) — PARTIALLY ADDRESSED (2 of 3).**
- (a) Monotonicity: **fixed** — "not strictly monotone… essentially flat (≈0.13–0.15) across the 2–5 Å bins" matches residue_spectrum.json (0.129, 0.147).
- (b) Distance claim: **fixed** — now correctly states the lowest ratio is in the 0–5 window (0.121) with no trend (range 0.121–0.325 ≈ "0.12–0.33").
- (c) ESMFold per-protein claims: **NOT FIXED.** §4.3 still says "ESMFold's per-protein direction cosines are also centred near zero and its per-protein magnitude correlations are uninformative (Table 2, Fig 2b)." From esmfold_vs_pdb_summary.json: ESM magnitude r = 0.99, 0.77, 0.72, 0.68, 0.68, 0.67, 0.63, 0.55, 0.52… **mean ≈ 0.50, slightly higher than AF3's 0.46**; ESM mean cosine ≈ +0.19 vs AF3's +0.22. The sentence contradicts the very table it cites, and it treats AF3's identical metric as signal elsewhere ("p38, r=+0.96; AURKA, r=+0.85"). This is the single most damaging remaining item.

**W3 (unapplied exclusion rule) — ADDRESSED.** §2.2 now states the actual rule (P1/P2 hard; P3–P5 flags), the n=13 sensitivity is shipped and in Table S1, and §4.1 reports it. *Residuals:* §2.2 promises "each cohort-level result is also reported after removing pairs that fail P3 **or P5** (Supplementary Table S1)" — Table S1 contains only the P3 exclusion (the P5 result lives in §4.1/Fig 1b; cross-reference it instead). And §4.1's "either compositionally inconsistent (P5; IKBKB, IRF3, CDK7, CDK1) or failing coverage/overlap (P3; AMPKα2, TBK1)" misclassifies AMPK/TBK1, which fail **both** P3 and P5 per the manifest.

**W4 (purity gate blind to its own confounds) — PARTIALLY ADDRESSED: artifact yes, manuscript no.** `audit_confounds.py` + `pair_confound_audit.json` are real and thorough (ligands, SEQADV mutations, MODRES sets, apo-input statement), and the README adds the caveat. But **main.tex never mentions any of it**: no ligand, mutation, extra-phospho-site, or "composition-matched" discussion anywhere; the title and all text still sell "pure." Worse, the new audit now *contradicts* two manuscript sentences: (i) §2.2 P2's "a structure that merely lacks the phosphate but is otherwise **wild type**" — the audit shows AURKA's unphos state carries T287D (engineered, target chain, one residue from the phosphosite T288) and IKBKB's unphos reference 4E3C carries **S177E/S181E phosphomimetic substitutions at exactly the two sites under study**; (ii) the README's "'pure' 应理解为 composition-matched" collides with the audit's own strict `composition_matched` flag, which is true for only **2 of 16** pairs (ERK2, p38) — 7/9 "pure" pairs differ in bound ligands. Also, the audit's `extra_phos_sites` conflates PDB and UniProt numbering: ERK2 "183" and PYGM "14" are numbering offsets (same site as modelled), while CDK7 pS164, AMPK pS108, IKBKB pS257 are genuine additional experimental sites AF3 did not model — the flag as computed will mislead future users.

**W5 (title/abstract overclaim) — PARTIALLY ADDRESSED.** "Can dominate / can govern" hedging is an improvement. But the pure-vs-impure sign flip (+0.22, n=9 vs −0.24, n=6) is still presented as a finding ("This is not a cosmetic filter"; Discussion: "pooling such pairs changes the sign") with **no uncertainty quantification** — both stratum CIs span zero by a wide margin. No permutation test or Fisher CI was added.

**W6 (noise floor) — PARTIALLY ADDRESSED.** The eight floor-less proteins are now flagged in §4.4 (good), and the predictor-floor residue restriction is documented (good). But the sentence "The experimentally observed change therefore exceeds the upper bound of same-state noise" survives unchanged and is false per-protein: 6 of 15 signals are below the maximum floor (1.309 Å) — PAK1 0.25, IRF3 0.33, CDK7 0.53, GSK3B 0.63, CDK2 0.78, PYGM 0.94 — and PAK1's 0.25 Å is below the cohort floor *median* (0.335). The defensible statement is per-protein: all 7 proteins with replicates exceed their own floor (IRF3 marginally, 0.33 vs 0.288; GSK3B modestly, 0.63 vs 0.495). No cross-crystal floor was added (still the only quasi-new-data ask; optional). The predictor-floor max (3.03 Å, AKT1 — 4× its own signal) is also unmentioned.

**W7 (two pipelines disagree) — FULLY ADDRESSED.** See verification above; all 15 proteins reconciled and the §4.6 correlations re-run on consistent inputs, with paper numbers matching to the third decimal in r and second in p.

**W8 (stale README) — FULLY ADDRESSED.** Regenerated from the frozen manifest; all counts (18/16/2/15/9/10-4-2/6) match; full script→artifact chain; headline table matches the JSONs; historical n=2/6/10/11 documents explicitly superseded.

**W9 (Table 1 contradictions) — FULLY ADDRESSED.** Single-state rows now "n/a" across the gate; ERK1 shows Pure=Y / In-cohort=N with "AF3 models missing"; manifest now distinguishes `gate_pass_paired` (10) from `pure_cohort` (9), matching the text.

**W10 (minor items) — mixed.** (i) Fixed ("one value per protein"). (ii) Fixed (caption no longer promises a phantom IRF3-removed panel). (viii) Fixed (funnel 24 + "24 screened → 18 curated" now consistent). **Not fixed:** (iii) the residue spectrum is still model-0-only and unstated, while §2.3 claims "aggregates all models unless otherwise stated"; (iv) within-protein model instability is now *computed and shipped* (`model_disp_cv`: CDK2 1.006 spanning 0.22–2.02 Å; AKT1 0.745; CDK1 0.933) but the paper still says "close to constant" with no qualification; (v) PKA Thr198 vs 6QJ7 MODRES TPO197 and CTD 1614/1621 vs MODRES 1621/1628 off-by-ones persist in Table 1; (vi) "Existing benchmarks rarely audit…" still has zero citations and the related-work section still does not engage the literature being critiqued; (vii) `nam2025alphafold3` key vs year=2026 persists; (ix) "effect-type strata… reported only descriptively" — still no strata reported anywhere.

## 3. Remaining critical weaknesses, ranked

**R1 (severe — contradicted claim).** §4.3's ESMFold per-protein sentence (W2c above). The shipped JSON and Table 2 show ESM's magnitude correlations average ≈0.50 (≥ AF3's 0.46) and its cosines average +0.19 (≈ AF3's +0.22). *Minimum fix:* one honest sentence, e.g. "ESMFold, like AF3, retains per-protein magnitude correlation (mean r≈0.5) while under-predicting amplitude by roughly an order of magnitude at the median (median ratio 0.04); the shared failure is amplitude and cohort-level tracking, not local shape." Do not fix by deleting the metric from the table.

**R2 (severe — silent analysis-scope error, newly found; I missed it in Round 1).** The §4.6 residue-level analysis ("binning **all 3,355 common residues**") covers only **10 of the 15 cohort proteins**. `residue_spectrum.py`'s stale `SHORT2KEY` map omits GSK3B, PAK1, AKT1, TBK1, IKBKB; the included proteins' n_common sum to exactly 3,355 (296+352+263+232+338+259+274+255+803+283), and the excluded five account for another 1,793 residues (cohort total 5,148). The exclusion is arbitrary, unstated, and presented as cohort-level ("all"). *Minimum fix:* extend the map to all 15 short names, re-run on shipped artifacts, and update §4.6/Fig 3 numbers (the amplitude-compression trend will likely survive; the distance-decomposition "no trend" claim must be re-checked). Also state the model-0-only convention (W10iii).

**R3 (major — W4 not propagated into the manuscript).** The paper's central concept, "pure," is contradicted by its own new audit: 7/9 pure pairs differ in bound ligands; AURKA (pure) has a state-asymmetric T287D adjacent to the phosphosite; IKBKB's "unphosphorylated" reference is an S177E/S181E phosphomimetic double mutant at the study sites; CDK7/AMPK/IKBKB experimental phospho states carry extra sites AF3 did not model; AF3 ran apo full-length against ligand-bearing crystals. *Minimum fix (text only):* one paragraph in §2.2/§4.1 disclosing these (all extractable from `pair_confound_audit.json`), reframe "pure" as "chain-composition-matched," add a caveat to Fig 1b, fix P2's "otherwise wild type," and reconcile the README's "composition-matched" wording with the audit's strict flag (2/16). Fix the audit's numbering conflation in `extra_phos_sites` (map PDB→UniProt before flagging).

**R4 (moderate).** §4.4's "exceeds the upper bound of same-state noise" (W6 residual). *Minimum fix:* per-protein statement for the 7 with floors, flag PAK1 (0.25 Å, no floor, below cohort floor median), soften to "for every protein with a measurable floor, the phospho displacement exceeds its own same-state floor, marginally for IRF3."

**R5 (moderate).** Pure-vs-impure sign flip without uncertainty (W5 residual). *Minimum fix:* Fisher CIs on both stratum r's or a permutation p-value for the difference; one sentence noting both stratum CIs span zero, keeping the "can" hedging.

**R6 (minor).** §2.4 "onto its experimental counterpart" (W1 residual) — one-clause correction.

**R7 (minor batch).** (a) §4.1 either/or misclassification of AMPK/TBK1 and "coverage/overlap" → "coverage"; (b) Table S1 "AF3 $r$ under direct superposition **(Å)**" — a dimensionless r with an Ångström unit (csv and tex); (c) §2.2's Table-S1 promise for the P5 exclusion — cross-reference §4.1 instead; (d) PKA/CTD site off-by-ones in Table 1; (e) cite at least one actual phospho-pairing benchmark for the "rarely audit" claim; (f) bib key/year; (g) either report effect-type strata descriptively or delete the sentence promising them; (h) one sentence qualifying "near-constant" with the shipped within-protein CVs (6/15 proteins CV>0.3; CDK2 spans 0.22–2.02 Å across models); (i) mention the predictor-floor max (3.03 Å, AKT1).

## 4. Ready for submission? **No — but close.**

Every remaining fix is a text edit, a one-line map extension plus re-run of `residue_spectrum.py`, or a caveat paragraph sourced from artifacts already shipped. No new predictions, no wet-lab work. The empirical core — a purity-audited 15-pair benchmark on which AF3's predicted displacement is uncorrelated with experiment (r=+0.11, robust to convention, P3-exclusion, and leverage removal), near-constant, at the predictor's own sampling floor, with both confidence signals reducible to one leverage point — is real, now better triangulated than in Round 1, and useful. If R1–R3 are fixed honestly (corrected rather than deleted), I would rate the next revision 7.5–8 and support acceptance at Bioinformatics.

## 5. Memory update (for Round 3)

1. **Track whether R1 (ESM "uninformative") is corrected or merely deleted** — watch for the metric quietly vanishing from Table 2/Fig 2b rather than being honestly restated. The correct statement is "shape correlation persists at ~10–25× amplitude under-prediction."
2. **Track the 3,355 number**: after the SHORT2KEY fix it should become ~5,148; verify §4.6 and Fig 3 are regenerated and the "all" claim updated; check whether the distance-decomposition "no trend" conclusion survives inclusion of GSK3B/PAK1/AKT1/TBK1/IKBKB.
3. **Track whether "pure" is reframed in main.tex itself** (not just README), and whether AURKA T287D and IKBKB S177E/S181E are disclosed. If a "composition-matched" claim appears, verify it against `pair_confound_audit.json` (strict flag = 2/16, not 9).
4. **Track whether stratified uncertainty (R5) and the per-protein noise-floor statement (R4) appear**, and whether any *numbers shift* when artifacts are re-run — a shifting headline r would be a new problem.
5. **Pattern to watch (executor behavior):** artifact-level fixes are thorough and faithful, but manuscript propagation lags — W4, W10(iii), W10(iv) were all "fixed in the JSON, absent from the paper." In Round 3, diff the *paper text* against the *newest* JSONs, not against the changelog.
6. **Check `extra_phos_sites` numbering fix** (PDB↔UniProt mapping) — ERK2 "183" and PYGM "14" should disappear as false extras; CDK7 "164", AMPK "108", IKBKB "257" should remain and ideally be disclosed in the paper.
7. Residual minor ledger to clear: Table S1 "(Å)" typo, AMPK/TBK1 P3/P5 classification, PKA/CTD off-by-ones, uncited "rarely audit" claim, `nam2025alphafold3` year, effect-strata promise, model-CV sentence, predictor-floor max.
8. My own Round 1 audit missed the 10/15 residue-spectrum coverage — in Round 3, explicitly check *which proteins enter each analysis*, not just whether quoted numbers match.

---

## Round 3 — 2026-11-13

- Transport: llm-chat (local HTTP client)
- Endpoint: https://api.siliconflow.cn/v1
- Reviewer model: `zai-org/GLM-5.3` (backend-reported) -> family `zhipu`
- Executor model: `deepseek-v4-pro` (caller-declared) -> family `deepseek`
- Independence verified: true (different family)
- Files supplied: 35 primary artifacts (verbatim, incl. full Round 1+2 memory)
- Usage: prompt 149,903 tok / completion 41,478 tok (reasoning 37,177)
- Elapsed: 509.5 s
- Gate (`review_gate.py --round-backend llm-chat --score 7.5 --verdict "almost"`): `{"decision":"stop","identity_assurance":"caller_declared","next_backend":null,"reason":"HTTP reviewer returned a positive cross-family verdict; executor identity remains caller-declared","requires_external_acquittal":false}`
- Reviewer verdict: **Score 7.5/10 — almost** (Almost ready; one text-only pass remains)

### Reviewer response (verbatim)

# Round 3 Review: "Data purity, not model capability, can dominate apparent AlphaFold3 performance on phosphorylation benchmarks"

## 1. Score: **7.5 / 10** — Almost ready; one text-only pass remains

This is the revision I asked for. R1 and R2 — the two items I called most damaging — are fixed *honestly* (corrected, not deleted) and verifiably. R4–R7 are cleared. Every number I checked in the manuscript this round reproduces from the shipped JSONs, including the regenerated residue-spectrum analysis. But the revision introduced **one new manuscript-vs-artifact contradiction** (IKBKB misattributed to the pure stratum, in three places) and **one backwards logical justification** in the very paragraph that was added to fix R3. Both are sentence-level fixes; neither touches a statistic.

## 2. What I verified holds (so the authors know the audit was real)

- **Headline unchanged and reproducible**: r=+0.1132, CI [−0.3433, +0.5744], slope +0.053, Spearman +0.1536/p=0.5848, mean 0.521±0.359, range [0.25, 2.51], Src leverage → +0.2202; ESM r=−0.3622, CI [−0.7175, +0.3079], PAK1 leverage → −0.4705. The headline r did **not** shift when artifacts were re-run — my Round 2 demand was met.
- **R2 fix verified end-to-end**: residue_spectrum.json now has n_residues=5148, n_proteins=15; the per-protein n_common values in af3_vs_pdb_summary.json sum to exactly 5,148 (262+255+263+283+296+259+352+344+803+617+232+338+291+279+274). The amplitude bins are now **strictly monotone** (0.498, 0.308, 0.273, 0.224, 0.168, 0.067), so the paper's new "monotonically across the six bins (0.50, 0.31, 0.27, 0.22, 0.17, 0.07)" is true. The distance claim re-verified: ratios 0.147–0.312, non-monotone, lowest in 0–5 (where median exp displacement, 4.34 Å, is indeed largest). The script now asserts full SHORT2KEY coverage — a good defensive fix.
- **R1 fix verified by recomputation**: ESM mean magnitude r = 7.562/15 = **+0.504** (paper +0.50) vs AF3 6.901/15 = **+0.460**; ESM mean cosine (14 non-null) = **+0.187** vs AF3 **+0.219**; ESM median amplitude ratio = **0.043** (paper 0.04) vs AF3 **0.39**. All match §4.3 exactly, and the metrics remain in Table 2 and Fig 2b.
- **R5 verified**: purity_strata pure r=+0.2213 CI [−0.3895, +0.9109]; impure −0.2385 CI [−1.0, +0.8369]; permutation p=0.4007 — in abstract, §4.1, and Table S1. The permutation implementation (20k label reshuffles, seed 7, two-sided |Δr|) is sound.
- **R4 verified against noise_floor.json**: all seven floored proteins exceed their own floor (IRF3 0.33 vs 0.288; GSK3B 0.63 vs 0.495; IKBKB floor 1.309 > CDK7/GSK3B/CDK2/PYGM signals); PAK1 (no replicates, 0.25 Å) correctly flagged; predictor floor 3.033 Å for AKT1 ≈ 3.8× its 0.805 Å signal ("about four times" ✓).
- **R7h verified**: model_disp_cv median = 0.332 ("0.33" ✓); 8 of 15 exceed 0.3 ✓; CDK2 spans 0.22–2.02 Å ✓.
- **Effect-type composition verified**: 11 local / 3 allosteric / 1 disorder-order in the manifest — matches the new Limitations sentence.
- **Audit counts verified**: 7/9 pure pairs with ligand differences; strict composition_matched = 2/16 (ERK2, p38α); extras CDK7 pS164, IKBKB pS257, AMPK pS108 on β-subunit chains B/D — all as §2.2 states.
- Table 1, Table 2 (all 15 rows, both predictors), Table S1, funnel counts (24/16/15/9), and §4.6 geometry correlations (−0.54/0.040, −0.58/0.023, −0.53/0.041, −0.38/0.16) all reproduce.

## 3. Disposition of Round 2 weaknesses

**R1 (ESM "uninformative") — ADDRESSED, honestly.** §4.3 now states ESM "does retain local shape agreement — its mean per-protein magnitude correlation is +0.50, essentially the same as AF3's +0.46... while under-predicting amplitude about an order of magnitude more severely (median ratio 0.04, against 0.39)". Every number verified above. Not deleted: Table 2 and Fig 2b retain the metrics. *Residual:* the section title "A second, sequence-only predictor is **equally uninformative**" is now in tension with the corrected body.

**R2 (silent 10/15 scope error) — FULLY ADDRESSED.** See §2. The "all" claim is now true (5,148/15), the model-0 convention is stated in §4.6, the Fig 3 caption, and the JSON, and the distance-decomposition conclusion was genuinely re-checked rather than assumed.

**R3 ("pure" not propagated into main.tex) — LARGELY ADDRESSED, with a new error.** The "What 'pure' does and does not mean" paragraph is in §2.2; P2 no longer claims "otherwise wild type"; Fig 1's caption carries the caveat; the README reconciles the strict 2/16 flag. **But** the revision misattributes IKBKB to the pure stratum (see W1 below) and omits PAK1's K299R, which the audit identifies as the actual second pure-pair mutation.

**R4 (noise-floor overclaim) — ADDRESSED.** The false cohort-wide sentence is gone; the per-protein statement is present and verified.

**R5 (stratified uncertainty) — ADDRESSED.** CIs and permutation p in abstract, Results, and Table S1; §4.1 explicitly says "evidence that the pooled mixture is unstable, not a demonstrated contrast".

**R6 (superposition clause) — ADDRESSED.** §2.4 now says "superposing each predicted state onto the common experimental non-phospho reference frame... independently for the two states" — exactly what the code does (`kabsch(P_au, Q)`, `kabsch(P_ap, Q)` with `Q = r_unph[0][ip_up]`).

**R7 (minor batch) — 8 of 9 cleared.** (a) AMPK/TBK1 now "additionally fail the P3 coverage criterion" ✓; (b) Table S1 "(Å)" removed from the dimensionless r ✓; (c) P5 exclusion now cross-references §4.1/Fig 1b ✓; (d) off-by-ones resolved by a Table 1 caption footnote (TPO197 = UniProt Thr198; SEP1621/1628 = Ser1614/Ser1621) — internally consistent with the manifest's MODRES records and plausible (PKA mature-chain numbering), though no SIFTS mapping ships to check it externally; (f) `nam2026alphafold3` key/year consistent ✓; (g) effect strata now reported descriptively (11/3/1, verified) ✓; (h) model-CV sentence ✓; (i) predictor-floor max ✓. **(e) partial:** Background now engages two actual AF3-phospho evaluations (nam2026, ramasamy2026 — a real improvement), but the Introduction's "such pairs are routinely treated as 'phosphorylation pairs'" / "Existing benchmarks... rarely audit" still carries no citation.

## 4. Round 2 memory ledger (§5, items 1–8)

1. **ESM sentence corrected, not deleted — CLEARED** (verified §2/§3 above; residual: section title).
2. **3,355 → 5,148 — CLEARED**; "all" claim now true; distance "no trend" re-verified and survives; monotonicity claim now *actually true* (it was false on the 10-protein subset).
3. **"pure" reframed in main.tex — CLEARED WITH A NEW ERROR**: IKBKB is presented as pure-pair heterogeneity in Fig 1's caption, the README, and §2.2's frame (see W1); AURKA T287D and IKBKB S177E/S181E *are* disclosed, PAK1 K299R is not.
4. **Stratified uncertainty + per-protein floor — CLEARED**; headline r did not shift.
5. **Propagation pattern — largely broken**, but the failure mode has *changed*: from omission (fixed-in-JSON, absent-from-paper) to **mis-transcription** (text present but factually wrong against the audit). See W1 and memory item 5 below.
6. **extra_phos_sites numbering — CLEARED**: ERK2 PDB 183/185 → UniProt 185/187 (extra = []); PYGM 14 → 15 (extra = []); CDK7 [164] and IKBKB [257] retained as genuine; AMPK pS108 correctly moved to `extra_phos_sites_other_chains` {B, D}. All three disclosed in §2.2.
7. **Minor ledger — cleared except (e)**, which is partial (two citations added in Background; the Intro claim itself still uncited).
8. **Analysis-scope check — done**: residue spectrum 15/15 (with an assert); all other analyses cover the 15; pae_metrics additionally contains the two single-state entries but they are correctly excluded from correlations (no amplitude_ratio).

## 5. Remaining weaknesses, ranked (N = newly introduced)

**W1 (N, moderate — a checkable factual error in the paper's core conceptual section).** Fig 1's caption states: "the nine pure pairs are not identical apart from the phosphate, and two carry engineered substitutions adjacent to (AURKA) or at (IKBKB) the site under study." The README repeats it ("9 对 pure 中 … 2 对含点突变（AURKA…；IKBKB…）"), and §2.2's paragraph opens "Auditing the nine pure pairs... shows how much heterogeneity remains:" and then lists IKBKB as the final example. **IKBKB is not one of the nine pure pairs** — the manifest shows P5_composition: false, verdict WARN (chains A–F vs A,B), pure_pair: false. Per pair_confound_audit.json, the two pure pairs with mutations are **AURKA (T287D) and PAK1 (K299R)** — `n_pure_with_point_mutation: 2` counts exactly those. PAK1's K299R (present in *both* states) is nowhere disclosed in the paper, even though it is a pure-pair confound and the paper's own audit ships it. A reader checking the audit against the caption will find the contradiction immediately — this is the same class of error that has gated every round.
*Minimum fix:* rewrite the caption, the README clause, and the §2.2 frame: "two of the nine carry engineered substitutions (AURKA T287D, state-asymmetric and one residue from the site; PAK1 K299R, present in both states); among impure pairs, the IKBKB unphosphorylated reference (4E3C) is an S177E/S181E phosphomimetic double mutant at the two sites under study."

**W2 (N, moderate-minor — the new §2.2 justification is backwards, and the Discussion contradicts it).** §2.2: "if residual confounding has any effect it should inflate the apparent response, so the finding is conservative." Two problems. (i) The direction is not uniform in the paper's *own* audit: the IKBKB S177E/S181E and AURKA T287D substitutions sit in the **unphosphorylated** references, i.e., partially phosphomimetic baselines that **deflate** the apparent response, not inflate it. (ii) Where confounding *does* inflate the experimental displacement, it attenuates the pred-vs-exp correlation toward the null — i.e., it could **manufacture** the negative result; that is the opposite of "conservative." The defensible argument is available and stronger: confounding acts in heterogeneous directions, so the pure-stratum sensitivity (r=+0.22, CI spanning zero) is the relevant robustness check — and the negative conclusion survives it. Separately, the Discussion says "pairing impurity **inflates apparent agreement**" while §2.2 says "inflate the apparent **response**"; as written, the Discussion's claim is contradicted by the paper's own stratification (impure r=−0.24 < pure r=+0.22), and the paragraph's "two artefacts work in the same direction" does not cohere under either reading.
*Minimum fix:* rewrite the two sentences as above; in the Discussion, "inflates the apparent response and can destabilise the pooled statistic (here, sign-flipping it)".

**W3 (minor).** Conclusion: "pairing purity and same-state noise calibration, not model capability, **determine** the apparent result." Unhedged verb, in a paper whose own permutation test on the purity contrast is p=0.40. *Fix:* "can dominate".

**W4 (minor).** §4.3 title "equally uninformative" vs the corrected body ("does retain local shape agreement"). *Fix:* retitle, e.g., "fails in the same way".

**W5 (minor).** The "routinely treated as phosphorylation pairs" / "rarely audit" claims remain uncited (Background engagement is new and welcome, but the Intro sentence itself still asserts a literature practice with no reference). *Fix:* attach the nam2026/ramasamy citations or soften to a first-person statement ("we could not find a published phospho-pairing benchmark that audits composition").

**W6 (minor, repo hygiene).** `model_stability.json` is an orphan (no shipped script generates it; absent from the README's script→artifact table) and **disagrees** with af3_vs_pdb_summary.json's model_disp_* fields on the same quantity (AKT1 CV 0.953 vs 0.745; CDK2 0.27–1.95 vs 0.22–2.02). The paper cites the summary values, so there is no manuscript contradiction, but this is a W7-class duplicate-artifact hazard. Relatedly, the audit's `has_mutation` for CDK1 counts SEQADV records on the partner chain (cyclin B1, P14635), not the target — harmless here (CDK1 is impure and not cited for mutations) but worth a one-line filter for future users. *Fix:* delete or regenerate model_stability.json; filter SEQADV by target chain.

**W7 (very minor).** (i) P2's caveat cross-reference points to §4.1 (`sec:res-purity`) but the caveats paragraph lives in §3.2 (`sec:purity`). (ii) §2.3's "aggregates all models unless otherwise stated" has two more unstated model-0 exceptions: the pLDDT analysis (seed1_model_0) and the PAE analysis (model 0 across seeds).

**Not done, acknowledged optional:** the cross-crystal same-state floor (my Round 1/2 ask). The per-protein statement now carries the honest version; I do not press this.

## 6. Ready for submission? **Almost.**

Every remaining fix is a text edit or repo hygiene — no new predictions, no re-analysis. The empirical core is now better triangulated than at any prior round: the null result is robust to projection convention, P3 exclusion, leverage removal, and purity stratification (with uncertainty now quantified); the residue-level analysis finally covers what it claims to cover; and the confound audit is disclosed in the manuscript rather than buried in a JSON. If W1 and W2 are fixed honestly — and W1 in particular must be, because it is a factual error against the paper's own shipped audit — I would rate the next revision **8/10** and support acceptance at Bioinformatics.

## 7. Memory update (for Round 4)

1. **Verify the IKBKB misattribution is corrected in all three locations** (Fig 1 caption, README, §2.2) and that PAK1 K299R is disclosed. Check against `pair_confound_audit.json`: `n_pure_with_point_mutation: 2` = AURKA + PAK1; IKBKB `pure_pair: false`.
2. **Verify the "conservative" rewrite.** Watch for a cosmetic fix (word-swap "agreement"→"response") that leaves "so the finding is conservative" intact — that inference is backwards regardless of wording. The correct sentence must acknowledge heterogeneous confound directions and lean on the pure-stratum sensitivity.
3. **Watch the strong verbs** ("determine" in the Conclusion, "can dominate/govern" elsewhere) for consistency with p=0.40 and the §4.1 "not a demonstrated contrast" statement.
4. **Check model_stability.json is removed or reconciled** with af3_vs_pdb_summary.json's model_disp_* fields.
5. **Pattern update:** the executor's failure mode has shifted from *omission* (fixed in JSON, absent from paper) to *mis-transcription* (text present but factually wrong vs the artifact — IKBKB/PAK1). In Round 4, verify every **named example** (protein, site, mutation, stratum membership) against the audit JSONs, not just the statistics.
6. If the authors assert the Table 1 numbering footnote (TPO197=Thr198; SEP1621/1628=S1614/S1621) is externally verified, ask for the mapping artifact; it is currently an assertion that is internally consistent with the manifest but not independently checkable from shipped files.
7. Residual minor ledger to clear: §4.3 title, "rarely audit" citation, P2 cross-reference, unstated model-0 conventions for pLDDT/PAE, CDK1 partner-chain SEQADV false positive in the audit.

---

## Round 3 — STOP GATE OUTCOME

Stop gate (`review_gate.py --round-backend llm-chat --score 7.5 --verdict "almost"
--executor-model deepseek-v4-pro --reviewer-model zai-org/GLM-5.3`) returned
`decision: stop` — positive threshold met (score 7.5 ≥ 6 AND verdict "almost").
Cross-family independence: executor `deepseek-v4-pro` (deepseek) vs reviewer
`zai-org/GLM-5.3` (zhipu) = different. `requires_external_acquittal: false`, so no
Codex/manual finalizer receipt is required.

The loop therefore terminates at Round 3 of MAX_ROUNDS=4. Remaining reviewer items
W1–W7 are recorded above as a post-stop residual ledger.