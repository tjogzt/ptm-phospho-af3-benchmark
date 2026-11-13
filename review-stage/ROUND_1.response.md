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