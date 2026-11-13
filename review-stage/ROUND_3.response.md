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