[Round 1/4 of autonomous review loop]

Review the work directly from its artifacts — executor notes are not evidence, so
read the files yourself rather than trusting my framing. The full verbatim contents
of the artifacts listed below are appended after this prompt.

- Claims / paper draft: paper/main.tex, paper/refs.bib
- Methods / analysis code under review: scripts/analysis/*.py
- Frozen cohort definition and inclusion criteria: results/manifest_frozen.json
- Raw analysis outputs (verbatim machine JSON, not a summary): results/analysis/*.json
- Result tables as shipped: results/tables/Table1_benchmark.{csv,tex}, Table2_metrics.{csv,tex}, TableS1_stats.{csv,tex}
- Reproducibility entry point: README.md

Please act as a senior computational-biology / ML reviewer (Bioinformatics /
NeurIPS level). Start from the assumption that the work is broken somewhere —
your job is to find where. Be adversarial. Trust nothing the author tells you —
verify everything yourself.

1. Score this work 1-10 for a top venue (Bioinformatics is the target).
2. List remaining critical weaknesses (ranked by severity).
3. For each weakness, specify the MINIMUM fix (experiment, analysis, or reframing).
4. State clearly: is this READY for submission? Yes / No / Almost.

Be brutally honest. If, after genuinely trying to break it, the work holds up and
is ready, say so clearly.

=== SCOPE LIMITS (these bound what you PROPOSE, never what you look for) ===
Report anything that is actually wrong here — including a rare-looking case, if
this work actually produces it. Then keep the fix in scope:
1. This is a small-sample structural-bioinformatics benchmarking paper, not a
   security paper and not a large-scale ML methods paper. Statistical realism
   matters more than defensive over-engineering; and a "fix" that would require
   a new wet-lab dataset or new prediction runs the authors cannot perform is out
   of scope. Say so if the honest fix is "narrow the claim" instead.
2. Check the paper's numbers against the raw JSON/CSV it ships. A claim that the
   manuscript asserts but the shipped result files do not support is the single
   most damaging finding you can report — look for it first.
3. NO speculative machinery: do not ask for extra feature flags, frameworks,
   compat layers or pipelines unless evidence shows a defect they fix.
4. NO corner-case obsession: exotic encodings and millisecond races are out of
   scope for this work.
5. Where a rubric or checklist is genuinely needed, do not over-mechanize
   judgement. A clear sentence a human reads beats a scored table nobody
   maintains.
Say plainly when something is correct. Do not manufacture findings.
