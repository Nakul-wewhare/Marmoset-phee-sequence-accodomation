# Results directory

This directory separates generated analysis output from values copied from the
manuscript for validation.

- `figures/` contains provenance-locked current manuscript snapshots for
  Figures 1–4 and S1–S8. Editable SVG sources are also present for Figure 2 and
  Figure S1. The snapshots are visual references, not regenerated evidence.
- `posterior_draws/` is the canonical location for portable posterior contrast
  draws exported by the R workflow when a complete model cache is installed.
  Every `.csv.gz` draw file requires its adjacent `.csv.gz.metadata.json`
  provenance sidecar.
- `tables/` contains summaries generated from those draws.
- `reference/` contains transcription-only targets from the current manuscript.
  Reference values are never used as model inputs or presented as regenerated
  results.

The notebooks and validators treat a missing generated output as missing. They
never replace it with a reference value, a locked PNG, or a legacy result.
