# Non-partner sequence robustness analyses

This directory contains the reference outputs for two reviewer-directed sensitivity analyses reported in Figure S9 and Tables S2–S3. Both analyses are restricted to Non-partner sequence structure because this was the outcome and context with evidence of convergence in the primary Bayesian analysis.

## Analyses

1. **Across- versus within-stage variation.** Model 1 compares distances between the same focal animal's session repertoires when the receiver is held constant. Its primary contrast is the Before–After distance minus the mean of the Before–Before and After–After distances. A positive value means that across-stage distance exceeds this average; the component table reports the Before and After comparisons separately.
2. **Bonded versus nonbonded change.** Model 2 compares the stage-related change in distance for three bonded dyads with that for six nonbonded opposite-sex dyads. The primary contrast is the bonded stage change minus the nonbonded stage change; negative values indicate a larger decrease among bonded dyads. All eligible session-level comparisons are retained, and multiple-membership terms account for reused animals and sessions.

Model 2 does not include same-sex dyads. Receiver composition is not balanced between bonded and nonbonded comparisons, so this analysis is interpreted as a sensitivity analysis rather than a perfectly matched control experiment.

## Inputs and code

The analysis reuses the repository's frozen inputs rather than duplicating them:

- `data/sequence/session_order_107.csv`
- `data/sequence/distances/transition_probability_107.npy`
- `data/sequence/distances/bigram_107.npy`
- `data/sequence/distances/phee_repeat_107.npy`
- `data/sequence/distances/local_alignment_107.npy`

The analysis script is `code/script_6_sequence_robustness_models.R`. From the repository root:

```bash
Rscript environment/install_R_packages.R
Rscript code/script_6_sequence_robustness_models.R --validate-only
Rscript code/script_6_sequence_robustness_models.R --refit
```

Create the Python environment described in the root README before running the script. Alternatively, set `RETICULATE_PYTHON` to a Python interpreter with NumPy installed. Validation checks input hashes, constructs both analysis datasets, and builds the Stan data without sampling. A full refit uses four chains, 4,000 iterations per chain (2,000 warmup), seed 123, `adapt_delta = 0.995`, and `max_treedepth = 15`; it may take several hours.

## Included reference outputs

- `figures/Fig_S9_sequence_robustness_posteriors.png` and `.pdf`: combined raster and vector versions of Figure S9.
- `tables/model1_individual_stage_shift_results.csv`: Model 1 primary contrast by metric and combined.
- `tables/model1_component_contrasts.csv`: Model 1 component contrasts.
- `tables/model2_pair_specificity_results.csv`: bonded-minus-nonbonded stage contrast by metric and combined.
- `tables/model2_group_stage_effects_by_metric.csv`: underlying stage effects for bonded and nonbonded dyads.
- `model_inputs/`: frozen-input hashes, support counts, and scaling checks.
- `diagnostics/diagnostics_summary.csv`: convergence and sampler diagnostics for the final fits.
- `analysis_provenance.json`: model specifications, sampling settings, input paths, and input hashes.

The final fits used R 4.5.2, brms 2.23.0, tidybayes 3.0.7, posterior 1.6.1, patchwork 1.3.2, reticulate 1.45.0, rstan 2.32.7, and NumPy 2.5.2. The large fitted-model and posterior-draw files are reproducible from the script and are intentionally omitted from this minimal update.
