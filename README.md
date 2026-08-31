# Vocal accommodation in marmoset phee-call acoustic and sequence structure

This repository contains the analysis code and processed data used for the manuscript:

**Vocal accommodation in marmoset phee-call sequence structure post pairing**

Wewhare et al. (manuscript)

The project tests whether common marmoset pairs become more similar in the acoustic structure and sequence structure of their phee calls after pair formation, and whether this differs between partner and non-partner social contexts.

Everything needed to reproduce the reported acoustic-space figures, sequence-space figures, Bayesian result tables, and final model figure is included. The notebooks use the saved distance matrices and fixed figure coordinates by default, and the R script uses the saved model objects and posterior draws. A normal run therefore reproduces the results without recalculating DTW distances, retraining the VAE, refitting embeddings, or sampling Bayesian models.

---

## Data

### Acoustic calls

`data/acoustic/` contains:

- `acoustic_calls_processed.csv`: 3,904 annotated call rows;
- `call_order_3612.csv`: the ordered table of 3,612 unique calls used by every acoustic matrix;
- `distances/`: the four accepted `3612 × 3612` distance matrices for STP, MFCC, DTW, and VAE representations;
- `features/`: PCA scores, VAE latent means, and the settings used to obtain the saved acoustic representations;
- `figure_inputs/`: the fixed two-dimensional coordinates, figure settings, subset counts, and selected example calls; and
- `example_audio/`: the 90 unique WAV files used in the final acoustic-space figures.

The order of `call_order_3612.csv` is the row and column order of every acoustic distance matrix.

### Call sequences

`data/sequence/` contains:

- `sequences_processed.csv`: 2,042 annotated sequences;
- `multi_call_sequences_1619.csv`: the 1,619 retained multi-call sequences;
- `session_order_107.csv`: the ordered table of 107 eligible session repertoires;
- `features/`: transition, bigram, phee-repeat, and local-alignment representations;
- `distances/`: the four accepted `107 × 107` sequence-distance matrices; and
- `figure_inputs/`: the fixed coordinates, example sequences, figure settings, and subset counts.

The order of `session_order_107.csv` is the row and column order of every sequence distance matrix.

### Bayesian models

`data/models/` contains:

- `acoustic_session_distances.csv`: 1,724 model rows, corresponding to 431 session comparisons across four acoustic metrics;
- `sequence_session_distances.csv`: 1,324 model rows, corresponding to 331 session comparisons across four sequence metrics;
- `cached_models/`: the two fitted `brms` model objects used for the reported results;
- `posterior_draws/`: the saved posterior stage effects and context contrasts; and
- `input_checks/`: summaries of filtering, scaling, and matrix checks.

The included model objects and posterior draws are the completed fits used for the manuscript. No additional data need to be supplied for the default analysis.

### Source audio

The ordinary repository contains the exact example WAVs required by the
acoustic figures. Larger source archives are provided with the project release
rather than stored in Git history:

- `extracted_calls_3612.tar`: all 3,612 extracted calls (approximately 1.1 GB);
- `source_sequence_recordings_part1.tar` and
  `source_sequence_recordings_part2.tar`: the 2,814 source sequence recordings
  (approximately 2.2 GB in total); and
- `vae_reconstruction_4to12khz.tar`: the accepted VAE model, training record,
  and spectrogram cache for optional reconstruction.

These archives are needed only for optional reconstruction from the recordings. They are not needed to reproduce any reported figure, table, or model result using the saved calculations in this repository.

---

## Analysis

### 1. Prepare and check the analysis data

`code/script_1_data_preparation.ipynb`

**Main tasks:**

- Load the processed acoustic and sequence tables.
- Check the expected row counts and unique identifiers.
- Summarise calls and sequences by individual, stage, and social context.
- Explain how the ordered call and session tables connect to later matrices.

**Default setting:**

```python
REBUILD_ANALYSIS_TABLES = False
```

The default run reads and checks the included files. Setting the switch to `True` rebuilds the ordered call and multi-call sequence tables and writes clearly named `_rebuilt.csv` files to `results/data_preparation/`; it does not overwrite the accepted inputs.

---

### 2. Acoustic spaces

`code/script_2_acoustic_spaces.ipynb`

**Acoustic representations:**

1. **Spectro-temporal parameters (STP):** Euclidean distances between the first five principal components of traditional spectral and temporal measurements.
2. **MFCC:** Euclidean distances between the first five MFCC principal components.
3. **Dynamic time warping (DTW):** exact constrained DTW distance between fixed-reference 4–12 kHz spectrograms.
4. **Variational autoencoder (VAE):** cosine distance between 32-dimensional Chatter latent means.

**Main tasks:**

- Load the four accepted `3612 × 3612` distance matrices.
- Check their call order, shape, symmetry, finite values, and zero diagonals.
- Load the fixed two-dimensional coordinates and exact example calls.
- Recreate the four final 50% acoustic-distribution figures.

Panel A of each figure shows all 3,612 calls and the selected spectrogram examples. Panel B separates the calls by stage and social context. The translucent coloured outlines enclose the central 50% of each animal's observed distribution and are descriptive ellipses, not confidence intervals.

**Default settings:**

```python
RECOMPUTE_ACOUSTIC_DISTANCES = False
RECOMPUTE_EMBEDDINGS = False
MAKE_FIGURES = True
```

**Outputs:**

- `results/acoustic_spaces/acoustic_space_stp_distribution50_reproduced.png`
- `results/acoustic_spaces/acoustic_space_mfcc_distribution50_reproduced.png`
- `results/acoustic_spaces/acoustic_space_dtw_distribution50_reproduced.png`
- `results/acoustic_spaces/acoustic_space_vae_distribution50_reproduced.png`

The matching files without `_reproduced` are the included reference exports.

---

### 3. Sequence spaces

`code/script_3_sequence_spaces.ipynb`

**Sequence representations:**

1. **Transition probabilities:** Euclidean distance between first-order transition-probability vectors.
2. **Bigram distributions:** Euclidean distance between relative bigram-frequency vectors.
3. **Phee-repeat distributions:** Euclidean distance between probabilities of phee repeat lengths two to five.
4. **Local alignment:** distance based on mean normalized Smith–Waterman similarity across all pairs of sequences in two repertoires.

**Main tasks:**

- Load the 107 eligible session repertoires in their exact matrix order.
- Load and validate the four `107 × 107` distance matrices.
- Summarise agreement among the sequence representations.
- Recreate the four final 50% sequence-distribution figures.

The translucent outlines in the distribution figures enclose the central 50% of the observed repertoire distribution for each animal within a stage and social context. They are descriptive ellipses rather than uncertainty intervals.

**Default settings:**

```python
RECOMPUTE_SEQUENCE_DISTANCES = False
RECOMPUTE_EMBEDDINGS = False
MAKE_FIGURES = True
```

**Outputs:**

- `results/sequence_spaces/sequence_space_transition_probability_distribution50_reproduced.png`
- `results/sequence_spaces/sequence_space_bigram_distribution50_reproduced.png`
- `results/sequence_spaces/sequence_space_phee_repeat_distribution50_reproduced.png`
- `results/sequence_spaces/sequence_space_local_alignment_distribution50_reproduced.png`

The accepted reference figures are supplied as PNG files without `_reproduced`.

---

### 4. Bayesian models and final results

`code/script_4_bayesian_models_and_results.R`

The sequence and acoustic analyses are fitted separately using the same Gaussian multilevel model:

```r
z_distance ~ stage * context * metric +
  (1 + stage * context || pair) +
  (1 | mm(session_1_id, session_2_id))
```

**Main tasks:**

- Load and check the two exact model-input tables.
- Load the two completed `brms` model objects and saved posterior draws.
- Calculate the posterior stage effect as `after - before` within every draw.
- Give each experimental pair equal weight and, for combined results, each metric equal weight.
- Recreate metric-specific, combined, and sensitivity result tables.
- Recreate the final interaction figure and model diagnostics.

Negative stage effects indicate convergence; positive stage effects indicate divergence. The stage-by-context interaction is calculated as the Partner stage effect minus the Non-partner stage effect.

**Default setting:**

```r
REFIT_MODELS <- FALSE
```

An ordinary run reads the supplied model objects and posterior cache. The only way to start model sampling is to call the script explicitly with `--refit`.

**Main outputs:**

- `results/bayesian_models/tables/manuscript_primary_results_reproduced.csv`
- `results/bayesian_models/tables/stage_effects_by_context_metric_reproduced.csv`
- `results/bayesian_models/figures/Fig_combined_interaction_panel_reproduced.png`
- `results/bayesian_models/figures/Fig_sequence_stage_by_context_metric_reproduced.png`
- `results/bayesian_models/figures/Fig_acoustic_stage_by_context_metric_reproduced.png`
- `results/bayesian_models/diagnostics/diagnostics_summary_reproduced.csv`
- regenerated posterior predictive checks and trace plots in `results/bayesian_models/diagnostics/`

The included manuscript reference outputs are in the same folders under concise names without a `_reproduced` suffix.

---

## Reproducing the results

Create the Python environment from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r environment/python-requirements.txt
```

Install the R packages once:

```bash
Rscript environment/install_R_packages.R
```

Then open Jupyter and run Scripts 1–3 in order using **Run All**:

```bash
jupyter notebook
```

Finally, recreate the Bayesian tables, figures, and diagnostics from the saved fits:

```bash
Rscript code/script_4_bayesian_models_and_results.R
```

The normal workflow uses only the saved calculations. It does not start any of the costly steps.

---

## Optional recalculation

The expensive operations are kept behind explicit switches:

- Script 2 can recalculate STP and MFCC distances and can refit all four two-dimensional embeddings.
- Exact DTW and VAE reconstruction use the complete-call and model-training materials supplied with the full project release.
- Script 3 can recalculate all four sequence-distance matrices and refit the sequence embeddings.
- Script 4 samples both Bayesian models only when run as:

```bash
Rscript code/script_4_bayesian_models_and_results.R --refit
```

Recalculated files receive `_rebuilt`, `_recomputed`, `refitted`, or `_reproduced` names so that the accepted inputs and reference results are not silently overwritten.

---

## Requirements

The cached workflow was tested with Python 3.13 and R 4.5.2. Exact Python
package versions are listed in `environment/python-requirements.txt`; the R
installer lists the packages used by Script 4. A working C++/Stan toolchain is
needed only for the optional `--refit` route.

---

## Citation

If you use this code or these data, please cite the accompanying manuscript and the archived release of this repository.

## Licence

Code and documentation are released under the MIT License (`LICENSE-CODE`).
Data, audio, figures, and derived outputs are released under CC BY 4.0
(`LICENSE-DATA`).
