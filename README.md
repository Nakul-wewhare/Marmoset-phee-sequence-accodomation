# Marmoset phee-call sequence accommodation

This repository is a readable record of the analysis used to study changes in
common-marmoset phee calls and call sequences before and after social pairing.
It keeps the analysis in the same notebook-first style in which it was
developed: the calculations, checks, and figures are shown together, with
comments written for someone encountering the project for the first time.

The aim is to make the analysis understandable and reusable. It is not a
one-command reconstruction of every manuscript file, and it is intentionally
not organized as a Python package.

## Analysis files

Read the files in this order:

1. `code/script_1_data_preprocess.ipynb` cleans the extracted call and sequence
   annotations and explains the columns used later.
2. `code/script_2a_spectral_measurements.ipynb` calculates and visualizes the
   acoustic call measurements.
3. `code/script_3a_sequence_measurements.ipynb` calculates and visualizes the
   sequence measurements.
4. `code/script_5a_R_Bayesian_interaction_models.R` fits the final interaction
   models and summarizes the before-to-after posterior contrasts.

The notebooks retain the direct, exploratory structure of the original
analysis, but duplicated cells and private working notes have been removed.
Markdown and code comments explain why each section exists, what it reads, and
what it produces.

## Expensive calculations are off by default

Some pairwise acoustic and sequence comparisons take a long time. Each
notebook has a clearly marked setting near the beginning that leaves those
steps disabled by default and loads the previously calculated tables instead.
Turn a recomputation switch on only when you have the complete source data and
intend to replace the cached result.

The saved call- and repertoire-distance matrices expected by Scripts 2a and 3a
are not included in the current public checkout. Their default branches report
that clearly instead of beginning the costly all-pairs calculations. Script 3a
can still produce its descriptive sequence summaries from the processed CSV.
Script 2a looks for the three `*_dist_matrix_2200calls.npy` files and
`labels_2200calls.csv`; Script 3a looks for
`sequence_distance_matrices.npz` and `sequence_session_labels.csv`. All belong
in a local `distance_matrices/` folder.

Script 4a is not presented as another public workflow step. Its model-table
outputs are treated as precomputed inputs to Script 5a:

- `glm data/glm_data_seq_interaction.csv`
- `glm data/glm_data_call_interaction.csv`

The public repository does not include the call interaction table or saved
`brms` model objects. They must be supplied before Script 5a can complete the
call-structure analysis. The R script checks for these files and stops with a
plain explanation when they are unavailable; it does not silently rebuild
them.

## Setup

Python 3.11 or newer is recommended. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r environment/python-requirements.txt
jupyter notebook
```

Open the notebooks from `code/` and run them in numerical order. Their default
settings favor the included processed or cached inputs.

Install the R packages separately:

```bash
Rscript environment/install_R_packages.R
```

Stan compilation also requires a working C++ toolchain. That toolchain is only
needed when model fitting is explicitly requested.

## Running Script 5a safely

The ordinary command is cache-only:

```bash
Rscript code/script_5a_R_Bayesian_interaction_models.R
```

It validates the two precomputed input tables and loads compatible saved model
objects when they exist. It never starts sampling. If an input or model cache
is missing, the script reports it and stops safely.

Model fitting is always an explicit choice:

```bash
# Fit only models whose saved model is missing
Rscript code/script_5a_R_Bayesian_interaction_models.R --fit

# Replace both saved models, even when caches exist
Rscript code/script_5a_R_Bayesian_interaction_models.R --refit
```

Both fitting modes can take hours. Saved models and generated summaries are
written below `glm data/brms_interaction_outputs/`.

## Data folders

- `extracted calls/` contains the individual WAV clips available with this
  repository.
- `spectral data/` contains processed acoustic measurements.
- `seqeunce data/` contains processed sequence annotations. The historical
  spelling is retained because the notebooks use it.
- `glm data/` contains model-ready tables.

The repository does not contain every intermediate object from the original
working computer. In particular, a reader should not enable an expensive
recalculation unless the inputs named in that notebook are available locally.

## Citation and licences

Citation metadata is in `CITATION.cff`. Code is released under the MIT License;
data, audio, figures, and derived outputs are released under CC BY 4.0. See
`LICENSE`, `LICENSE-CODE`, and `LICENSE-DATA` for the licence split.
