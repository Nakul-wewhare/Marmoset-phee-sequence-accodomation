# Sequence-level vocal convergence in common marmosets

This repository contains the essential code, data, and manuscript outputs needed to reproduce the analyses for the current version of the manuscript:

*Sequence-level vocal convergence in common marmosets*  
Wewhare et al., 2026 (manuscript)

Repository URL: [github.com/Nakul-wewhare/Marmoset-phee-sequence-accodomation](https://github.com/Nakul-wewhare/Marmoset-phee-sequence-accodomation)  
Manuscript-ready repository prepared on March 30, 2026

The repository includes the manuscript-essential reproducibility materials:

- extracted single-call WAV files in `extracted calls/`
- processed sequence and spectral tables in `seqeunce data/` and `spectral data/`
- model-ready tables in `glm data/`
- final manuscript summary tables and figures in `glm data/brms_manuscript_outputs-newsession/`

The code here reproduces the main data processing, distance metrics, and Bayesian multilevel models reported in the paper.

The working tree has been trimmed to the files required for the current manuscript workflow. Historical folder names such as `seqeunce data` are intentionally retained because the notebooks refer to them directly.

## Main workflow

The main reproducible workflow follows the manuscript Methods section:

1. `code/script_2_spectral_measurements.ipynb`
2. `code/script_3_seq_measurements.ipynb`
3. `code/script_4_GLM preprocessing.ipynb`
4. `code/script_5_R_Bayesian_models.R`

The repository now starts from the processed sequence table, processed spectral table, extracted call WAV files, and model-ready analysis tables required by the current manuscript. Earlier manual extraction and annotation materials, legacy code, duplicate outputs, and bulky saved model objects were removed because they are not required for the current manuscript version.

## Citation and licensing

- Citation metadata is provided in `CITATION.cff`.
- Code in `code/`, `scripts/`, `environment/`, and repository documentation is released under the MIT License.
- Data, audio, figures, and derived analysis outputs are released under CC BY 4.0.
- The licence split is described in `LICENSE`, with full text or pointers in `LICENSE-CODE` and `LICENSE-DATA`.
- Copy-ready manuscript and submission text is provided in `docs/MANUSCRIPT_TEXT.md`.
- Submission preparation notes are provided in `docs/PROCEEDINGS_B_DATA_CODE_SHARING.md` and `docs/SUBMISSION_CHECKLIST.md`.

## Directory guide

- `extracted calls/`: extracted phee-call WAV files used for acoustic comparisons.
- `seqeunce data/`: processed sequence data tables.
- `spectral data/`: processed spectral tables used by the current manuscript workflow.
- `glm data/`: model-ready analysis tables and current manuscript summary outputs.
- `figures/manuscript_exports/`: manuscript figure exports copied from the working analysis folder.

## Reproducing the results

### Fastest route: rerun only the final Bayesian models

The model-ready tables are already included in `glm data/`.

```bash
Rscript environment/install_R_packages.R
cd code
Rscript script_5_R_Bayesian_models.R
```

This writes the final outputs to `glm data/brms_manuscript_outputs-newsession/`.

Set `MANUSCRIPT_REPO_SAVE_AUXILIARY=1` to also save per-model diagnostic files, and set `MANUSCRIPT_REPO_SAVE_MODELS=1` to save full fitted model objects.

### Full scripted route from the included processed inputs

Install the Python and R dependencies, then run the notebooks and final model script in order:

```bash
python -m pip install -r environment/python-requirements.txt
Rscript environment/install_R_packages.R
bash scripts/run_repro_pipeline.sh
```

### Scope of the trimmed repository

This trimmed repository is designed to reproduce the current manuscript from the inputs actually required by the automated workflow: extracted call WAV files, processed sequence and spectral tables, and model-ready analysis tables. It does not include earlier manual preprocessing artifacts, duplicate historical outputs, or large saved model objects that are not required for the current manuscript version.

## Methods crosswalk

- Sequence metrics (transition probabilities, bigrams, repeat distribution, local alignment): `code/script_3_seq_measurements.ipynb`
- Acoustic metrics (STP, MFCC, DTW): `code/script_2_spectral_measurements.ipynb`
- Model-table assembly: `code/script_4_GLM preprocessing.ipynb`
- Bayesian multilevel models and manuscript-ready outputs: `code/script_5_R_Bayesian_models.R`
