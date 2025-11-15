# Marmoset-phee-sequence-accodomation
Here’s a draft README.md you can drop straight into the repo and tweak as needed:

⸻

Vocal accommodation in marmoset phee call sequence structure

This repository contains the analysis code and processed data used for the manuscript:

Vocal accommodation in marmoset phee call sequence structure post pairing
Wewhare et al., 2025 (manuscript)

The project tests whether common marmoset pairs converge in phee call sequence structure (syntax) and/or call structure (spectral features) after pair formation, and whether convergence depends on audience (partner vs non-partner).

The code here reproduces the main data processing, distance metrics, and linear mixed-effects models reported in the paper.

⸻

Repository structure

You can organise (or check) the repo roughly like this:

.
├── README.md
├── manuscript/
│   └── marmoset_syntax_manuscript_V3.pdf
├── notebooks/
│   ├── script_1_data_preprocess.ipynb
│   ├── script_2_spectral_measurements.ipynb
│   ├── script_3_seq_measurements.ipynb
│   └── script_4_GLM_preprocessing.ipynb
├── R/
│   ├── script_5_R_LMM_model.R
│   ├── dfbeta.R
│   ├── diagnostic_fcns2.r
│   └── glmm_stability.r
├── data/
│   ├── raw/              # (not in repo) original wav files, Raven selections, etc.
│   ├── processed/        # cleaned call / sequence tables, MFCCs, distance matrices
│   └── glm_input/        # per-metric distance tables used by the LMMs
└── figures/              # optional: exported plots / figure panels

If your layout differs a bit that’s fine – the scripts don’t hard-require this exact structure as long as file paths are adjusted accordingly.

⸻

Dependencies

Python

The notebooks were run with Python 3.x. You’ll need the usual scientific stack; something like:
	•	numpy
	•	pandas
	•	scipy
	•	scikit-learn
	•	matplotlib
	•	tqdm
	•	jupyter / jupyterlab

Optionally (depending on your edits):
	•	numba
	•	statsmodels

A simple conda env (example):

conda create -n marmoset-syntax python=3.10
conda activate marmoset-syntax
conda install numpy pandas scipy scikit-learn matplotlib jupyterlab
pip install tqdm

R

Analyses were run in R 4.2.2. Required packages (see script_5_R_LMM_model.R for the exact list):
	•	lme4
	•	lmerTest
	•	emmeans
	•	multcomp
	•	DHARMa
	•	tidyverse (or at least dplyr, ggplot2)
	•	data.table (optional but convenient)

Install from within R:

install.packages(c(
  "lme4", "lmerTest", "emmeans", "multcomp",
  "DHARMa", "tidyverse", "data.table"
))

The files dfbeta.R, diagnostic_fcns2.r, and glmm_stability.r are helper functions for model diagnostics (originally from Roger Mundry); these are sourced inside the main R script.

⸻

Data
	•	Raw audio (.wav) files are not included in this repository (file sizes + animal welfare regulations).
	•	The analysis starts from:
	•	call/sequence tables exported from Avisoft-SASLab Pro,
	•	spectral parameter tables exported from Raven Pro, and
	•	MFCCs computed in R using the BehaviouR pipeline (see manuscript Methods).

If you want to rerun everything from scratch you will need:
	1.	A table of phee calls and sequences with:
	•	individual ID,
	•	session ID,
	•	stage (before / after),
	•	conspecific ID,
	•	recipient type (partner / non-partner),
	•	element type (phee, peep, ek/egg, tsk, …),
	•	within-sequence position.
	2.	Per-call spectral parameters (Raven output), and
	3.	Per-call MFCCs (133 coefficients per call in the current implementation).

If you are only interested in reproducing the models and plots and the repository already includes data/processed/ CSVs, you can skip straight to the R and GLM-preprocessing steps below.

⸻

Analysis pipeline

1. Data preprocessing (script_1_data_preprocess.ipynb)

Main tasks:
	•	Read in annotation tables from Avisoft / Raven.
	•	Apply sequence definition rules (≤0.5 s gaps; optional ±1 s inclusion of flanking elements, as described in the Methods).
	•	Clean and filter:
	•	exclude low SNR calls and overlaps,
	•	keep only calls that can be assigned unambiguously to an individual,
	•	drop sessions with <5 units per individual (for a given metric).
	•	Construct per-session / per-individual repertoire tables used by later scripts.

Outputs (typical):
	•	data/processed/phee_calls_clean.csv
	•	data/processed/phee_sequences_clean.csv

You should run this notebook first.

⸻

2. Spectral measurements & distances (script_2_spectral_measurements.ipynb)

Main tasks:
	•	Import spectral parameters (Raven) and MFCCs for all phee elements.
	•	Do separate PCAs for:
	•	traditional spectral parameters,
	•	MFCCs.
	•	Retain first 5 PCs for each representation (matching the manuscript).
	•	Compute pairwise repertoire distances between individuals for each session group and context.

Metrics implemented (see manuscript “Call structure analysis”):
	1.	Mean pairwise Euclidean distance between calls.
	2.	Euclidean distance between repertoire centroids.
	3.	Bhattacharyya / Mahalanobis-style centroid distance (depending on the exact final version of the script).

Outputs:
	•	data/processed/spectral_distances_by_metric.csv
(one row per comparison × metric, with variables like focal_ID, partner_ID, stage, recipient_context, metric, distance).

⸻

3. Sequence metrics & distances (script_3_seq_measurements.ipynb)

Main tasks:

Starting from the cleaned sequence tables, represent each repertoire as an ordered string of call-type symbols, then compute four complementary sequence-structure distances:
	1.	Transition probabilities
	•	Estimate first-order transition matrix T, vectorise, and take Euclidean distance between matrices.
	2.	Bigram distribution
	•	Count all 2-grams, normalise by total bigram count, then compute Euclidean distance between bigram vectors.
	3.	Repeat distribution
	•	Compute distribution of run lengths for phee elements and use Euclidean distance between the resulting probability distributions.
	4.	Local alignment (Smith–Waterman)
	•	Implement a local alignment with match = 1, mismatch = −2, gap = −2.
	•	Normalise the best alignment score by log(sequence length), symmetrise across repertoire pairs, and convert similarity to distance.

Outputs:
	•	data/processed/sequence_distances_by_metric.csv

Again, one row per comparison × metric, to be combined with the spectral distances in the next step.

⸻

4. GLM / LMM preprocessing (script_4_GLM_preprocessing.ipynb)

Main tasks:
	•	Combine spectral and sequence distance tables into a single long-format dataset.
	•	Create variables:
	•	stage (before / after),
	•	Partner_Non_partner (partner vs non-partner),
	•	Conspecific / comparison type (direct vs indirect, or equivalent),
	•	metric (transition, bigram, repeat, alignment, spectral_PC, etc.),
	•	focal_ID, pair_ID, session numbers.
	•	Scale distances:
	•	distance scaled to 0–1 within each metric (for the LMMs reported in the manuscript), or
	•	z-score within metric (earlier report version).
	•	Compute distance change between stages for each combination (after − before).

Outputs (examples):
	•	data/glm_input/sequence_distance_stage_scaled.csv
	•	data/glm_input/spectral_distance_stage_scaled.csv
	•	data/glm_input/distance_change_by_metric.csv

These files are what the R script mainly consumes.

⸻

5. Linear mixed-effects models & plots (R/script_5_R_LMM_model.R)

Main tasks:
	•	Load the prepared GLM input tables.
	•	Fit the LMMs described in the manuscript, separately for:
	•	Sequence structure — partner context
	•	Sequence structure — non-partner context
	•	Call structure — partner context
	•	Call structure — non-partner context

A typical model (final manuscript version) is of the form:

distance_scaled ~ stage +
  (1 + stage | metric) +
  (1 | pair_ID) +
  (1 | session_ID_ind1) +
  (1 | session_ID_ind2)

Key points:
	•	stage (before vs after) is the fixed effect of interest.
	•	Random slopes for stage within metric, plus random intercepts for pair and session, to keep the random-effects structure “maximal” (Barr et al. 2013).
	•	Models are fitted with lme4; p-values from likelihood-ratio tests (full vs reduced model without stage) and lmerTest.
	•	Diagnostic helpers:
	•	dfbeta.R, glmm_stability.r, and diagnostic_fcns2.r provide functions for:
	•	influence diagnostics (DFBETAs),
	•	stability checks for random-effects structures,
	•	residual checks.
	•	Figures:
	•	The script produces the coefficients for Figure 3 (effect of stage on repertoire distance, with 95% CIs) and any other summary plots.
	•	Code comments in the script indicate which chunks correspond to which figure panels.

⸻

How to reproduce the main results

Very short version:
	1.	Clone the repo

git clone https://github.com/<your-username>/Marmoset-phee-sequence-accommodation.git
cd Marmoset-phee-sequence-accommodation


	2.	Set up Python environment and install dependencies.
	3.	Run notebooks in order:
	1.	notebooks/script_1_data_preprocess.ipynb
	2.	notebooks/script_2_spectral_measurements.ipynb
	3.	notebooks/script_3_seq_measurements.ipynb
	4.	notebooks/script_4_GLM_preprocessing.ipynb
Make sure the input paths at the top of each notebook point to your local data/ folder.
	4.	Run the R analysis
Open R/script_5_R_LMM_model.R in R / RStudio and source it, or run sections interactively. This will:
	•	fit the LMMs for call and sequence structure,
	•	print summaries & likelihood-ratio tests,
	•	generate the model-based plots.

If something fails, it is most likely a path issue (here() vs relative paths) or a package version difference. Feel free to adapt path definitions at the top of each script.

⸻

Citation

If you use this code or pipeline in your own work, please cite the manuscript:

Wewhare, N., Burkart, J. M., Wierucka, K. (2025).
Vocal accommodation in marmoset phee call sequence structure post pairing.
[journal / preprint details to be added]

You’re also welcome to adapt the sequence-metric code (transition matrices, n-grams, local alignment) for other species; if you do, a brief acknowledgement or citation is appreciated.

⸻

Contact

If you have questions, spot bugs, or want to reuse parts of the code in a slightly different context:
	•	Maintainer: Nakul Wewhare
	•	Issues and pull requests via the GitHub repo are very welcome.

I’ve tried to keep the scripts readable and the pipeline transparent; if any step is unclear, open an issue and I’ll add more comments or a small example.
