# Manuscript-to-code crosswalk

This guide maps each analytical claim and figure in the manuscript to the public, reader-facing workflow. The numbered notebooks are intended to be read and run in order. They load versioned cached artifacts by default; expensive distance estimation, neural-network training, and Bayesian sampling are never started implicitly.

## Workflow at a glance

| Manuscript component | Primary notebook | Canonical inputs | Validated artifacts or linked figures |
|---|---|---|---|
| Study design, data processing, session repertoires, and eligibility | `notebooks/01_data_and_session_repertoires.ipynb` | `data/processed/sequences.csv`, `data/processed/calls.csv` | Validated sequence and call session inventories; sample-size audit; locked Figure 1/S1 references |
| Call-structure methods and descriptive results | `notebooks/02_call_structure.ipynb` | STP and MFCC PCA scores, exact-DTW distances, VAE latent means, cached embeddings | Call-distance audit; regenerated analytical Figures S5–S8 plus locked references |
| Sequence-structure methods and descriptive results | `notebooks/03_sequence_structure.ipynb` | Processed sequences, sequence representations and distances, cached embeddings | Sequence-distance audit; regenerated analytical Figures 2, 4, and S2–S4 plus locked references |
| Session-pair construction, distance standardization, and model tables | `notebooks/04_model_tables.ipynb` | Cached repertoire distances | Validated `data/derived/model_sequence.csv` and `data/derived/model_call.csv`; model-input audit |
| Bayesian multilevel models, posterior contrasts, diagnostics, and common support | `notebooks/05_bayesian_models_and_results.ipynb` | Model tables and saved posterior draws | Reconstructed contrast summaries; regenerated Figure 3 when draws are present; diagnostic audit |

## Methods crosswalk

### Shared implementations

The notebooks explain and validate the workflow; shared modules hold the reusable implementation:

| Responsibility | Implementation |
|---|---|
| Unique-call cleaning and full-audio preflight | `src/marmoset_convergence/calls.py` |
| Stable identifiers, eligible repertoires, and Cartesian session pairs | `src/marmoset_convergence/identifiers.py`, `src/marmoset_convergence/sessions.py` |
| Sequence representations and canonical local alignment | `src/marmoset_convergence/sequence.py` |
| Exact constrained DTW and condensed-matrix handling | `src/marmoset_convergence/dtw.py`, `src/marmoset_convergence/matrices.py` |
| Repertoire aggregation, z-scoring, and long-form model tables | `src/marmoset_convergence/model_tables.py` |
| Cache parameter contracts, checksums, and cached/full validation | `src/marmoset_convergence/provenance.py`, `src/marmoset_convergence/cli.py` |
| Cheap analytical plots from validated tables/coordinates | `src/marmoset_convergence/figures.py`, `scripts/regenerate_cached_figures.py` |
| Unified model formula, posterior contrasts, diagnostics, and Figure 3 | `R/model_spec.R`, `R/model_engine.R`, `R/posterior_results.R`, `R/run_models.R` |
| Cached-default orchestration | `scripts/run_repro_pipeline.sh` |

### Study subjects and experimental design

The design is documented in Notebook 01. Six adults (three males and three females) were recorded in all nine possible male–female combinations before pair formation and four months after three permanent pairs were formed. Each combination was scheduled for six sessions at each stage. Analyses compare the same eventual pair members in two audience contexts: Partner sessions, in which they called with one another, and Non-partner sessions, in which each member called with another opposite-sex individual.

Figure 1 is a design schematic rather than a computed statistical result. Its editable source and final export belong in `results/figures/main/`, with provenance recorded alongside the export.

### Data processing and session repertoires

Notebook 01 covers bout construction, call attribution, the 2–10 kHz filtering step, and the definition of a receiver-specific session repertoire. It verifies the manuscript’s sample sizes:

- 1,619 analysed phee sequences and 3,612 unique analysed phee calls;
- 107 eligible sequence repertoires containing 1,496 sequences;
- 130 eligible call repertoires containing 3,527 calls;
- 331 sequence-repertoire comparisons (62 Partner, 269 Non-partner);
- 431 call-repertoire comparisons (74 Partner, 357 Non-partner).

The minimum inclusion threshold is five phee sequences for a sequence repertoire and five individual phee calls for a call repertoire.

### Sequence structure

Notebook 03 covers all four prespecified sequence representations:

1. first-order transition probabilities;
2. normalized bigram frequencies;
3. phee-repeat-length distributions;
4. Smith–Waterman local alignment.

The first three use Euclidean distance between repertoire-level feature vectors. Local alignment is calculated between complete sequences and aggregated to the repertoire level. The exact scoring and normalization rules are fixed in `docs/analysis_specification.md` and in the cache manifest.

Figure 2 summarizes the observed call-type and sequence-length distributions, transition network, and bigram counts. Figure 4 displays the transition-probability repertoire space. Figures S2–S4 display the bigram, local-alignment, and phee-repeat spaces. PaCMAP is used for transition probability and bigram; UMAP is used for local alignment and phee repeat, with the alignment embedding consuming its precomputed distance matrix. Embeddings are descriptive: each sequence metric is embedded once using the 107 eligible repertoires, and the same cached coordinates and axis limits are reused in every stage-by-context panel.

### Call structure

Notebook 02 covers all four prespecified call representations:

1. 11 spectro-temporal parameters (STP), reduced to five principal components;
2. 132 mel-frequency cepstral coefficients (MFCC), reduced to five principal components;
3. exact, constrained spectrogram dynamic time warping (DTW);
4. 32-dimensional variational-autoencoder (VAE) latent means.

The session-pair distance for each call metric is the mean of every call-to-call distance across the two eligible repertoires. The STP and MFCC caches preserve the historical standardized PCA basis used for the reported results (fitted on 3,907 merged rows before audio availability filtering and exact-filename deduplication), while exposing one score row for each of the 3,612 canonical calls. Figures S5–S8 display the MFCC, STP, DTW, and VAE spaces, respectively. Each call embedding contains all 3,612 unique analysed calls; repertoire aggregation later filters to the 3,527 calls in eligible repertoires. PaCMAP is used for vector representations and UMAP with a precomputed distance matrix is used for DTW.

### Model-table preparation

Notebook 04 documents the Cartesian pairing of eligible session repertoires within focal pair, stage, and audience context. Every repertoire-pair comparison contributes one row per metric. Raw distances are standardized separately within each metric, jointly across pairs, stages, and contexts. The expected long-form model tables contain 1,324 sequence rows and 1,724 call rows.

### Bayesian modelling and posterior contrasts

Notebook 05 documents the two unified Gaussian multilevel models, one for sequence structure and one for call structure. In each model, metric interacts with stage and context; the same model produces the combined and metric-specific posterior contrasts. Pair-level varying effects and a multi-membership session term account for repeated dyads and repeated use of session repertoires.

Figure 3 is generated from saved posterior draws. It shows the stage contrast

\[
\Delta = E(d_\mathrm{after}) - E(d_\mathrm{before})
\]

for Partner and Non-partner contexts in separate rows, with one combined and
four metric-specific columns for each structural level. The direct audience
contrast is intentionally reported in the text and posterior tables—not as an
additional Figure 3 panel—and is

\[
\Psi = \Delta_\mathrm{Partner} - \Delta_\mathrm{Non\mbox{-}partner}.
\]

Negative values of \(\Delta\) indicate convergence. Positive values of \(\Psi\) indicate a larger reduction in distance in the Non-partner context. The common-support sensitivity estimate averages over the two focal pairs represented in every stage-by-context cell.

## Results crosswalk

| Reported result | Reproducibility location |
|---|---|
| Eligible-repertoire and comparison counts | Notebooks 01 and 04 |
| Sequence convergence in the Non-partner context | Notebook 05; sequence posterior-draw cache; Figure 3 |
| No resolved sequence change in the Partner context | Notebook 05; sequence posterior-draw cache; Figure 3 |
| Positive stage-by-context sequence contrast and common-support result | Notebook 05; contrast table and posterior draws |
| Three of four sequence metrics show Non-partner convergence | Notebook 05; metric-specific contrasts; Figure 3 |
| No resolved call-structure change in either context | Notebook 05; call posterior-draw cache; Figure 3 |
| Sequence-space descriptive patterns | Notebook 03; Figures 4 and S2–S4 |
| Call-space descriptive patterns | Notebook 02; Figures S5–S8 |
| Sampling and posterior-predictive diagnostics | Notebook 05; `results/tables/` |

## Figure provenance

The PNGs named `figure_1.png` through `figure_4.png` and `figure_s1_*.png` through `figure_s8_*.png` are provenance-locked snapshots of the current manuscript figures. `results/manifest.json` and `results/figures/figure_provenance.csv` record their document revision, source object, and SHA-256. They are visual reference targets, not numerical inputs to the analysis.

Displaying one of these locked snapshots is not, by itself, a reproduction. If the quantitative cache needed by a notebook is unavailable, the notebook reports that absence and does not treat the snapshot as computed evidence.

Regenerated analytical exports read processed tables, derived distances, cached coordinates, or saved posterior draws. They do not fit a model, train a VAE, or calculate DTW as a side effect. Figure panels based on an embedding always reuse the versioned pooled coordinates in `data/cache/embeddings/` and the same full-data limits in every facet; they are not embedded independently. The Python figures can be written with `python scripts/regenerate_cached_figures.py --root .` (or selected with repeated `--figure` flags). Regenerated files use a distinct name, such as `figure_4_regenerated.png`, so they cannot silently overwrite a locked snapshot.

The compact regenerated call-space panels honestly omit the manuscript's linked spectrogram cards: complete WAV coverage is not part of the cached figure contract. Sequence-space panels similarly omit the locked repertoire cards. These `_regenerated.png` files are disposable reader renders, not quantitative inputs or archival snapshots, so they are displayed directly after creation and are not accepted by `registered_artifact()` as substitutes for a manifest-registered cache or locked PNG.

The quantitative source for every analytical figure must be recoverable from the paths named in this crosswalk and the relevant artifact manifest. If a required artifact is absent or fails validation, the notebooks stop with an explanatory error rather than producing a substitute result.
