# Reproducibility guide

The public workflow is cache-first. The normal command verifies artifact
checksums, ordered identifiers, schemas, method parameters, and notebook
hygiene. It does not calculate a missing distance, fit an embedding, train a
network, or sample a model as a fallback.

## Environment

Run commands from the repository root. The Python file pins the direct
dependencies used by cached validation and the notebooks; the repository must
also be installed so scripts can import the shared package:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r environment/python-requirements.txt
python -m pip install -e .
Rscript environment/install_R_packages.R
```

`environment/install_R_packages.R` records the direct R packages but is not a
binary lockfile. An archival release should preserve its resolved R session
information alongside the model cache. Full acoustic or embedding work uses
the separate, optional `environment/python-full-rebuild-requirements.txt`.
Approximate FastDTW is intentionally absent from both environments.

## Cached validation: default

```bash
bash scripts/run_repro_pipeline.sh
```

This is the same as `--mode cached`. It is read-only and validates every
artifact registered in the checkout. Because this development branch has a
partial cache bundle, the command succeeds while printing the precise list of
current-analysis assets that are still missing.

To audit an archival bundle and fail unless every required current artifact is
installed, use:

```bash
bash scripts/run_repro_pipeline.sh --mode cached --strict
```

The five notebooks apply the same rule. Their expensive-step flags are all
`False`, and no notebook contains a fallback that enables a costly operation.
Notebook 01 can be executed with the included bundle; Notebooks 02–05 require
additional quantitative release caches and stop before analysis if one is
missing or unregistered.

## Reader-facing notebook order

```text
notebooks/01_data_and_session_repertoires.ipynb
notebooks/02_call_structure.ipynb
notebooks/03_sequence_structure.ipynb
notebooks/04_model_tables.ipynb
notebooks/05_bayesian_models_and_results.ipynb
```

The notebooks find the repository root without a user-specific path. Outputs
and execution counts are stripped in version control. They explain and audit
the shared implementation rather than carrying separate copies of the
calculation code.

## Included inexpensive artifacts

The following files are generated deterministically by
`scripts/prepare_cached_artifacts.py` without reading audio or refitting PCA:

```text
data/cache/sequence_session_inventory.csv
data/cache/call_session_inventory.csv
data/cache/sequence_session_pairs.csv
data/cache/call_session_pairs.csv
data/cache/sequence_representations.npz
data/cache/stp_pca_scores.csv
data/cache/mfcc_pca_scores.csv
data/cache/manifest.json
```

The STP and MFCC files reuse the five score columns present in the historical
processed call table after canonical filename deduplication. Their basis was
originally fitted to 3,907 merged rows, before three unavailable-audio rows
were omitted and before the 3,904 exported rows were collapsed to 3,612 unique
calls. The script records that fit population, component variances, source
columns, and `pca_refit = false` in the manifest. Its default action never
computes local alignment, DTW, VAE features, embeddings, or models.

To rebuild only these inexpensive artifacts explicitly:

```bash
python scripts/prepare_cached_artifacts.py
```

Once the released exact-DTW condensed matrix and VAE latent means are present
and registered, assemble the four-metric call tables with:

```bash
python scripts/assemble_call_model_table.py --root .
# Equivalent runner opt-in:
bash scripts/run_repro_pipeline.sh --assemble-call-model
```

The assembler uses STP/MFCC Euclidean distance, VAE cosine distance, and direct
lookup from the exact-DTW cache. It averages every cross-repertoire call pair,
records that Cartesian-product count, applies population-SD standardization,
and registers both outputs. It never reads WAVs, refits PCA, calculates DTW, or
trains a VAE. A missing, unregistered, tampered, misordered, or
parameter-incompatible cache stops the command before either table is written.
The default runner invocation does not call the assembler and remains
read-only.

## Complete cached-release contract

Running all five notebooks and reconstructing the statistical summaries also
requires these versioned assets:

```text
data/cache/sequence_distances.npz
data/cache/dtw_distances_condensed.npz
data/cache/vae_latent_means.csv
data/cache/embeddings/transition_probability.csv
data/cache/embeddings/bigram.csv
data/cache/embeddings/phee_repeat.csv
data/cache/embeddings/local_alignment.csv
data/cache/embeddings/stp.csv
data/cache/embeddings/mfcc.csv
data/cache/embeddings/dtw.csv
data/cache/embeddings/vae.csv
data/derived/sequence_repertoire_distances.csv
data/derived/call_repertoire_distances.csv
data/derived/model_sequence.csv
data/derived/model_call.csv
results/posterior_draws/sequence_stage_effect_draws.csv.gz
results/posterior_draws/sequence_stage_effect_draws.csv.gz.metadata.json
results/posterior_draws/call_stage_effect_draws.csv.gz
results/posterior_draws/call_stage_effect_draws.csv.gz.metadata.json
results/tables/stage_effects_by_context_metric.csv
results/tables/stage_effects_combined.csv
results/tables/psi_by_metric.csv
results/tables/psi_combined.csv
results/tables/model_diagnostics.csv
results/figures/main/diagnostics/sequence/sequence_interaction_model_pp_check.png
results/figures/main/diagnostics/sequence/sequence_interaction_model_trace.png
results/figures/main/diagnostics/call/call_interaction_model_pp_check.png
results/figures/main/diagnostics/call/call_interaction_model_trace.png
results/model_manifest.json
```

The RDS draw caches and `data/cache/posterior_draws/model_diagnostics.csv` may
also be retained as faster local alternatives, but they are not required when
the portable draw CSVs, their metadata sidecars, and the registered results
diagnostics are complete.

The cache and result manifests must register each artifact’s portable path and
SHA-256. Matrix records additionally bind the exact method parameters and
ordered identifier hash. Posterior metadata binds the draws to the exact model
table content, formula, factor levels, pair grid, sampling settings, software,
and draw checksum. A rounded manuscript table is never a posterior cache.

The current development checkout intentionally lacks several files in this
list. The older 2,200-call FastDTW matrix and older context-separated models do
not satisfy the contract.

## Explicit full mode

Full mode never selects a costly component implicitly. Sequence local
alignment and Bayesian refitting validate their own non-audio inputs and do
not require WAV availability:

```bash
bash scripts/run_repro_pipeline.sh --mode full \
  --recompute-sequence-alignment

bash scripts/run_repro_pipeline.sh --mode full --refit-models
```

Passing full mode alone still performs validation only. Before a future
audio-dependent DTW or VAE rebuild, separately require one WAV basename for
every canonical `call_id`:

```bash
marmoset-repro preflight-audio --root . \
  --audio-dir=/path/to/complete/extracted-calls
```

The tracked 728-file subset fails this identity-only check before any signal
processing starts.

It uses four chains, 4,000 iterations per chain, 2,000 warm-up iterations, the
fixed priors and interaction structure, mandatory trace/posterior-predictive
diagnostics, and minimum bulk and tail effective sample sizes of 1,000. A
release refit stops instead of publishing a cache when any required diagnostic
fails. The
safe default `Rscript R/run_models.R --report-only` loads validated draws and
does not sample.

The exact constrained alignment kernel and all-pairs frame-array driver are implemented by
`marmoset_convergence.dtw`, and all-pairs execution is disabled unless
`allow_expensive=True`. The current repository validates and reuses the
released DTW cache; it does not yet expose an end-to-end WAV-to-spectrogram DTW
command. Chatter and embedding parameters are fixed in
`marmoset_convergence.provenance`. A complete archival release must pair every
large artifact with its producing driver and resolved environment; this branch
does not pretend the legacy approximate matrix is a current precomputation.

## Expected validation totals

| Quantity | Expected value |
|---|---:|
| Analysed phee sequences | 1,619 |
| Unique analysed phee calls | 3,612 |
| Eligible sequence repertoires | 107 |
| Sequences in eligible repertoires | 1,496 |
| Eligible call repertoires | 130 |
| Calls in eligible repertoires | 3,527 |
| Unique sequence comparisons | 331 |
| Partner sequence comparisons | 62 |
| Non-partner sequence comparisons | 269 |
| Sequence model rows | 1,324 |
| Unique call comparisons | 431 |
| Partner call comparisons | 74 |
| Non-partner call comparisons | 357 |
| Call model rows | 1,724 |
| Exact-DTW condensed distances | 6,521,466 |

## Figure snapshots and regeneration

The unsuffixed PNGs in `results/figures/` are checksum-locked exports from the
identified current manuscript revision. They are comparison targets, not
computational inputs. Regenerated exports use a distinct `_regenerated` suffix
and never overwrite them.

After installing the cached-mode environment, regenerate any Python-owned
analytical figure without executing a notebook:

```bash
python scripts/regenerate_cached_figures.py --root . --figure 2
python scripts/regenerate_cached_figures.py --root . --figure 4 --figure s2 --figure s3 --figure s4
python scripts/regenerate_cached_figures.py --root . --figure s5 --figure s6 --figure s7 --figure s8
```

Omitting `--figure` requests all Python-owned figures. The command validates
every requested input before its first write and stops if a checksum, exact ID
order, seven-column embedding schema, coordinate value, grouping value,
method, package version, random seed, pooled-fit scope, or source fingerprint
does not match the manifest. Existing regenerated files are replaced only with
the explicit `--overwrite-regenerated` flag; unsuffixed locked PNGs are never
writable targets.

- Figure 2 can be reconstructed from the processed sequence table.
- Figure 3 requires validated posterior contrast draws and diagnostics.
- Figure 4 requires the cached transition-probability embedding.
- Figures S2–S4 require the cached sequence embeddings.
- Figures S5–S8 require the cached call embeddings.

Each sequence embedding is fitted once to the 107 eligible repertoires, and
each call embedding is fitted once to all 3,612 analysed unique calls. All
pooled and stage-by-context panels reuse the same coordinates and limits.
Statistical conclusions use original-space distances, never two-dimensional
display distances.

The regenerated call panels omit the locked spectrogram-card annotations
because the complete WAV collection is not in the cached release. The compact
sequence panels omit the locked repertoire cards. These `_regenerated.png`
files are disposable reader-facing renders rather than release inputs or
archival manuscript snapshots; notebooks display a file immediately after
creating it, but later cached analyses do not accept it through
`registered_artifact()` unless a release curator explicitly registers a new
presentation artifact.

## Verification checklist for an archival release

1. Start from a clean clone and install the pinned Python requirements plus the
   local package.
2. Install the large release assets at their canonical relative paths.
3. Run `bash scripts/run_repro_pipeline.sh --mode cached --strict`.
4. Execute the five notebooks in numeric order; their figure cells write only
   `_regenerated` exports.
5. Confirm every expected count and within-metric population z-score audit.
6. Confirm the model diagnostics, input fingerprints, 8,000 post-warm-up draws,
   equal-pair/equal-metric contrasts, and common-support calculation.
7. Compare regenerated figures with the locked snapshots without overwriting
   either version.
8. Archive the repository commit, manifests, full Python/R environments, and
   platform metadata.

## Troubleshooting

**A cache is missing.** Obtain the complete release asset. Do not substitute an
approximation, an old matrix, a figure pixel value, or a manuscript number.

**An identifier-order check fails.** Rebuild every dependent artifact from the
same ordered index. Reordering only a matrix invalidates pairwise lookup.

**Fewer than 3,612 WAVs are available.** Cached analysis remains possible, but
DTW and VAE caches cannot be rebuilt from that directory.

**Posterior metadata does not match a model table.** Treat the draw cache as
stale and obtain or explicitly refit the model from the exact table. Matching
pair names alone is not sufficient provenance.
