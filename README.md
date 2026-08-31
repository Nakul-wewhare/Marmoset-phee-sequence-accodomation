# Sequence-level vocal convergence in common marmosets

This repository is the reader-facing analysis companion to *Sequence-level
vocal convergence in common marmosets* (Wewhare, Burkart, and Wierucka). It
organizes the current methods into five numbered notebooks, keeps reusable
calculations in tested Python and R modules, and treats costly results as
versioned caches.

The default workflow is read-only and inexpensive. It verifies the processed
tables, available caches, method parameters, identifier order, and SHA-256
checksums. It never calculates the all-call DTW matrix, trains the VAE, fits an
embedding, or samples a Bayesian model as a fallback. The implemented costly
components require explicit full-mode flags; unavailable drivers are named as
release limitations rather than simulated.

## Start here

Create an environment and install the local package:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r environment/python-requirements.txt
python -m pip install -e .
```

Then validate every artifact currently included in the checkout:

```bash
bash scripts/run_repro_pipeline.sh
```

Open the notebooks in this order:

1. `notebooks/01_data_and_session_repertoires.ipynb`
2. `notebooks/02_call_structure.ipynb`
3. `notebooks/03_sequence_structure.ipynb`
4. `notebooks/04_model_tables.ipynb`
5. `notebooks/05_bayesian_models_and_results.ipynb`

Each notebook begins with the corresponding Methods/Results/Figure crosswalk.
Comments explain the analytical decisions and cache contracts needed by a new
reader. Notebook outputs are stripped, so readers see only results produced
from their own validated local artifacts.

## What is included in this branch

- Canonical processed tables: 1,619 sequences and 3,612 deduplicated calls.
- Checksummed eligible-session inventories: 107 sequence repertoires and 130
  call repertoires.
- Canonical session-pair indexes: 331 sequence comparisons and 431 call
  comparisons.
- Inexpensive sequence representations and the five previously calculated STP
  and MFCC principal-component scores for every unique call.
- Exact method implementations for sequence metrics, the constrained-DTW
  alignment kernel, model-table construction, cache validation, and the
  unified Bayesian models.
- Provenance-locked snapshots of Figures 1–4 and S1–S8 from the current
  manuscript revision. These are visual reference targets, not numerical
  substitutes for missing caches.

This development branch is intentionally honest about incomplete release
assets. The current exact-DTW matrix, VAE latent means, descriptive embedding
coordinates, complete call model table, and posterior draws are not present in
the source material available for this refresh. Older 2,200-call FastDTW files
and context-separated three-metric models are methodologically incompatible and
were not relabeled as current results. Their required canonical paths are listed
in `data/cache/manifest.json` and `docs/reproducibility.md`.

Notebook 01 can run from the included assets. Later notebooks stop with a clear
missing-cache message until the corresponding versioned release assets are
installed. A manuscript figure snapshot is never accepted in place of its
quantitative source.

## Analysis map

| Stage | Reader notebook | Shared implementation |
|---|---|---|
| Data cleaning, identity, eligibility, session pairs | `01_data_and_session_repertoires.ipynb` | `scripts/import_legacy_processed_data.py`, `marmoset_convergence.sessions` |
| STP, MFCC, exact DTW, VAE call structure | `02_call_structure.ipynb` | `marmoset_convergence.calls`, `marmoset_convergence.dtw` |
| Transition, bigram, phee repeat, local alignment | `03_sequence_structure.ipynb` | `marmoset_convergence.sequence` |
| Repertoire aggregation and within-metric standardization | `04_model_tables.ipynb` | `marmoset_convergence.model_tables` |
| Unified models, contrasts, diagnostics, Figure 3 | `05_bayesian_models_and_results.ipynb` | `R/model_spec.R`, `R/model_engine.R`, `R/posterior_results.R` |

The full manuscript-to-code mapping is in
`docs/manuscript_code_crosswalk.md`; exact mathematical and computational rules
are in `docs/analysis_specification.md`.

## Cached and full modes

Cached mode validates only; it does not silently repair or recompute an
artifact:

```bash
bash scripts/run_repro_pipeline.sh --mode cached
```

To require a complete current release rather than validate the available
partial bundle:

```bash
bash scripts/run_repro_pipeline.sh --mode cached --strict
```

After the exact-DTW and VAE caches have been installed and registered, the two
call tables can be assembled explicitly without audio or feature recomputation:

```bash
bash scripts/run_repro_pipeline.sh --assemble-call-model
```

This opt-in action verifies every input checksum, method contract, and ordered
call ID, writes the two derived CSVs, and registers their checksums. Omitting
the flag leaves cached mode read-only.

Full mode is opt-in, and every expensive component still requires its
individual flag. Sequence alignment and Bayesian refitting do not depend on
the audio checkout:

```bash
bash scripts/run_repro_pipeline.sh --mode full \
  --recompute-sequence-alignment

# Bayesian sampling is also explicit:
bash scripts/run_repro_pipeline.sh --mode full --refit-models
```

Do not use the 728 WAV files currently tracked in `extracted calls/` for a full
acoustic rebuild; they are an incomplete historical subset. Cached validation
does not read audio. Before any audio-dependent rebuild, run the separate
identity-only preflight:

```bash
marmoset-repro preflight-audio --root . \
  --audio-dir=/path/to/complete/extracted-calls
```

## Repository layout

- `data/source/`: historical processed inputs retained for provenance.
- `data/processed/`: canonical, portable sequence and call tables.
- `data/cache/`: versioned representations, matrices, scores, and manifest.
- `data/derived/`: repertoire-distance and model-ready tables.
- `src/marmoset_convergence/`: reusable Python implementation.
- `R/`: unified Bayesian model and posterior-reporting implementation.
- `notebooks/`: five reader-facing, cached-only narratives.
- `results/`: locked manuscript snapshots, generated outputs, and reference
  values kept separate by role.
- `docs/`: analysis specification, data dictionary, crosswalk, and detailed
  reproduction guide.
- `tests/`: lightweight Python and R tests; no test runs costly components.

## Citation and licensing

Citation metadata is in `CITATION.cff`. Code is released under the MIT License;
data, audio, figures, and derived outputs are released under CC BY 4.0, as
described in `LICENSE`, `LICENSE-CODE`, and `LICENSE-DATA`. The archival DOI can
be added to the citation metadata after a complete, checksummed bundle is
deposited.
