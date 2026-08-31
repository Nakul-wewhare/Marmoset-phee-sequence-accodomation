# Data dictionary

This document describes the public analytical data model. Column order is not semantically meaningful; identifiers, units, categories, and relationships are. CSV files use UTF-8 encoding, a header row, and explicit missing values. Identifier columns should be read as strings even when their values appear numeric.

## Controlled vocabulary

| Concept | Allowed values | Meaning |
|---|---|---|
| `stage` | `before`, `after` | Recording stage before pair formation or four months after pair formation |
| `context` | `partner`, `non-partner` | Identity of the receiver relative to the caller’s eventual pair assignment |
| `structure` | `sequence`, `call` | Level of vocal organization represented by a distance |
| sequence `metric` | `transition_probability`, `bigram`, `phee_repeat`, `local_alignment` | Sequence representation used to calculate distance |
| call `metric` | `stp`, `mfcc`, `dtw`, `vae` | Acoustic representation used to calculate distance |

`Non-partner` is the preferred manuscript label. Historical labels such as `stranger`, `non_partner`, or `non-paired` are normalized to `non-partner` when legacy inputs are imported.

## Identifier relationships

- `individual_id` identifies a caller.
- `receiver_id` identifies the conspecific present across the barrier.
- `recording_id` identifies a dyadic recording event.
- `repertoire_id` identifies the calls or sequences produced by one caller within one recording event. It must therefore include caller identity, not only a session number.
- `sequence_id` identifies one annotated bout within a repertoire.
- `call_id` identifies one unique phee-call selection. The canonical call table contains one row per `call_id`.
- `pair_id` identifies one of the three eventual bonded pairs.
- `comparison_id` identifies an unordered pair of eligible repertoires belonging to the two members of one focal pair in the same stage and audience context.

Every derived comparison records both `repertoire_a_id` and `repertoire_b_id`. These identifiers are also used by the multi-membership random effect in the Bayesian models.

## `data/processed/sequences.csv`

One row represents one annotated phee sequence. The canonical table contains 1,619 rows.

| Field | Type | Description |
|---|---|---|
| `sequence_id` | string | Stable unique sequence identifier |
| `repertoire_id` | string | Receiver-specific caller-session repertoire containing the sequence |
| `recording_id` | string | Dyadic recording-event identifier |
| `individual_id` | string | Caller identity |
| `receiver_id` | string | Audience individual |
| `pair_id` | string | Eventual bonded pair to which the caller belongs |
| `stage` | category | `before` or `after` |
| `context` | category | `partner` or `non-partner` |
| `sequence` | string | Ordered call-type symbols after cleaning |
| `n_elements` | integer | Number of vocal elements in `sequence` |
| `session_number` | string | Source session label; not globally unique by itself |
| `recording_datetime` | datetime | Recording date and time when available |
| `source_file` | string | Repository-relative provenance for the annotation source |

Additional provenance columns may be retained when they do not expose local machine paths. The ordered `sequence` is the analytical value; wide legacy element columns are not required for new analysis code.

## `data/processed/calls.csv`

One row represents one unique phee call. Duplicate legacy rows are removed by stable call identity before this table is written. The canonical table contains 3,612 rows.

| Field | Type | Unit / description |
|---|---|---|
| `call_id` | string | Stable unique call identifier |
| `filename` | string | Portable WAV basename used to derive and audit `call_id` |
| `repertoire_id` | string | Receiver-specific caller-session repertoire containing the call |
| `recording_id` | string | Dyadic recording-event identifier |
| `individual_id` | string | Caller identity |
| `receiver_id` | string | Audience individual |
| `pair_id` | string | Eventual bonded pair to which the caller belongs |
| `stage` | category | `before` or `after` |
| `context` | category | `partner` or `non-partner` |
| `session_number` | string | Source session label; not globally unique by itself |
| `recording_datetime` | datetime | Recording date and time when available |
| `audio_path` | string | Repository-relative WAV path |
| `duration_90_s`, `duration_50_s` | numeric | Energy-duration measures in seconds |
| `center_frequency_hz` | numeric | Center frequency in hertz |
| `frequency_5_hz`, `frequency_25_hz`, `frequency_75_hz`, `frequency_95_hz` | numeric | Energy-frequency percentiles in hertz |
| `bandwidth_50_hz`, `bandwidth_90_hz` | numeric | Energy bandwidths in hertz |
| `average_entropy_bits`, `aggregate_entropy_bits` | numeric | Spectral entropy measures in bits |
| `mfcc_001` … `mfcc_132` | numeric | MFCC features before standardization and PCA |

## Session inventories

`data/cache/sequence_session_inventory.csv` contains one row per eligible sequence repertoire; `data/cache/call_session_inventory.csv` contains one row per eligible call repertoire. Both include `repertoire_id`, `individual_id`, `receiver_id`, `pair_id`, `stage`, `context`, and `session_number`. Eligibility has already been applied before these caches are written, so there is no separate `eligible` flag.

| File-specific field | Type | Description |
|---|---|---|
| `session_number` | string | Source session label within the caller/receiver/stage key |
| `n_sequences` | integer | Number of retained phee sequences in a sequence repertoire |
| `n_calls` | integer | Number of retained phee calls in a call repertoire |

Expected inventory totals are 107 sequence repertoires containing 1,496 sequences and 130 call repertoires containing 3,527 calls.

## Representation caches

| Path | Contents |
|---|---|
| `data/cache/sequence_representations.npz` | Ordered repertoire identifiers and feature matrices for transition probabilities, bigram proportions, and phee-repeat probabilities |
| `data/cache/sequence_distances.npz` | Ordered repertoire identifiers and the 5,671-value SciPy-order condensed local-alignment distance vector |
| `data/cache/stp_pca_scores.csv` | `call_id` plus five STP principal-component scores |
| `data/cache/mfcc_pca_scores.csv` | `call_id` plus five MFCC principal-component scores |
| `data/cache/dtw_distances_condensed.npz` | Ordered `call_id` values and the condensed upper triangle of the exact-DTW matrix |
| `data/cache/vae_latent_means.csv` | `call_id` plus 32 latent-mean coordinates |
| `data/cache/embeddings/<metric>.csv` | Exactly the observation identifier, `embedding_x`, `embedding_y`, `individual_id`, `pair_id`, `stage`, and `context`, in canonical source order |

The condensed DTW vector contains \(3,612 \times 3,611 / 2 = 6,521,466\) distances. A condensed vector excludes the zero diagonal and stores each unordered call pair once.

The STP/MFCC score caches preserve the PCA bases used by the reported results.
Those bases were fitted to the 3,907-row historical merged selection table,
before three rows lacking exported audio were omitted and before 3,904 exported
rows were collapsed to 3,612 unique filenames. The manifest records this
provenance; the canonical cache preparation step does not refit PCA.

## Repertoire-distance tables

`data/derived/sequence_repertoire_distances.csv` and `data/derived/call_repertoire_distances.csv` record one repertoire comparison per metric, preferably in long form.

| Field | Type | Description |
|---|---|---|
| `comparison_id` | string | Stable unique repertoire-pair identifier |
| `pair_id` | string | Eventual bonded pair |
| `stage` | category | `before` or `after` |
| `context` | category | `partner` or `non-partner` |
| `repertoire_a_id`, `repertoire_b_id` | string | The two eligible repertoires compared |
| `metric` | category | One of the four metrics for the structural level |
| `distance` | numeric | Raw non-negative distance; smaller values denote greater similarity |
| `n_lower_level_pairs` | integer | Number of sequence-to-sequence or call-to-call comparisons contributing to an aggregate distance, when applicable |

There are 331 unique sequence comparisons (62 Partner and 269 Non-partner) and 431 unique call comparisons (74 Partner and 357 Non-partner).

## Model tables

`data/derived/model_sequence.csv` and `data/derived/model_call.csv` are long-form analysis tables with one row per repertoire comparison and metric. They retain the identifiers above and add:

| Field | Type | Description |
|---|---|---|
| `standardized_distance` | numeric | Distance z-score calculated within metric across all eligible pairs, stages, and contexts |
| `stage_centered` | numeric | `-0.5` before pairing and `+0.5` after pairing |

The sequence model table has 1,324 rows (331 comparisons × four metrics). The call model table has 1,724 rows (431 comparisons × four metrics).

## Posterior draws and summaries

`results/posterior_draws/sequence_stage_effect_draws.csv.gz` and `results/posterior_draws/call_stage_effect_draws.csv.gz` contain pair-level posterior expected values. Each row is one posterior draw × focal pair × audience context × metric combination.

Each portable draw file has an adjacent `.csv.gz.metadata.json` sidecar. The
sidecar records and fingerprints the exact model-table file/content SHA-256,
formula, metric and pair order, common-support pairs, sampling seed/settings,
software versions, and draw checksum. Report mode rejects a missing, edited, or
stale sidecar.

| Field | Type | Description |
|---|---|---|
| `.draw` | integer | Posterior draw identifier; 8,000 post-warm-up draws are expected |
| `analysis` | category | Sequence structure or call structure |
| `pair` | string | Focal bonded pair included in the prediction |
| `context` | category | Partner or Non-partner |
| `metric` | category | Canonical metric code |
| `metric_label` | string | Reader-facing metric label |
| `expected_before` | numeric | Posterior expected standardized distance before pairing, retaining the pair effect and marginalizing the session multi-membership effect |
| `expected_after` | numeric | Corresponding posterior expectation after pairing |
| `delta` | numeric | `expected_after - expected_before` |
| `common_support_pair` | boolean | Whether this pair has observed comparisons in all four stage-by-context cells |

The raw portable draw files do not contain pre-collapsed support sets, combined metrics, or Ψ rows. Report mode derives `all_pairs` and `common_support_A_B`, averages pairs and metrics with the specified equal weights, and calculates Ψ from these pair-level rows.

Tables in `results/tables/` summarize posterior medians, 95% credible intervals, sample sizes, and diagnostics. The diagnostics table includes registered `ppc_path` and `trace_path` artifacts for both models and requires bulk and tail ESS of at least 1,000. Posterior draw files—not rounded summary tables—are the source for Figure 3.

## Artifact manifests

Every cache is described either by an adjacent JSON manifest or by
`data/cache/manifest.json`. The core schema requires a portable path, artifact
SHA-256, byte size, and parameter object, with an ordered-identifier hash and
count for order-sensitive arrays. File-specific parameters record available
dimensions and input provenance. Costly or stochastic release artifacts must
add the algorithm/package versions, random seed, source hashes, and creation
provenance needed to reproduce that calculation. Array order is meaningful and
must agree exactly across identifiers, matrices, and embedding tables.
