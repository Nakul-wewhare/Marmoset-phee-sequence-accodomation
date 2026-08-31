# Canonical analysis specification

This file is the human-readable specification for the manuscript analysis. Machine-readable parameters and exact ordered identifiers are recorded in the cache manifest. Changing a rule below creates a new analysis version and requires rebuilding every dependent artifact.

## 1. Experimental units and eligibility

The six focal adults comprise three eventual male–female pairs. Recordings occurred before pair formation and four months after pair formation. Audience context is defined relative to eventual pair membership at both stages.

A phee sequence is a bout containing at least one phee element. Consecutive vocal elements separated by no more than 0.5 s belong to the same bout. The documented edge rule may attach an isolated element within 1 s of a bout boundary when the separating silence is comparable to within-bout intervals. Only confidently attributed, non-overlapping calls with adequate signal-to-noise ratio are retained. Audio preprocessing applies a 2–10 kHz bandpass filter.

A session repertoire is every retained unit produced by one caller during one recording session with one receiver. Sequence analyses retain repertoires containing at least five phee sequences. Call analyses retain repertoires containing at least five individual phee calls.

For a focal bonded pair, stage, and context, every eligible repertoire from one partner is compared with every eligible repertoire from the other partner. In the Non-partner context, the two receivers need not be the same individual. A repertoire is never compared across stages or contexts.

## 2. Sequence metrics

All call types within an eligible phee sequence are retained. The manifest fixes the call-type alphabet and feature order.

### 2.1 Transition probability

For every repertoire, count adjacent transitions from call type \(i\) to call type \(j\). Normalize each source-call row to obtain the first-order conditional transition matrix \(T\), where \(T_{ij}=P(j\mid i)\). Flatten the matrix in the manifest’s fixed call-type order. Repertoire distance is Euclidean distance between flattened matrices.

### 2.2 Bigram distribution

Count every adjacent ordered call pair in the repertoire and divide by the total number of observed bigrams. With six retained call-type categories, the fixed representation has 36 entries. Repertoire distance is Euclidean distance between bigram-proportion vectors.

### 2.3 Phee-repeat distribution

Identify contiguous runs of phee elements and calculate the relative-frequency vector for the prespecified exact run-length categories 2, 3, 4, and 5. Runs longer than five are excluded rather than truncated or pooled into the final bin. A repertoire with no eligible run is represented by an all-zero vector. Repertoire distance is Euclidean distance between repeat-probability vectors.

### 2.4 Local alignment

Use Smith–Waterman local alignment for every cross-repertoire pair of complete call-type sequences. The canonical scores are:

| Operation | Score |
|---|---:|
| Match | +2 |
| Mismatch | −1 |
| Gap | −1 |

For sequences \(x\) and \(y\), divide the optimal local-alignment score by \(\log(|x|+|y|)\). For a pair of session repertoires, average this normalized similarity over **all** sequence pairs in their Cartesian product. Let \(s_{ab}\) denote that repertoire similarity and let \(s_{\max}\) be the maximum across all 5,671 unordered pairs among the 107 eligible repertoires, including pairs that are not among the 331 selected model comparisons. The stored non-negative distance is

\[
d_{ab}=s_{\max}-s_{ab}.
\]

The value of \(s_{\max}\), the eligible comparison set over which it was calculated, sequence counts, and implementation version are recorded in the cache manifest. The symmetric 107-repertoire matrix is stored as its 5,671-value upper triangle in SciPy condensed order (`storage = "condensed_scipy"`), excluding the zero diagonal. This definition does not use directional best matches and does not substitute an approximate alignment algorithm.

## 3. Call metrics

Call analyses use only individual phee calls occurring within sequences. The canonical processed table has 3,612 unique calls.

### 3.1 Spectro-temporal parameters (STP)

Extract 11 Raven energy and duration measures using a 400 ms Hann window, 200
ms hop, and 50% overlap: two duration measures, center frequency, four
frequency percentiles, two bandwidth measures, and two entropy measures. The
frozen PCA basis underlying the reported results was fitted during historical
preprocessing: `StandardScaler` and five-component PCA were fitted to the 3,907
rows in the merged retained-selection table. Three rows without exported audio
were then omitted, producing 3,904 exported rows; exact repeated filenames were
subsequently collapsed to the 3,612 canonical unique calls. The five retained
components explained 97.88% of variance. The public cache reuses those saved
scores rather than refitting a slightly different basis after deduplication.
Call distance is Euclidean distance in this five-dimensional score space.

### 3.2 Mel-frequency cepstral coefficients (MFCC)

Extract the canonical 132-coefficient MFCC representation and fit a separate
standardized PCA using the same frozen 3,907-row historical fit population and
the same subsequent 3,904-to-3,612 export/deduplication sequence described
above. The five retained components explained 79.12% of variance. Call
distance is Euclidean distance in this five-dimensional score space. MFCC
feature order and the PCA-score provenance are fixed in the manifest. Refitting
either PCA on only 3,612 unique calls defines a new analysis version and
invalidates dependent distances, embeddings, and posterior results.

### 3.3 Exact dynamic time warping (DTW)

Construct a fixed-reference spectrogram for each call with these parameters:

| Parameter | Value |
|---|---:|
| Frequency range | 4–12 kHz |
| Analysis window | 2,048-sample Hann |
| Hop | 512 samples |
| Display/amplitude range | −80 to 0 dBFS |
| Frame-to-frame cost | Cosine distance |
| Algorithm | Exact dynamic programming |
| Allowed normalized-time deviation | 0.20 |
| Path normalization | Cumulative cost divided by alignment-step count |

For spectrograms with \(n\) and \(m\) frames, a frame pair \((i,j)\) is admissible only when

\[
\left|\frac{i}{n-1}-\frac{j}{m-1}\right|\leq 0.20.
\]

Use the standard monotone DTW predecessor steps (advance the first call, the second call, or both) within this band. The result is symmetric and stored once in SciPy condensed order (`storage = "condensed_scipy"`) following the manifest’s ordered `call_id` list. The vector has 6,521,466 entries for 3,612 calls. Approximate DTW implementations are not valid substitutes.

### 3.4 Variational autoencoder (VAE)

Train the Chatter 0.1.6 convolutional VAE once on all 3,612 unique calls. The fixed training configuration is:

| Parameter | Value |
|---|---:|
| Input | 128 × 128 spectrogram image |
| Time canvas | 2.5 s, centred |
| Latent dimension | 32 |
| Epochs | 100 |
| Batch size | 32 |
| Learning rate | \(10^{-4}\) |
| \(\beta\) | 0.5 |
| Random seed | 42 |

The cached representation is the 32-dimensional latent mean for each call. Call distance is cosine distance between latent-mean vectors. Training labels do not include caller, stage, or audience context.

### 3.5 Repertoire aggregation

For each call metric, calculate every call-to-call distance across the Cartesian product of the two session repertoires and take the arithmetic mean. Retain `n_lower_level_pairs` to make the aggregation auditable. Smaller distances always mean greater acoustic similarity.

## 4. Descriptive embeddings

Fit each embedding once, pooled across stage and context. Sequence embeddings
contain the 107 eligible sequence repertoires. Call embeddings contain all
3,612 unique analysed calls, matching the PCA/VAE/DTW representation set and
the locked supplementary figures; only subsequent repertoire aggregation and
modelling filter to the 3,527 calls belonging to eligible call repertoires. Use
PaCMAP for transition probability, bigram, STP, MFCC, and VAE feature vectors.
Use UMAP on phee-repeat vectors, and use UMAP with `metric="precomputed"` for
local-alignment and DTW matrices. Algorithm hyperparameters, package versions,
seed, ordered identifiers, input and coordinate shapes, and the source-artifact
checksum are stored in the manifest. Coordinate CSVs use the fixed columns
`embedding_x`, `embedding_y`, the observation ID, and canonical
`individual_id`, `pair_id`, `stage`, and `context` values in source order.

Pooled and faceted panels must use the same cached coordinates and axis limits. Embeddings are descriptive and are never used as model responses or as substitutes for original-space distance matrices.

## 5. Model tables and standardization

Write one long-form row per repertoire comparison and metric. Standardize raw distance independently for each metric, using the complete eligible set for that structural level jointly across all pairs, stages, and contexts:

\[
z_{im}=\frac{d_{im}-\bar d_m}{s_m}.
\]

Here \(s_m\) is the population standard deviation over all eligible rows for
metric \(m\) (`ddof = 0`), not the sample standard deviation. Record the
standardization convention, means, and scale values in the derived-data
manifest. Stage is coded `before = -0.5` and `after = +0.5`.

Expected table sizes are 331 × 4 = 1,324 rows for sequence structure and 431 × 4 = 1,724 rows for call structure.

## 6. Bayesian multilevel models

Fit one unified model for sequence structure and one for call structure. For observation \(i\), the Gaussian identity-link model is represented by the brms-style formula

```r
standardized_distance ~ stage * context * metric +
  (1 + stage + context + stage:context || pair_id) +
  (1 | mm(repertoire_a_id, repertoire_b_id))
```

The full fixed-effect interaction permits stage and context effects to differ by metric. The uncorrelated pair term supplies a varying intercept and varying slopes for stage, context, and their interaction. The equal-weight multi-membership intercept accounts for the repeated use of both session repertoires across comparisons. Repertoire identifiers must be globally unique.

Use these priors:

| Parameter class | Prior |
|---|---|
| Intercept and population-level coefficients | Normal(0, 1) |
| Group-level standard deviations | Exponential(1) |
| Residual standard deviation | Exponential(1) |

Run four chains with 4,000 iterations per chain, including 2,000 warm-up iterations. Record the sampler seed and control settings. Release caches require trace plots and posterior predictive checks, maximum \(\hat R \le 1.01\), minimum bulk and tail effective sample sizes of 1,000, no divergent transitions, and no maximum-tree-depth exceedances. A refit that fails this contract is not a reportable release cache.

## 7. Posterior estimands

Obtain manuscript estimates from posterior expected standardized distances generated by the unified model. For each draw, pair, metric, and context:

\[
\Delta_{pmc}=E(d\mid\mathrm{after},p,m,c)-E(d\mid\mathrm{before},p,m,c).
\]

Average metric-specific \(\Delta\) values equally over the three observed focal pairs. Obtain the combined structural-level estimate by then averaging equally over the four standardized metrics. Every observed repertoire comparison still contributes separately to fitting; equal weighting applies to the reported posterior contrast, not to row replication.

The audience contrast is

\[
\Psi=\Delta_\mathrm{Partner}-\Delta_\mathrm{Non\mbox{-}partner}.
\]

Positive \(\Psi\) indicates that distance decreased more in the Non-partner context. The primary estimate averages over all three focal pairs and necessarily uses partial pooling for an unobserved cell when a pair lacks Partner-after data. The common-support sensitivity estimate repeats the posterior marginalization over only the two pairs observed in all four stage-by-context cells. Report posterior medians and equal-tail 95% credible intervals.

## 8. Cache contract

All expensive computations are opt-in. Cached mode must:

1. require every artifact needed for the requested result;
2. validate schema version, source hashes, artifact checksum, shape, ordered identifiers, and parameters;
3. stop on missing or incompatible input;
4. never calculate an approximate replacement;
5. never train or fit a model as a fallback.

An explicit full rebuild creates or replaces an artifact only after validating complete source inputs. Every cache record stores a portable path, byte size, SHA-256, method parameters, and—where order matters—an ordered-ID hash and count. Method-specific exporters also record their applicable package versions, seeds, shapes, and source fingerprints. Posterior sidecars bind model input content, formula, sampling settings, software, and draw checksum. Downstream artifacts must be invalidated whenever an upstream ordered index or parameter changes.
