# Versioned analysis caches

The normal workflow is cache-first. It validates checksums, shapes, ordered
identifiers, and method parameters before reading an artifact, and it never
recomputes a missing cache automatically.

The revised analysis requires exact constrained DTW distances for 3,612 calls,
32-dimensional Chatter VAE latent means, pooled embedding coordinates, and
posterior draws from the unified interaction models. The older 2,200-call
FastDTW subset and the older context-separated three-metric models are not
compatible and must not be placed here.

Release assets that are not committed directly to Git should be downloaded to
the canonical paths listed in `docs/reproducibility.md` and registered in
`manifest.json`. Default validation checks every registered artifact that is
present and reports the incomplete release. Strict validation stops and names
every missing or unregistered revised-analysis artifact.

With the exact-DTW and VAE files installed, run
`python scripts/assemble_call_model_table.py --root .` from the repository root
to validate all four call representations and build the derived call tables.
The command is cache-only; it performs no audio, DTW, PCA, or VAE recomputation.
