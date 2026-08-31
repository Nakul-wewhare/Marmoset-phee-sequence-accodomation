# Archival release checklist

This checkout is a transparent partial reproducibility bundle. It includes the
canonical processed tables, inexpensive caches, tested analysis code, and
provenance-locked manuscript figure snapshots. The exact-DTW matrix, VAE latent
means, pooled embeddings, complete derived tables, posterior draws, and model
diagnostics must be installed and registered before describing the repository
as a complete archival release.

Before assigning an archive DOI:

1. Deposit or link the complete 3,612-call audio collection, with its access
   conditions and canonical identifier mapping.
2. Install the current exact-DTW, Chatter 0.1.6, embedding, derived-table,
   posterior-draw, and diagnostic artifacts at the paths in
   `docs/reproducibility.md`.
3. Register every artifact with portable paths, SHA-256 checksums, ordered-ID
   fingerprints, method parameters, and applicable software/seed metadata.
4. Run `bash scripts/run_repro_pipeline.sh --mode cached --strict` from a clean
   clone.
5. Execute the five numbered notebooks and confirm that regenerated figures
   are written only to `_regenerated` paths.
6. Archive the repository commit and resolved Python and R environments with
   the data bundle.
7. Add the verified archive and audio-access identifiers to `CITATION.cff`,
   `.zenodo.json`, the manuscript data-availability statement, and the
   repository release notes.

Until those checks pass, cite the source repository and associated manuscript
without claiming that the partial Git checkout contains every costly result.

Code is licensed under MIT. Data, audio, figures, and derived outputs are
licensed under CC BY 4.0 as described in the repository license files.
