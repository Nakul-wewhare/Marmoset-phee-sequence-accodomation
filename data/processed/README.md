# Canonical processed tables

`sequences.csv` contains 1,619 unique sequence observations and `calls.csv`
contains 3,612 unique phee calls. These are the canonical analytical tables
consumed by the rebuilt workflow. The inexpensive cache-preparation migration
also reads the frozen historical PCA-score columns from
`data/source/legacy_processed_calls.csv`; it validates and reorders those
scores against `calls.csv` and never refits the PCA basis. Schemas are documented in
`docs/data_dictionary.md`; their source and output checksums are in
`manifest.json`.

This layer is inexpensive to recreate. It contains no local absolute paths and
does not depend on the availability of the WAV files.
