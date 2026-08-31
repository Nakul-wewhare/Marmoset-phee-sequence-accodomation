# Analysis data

The public analysis layout separates historical inputs from canonical tables,
versioned caches, and derived model data:

- `source/` preserves the two historical wide CSV inputs.
- `processed/` contains portable, one-row-per-unit sequence and call tables.
- `cache/` is reserved for versioned representations and expensive results.
- `derived/` contains repertoire distances and model-ready long tables.

Every analysis artifact must be registered with a SHA-256 checksum in a
manifest. A file that is present but unregistered is not accepted by cached
mode. See `docs/data_dictionary.md` for schemas and
`docs/reproducibility.md` for the cached/full distinction.
