# Historical processed inputs

These wide CSV files are retained solely to make the inexpensive migration to
the canonical schema auditable. They preserve historical column names,
duplicate call rows, and derived PCA columns. Author-machine audio paths have
been replaced with portable `extracted calls/<filename>` values; this metadata
change does not alter any analytical value.

Run `python scripts/import_legacy_processed_data.py` to deduplicate calls,
normalize audio paths as repository-relative paths, construct collision-safe
receiver-specific repertoire identifiers, and write `data/processed/`.

Downstream analysis reads the canonical processed tables. The sole documented
exception is `scripts/prepare_cached_artifacts.py`, which extracts and validates
the frozen historical `trad_PC*` and `mfcc_PC*` scores without refitting PCA.

The migration removes 292 repeated call rows by normalized WAV filename. It
also repairs one documented legacy string-concatenation sentinel (`nanAA` to
`AA`) before validating the six-symbol sequence alphabet. Both transformations
are deterministic and recorded in the processed-data manifest.
