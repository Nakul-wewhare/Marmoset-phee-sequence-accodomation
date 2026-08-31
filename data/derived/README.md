# Derived analysis tables

This directory contains session-pair repertoire distances and the two
long-form model tables. Each comparison occurs once per metric. Raw distances
are standardized within metric across every pair, stage, and context before the
sequence and call tables are passed to R.

The complete expected sizes are 1,324 sequence rows and 1,724 call rows. A
partial table is never treated as a valid cached result.

`scripts/assemble_call_model_table.py` creates and registers the two call CSVs
only when every canonical call input and precomputed DTW/VAE cache passes its
manifest checksum, parameter, schema, and ordered-ID checks.
