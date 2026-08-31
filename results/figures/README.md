# Figure exports

The PNG files listed in `figure_provenance.csv` are the inline figure exports
from the identified current manuscript document revision. The provenance table
records the document and revision identifiers, inline-object identifier,
retrieval date, and SHA-256 digest for every export.

These are versioned presentation artifacts. They connect notebook sections to
the manuscript and provide visual comparison targets, but they are not
substitutes for distance matrices, embeddings, posterior draws, or model
diagnostics. Regeneration code must use the corresponding analysis cache.

The two SVG files are editable sources retained from the project workspace for
Figure 2 and Figure S1. The provenance-listed PNGs are the authoritative
manuscript snapshots for this refresh.

Files ending in `_regenerated.png` are disposable reader-facing analytical
renders. The cached plotting script creates them from manifest-validated
tables and coordinates, never overwrites an unsuffixed snapshot, and does not
register them as quantitative inputs. A future archival release may register a
regenerated presentation artifact explicitly; mere file presence is not
provenance.
