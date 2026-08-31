#!/usr/bin/env python3
"""Replace author-machine audio paths in the archival call CSV.

The historical table is retained because it contains the frozen STP and MFCC
principal-component scores used by the reported analysis.  Absolute paths are
not analytical inputs, so this one-time release-preparation step rewrites only
``call_audio_path`` as ``extracted calls/<filename>`` while preserving row
order, duplicate rows, and every numerical value.
"""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path


def sanitize(path: Path) -> int:
    source = path.resolve()
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Historical call table has no header: {source}")
        required = {"filename", "call_audio_path"}
        missing = sorted(required.difference(reader.fieldnames))
        if missing:
            raise ValueError("Historical call table lacks: " + ", ".join(missing))

        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{source.name}.", suffix=".tmp", dir=source.parent
        )
        rows = 0
        try:
            with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as output:
                writer = csv.DictWriter(
                    output,
                    fieldnames=reader.fieldnames,
                    lineterminator="\n",
                )
                writer.writeheader()
                for row in reader:
                    filename = Path(str(row["filename"])).name
                    if not filename:
                        raise ValueError(f"Row {rows + 2} has an empty filename")
                    row["call_audio_path"] = f"extracted calls/{filename}"
                    writer.writerow(row)
                    rows += 1
            os.replace(temporary_name, source)
        except BaseException:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass
            raise
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "path",
        type=Path,
        nargs="?",
        default=Path("data/source/legacy_processed_calls.csv"),
    )
    args = parser.parse_args()
    rows = sanitize(args.path)
    print(f"Sanitized {rows:,} historical call rows in {args.path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
