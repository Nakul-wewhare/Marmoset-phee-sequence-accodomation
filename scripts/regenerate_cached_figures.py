#!/usr/bin/env python3
"""Regenerate analytical PNGs from validated, already-computed caches.

The script performs a complete read-only preflight before writing any figure.
It never reads audio, computes distances, or fits an embedding, and every output
is forced to a distinct ``*_regenerated.png`` path beside the locked snapshot.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd

from marmoset_convergence.cli import validate_cache
from marmoset_convergence.figures import (
    EMBEDDING_FIGURES,
    prepare_embedding_table,
    regenerate_embedding_space,
    regenerate_sequence_summary,
    regenerated_figure_path,
)
from marmoset_convergence.paths import ProjectPaths
from marmoset_convergence.provenance import CacheManifest

FIGURE_2_LOCKED_PATH = "results/figures/main/figure_2.png"
FIGURE_CODES = ("2", "4", "s2", "s3", "s4", "s5", "s6", "s7", "s8")


def _artifact_name_for_path(manifest: CacheManifest, relative_path: str) -> str:
    matches = [
        name
        for name, record in manifest.artifacts.items()
        if record.path == Path(relative_path).as_posix()
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one manifest record for {relative_path!r}; "
            f"found {len(matches)}"
        )
    return matches[0]


def _required_paths(codes: Iterable[str]) -> tuple[str, ...]:
    paths: list[str] = []
    for code in codes:
        if code == "2":
            paths.append("data/processed/sequences.csv")
            continue
        spec = EMBEDDING_FIGURES[code]
        metadata_path = (
            "data/cache/sequence_session_inventory.csv"
            if spec.level == "sequence"
            else "data/processed/calls.csv"
        )
        paths.extend((metadata_path, f"data/cache/embeddings/{spec.metric}.csv"))
    return tuple(dict.fromkeys(paths))


def regenerate_cached_figures(
    root: Path,
    codes: Iterable[str],
    *,
    overwrite_regenerated: bool = False,
) -> tuple[Path, ...]:
    """Validate all requested inputs, then write only regenerated exports."""

    paths = ProjectPaths.from_root(root)
    selected = tuple(dict.fromkeys(code.casefold() for code in codes))
    unknown = sorted(set(selected).difference(FIGURE_CODES))
    if unknown:
        raise ValueError("Unknown figure codes: " + ", ".join(unknown))
    if not selected:
        raise ValueError("Select at least one figure")

    locked_targets = tuple(
        paths.root
        / (
            FIGURE_2_LOCKED_PATH
            if code == "2"
            else EMBEDDING_FIGURES[code].locked_relative_path
        )
        for code in selected
    )
    destinations = tuple(regenerated_figure_path(target) for target in locked_targets)
    if len(set(destinations)) != len(destinations):
        raise ValueError("Multiple requested figures resolve to the same regenerated target")
    for locked, destination in zip(locked_targets, destinations):
        if locked.exists() and destination.exists() and locked.resolve() == destination.resolve():
            raise ValueError(f"Regenerated target resolves to a locked snapshot: {destination}")
        if destination.exists() and not overwrite_regenerated:
            raise FileExistsError(
                f"Regenerated export already exists: {destination}. Pass "
                "--overwrite-regenerated to replace only derived outputs."
            )

    manifest = CacheManifest.load(paths.manifest)
    required_paths = _required_paths(selected)
    artifact_names = tuple(
        _artifact_name_for_path(manifest, relative_path)
        for relative_path in required_paths
    )
    # Explicit artifact names keep this cheap figure command independent from
    # the completeness profile used by the main cache-validation CLI.
    validate_cache(paths.root, artifact_names)

    loaded: dict[str, pd.DataFrame] = {
        relative_path: pd.read_csv(paths.root / relative_path)
        for relative_path in required_paths
    }
    # Complete all schema/order/grouping checks before the first PNG write.
    for code in selected:
        if code == "2":
            continue
        spec = EMBEDDING_FIGURES[code]
        metadata_path = (
            "data/cache/sequence_session_inventory.csv"
            if spec.level == "sequence"
            else "data/processed/calls.csv"
        )
        prepare_embedding_table(
            loaded[f"data/cache/embeddings/{spec.metric}.csv"],
            loaded[metadata_path],
            id_column=spec.id_column,
        )

    outputs: list[Path] = []
    for code in selected:
        if code == "2":
            outputs.append(
                regenerate_sequence_summary(
                    loaded["data/processed/sequences.csv"],
                    paths.root / FIGURE_2_LOCKED_PATH,
                    overwrite_regenerated=overwrite_regenerated,
                )
            )
            continue
        spec = EMBEDDING_FIGURES[code]
        metadata_path = (
            "data/cache/sequence_session_inventory.csv"
            if spec.level == "sequence"
            else "data/processed/calls.csv"
        )
        outputs.append(
            regenerate_embedding_space(
                loaded[f"data/cache/embeddings/{spec.metric}.csv"],
                loaded[metadata_path],
                spec,
                paths.root / spec.locked_relative_path,
                overwrite_regenerated=overwrite_regenerated,
            )
        )
    return tuple(outputs)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--figure",
        action="append",
        type=str.casefold,
        choices=FIGURE_CODES,
        dest="figures",
        help="Figure code to regenerate (repeatable); default: all Python figures",
    )
    parser.add_argument(
        "--overwrite-regenerated",
        action="store_true",
        help="Replace existing *_regenerated.png files; locked PNGs remain protected",
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    outputs = regenerate_cached_figures(
        args.root,
        args.figures or FIGURE_CODES,
        overwrite_regenerated=args.overwrite_regenerated,
    )
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
