#!/usr/bin/env python3
"""Prepare the inexpensive, versioned artifacts used by the public notebooks.

The default action builds session inventories, comparison indexes, sequence
feature representations, and PCA-score tables by reusing the scores already
stored in the historical processed call table.  It does not read audio,
calculate DTW, train a VAE, fit an embedding, or sample a Bayesian model.

Local alignment is deterministic and much smaller than the all-call acoustic
steps, but it is still opt-in.  Pass ``--recompute-sequence-alignment`` to
calculate it and complete the four-metric sequence model table.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from marmoset_convergence import __version__
from marmoset_convergence.calls import call_basename
from marmoset_convergence.matrices import (
    CondensedDistanceMatrix,
    SequenceRepresentationCache,
    feature_distance_condensed,
    load_sequence_distance_cache,
    save_sequence_distance_cache,
    save_sequence_representation_cache,
    square_to_condensed,
)
from marmoset_convergence.model_tables import (
    SEQUENCE_METRICS,
    attach_repertoire_distances,
    build_model_table,
)
from marmoset_convergence.paths import ProjectPaths
from marmoset_convergence.provenance import (
    CURRENT_RELEASE_ARTIFACT_PATHS,
    LOCAL_ALIGNMENT_PARAMETERS,
    MODEL_TABLE_PARAMETERS,
    SEQUENCE_REPRESENTATION_PARAMETERS,
    CacheManifest,
    ordered_ids_sha256,
    sha256_file,
)
from marmoset_convergence.sequence import (
    build_sequence_representations,
    local_alignment_distance_matrix,
)
from marmoset_convergence.sessions import (
    build_session_inventory,
    build_session_pairs,
)


BONDED_PAIRS = {
    "A": ("Tabor", "Lola"),
    "B": ("Wuschel", "Olympia"),
    "C": ("Odin", "Nougatti"),
}
SEQUENCE_ALPHABET = tuple("AKPRST")
EXPECTED = {
    "sequence_inventory_rows": 107,
    "eligible_sequences": 1496,
    "sequence_comparisons": 331,
    "call_inventory_rows": 130,
    "eligible_calls": 3527,
    "call_comparisons": 431,
}


def _write_csv(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False, lineterminator="\n")


def _require_equal(observed: int, expected: int, label: str) -> None:
    if observed != expected:
        raise ValueError(f"Expected {expected:,} {label}, observed {observed:,}")


def _extract_precomputed_pca_scores(
    source: Path,
    calls: pd.DataFrame,
    *,
    source_prefix: str,
    output_prefix: str,
) -> pd.DataFrame:
    """Extract legacy PC scores in canonical call order without refitting PCA."""

    source_columns = [f"{source_prefix}_PC{index}" for index in range(1, 6)]
    raw = pd.read_csv(source, usecols=["filename", *source_columns])
    raw["call_id"] = raw["filename"].map(call_basename)

    conflicts = []
    for call_id, group in raw.groupby("call_id", sort=False):
        numeric = group[source_columns].apply(pd.to_numeric, errors="coerce")
        if numeric.isna().any().any():
            conflicts.append(f"{call_id}: non-numeric score")
        elif len(group) > 1 and not np.allclose(
            numeric.to_numpy(),
            np.repeat(numeric.iloc[[0]].to_numpy(), len(group), axis=0),
            rtol=0.0,
            atol=0.0,
        ):
            conflicts.append(f"{call_id}: duplicate rows disagree")
    if conflicts:
        raise ValueError("Invalid precomputed PCA scores; " + "; ".join(conflicts[:5]))

    scores = raw.drop_duplicates("call_id", keep="first").set_index("call_id")
    expected_ids = calls["call_id"].astype(str)
    missing = sorted(set(expected_ids) - set(scores.index.astype(str)))
    extra = sorted(set(scores.index.astype(str)) - set(expected_ids))
    if missing or extra:
        raise ValueError(
            f"Precomputed {source_prefix} IDs disagree with canonical calls "
            f"(missing={len(missing)}, extra={len(extra)})"
        )
    scores = scores.loc[expected_ids, source_columns].reset_index()
    scores = scores.rename(
        columns={
            source_column: f"{output_prefix}_pc_{index}"
            for index, source_column in enumerate(source_columns, start=1)
        }
    )
    values = scores.drop(columns="call_id").to_numpy(dtype=float)
    if not np.isfinite(values).all() or (np.std(values, axis=0, ddof=0) == 0).any():
        raise ValueError(f"Precomputed {source_prefix} scores are invalid")
    return scores


def _sequence_repertoires(
    sequences: pd.DataFrame, inventory: pd.DataFrame
) -> dict[str, tuple[tuple[str, ...], ...]]:
    eligible = set(inventory["repertoire_id"].astype(str))
    filtered = sequences.loc[sequences["repertoire_id"].astype(str).isin(eligible)]
    grouped = filtered.groupby("repertoire_id", sort=False)["sequence"].apply(
        lambda values: tuple(tuple(str(value)) for value in values)
    )
    ordered_ids = tuple(inventory["repertoire_id"].astype(str))
    missing = [repertoire_id for repertoire_id in ordered_ids if repertoire_id not in grouped]
    if missing:
        raise ValueError(f"Eligible repertoires lack sequences: {missing[:5]}")
    return {repertoire_id: grouped[repertoire_id] for repertoire_id in ordered_ids}


def _manifest(
    paths: ProjectPaths,
    *,
    sequence_inventory: pd.DataFrame,
    call_inventory: pd.DataFrame,
    sequence_pairs: pd.DataFrame,
    call_pairs: pd.DataFrame,
    alignment_parameters: dict[str, object] | None,
) -> CacheManifest:
    # Exact DTW, VAE latents, and the derived call tables are owned by their
    # release/assembly workflows. They are intentionally absent here, so valid
    # existing records for those paths are preserved by the loop below.
    candidates = {
        "processed_manifest": paths.processed / "manifest.json",
        "processed_sequences": paths.sequences,
        "processed_calls": paths.calls,
        "sequence_session_inventory": paths.sequence_inventory,
        "call_session_inventory": paths.call_inventory,
        "sequence_session_pairs": paths.cache / "sequence_session_pairs.csv",
        "call_session_pairs": paths.cache / "call_session_pairs.csv",
        "sequence_representations": paths.sequence_representations,
        "stp_pca_scores": paths.stp_pca_scores,
        "mfcc_pca_scores": paths.mfcc_pca_scores,
        "sequence_distances": paths.sequence_distances,
        "sequence_repertoire_distances": paths.sequence_repertoire_distances,
        "model_sequence": paths.model_sequence,
    }
    artifacts = {name: path for name, path in candidates.items() if path.is_file()}
    owned_paths = {
        path.resolve().relative_to(paths.root).as_posix() for path in candidates.values()
    }
    preserved_artifacts = {}
    preserved_metadata: dict[str, object] = {}
    if paths.manifest.is_file():
        existing_manifest = CacheManifest.load(paths.manifest)
        preserved_metadata.update(existing_manifest.metadata)
        for name, record in existing_manifest.artifacts.items():
            if record.path in owned_paths:
                continue
            existing_manifest.validate(paths.root, required_artifacts=[name])
            preserved_artifacts[name] = record

    parameters: dict[str, dict[str, object]] = {
        "processed_manifest": {"schema": "processed-input-manifest-1.0"},
        "processed_sequences": {"rows": 1619, "schema": "canonical-sequences-1.0"},
        "processed_calls": {"rows": 3612, "schema": "canonical-calls-1.0"},
        "sequence_session_inventory": {
            "minimum_sequences": 5,
            "rows": EXPECTED["sequence_inventory_rows"],
            "eligible_units": EXPECTED["eligible_sequences"],
        },
        "call_session_inventory": {
            "minimum_calls": 5,
            "rows": EXPECTED["call_inventory_rows"],
            "eligible_units": EXPECTED["eligible_calls"],
        },
        "sequence_session_pairs": {
            "comparison_rule": "all_cross_individual_sessions_within_pair_stage_context",
            "rows": EXPECTED["sequence_comparisons"],
        },
        "call_session_pairs": {
            "comparison_rule": "all_cross_individual_sessions_within_pair_stage_context",
            "rows": EXPECTED["call_comparisons"],
        },
        "sequence_representations": {
            **SEQUENCE_REPRESENTATION_PARAMETERS,
            "shapes": {
                "transition_probability": [107, 36],
                "bigram": [107, 36],
                "phee_repeat": [107, 4],
            },
        },
        "stp_pca_scores": {
            "components": 5,
            "source": "legacy precomputed trad_PC1..trad_PC5",
            "pca_refit": False,
            "duplicate_rule": "one canonical row per call_id",
            "basis_fit_rows": 3907,
            "exported_rows_before_deduplication": 3904,
            "canonical_unique_calls": 3612,
            "basis_fit_stage": "before three unavailable-audio rows were omitted and before filename deduplication",
            "standardization": "sklearn StandardScaler on the 3,907-row historical merged table",
            "pca": "sklearn PCA(n_components=5, random_state=42)",
            "explained_variance_percent": [45.97, 33.45, 12.16, 4.06, 2.24],
        },
        "mfcc_pca_scores": {
            "components": 5,
            "source": "legacy precomputed mfcc_PC1..mfcc_PC5",
            "pca_refit": False,
            "duplicate_rule": "one canonical row per call_id",
            "basis_fit_rows": 3907,
            "exported_rows_before_deduplication": 3904,
            "canonical_unique_calls": 3612,
            "basis_fit_stage": "before three unavailable-audio rows were omitted and before filename deduplication",
            "standardization": "sklearn StandardScaler on the 3,907-row historical merged table",
            "pca": "sklearn PCA(n_components=5, random_state=42)",
            "explained_variance_percent": [34.19, 21.91, 11.32, 6.29, 5.41],
        },
        "sequence_repertoire_distances": {
            "metrics": list(SEQUENCE_METRICS),
            "rows": EXPECTED["sequence_comparisons"] * len(SEQUENCE_METRICS),
        },
        "model_sequence": {
            **MODEL_TABLE_PARAMETERS,
            "metrics": list(SEQUENCE_METRICS),
            "rows": EXPECTED["sequence_comparisons"] * len(SEQUENCE_METRICS),
        },
    }
    if alignment_parameters is not None:
        parameters["sequence_distances"] = alignment_parameters
    if paths.sequence_repertoire_distances.is_file():
        distances = pd.read_csv(paths.sequence_repertoire_distances)
        standardization = {}
        for metric, values in distances.groupby("metric", sort=False)["distance"]:
            numeric = values.to_numpy(dtype=float)
            standardization[str(metric)] = {
                "mean": float(np.mean(numeric)),
                "standard_deviation": float(np.std(numeric, ddof=0)),
                "n": int(len(numeric)),
            }
        parameters["model_sequence"]["standardization_by_metric"] = standardization

    ordered_ids: dict[str, Iterable[str]] = {
        "processed_sequences": tuple(pd.read_csv(paths.sequences, usecols=["sequence_id"])["sequence_id"].astype(str)),
        "processed_calls": tuple(pd.read_csv(paths.calls, usecols=["call_id"])["call_id"].astype(str)),
        "sequence_session_inventory": tuple(sequence_inventory["repertoire_id"].astype(str)),
        "call_session_inventory": tuple(call_inventory["repertoire_id"].astype(str)),
        "sequence_session_pairs": tuple(sequence_pairs["comparison_id"].astype(str)),
        "call_session_pairs": tuple(call_pairs["comparison_id"].astype(str)),
        "sequence_representations": tuple(sequence_inventory["repertoire_id"].astype(str)),
        "stp_pca_scores": tuple(pd.read_csv(paths.stp_pca_scores, usecols=["call_id"])["call_id"].astype(str)),
        "mfcc_pca_scores": tuple(pd.read_csv(paths.mfcc_pca_scores, usecols=["call_id"])["call_id"].astype(str)),
    }
    if "sequence_distances" in artifacts:
        ordered_ids["sequence_distances"] = tuple(
            sequence_inventory["repertoire_id"].astype(str)
        )

    missing_current_release = [
        relative
        for relative in CURRENT_RELEASE_ARTIFACT_PATHS
        if not (paths.root / relative).is_file()
    ]
    generated_manifest = CacheManifest.from_files(
        paths.root,
        artifacts,
        parameters=parameters,
        ordered_ids=ordered_ids,
        metadata={
            **preserved_metadata,
            "workflow": "cache-first; expensive components require explicit flags",
            "producer": "scripts/prepare_cached_artifacts.py",
            "software": {
                "marmoset_convergence": __version__,
                "numpy": np.__version__,
                "pandas": pd.__version__,
            },
            "release_completeness": (
                "complete" if not missing_current_release else "partial"
            ),
            "missing_current_release_artifacts": missing_current_release,
            "prohibited_substitutions": [
                "legacy 2,200-call FastDTW matrix",
                "legacy context-separated three-metric models",
                "synthetic posterior draws",
            ],
        },
    )
    merged_artifacts = {**preserved_artifacts, **generated_manifest.artifacts}
    return CacheManifest(
        artifacts=merged_artifacts,
        metadata=generated_manifest.metadata,
    )


def prepare(root: Path, *, recompute_sequence_alignment: bool) -> dict[str, int | float | str]:
    paths = ProjectPaths.from_root(root)
    paths.create_output_directories()

    sequences = pd.read_csv(paths.sequences)
    calls = pd.read_csv(paths.calls)
    _require_equal(len(sequences), 1619, "processed sequences")
    _require_equal(len(calls), 3612, "unique processed calls")

    sequence_inventory = build_session_inventory(
        sequences,
        unit_id_col="sequence_id",
        count_col="n_sequences",
        minimum_units=5,
        bonded_pairs=BONDED_PAIRS,
    )
    call_inventory = build_session_inventory(
        calls,
        unit_id_col="call_id",
        count_col="n_calls",
        minimum_units=5,
        bonded_pairs=BONDED_PAIRS,
    )
    sequence_pairs = build_session_pairs(sequence_inventory, BONDED_PAIRS)
    call_pairs = build_session_pairs(call_inventory, BONDED_PAIRS)

    _require_equal(len(sequence_inventory), EXPECTED["sequence_inventory_rows"], "eligible sequence repertoires")
    _require_equal(int(sequence_inventory["n_sequences"].sum()), EXPECTED["eligible_sequences"], "eligible sequences")
    _require_equal(len(sequence_pairs), EXPECTED["sequence_comparisons"], "sequence comparisons")
    _require_equal(len(call_inventory), EXPECTED["call_inventory_rows"], "eligible call repertoires")
    _require_equal(int(call_inventory["n_calls"].sum()), EXPECTED["eligible_calls"], "eligible calls")
    _require_equal(len(call_pairs), EXPECTED["call_comparisons"], "call comparisons")

    _write_csv(sequence_inventory, paths.sequence_inventory)
    _write_csv(call_inventory, paths.call_inventory)
    _write_csv(sequence_pairs, paths.cache / "sequence_session_pairs.csv")
    _write_csv(call_pairs, paths.cache / "call_session_pairs.csv")

    legacy_calls = paths.root / "data/source/legacy_processed_calls.csv"
    stp_scores = _extract_precomputed_pca_scores(
        legacy_calls, calls, source_prefix="trad", output_prefix="stp"
    )
    mfcc_scores = _extract_precomputed_pca_scores(
        legacy_calls, calls, source_prefix="mfcc", output_prefix="mfcc"
    )
    _write_csv(stp_scores, paths.stp_pca_scores)
    _write_csv(mfcc_scores, paths.mfcc_pca_scores)

    repertoires = _sequence_repertoires(sequences, sequence_inventory)
    representations = build_sequence_representations(
        repertoires,
        alphabet=SEQUENCE_ALPHABET,
        phee_token="A",
        repeat_minimum=2,
        repeat_maximum=5,
    )
    representation_cache = SequenceRepresentationCache(
        representations.repertoire_ids,
        representations.transition_probability,
        representations.bigram,
        representations.phee_repeat,
    )
    save_sequence_representation_cache(paths.sequence_representations, representation_cache)

    alignment_parameters: dict[str, object] | None = None
    if recompute_sequence_alignment:
        alignment = local_alignment_distance_matrix(repertoires, allow_expensive=True)
        save_sequence_distance_cache(
            paths.sequence_distances,
            alignment.repertoire_ids,
            square_to_condensed(alignment.distances),
        )
        sequence_counts = tuple(
            f"{repertoire_id}:{len(repertoires[repertoire_id])}"
            for repertoire_id in alignment.repertoire_ids
        )
        alignment_parameters = {
            **LOCAL_ALIGNMENT_PARAMETERS,
            "global_max_similarity": alignment.global_max_similarity,
            "eligible_set_version": ordered_ids_sha256(alignment.repertoire_ids),
            "eligible_repertoires": len(alignment.repertoire_ids),
            "eligible_sequences": sum(len(value) for value in repertoires.values()),
            "sequence_counts_sha256": ordered_ids_sha256(sequence_counts),
            "processed_sequences_sha256": sha256_file(paths.sequences),
            "implementation": "marmoset_convergence.sequence.local_alignment_distance_matrix",
            "implementation_version": __version__,
            "implementation_sha256": sha256_file(
                paths.root / "src/marmoset_convergence/sequence.py"
            ),
        }
    elif paths.sequence_distances.is_file():
        existing_manifest = CacheManifest.load(paths.manifest)
        matching_records = [
            (name, record)
            for name, record in existing_manifest.artifacts.items()
            if record.path == paths.sequence_distances.relative_to(paths.root).as_posix()
        ]
        if len(matching_records) != 1:
            raise ValueError(
                "The existing sequence alignment cache must have exactly one record "
                "in data/cache/manifest.json"
            )
        alignment_name, alignment_record = matching_records[0]
        existing_manifest.validate(
            paths.root,
            required_artifacts=[alignment_name],
            expected_parameters={alignment_name: LOCAL_ALIGNMENT_PARAMETERS},
        )
        required_dynamic = {
            "global_max_similarity",
            "eligible_set_version",
            "eligible_repertoires",
            "eligible_sequences",
            "sequence_counts_sha256",
            "processed_sequences_sha256",
            "implementation",
            "implementation_version",
            "implementation_sha256",
        }
        missing_dynamic = sorted(required_dynamic.difference(alignment_record.parameters))
        if missing_dynamic:
            raise ValueError(
                "Existing sequence alignment cache lacks current provenance: "
                + ", ".join(missing_dynamic)
            )
        alignment_parameters = dict(alignment_record.parameters)

    if paths.sequence_distances.is_file():
        alignment_matrix = load_sequence_distance_cache(paths.sequence_distances)
        if alignment_matrix.ids != representations.repertoire_ids:
            raise ValueError("Sequence alignment and representation ID orders disagree")
        matrices = {
            "transition_probability": CondensedDistanceMatrix(
                representations.repertoire_ids,
                feature_distance_condensed(representations.transition_probability),
                "transition_probability",
            ),
            "bigram": CondensedDistanceMatrix(
                representations.repertoire_ids,
                feature_distance_condensed(representations.bigram),
                "bigram",
            ),
            "phee_repeat": CondensedDistanceMatrix(
                representations.repertoire_ids,
                feature_distance_condensed(representations.phee_repeat),
                "phee_repeat",
            ),
            "local_alignment": alignment_matrix,
        }
        sequence_distances = attach_repertoire_distances(sequence_pairs, matrices)
        _require_equal(len(sequence_distances), 1324, "sequence model-distance rows")
        _write_csv(sequence_distances, paths.sequence_repertoire_distances)
        model_sequence = build_model_table(
            sequence_distances, expected_metrics=SEQUENCE_METRICS, ddof=0
        )
        _require_equal(len(model_sequence), 1324, "sequence model rows")
        _write_csv(model_sequence, paths.model_sequence)

    manifest = _manifest(
        paths,
        sequence_inventory=sequence_inventory,
        call_inventory=call_inventory,
        sequence_pairs=sequence_pairs,
        call_pairs=call_pairs,
        alignment_parameters=alignment_parameters,
    )
    manifest.write(paths.manifest)
    return {
        "status": "ok",
        "sequence_repertoires": len(sequence_inventory),
        "sequence_comparisons": len(sequence_pairs),
        "call_repertoires": len(call_inventory),
        "call_comparisons": len(call_pairs),
        "local_alignment": "recomputed" if recompute_sequence_alignment else "not requested",
        "registered_artifacts": len(manifest.artifacts),
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--recompute-sequence-alignment",
        action="store_true",
        help="Explicitly calculate the complete eligible-set local-alignment matrix.",
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    result = prepare(
        args.root.resolve(),
        recompute_sequence_alignment=args.recompute_sequence_alignment,
    )
    for key, value in result.items():
        print(f"{key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
