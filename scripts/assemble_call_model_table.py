#!/usr/bin/env python3
"""Assemble the call-distance and model tables from registered caches only.

This command reads no audio and never calculates DTW, trains a VAE, or refits
either PCA.  It requires the canonical tables and all four call representations
to be present in ``data/cache/manifest.json`` with valid checksums, parameters,
and ordered-call provenance before it writes either derived table.
"""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from marmoset_convergence.exceptions import (
    CacheValidationError,
    DataValidationError,
    MarmosetConvergenceError,
)
from marmoset_convergence.identifiers import normalize_stage
from marmoset_convergence.matrices import (
    FeatureDistanceMatrix,
    load_dtw_distance_cache,
)
from marmoset_convergence.model_tables import (
    CALL_METRICS,
    aggregate_item_distances,
    build_model_table,
)
from marmoset_convergence.paths import ProjectPaths
from marmoset_convergence.provenance import (
    CURRENT_RELEASE_ARTIFACT_PATHS,
    DTW_PARAMETERS,
    MODEL_TABLE_PARAMETERS,
    VAE_PARAMETERS,
    CacheManifest,
    ordered_ids_sha256,
)
from marmoset_convergence.sessions import normalize_context


INPUT_PATHS: Mapping[str, str] = {
    "processed_calls": "data/processed/calls.csv",
    "call_session_inventory": "data/cache/call_session_inventory.csv",
    "call_session_pairs": "data/cache/call_session_pairs.csv",
    "stp_pca_scores": "data/cache/stp_pca_scores.csv",
    "mfcc_pca_scores": "data/cache/mfcc_pca_scores.csv",
    "dtw_distances_condensed": "data/cache/dtw_distances_condensed.npz",
    "vae_latent_means": "data/cache/vae_latent_means.csv",
}

STATIC_PARAMETERS: Mapping[str, Mapping[str, object]] = {
    "processed_calls": {"schema": "canonical-calls-1.0"},
    "call_session_inventory": {"minimum_calls": 5},
    "call_session_pairs": {
        "comparison_rule": (
            "all_cross_individual_sessions_within_pair_stage_context"
        )
    },
    "stp_pca_scores": {"components": 5, "pca_refit": False},
    "mfcc_pca_scores": {"components": 5, "pca_refit": False},
    "dtw_distances_condensed": DTW_PARAMETERS,
    "vae_latent_means": VAE_PARAMETERS,
}

STP_COLUMNS = tuple(f"stp_pc_{index}" for index in range(1, 6))
MFCC_COLUMNS = tuple(f"mfcc_pc_{index}" for index in range(1, 6))


def _registered_inputs(root: Path, manifest: CacheManifest) -> dict[str, str]:
    """Resolve logical inputs to unique manifest records at canonical paths."""

    names: dict[str, str] = {}
    missing_files = []
    unregistered = []
    duplicate_records = []
    for logical_name, relative in INPUT_PATHS.items():
        if not (root / relative).is_file():
            missing_files.append(relative)
        matches = [
            name
            for name, record in manifest.artifacts.items()
            if record.path == relative
        ]
        if not matches:
            unregistered.append(relative)
        elif len(matches) > 1:
            duplicate_records.append(relative)
        else:
            names[logical_name] = matches[0]

    if missing_files or unregistered or duplicate_records:
        details = []
        if missing_files:
            details.append("missing files: " + ", ".join(missing_files))
        if unregistered:
            details.append(
                "not registered in data/cache/manifest.json: "
                + ", ".join(unregistered)
            )
        if duplicate_records:
            details.append(
                "registered more than once: " + ", ".join(duplicate_records)
            )
        raise CacheValidationError(
            "Cannot assemble the call model table from caches; "
            + "; ".join(details)
            + ". Install and register the exact DTW and VAE release caches; "
            "audio/DTW/VAE recomputation is never attempted by this command."
        )
    return names


def _validate_manifest_inputs(
    root: Path, manifest: CacheManifest, names: Mapping[str, str]
) -> None:
    expected_parameters = {
        names[logical_name]: parameters
        for logical_name, parameters in STATIC_PARAMETERS.items()
    }
    manifest.validate(
        root,
        required_artifacts=names.values(),
        expected_parameters=expected_parameters,
    )


def _read_csv(path: Path, label: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except (OSError, UnicodeError, pd.errors.ParserError) as exc:
        raise CacheValidationError(f"Could not read {label} cache {path}: {exc}") from exc


def _require_columns(
    table: pd.DataFrame, required: Iterable[str], label: str
) -> None:
    missing = [column for column in required if column not in table]
    if missing:
        raise DataValidationError(
            f"{label} is missing required columns: " + ", ".join(missing)
        )


def _string_ids(table: pd.DataFrame, column: str, label: str) -> tuple[str, ...]:
    _require_columns(table, [column], label)
    raw = table[column]
    if raw.isna().any() or raw.astype(str).str.strip().eq("").any():
        raise DataValidationError(f"{label} contains missing or empty {column} values")
    ids = tuple(raw.astype(str))
    if len(set(ids)) != len(ids):
        raise DataValidationError(f"{label} contains duplicate {column} values")
    return ids


def _feature_cache(
    table: pd.DataFrame,
    *,
    label: str,
    metric: str,
    distance_kind: str,
    feature_columns: Sequence[str] | None = None,
    dimensions: int | None = None,
) -> FeatureDistanceMatrix:
    ids = _string_ids(table, "call_id", label)
    if feature_columns is None:
        columns = tuple(column for column in table if column != "call_id")
        if dimensions is None:
            raise ValueError("dimensions is required for an inferred feature schema")
        if len(table.columns) != dimensions + 1 or len(columns) != dimensions:
            raise DataValidationError(
                f"{label} must contain call_id plus exactly {dimensions} coordinates"
            )
    else:
        columns = tuple(feature_columns)
        expected = ("call_id", *columns)
        if tuple(table.columns) != expected:
            raise DataValidationError(
                f"{label} columns must be exactly: " + ", ".join(expected)
            )

    numeric = table.loc[:, columns].apply(pd.to_numeric, errors="coerce")
    values = numeric.to_numpy(dtype=np.float64)
    if not np.isfinite(values).all():
        raise DataValidationError(
            f"{label} feature coordinates must be finite numeric values"
        )
    return FeatureDistanceMatrix(ids, values, metric, distance_kind)


def _require_same_ids(
    label: str, observed: Sequence[str], expected: Sequence[str]
) -> None:
    if tuple(observed) == tuple(expected):
        return
    observed_set = set(observed)
    expected_set = set(expected)
    missing = sorted(expected_set - observed_set)
    extra = sorted(observed_set - expected_set)
    if missing or extra:
        raise DataValidationError(
            f"{label} call IDs disagree with canonical calls "
            f"(missing={len(missing)}, extra={len(extra)})"
        )
    raise DataValidationError(f"{label} call-ID order disagrees with canonical calls")


def _require_record_ids(
    manifest: CacheManifest,
    record_name: str,
    ids: Sequence[str],
) -> None:
    record = manifest.artifacts[record_name]
    if record.n_ids is None or record.ordered_ids_sha256 is None:
        raise CacheValidationError(
            f"Artifact {record_name!r} does not record ordered-ID provenance"
        )
    if record.n_ids != len(ids):
        raise CacheValidationError(
            f"Artifact {record_name!r} records {record.n_ids} IDs; observed {len(ids)}"
        )
    if record.ordered_ids_sha256 != ordered_ids_sha256(ids):
        raise CacheValidationError(
            f"Artifact {record_name!r} ordered-ID checksum mismatch"
        )


def _require_record_count(
    manifest: CacheManifest,
    record_name: str,
    key: str,
    observed: int,
) -> None:
    recorded = manifest.artifacts[record_name].parameters.get(key)
    if recorded != observed:
        raise CacheValidationError(
            f"Artifact {record_name!r} parameter {key!r} is {recorded!r}; "
            f"observed {observed}"
        )


def _validate_inventory_and_membership(
    calls: pd.DataFrame,
    inventory: pd.DataFrame,
) -> pd.DataFrame:
    _require_columns(calls, ["call_id", "repertoire_id"], "processed calls")
    _require_columns(
        inventory,
        ["repertoire_id", "pair_id", "stage", "context", "n_calls"],
        "call session inventory",
    )
    call_ids = _string_ids(calls, "call_id", "processed calls")
    repertoire_ids = _string_ids(
        inventory, "repertoire_id", "call session inventory"
    )
    if calls["repertoire_id"].isna().any():
        raise DataValidationError("Processed calls contain missing repertoire_id values")

    counts = pd.to_numeric(inventory["n_calls"], errors="coerce")
    if (
        counts.isna().any()
        or not np.equal(counts, np.floor(counts)).all()
        or (counts < 5).any()
    ):
        raise DataValidationError(
            "Call inventory n_calls must contain integers greater than or equal to five"
        )

    inventory_ids = set(repertoire_ids)
    membership = calls.loc[
        calls["repertoire_id"].astype(str).isin(inventory_ids),
        ["call_id", "repertoire_id"],
    ].copy()
    membership["call_id"] = membership["call_id"].astype(str)
    membership["repertoire_id"] = membership["repertoire_id"].astype(str)
    observed_counts = (
        membership.groupby("repertoire_id", sort=False)["call_id"]
        .size()
        .reindex(repertoire_ids, fill_value=0)
        .to_numpy()
    )
    if not np.array_equal(observed_counts, counts.to_numpy(dtype=int)):
        mismatched = [
            repertoire_id
            for repertoire_id, observed, expected in zip(
                repertoire_ids, observed_counts, counts.to_numpy(dtype=int)
            )
            if observed != expected
        ]
        raise DataValidationError(
            "Call inventory counts disagree with canonical membership: "
            + ", ".join(mismatched[:5])
        )
    if len(membership) != int(counts.sum()):
        raise DataValidationError("Eligible-call total disagrees with call inventory")
    if len(call_ids) != len(calls):  # Defensive; _string_ids already establishes this.
        raise DataValidationError("Processed call IDs are invalid")
    return membership


def _validate_comparison_inventory(
    comparisons: pd.DataFrame,
    inventory: pd.DataFrame,
) -> None:
    required = (
        "comparison_id",
        "pair_id",
        "stage",
        "context",
        "repertoire_a_id",
        "repertoire_b_id",
    )
    _require_columns(comparisons, required, "call session pairs")
    _string_ids(comparisons, "comparison_id", "call session pairs")

    normalized_inventory = inventory.copy()
    normalized_inventory["repertoire_id"] = normalized_inventory[
        "repertoire_id"
    ].astype(str)
    normalized_inventory["stage"] = normalized_inventory["stage"].map(normalize_stage)
    normalized_inventory["context"] = normalized_inventory["context"].map(
        normalize_context
    )
    by_id = normalized_inventory.set_index("repertoire_id", drop=False)
    for comparison in comparisons.itertuples(index=False):
        stage = normalize_stage(comparison.stage)
        context = normalize_context(comparison.context)
        left_id = str(comparison.repertoire_a_id)
        right_id = str(comparison.repertoire_b_id)
        if left_id == right_id:
            raise DataValidationError(
                f"Comparison {comparison.comparison_id!r} repeats one repertoire"
            )
        for repertoire_id in (left_id, right_id):
            if repertoire_id not in by_id.index:
                raise DataValidationError(
                    f"Comparison {comparison.comparison_id!r} references ineligible "
                    f"repertoire {repertoire_id!r}"
                )
            repertoire = by_id.loc[repertoire_id]
            expected = (
                str(repertoire["pair_id"]),
                repertoire["stage"],
                repertoire["context"],
            )
            observed = (str(comparison.pair_id), stage, context)
            if observed != expected:
                raise DataValidationError(
                    f"Comparison {comparison.comparison_id!r} metadata disagrees with "
                    f"repertoire {repertoire_id!r}"
                )


def _atomic_csv(table: pd.DataFrame, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            prefix=f".{destination.name}.",
            suffix=".tmp",
            dir=destination.parent,
            delete=False,
        ) as handle:
            temporary_name = Path(handle.name)
            table.to_csv(handle, index=False, lineterminator="\n")
        temporary_name.replace(destination)
    finally:
        if temporary_name is not None and temporary_name.exists():
            temporary_name.unlink()


def _row_ids(table: pd.DataFrame) -> tuple[str, ...]:
    return tuple(
        f"{comparison_id}|metric={metric}"
        for comparison_id, metric in zip(table["comparison_id"], table["metric"])
    )


def _standardization_parameters(
    distances: pd.DataFrame,
) -> dict[str, dict[str, int | float]]:
    result = {}
    for metric, values in distances.groupby("metric", sort=False)["distance"]:
        numeric = values.to_numpy(dtype=np.float64)
        result[str(metric)] = {
            "mean": float(np.mean(numeric)),
            "standard_deviation": float(np.std(numeric, ddof=0)),
            "n": int(len(numeric)),
        }
    return result


def _register_outputs(
    paths: ProjectPaths,
    manifest: CacheManifest,
    input_names: Mapping[str, str],
    distances: pd.DataFrame,
    model: pd.DataFrame,
) -> CacheManifest:
    fingerprints = {
        logical_name: manifest.artifacts[record_name].sha256
        for logical_name, record_name in input_names.items()
    }
    parameters: dict[str, Mapping[str, object]] = {
        "call_repertoire_distances": {
            "aggregation": "mean_all_cross_repertoire_pairs",
            "distance_by_metric": {
                "stp": "euclidean",
                "mfcc": "euclidean",
                "dtw": "exact_precomputed_condensed_lookup",
                "vae": "cosine",
            },
            "metrics": list(CALL_METRICS),
            "comparisons": int(distances["comparison_id"].nunique()),
            "rows": int(len(distances)),
            "n_lower_level_pairs": "n_calls_a_times_n_calls_b",
            "input_artifact_sha256": fingerprints,
        },
        "model_call": {
            **MODEL_TABLE_PARAMETERS,
            "metrics": list(CALL_METRICS),
            "rows": int(len(model)),
            "source": "data/derived/call_repertoire_distances.csv",
            "standardization_by_metric": _standardization_parameters(distances),
            "input_artifact_sha256": fingerprints,
        },
    }
    generated = CacheManifest.from_files(
        paths.root,
        {
            "call_repertoire_distances": paths.call_repertoire_distances,
            "model_call": paths.model_call,
        },
        parameters=parameters,
        ordered_ids={
            "call_repertoire_distances": _row_ids(distances),
            "model_call": _row_ids(model),
        },
    )

    output_paths = {record.path for record in generated.artifacts.values()}
    retained = {
        name: record
        for name, record in manifest.artifacts.items()
        if record.path not in output_paths
    }
    missing_release = [
        relative
        for relative in CURRENT_RELEASE_ARTIFACT_PATHS
        if not (paths.root / relative).is_file()
    ]
    metadata = {
        **manifest.metadata,
        "release_completeness": "complete" if not missing_release else "partial",
        "missing_current_release_artifacts": missing_release,
        "derived_call_tables_producer": "scripts/assemble_call_model_table.py",
        "derived_call_tables_cache_only": True,
    }
    return CacheManifest(
        artifacts={**retained, **generated.artifacts},
        metadata=metadata,
    )


def assemble(root: Path) -> dict[str, int | str]:
    """Validate all caches, write both call tables, and register their checksums."""

    paths = ProjectPaths.from_root(root)
    manifest = CacheManifest.load(paths.manifest)
    input_names = _registered_inputs(paths.root, manifest)
    _validate_manifest_inputs(paths.root, manifest, input_names)

    calls = _read_csv(paths.calls, "processed calls")
    inventory = _read_csv(paths.call_inventory, "call session inventory")
    comparisons = _read_csv(
        paths.cache / "call_session_pairs.csv", "call session pairs"
    )
    stp_table = _read_csv(paths.stp_pca_scores, "STP PCA")
    mfcc_table = _read_csv(paths.mfcc_pca_scores, "MFCC PCA")
    vae_table = _read_csv(paths.vae_latent_means, "VAE latent means")

    canonical_call_ids = _string_ids(calls, "call_id", "processed calls")
    inventory_ids = _string_ids(
        inventory, "repertoire_id", "call session inventory"
    )
    comparison_ids = _string_ids(
        comparisons, "comparison_id", "call session pairs"
    )
    membership = _validate_inventory_and_membership(calls, inventory)
    _validate_comparison_inventory(comparisons, inventory)

    stp = _feature_cache(
        stp_table,
        label="STP PCA cache",
        metric="stp",
        distance_kind="euclidean",
        feature_columns=STP_COLUMNS,
    )
    mfcc = _feature_cache(
        mfcc_table,
        label="MFCC PCA cache",
        metric="mfcc",
        distance_kind="euclidean",
        feature_columns=MFCC_COLUMNS,
    )
    vae = _feature_cache(
        vae_table,
        label="VAE latent cache",
        metric="vae",
        distance_kind="cosine",
        dimensions=32,
    )
    dtw = load_dtw_distance_cache(paths.dtw_distances)

    for label, observed in (
        ("STP PCA cache", stp.ids),
        ("MFCC PCA cache", mfcc.ids),
        ("exact-DTW cache", dtw.ids),
        ("VAE latent cache", vae.ids),
    ):
        _require_same_ids(label, observed, canonical_call_ids)

    ids_by_logical_name: Mapping[str, Sequence[str]] = {
        "processed_calls": canonical_call_ids,
        "call_session_inventory": inventory_ids,
        "call_session_pairs": comparison_ids,
        "stp_pca_scores": stp.ids,
        "mfcc_pca_scores": mfcc.ids,
        "dtw_distances_condensed": dtw.ids,
        "vae_latent_means": vae.ids,
    }
    for logical_name, ids in ids_by_logical_name.items():
        _require_record_ids(manifest, input_names[logical_name], ids)

    _require_record_count(
        manifest, input_names["processed_calls"], "rows", len(calls)
    )
    _require_record_count(
        manifest, input_names["call_session_inventory"], "rows", len(inventory)
    )
    _require_record_count(
        manifest,
        input_names["call_session_inventory"],
        "eligible_units",
        len(membership),
    )
    _require_record_count(
        manifest, input_names["call_session_pairs"], "rows", len(comparisons)
    )
    for logical_name, table in (
        ("stp_pca_scores", stp_table),
        ("mfcc_pca_scores", mfcc_table),
    ):
        _require_record_count(
            manifest,
            input_names[logical_name],
            "canonical_unique_calls",
            len(table),
        )

    matrices = {"stp": stp, "mfcc": mfcc, "dtw": dtw, "vae": vae}
    distances = aggregate_item_distances(comparisons, membership, matrices)
    expected_rows = len(comparisons) * len(CALL_METRICS)
    if len(distances) != expected_rows:
        raise DataValidationError(
            f"Expected {expected_rows} call-distance rows, observed {len(distances)}"
        )
    model = build_model_table(distances, expected_metrics=CALL_METRICS, ddof=0)
    if len(model) != expected_rows:
        raise DataValidationError(
            f"Expected {expected_rows} call-model rows, observed {len(model)}"
        )

    _atomic_csv(distances, paths.call_repertoire_distances)
    _atomic_csv(model, paths.model_call)
    updated_manifest = _register_outputs(
        paths, manifest, input_names, distances, model
    )
    temporary_manifest = paths.manifest.with_name(f".{paths.manifest.name}.tmp")
    try:
        updated_manifest.write(temporary_manifest)
        temporary_manifest.replace(paths.manifest)
    finally:
        if temporary_manifest.exists():
            temporary_manifest.unlink()

    return {
        "status": "ok",
        "call_comparisons": int(len(comparisons)),
        "call_distance_rows": int(len(distances)),
        "call_model_rows": int(len(model)),
        "registered_artifacts": int(len(updated_manifest.artifacts)),
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        payload = assemble(args.root.resolve())
    except (MarmosetConvergenceError, OSError, ValueError, KeyError) as exc:
        print(json.dumps({"status": "error", "message": str(exc)}, indent=2))
        return 2
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
