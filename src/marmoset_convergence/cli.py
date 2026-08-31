"""Read-only validation entry points for cached and full workflows."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterable, Mapping, Optional

import numpy as np
import pandas as pd

from .calls import preflight_audio
from .exceptions import MarmosetConvergenceError
from .matrices import (
    load_dtw_distance_cache,
    load_sequence_distance_cache,
    load_sequence_representation_cache,
)
from .paths import ProjectPaths
from .provenance import (
    DTW_PARAMETERS,
    EMBEDDING_PARAMETERS,
    LOCAL_ALIGNMENT_PARAMETERS,
    SEQUENCE_REPRESENTATION_PARAMETERS,
    VAE_PARAMETERS,
    CacheManifest,
    ordered_ids_sha256,
    sha256_file,
)

_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _artifact_kind(name: str, path: str) -> Optional[str]:
    combined = f"{name} {path}".casefold()
    if "dtw_distances_condensed" in combined:
        return "dtw"
    if "sequence_distances" in combined:
        return "sequence_distances"
    if "sequence_representations" in combined:
        return "sequence_representations"
    if "vae_latent_means" in combined:
        return "vae"
    if "stp_pca_scores" in combined:
        return "stp_pca"
    if "mfcc_pca_scores" in combined:
        return "mfcc_pca"
    normalized_path = path.replace("\\", "/")
    if "/embeddings/" in f"/{normalized_path}" and normalized_path.endswith(".csv"):
        metric = Path(normalized_path).stem
        if metric in EMBEDDING_PARAMETERS:
            return f"embedding:{metric}"
    return None


def _expected_parameters(
    manifest: CacheManifest, names: Iterable[str]
) -> Mapping[str, Mapping[str, object]]:
    expected = {}
    for name in names:
        kind = _artifact_kind(name, manifest.artifacts[name].path)
        if kind == "dtw":
            expected[name] = DTW_PARAMETERS
        elif kind == "sequence_distances":
            expected[name] = LOCAL_ALIGNMENT_PARAMETERS
        elif kind == "sequence_representations":
            expected[name] = SEQUENCE_REPRESENTATION_PARAMETERS
        elif kind == "vae":
            expected[name] = VAE_PARAMETERS
        elif kind in {"stp_pca", "mfcc_pca"}:
            expected[name] = {"components": 5, "pca_refit": False}
        elif kind is not None and kind.startswith("embedding:"):
            expected[name] = EMBEDDING_PARAMETERS[kind.split(":", 1)[1]]
    return expected


def _validate_embedding_contract(
    root: Path, record_name: str, manifest: CacheManifest
) -> None:
    record = manifest.artifacts[record_name]
    kind = _artifact_kind(record_name, record.path)
    if kind is None or not kind.startswith("embedding:"):
        return
    metric = kind.split(":", 1)[1]
    is_sequence = metric in {
        "transition_probability",
        "bigram",
        "phee_repeat",
        "local_alignment",
    }
    id_column = "repertoire_id" if is_sequence else "call_id"
    canonical_path = (
        root / "data/cache/sequence_session_inventory.csv"
        if is_sequence
        else root / "data/processed/calls.csv"
    )
    grouping_columns = ["individual_id", "pair_id", "stage", "context"]
    canonical_table = pd.read_csv(
        canonical_path, usecols=[id_column, *grouping_columns]
    )
    canonical = canonical_table[id_column].astype(str)
    table = pd.read_csv(root / record.path)
    expected_columns = {
        id_column,
        "embedding_x",
        "embedding_y",
        *grouping_columns,
    }
    if set(table.columns) != expected_columns:
        raise ValueError(
            f"Embedding {record_name!r} must have exactly the canonical columns "
            + ", ".join(sorted(expected_columns))
        )
    observed = table[id_column].astype(str)
    if observed.duplicated().any():
        raise ValueError(f"Embedding {record_name!r} contains duplicate identifiers")
    if tuple(observed) != tuple(canonical):
        raise ValueError(
            f"Embedding {record_name!r} identifier order differs from the canonical input"
        )
    for column in grouping_columns:
        if tuple(table[column].astype(str)) != tuple(canonical_table[column].astype(str)):
            raise ValueError(
                f"Embedding {record_name!r} {column!r} values differ from canonical metadata"
            )
    coordinates = table[["embedding_x", "embedding_y"]].apply(
        pd.to_numeric, errors="coerce"
    )
    if not pd.notna(coordinates).all().all():
        raise ValueError(f"Embedding {record_name!r} coordinates must be finite numeric values")
    values = coordinates.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"Embedding {record_name!r} coordinates must be finite numeric values")
    if record.ordered_ids_sha256 != ordered_ids_sha256(observed) or record.n_ids != len(observed):
        raise ValueError(f"Embedding {record_name!r} ordered-ID provenance mismatch")
    parameters = record.parameters
    hyperparameters = parameters.get("hyperparameters")
    if not isinstance(hyperparameters, Mapping) or not hyperparameters:
        raise ValueError(f"Embedding {record_name!r} lacks a nonempty hyperparameters object")
    source_sha = str(parameters.get("source_artifact_sha256", "")).casefold()
    if not _SHA256.fullmatch(source_sha):
        raise ValueError(f"Embedding {record_name!r} lacks a valid source_artifact_sha256")
    source_by_metric = {
        "transition_probability": root / "data/cache/sequence_representations.npz",
        "bigram": root / "data/cache/sequence_representations.npz",
        "phee_repeat": root / "data/cache/sequence_representations.npz",
        "local_alignment": root / "data/cache/sequence_distances.npz",
        "stp": root / "data/cache/stp_pca_scores.csv",
        "mfcc": root / "data/cache/mfcc_pca_scores.csv",
        "dtw": root / "data/cache/dtw_distances_condensed.npz",
        "vae": root / "data/cache/vae_latent_means.csv",
    }
    source_path = source_by_metric[metric]
    if not source_path.is_file() or source_sha != sha256_file(source_path):
        raise ValueError(f"Embedding {record_name!r} source-artifact checksum mismatch")
    expected_width = {
        "transition_probability": 36,
        "bigram": 36,
        "phee_repeat": 4,
        "local_alignment": len(table),
        "stp": 5,
        "mfcc": 5,
        "dtw": len(table),
        "vae": 32,
    }[metric]
    if parameters.get("input_shape") != [len(table), expected_width] or parameters.get(
        "coordinate_shape"
    ) != [len(table), 2]:
        raise ValueError(f"Embedding {record_name!r} lacks valid input/coordinate shapes")


def _validate_npz_contracts(root: Path, manifest: CacheManifest, names: Iterable[str]) -> None:
    ids_by_kind = {}
    for name in names:
        record = manifest.artifacts[name]
        kind = _artifact_kind(name, record.path)
        path = root / record.path
        if kind == "dtw":
            matrix = load_dtw_distance_cache(path)
            ids_by_kind[kind] = matrix.ids
        elif kind == "sequence_distances":
            matrix = load_sequence_distance_cache(path)
            ids_by_kind[kind] = matrix.ids
            dynamic = {
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
            missing = sorted(dynamic.difference(record.parameters))
            if missing:
                raise ValueError(
                    f"Artifact {name!r} lacks dynamic local-alignment parameters: "
                    + ", ".join(missing)
                )
        elif kind == "sequence_representations":
            cache = load_sequence_representation_cache(path)
            ids_by_kind[kind] = cache.repertoire_ids
        elif kind in {"stp_pca", "mfcc_pca", "vae"}:
            table = pd.read_csv(path)
            if "call_id" not in table:
                raise ValueError(f"Artifact {name!r} lacks 'call_id'")
            ids = tuple(table["call_id"].astype(str))
            if len(set(ids)) != len(ids):
                raise ValueError(f"Artifact {name!r} contains duplicate call IDs")
            if kind in {"stp_pca", "mfcc_pca"}:
                prefix = "stp" if kind == "stp_pca" else "mfcc"
                feature_columns = [f"{prefix}_pc_{index}" for index in range(1, 6)]
                expected_columns = ["call_id", *feature_columns]
                if list(table.columns) != expected_columns:
                    raise ValueError(
                        f"Artifact {name!r} columns must be exactly "
                        + ", ".join(expected_columns)
                    )
            else:
                feature_columns = [column for column in table if column != "call_id"]
                if len(feature_columns) != 32 or len(table.columns) != 33:
                    raise ValueError(
                        f"Artifact {name!r} must contain call_id plus 32 latent means"
                    )
                recorded_columns = record.parameters.get("feature_columns")
                if recorded_columns != feature_columns:
                    raise ValueError(
                        f"Artifact {name!r} feature columns disagree with its manifest"
                    )
            numeric = table[feature_columns].apply(pd.to_numeric, errors="coerce")
            if not np.isfinite(numeric.to_numpy(dtype=float)).all():
                raise ValueError(f"Artifact {name!r} feature values must be finite")
            ids_by_kind[kind] = ids

        ids = ids_by_kind.get(kind)
        if ids is not None and (
            record.ordered_ids_sha256 != ordered_ids_sha256(ids)
            or record.n_ids != len(ids)
        ):
            raise ValueError(f"Artifact {name!r} ordered-ID provenance mismatch")
        if kind is not None and kind.startswith("embedding:"):
            _validate_embedding_contract(root, name, manifest)

    sequence_ids = ids_by_kind.get("sequence_distances")
    representation_ids = ids_by_kind.get("sequence_representations")
    if (
        sequence_ids is not None
        and representation_ids is not None
        and sequence_ids != representation_ids
    ):
        raise ValueError(
            "Sequence distance and representation caches use different repertoire order"
        )
    inventory_path = root / "data/cache/sequence_session_inventory.csv"
    if inventory_path.is_file():
        expected_repertoires = tuple(
            pd.read_csv(inventory_path, usecols=["repertoire_id"])[
                "repertoire_id"
            ].astype(str)
        )
        for kind in ("sequence_distances", "sequence_representations"):
            if kind in ids_by_kind and ids_by_kind[kind] != expected_repertoires:
                raise ValueError(
                    f"{kind.replace('_', ' ').title()} IDs differ from the canonical inventory"
                )
    calls_path = root / "data/processed/calls.csv"
    if calls_path.is_file():
        expected_calls = tuple(
            pd.read_csv(calls_path, usecols=["call_id"])["call_id"].astype(str)
        )
        for kind in ("dtw", "stp_pca", "mfcc_pca", "vae"):
            if kind in ids_by_kind and ids_by_kind[kind] != expected_calls:
                raise ValueError(
                    f"{kind.replace('_', ' ').upper()} IDs differ from canonical calls"
                )


def validate_cache(root: Path, required: Optional[Iterable[str]] = None) -> dict:
    paths = ProjectPaths.from_root(root)
    manifest = CacheManifest.load(paths.manifest)
    names = tuple(manifest.artifacts) if required is None else tuple(required)
    manifest.validate(
        paths.root,
        required_artifacts=names,
        expected_parameters=_expected_parameters(manifest, names),
    )
    _validate_npz_contracts(paths.root, manifest, names)
    return {"status": "ok", "validated_artifacts": list(names)}


def run_audio_preflight(
    root: Path,
    *,
    calls_path: Optional[Path] = None,
    audio_dir: Optional[Path] = None,
) -> dict:
    paths = ProjectPaths.from_root(root)
    calls_file = (calls_path or paths.calls).resolve()
    audio = (audio_dir or paths.audio).resolve()
    calls = pd.read_csv(calls_file)

    if paths.manifest.is_file():
        manifest = CacheManifest.load(paths.manifest)
        relative_calls = calls_file.relative_to(paths.root).as_posix()
        matching = [
            name
            for name, record in manifest.artifacts.items()
            if record.path == relative_calls
        ]
        if not matching:
            raise ValueError(f"Canonical call table {relative_calls!r} is not in manifest")
        name = matching[0]
        manifest.validate(
            paths.root,
            required_artifacts=[name],
            expected_ordered_ids={name: tuple(calls["call_id"].astype(str))},
        )

    report = preflight_audio(calls, audio)
    report.require_complete()
    return {
        "status": "ok",
        "expected_wavs": report.expected_count,
        "available_wavs": report.available_count,
        "unexpected_wavs": len(report.unexpected),
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="marmoset-repro",
        description="Validate cached artifacts or preflight a full recomputation.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    cached = subparsers.add_parser(
        "validate-cache", help="Validate manifest, checksums, parameters, IDs, and NPZ schemas"
    )
    cached.add_argument("--root", type=Path, default=Path.cwd())
    cached.add_argument(
        "--require",
        action="append",
        default=None,
        metavar="ARTIFACT",
        help="Require one manifest artifact (repeatable); default validates all",
    )

    audio = subparsers.add_parser(
        "preflight-audio", help="Fail unless every canonical call has one available WAV"
    )
    audio.add_argument("--root", type=Path, default=Path.cwd())
    audio.add_argument("--calls", type=Path, default=None)
    audio.add_argument("--audio-dir", type=Path, default=None)
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "validate-cache":
            payload = validate_cache(args.root, args.require)
        else:
            payload = run_audio_preflight(
                args.root, calls_path=args.calls, audio_dir=args.audio_dir
            )
    except (MarmosetConvergenceError, OSError, ValueError, KeyError) as exc:
        print(json.dumps({"status": "error", "message": str(exc)}, indent=2))
        return 2
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
