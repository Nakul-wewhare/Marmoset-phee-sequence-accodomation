"""Checksummed cache manifests and fixed method-parameter contracts."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple, Union

from .exceptions import CacheValidationError

PathLike = Union[str, Path]
MANIFEST_SCHEMA_VERSION = 1

# Numerical artifacts required to execute all five reader-facing notebooks.
# Locked manuscript PNGs are validated separately and are not numerical inputs.
CURRENT_RELEASE_ARTIFACT_PATHS: Tuple[str, ...] = (
    "data/processed/sequences.csv",
    "data/processed/calls.csv",
    "data/cache/sequence_session_inventory.csv",
    "data/cache/call_session_inventory.csv",
    "data/cache/sequence_session_pairs.csv",
    "data/cache/call_session_pairs.csv",
    "data/cache/sequence_representations.npz",
    "data/cache/sequence_distances.npz",
    "data/cache/stp_pca_scores.csv",
    "data/cache/mfcc_pca_scores.csv",
    "data/cache/dtw_distances_condensed.npz",
    "data/cache/vae_latent_means.csv",
    "data/cache/embeddings/transition_probability.csv",
    "data/cache/embeddings/bigram.csv",
    "data/cache/embeddings/phee_repeat.csv",
    "data/cache/embeddings/local_alignment.csv",
    "data/cache/embeddings/stp.csv",
    "data/cache/embeddings/mfcc.csv",
    "data/cache/embeddings/dtw.csv",
    "data/cache/embeddings/vae.csv",
    "data/derived/sequence_repertoire_distances.csv",
    "data/derived/call_repertoire_distances.csv",
    "data/derived/model_sequence.csv",
    "data/derived/model_call.csv",
    "results/posterior_draws/sequence_stage_effect_draws.csv.gz",
    "results/posterior_draws/sequence_stage_effect_draws.csv.gz.metadata.json",
    "results/posterior_draws/call_stage_effect_draws.csv.gz",
    "results/posterior_draws/call_stage_effect_draws.csv.gz.metadata.json",
    "results/tables/stage_effects_by_context_metric.csv",
    "results/tables/stage_effects_combined.csv",
    "results/tables/psi_by_metric.csv",
    "results/tables/psi_combined.csv",
    "results/tables/model_diagnostics.csv",
    "results/figures/main/diagnostics/sequence/sequence_interaction_model_pp_check.png",
    "results/figures/main/diagnostics/sequence/sequence_interaction_model_trace.png",
    "results/figures/main/diagnostics/call/call_interaction_model_pp_check.png",
    "results/figures/main/diagnostics/call/call_interaction_model_trace.png",
    "results/model_manifest.json",
)

DTW_PARAMETERS: Mapping[str, Any] = {
    "algorithm": "exact_dynamic_programming",
    "frequency_min_hz": 4000,
    "frequency_max_hz": 12000,
    "window_samples": 2048,
    "hop_samples": 512,
    "db_min": -80,
    "db_max": 0,
    "local_cost": "cosine",
    "normalized_time_constraint": 0.20,
    "path_normalization": "alignment_steps",
    "storage": "condensed_scipy",
}

VAE_PARAMETERS: Mapping[str, Any] = {
    "implementation": "chatter",
    "version": "0.1.6",
    "image_size": [128, 128],
    "canvas_seconds": 2.5,
    "latent_dim": 32,
    "epochs": 100,
    "batch_size": 32,
    "learning_rate": 0.0001,
    "beta": 0.5,
    "seed": 42,
}

# Coordinate caches are descriptive outputs, not model inputs.  Their records
# still bind the pooled fit, implementation, seed, and source population so a
# context-specific or differently ordered fit cannot masquerade as the locked
# analysis.  The full algorithm-specific tuning dictionary and source checksum
# are required by the semantic validator in ``cli.py``.
EMBEDDING_PARAMETERS: Mapping[str, Mapping[str, Any]] = {
    "transition_probability": {
        "algorithm": "pacmap",
        "package": "pacmap",
        "package_version": "0.8.0",
        "random_seed": 42,
        "n_components": 2,
        "fit_scope": "pooled_across_stage_and_context",
        "input_scope": "107_eligible_sequence_repertoires",
        "input_metric": "euclidean",
    },
    "bigram": {
        "algorithm": "pacmap",
        "package": "pacmap",
        "package_version": "0.8.0",
        "random_seed": 42,
        "n_components": 2,
        "fit_scope": "pooled_across_stage_and_context",
        "input_scope": "107_eligible_sequence_repertoires",
        "input_metric": "euclidean",
    },
    "phee_repeat": {
        "algorithm": "umap",
        "package": "umap-learn",
        "package_version": "0.5.9.post2",
        "random_seed": 42,
        "n_components": 2,
        "fit_scope": "pooled_across_stage_and_context",
        "input_scope": "107_eligible_sequence_repertoires",
        "input_metric": "euclidean",
    },
    "local_alignment": {
        "algorithm": "umap",
        "package": "umap-learn",
        "package_version": "0.5.9.post2",
        "random_seed": 42,
        "n_components": 2,
        "fit_scope": "pooled_across_stage_and_context",
        "input_scope": "107_eligible_sequence_repertoires",
        "input_metric": "precomputed",
    },
    "stp": {
        "algorithm": "pacmap",
        "package": "pacmap",
        "package_version": "0.8.0",
        "random_seed": 42,
        "n_components": 2,
        "fit_scope": "pooled_across_stage_and_context",
        "input_scope": "3612_unique_calls",
        "input_metric": "euclidean",
    },
    "mfcc": {
        "algorithm": "pacmap",
        "package": "pacmap",
        "package_version": "0.8.0",
        "random_seed": 42,
        "n_components": 2,
        "fit_scope": "pooled_across_stage_and_context",
        "input_scope": "3612_unique_calls",
        "input_metric": "euclidean",
    },
    "dtw": {
        "algorithm": "umap",
        "package": "umap-learn",
        "package_version": "0.5.9.post2",
        "random_seed": 42,
        "n_components": 2,
        "fit_scope": "pooled_across_stage_and_context",
        "input_scope": "3612_unique_calls",
        "input_metric": "precomputed",
    },
    "vae": {
        "algorithm": "pacmap",
        "package": "pacmap",
        "package_version": "0.8.0",
        "random_seed": 42,
        "n_components": 2,
        "fit_scope": "pooled_across_stage_and_context",
        "input_scope": "3612_unique_calls",
        "input_metric": "euclidean",
    },
}

LOCAL_ALIGNMENT_PARAMETERS: Mapping[str, Any] = {
    "algorithm": "smith_waterman",
    "match": 2,
    "mismatch": -1,
    "gap": -1,
    "normalization": "log_combined_lengths",
    "aggregation": "mean_all_cross_repertoire_pairs",
    "distance_transform": "global_max_minus_similarity",
    "storage": "condensed_scipy",
}

SEQUENCE_REPRESENTATION_PARAMETERS: Mapping[str, Any] = {
    "alphabet": ["A", "K", "P", "R", "S", "T"],
    "transition_normalization": "row_conditional",
    "bigram_normalization": "global_total",
    "phee_token": "A",
    "repeat_min": 2,
    "repeat_max": 5,
    "longer_runs": "exclude",
}

MODEL_TABLE_PARAMETERS: Mapping[str, Any] = {
    "standardization": "within_metric_joint_all_pairs_stages_contexts",
    "standard_deviation_ddof": 0,
}

_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def sha256_file(path: PathLike, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def ordered_ids_sha256(ids: Iterable[Any]) -> str:
    """Hash ordered IDs using an unambiguous UTF-8 JSON encoding."""

    payload = json.dumps(
        [str(identifier) for identifier in ids],
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _portable_relative_path(value: Any) -> str:
    text = str(value).strip().replace("\\", "/")
    path = PurePosixPath(text)
    if not text or path.is_absolute() or ".." in path.parts:
        raise CacheValidationError(
            f"Manifest artifact paths must be project-relative: {value!r}"
        )
    return path.as_posix()


def _require_parameter_subset(
    actual: Mapping[str, Any], expected: Mapping[str, Any], *, artifact: str
) -> None:
    for key, expected_value in expected.items():
        if key not in actual:
            raise CacheValidationError(
                f"Artifact {artifact!r} lacks required parameter {key!r}"
            )
        actual_value = actual[key]
        if isinstance(expected_value, Mapping):
            if not isinstance(actual_value, Mapping):
                raise CacheValidationError(
                    f"Artifact {artifact!r} parameter {key!r} must be an object"
                )
            _require_parameter_subset(
                actual_value, expected_value, artifact=f"{artifact}.{key}"
            )
        elif actual_value != expected_value:
            raise CacheValidationError(
                f"Artifact {artifact!r} parameter {key!r} is {actual_value!r}; "
                f"expected {expected_value!r}"
            )


@dataclass(frozen=True)
class ArtifactRecord:
    path: str
    sha256: str
    size_bytes: int
    parameters: Mapping[str, Any] = field(default_factory=dict)
    ordered_ids_sha256: Optional[str] = None
    n_ids: Optional[int] = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", _portable_relative_path(self.path))
        checksum = str(self.sha256).casefold()
        if not _SHA256.fullmatch(checksum):
            raise CacheValidationError(f"Invalid SHA-256 for artifact {self.path!r}")
        object.__setattr__(self, "sha256", checksum)
        if not isinstance(self.size_bytes, int) or self.size_bytes < 0:
            raise CacheValidationError(f"Invalid size_bytes for artifact {self.path!r}")
        if not isinstance(self.parameters, Mapping):
            raise CacheValidationError(f"Artifact {self.path!r} parameters must be an object")
        if self.ordered_ids_sha256 is not None:
            ids_checksum = str(self.ordered_ids_sha256).casefold()
            if not _SHA256.fullmatch(ids_checksum):
                raise CacheValidationError(
                    f"Invalid ordered_ids_sha256 for artifact {self.path!r}"
                )
            object.__setattr__(self, "ordered_ids_sha256", ids_checksum)
        if self.n_ids is not None and (
            not isinstance(self.n_ids, int) or self.n_ids < 0
        ):
            raise CacheValidationError(f"Invalid n_ids for artifact {self.path!r}")

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ArtifactRecord":
        required = {"path", "sha256", "size_bytes", "parameters"}
        missing = sorted(required.difference(data))
        unknown = sorted(
            set(data).difference(required | {"ordered_ids_sha256", "n_ids"})
        )
        if missing or unknown:
            detail = []
            if missing:
                detail.append("missing " + ", ".join(missing))
            if unknown:
                detail.append("unknown " + ", ".join(unknown))
            raise CacheValidationError("Invalid artifact record: " + "; ".join(detail))
        return cls(
            path=data["path"],
            sha256=data["sha256"],
            size_bytes=data["size_bytes"],
            parameters=data["parameters"],
            ordered_ids_sha256=data.get("ordered_ids_sha256"),
            n_ids=data.get("n_ids"),
        )

    def to_dict(self) -> Dict[str, Any]:
        data: Dict[str, Any] = {
            "path": self.path,
            "sha256": self.sha256,
            "size_bytes": self.size_bytes,
            "parameters": dict(self.parameters),
        }
        if self.ordered_ids_sha256 is not None:
            data["ordered_ids_sha256"] = self.ordered_ids_sha256
        if self.n_ids is not None:
            data["n_ids"] = self.n_ids
        return data


@dataclass(frozen=True)
class CacheManifest:
    artifacts: Mapping[str, ArtifactRecord]
    schema_version: int = MANIFEST_SCHEMA_VERSION
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.schema_version != MANIFEST_SCHEMA_VERSION:
            raise CacheValidationError(
                f"Unsupported manifest schema_version {self.schema_version}; "
                f"expected {MANIFEST_SCHEMA_VERSION}"
            )
        if not self.artifacts:
            raise CacheValidationError("Manifest must register at least one artifact")
        if not isinstance(self.metadata, Mapping):
            raise CacheValidationError("Manifest metadata must be an object")

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "CacheManifest":
        allowed = {"schema_version", "artifacts", "metadata"}
        unknown = sorted(set(data).difference(allowed))
        if unknown:
            raise CacheValidationError("Unknown manifest fields: " + ", ".join(unknown))
        if "schema_version" not in data or "artifacts" not in data:
            raise CacheValidationError("Manifest requires schema_version and artifacts")
        raw_artifacts = data["artifacts"]
        if not isinstance(raw_artifacts, Mapping):
            raise CacheValidationError("Manifest artifacts must be an object")
        artifacts = {
            str(name): ArtifactRecord.from_dict(record)
            for name, record in raw_artifacts.items()
        }
        return cls(
            artifacts=artifacts,
            schema_version=data["schema_version"],
            metadata=data.get("metadata", {}),
        )

    @classmethod
    def load(cls, path: PathLike) -> "CacheManifest":
        manifest_path = Path(path)
        if not manifest_path.is_file():
            raise CacheValidationError(f"Cache manifest is missing: {manifest_path}")
        try:
            with manifest_path.open("r", encoding="utf-8") as handle:
                payload = json.load(handle)
        except (OSError, json.JSONDecodeError) as exc:
            raise CacheValidationError(
                f"Could not parse cache manifest {manifest_path}: {exc}"
            ) from exc
        if not isinstance(payload, Mapping):
            raise CacheValidationError("Cache manifest root must be a JSON object")
        return cls.from_dict(payload)

    @classmethod
    def from_files(
        cls,
        root: PathLike,
        artifacts: Mapping[str, PathLike],
        *,
        parameters: Optional[Mapping[str, Mapping[str, Any]]] = None,
        ordered_ids: Optional[Mapping[str, Sequence[Any]]] = None,
        metadata: Optional[Mapping[str, Any]] = None,
    ) -> "CacheManifest":
        """Construct records from files; writing remains a separate explicit step."""

        root_path = Path(root).expanduser().resolve()
        records = {}
        for name, raw_path in artifacts.items():
            artifact_path = Path(raw_path)
            if not artifact_path.is_absolute():
                artifact_path = root_path / artifact_path
            artifact_path = artifact_path.resolve()
            try:
                relative = artifact_path.relative_to(root_path).as_posix()
            except ValueError as exc:
                raise CacheValidationError(
                    f"Artifact {artifact_path} is outside project root {root_path}"
                ) from exc
            if not artifact_path.is_file():
                raise CacheValidationError(f"Cannot register missing artifact {artifact_path}")
            ids = None if ordered_ids is None else ordered_ids.get(name)
            records[str(name)] = ArtifactRecord(
                path=relative,
                sha256=sha256_file(artifact_path),
                size_bytes=artifact_path.stat().st_size,
                parameters={} if parameters is None else parameters.get(name, {}),
                ordered_ids_sha256=None if ids is None else ordered_ids_sha256(ids),
                n_ids=None if ids is None else len(ids),
            )
        return cls(records, metadata={} if metadata is None else metadata)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "metadata": dict(self.metadata),
            "artifacts": {
                name: record.to_dict() for name, record in sorted(self.artifacts.items())
            },
        }

    def write(self, path: PathLike) -> None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("w", encoding="utf-8") as handle:
            json.dump(self.to_dict(), handle, indent=2, sort_keys=True)
            handle.write("\n")

    def validate(
        self,
        root: PathLike,
        *,
        required_artifacts: Optional[Iterable[str]] = None,
        expected_parameters: Optional[Mapping[str, Mapping[str, Any]]] = None,
        expected_ordered_ids: Optional[Mapping[str, Sequence[Any]]] = None,
    ) -> "CacheManifest":
        """Verify registration, path confinement, size, checksum, IDs, and parameters."""

        root_path = Path(root).expanduser().resolve()
        names = (
            tuple(self.artifacts)
            if required_artifacts is None
            else tuple(str(name) for name in required_artifacts)
        )
        missing_records = sorted(set(names).difference(self.artifacts))
        if missing_records:
            raise CacheValidationError(
                "Artifacts are not registered in manifest: " + ", ".join(missing_records)
            )
        for name in names:
            record = self.artifacts[name]
            path = (root_path / record.path).resolve()
            try:
                path.relative_to(root_path)
            except ValueError as exc:
                raise CacheValidationError(
                    f"Artifact {name!r} resolves outside project root"
                ) from exc
            if not path.is_file():
                raise CacheValidationError(f"Artifact {name!r} is missing: {path}")
            observed_size = path.stat().st_size
            if observed_size != record.size_bytes:
                raise CacheValidationError(
                    f"Artifact {name!r} size is {observed_size}; expected {record.size_bytes}"
                )
            observed_checksum = sha256_file(path)
            if observed_checksum != record.sha256:
                raise CacheValidationError(f"Artifact {name!r} SHA-256 mismatch")
            if expected_parameters is not None and name in expected_parameters:
                _require_parameter_subset(
                    record.parameters, expected_parameters[name], artifact=name
                )
            if expected_ordered_ids is not None and name in expected_ordered_ids:
                if record.ordered_ids_sha256 is None or record.n_ids is None:
                    raise CacheValidationError(
                        f"Artifact {name!r} does not record ordered-ID provenance"
                    )
                ids = expected_ordered_ids[name]
                if record.n_ids != len(ids):
                    raise CacheValidationError(
                        f"Artifact {name!r} records {record.n_ids} IDs; received {len(ids)}"
                    )
                if record.ordered_ids_sha256 != ordered_ids_sha256(ids):
                    raise CacheValidationError(
                        f"Artifact {name!r} ordered-ID checksum mismatch"
                    )
        return self


def load_and_validate_manifest(
    manifest_path: PathLike,
    *,
    root: Optional[PathLike] = None,
    required_artifacts: Optional[Iterable[str]] = None,
    expected_parameters: Optional[Mapping[str, Mapping[str, Any]]] = None,
    expected_ordered_ids: Optional[Mapping[str, Sequence[Any]]] = None,
) -> CacheManifest:
    """Reader-facing one-call validation used before loading any cache."""

    path = Path(manifest_path).expanduser().resolve()
    project_root = path.parents[2] if root is None else Path(root)
    manifest = CacheManifest.load(path)
    return manifest.validate(
        project_root,
        required_artifacts=required_artifacts,
        expected_parameters=expected_parameters,
        expected_ordered_ids=expected_ordered_ids,
    )
