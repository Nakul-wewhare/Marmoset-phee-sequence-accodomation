#!/usr/bin/env python3
"""Validate the artifacts available in a checkout without changing them."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path, PurePosixPath

from marmoset_convergence.cli import validate_cache
from marmoset_convergence.provenance import (
    CURRENT_RELEASE_ARTIFACT_PATHS,
    CacheManifest,
    sha256_file,
)

LOCKED_RESULTS_ARTIFACTS = (
    "results/figures/main/figure_1.png",
    "results/figures/main/figure_2.png",
    "results/figures/main/figure_3.png",
    "results/figures/main/figure_4.png",
    "results/figures/supplement/figure_s1_workflow.png",
    "results/figures/supplement/figure_s2_bigram.png",
    "results/figures/supplement/figure_s3_local_alignment.png",
    "results/figures/supplement/figure_s4_phee_repeat.png",
    "results/figures/supplement/figure_s5_mfcc.png",
    "results/figures/supplement/figure_s6_stp.png",
    "results/figures/supplement/figure_s7_dtw.png",
    "results/figures/supplement/figure_s8_vae.png",
    "results/figures/figure_provenance.csv",
    "results/reference/manuscript_reported_stage_contrasts.csv",
)


def _relative_file(root: Path, value: object) -> Path:
    text = str(value).strip().replace("\\", "/")
    path = PurePosixPath(text)
    if not text or path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Manifest path is not repository-relative: {value!r}")
    resolved = (root / path.as_posix()).resolve()
    resolved.relative_to(root)
    if not resolved.is_file():
        raise FileNotFoundError(f"Registered artifact is missing: {path.as_posix()}")
    return resolved


def validate_processed_manifest(root: Path) -> int:
    path = root / "data/processed/manifest.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("kind") != "processed-input-manifest":
        raise ValueError("Unexpected processed manifest kind")
    records = [*payload.get("sources", []), *payload.get("artifacts", [])]
    if not records:
        raise ValueError("Processed manifest contains no file records")
    for record in records:
        artifact = _relative_file(root, record["path"])
        if sha256_file(artifact) != record["sha256"]:
            raise ValueError(f"SHA-256 mismatch: {record['path']}")
        if "rows" in record:
            with artifact.open("r", encoding="utf-8", newline="") as handle:
                observed = sum(1 for _ in csv.reader(handle)) - 1
            if observed != int(record["rows"]):
                raise ValueError(
                    f"Row-count mismatch for {record['path']}: {observed} != {record['rows']}"
                )
    return len(records)


def validate_results_manifest(root: Path) -> int:
    path = root / "results/manifest.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("kind") != "manuscript-artifact-manifest":
        raise ValueError("Unexpected results manifest kind")
    records = payload.get("artifacts", [])
    if not records:
        raise ValueError("Results manifest contains no artifacts")
    seen = set()
    checksums = {}
    for record in records:
        relative = str(record["path"])
        if relative in seen:
            raise ValueError(f"Duplicate results-manifest path: {relative}")
        seen.add(relative)
        artifact = _relative_file(root, relative)
        if sha256_file(artifact) != record["sha256"]:
            raise ValueError(f"SHA-256 mismatch: {relative}")
        checksums[relative] = str(record["sha256"])
    missing = sorted(set(LOCKED_RESULTS_ARTIFACTS).difference(seen))
    if missing:
        raise ValueError(
            "Results manifest lacks required locked artifacts: " + ", ".join(missing)
        )

    provenance_path = root / "results/figures/figure_provenance.csv"
    with provenance_path.open("r", encoding="utf-8", newline="") as handle:
        provenance_rows = list(csv.DictReader(handle))
    expected_pngs = {
        path for path in LOCKED_RESULTS_ARTIFACTS if path.casefold().endswith(".png")
    }
    observed_pngs = {row["path"] for row in provenance_rows}
    if observed_pngs != expected_pngs:
        raise ValueError("Figure provenance index does not list exactly Figures 1–4 and S1–S8")
    for row in provenance_rows:
        if row["sha256"] != checksums[row["path"]]:
            raise ValueError(f"Figure provenance checksum disagrees for {row['path']}")
    return len(records)


def _iter_manifest_records(payload: dict):
    for section_name in ("artifacts", "sources", "inputs", "files", "outputs"):
        section = payload.get(section_name, {})
        if isinstance(section, dict):
            for key, value in section.items():
                record = value if isinstance(value, dict) else {"sha256": value}
                yield str(record.get("path", key)), record
        elif isinstance(section, list):
            for record in section:
                if isinstance(record, dict) and record.get("path"):
                    yield str(record["path"]), record


def _record_checksum(record: dict) -> str | None:
    checksum = record.get("sha256") or record.get("checksum_sha256")
    if checksum:
        return str(checksum).removeprefix("sha256:").casefold()
    value = record.get("checksum")
    if isinstance(value, str):
        return value.removeprefix("sha256:").casefold()
    if isinstance(value, dict) and str(value.get("algorithm", "")).casefold() == "sha256":
        return str(value.get("value", "")).casefold()
    return None


def validate_generic_manifest(root: Path, path: Path) -> set[str]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    registered = set()
    for relative, record in _iter_manifest_records(payload):
        artifact = _relative_file(root, relative)
        expected = _record_checksum(record)
        if not expected:
            raise ValueError(f"No SHA-256 for {relative} in {path.relative_to(root)}")
        if sha256_file(artifact) != expected:
            raise ValueError(f"SHA-256 mismatch: {relative}")
        registered.add(PurePosixPath(relative.replace("\\", "/")).as_posix())
    if not registered:
        raise ValueError(f"Manifest contains no file records: {path.relative_to(root)}")
    return registered


def registered_artifact_paths(root: Path) -> set[str]:
    registered = set()
    for manifest_path in (
        root / "data/processed/manifest.json",
        root / "data/derived/manifest.json",
        root / "results/manifest.json",
        root / "results/model_manifest.json",
    ):
        if manifest_path.is_file():
            registered.update(validate_generic_manifest(root, manifest_path))
    cache_manifest = CacheManifest.load(root / "data/cache/manifest.json")
    registered.update(record.path for record in cache_manifest.artifacts.values())
    return registered


def validate_notebook_hygiene(root: Path) -> int:
    notebooks = sorted((root / "notebooks").glob("[0-9][0-9]_*.ipynb"))
    if len(notebooks) != 5:
        raise ValueError(f"Expected five numbered notebooks, found {len(notebooks)}")
    for path in notebooks:
        payload = json.loads(path.read_text(encoding="utf-8"))
        for index, cell in enumerate(payload.get("cells", [])):
            if cell.get("cell_type") == "code":
                if cell.get("outputs"):
                    raise ValueError(f"Notebook output is not stripped: {path}:{index}")
                if cell.get("execution_count") is not None:
                    raise ValueError(
                        f"Notebook execution count is not stripped: {path}:{index}"
                    )
    return len(notebooks)


def validate(root: Path, *, strict: bool) -> dict[str, object]:
    root = root.resolve()
    processed_records = validate_processed_manifest(root)
    cache_result = validate_cache(root)
    results_records = validate_results_manifest(root)
    notebooks = validate_notebook_hygiene(root)

    cache_manifest = CacheManifest.load(root / "data/cache/manifest.json")
    completeness = str(cache_manifest.metadata.get("release_completeness", "unknown"))
    missing = list(cache_manifest.metadata.get("missing_current_release_artifacts", []))
    registered = registered_artifact_paths(root)
    if strict:
        missing_files = [
            relative
            for relative in CURRENT_RELEASE_ARTIFACT_PATHS
            if not (root / relative).is_file()
        ]
        unregistered = [
            relative
            for relative in CURRENT_RELEASE_ARTIFACT_PATHS
            if (root / relative).is_file()
            if relative != "results/model_manifest.json"
            and relative not in registered
        ]
        problems = []
        if missing_files:
            problems.append("missing files:\n  - " + "\n  - ".join(missing_files))
        if unregistered:
            problems.append(
                "files absent from all provenance manifests:\n  - "
                + "\n  - ".join(unregistered)
            )
        if completeness != "complete" or missing:
            problems.append(
                "data/cache/manifest.json does not declare a complete release; "
                "regenerate it after installing and registering every asset"
            )
        if problems:
            raise ValueError(
                "A complete current-analysis release was requested, but the "
                "explicit release contract failed:\n" + "\n".join(problems)
            )
    return {
        "status": "ok",
        "validation_scope": "complete-current-release" if strict else "available-artifacts",
        "processed_manifest_records": processed_records,
        "cache_artifacts": len(cache_result["validated_artifacts"]),
        "results_manifest_records": results_records,
        "numbered_notebooks": notebooks,
        "release_completeness": completeness,
        "missing_current_release_artifacts": missing,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail unless the manifest declares a complete current-analysis release.",
    )
    args = parser.parse_args()
    try:
        payload = validate(args.root, strict=args.strict)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(json.dumps({"status": "error", "message": str(exc)}, indent=2))
        return 2
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
