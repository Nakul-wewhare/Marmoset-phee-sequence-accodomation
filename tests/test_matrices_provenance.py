import json
from pathlib import Path

import numpy as np
import pytest

from marmoset_convergence import (
    CacheManifest,
    CacheValidationError,
    CondensedDistanceMatrix,
    DataValidationError,
    FeatureDistanceMatrix,
    condensed_index,
    condensed_to_square,
    load_dtw_distance_cache,
    load_sequence_distance_cache,
    ordered_ids_sha256,
    square_to_condensed,
)
from marmoset_convergence.matrices import (
    save_dtw_distance_cache,
    save_sequence_distance_cache,
)


def test_condensed_helpers_match_scipy_upper_triangle_order():
    square = np.array(
        [
            [0, 1, 2, 3],
            [1, 0, 4, 5],
            [2, 4, 0, 6],
            [3, 5, 6, 0],
        ],
        dtype=float,
    )

    condensed = square_to_condensed(square)

    np.testing.assert_array_equal(condensed, [1, 2, 3, 4, 5, 6])
    np.testing.assert_array_equal(condensed_to_square(condensed), square)
    assert condensed_index(4, 2, 0) == 1
    assert condensed_index(4, 1, 3) == 4


def test_condensed_matrix_validates_length_and_ordered_ids():
    matrix = CondensedDistanceMatrix(("a", "b", "c"), np.array([1, 2, 3]), "x")
    assert matrix.distance("c", "a") == 2
    assert matrix.distance("b", "b") == 0
    assert matrix.mean_between(["a"], ["b", "c"]) == pytest.approx(1.5)

    with pytest.raises(DataValidationError, match="has 2 distances"):
        CondensedDistanceMatrix(("a", "b", "c"), np.array([1, 2]), "x")


def test_feature_distance_matrix_averages_only_requested_cross_blocks():
    euclidean = FeatureDistanceMatrix(
        ("a1", "a2", "b1", "b2"),
        np.array([[0.0], [2.0], [4.0], [8.0]]),
        "stp",
        "euclidean",
    )
    assert euclidean.mean_between(["a1", "a2"], ["b1", "b2"]) == pytest.approx(
        (4 + 8 + 2 + 6) / 4
    )

    cosine = FeatureDistanceMatrix(
        ("x", "y", "z"),
        np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]]),
        "vae",
        "cosine",
    )
    assert cosine.mean_between(["x"], ["y", "z"]) == pytest.approx(1.5)

    with pytest.raises(DataValidationError, match="zero-norm"):
        FeatureDistanceMatrix(("x", "zero"), np.array([[1.0], [0.0]]), "vae", "cosine")


def test_dtw_npz_contract_round_trip_and_rejects_extra_keys(tmp_path: Path):
    path = tmp_path / "dtw.npz"
    save_dtw_distance_cache(path, ["a", "b", "c"], np.array([1.0, 2.0, 3.0]))
    loaded = load_dtw_distance_cache(path)
    assert loaded.ids == ("a", "b", "c")

    np.savez(path, call_ids=np.array(["a", "b"]), distances=np.array([1.0]), extra=1)
    with pytest.raises(CacheValidationError, match="cache keys"):
        load_dtw_distance_cache(path)


def test_sequence_distance_npz_uses_documented_keys(tmp_path: Path):
    path = tmp_path / "sequence.npz"
    save_sequence_distance_cache(path, ["r1", "r2"], np.array([0.25]))

    with np.load(path, allow_pickle=False) as archive:
        assert set(archive.files) == {"repertoire_ids", "local_alignment"}
    assert load_sequence_distance_cache(path).distance("r1", "r2") == pytest.approx(0.25)


def test_manifest_validates_checksum_parameters_and_id_order(tmp_path: Path):
    artifact = tmp_path / "data" / "cache" / "tiny.bin"
    artifact.parent.mkdir(parents=True)
    artifact.write_bytes(b"cached result")
    ids = ["a", "b"]
    manifest = CacheManifest.from_files(
        tmp_path,
        {"tiny": artifact},
        parameters={"tiny": {"algorithm": "test", "nested": {"x": 1}}},
        ordered_ids={"tiny": ids},
    )
    manifest_path = tmp_path / "data" / "cache" / "manifest.json"
    manifest.write(manifest_path)

    loaded = CacheManifest.load(manifest_path)
    loaded.validate(
        tmp_path,
        required_artifacts=["tiny"],
        expected_parameters={"tiny": {"algorithm": "test", "nested": {"x": 1}}},
        expected_ordered_ids={"tiny": ids},
    )
    assert loaded.artifacts["tiny"].ordered_ids_sha256 == ordered_ids_sha256(ids)

    with pytest.raises(CacheValidationError, match="ordered-ID checksum"):
        loaded.validate(
            tmp_path,
            required_artifacts=["tiny"],
            expected_ordered_ids={"tiny": list(reversed(ids))},
        )


def test_manifest_detects_tampered_artifact(tmp_path: Path):
    artifact = tmp_path / "cache.bin"
    artifact.write_bytes(b"first")
    manifest = CacheManifest.from_files(tmp_path, {"cache": artifact})
    artifact.write_bytes(b"other")

    with pytest.raises(CacheValidationError, match="SHA-256 mismatch"):
        manifest.validate(tmp_path)


def test_manifest_rejects_paths_outside_project():
    payload = {
        "schema_version": 1,
        "artifacts": {
            "bad": {
                "path": "../outside.bin",
                "sha256": "0" * 64,
                "size_bytes": 0,
                "parameters": {},
            }
        },
    }
    with pytest.raises(CacheValidationError, match="project-relative"):
        CacheManifest.from_dict(payload)
