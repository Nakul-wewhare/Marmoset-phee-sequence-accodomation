import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from marmoset_convergence import (
    CacheManifest,
    CacheValidationError,
    DTW_PARAMETERS,
    VAE_PARAMETERS,
    square_to_condensed,
)
from marmoset_convergence.matrices import save_dtw_distance_cache


SCRIPT = Path(__file__).parents[1] / "scripts" / "assemble_call_model_table.py"
SPEC = importlib.util.spec_from_file_location("assemble_call_model_table", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
ASSEMBLER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ASSEMBLER)


def _synthetic_cache_tree(root: Path) -> None:
    processed = root / "data" / "processed"
    cache = root / "data" / "cache"
    processed.mkdir(parents=True)
    cache.mkdir(parents=True)

    repertoire_ids = ("r1", "r2", "r3", "r4")
    calls = []
    for repertoire_id in repertoire_ids:
        for index in range(5):
            calls.append(
                {
                    "call_id": f"{repertoire_id}_call_{index}",
                    "repertoire_id": repertoire_id,
                }
            )
    calls = pd.DataFrame(calls)
    call_ids = tuple(calls["call_id"])
    calls.to_csv(processed / "calls.csv", index=False)

    inventory = pd.DataFrame(
        {
            "repertoire_id": repertoire_ids,
            "pair_id": ["A"] * 4,
            "stage": ["before", "before", "after", "after"],
            "context": ["partner"] * 4,
            "n_calls": [5] * 4,
        }
    )
    inventory.to_csv(cache / "call_session_inventory.csv", index=False)
    comparisons = pd.DataFrame(
        [
            {
                "comparison_id": "before_pair",
                "pair_id": "A",
                "stage": "before",
                "context": "partner",
                "repertoire_a_id": "r1",
                "repertoire_b_id": "r2",
            },
            {
                "comparison_id": "after_pair",
                "pair_id": "A",
                "stage": "after",
                "context": "partner",
                "repertoire_a_id": "r3",
                "repertoire_b_id": "r4",
            },
        ]
    )
    comparisons.to_csv(cache / "call_session_pairs.csv", index=False)

    repertoire_for_call = calls["repertoire_id"].map(
        {"r1": 0.0, "r2": 1.0, "r3": 0.0, "r4": 3.0}
    )
    stp = pd.DataFrame({"call_id": call_ids})
    mfcc = pd.DataFrame({"call_id": call_ids})
    for index in range(1, 6):
        stp[f"stp_pc_{index}"] = repertoire_for_call if index == 1 else 0.0
        mfcc[f"mfcc_pc_{index}"] = (
            calls["repertoire_id"].map(
                {"r1": 0.0, "r2": 2.0, "r3": 0.0, "r4": 5.0}
            )
            if index == 1
            else 0.0
        )
    stp.to_csv(cache / "stp_pca_scores.csv", index=False)
    mfcc.to_csv(cache / "mfcc_pca_scores.csv", index=False)

    vae = pd.DataFrame({"call_id": call_ids})
    latent_centers = {
        "r1": (1.0, 0.0),
        "r2": (0.0, 1.0),
        "r3": (1.0, 0.0),
        "r4": (-1.0, 0.0),
    }
    for index in range(1, 33):
        if index <= 2:
            vae[f"latent_{index}"] = [
                latent_centers[repertoire_id][index - 1]
                for repertoire_id in calls["repertoire_id"]
            ]
        else:
            vae[f"latent_{index}"] = 0.0
    vae.to_csv(cache / "vae_latent_means.csv", index=False)

    dtw_square = np.ones((len(call_ids), len(call_ids)), dtype=float)
    np.fill_diagonal(dtw_square, 0.0)
    repertoire_positions = {
        repertoire_id: np.flatnonzero(
            calls["repertoire_id"].to_numpy() == repertoire_id
        )
        for repertoire_id in repertoire_ids
    }
    for left, right, value in (("r1", "r2", 2.0), ("r3", "r4", 4.0)):
        left_positions = repertoire_positions[left]
        right_positions = repertoire_positions[right]
        dtw_square[np.ix_(left_positions, right_positions)] = value
        dtw_square[np.ix_(right_positions, left_positions)] = value
    save_dtw_distance_cache(
        cache / "dtw_distances_condensed.npz",
        call_ids,
        square_to_condensed(dtw_square),
    )

    artifacts = {
        "processed_calls": processed / "calls.csv",
        "call_session_inventory": cache / "call_session_inventory.csv",
        "call_session_pairs": cache / "call_session_pairs.csv",
        "stp_pca_scores": cache / "stp_pca_scores.csv",
        "mfcc_pca_scores": cache / "mfcc_pca_scores.csv",
        "dtw_distances_condensed": cache / "dtw_distances_condensed.npz",
        "vae_latent_means": cache / "vae_latent_means.csv",
    }
    manifest = CacheManifest.from_files(
        root,
        artifacts,
        parameters={
            "processed_calls": {"schema": "canonical-calls-1.0", "rows": 20},
            "call_session_inventory": {
                "minimum_calls": 5,
                "rows": 4,
                "eligible_units": 20,
            },
            "call_session_pairs": {
                "comparison_rule": (
                    "all_cross_individual_sessions_within_pair_stage_context"
                ),
                "rows": 2,
            },
            "stp_pca_scores": {
                "components": 5,
                "pca_refit": False,
                "canonical_unique_calls": 20,
            },
            "mfcc_pca_scores": {
                "components": 5,
                "pca_refit": False,
                "canonical_unique_calls": 20,
            },
            "dtw_distances_condensed": DTW_PARAMETERS,
            "vae_latent_means": VAE_PARAMETERS,
        },
        ordered_ids={
            "processed_calls": call_ids,
            "call_session_inventory": repertoire_ids,
            "call_session_pairs": tuple(comparisons["comparison_id"]),
            "stp_pca_scores": call_ids,
            "mfcc_pca_scores": call_ids,
            "dtw_distances_condensed": call_ids,
            "vae_latent_means": call_ids,
        },
    )
    manifest.write(cache / "manifest.json")


def test_cache_only_assembly_writes_four_metric_tables_and_registers_them(
    tmp_path: Path,
):
    _synthetic_cache_tree(tmp_path)

    result = ASSEMBLER.assemble(tmp_path)

    assert result["call_comparisons"] == 2
    assert result["call_distance_rows"] == 8
    assert result["call_model_rows"] == 8
    distances = pd.read_csv(
        tmp_path / "data" / "derived" / "call_repertoire_distances.csv"
    )
    assert tuple(distances["metric"].unique()) == ("stp", "mfcc", "dtw", "vae")
    assert set(distances["n_lower_level_pairs"]) == {25}
    observed = distances.set_index(["comparison_id", "metric"])["distance"]
    assert observed[("before_pair", "stp")] == pytest.approx(1.0)
    assert observed[("after_pair", "stp")] == pytest.approx(3.0)
    assert observed[("before_pair", "mfcc")] == pytest.approx(2.0)
    assert observed[("after_pair", "mfcc")] == pytest.approx(5.0)
    assert observed[("before_pair", "dtw")] == pytest.approx(2.0)
    assert observed[("after_pair", "dtw")] == pytest.approx(4.0)
    assert observed[("before_pair", "vae")] == pytest.approx(1.0)
    assert observed[("after_pair", "vae")] == pytest.approx(2.0)

    model = pd.read_csv(tmp_path / "data" / "derived" / "model_call.csv")
    assert set(model["stage_centered"]) == {-0.5, 0.5}
    for _, values in model.groupby("metric")["standardized_distance"]:
        np.testing.assert_allclose(sorted(values), [-1.0, 1.0])

    manifest = CacheManifest.load(tmp_path / "data" / "cache" / "manifest.json")
    manifest.validate(
        tmp_path,
        required_artifacts=["call_repertoire_distances", "model_call"],
    )
    assert manifest.artifacts["model_call"].parameters[
        "standard_deviation_ddof"
    ] == 0
    assert manifest.metadata["derived_call_tables_cache_only"] is True


def test_assembly_stops_before_writing_when_a_required_cache_is_missing(
    tmp_path: Path,
):
    _synthetic_cache_tree(tmp_path)
    (tmp_path / "data" / "cache" / "vae_latent_means.csv").unlink()

    with pytest.raises(CacheValidationError, match="missing files:.*vae_latent_means"):
        ASSEMBLER.assemble(tmp_path)

    assert not (tmp_path / "data" / "derived" / "model_call.csv").exists()


def test_assembly_rejects_a_tampered_cache_before_writing(tmp_path: Path):
    _synthetic_cache_tree(tmp_path)
    path = tmp_path / "data" / "cache" / "stp_pca_scores.csv"
    path.write_text(path.read_text(encoding="utf-8") + "\n", encoding="utf-8")

    with pytest.raises(CacheValidationError, match="size is|SHA-256 mismatch"):
        ASSEMBLER.assemble(tmp_path)

    assert not (tmp_path / "data" / "derived" / "model_call.csv").exists()


def test_assembly_rejects_an_unregistered_cache_before_writing(tmp_path: Path):
    _synthetic_cache_tree(tmp_path)
    manifest_path = tmp_path / "data" / "cache" / "manifest.json"
    manifest = CacheManifest.load(manifest_path)
    CacheManifest(
        artifacts={
            name: record
            for name, record in manifest.artifacts.items()
            if name != "vae_latent_means"
        },
        metadata=manifest.metadata,
    ).write(manifest_path)

    with pytest.raises(CacheValidationError, match="not registered.*vae_latent_means"):
        ASSEMBLER.assemble(tmp_path)

    assert not (tmp_path / "data" / "derived" / "model_call.csv").exists()
