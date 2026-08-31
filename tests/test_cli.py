from pathlib import Path

import numpy as np
import pandas as pd

from marmoset_convergence import CacheManifest, DTW_PARAMETERS
from marmoset_convergence.cli import main, run_audio_preflight, validate_cache
from marmoset_convergence.matrices import save_dtw_distance_cache


def test_validate_cache_cli_checks_manifest_parameters_ids_and_npz(tmp_path: Path):
    cache_dir = tmp_path / "data" / "cache"
    cache_dir.mkdir(parents=True)
    artifact = cache_dir / "dtw_distances_condensed.npz"
    ids = ["a.wav", "b.wav"]
    save_dtw_distance_cache(artifact, ids, np.array([0.25]))
    manifest = CacheManifest.from_files(
        tmp_path,
        {"dtw_distances_condensed": artifact},
        parameters={"dtw_distances_condensed": DTW_PARAMETERS},
        ordered_ids={"dtw_distances_condensed": ids},
    )
    manifest.write(cache_dir / "manifest.json")

    result = validate_cache(tmp_path)

    assert result == {
        "status": "ok",
        "validated_artifacts": ["dtw_distances_condensed"],
    }
    assert main(["validate-cache", "--root", str(tmp_path)]) == 0


def test_audio_preflight_cli_defaults_to_extracted_calls(tmp_path: Path):
    processed = tmp_path / "data" / "processed"
    processed.mkdir(parents=True)
    pd.DataFrame(
        {"call_id": ["a.wav"], "filename": ["a.wav"]}
    ).to_csv(processed / "calls.csv", index=False)
    audio = tmp_path / "extracted calls"
    audio.mkdir()
    (audio / "a.wav").write_bytes(b"header not decoded")

    result = run_audio_preflight(tmp_path)

    assert result["status"] == "ok"
    assert result["expected_wavs"] == 1
    assert main(["preflight-audio", "--root", str(tmp_path)]) == 0
