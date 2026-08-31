from pathlib import Path

from marmoset_convergence import ProjectPaths


def test_project_paths_match_public_layout(tmp_path: Path):
    paths = ProjectPaths.from_root(tmp_path)

    assert paths.sequences == tmp_path / "data" / "processed" / "sequences.csv"
    assert paths.calls == tmp_path / "data" / "processed" / "calls.csv"
    assert paths.sequence_distances == tmp_path / "data" / "cache" / "sequence_distances.npz"
    assert paths.dtw_distances == tmp_path / "data" / "cache" / "dtw_distances_condensed.npz"
    assert paths.model_sequence == tmp_path / "data" / "derived" / "model_sequence.csv"
    assert paths.model_call == tmp_path / "data" / "derived" / "model_call.csv"
    assert paths.audio == tmp_path / "extracted calls"
