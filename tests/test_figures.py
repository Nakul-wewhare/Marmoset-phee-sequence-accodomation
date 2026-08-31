from pathlib import Path
import importlib.util

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from marmoset_convergence.exceptions import DataValidationError
from marmoset_convergence.figures import (
    EMBEDDING_FIGURES,
    plot_embedding_space,
    prepare_embedding_table,
    regenerate_sequence_summary,
    regenerated_figure_path,
)


def _embedding_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    position = 0
    for context in ("partner", "non-partner"):
        for stage in ("before", "after"):
            for individual, offset in (("Tabor", 0.0), ("Lola", 0.7)):
                for replicate in range(2):
                    rows.append(
                        {
                            "repertoire_id": f"r{position:02d}",
                            "embedding_x": position / 3,
                            "embedding_y": offset + replicate + (stage == "after"),
                            "individual_id": individual,
                            "pair_id": "A",
                            "stage": stage,
                            "context": context,
                        }
                    )
                    position += 1
    embedding = pd.DataFrame(rows)
    metadata = embedding[
        ["repertoire_id", "individual_id", "pair_id", "stage", "context"]
    ].copy()
    return embedding, metadata


def test_embedding_contract_requires_exact_order_schema_and_finite_coordinates():
    embedding, metadata = _embedding_inputs()
    table, x_name, y_name = prepare_embedding_table(
        embedding, metadata, id_column="repertoire_id"
    )
    assert (x_name, y_name) == ("embedding_x", "embedding_y")
    assert tuple(table["repertoire_id"]) == tuple(metadata["repertoire_id"])

    with pytest.raises(DataValidationError, match="order differs"):
        prepare_embedding_table(
            embedding.iloc[::-1].reset_index(drop=True),
            metadata,
            id_column="repertoire_id",
        )
    with pytest.raises(DataValidationError, match="seven-column schema"):
        prepare_embedding_table(
            embedding.assign(extra_coordinate=0.0),
            metadata,
            id_column="repertoire_id",
        )
    invalid = embedding.copy()
    invalid.loc[0, "embedding_x"] = np.inf
    with pytest.raises(DataValidationError, match="finite numeric"):
        prepare_embedding_table(invalid, metadata, id_column="repertoire_id")


def test_embedding_panels_reuse_pooled_axis_limits():
    embedding, metadata = _embedding_inputs()
    figure = plot_embedding_space(embedding, metadata, EMBEDDING_FIGURES["4"])
    try:
        assert len(figure.axes) == 5
        x_limits = [axis.get_xlim() for axis in figure.axes]
        y_limits = [axis.get_ylim() for axis in figure.axes]
        assert all(limits == x_limits[0] for limits in x_limits[1:])
        assert all(limits == y_limits[0] for limits in y_limits[1:])
    finally:
        plt.close(figure)


def test_regenerated_export_never_replaces_locked_snapshot(tmp_path: Path):
    locked = tmp_path / "results" / "figures" / "main" / "figure_2.png"
    locked.parent.mkdir(parents=True)
    locked.write_bytes(b"locked manuscript snapshot")
    sequences = pd.DataFrame({"sequence": ["AA", "APK", "RRT", "SAP"]})

    output = regenerate_sequence_summary(sequences, locked)

    assert output == locked.with_name("figure_2_regenerated.png")
    assert output.is_file()
    assert locked.read_bytes() == b"locked manuscript snapshot"
    with pytest.raises(FileExistsError, match="already exists"):
        regenerate_sequence_summary(sequences, locked)
    regenerate_sequence_summary(sequences, locked, overwrite_regenerated=True)
    assert locked.read_bytes() == b"locked manuscript snapshot"


def test_regenerated_path_and_missing_audio_annotation_are_explicit(tmp_path: Path):
    locked = tmp_path / "figure_s5_mfcc.png"
    assert regenerated_figure_path(locked).name == "figure_s5_mfcc_regenerated.png"
    with pytest.raises(DataValidationError, match="locked manuscript path"):
        regenerated_figure_path(tmp_path / "figure_s5_mfcc_regenerated.png")
    note = EMBEDDING_FIGURES["s5"].annotation_note.casefold()
    assert "spectrogram cards" in note
    assert "complete wav collection" in note


def test_figure_script_checks_all_destination_conflicts_before_manifest_io(tmp_path: Path):
    script_path = Path(__file__).resolve().parents[1] / "scripts" / "regenerate_cached_figures.py"
    module_spec = importlib.util.spec_from_file_location(
        "regenerate_cached_figures_test_module", script_path
    )
    assert module_spec is not None and module_spec.loader is not None
    module = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)

    destination = (
        tmp_path / "results" / "figures" / "main" / "figure_2_regenerated.png"
    )
    destination.parent.mkdir(parents=True)
    destination.write_bytes(b"existing reader render")
    with pytest.raises(FileExistsError, match="already exists"):
        module.regenerate_cached_figures(tmp_path, ["2"])
