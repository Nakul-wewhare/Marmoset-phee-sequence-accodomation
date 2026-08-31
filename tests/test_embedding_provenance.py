from pathlib import Path

import pandas as pd
import pytest

from marmoset_convergence.cli import validate_cache
from marmoset_convergence.provenance import (
    EMBEDDING_PARAMETERS,
    CacheManifest,
    sha256_file,
)


def _write_embedding_fixture(root: Path, *, reverse: bool = False) -> None:
    processed = root / "data/processed"
    cache = root / "data/cache"
    embeddings = cache / "embeddings"
    processed.mkdir(parents=True)
    embeddings.mkdir(parents=True)

    calls = pd.DataFrame(
        {
            "call_id": ["c1", "c2"],
            "individual_id": ["Tabor", "Lola"],
            "pair_id": ["A", "A"],
            "stage": ["before", "after"],
            "context": ["non-partner", "partner"],
        }
    )
    calls.to_csv(processed / "calls.csv", index=False)
    scores = pd.DataFrame(
        {
            "call_id": ["c1", "c2"],
            **{f"stp_pc_{index}": [0.0, float(index)] for index in range(1, 6)},
        }
    )
    scores.to_csv(cache / "stp_pca_scores.csv", index=False)
    embedding = calls.copy()
    embedding.insert(1, "embedding_x", [0.0, 1.0])
    embedding.insert(2, "embedding_y", [1.0, 0.0])
    if reverse:
        embedding = embedding.iloc[::-1].reset_index(drop=True)
    embedding_path = embeddings / "stp.csv"
    embedding.to_csv(embedding_path, index=False)

    parameters = {
        **EMBEDDING_PARAMETERS["stp"],
        "hyperparameters": {"n_neighbors": 10, "MN_ratio": 0.5, "FP_ratio": 2.0},
        "source_artifact_sha256": sha256_file(cache / "stp_pca_scores.csv"),
        "input_shape": [2, 5],
        "coordinate_shape": [2, 2],
    }
    manifest = CacheManifest.from_files(
        root,
        {"embedding_stp": embedding_path},
        parameters={"embedding_stp": parameters},
        ordered_ids={"embedding_stp": embedding["call_id"].astype(str)},
    )
    manifest.write(cache / "manifest.json")


def test_embedding_contract_binds_order_metadata_source_and_parameters(tmp_path: Path):
    _write_embedding_fixture(tmp_path)

    result = validate_cache(tmp_path)

    assert result["validated_artifacts"] == ["embedding_stp"]


def test_embedding_contract_rejects_noncanonical_identifier_order(tmp_path: Path):
    _write_embedding_fixture(tmp_path, reverse=True)

    with pytest.raises(ValueError, match="identifier order"):
        validate_cache(tmp_path)
