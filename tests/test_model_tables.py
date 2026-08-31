import numpy as np
import pandas as pd
import pytest

from marmoset_convergence import (
    CondensedDistanceMatrix,
    DataValidationError,
    aggregate_item_distances,
    attach_repertoire_distances,
    build_model_table,
    joint_zscore_within_metric,
)


def _comparisons():
    return pd.DataFrame(
        [
            {
                "comparison_id": "c1",
                "pair_id": "Pair A",
                "stage": "before",
                "context": "stranger",
                "repertoire_a_id": "r1",
                "repertoire_b_id": "r2",
            },
            {
                "comparison_id": "c2",
                "pair_id": "Pair A",
                "stage": "post",
                "context": "non_partner",
                "repertoire_a_id": "r1",
                "repertoire_b_id": "r3",
            },
        ]
    )


def test_attach_repertoire_distances_builds_canonical_long_table():
    matrix = CondensedDistanceMatrix(
        ("r1", "r2", "r3"), np.array([1.0, 2.0, 3.0]), "bigram"
    )
    result = attach_repertoire_distances(_comparisons(), {"bigram": matrix})

    assert result["distance"].tolist() == [1.0, 2.0]
    assert set(result["context"]) == {"non-partner"}


def test_aggregate_item_distances_means_all_cross_call_pairs():
    # Order: a1-a2, a1-b1, a1-b2, a2-b1, a2-b2, b1-b2.
    matrix = CondensedDistanceMatrix(
        ("a1", "a2", "b1", "b2"),
        np.array([9.0, 1.0, 2.0, 3.0, 4.0, 9.0]),
        "dtw",
    )
    membership = pd.DataFrame(
        {
            "call_id": ["a1", "a2", "b1", "b2"],
            "repertoire_id": ["r1", "r1", "r2", "r2"],
        }
    )
    result = aggregate_item_distances(
        _comparisons().iloc[[0]], membership, {"dtw": matrix}
    )

    assert result.loc[0, "distance"] == pytest.approx((1 + 2 + 3 + 4) / 4)
    assert result.loc[0, "n_lower_level_pairs"] == 4


def test_joint_zscore_defaults_to_population_sd_and_is_joint_across_contexts():
    table = pd.DataFrame(
        {
            "metric": ["x", "x", "y", "y"],
            "distance": [1.0, 3.0, 10.0, 14.0],
            "context": ["partner", "non-partner", "partner", "non-partner"],
        }
    )
    result = joint_zscore_within_metric(table)

    np.testing.assert_allclose(result["standardized_distance"], [-1, 1, -1, 1])


def test_model_table_requires_every_metric_for_every_comparison():
    rows = []
    for comparison in _comparisons().to_dict("records"):
        for metric, offset in (("bigram", 0.0), ("local_alignment", 1.0)):
            rows.append({**comparison, "metric": metric, "distance": len(rows) + offset})
    table = build_model_table(
        pd.DataFrame(rows), expected_metrics=("bigram", "local_alignment")
    )

    assert set(table["stage_centered"]) == {-0.5, 0.5}
    assert set(table["context"]) == {"non-partner"}
    assert table.groupby("metric")["standardized_distance"].mean().abs().max() < 1e-12

    incomplete = pd.DataFrame(rows[:-1])
    with pytest.raises(DataValidationError, match="do not contain every"):
        build_model_table(
            incomplete, expected_metrics=("bigram", "local_alignment")
        )


def test_joint_zscore_rejects_constant_metric():
    table = pd.DataFrame({"metric": ["x", "x"], "distance": [1.0, 1.0]})
    with pytest.raises(DataValidationError, match="fewer than two varying"):
        joint_zscore_within_metric(table)
