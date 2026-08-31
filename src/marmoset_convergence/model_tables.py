"""Repertoire-distance and model-table assembly from validated cached inputs."""

from __future__ import annotations

from typing import Mapping, Optional, Protocol, Sequence

import numpy as np
import pandas as pd

from .exceptions import DataValidationError
from .matrices import CondensedDistanceMatrix
from .sessions import normalize_context
from .identifiers import normalize_stage

SEQUENCE_METRICS = (
    "transition_probability",
    "bigram",
    "phee_repeat",
    "local_alignment",
)
CALL_METRICS = ("stp", "mfcc", "dtw", "vae")

COMPARISON_COLUMNS = (
    "comparison_id",
    "pair_id",
    "stage",
    "context",
    "repertoire_a_id",
    "repertoire_b_id",
)


class CrossSetDistance(Protocol):
    """Distance source that can average one Cartesian product of item IDs."""

    metric: str

    def mean_between(
        self, left_ids: Sequence[str], right_ids: Sequence[str]
    ) -> float: ...


def _validate_comparisons(comparisons: pd.DataFrame) -> pd.DataFrame:
    missing = [column for column in COMPARISON_COLUMNS if column not in comparisons]
    if missing:
        raise DataValidationError("Missing comparison columns: " + ", ".join(missing))
    if comparisons["comparison_id"].duplicated().any():
        raise DataValidationError("comparison_id values must be unique")
    result = comparisons.copy()
    result["stage"] = result["stage"].map(normalize_stage)
    result["context"] = result["context"].map(normalize_context)
    return result


def attach_repertoire_distances(
    comparisons: pd.DataFrame,
    matrices: Mapping[str, CondensedDistanceMatrix],
) -> pd.DataFrame:
    """Attach direct repertoire-to-repertoire distances in long form."""

    work = _validate_comparisons(comparisons)
    if not matrices:
        raise DataValidationError("At least one repertoire distance matrix is required")
    rows = []
    for comparison in work.itertuples(index=False):
        base = {column: getattr(comparison, column) for column in COMPARISON_COLUMNS}
        for metric, matrix in matrices.items():
            rows.append(
                {
                    **base,
                    "metric": str(metric),
                    "distance": matrix.distance(
                        comparison.repertoire_a_id, comparison.repertoire_b_id
                    ),
                }
            )
    return pd.DataFrame(rows, columns=[*COMPARISON_COLUMNS, "metric", "distance"])


def aggregate_item_distances(
    comparisons: pd.DataFrame,
    membership: pd.DataFrame,
    matrices: Mapping[str, CrossSetDistance],
    *,
    item_id_col: str = "call_id",
    repertoire_id_col: str = "repertoire_id",
) -> pd.DataFrame:
    """Average every cross-repertoire call distance for each session pair.

    This consumes precomputed item-level matrices.  It does not extract audio or
    calculate DTW, and is therefore safe in the default cached workflow.
    """

    work = _validate_comparisons(comparisons)
    for column in (item_id_col, repertoire_id_col):
        if column not in membership:
            raise DataValidationError(f"Missing membership column {column!r}")
        if membership[column].isna().any():
            raise DataValidationError(f"Membership column {column!r} cannot be missing")
    if membership[item_id_col].duplicated().any():
        raise DataValidationError(f"Each {item_id_col} must belong to one repertoire")
    if not matrices:
        raise DataValidationError("At least one item distance matrix is required")
    mismatched_metrics = [
        str(metric)
        for metric, matrix in matrices.items()
        if str(matrix.metric) != str(metric)
    ]
    if mismatched_metrics:
        raise DataValidationError(
            "Distance-source metric labels disagree with mapping keys: "
            + ", ".join(mismatched_metrics)
        )
    ids_by_repertoire = (
        membership.groupby(repertoire_id_col, sort=False)[item_id_col]
        .apply(lambda values: tuple(values.astype(str)))
        .to_dict()
    )
    rows = []
    for comparison in work.itertuples(index=False):
        try:
            ids_a = ids_by_repertoire[comparison.repertoire_a_id]
            ids_b = ids_by_repertoire[comparison.repertoire_b_id]
        except KeyError as exc:
            raise DataValidationError(
                f"No item membership for repertoire {exc.args[0]!r}"
            ) from exc
        base = {column: getattr(comparison, column) for column in COMPARISON_COLUMNS}
        for metric, matrix in matrices.items():
            distance = matrix.mean_between(ids_a, ids_b)
            if not np.isfinite(distance) or distance < 0:
                raise DataValidationError(
                    f"Metric {metric!r} produced an invalid repertoire distance"
                )
            rows.append(
                {
                    **base,
                    "metric": str(metric),
                    "distance": distance,
                    "n_lower_level_pairs": len(ids_a) * len(ids_b),
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            *COMPARISON_COLUMNS,
            "metric",
            "distance",
            "n_lower_level_pairs",
        ],
    )


def joint_zscore_within_metric(
    table: pd.DataFrame,
    *,
    metric_col: str = "metric",
    value_col: str = "distance",
    output_col: str = "standardized_distance",
    ddof: int = 0,
) -> pd.DataFrame:
    """Z-score jointly using the population SD across pairs/stages/contexts."""

    if metric_col not in table or value_col not in table:
        raise DataValidationError(f"Table requires {metric_col!r} and {value_col!r}")
    if ddof < 0:
        raise DataValidationError("ddof must be non-negative")
    result = table.copy()
    numeric = pd.to_numeric(result[value_col], errors="coerce")
    if numeric.isna().any() or not np.isfinite(numeric.to_numpy()).all():
        raise DataValidationError(f"{value_col} must contain finite numeric distances")
    result[value_col] = numeric
    grouped = result.groupby(metric_col, sort=False, dropna=False)[value_col]
    means = grouped.transform("mean")
    standard_deviations = grouped.transform(lambda values: values.std(ddof=ddof))
    invalid_metrics = sorted(
        result.loc[
            standard_deviations.isna() | standard_deviations.eq(0), metric_col
        ]
        .astype(str)
        .unique()
    )
    if invalid_metrics:
        raise DataValidationError(
            "Cannot standardize metrics with fewer than two varying values: "
            + ", ".join(invalid_metrics)
        )
    result[output_col] = (result[value_col] - means) / standard_deviations
    return result


def build_model_table(
    repertoire_distances: pd.DataFrame,
    *,
    expected_metrics: Optional[Sequence[str]] = None,
    ddof: int = 0,
) -> pd.DataFrame:
    """Validate a complete long table and add model-ready transformations."""

    required = [*COMPARISON_COLUMNS, "metric", "distance"]
    missing = [column for column in required if column not in repertoire_distances]
    if missing:
        raise DataValidationError("Missing distance-table columns: " + ", ".join(missing))
    retained = [
        *required,
        *(
            ["n_lower_level_pairs"]
            if "n_lower_level_pairs" in repertoire_distances
            else []
        ),
    ]
    result = repertoire_distances[retained].copy()
    result["stage"] = result["stage"].map(normalize_stage)
    result["context"] = result["context"].map(normalize_context)
    result["metric"] = result["metric"].astype(str)
    if result.duplicated(["comparison_id", "metric"]).any():
        raise DataValidationError("Each comparison_id/metric combination must be unique")
    if "n_lower_level_pairs" in result:
        lower_counts = pd.to_numeric(result["n_lower_level_pairs"], errors="coerce")
        if lower_counts.isna().any() or (lower_counts < 1).any() or not np.equal(
            lower_counts, np.floor(lower_counts)
        ).all():
            raise DataValidationError("n_lower_level_pairs must contain positive integers")
        result["n_lower_level_pairs"] = lower_counts.astype(int)

    expected = tuple(expected_metrics or sorted(result["metric"].unique()))
    if not expected or len(set(expected)) != len(expected):
        raise DataValidationError("expected_metrics must contain distinct metric codes")
    observed_global = set(result["metric"])
    if observed_global != set(expected):
        raise DataValidationError(
            f"Metric set is {sorted(observed_global)}, expected {sorted(expected)}"
        )
    counts = result.groupby("comparison_id", sort=False)["metric"].agg(
        lambda values: frozenset(values)
    )
    incomplete = counts[counts != frozenset(expected)]
    if not incomplete.empty:
        raise DataValidationError(
            f"{len(incomplete)} comparison(s) do not contain every expected metric"
        )

    result = joint_zscore_within_metric(result, ddof=ddof)
    result["stage_centered"] = result["stage"].map(
        {"before": -0.5, "after": 0.5}
    )
    if result["stage_centered"].isna().any():  # normalize_stage should prevent this.
        raise DataValidationError("Could not construct centered stage coding")
    order = {metric: position for position, metric in enumerate(expected)}
    result["_metric_order"] = result["metric"].map(order)
    result = result.sort_values(
        ["pair_id", "stage", "context", "comparison_id", "_metric_order"],
        kind="stable",
    ).drop(columns="_metric_order")
    return result.reset_index(drop=True)
