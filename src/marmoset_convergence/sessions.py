"""Construction of eligible receiver-specific session repertoires."""

from __future__ import annotations

from itertools import product
from typing import Any, Mapping, Optional

import pandas as pd

from .exceptions import DataValidationError
from .identifiers import (
    BondedPairs,
    canonical_comparison_id,
    canonical_session_id,
    normalize_bonded_pairs,
    normalize_stage,
)

_CONTEXT_ALIASES = {
    "partner": "partner",
    "paired": "partner",
    "non-partner": "non-partner",
    "non_partner": "non-partner",
    "nonpartner": "non-partner",
    "stranger": "non-partner",
}


def normalize_context(value: Any) -> str:
    if value is None:
        raise DataValidationError("context cannot be missing")
    raw = str(value).strip().casefold()
    try:
        return _CONTEXT_ALIASES[raw]
    except KeyError as exc:
        raise DataValidationError(
            f"Unknown context {value!r}; expected partner or non-partner"
        ) from exc


def _individual_to_pair(
    bonded_pairs: BondedPairs,
) -> tuple[Mapping[str, str], Mapping[str, str], Mapping[str, tuple[str, str]]]:
    normalized = normalize_bonded_pairs(bonded_pairs)
    pair_by_individual = {}
    mate_by_individual = {}
    for pair_id, (a, b) in normalized.items():
        pair_by_individual[a] = pair_id
        pair_by_individual[b] = pair_id
        mate_by_individual[a] = b
        mate_by_individual[b] = a
    return pair_by_individual, mate_by_individual, normalized


def add_repertoire_ids(
    table: pd.DataFrame,
    *,
    individual_col: str = "individual_id",
    receiver_col: str = "receiver_id",
    stage_col: str = "stage",
    session_col: str = "session_number",
    output_col: str = "repertoire_id",
) -> pd.DataFrame:
    """Add canonical IDs containing individual, receiver, stage, and session."""

    required = (individual_col, receiver_col, stage_col, session_col)
    missing = [column for column in required if column not in table]
    if missing:
        raise DataValidationError("Missing repertoire key columns: " + ", ".join(missing))
    result = table.copy()
    result[stage_col] = result[stage_col].map(normalize_stage)
    identifiers = [
        canonical_session_id(individual, receiver, stage, session)
        for individual, receiver, stage, session in result[list(required)].itertuples(
            index=False, name=None
        )
    ]
    if output_col in result:
        existing = result[output_col].astype(str)
        mismatch = existing != pd.Series(identifiers, index=result.index)
        if mismatch.any():
            raise DataValidationError(
                f"Existing {output_col!r} disagrees with canonical session keys"
            )
        result[output_col] = identifiers
    else:
        result.insert(0, output_col, identifiers)
    return result


def build_session_inventory(
    units: pd.DataFrame,
    *,
    unit_id_col: str,
    count_col: str,
    minimum_units: int = 5,
    bonded_pairs: Optional[BondedPairs] = None,
    individual_col: str = "individual_id",
    receiver_col: str = "receiver_id",
    stage_col: str = "stage",
    session_col: str = "session_number",
    context_col: str = "context",
    pair_col: str = "pair_id",
) -> pd.DataFrame:
    """Create the eligible sequence or call session inventory.

    ``unit_id_col`` should be ``sequence_id`` or ``call_id`` and must already be
    unique.  Eligibility is evaluated after grouping by the full receiver-
    specific session key.
    """

    if not isinstance(minimum_units, int) or minimum_units < 1:
        raise DataValidationError("minimum_units must be a positive integer")
    required = [unit_id_col, individual_col, receiver_col, stage_col, session_col]
    missing = [column for column in required if column not in units]
    if missing:
        raise DataValidationError("Missing unit-table columns: " + ", ".join(missing))
    if units[unit_id_col].isna().any() or units[unit_id_col].astype(str).str.strip().eq("").any():
        raise DataValidationError(f"{unit_id_col} cannot contain missing/empty IDs")
    if units[unit_id_col].duplicated().any():
        duplicated = units.loc[
            units[unit_id_col].duplicated(keep=False), unit_id_col
        ].astype(str)
        raise DataValidationError(
            f"{unit_id_col} must be unique before inventory construction; duplicate "
            + ", ".join(sorted(duplicated.unique())[:5])
        )

    work = units.copy()
    work[stage_col] = work[stage_col].map(normalize_stage)

    pair_by_individual = mate_by_individual = None
    if bonded_pairs is not None:
        pair_by_individual, mate_by_individual, _ = _individual_to_pair(bonded_pairs)
    if pair_col not in work:
        if pair_by_individual is None:
            raise DataValidationError(
                f"{pair_col!r} is absent; provide bonded_pairs to derive it"
            )
        work[pair_col] = work[individual_col].map(pair_by_individual)
    if work[pair_col].isna().any():
        unknown = sorted(work.loc[work[pair_col].isna(), individual_col].astype(str).unique())
        raise DataValidationError("Individuals absent from bonded_pairs: " + ", ".join(unknown))

    if context_col not in work:
        if mate_by_individual is None:
            raise DataValidationError(
                f"{context_col!r} is absent; provide bonded_pairs to derive it"
            )
        work[context_col] = [
            "partner" if str(receiver) == mate_by_individual.get(str(individual)) else "non-partner"
            for individual, receiver in work[[individual_col, receiver_col]].itertuples(
                index=False, name=None
            )
        ]
    work[context_col] = work[context_col].map(normalize_context)
    work = add_repertoire_ids(
        work,
        individual_col=individual_col,
        receiver_col=receiver_col,
        stage_col=stage_col,
        session_col=session_col,
    )

    group_columns = [
        "repertoire_id",
        individual_col,
        receiver_col,
        pair_col,
        stage_col,
        context_col,
        session_col,
    ]
    inventory = (
        work.groupby(group_columns, sort=False, dropna=False)[unit_id_col]
        .size()
        .rename(count_col)
        .reset_index()
    )
    if inventory["repertoire_id"].duplicated().any():
        raise DataValidationError(
            "A canonical repertoire_id maps to inconsistent pair/context metadata"
        )
    inventory = inventory.loc[inventory[count_col] >= minimum_units].copy()
    inventory = inventory.rename(
        columns={
            individual_col: "individual_id",
            receiver_col: "receiver_id",
            pair_col: "pair_id",
            stage_col: "stage",
            context_col: "context",
            session_col: "session_number",
        }
    )
    columns = [
        "repertoire_id",
        "individual_id",
        "receiver_id",
        "pair_id",
        "stage",
        "context",
        "session_number",
        count_col,
    ]
    return inventory[columns].sort_values(
        ["pair_id", "stage", "context", "individual_id", "repertoire_id"],
        kind="stable",
    ).reset_index(drop=True)


def build_session_pairs(
    inventory: pd.DataFrame,
    bonded_pairs: BondedPairs,
    *,
    validate_context_receivers: bool = True,
) -> pd.DataFrame:
    """Generate all cross-individual repertoire pairs within stage and context."""

    required = {
        "repertoire_id",
        "individual_id",
        "receiver_id",
        "pair_id",
        "stage",
        "context",
        "session_number",
    }
    missing = sorted(required.difference(inventory.columns))
    if missing:
        raise DataValidationError("Missing inventory columns: " + ", ".join(missing))
    if inventory["repertoire_id"].duplicated().any():
        raise DataValidationError("Inventory repertoire_id values must be unique")

    pair_by_individual, mate_by_individual, normalized_pairs = _individual_to_pair(
        bonded_pairs
    )
    work = inventory.copy()
    work["stage"] = work["stage"].map(normalize_stage)
    work["context"] = work["context"].map(normalize_context)

    expected_pair = work["individual_id"].astype(str).map(pair_by_individual)
    if expected_pair.isna().any():
        unknown = sorted(work.loc[expected_pair.isna(), "individual_id"].astype(str).unique())
        raise DataValidationError("Inventory contains unknown individuals: " + ", ".join(unknown))
    if not expected_pair.equals(work["pair_id"].astype(str)):
        raise DataValidationError("Inventory pair_id is inconsistent with bonded_pairs")

    if validate_context_receivers:
        for row in work.itertuples(index=False):
            mate = mate_by_individual[str(row.individual_id)]
            is_partner_receiver = str(row.receiver_id) == mate
            if (row.context == "partner") != is_partner_receiver:
                raise DataValidationError(
                    f"Context/receiver mismatch for repertoire {row.repertoire_id!r}"
                )

    rows = []
    for pair_id, (individual_a, individual_b) in normalized_pairs.items():
        pair_rows = work.loc[work["pair_id"].astype(str) == pair_id]
        conditions = pair_rows[["stage", "context"]].drop_duplicates()
        for stage, context in conditions.itertuples(index=False, name=None):
            subset = pair_rows.loc[
                (pair_rows["stage"] == stage) & (pair_rows["context"] == context)
            ]
            a_rows = subset.loc[subset["individual_id"].astype(str) == individual_a]
            b_rows = subset.loc[subset["individual_id"].astype(str) == individual_b]
            for (_, a), (_, b) in product(a_rows.iterrows(), b_rows.iterrows()):
                comparison_id = canonical_comparison_id(
                    pair_id, stage, context, a["repertoire_id"], b["repertoire_id"]
                )
                rows.append(
                    {
                        "comparison_id": comparison_id,
                        "pair_id": pair_id,
                        "stage": stage,
                        "context": context,
                        "repertoire_a_id": a["repertoire_id"],
                        "repertoire_b_id": b["repertoire_id"],
                        "individual_a": individual_a,
                        "individual_b": individual_b,
                        "receiver_a": a["receiver_id"],
                        "receiver_b": b["receiver_id"],
                        "session_a": a["session_number"],
                        "session_b": b["session_number"],
                    }
                )
    columns = [
        "comparison_id",
        "pair_id",
        "stage",
        "context",
        "repertoire_a_id",
        "repertoire_b_id",
        "individual_a",
        "individual_b",
        "receiver_a",
        "receiver_b",
        "session_a",
        "session_b",
    ]
    result = pd.DataFrame(rows, columns=columns)
    if result["comparison_id"].duplicated().any():
        raise DataValidationError("Generated comparison_id values are not unique")
    return result.sort_values(
        ["pair_id", "stage", "context", "repertoire_a_id", "repertoire_b_id"],
        kind="stable",
    ).reset_index(drop=True)
