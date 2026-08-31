import pandas as pd
import pytest

from marmoset_convergence import (
    DataValidationError,
    SessionKey,
    build_session_inventory,
    build_session_pairs,
    canonical_session_id,
    canonical_comparison_id,
    normalize_context,
)


def test_session_id_contains_every_identity_component_and_round_trips():
    identifier = canonical_session_id("Odin | focal", "Nougatti/a", "pre", 1.0)

    assert identifier == (
        "focal=Odin%20%7C%20focal|conspecific=Nougatti%2Fa|stage=before|session=1"
    )
    parsed = SessionKey.from_id(identifier)
    assert parsed.focal == "Odin | focal"
    assert parsed.conspecific == "Nougatti/a"
    assert parsed.stage == "before"
    assert parsed.session == "1"


def test_receiver_and_stage_prevent_session_id_collisions():
    base = canonical_session_id("A", "B", "before", 1)
    assert base != canonical_session_id("A", "C", "before", 1)
    assert base != canonical_session_id("A", "B", "after", 1)


def test_comparison_id_is_unordered_in_its_two_repertoires():
    forward = canonical_comparison_id("Pair A", "before", "partner", "r1", "r2")
    reverse = canonical_comparison_id("Pair A", "before", "partner", "r2", "r1")
    assert forward == reverse


def test_context_aliases_emit_public_hyphenated_label():
    assert normalize_context("non_partner") == "non-partner"
    assert normalize_context("stranger") == "non-partner"
    assert normalize_context("partner") == "partner"


def _unit_table():
    rows = []
    for unit in range(5):
        rows.append(
            {
                "sequence_id": f"a-partner-{unit}",
                "individual_id": "A",
                "receiver_id": "B",
                "stage": "pre",
                "session_number": 1.0,
            }
        )
        rows.append(
            {
                "sequence_id": f"b-partner-{unit}",
                "individual_id": "B",
                "receiver_id": "A",
                "stage": "before",
                "session_number": 1,
            }
        )
    for unit in range(4):
        rows.append(
            {
                "sequence_id": f"ineligible-{unit}",
                "individual_id": "A",
                "receiver_id": "C",
                "stage": "before",
                "session_number": 2,
            }
        )
    return pd.DataFrame(rows)


def test_inventory_derives_pair_context_and_applies_minimum_per_full_session_key():
    inventory = build_session_inventory(
        _unit_table(),
        unit_id_col="sequence_id",
        count_col="n_sequences",
        minimum_units=5,
        bonded_pairs={"Pair A": ("A", "B")},
    )

    assert len(inventory) == 2
    assert set(inventory["context"]) == {"partner"}
    assert set(inventory["n_sequences"]) == {5}
    assert all("conspecific=" in value for value in inventory["repertoire_id"])


def test_inventory_refuses_duplicate_unit_ids_before_counting():
    units = _unit_table()
    units.loc[1, "sequence_id"] = units.loc[0, "sequence_id"]
    with pytest.raises(DataValidationError, match="must be unique"):
        build_session_inventory(
            units,
            unit_id_col="sequence_id",
            count_col="n_sequences",
            bonded_pairs={"Pair A": ("A", "B")},
        )


def test_session_pairs_use_cartesian_products_within_stage_and_context():
    records = []
    for individual, mate, other, sessions in (
        ("A", "B", "C", (1, 2)),
        ("B", "A", "D", (3, 4, 5)),
    ):
        for session in sessions:
            repertoire_id = canonical_session_id(
                individual, other, "after", session
            )
            records.append(
                {
                    "repertoire_id": repertoire_id,
                    "individual_id": individual,
                    "receiver_id": other,
                    "pair_id": "Pair A",
                    "stage": "after",
                    "context": "non_partner",
                    "session_number": session,
                    "n_sequences": 5,
                }
            )
    inventory = pd.DataFrame(records)

    comparisons = build_session_pairs(
        inventory, {"Pair A": ("A", "B")}
    )

    assert len(comparisons) == 2 * 3
    assert set(comparisons["context"]) == {"non-partner"}
    assert comparisons["comparison_id"].is_unique
    assert set(comparisons["receiver_a"]) == {"C"}
    assert set(comparisons["receiver_b"]) == {"D"}


def test_session_pairs_validate_partner_receiver_identity():
    inventory = pd.DataFrame(
        [
            {
                "repertoire_id": canonical_session_id("A", "C", "before", 1),
                "individual_id": "A",
                "receiver_id": "C",
                "pair_id": "Pair A",
                "stage": "before",
                "context": "partner",
                "session_number": 1,
            },
            {
                "repertoire_id": canonical_session_id("B", "A", "before", 1),
                "individual_id": "B",
                "receiver_id": "A",
                "pair_id": "Pair A",
                "stage": "before",
                "context": "partner",
                "session_number": 1,
            },
        ]
    )
    with pytest.raises(DataValidationError, match="Context/receiver mismatch"):
        build_session_pairs(inventory, {"Pair A": ("A", "B")})
