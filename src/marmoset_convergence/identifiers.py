"""Stable identifiers for session repertoires and their comparisons."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from typing import Any, Mapping, Sequence, Tuple, Union
from urllib.parse import quote, unquote

from .exceptions import DataValidationError

_INTEGER_LIKE = re.compile(r"^[+-]?\d+(?:\.0+)?$")
_STAGE_ALIASES = {
    "before": "before",
    "pre": "before",
    "pre-pairing": "before",
    "pre_pairing": "before",
    "prepairing": "before",
    "after": "after",
    "post": "after",
    "post-pairing": "after",
    "post_pairing": "after",
    "postpairing": "after",
}


def _nonempty_text(value: Any, field: str) -> str:
    if value is None:
        raise DataValidationError(f"{field} cannot be missing")
    text = str(value).strip()
    if not text or text.casefold() in {"nan", "none", "<na>"}:
        raise DataValidationError(f"{field} cannot be empty")
    return text


def normalize_stage(stage: Any) -> str:
    """Return canonical ``before``/``after`` stage labels."""

    raw = _nonempty_text(stage, "stage").casefold()
    try:
        return _STAGE_ALIASES[raw]
    except KeyError as exc:
        raise DataValidationError(
            f"Unknown stage {stage!r}; expected before/pre or after/post"
        ) from exc


def normalize_session_number(session: Any) -> str:
    """Normalize equivalent numeric labels (for example, 1 and 1.0)."""

    value = _nonempty_text(session, "session")
    if _INTEGER_LIKE.fullmatch(value):
        return str(int(float(value)))
    return value


def _encode(value: str) -> str:
    return quote(value, safe="-._~")


@dataclass(frozen=True, order=True)
class SessionKey:
    """Identity of one individual's repertoire in one recording session.

    The receiver/conspecific is part of the identity.  Omitting it conflates
    different audience-specific repertoires when session numbers are reused.
    """

    focal: str
    conspecific: str
    stage: str
    session: str

    @classmethod
    def create(
        cls, focal: Any, conspecific: Any, stage: Any, session: Any
    ) -> "SessionKey":
        return cls(
            focal=_nonempty_text(focal, "focal"),
            conspecific=_nonempty_text(conspecific, "conspecific"),
            stage=normalize_stage(stage),
            session=normalize_session_number(session),
        )

    @property
    def id(self) -> str:
        return (
            f"focal={_encode(self.focal)}|"
            f"conspecific={_encode(self.conspecific)}|"
            f"stage={_encode(self.stage)}|"
            f"session={_encode(self.session)}"
        )

    @classmethod
    def from_id(cls, identifier: str) -> "SessionKey":
        fields = {}
        for component in _nonempty_text(identifier, "repertoire_id").split("|"):
            if "=" not in component:
                raise DataValidationError(f"Malformed repertoire_id: {identifier!r}")
            key, value = component.split("=", 1)
            if key in fields:
                raise DataValidationError(f"Repeated field {key!r} in repertoire_id")
            fields[key] = unquote(value)
        expected = {"focal", "conspecific", "stage", "session"}
        if set(fields) != expected:
            raise DataValidationError(
                f"Malformed repertoire_id fields: expected {sorted(expected)}, "
                f"found {sorted(fields)}"
            )
        parsed = cls.create(**fields)
        if parsed.id != identifier:
            raise DataValidationError(
                f"repertoire_id is not canonical; expected {parsed.id!r}"
            )
        return parsed


def canonical_session_id(
    focal: Any, conspecific: Any, stage: Any, session: Any
) -> str:
    """Construct a collision-safe, human-readable repertoire identifier."""

    return SessionKey.create(focal, conspecific, stage, session).id


def canonical_pair_id(individual_a: Any, individual_b: Any) -> str:
    """Construct an order-independent bonded-pair identifier."""

    a = _nonempty_text(individual_a, "individual_a")
    b = _nonempty_text(individual_b, "individual_b")
    if a == b:
        raise DataValidationError("A bonded pair must contain two individuals")
    first, second = sorted((a, b), key=lambda value: (value.casefold(), value))
    return f"pair={_encode(first)}|{_encode(second)}"


BondedPairs = Union[
    Mapping[str, Sequence[Any]],
    Sequence[Tuple[Any, Any]],
]


def normalize_bonded_pairs(bonded_pairs: BondedPairs) -> Mapping[str, Tuple[str, str]]:
    """Normalize either ``{pair_id: (a, b)}`` or a sequence of pairs."""

    normalized = {}
    items = bonded_pairs.items() if isinstance(bonded_pairs, Mapping) else (
        (canonical_pair_id(pair[0], pair[1]), pair) for pair in bonded_pairs
    )
    seen_individuals = set()
    for raw_pair_id, pair in items:
        if len(pair) != 2:
            raise DataValidationError(f"Pair {raw_pair_id!r} does not have two members")
        a = _nonempty_text(pair[0], "individual_a")
        b = _nonempty_text(pair[1], "individual_b")
        if a == b:
            raise DataValidationError(f"Pair {raw_pair_id!r} repeats one individual")
        pair_id = _nonempty_text(raw_pair_id, "pair_id")
        if pair_id in normalized:
            raise DataValidationError(f"Duplicate pair_id {pair_id!r}")
        overlap = seen_individuals.intersection((a, b))
        if overlap:
            raise DataValidationError(
                "Each focal individual must belong to one bonded pair; repeated: "
                + ", ".join(sorted(overlap))
            )
        normalized[pair_id] = (a, b)
        seen_individuals.update((a, b))
    if not normalized:
        raise DataValidationError("At least one bonded pair is required")
    return normalized


def canonical_comparison_id(
    pair_id: Any,
    stage: Any,
    context: Any,
    repertoire_a_id: Any,
    repertoire_b_id: Any,
) -> str:
    """Return a short deterministic ID for one repertoire comparison."""

    repertoire_ids = sorted(
        (
            _nonempty_text(repertoire_a_id, "repertoire_a_id"),
            _nonempty_text(repertoire_b_id, "repertoire_b_id"),
        )
    )
    if repertoire_ids[0] == repertoire_ids[1]:
        raise DataValidationError("A comparison requires two distinct repertoires")
    payload = {
        "context": _nonempty_text(context, "context").casefold(),
        "pair_id": _nonempty_text(pair_id, "pair_id"),
        "repertoire_ids": repertoire_ids,
        "stage": normalize_stage(stage),
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return f"comparison={digest[:24]}"
