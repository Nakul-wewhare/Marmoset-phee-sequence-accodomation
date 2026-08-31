"""Canonical sequence representations and distances from the manuscript Methods."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import log
from typing import Hashable, Mapping, Sequence, Tuple

import numpy as np

from .exceptions import DataValidationError, ExpensiveComputationDisabled

Token = Hashable
VocalSequence = Sequence[Token]
Repertoire = Sequence[VocalSequence]


def _tokens(sequence: VocalSequence) -> Tuple[Token, ...]:
    if isinstance(sequence, (bytes, bytearray)):
        raise DataValidationError("Sequences must contain tokens, not raw bytes")
    values = tuple(sequence)
    if not values:
        raise DataValidationError("Vocal sequences cannot be empty")
    return values


def _repertoire(repertoire: Repertoire) -> Tuple[Tuple[Token, ...], ...]:
    sequences = tuple(_tokens(sequence) for sequence in repertoire)
    if not sequences:
        raise DataValidationError("Session repertoires cannot be empty")
    return sequences


def discover_alphabet(repertoires: Sequence[Repertoire]) -> Tuple[Token, ...]:
    """Return the sorted union of tokens used by the supplied repertoires."""

    tokens = {token for repertoire in repertoires for sequence in repertoire for token in sequence}
    if not tokens:
        raise DataValidationError("Cannot discover an alphabet from empty repertoires")
    try:
        return tuple(sorted(tokens))
    except TypeError as exc:
        raise DataValidationError("Alphabet tokens must have a consistent sortable type") from exc


def _alphabet_index(alphabet: Sequence[Token]) -> tuple[Tuple[Token, ...], Mapping[Token, int]]:
    ordered = tuple(alphabet)
    if not ordered or len(set(ordered)) != len(ordered):
        raise DataValidationError("alphabet must contain distinct tokens")
    return ordered, {token: index for index, token in enumerate(ordered)}


def transition_probability_vector(
    repertoire: Repertoire, alphabet: Sequence[Token]
) -> np.ndarray:
    """Flatten the row-normalized first-order transition table."""

    sequences = _repertoire(repertoire)
    ordered, index = _alphabet_index(alphabet)
    counts = np.zeros((len(ordered), len(ordered)), dtype=np.float64)
    for sequence in sequences:
        for preceding, following in zip(sequence[:-1], sequence[1:]):
            try:
                counts[index[preceding], index[following]] += 1.0
            except KeyError as exc:
                raise DataValidationError(
                    f"Sequence token {exc.args[0]!r} is absent from the fixed alphabet"
                ) from exc
    row_totals = counts.sum(axis=1, keepdims=True)
    probabilities = np.divide(
        counts,
        row_totals,
        out=np.zeros_like(counts),
        where=row_totals != 0,
    )
    return probabilities.reshape(-1)


def bigram_frequency_vector(
    repertoire: Repertoire, alphabet: Sequence[Token]
) -> np.ndarray:
    """Return global bigram proportions in fixed alphabet-product order."""

    sequences = _repertoire(repertoire)
    ordered, index = _alphabet_index(alphabet)
    counts = np.zeros((len(ordered), len(ordered)), dtype=np.float64)
    for sequence in sequences:
        for preceding, following in zip(sequence[:-1], sequence[1:]):
            try:
                counts[index[preceding], index[following]] += 1.0
            except KeyError as exc:
                raise DataValidationError(
                    f"Sequence token {exc.args[0]!r} is absent from the fixed alphabet"
                ) from exc
    total = counts.sum()
    if total:
        counts /= total
    return counts.reshape(-1)


def phee_repeat_vector(
    repertoire: Repertoire,
    *,
    target: Token = "A",
    minimum_length: int = 2,
    maximum_length: int = 5,
) -> np.ndarray:
    """Return proportions of exact target-token run lengths 2 through 5.

    Runs longer than ``maximum_length`` are outside the prespecified feature set
    and are not truncated into its final bin, matching the original analysis.
    A repertoire without an eligible run is represented by an all-zero vector.
    """

    if minimum_length < 1 or maximum_length < minimum_length:
        raise DataValidationError("Invalid repeat-length range")
    sequences = _repertoire(repertoire)
    counts = np.zeros(maximum_length - minimum_length + 1, dtype=np.float64)
    for sequence in sequences:
        cursor = 0
        while cursor < len(sequence):
            end = cursor + 1
            while end < len(sequence) and sequence[end] == sequence[cursor]:
                end += 1
            run_length = end - cursor
            if (
                sequence[cursor] == target
                and minimum_length <= run_length <= maximum_length
            ):
                counts[run_length - minimum_length] += 1.0
            cursor = end
    total = counts.sum()
    return counts / total if total else counts


def euclidean_feature_distance(left: np.ndarray, right: np.ndarray) -> float:
    left_array = np.asarray(left, dtype=np.float64)
    right_array = np.asarray(right, dtype=np.float64)
    if left_array.shape != right_array.shape:
        raise DataValidationError(
            f"Feature shapes differ: {left_array.shape} versus {right_array.shape}"
        )
    if not np.isfinite(left_array).all() or not np.isfinite(right_array).all():
        raise DataValidationError("Feature vectors must be finite")
    return float(np.linalg.norm(left_array - right_array))


def transition_probability_distance(
    repertoire_a: Repertoire, repertoire_b: Repertoire, alphabet: Sequence[Token]
) -> float:
    return euclidean_feature_distance(
        transition_probability_vector(repertoire_a, alphabet),
        transition_probability_vector(repertoire_b, alphabet),
    )


def bigram_distance(
    repertoire_a: Repertoire, repertoire_b: Repertoire, alphabet: Sequence[Token]
) -> float:
    return euclidean_feature_distance(
        bigram_frequency_vector(repertoire_a, alphabet),
        bigram_frequency_vector(repertoire_b, alphabet),
    )


def phee_repeat_distance(
    repertoire_a: Repertoire,
    repertoire_b: Repertoire,
    *,
    target: Token = "A",
    minimum_length: int = 2,
    maximum_length: int = 5,
) -> float:
    return euclidean_feature_distance(
        phee_repeat_vector(
            repertoire_a,
            target=target,
            minimum_length=minimum_length,
            maximum_length=maximum_length,
        ),
        phee_repeat_vector(
            repertoire_b,
            target=target,
            minimum_length=minimum_length,
            maximum_length=maximum_length,
        ),
    )


@dataclass(frozen=True)
class SmithWatermanConfig:
    """Scoring parameters fixed by the current manuscript Methods."""

    match: float = 2.0
    mismatch: float = -1.0
    gap: float = -1.0

    def __post_init__(self) -> None:
        if self.match <= 0:
            raise DataValidationError("Smith-Waterman match score must be positive")
        if self.mismatch > 0 or self.gap > 0:
            raise DataValidationError("Mismatch and gap scores must be non-positive")


DEFAULT_SMITH_WATERMAN = SmithWatermanConfig()


def smith_waterman_score(
    sequence_a: VocalSequence,
    sequence_b: VocalSequence,
    *,
    config: SmithWatermanConfig = DEFAULT_SMITH_WATERMAN,
) -> float:
    """Compute the highest local-alignment score with one canonical algorithm."""

    a = _tokens(sequence_a)
    b = _tokens(sequence_b)
    previous = np.zeros(len(b) + 1, dtype=np.float64)
    best = 0.0
    for token_a in a:
        current = np.zeros(len(b) + 1, dtype=np.float64)
        for j, token_b in enumerate(b, start=1):
            substitution = config.match if token_a == token_b else config.mismatch
            current[j] = max(
                0.0,
                previous[j - 1] + substitution,
                previous[j] + config.gap,
                current[j - 1] + config.gap,
            )
            if current[j] > best:
                best = float(current[j])
        previous = current
    return best


def normalized_local_alignment_similarity(
    sequence_a: VocalSequence,
    sequence_b: VocalSequence,
    *,
    config: SmithWatermanConfig = DEFAULT_SMITH_WATERMAN,
) -> float:
    """Normalize alignment score by log of the two complete sequence lengths."""

    a = _tokens(sequence_a)
    b = _tokens(sequence_b)
    denominator = log(len(a) + len(b))
    if denominator <= 0:
        raise DataValidationError("Combined sequence length must exceed one")
    return smith_waterman_score(a, b, config=config) / denominator


def repertoire_local_alignment_similarity(
    repertoire_a: Repertoire,
    repertoire_b: Repertoire,
    *,
    config: SmithWatermanConfig = DEFAULT_SMITH_WATERMAN,
) -> float:
    """Mean normalized score across every cross-repertoire sequence pair."""

    a = _repertoire(repertoire_a)
    b = _repertoire(repertoire_b)
    total = sum(
        normalized_local_alignment_similarity(seq_a, seq_b, config=config)
        for seq_a in a
        for seq_b in b
    )
    return float(total / (len(a) * len(b)))


@dataclass(frozen=True)
class LocalAlignmentDistances:
    repertoire_ids: Tuple[str, ...]
    similarities: np.ndarray
    distances: np.ndarray
    global_max_similarity: float


def similarities_to_distances(
    similarities: np.ndarray,
    *,
    global_max_similarity: float,
) -> np.ndarray:
    """Apply one fixed global ``maximum - similarity`` transformation."""

    values = np.asarray(similarities, dtype=np.float64)
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
        raise DataValidationError("Similarity matrix must be square")
    if not np.isfinite(values).all() or not np.isfinite(global_max_similarity):
        raise DataValidationError("Similarities and their global reference must be finite")
    if not np.allclose(values, values.T, rtol=1e-10, atol=1e-12):
        raise DataValidationError("Similarity matrix must be symmetric")
    off_diagonal = values[~np.eye(len(values), dtype=bool)]
    if off_diagonal.size and off_diagonal.max() > global_max_similarity + 1e-12:
        raise DataValidationError(
            "global_max_similarity is below an observed repertoire-pair similarity"
        )
    distances = global_max_similarity - values
    distances[np.abs(distances) < 1e-12] = 0.0
    np.fill_diagonal(distances, 0.0)
    return distances


def local_alignment_distance_matrix(
    repertoires: Mapping[str, Repertoire],
    *,
    config: SmithWatermanConfig = DEFAULT_SMITH_WATERMAN,
    global_max_similarity: float | None = None,
    allow_expensive: bool = False,
) -> LocalAlignmentDistances:
    """Build the all-repertoire matrix only after explicit opt-in.

    The automatically selected reference is the maximum across distinct
    repertoire pairs in this complete eligible set.  Pass the saved global
    reference when transforming any later subset.
    """

    if not allow_expensive:
        raise ExpensiveComputationDisabled(
            "Batch local alignment is disabled by default; load "
            "data/cache/sequence_distances.npz or pass allow_expensive=True"
        )
    ids = tuple(str(identifier) for identifier in repertoires)
    if len(ids) < 2:
        raise DataValidationError("At least two repertoires are required")
    if len(set(ids)) != len(ids):
        raise DataValidationError("Repertoire IDs must be unique")
    normalized = [_repertoire(repertoires[identifier]) for identifier in repertoires]
    similarities = np.zeros((len(ids), len(ids)), dtype=np.float64)
    for i, j in combinations(range(len(ids)), 2):
        value = repertoire_local_alignment_similarity(
            normalized[i], normalized[j], config=config
        )
        similarities[i, j] = similarities[j, i] = value
    observed_max = float(similarities[np.triu_indices(len(ids), k=1)].max())
    reference = observed_max if global_max_similarity is None else float(global_max_similarity)
    if reference < observed_max - 1e-12:
        raise DataValidationError(
            "Provided global_max_similarity is below this eligible set's maximum"
        )
    np.fill_diagonal(similarities, reference)
    distances = similarities_to_distances(
        similarities, global_max_similarity=reference
    )
    return LocalAlignmentDistances(ids, similarities, distances, reference)


@dataclass(frozen=True)
class SequenceRepresentations:
    repertoire_ids: Tuple[str, ...]
    alphabet: Tuple[Token, ...]
    transition_probability: np.ndarray
    bigram: np.ndarray
    phee_repeat: np.ndarray


def build_sequence_representations(
    repertoires: Mapping[str, Repertoire],
    *,
    alphabet: Sequence[Token] | None = None,
    phee_token: Token = "A",
    repeat_minimum: int = 2,
    repeat_maximum: int = 5,
) -> SequenceRepresentations:
    """Build the three inexpensive repertoire-level feature matrices."""

    ids = tuple(str(identifier) for identifier in repertoires)
    if not ids or len(set(ids)) != len(ids):
        raise DataValidationError("Repertoire IDs must be nonempty and unique")
    normalized = [_repertoire(repertoire) for repertoire in repertoires.values()]
    ordered_alphabet = tuple(alphabet) if alphabet is not None else discover_alphabet(normalized)
    transition = np.vstack(
        [transition_probability_vector(repertoire, ordered_alphabet) for repertoire in normalized]
    )
    bigram = np.vstack(
        [bigram_frequency_vector(repertoire, ordered_alphabet) for repertoire in normalized]
    )
    repeat = np.vstack(
        [
            phee_repeat_vector(
                repertoire,
                target=phee_token,
                minimum_length=repeat_minimum,
                maximum_length=repeat_maximum,
            )
            for repertoire in normalized
        ]
    )
    return SequenceRepresentations(ids, ordered_alphabet, transition, bigram, repeat)
