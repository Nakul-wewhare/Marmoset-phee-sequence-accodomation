"""Validated condensed symmetric matrices and stable NPZ cache contracts."""

from __future__ import annotations

from dataclasses import dataclass, field
from math import isqrt
from pathlib import Path
from typing import Dict, Iterable, Literal, Mapping, Sequence, Tuple, Union

import numpy as np

from .exceptions import CacheValidationError, DataValidationError

PathLike = Union[str, Path]


def condensed_length(size: int) -> int:
    if not isinstance(size, int) or size < 0:
        raise DataValidationError("Matrix size must be a non-negative integer")
    return size * (size - 1) // 2


def infer_square_size(length: int) -> int:
    """Infer ``n`` from a SciPy-order condensed length ``n * (n - 1) / 2``."""

    if not isinstance(length, int) or length < 0:
        raise DataValidationError("Condensed length must be a non-negative integer")
    discriminant = 1 + 8 * length
    root = isqrt(discriminant)
    if root * root != discriminant or (1 + root) % 2:
        raise DataValidationError(f"{length} is not a valid condensed-matrix length")
    return (1 + root) // 2


def condensed_index(size: int, row: int, column: int) -> int:
    """Return the SciPy ``squareform`` index for one off-diagonal pair."""

    condensed_length(size)
    if not (0 <= row < size and 0 <= column < size):
        raise IndexError(f"Matrix indices {(row, column)} are outside size {size}")
    if row == column:
        raise IndexError("The diagonal is not stored in condensed form")
    if row > column:
        row, column = column, row
    return size * row - row * (row + 1) // 2 + column - row - 1


def square_to_condensed(matrix: np.ndarray) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64)
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
        raise DataValidationError("Distance matrix must be square")
    if not np.isfinite(values).all():
        raise DataValidationError("Distance matrix must be finite")
    if not np.allclose(values, values.T, rtol=1e-10, atol=1e-12):
        raise DataValidationError("Distance matrix must be symmetric")
    if not np.allclose(np.diag(values), 0.0, rtol=0.0, atol=1e-12):
        raise DataValidationError("Distance matrix diagonal must be zero")
    if (values < -1e-12).any():
        raise DataValidationError("Distances cannot be negative")
    upper = np.triu_indices(values.shape[0], k=1)
    result = values[upper].copy()
    result[np.abs(result) < 1e-12] = 0.0
    return result


def condensed_to_square(values: np.ndarray, *, size: int | None = None) -> np.ndarray:
    vector = np.asarray(values, dtype=np.float64)
    if vector.ndim != 1:
        raise DataValidationError("Condensed distances must be a 1-D array")
    inferred = infer_square_size(len(vector))
    if size is not None and size != inferred:
        raise DataValidationError(
            f"Condensed length implies size {inferred}, not requested size {size}"
        )
    if not np.isfinite(vector).all() or (vector < -1e-12).any():
        raise DataValidationError("Condensed distances must be finite and non-negative")
    matrix = np.zeros((inferred, inferred), dtype=np.float64)
    upper = np.triu_indices(inferred, k=1)
    matrix[upper] = vector
    matrix[(upper[1], upper[0])] = vector
    return matrix


@dataclass(frozen=True)
class CondensedDistanceMatrix:
    ids: Tuple[str, ...]
    distances: np.ndarray
    metric: str
    _index: Mapping[str, int] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        identifiers = tuple(str(identifier) for identifier in self.ids)
        if len(set(identifiers)) != len(identifiers):
            raise DataValidationError("Distance-matrix IDs must be unique")
        values = np.asarray(self.distances, dtype=np.float64)
        if len(values) != condensed_length(len(identifiers)):
            raise DataValidationError(
                f"Metric {self.metric!r} has {len(values)} distances for "
                f"{len(identifiers)} IDs"
            )
        if values.ndim != 1 or not np.isfinite(values).all() or (values < -1e-12).any():
            raise DataValidationError(
                f"Metric {self.metric!r} distances must be a finite non-negative vector"
            )
        object.__setattr__(self, "ids", identifiers)
        object.__setattr__(self, "distances", values)
        object.__setattr__(self, "metric", str(self.metric))
        object.__setattr__(
            self,
            "_index",
            {identifier: position for position, identifier in enumerate(identifiers)},
        )

    @property
    def size(self) -> int:
        return len(self.ids)

    def square(self) -> np.ndarray:
        return condensed_to_square(self.distances, size=self.size)

    def distance(self, left_id: str, right_id: str) -> float:
        try:
            left = self._index[str(left_id)]
            right = self._index[str(right_id)]
        except KeyError as exc:
            raise DataValidationError(
                f"ID {exc.args[0]!r} is absent from the {self.metric!r} matrix"
            ) from exc
        if left == right:
            return 0.0
        return float(self.distances[condensed_index(self.size, left, right)])

    def mean_between(
        self, left_ids: Iterable[str], right_ids: Iterable[str]
    ) -> float:
        """Mean of every cross-set item distance, including repeats if supplied."""

        left = tuple(str(identifier) for identifier in left_ids)
        right = tuple(str(identifier) for identifier in right_ids)
        if not left or not right:
            raise DataValidationError("Cannot average distances over an empty item set")
        try:
            left_positions = np.fromiter(
                (self._index[identifier] for identifier in left),
                dtype=np.int64,
                count=len(left),
            )
            right_positions = np.fromiter(
                (self._index[identifier] for identifier in right),
                dtype=np.int64,
                count=len(right),
            )
        except KeyError as exc:
            raise DataValidationError(
                f"ID {exc.args[0]!r} is absent from the {self.metric!r} matrix"
            ) from exc

        rows = left_positions[:, np.newaxis]
        columns = right_positions[np.newaxis, :]
        diagonal = rows == columns
        lower = np.minimum(rows, columns)
        upper = np.maximum(rows, columns)
        indexes = (
            self.size * lower
            - lower * (lower + 1) // 2
            + upper
            - lower
            - 1
        )
        values = np.zeros(indexes.shape, dtype=np.float64)
        values[~diagonal] = self.distances[indexes[~diagonal]]
        return float(values.mean())


@dataclass(frozen=True)
class FeatureDistanceMatrix:
    """Lazy Euclidean or cosine distances over an ordered feature table.

    Only requested cross-set blocks are evaluated.  This avoids constructing a
    full condensed matrix for vector representations when repertoire means are
    the only downstream quantity required.
    """

    ids: Tuple[str, ...]
    features: np.ndarray
    metric: str
    distance_kind: Literal["euclidean", "cosine"] = "euclidean"
    _index: Mapping[str, int] = field(init=False, repr=False, compare=False)
    _normalized: np.ndarray | None = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        identifiers = tuple(str(identifier) for identifier in self.ids)
        if not identifiers or len(set(identifiers)) != len(identifiers):
            raise DataValidationError("Feature-matrix IDs must be non-empty and unique")
        values = np.asarray(self.features, dtype=np.float64)
        if (
            values.ndim != 2
            or values.shape[0] != len(identifiers)
            or values.shape[1] < 1
        ):
            raise DataValidationError(
                "Feature matrix must have one non-empty feature row per ID"
            )
        if not np.isfinite(values).all():
            raise DataValidationError("Feature matrix contains non-finite values")
        kind = str(self.distance_kind).casefold()
        if kind not in {"euclidean", "cosine"}:
            raise DataValidationError(
                "Feature distance_kind must be 'euclidean' or 'cosine'"
            )
        normalized = None
        if kind == "cosine":
            norms = np.linalg.norm(values, axis=1)
            invalid_rows = np.flatnonzero(~np.isfinite(norms) | (norms == 0))
            if invalid_rows.size:
                preview = ", ".join(identifiers[index] for index in invalid_rows[:5])
                raise DataValidationError(
                    "Cosine-distance features contain invalid or zero-norm row(s): "
                    + preview
                )
            normalized = values / norms[:, np.newaxis]
        object.__setattr__(self, "ids", identifiers)
        object.__setattr__(self, "features", values)
        object.__setattr__(self, "metric", str(self.metric))
        object.__setattr__(self, "distance_kind", kind)
        object.__setattr__(
            self,
            "_index",
            {identifier: position for position, identifier in enumerate(identifiers)},
        )
        object.__setattr__(self, "_normalized", normalized)

    def _positions(self, ids: Iterable[str]) -> np.ndarray:
        identifiers = tuple(str(identifier) for identifier in ids)
        if not identifiers:
            raise DataValidationError("Cannot average distances over an empty item set")
        try:
            return np.fromiter(
                (self._index[identifier] for identifier in identifiers),
                dtype=np.int64,
                count=len(identifiers),
            )
        except KeyError as exc:
            raise DataValidationError(
                f"ID {exc.args[0]!r} is absent from the {self.metric!r} feature matrix"
            ) from exc

    def mean_between(
        self, left_ids: Iterable[str], right_ids: Iterable[str]
    ) -> float:
        """Mean distance across one Cartesian product of feature rows."""

        left = self._positions(left_ids)
        right = self._positions(right_ids)
        if self.distance_kind == "euclidean":
            differences = (
                self.features[left, np.newaxis, :]
                - self.features[np.newaxis, right, :]
            )
            return float(np.linalg.norm(differences, axis=2).mean())

        assert self._normalized is not None  # Established in __post_init__.
        similarities = self._normalized[left] @ self._normalized[right].T
        distances = 1.0 - np.clip(similarities, -1.0, 1.0)
        return float(distances.mean())


@dataclass(frozen=True)
class DistanceBundle:
    ids: Tuple[str, ...]
    matrices: Mapping[str, CondensedDistanceMatrix]

    def __post_init__(self) -> None:
        ids = tuple(str(identifier) for identifier in self.ids)
        if not self.matrices:
            raise DataValidationError("A distance bundle must contain at least one metric")
        for metric, matrix in self.matrices.items():
            if matrix.ids != ids:
                raise DataValidationError(
                    f"Metric {metric!r} does not use the bundle's ordered IDs"
                )
        object.__setattr__(self, "ids", ids)


@dataclass(frozen=True)
class SequenceRepresentationCache:
    repertoire_ids: Tuple[str, ...]
    transition_probability: np.ndarray
    bigram: np.ndarray
    phee_repeat: np.ndarray

    def __post_init__(self) -> None:
        ids = tuple(str(identifier) for identifier in self.repertoire_ids)
        if not ids or len(set(ids)) != len(ids):
            raise CacheValidationError("Representation repertoire_ids must be unique")
        for name in ("transition_probability", "bigram", "phee_repeat"):
            values = np.asarray(getattr(self, name), dtype=np.float64)
            if values.ndim != 2 or values.shape[0] != len(ids):
                raise CacheValidationError(
                    f"{name} must have one feature row per repertoire_id"
                )
            if not np.isfinite(values).all():
                raise CacheValidationError(f"{name} contains non-finite values")
            object.__setattr__(self, name, values)
        object.__setattr__(self, "repertoire_ids", ids)


def feature_distance_condensed(features: np.ndarray) -> np.ndarray:
    """Euclidean distances between all rows, in condensed order."""

    values = np.asarray(features, dtype=np.float64)
    if values.ndim != 2 or values.shape[0] < 2 or not np.isfinite(values).all():
        raise DataValidationError("Features must be a finite 2-D array with at least two rows")
    distances = np.empty(condensed_length(values.shape[0]), dtype=np.float64)
    cursor = 0
    for i in range(values.shape[0] - 1):
        block = np.linalg.norm(values[i + 1 :] - values[i], axis=1)
        distances[cursor : cursor + len(block)] = block
        cursor += len(block)
    return distances


def _read_npz(path: PathLike) -> Dict[str, np.ndarray]:
    cache_path = Path(path)
    if not cache_path.is_file():
        raise CacheValidationError(f"Required cache is missing: {cache_path}")
    try:
        with np.load(cache_path, allow_pickle=False) as archive:
            return {key: archive[key] for key in archive.files}
    except (OSError, ValueError) as exc:
        raise CacheValidationError(f"Could not read NPZ cache {cache_path}: {exc}") from exc


def _ids_from_array(values: np.ndarray, key: str) -> Tuple[str, ...]:
    if values.ndim != 1:
        raise CacheValidationError(f"NPZ key {key!r} must be a 1-D ID array")
    ids = tuple(str(value) for value in values.tolist())
    if len(set(ids)) != len(ids):
        raise CacheValidationError(f"NPZ key {key!r} contains duplicate IDs")
    return ids


def load_dtw_distance_cache(path: PathLike) -> CondensedDistanceMatrix:
    """Load the exact ``call_ids``/``distances`` DTW NPZ contract."""

    archive = _read_npz(path)
    expected = {"call_ids", "distances"}
    if set(archive) != expected:
        raise CacheValidationError(
            f"DTW cache keys must be {sorted(expected)}, found {sorted(archive)}"
        )
    ids = _ids_from_array(archive["call_ids"], "call_ids")
    try:
        return CondensedDistanceMatrix(ids, archive["distances"], "dtw")
    except DataValidationError as exc:
        raise CacheValidationError(str(exc)) from exc


def load_sequence_distance_cache(path: PathLike) -> CondensedDistanceMatrix:
    """Load the exact ``repertoire_ids``/``local_alignment`` NPZ contract."""

    archive = _read_npz(path)
    expected = {"repertoire_ids", "local_alignment"}
    if set(archive) != expected:
        raise CacheValidationError(
            f"Sequence-distance cache keys must be {sorted(expected)}, "
            f"found {sorted(archive)}"
        )
    ids = _ids_from_array(archive["repertoire_ids"], "repertoire_ids")
    try:
        return CondensedDistanceMatrix(ids, archive["local_alignment"], "local_alignment")
    except DataValidationError as exc:
        raise CacheValidationError(str(exc)) from exc


def load_sequence_representation_cache(path: PathLike) -> SequenceRepresentationCache:
    """Load the fixed four-key sequence-representation NPZ contract."""

    archive = _read_npz(path)
    expected = {"repertoire_ids", "transition_probability", "bigram", "phee_repeat"}
    if set(archive) != expected:
        raise CacheValidationError(
            f"Sequence-representation cache keys must be {sorted(expected)}, "
            f"found {sorted(archive)}"
        )
    return SequenceRepresentationCache(
        _ids_from_array(archive["repertoire_ids"], "repertoire_ids"),
        archive["transition_probability"],
        archive["bigram"],
        archive["phee_repeat"],
    )


def save_dtw_distance_cache(
    path: PathLike, call_ids: Sequence[str], distances: np.ndarray
) -> None:
    """Explicitly write the stable DTW NPZ contract."""

    matrix = CondensedDistanceMatrix(tuple(call_ids), distances, "dtw")
    np.savez_compressed(
        Path(path), call_ids=np.asarray(matrix.ids), distances=matrix.distances
    )


def save_sequence_distance_cache(
    path: PathLike, repertoire_ids: Sequence[str], local_alignment: np.ndarray
) -> None:
    matrix = CondensedDistanceMatrix(
        tuple(repertoire_ids), local_alignment, "local_alignment"
    )
    np.savez_compressed(
        Path(path),
        repertoire_ids=np.asarray(matrix.ids),
        local_alignment=matrix.distances,
    )


def save_sequence_representation_cache(
    path: PathLike,
    cache: SequenceRepresentationCache,
) -> None:
    np.savez_compressed(
        Path(path),
        repertoire_ids=np.asarray(cache.repertoire_ids),
        transition_probability=cache.transition_probability,
        bigram=cache.bigram,
        phee_repeat=cache.phee_repeat,
    )
