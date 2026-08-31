"""Exact, constrained dynamic time warping for cached acoustic analysis.

The pairwise reference implementation is intentionally separate from the
full-matrix driver.  Building all call-to-call distances requires an explicit
``allow_expensive=True`` flag and is never triggered by importing the package.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Mapping, Optional, Tuple

import numpy as np

from .exceptions import DataValidationError, ExpensiveComputationDisabled


@dataclass(frozen=True)
class DTWConfig:
    normalized_time_constraint: float = 0.20

    def __post_init__(self) -> None:
        if not 0 <= self.normalized_time_constraint <= 1:
            raise DataValidationError("normalized_time_constraint must be in [0, 1]")


DEFAULT_DTW = DTWConfig()


@dataclass(frozen=True)
class DTWResult:
    cumulative_cost: float
    alignment_steps: int
    distance: float
    path: Optional[Tuple[Tuple[int, int], ...]] = None


def cosine_frame_costs(frames_a: np.ndarray, frames_b: np.ndarray) -> np.ndarray:
    """Return all cosine distances for arrays shaped ``frames × features``."""

    a = np.asarray(frames_a, dtype=np.float64)
    b = np.asarray(frames_b, dtype=np.float64)
    if a.ndim != 2 or b.ndim != 2:
        raise DataValidationError("DTW inputs must be 2-D frames-by-features arrays")
    if a.shape[0] == 0 or b.shape[0] == 0:
        raise DataValidationError("DTW inputs must each contain at least one frame")
    if a.shape[1] != b.shape[1]:
        raise DataValidationError(
            f"DTW feature dimensions differ: {a.shape[1]} versus {b.shape[1]}"
        )
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise DataValidationError("DTW input frames must be finite")

    dot = a @ b.T
    norm_a = np.linalg.norm(a, axis=1)
    norm_b = np.linalg.norm(b, axis=1)
    denominator = norm_a[:, None] * norm_b[None, :]
    similarity = np.divide(
        dot,
        denominator,
        out=np.zeros_like(dot),
        where=denominator != 0,
    )
    both_zero = (norm_a[:, None] == 0) & (norm_b[None, :] == 0)
    similarity[both_zero] = 1.0
    distance = np.clip(1.0 - similarity, 0.0, 2.0)
    distance[np.abs(distance) < 1e-12] = 0.0
    return distance


def _relative_positions(length: int) -> np.ndarray:
    return np.array([0.0]) if length == 1 else np.linspace(0.0, 1.0, length)


def exact_dtw(
    frames_a: np.ndarray,
    frames_b: np.ndarray,
    *,
    config: DTWConfig = DEFAULT_DTW,
    return_path: bool = False,
) -> DTWResult:
    """Find the minimum cumulative-cost path inside the normalized-time band.

    Local frame cost is cosine distance.  The reported distance is the optimal
    path's cumulative cost divided by its number of alignment cells.  Ties are
    resolved deterministically in diagonal, vertical, horizontal order.
    """

    local = cosine_frame_costs(frames_a, frames_b)
    n, m = local.shape
    positions_a = _relative_positions(n)
    positions_b = _relative_positions(m)
    allowed = (
        np.abs(positions_a[:, None] - positions_b[None, :])
        <= config.normalized_time_constraint + 1e-12
    )

    cumulative = np.full((n, m), np.inf, dtype=np.float64)
    steps = np.zeros((n, m), dtype=np.int64)
    predecessor = np.full((n, m), -1, dtype=np.int8) if return_path else None
    if allowed[0, 0]:
        cumulative[0, 0] = local[0, 0]
        steps[0, 0] = 1

    # Codes and ordering: 0 diagonal, 1 vertical, 2 horizontal.
    moves = ((-1, -1, 0), (-1, 0, 1), (0, -1, 2))
    for i in range(n):
        for j in range(m):
            if (i == 0 and j == 0) or not allowed[i, j]:
                continue
            candidates = []
            for di, dj, code in moves:
                previous_i, previous_j = i + di, j + dj
                if previous_i < 0 or previous_j < 0:
                    continue
                previous_cost = cumulative[previous_i, previous_j]
                if np.isfinite(previous_cost):
                    candidates.append(
                        (previous_cost, code, previous_i, previous_j)
                    )
            if not candidates:
                continue
            previous_cost, code, previous_i, previous_j = min(
                candidates, key=lambda item: (item[0], item[1])
            )
            cumulative[i, j] = previous_cost + local[i, j]
            steps[i, j] = steps[previous_i, previous_j] + 1
            if predecessor is not None:
                predecessor[i, j] = code

    final_cost = float(cumulative[-1, -1])
    final_steps = int(steps[-1, -1])
    if not np.isfinite(final_cost) or final_steps == 0:
        raise DataValidationError(
            "No DTW path satisfies the normalized-time constraint; verify frame "
            "counts or explicitly choose a wider prespecified band"
        )

    path = None
    if predecessor is not None:
        reversed_path = [(n - 1, m - 1)]
        i, j = n - 1, m - 1
        while i or j:
            code = int(predecessor[i, j])
            if code == 0:
                i, j = i - 1, j - 1
            elif code == 1:
                i -= 1
            elif code == 2:
                j -= 1
            else:  # Defensive: an unreachable cell would have failed above.
                raise DataValidationError("Invalid DTW traceback state")
            reversed_path.append((i, j))
        path = tuple(reversed(reversed_path))
    return DTWResult(final_cost, final_steps, final_cost / final_steps, path)


def exact_dtw_distance(
    frames_a: np.ndarray,
    frames_b: np.ndarray,
    *,
    config: DTWConfig = DEFAULT_DTW,
) -> float:
    """Convenience wrapper returning path-length-normalized distance only."""

    return exact_dtw(frames_a, frames_b, config=config).distance


def pairwise_exact_dtw_condensed(
    frames_by_call: Mapping[str, np.ndarray],
    *,
    config: DTWConfig = DEFAULT_DTW,
    allow_expensive: bool = False,
) -> tuple[Tuple[str, ...], np.ndarray]:
    """Build an all-call condensed matrix only after explicit full-mode opt-in."""

    if not allow_expensive:
        raise ExpensiveComputationDisabled(
            "All-call DTW is disabled by default; load "
            "data/cache/dtw_distances_condensed.npz or pass allow_expensive=True"
        )
    call_ids = tuple(str(call_id) for call_id in frames_by_call)
    if len(call_ids) < 2 or len(set(call_ids)) != len(call_ids):
        raise DataValidationError("At least two unique call IDs are required")
    values = np.empty(len(call_ids) * (len(call_ids) - 1) // 2, dtype=np.float64)
    cursor = 0
    frame_arrays = tuple(frames_by_call.values())
    for i, j in combinations(range(len(call_ids)), 2):
        values[cursor] = exact_dtw_distance(
            frame_arrays[i], frame_arrays[j], config=config
        )
        cursor += 1
    return call_ids, values
