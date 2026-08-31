import numpy as np
import pytest

from marmoset_convergence import (
    DTWConfig,
    DataValidationError,
    ExpensiveComputationDisabled,
    exact_dtw,
    exact_dtw_distance,
    pairwise_exact_dtw_condensed,
)


def test_exact_dtw_identical_trajectories_have_zero_distance_and_diagonal_path():
    frames = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])

    result = exact_dtw(frames, frames, return_path=True)

    assert result.cumulative_cost == pytest.approx(0)
    assert result.alignment_steps == 3
    assert result.distance == pytest.approx(0)
    assert result.path == ((0, 0), (1, 1), (2, 2))


def test_exact_dtw_uses_cosine_frame_cost():
    left = np.array([[1.0, 0.0], [1.0, 0.0]])
    right = np.array([[0.0, 1.0], [0.0, 1.0]])

    assert exact_dtw_distance(left, right) == pytest.approx(1.0)
    assert exact_dtw_distance(right, left) == pytest.approx(1.0)


def test_dtw_normalizes_cumulative_cost_by_alignment_steps():
    left = np.array([[1.0, 0.0], [1.0, 0.0]])
    right = np.array([[0.0, 1.0], [0.0, 1.0]])
    result = exact_dtw(left, right)

    assert result.cumulative_cost == pytest.approx(2.0)
    assert result.alignment_steps == 2
    assert result.distance == pytest.approx(1.0)


def test_normalized_time_band_can_reject_implausible_length_warp():
    short = np.ones((2, 3))
    long = np.ones((100, 3))

    with pytest.raises(DataValidationError, match="No DTW path"):
        exact_dtw(short, long, config=DTWConfig(0.20))
    assert exact_dtw_distance(short, long, config=DTWConfig(1.0)) == pytest.approx(0)


def test_zero_frames_have_defined_cosine_behavior():
    zeros = np.zeros((2, 2))
    nonzero = np.ones((2, 2))

    assert exact_dtw_distance(zeros, zeros) == pytest.approx(0)
    assert exact_dtw_distance(zeros, nonzero) == pytest.approx(1)


def test_all_call_dtw_is_disabled_without_explicit_opt_in():
    frames = {"a": np.ones((2, 2)), "b": np.ones((2, 2))}
    with pytest.raises(ExpensiveComputationDisabled):
        pairwise_exact_dtw_condensed(frames)

    ids, distances = pairwise_exact_dtw_condensed(frames, allow_expensive=True)
    assert ids == ("a", "b")
    np.testing.assert_allclose(distances, [0])
