from math import log, sqrt

import numpy as np
import pytest

from marmoset_convergence import (
    DataValidationError,
    ExpensiveComputationDisabled,
    bigram_frequency_vector,
    local_alignment_distance_matrix,
    normalized_local_alignment_similarity,
    phee_repeat_vector,
    repertoire_local_alignment_similarity,
    similarities_to_distances,
    smith_waterman_score,
    transition_probability_vector,
)


def test_transition_probabilities_are_row_conditional():
    alphabet = ("A", "B")
    vector = transition_probability_vector(["AAB", "ABB"], alphabet)
    matrix = vector.reshape(2, 2)

    np.testing.assert_allclose(matrix, [[1 / 3, 2 / 3], [0, 1]])


def test_bigrams_are_normalized_by_global_total():
    alphabet = ("A", "B")
    vector = bigram_frequency_vector(["AAB", "ABB"], alphabet)

    np.testing.assert_allclose(vector.reshape(2, 2), [[0.25, 0.5], [0, 0.25]])


def test_unknown_token_is_not_silently_dropped():
    with pytest.raises(DataValidationError, match="fixed alphabet"):
        transition_probability_vector(["AC"], ("A", "B"))


def test_repeat_profile_counts_exact_lengths_and_excludes_long_runs():
    vector = phee_repeat_vector(["AA", "AAA", "AAAAAA", "BAAAB"], target="A")

    # One length-2 run and two length-3 runs; the length-6 run is excluded.
    np.testing.assert_allclose(vector, [1 / 3, 2 / 3, 0, 0])


def test_repeat_profile_without_eligible_runs_is_zero_not_nan():
    np.testing.assert_array_equal(phee_repeat_vector(["AB", "BA"]), np.zeros(4))


def test_smith_waterman_uses_manuscript_scores():
    assert smith_waterman_score("AAB", "ACB") == 3
    assert smith_waterman_score("AAA", "AAA") == 6
    assert smith_waterman_score("A", "B") == 0


def test_local_alignment_normalizes_by_log_combined_lengths():
    observed = normalized_local_alignment_similarity("AA", "AA")
    assert observed == pytest.approx(4 / log(4))


def test_repertoire_similarity_means_all_pairs_instead_of_best_matches():
    observed = repertoire_local_alignment_similarity(["AA", "BB"], ["AA"])
    expected = ((4 / log(4)) + 0) / 2
    assert observed == pytest.approx(expected)


def test_similarity_to_distance_uses_supplied_fixed_global_reference():
    similarities = np.array([[9.0, 2.0], [2.0, 9.0]])
    distances = similarities_to_distances(
        similarities, global_max_similarity=5.0
    )
    np.testing.assert_allclose(distances, [[0, 3], [3, 0]])


def test_similarity_to_distance_rejects_too_small_global_reference():
    similarities = np.array([[0.0, 2.0], [2.0, 0.0]])
    with pytest.raises(DataValidationError, match="below an observed"):
        similarities_to_distances(similarities, global_max_similarity=1.0)


def test_local_alignment_matrix_is_opt_in_and_records_global_reference():
    repertoires = {"r1": ["AA"], "r2": ["AA"], "r3": ["BB"]}
    with pytest.raises(ExpensiveComputationDisabled):
        local_alignment_distance_matrix(repertoires)

    result = local_alignment_distance_matrix(repertoires, allow_expensive=True)

    assert result.global_max_similarity == pytest.approx(4 / log(4))
    np.testing.assert_allclose(np.diag(result.distances), 0)
    assert result.distances[0, 1] == 0
    assert result.distances[0, 2] == pytest.approx(result.global_max_similarity)
