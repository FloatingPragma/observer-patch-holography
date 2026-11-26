"""Regression checks on the actual charged source and neutral release algebra."""

from itertools import combinations

import numpy as np
import pytest
from scipy.linalg import eigh

from model.candidate_laws.interval_echo.openwave.row_10_oscillations.audit import (
    environment,
    release_sectors,
    sector_from_leaf_read,
    source_reflections,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release import (
    flight,
    geometry,
)


def test_all_sharp_tip_words_have_reflections_on_the_native_source():
    source = source_reflections()
    assert source["charged_chart_dimension"] == 3584
    assert source["exact_commuting_arm_permutations"] == 6
    assert set(source["sharp_tip_word_stabilizers"]) == set(map(str, range(8)))


def test_remote_ambiguity_uses_both_leaf_loads_and_the_retained_directions():
    leaves = (165, 158)
    launch = (160, 167, 1)
    moved = (165, 167, 0)
    assert environment(launch, leaves)[1] != environment(moved, leaves)[1]
    assert (
        sector_from_leaf_read(launch, leaves)
        == sector_from_leaf_read(moved, leaves)
        == 1
    )
    opposite_launch = (160, 167, 0)
    assert environment(launch, leaves)[0] == environment(opposite_launch, leaves)[0]
    assert sector_from_leaf_read(opposite_launch, leaves) != sector_from_leaf_read(
        launch, leaves
    )


def test_invariant_on_every_native_word_and_every_reachable_word_for_both_signs():
    release = release_sectors()
    assert (release["full_release_words"], release["directed_nonzero_exchanges"]) == (
        960,
        3472,
    )
    assert release["off_diagonal_contact_entries"] == [6, 4]
    for classes in release["launch_sectors"].values():
        assert classes == dict(
            reachable_words_by_mask=[210, 211, 211, 210], environment_keys=14
        )


@pytest.mark.parametrize("launch", [(160, 167), (167, 160)])
def test_partial_trace_erases_every_launch_sector_coherence(launch):
    leaves = geometry.layout()[2]
    waits = (0, 256, 2048)
    residuals, guide_words = {}, {}
    columns = []
    for sector in range(4):
        start = (*launch, sector)
        states, hamiltonian = flight.chart(start)
        energies, vectors = eigh(hamiltonian)
        amplitudes = (
            np.exp(-1j * np.outer(waits, energies)) * vectors[states.index(start)].conj()
        ) @ vectors.T
        residual_index, guide_index = [], []
        for vacancy, excess, head in states:
            residual = environment((vacancy, excess, head), leaves)
            # A defect on a traced leaf is absent from the accessible guide.
            guide = tuple(
                port if port not in leaves else -1 for port in (vacancy, excess)
            )
            residual_index.append(residuals.setdefault(residual, len(residuals)))
            guide_index.append(guide_words.setdefault(guide, len(guide_words)))
        assert len(set(zip(residual_index, guide_index, strict=True))) == len(states)
        columns.append((amplitudes, residual_index, guide_index))

    assert (len(residuals), len(guide_words)) == (14, 211)
    for time_index in range(len(waits)):
        field = np.zeros((4, len(residuals), len(guide_words)), dtype=complex)
        for sector, (amplitudes, residual_index, guide_index) in enumerate(columns):
            field[sector, residual_index, guide_index] = amplitudes[time_index]
        np.testing.assert_allclose(
            np.sum(abs(field) ** 2, axis=(1, 2)), 1, atol=1e-12, rtol=0
        )
        for first, second in combinations(field, 2):
            cross_term = first.T @ second.conj()
            np.testing.assert_array_equal(cross_term, np.zeros_like(cross_term))
            # A coherent projection at the source onto its uniform record state
            # removes the tag and gives a positive control for the cross term.
            projected_first = first.sum(axis=0) / np.sqrt(len(residuals))
            projected_second = second.sum(axis=0) / np.sqrt(len(residuals))
            projected_cross = np.outer(projected_first, projected_second.conj())
            assert np.linalg.norm(projected_cross) > 1e-4
