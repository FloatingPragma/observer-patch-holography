"""Complete total-two detector and independent block checks for contrast repair."""
from fractions import Fraction

import numpy as np
import pytest
from scipy.linalg import expm

from . import contrast_replay
from .contrast_replay import cell, family_pair_bound, primary_join, run, seam, seam_orientation_chart
from .scalar_channel.carrier import load_source


def test_total_two_sector_retains_double_occupancy_complement():
    h, states, adjacent = cell(8)
    kinetic, other_states, _ = cell(0)
    assert states == other_states
    assert h.shape == (78, 78) and len(set(states)) == 78
    doubles = np.array([a == b for a, b in states])
    assert doubles.sum() == 12 and adjacent.sum() == 30
    np.testing.assert_array_equal(h[np.ix_(doubles, ~doubles)], 0)
    np.testing.assert_array_equal(h[np.ix_(~doubles, doubles)], 0)
    expected_cost = np.where(doubles, 24, np.where(adjacent, 10, 12))
    np.testing.assert_array_equal(h - kinetic, np.diag(8 * expected_cost))

    # On double occupancy a site SWAP transports the entire load of two.
    # Its independent twelve-state Hamiltonian is L/2 + 24*kappa*I.
    adjacency = np.zeros((12, 12))
    for a, b in load_source('cell')['edges']:
        adjacency[a, b] = adjacency[b, a] = 1
    laplacian = np.diag(adjacency.sum(axis=1)) - adjacency
    np.testing.assert_array_equal(h[np.ix_(doubles, doubles)],
                                   laplacian / 2 + 192 * np.eye(12))


def test_full_cell_blocks_obey_the_ordered_theorem_hypotheses():
    h, _, adjacent = cell(8)
    A = h[np.ix_(adjacent, adjacent)]
    B = h[np.ix_(~adjacent, adjacent)]
    C = h[np.ix_(~adjacent, ~adjacent)]
    assert np.linalg.eigvalsh(A).max() <= 86 + 1e-12
    assert np.linalg.eigvalsh(C).min() >= 96 - 1e-12
    assert abs(B).sum(axis=0).max() <= 5
    assert abs(B).sum(axis=1).max() <= 2
    assert np.linalg.norm(B, 2)**2 <= 10
    assert family_pair_bound(8) == Fraction(2, 5)
    assert family_pair_bound(Fraction(7, 2)) == 1
    with pytest.raises(ValueError, match='kappa > 3'):
        family_pair_bound(3)


def test_complete_law_retains_zero_probability_double_outcomes():
    result = run()
    for experiment in result['cell']:
        assert len(experiment['outcomes']) == len(experiment['eigenvalues']) == 78
        doubles = np.array([a == b for a, b in experiment['outcomes']])
        for record in experiment['records']:
            law = np.array(record['probabilities'])
            weights = np.array(record['unnormalized_weights'])
            assert law.shape == (78,)
            assert weights.shape == (78,)
            assert np.isfinite(weights).all() and (weights >= 0).all()
            assert record['numerical_normalization'] == float(weights.sum())
            assert record['state_norm_error'] == abs(float(weights.sum()) - 1)
            assert record['state_norm_error'] < 1e-12
            np.testing.assert_array_equal(law, weights / record['numerical_normalization'])
            assert (law >= 0).all() and (law <= 1).all()
            assert law.sum() == pytest.approx(1, abs=1e-12)
            np.testing.assert_allclose(law[doubles], 0, atol=1e-25)
    assert result['all_time_bound']['numerator'] == 2
    assert result['all_time_bound']['denominator'] == 5
    assert max(r['no_adjacent_pair_probability']
               for r in result['cell'][1]['records']) > .4


def test_normalization_rejects_a_nonunitary_evolution(monkeypatch):
    original = contrast_replay.eigh
    def scaled_frame(*args, **kwargs):
        result = original(*args, **kwargs)
        if kwargs.get('eigvals_only', False):
            return result
        values, vectors = result
        return values, 1.001 * vectors
    monkeypatch.setattr(contrast_replay, 'eigh', scaled_frame)
    with pytest.raises(ValueError, match='Incomplete two-load instrument'):
        run()


def test_record_size_is_conserved_by_echo_but_changed_by_write():
    h, _, write = seam(8)
    # Independent memory labels, with each record carrying all nine load pairs.
    sizes = np.array([0] * 9 + [d * d for _coin in (0, 1)
                                for d in range(-2, 3) for _load in range(9)])
    size_operator = np.diag(sizes)
    np.testing.assert_allclose(size_operator @ h, h @ size_operator,
                               atol=1e-12, rtol=0)
    # The blank |2,0> input writes a record of contrast magnitude two.
    blank = np.eye(99)[:, 6]
    written = write @ blank
    np.testing.assert_array_equal(size_operator @ blank, np.zeros(99))
    np.testing.assert_allclose(size_operator @ written, 4 * written,
                               atol=1e-15, rtol=0)
    assert np.linalg.norm(size_operator @ write - write @ size_operator, 2) == pytest.approx(4)
    # Blank and recorded zero are separate invariant sectors of the echo.
    labels = np.array([-1] * 9 + [abs(d) for _coin in (0, 1)
                                  for d in range(-2, 3) for _load in range(9)])
    for label in (-1, 0, 1, 2):
        projection = np.diag((labels == label).astype(float))
        np.testing.assert_allclose(projection @ h, h @ projection,
                                   atol=1e-12, rtol=0)


def test_reversed_seam_chart_transports_every_memory_state_and_reference():
    chart = seam_orientation_chart(True)
    np.testing.assert_array_equal(chart[chart], np.arange(99))
    np.testing.assert_array_equal(seam_orientation_chart(False), np.arange(99))
    # Independently form the physical endpoint and memory label permutations.
    loads = [(a, b) for a in range(3) for b in range(3)]
    labels = [(c, d) for c in (0, 1) for d in range(-2, 3)]
    load_swap = np.eye(9)[:, [loads.index((b, a)) for a, b in loads]]
    record_swap = np.eye(10)[:, [labels.index((1-c, -d)) for c, d in labels]]
    expected = np.zeros((99, 99))
    expected[:9, :9] = load_swap
    expected[9:, 9:] = np.kron(record_swap, load_swap)
    actual = np.eye(99)[:, chart]
    np.testing.assert_array_equal(actual, expected)

    h, _, _ = seam(8)
    propagation = expm(-.13j * h)
    rng = np.random.default_rng(731)
    joint = rng.normal(size=(99, 2)) + 1j*rng.normal(size=(99, 2))
    joint /= np.linalg.norm(joint)
    np.testing.assert_allclose(propagation @ actual @ joint,
                               actual @ propagation @ joint, atol=3e-14, rtol=0)
    incomplete = np.kron(np.eye(11), load_swap)
    assert np.linalg.norm(propagation @ incomplete @ joint
                          - incomplete @ propagation @ joint) > .05

    result = primary_join()
    assert len(result['seam_charts']) == 30
    edge = next(item for item in result['seam_charts'] if item['source_edge'] == [4, 6])
    assert edge == dict(source_edge=[4, 6], primary_edge=[10, 11], orientation='reversed')
    for check in result['checks']:
        assert check['full_99_reversed_seam_error'] < 1e-12
        assert check['load_only_orientation_control_error'] > .5


def test_the_port_map_addresses_every_cell_of_every_federation():
    """The labels of the construction manual reach the local ports of every cell of the three federations."""
    from .contrast_memory_spherical_v1 import carrier
    from .scalar_channel.carrier import EDGES
    port_map = (0, 1, 5, 7, 11, 9, 10, 8, 4, 6, 2, 3)
    cells = 0
    for level in (0, 1, 2):
        fed = carrier.federation(level)
        seams = {tuple(sorted(map(int, e))) for e in fed.seams}
        for k in range(fed.n // 12):
            for a, b in EDGES:
                assert tuple(sorted((12 * k + port_map[a], 12 * k + port_map[b]))) in seams
            cells += 1
    assert cells == 420
    # the port opposite to manual port p is 11 - p; in a federation cell local port 0 is opposite local port 3
    assert port_map[11 - 0] == 3
