from fractions import Fraction
from itertools import permutations

import numpy as np
import pytest

from . import optimality


@pytest.mark.parametrize('bad', [None, [], (True,), (1.0,), (2, 1), (1, 1), (-1,), (12,)])
def test_invalid_locality_domain_refused(bad):
    with pytest.raises(ValueError):
        optimality.bounds(bad)


def _partial_trace(matrix, core):
    # Independent sparse embedding into twelve tensor factors; basis here is
    # lexicographic ordered positions, deliberately different from core.py.
    states = tuple(permutations(range(12), 2))
    configs = [tuple(1 if x == a else 2 if x == b else 0 for x in range(12))
               for a, b in states]
    kept = sorted({tuple(c[x] for x in core) for c in configs})
    index = {value: i for i, value in enumerate(kept)}
    outside = tuple(x for x in range(12) if x not in core)
    result = np.zeros((len(kept), len(kept)), complex)
    for i, ci in enumerate(configs):
        for j, cj in enumerate(configs):
            if all(ci[x] == cj[x] for x in outside):
                result[index[tuple(ci[x] for x in core)], index[tuple(cj[x] for x in core)]] += matrix[i, j]
    return result


@pytest.mark.parametrize('core', [(0, 1), (0, 5, 11), tuple(range(11)), tuple(range(12))])
def test_sharp_state_attains_universal_exterior_bound(core):
    states = tuple(permutations(range(12), 2))
    index = {pair: i for i, pair in enumerate(states)}
    identity = np.eye(132)
    swap = np.zeros_like(identity)
    for i, (a, b) in enumerate(states):
        swap[index[b, a], i] = 1
    sigma = identity/132
    sharp = (identity+swap)/132
    assert np.allclose(swap@sharp, sharp)
    exterior = tuple(x for x in range(12) if x not in core)
    difference = _partial_trace(sharp-sigma, exterior)
    measured = np.linalg.svd(difference, compute_uv=False).sum()/2
    bound = Fraction(optimality.bounds(core)['sharp_sector_min_exterior_trace_distance'])
    assert measured == pytest.approx(float(bound), abs=1e-14)
    if len(exterior) >= 2:
        assert measured > 0


def test_local_preparation_attains_ceiling_and_has_unchanged_full_exterior():
    states = tuple(permutations(range(12), 2))
    index = {pair: i for i, pair in enumerate(states)}
    rho = np.eye(132)/132
    core = (0, 5, 11)
    for i, (a, b) in enumerate(states):
        if a in core and b in core:
            rho[i, index[b, a]] = 1/132
    charge = sum(rho[index[b, a], i] for i, (a, b) in enumerate(states))
    assert charge == pytest.approx(float(Fraction(optimality.bounds(core)['local_preparation_max_abs_charge'])))
    assert np.linalg.eigvalsh(rho).min() >= -1e-15
    exterior = tuple(x for x in range(12) if x not in core)
    assert np.max(np.abs(_partial_trace(rho-np.eye(132)/132, exterior))) == 0
    # This saturating local preparation is still far from a sharp sector.
    assert 0 < charge < 1


def test_sharp_exterior_obstruction_requires_the_one_A_one_B_carrier():
    # Split the literal tensor space into port 0 and its eleven-port exterior.
    # Its reference marginal has 132 distinct basis states, each with weight
    # 1/132: 110 both-outside configurations and 22 single-label ones.
    exterior = tuple(range(1, 12))
    states = tuple(permutations(range(12), 2))
    configs = [tuple(1 if x == a else 2 if x == b else 0 for x in exterior)
               for a, b in states]
    basis = sorted(set(configs))
    assert len(basis) == 132
    index = {value: i for i, value in enumerate(basis)}
    flip = [index[tuple(3 - x if x else 0 for x in value)] for value in basis]
    assert [flip[flip[i]] for i in range(132)] == list(range(132))
    swap = np.eye(132)[flip]
    plus, minus = (np.eye(132) + swap) / 2, (np.eye(132) - swap) / 2
    vacuum = np.diag([1., 0., 0.])
    core_minus = np.array([[0., 0., 0.], [0., .5, -.5], [0., -.5, .5]])
    # Matching core/exterior character projectors are paired to make Q=+1.
    # Unlike H132, this ambient state permits extra and missing labels.
    rho = (np.kron(vacuum, plus) + np.kron(core_minus, minus)) / 132
    global_flip = [132 * (0, 2, 1)[a] + flip[b]
                   for a in range(3) for b in range(132)]
    np.testing.assert_array_equal(rho[global_flip], rho)
    assert np.trace(rho) == pytest.approx(1)
    reduced = np.trace(rho.reshape(3, 132, 3, 132), axis1=0, axis2=2)
    sigma_exterior = _partial_trace(np.eye(132) / 132, exterior)
    np.testing.assert_array_equal(reduced, sigma_exterior)
    # Thus its exterior distance is zero although the H132 lower bound is
    # 55/132. Check explicitly that the state violates that bound's carrier.
    admitted = [132 * a + i for a in range(3) for i, value in enumerate(basis)
                if int(a == 1) + value.count(1) == 1
                and int(a == 2) + value.count(2) == 1]
    sector_probability = np.diag(rho)[admitted].sum()
    assert sector_probability == pytest.approx(121 / 264)
    assert Fraction(optimality.bounds((0,))['sharp_sector_min_exterior_trace_distance']) == Fraction(55, 132)


def test_all_subsets_and_extreme_cores():
    value = optimality.certificate()
    assert value['all_port_subsets_checked'] == 4096
    rows = value['bounds_by_core_size']
    assert rows[0]['local_preparation_max_abs_charge'] == '0'
    assert rows[0]['sharp_sector_min_exterior_trace_distance'] == '1/2'
    assert rows[2]['local_preparation_max_abs_charge'] == '1/66'
    assert rows[2]['sharp_sector_min_exterior_trace_distance'] == '15/44'
    assert rows[12]['local_preparation_max_abs_charge'] == '1'
    assert rows[12]['sharp_sector_min_exterior_trace_distance'] == '0'
    assert 'one-A/one-B H132 carrier' in value['sharp_domain']
