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
