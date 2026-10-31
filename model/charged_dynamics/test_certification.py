"""Exact independent inertia oracles and hostile interval-certificate controls."""
from copy import deepcopy
from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import certification as c
from . import chain


@pytest.mark.parametrize('matrix,expected', (
    (((0,),), (0, 1, 0)), (((-2,),), (1, 0, 0)),
    (((0, 1), (1, 0)), (1, 0, 1)),
    (((0, 1, 1), (1, 0, 1), (1, 1, 0)), (2, 0, 1)),
    (((0, 0, 0), (0, -1, 0), (0, 0, 2)), (1, 1, 1)),
))
def test_zero_diagonal_and_singular_exact_inertia(matrix, expected):
    assert c.inertia(matrix) == expected


def test_inertia_against_known_integer_congruence_including_large_entries():
    random = np.random.default_rng(1157)
    for size in range(1, 13):
        for trial in range(5):
            transform = np.triu(random.integers(-3, 4, size=(size, size))).astype(object)
            for i in range(size): transform[i, i] = 1
            diagonal = [(-1, 0, 1)[(i+trial) % 3]*(i+1) for i in range(size)]
            matrix = transform.T@np.diag(np.array(diagonal, dtype=object))@transform
            if trial == 4: matrix *= 1 << 200
            literal = [[int(value) for value in row] for row in matrix]
            expected = (sum(x < 0 for x in diagonal), sum(x == 0 for x in diagonal), sum(x > 0 for x in diagonal))
            assert c.inertia(literal) == expected


@pytest.fixture(scope='module')
def grid():
    return c.build_grid()


def test_exact_grid_proves_the_whole_path_and_complete_reference(grid):
    result = c.replay_grid(grid)
    assert result['charged_complement_gap_uniform_lower'] == '1/256'
    assert result['half_filled_excitation_gap_uniform_lower'] == '1/128'
    assert result['other_number_sector_separation_lower'] == '1/8'
    assert result['charged_rank_in_each_of_N3_N5'] == 4
    assert c.structural_certificate()['full_Fock_columns'] == 256


@pytest.mark.parametrize('damage', ('missing', 'duplicate', 'bool', 'float', 'false_interval', 'same_endpoint'))
def test_grid_replay_refuses_incomplete_or_forged_evidence(grid, damage):
    changed = deepcopy(grid)
    if damage == 'missing': changed.pop()
    elif damage == 'duplicate': changed[2][0] = changed[1][0]
    elif damage == 'bool': changed[0][0] = False
    elif damage == 'float': changed[3][1] = float(changed[3][1])
    elif damage == 'false_interval': changed[32][2] += 1000
    else:
        for row in changed: row[1:] = changed[0][1:]
    with pytest.raises((ValueError, ArithmeticError)): c.replay_grid(changed)


def test_floating_proposals_cannot_change_the_exact_grid(grid, monkeypatch):
    original = np.linalg.eigvalsh
    def displaced(matrix):
        return original(matrix)+2/c.DENOMINATOR
    monkeypatch.setattr(np.linalg, 'eigvalsh', displaced)
    assert c.build_grid() == grid


def test_exact_spin_reduction_rejects_an_altered_generator(monkeypatch):
    original = chain.hubbard_integer_parts
    def changed(up, down):
        addresses, base, link = original(up, down)
        if (up, down) == (2, 1):
            base = [list(row) for row in base]; base[0][0] += 1
            base = tuple(map(tuple, base))
        return addresses, base, link
    monkeypatch.setattr(chain, 'hubbard_integer_parts', changed)
    with pytest.raises(ArithmeticError, match='spin|particle-hole|gauge'):
        c.structural_certificate()


def test_short_time_positive_transfer_has_exact_scalar_controls():
    result = c.short_time_transfer()
    remainder = F(result['odd_Dyson_remainder_upper'])
    # A much longer independently summed rational series lies below the
    # geometric upper enclosure, and above its first omitted-order prefix.
    long_sum = sum((F(1, 2)**order/math.factorial(order) for order in range(3, 101, 2)), F())
    assert remainder >= long_sum
    assert remainder-long_sum < F(1, 10**27)
    amplitude = F(result['opposite_bare_code_amplitude_lower'])
    probability = F(result['opposite_bare_code_probability_lower'])
    assert amplitude > F(1, 20) and probability == amplitude*amplitude
    ideal_source_error = F(841, 1572864)
    assert amplitude-ideal_source_error > F(1, 20)


@pytest.mark.parametrize('matrix', ([[True]], [[1.]], [[1, 2], [3, 4]], []))
def test_inertia_strict_integer_domain(matrix):
    with pytest.raises(ValueError): c.inertia(matrix)
