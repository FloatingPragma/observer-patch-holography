import math

import numpy as np
import pytest
from scipy.linalg import expm

from ..cube_eta_protection import cube, spectrum
from ..full_fock_field.field import FullFock, evolve
from . import factorization as f, certified_power as cp, budgets


def test_exact_original_source_grouping_eta_phase_and_singlet_chart():
    result = f.certificate()
    assert result['consecutive_complete_edges'] == 12
    assert result['orthonormal_singlet_chart_dimension'] == 36
    with pytest.raises(ValueError): f.ordered_edges([x for group in cube.GEOMETRY.groups for x in group][:-1])
    terms = [x for group in cube.GEOMETRY.groups for x in group]
    terms[1], terms[2] = terms[2], terms[1]
    with pytest.raises(ValueError): f.ordered_edges(terms)


def test_interval_symmetric_square_step_matches_independent_generator_exponentials():
    gram = spectrum.singlet_data()[2]
    root = np.sqrt(gram)
    direct = np.eye(36, dtype=complex)
    gates = []
    for x, y, _, sign in f.ordered_edges():
        h = root[:, None]*f.symmetric_generator(x, y, sign)/root[None, :]/16
        gate = expm(-1j*float(f.STEP/2)*h); gates.append(gate)
        direct = gate@direct
    contact = np.diag([np.exp(-1j*float(f.STEP)) if x == y else 1 for x, y in f.PAIRS])
    direct = contact@direct
    for gate in reversed(gates): direct = gate@direct
    represented = cp.from_intervals(*f.step_intervals()).diagnostic_matrix()
    np.testing.assert_allclose(represented, direct, atol=2e-15, rtol=0)
    assert np.linalg.norm(represented.conj().T@represented-np.eye(36)) < 2e-14


def test_lifted_step_matches_all_65536_row_native_field_evaluator():
    _, singlets, gram, _, _, q = spectrum.singlet_data()
    t = cube.parent_lift()[0]
    frame = cube.restrict(cube.operators().frame.matrix(), 2)
    chart = frame.conj().T@singlets/np.sqrt(gram)[None, :]
    initial = np.zeros(65536, complex)
    initial[cube.weights(6)] = t@(chart@q)/math.sqrt(60*56)
    space = FullFock(16)
    actual = evolve(cube.GEOMETRY, initial, float(f.STEP), 1, space=space)
    step = cp.from_intervals(*f.step_intervals()).diagnostic_matrix()
    expected = np.zeros(65536, complex)
    expected[cube.weights(6)] = np.exp(-2j*float(f.STEP))*(t@(chart@(step@q)))/math.sqrt(60*56)
    np.testing.assert_allclose(actual, expected, atol=3e-12, rtol=0)


def test_all_eight_times_meet_the_full_native_probability_budget():
    result = budgets.certificate()
    assert len(result['schedule']) == 8
    assert result['schedule'][-1]['total_steps'] == 3670016
    from fractions import Fraction as F
    assert all(F(x) < F(1, 65536) for row in result['schedule'] for x in row['native_systematic_probability_bounds'])
