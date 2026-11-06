"""Mutation-sensitive source-record and quantum-history checks."""
from fractions import Fraction

import numpy as np
import pytest
import sympy as s

from . import certificate as c
from ..common_source import native_family_record_memory as memory


def test_recorded_recovery_all_native_seams_and_zero_slots():
    result = c.recovery_certificate()
    assert result['branches_checked_including_zero'] == 1440
    assert result['positive_branches'] == 600
    assert result['exact_branch_weights'] == ['0', '1/12', '1/3', '5/12']


@pytest.mark.parametrize('images', [
    ((0,), (1,), (2,), (3,), (4,)),
    ((0,), (1,), (2,), (3,), (4,), (4,)),
    ((0,), (1,), (2,), (3,), (4,), (5, 6)),
    ((0,), (1,), (2,), (3,), (4,), (512,)),
])
def test_input_dependent_and_nonisometric_recorded_error_rejected(images):
    with pytest.raises(ValueError):
        c.validate_images(images)


def test_memory_capacity_identity_keeps_every_outcome():
    result = c.compression_certificate()
    assert result['positive_branch_intertwinings'] + result['zero_outcomes_on_admitted_code'] == 1536
    assert result['complete_history_identity']


def test_mutated_native_destination_breaks_capacity_comparison():
    rows = c.branch_destinations()
    for i, (label, w, destinations) in enumerate(rows):
        if w:
            changed = list(destinations); changed[2] = 5
            rows[i] = label, w, tuple(changed)
            break
    with pytest.raises(ValueError, match='source branch'):
        c.compression_certificate(destinations=rows)


def test_deleting_zero_outcomes_is_rejected():
    rows = c.branch_destinations()
    rows = [row for row in rows if row[1]]
    with pytest.raises(ValueError, match='source branch'):
        c.compression_certificate(destinations=rows)


def test_quantum_branches_preserve_entangled_work_reference():
    # Direct matrix test independent of the symbolic classical branch join.
    # Check a complex entangled state at three unrelated pulse values,
    # including one beyond the neutral small-pulse bound.
    rng = np.random.default_rng(24681)
    vec = rng.normal(size=(4, 6, 2)) + 1j*rng.normal(size=(4, 6, 2))
    vec /= np.linalg.norm(vec)
    full = np.zeros((12, 6, 2), complex)
    full[list(c.SUPPORT)] = vec
    data = memory.architecture()
    for h in (.13, .71, 1.27):
        U = memory.feedback_unitaries(h)
        total = 0.
        for _, w, destinations in c.branch_destinations():
            for j, old in enumerate(c.SUPPORT):
                if not w:
                    continue
                dest = destinations[old]
                large_output = np.sqrt(float(w))*U[dest]@full[old]
                small_output = np.sqrt(float(w))*U[c.SUPPORT[c.SUPPORT.index(dest)]]@vec[j]
                assert np.max(abs(large_output-small_output)) < 1e-15
                total += np.linalg.norm(small_output)**2
        assert abs(total-1) < 4e-14


def test_duty_is_not_a_clock_or_coordinate_correction():
    result = c.duty_certificate()
    assert result['policies']['always']['rank_exact'] == 3
    assert result['policies']['C']['rank_exact'] == 3
    assert result['policies']['E']['rank_exact'] == 2
    assert result['policies']['CE']['rank_exact'] == 2
    assert result['policies']['E']['ratios_diagnostic'] is None
    assert result['policies']['CE']['ratios_diagnostic'] is None
    assert (result['policies']['always']['normalized_trace_square_exact'] !=
            result['policies']['C']['normalized_trace_square_exact'])


def test_conditional_normalization_does_not_change_ratios():
    result = c.duty_certificate()['policies']['C']
    matrix = s.Matrix(result['matrix_exact']).applyfunc(s.sympify)
    p = s.Rational(result['trace_exact'])
    v = np.linalg.eigvalsh(np.asarray(matrix, float))
    normalized = np.linalg.eigvalsh(np.asarray(matrix/p, float))
    assert np.max(abs(v/v[0]-normalized/normalized[0])) < 1e-9


def test_complete_observer_inverse_and_mixing_bound():
    result = c.readout_certificate()
    assert Fraction(result['exact_memory_l1_error']) <= Fraction(result['exact_memory_l1_upper'])
    assert result['two_clock_inverse_error_diagnostic'] < result['two_clock_numerical_allowance']
    for value, interval in zip(result['finite_ratios_diagnostic'][1:],
                               result['mixing_only_stationary_ratio_diagnostic_intervals']):
        assert interval[0] <= value <= interval[1]
    assert result['mixing_only_stationary_ratio_diagnostic_intervals'][0][1] < 196


def test_recovery_preserves_visible_poles_and_detector_can_hide_one():
    result = c.pole_certificate()
    assert result['masked_poles'] == ['3']
    assert result['unreachable_carrier_pole'] == '5'
    assert not result['faithful_recovery_moves_logical_poles']


def test_full_source_quantum_channel_intertwines_with_encoded_channel():
    result = c.compression_certificate()
    T4 = np.asarray(s.Matrix(result['transition_exact']).applyfunc(s.sympify), float)
    E = np.zeros((12, 4))
    E[list(c.SUPPORT), np.arange(4)] = 1
    embed = np.kron(E, np.eye(36))
    for h in (.07, .58):
        source = memory.joint_channel(h)
        U = memory.feedback_unitaries(h)
        reduced = np.block([[T4[i, j]*np.kron(U[old], U[old].conj())
                             for j in range(4)] for i, old in enumerate(c.SUPPORT)])
        assert np.max(abs(source@embed-embed@reduced)) < 3e-15
        # A changed pulse keeps classical records but violates quantum equality.
        altered = reduced.copy()
        altered[:36, :36] = T4[0, 0]*np.eye(36)
        assert np.linalg.norm(source@embed-embed@altered) > 1e-3


def test_transient_initial_memory_is_not_an_admitted_encoding():
    result = c.compression_certificate()
    E = s.zeros(12, 4)
    for j, old in enumerate(c.SUPPORT):
        E[old, j] = 1
    transient = s.zeros(12, 1); transient[0] = 1
    assert E*E.T*transient != transient
    assert 'Transient initial memory' in result['scope']


def test_memory_zero_startup_keeps_original_first_record():
    result = c.compression_certificate()
    assert result['one_step_startup_domain'] == [0, 1, 2, 3, 9, 10]
    values = s.Matrix(result['zero_start_memory_after_first_event']).applyfunc(s.sympify)
    assert values == s.Matrix([s.Rational(71, 72), 0, s.Rational(1, 72), 0])
    assert values != s.Matrix(result['stationary_exact']).applyfunc(s.sympify)
    assert values != s.Matrix([1, 0, 0, 0])


def test_raw_probability_mutation_with_same_unread_channel_rejected():
    rows = c.branch_destinations()
    i, j = 0, 4
    assert rows[i][2] == rows[j][2]
    change = s.Rational(5, 288)
    rows[i] = rows[i][0], rows[i][1]+change, rows[i][2]
    rows[j] = rows[j][0], rows[j][1]-change, rows[j][2]
    assert sum(w for _, w, _ in rows) == 1
    with pytest.raises(ValueError, match='source branch'):
        c.compression_certificate(destinations=rows)


def test_reused_inertia_brackets_certify_three_distinct_eigenvalues():
    result = c.readout_certificate()
    roots = [[s.Rational(a), s.Rational(b)] for a, b in result['family_root_enclosures_exact']]
    assert len(roots) == 3 and all(lo < hi for lo, hi in roots)
    assert roots[0][1] < roots[1][0] < roots[1][1] < roots[2][0]
    ratio = result['finite_preparation_ratio_enclosures_exact'][0]
    assert s.Rational(ratio[0]) > 195 and s.Rational(ratio[1]) < 196
