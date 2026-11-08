from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from model.native_spin_observer import circuits
from model.spin_joint_field_limit import clock_comparison as clock


def _source():
    ops = clock._spin_operators(2)
    signs = (1, -1)
    h0 = sum((ops[0][a]@ops[1][a] for a in range(3)), np.zeros((4, 4), complex))-np.eye(4)/4
    b, p = .25, .25
    d = -b*sum((signs[x]*ops[x][0] for x in range(2)), np.zeros((4, 4), complex))/p
    q = p*p*sum((signs[x]*(x+1)*ops[x][2] for x in range(2)), np.zeros((4, 4), complex))
    momentum = -p*sum(((x+1)*ops[x][1] for x in range(2)), np.zeros((4, 4), complex))
    return h0/p, d, q, momentum, b, p


def test_full_actual_source_drift_and_wrong_sign_control():
    rows = clock.source_rotation_certificate()['finite_source_checks']
    assert [row['all_rows'] for row in rows] == [2, 4, 8]
    _, d, q, momentum, b, p = _source()
    assert np.array_equal(1j*(d@momentum-momentum@d), -b/p**2*q)
    assert np.linalg.norm(1j*(d@momentum-momentum@d)-b/p**2*q) > 0


@pytest.mark.parametrize('time', [F(-3, 2), F(1, 7), F(3)])
def test_exact_selection_rotation_and_input_vector_bounds(time):
    _, d, q, momentum, b, p = _source()
    u = expm(1j*float(time)*d)
    theta = b*float(time)/p
    qrot = np.cos(theta)*q+p*np.sin(theta)*momentum
    prot = np.cos(theta)*momentum-np.sin(theta)/p*q
    assert np.linalg.norm(u@q@u.conj().T-qrot) < 1e-12
    assert np.linalg.norm(u@momentum@u.conj().T-prot) < 1e-12
    vector = np.array([1, 1j, 2, -1j], complex)/np.sqrt(7)
    qnorm, pnorm = np.linalg.norm(q@vector), np.linalg.norm(momentum@vector)
    c = p**3*(1+4)/2
    qbound = (1-np.cos(theta))*qnorm+p*abs(np.sin(theta))*(pnorm+c/2)
    pbound = (1-np.cos(theta))*pnorm+abs(np.sin(theta))/p*(qnorm+c/2)
    assert np.linalg.norm((expm(1j*qrot)-expm(1j*q))@vector) <= qbound+1e-14
    assert np.linalg.norm((expm(1j*prot)-expm(1j*momentum))@vector) <= pbound+1e-14


def test_joint_rotation_error_rates_do_not_assert_mixed_history_equivalence():
    rows = [clock.selection_bounds(n, F(1), F(1), F(2), F(3)) for n in (2, 4, 8, 16)]
    assert all(a['P_gate_input_vector_error_upper'] > b['P_gate_input_vector_error_upper']
               for a, b in zip(rows, rows[1:]))
    assert all(a['Q_gate_input_vector_error_upper'] > b['Q_gate_input_vector_error_upper']
               for a, b in zip(rows, rows[1:]))
    assert all(not row['mixed_time_history_comparison_claimed'] for row in rows)
    for n, row in zip((2, 4, 8, 16), rows):
        assert row['selection_angle_absolute_upper'] == F(1, n**3)


def test_actual_mixed_history_duhamel_identity_and_nonzero_torque():
    h0, d, q, _, _, _ = _source()
    hb = h0+d
    a, t = expm(1j*q), .4
    def alpha(h, duration, value):
        u = expm(1j*duration*h)
        return u@value@u.conj().T
    actual = alpha(hb, t, a)-alpha(h0, t, a)
    nodes, weights = np.polynomial.legendre.leggauss(24)
    integral = np.zeros_like(actual)
    for node, weight in zip(nodes, weights):
        s = t*(node+1)/2
        evolved = alpha(hb, s, a)
        integral += weight*t/2*alpha(h0, t-s, 1j*(d@evolved-evolved@d))
    assert np.linalg.norm(actual-integral) < 1e-12
    assert np.linalg.norm(actual) > 1e-3
    assert np.linalg.norm(d@alpha(hb, t, a)-alpha(hb, t, a)@d
                          -(d@a-a@d)) > 1e-3


def test_count_only_native_clock_budget_on_huge_cube():
    row = clock.cube_clock_plan(10**8, 32, (F(3), F(-2), F()), F(1, 1000))
    assert row['sites'] == (2*10**8+1)**3
    assert row['edges'] == 3*(2*10**8)*(2*10**8+1)**2
    assert row['scale'] == F(1, 32**2)
    assert row['selection_field'] == F(1, 32**5)
    assert row['nonzero_clock_slots'] == 2
    assert row['product_operator_error_upper'] <= F(1, 2000)
    assert row['complete_native_isometry_error_upper'] <= F(1, 2000)
    assert row['complete_clock_half_trace_error_upper'] <= F(1, 1000)
    assert row['source_scalar_exchange_phase_retained']
    assert not row['H_B_and_H_0_histories_identified']
    assert not row['clock_sequence_executed']


def test_zero_time_and_invalid_count_controls():
    row = clock.clock_plan(1, 0, F(1), F(1, 4), F(1, 2), (F(),), F(1, 100))
    assert row['native_pulse_count_upper'] == row['product_slices_per_nonzero_slot'] == 0
    assert row['complete_clock_half_trace_error_upper'] == 0
    for arguments in ((0, 0), (1, -1), (1, 4)):
        with pytest.raises(ValueError):
            clock.clock_plan(*arguments, F(1), F(1, 4), F(1, 2), (F(1),), F(1, 100))
    with pytest.raises(ValueError):
        clock.selection_bounds(0, F(1), F(1), F(1), F(1))
    with pytest.raises(ValueError):
        clock.selection_bounds(2, F(1), F(-1), F(1), F(1))
    with pytest.raises(ValueError):
        clock.clock_plan(1, 0, F(1), F(1, 4), F(1, 2), (.5,), F(1, 100))


def test_literal_clock_words_include_each_phase_and_inverse_slot():
    arguments = (2, ((0, 1),), F(1), F(1, 4), F(1, 2))
    positive = list(clock.clock_word(*arguments, F(1, 16), 2))
    negative = list(clock.clock_word(*arguments, F(-1, 16), 2))
    assert len(positive) == len(negative) == 2*(62+11*2)
    assert len(list(circuits.emit(positive))) == len(positive)
    assert list(clock.clock_word(*arguments, F(), 0)) == []
    assert len(list(clock.dyadic_clock_word(*arguments, F(1, 16), 2, 60))) == len(positive)
    with pytest.raises(ValueError):
        list(clock.clock_word(*arguments, F(1, 16), 0))
    with pytest.raises(ValueError):
        list(clock.clock_word(3, ((0, 1), (1, 2), (0, 2)), F(1), F(1, 4), F(1, 2), F(1), 2))


def test_complete_native_clock_reference_and_failure_records():
    row = clock.small_native_clock()
    assert row['pulses'] == 168
    assert row['failure_slots'] == 600
    assert row['all_work_rows'] == 16
    assert row['reference_dimension'] == 4
    assert row['all_failure_records_retained']
    assert row['exact_H_B_clock_distance_within_sum_bound']
    assert not row['large_joint_family_schedule_executed']
