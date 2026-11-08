"""Complete charged records, exact source frame and analytic error budgets."""
from fractions import Fraction as F

import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import expm

from . import arctan_bands as a


def test_joint_band_preserves_original_weights_and_both_complements():
    row = a.continuum_band()
    assert row['weak_commutator_lower'] == F(1, 440)
    assert row['scaled_energy_band'] == [F(3, 387200), F(387200)]
    assert row['each_original_charge_probability_lower'] == F(1, 387200)
    lo, hi = row['scaled_energy_band']
    d = row['weak_commutator_lower']
    assert lo/12+F(1, 4)/hi == d*d/4
    assert row['complete_preparation_records'] == ['charge_minus', 'charge_plus', 'complement']
    assert row['joint_positive_field_clock_used']
    assert not row['old_zero_field_clock_substituted']
    assert not row['original_odd_CAR_identification']
    assert not row['sharp_shell_or_physical_mass_inferred']


@pytest.mark.parametrize('bad', [0, -1, F(3, 4), 0.1, True])
def test_order_bound_rejects_unphysical_or_inexact_inputs(bad):
    with pytest.raises((TypeError, ValueError)):
        a.continuum_band(bad)


def test_rational_scalar_enclosures_have_proved_endpoints():
    with mp.workdps(100):
        number = lambda f: mp.mpf(f.numerator)/f.denominator
        for x in [F(), F(1, 239), F(1, 16), F(1, 5), F(1)]:
            for terms in [1, 2, 9, 32]:
                lo, hi = a.atan_interval(x, terms)
                assert number(lo) <= mp.atan(number(x)) <= number(hi)
        lo, hi = a.pi_interval()
        assert number(lo) < mp.pi < number(hi)
    for q in [F(), F(1, 3), F(257), F(4)]:
        lo, hi = a.sqrt_interval(q)
        assert lo*lo <= q <= hi*hi


def test_finite_longitudinal_error_uses_actual_joint_scaling():
    n = 16
    row = a.longitudinal_errors(n, n**6, n**6)
    assert row['longitudinal_variance_upper'] == F(2, n)
    assert row['longitudinal_configuration_commutator_upper'] == F(1, 2*n**4)
    better = a.longitudinal_errors(2*n, (2*n)**6, (2*n)**6)
    assert better['longitudinal_variance_upper'] == row['longitudinal_variance_upper']/2
    assert better['longitudinal_configuration_commutator_upper'] == row['longitudinal_configuration_commutator_upper']/16
    floor = a.finite_commutator_lower(F(1, 10), F(1, 10**8), F(1, 10**8), 6)
    assert 0 < floor['weak_commutator_lower'] < F(1, 440)
    assert a.finite_commutator_lower(F(1, 10), 1, 1, 6)['weak_commutator_lower'] == 0


@pytest.mark.parametrize('call', [
    lambda: a.longitudinal_errors(1, 1, 1),
    lambda: a.longitudinal_errors(2, -1, 1),
    lambda: a.longitudinal_errors(2, 1, 1, 0),
    lambda: a.atan_interval(F(2)),
    lambda: a.atan_interval(F(1, 2), 0),
    lambda: a.sqrt_interval(-1),
    lambda: a.sqrt_interval(1, True),
])
def test_error_inputs_do_not_accept_invalid_domains(call):
    with pytest.raises((TypeError, ValueError)):
        call()


def test_actual_source_clock_charge_and_nine_record_identity():
    row = a.finite_source_control()
    assert row['complete_record_count'] == 9
    assert row['allocated_zero_probability_records'] == 5
    assert row['all_input_effect_sum_identity']
    assert row['whole_source_and_charge_rotated_together']
    assert not row['finite_control_has_continuum_momentum_assignment']
    assert not row['additional_native_pulse_tape_executed']


def test_native_source_mutation_is_rejected(monkeypatch):
    original = a.source.matrix
    def changed(*args, **kwargs):
        matrix, scale = original(*args, **kwargs)
        matrix = matrix.copy()
        matrix[0, 0] += 1
        return matrix, scale
    monkeypatch.setattr(a.source, 'matrix', changed)
    with pytest.raises(ArithmeticError, match='actual source'):
        a.finite_source_control()


def test_all_nine_rational_intervals_cover_independent_high_precision_probabilities():
    row = a.dimer_intervals()
    with mp.workdps(100):
        number = lambda f: mp.mpf(f.numerator)/f.denominator
        root = mp.sqrt(257)
        k = mp.atan(mp.mpf(1)/16)/(2*mp.pi)
        s = (1+16/root)/2
        c = mp.sqrt(1-4*k*k)
        exact = {
            ('charge_minus', 'middle'): 2*k*k*s,
            ('charge_plus', 'middle'): 2*k*k*s,
            ('complement', 'ground'): (1-s*(1-c))**2,
            ('complement', 'upper'): s*(1-s)*(1-c)**2}
        assert mp.almosteq(sum(exact.values()), 1)
        for item in row['records']:
            truth = exact.get((item['preparation'], item['energy_label']), mp.mpf(0))
            lo, hi = item['original_probability_interval']
            assert number(lo) <= truth <= number(hi)
        assert sum(item['original_probability_interval'] == [0, 0]
                   for item in row['records']) == 5
        lo, hi = row['charged_scaled_energy_interval']
        assert number(lo) <= (16+root)/8 <= number(hi)
    assert not row['preparation_success_probability_divided_out']


def test_finite_source_duhamel_error_on_actual_ground_and_all_angle_cases():
    h, scale = a.source.matrix(2, ((0, 1),), F(1), F(1, 32))
    pauli_x = np.array([[0, 1], [1, 0]], complex)
    pauli_y = np.array([[0, -1j], [1j, 0]], complex)
    pauli_z = np.diag([1, -1]).astype(complex)
    eye = np.eye(2)
    marshall = np.diag([1, 1, -1, -1])
    hh = np.array([[1, 1], [1, -1]], complex)/np.sqrt(2)
    rotate = np.kron(hh, hh)
    hb = rotate@marshall@h@marshall@rotate/scale
    _, vectors = np.linalg.eigh(hb)
    ground = vectors[:, 0]
    p = 0.25
    sx0, sx1 = np.kron(eye, pauli_x)/2, np.kron(pauli_x, eye)/2
    sy0, sy1 = np.kron(eye, pauli_y)/2, np.kron(pauli_y, eye)/2
    sz0, sz1 = np.kron(eye, pauli_z)/2, np.kron(pauli_z, eye)/2
    x, y, momentum, m = p*p*(sx0-sx1), -p*p*(sy0-sy1), p*(sy0+sy1), p**3*(sz0-sz1)
    expect = lambda op: np.vdot(ground, op@ground)
    mean = expect(m).real
    sigma = np.linalg.norm((m-mean*np.eye(4))@ground)
    c = p**5
    assert abs(mean-p**3/np.sqrt(257)) < 1e-15
    assert sigma*sigma <= 8*p**6/(1/32)
    for theta in [0, .13, .8, 1.7]:
        xr = np.cos(theta)*x+np.sin(theta)*y
        assert np.linalg.norm(m@xr-xr@m, 2) <= c+1e-15
        for u in [-5, -.4, .2, 3]:
            unitary = expm(1j*u*xr)
            difference = expect(momentum@unitary-unitary@momentum)-u*np.cos(theta)*mean*expect(unitary)
            bound = abs(u*np.cos(theta))*(sigma+abs(u)*c/2)
            assert abs(difference) <= bound+1e-15


def test_result_mutation_does_not_change_fresh_analytic_certificate():
    first = a.continuum_band()
    first['order_floor']._numerator = 999
    first['complete_preparation_records'].clear()
    second = a.continuum_band()
    assert second['order_floor'] == F(1, 10)
    assert len(second['complete_preparation_records']) == 3
