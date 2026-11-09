"""Complete-sector, covariance, cutoff and recorded-current controls."""
from fractions import Fraction as F

import mpmath as mp
import numpy as np
import pytest
import sympy as s

from model.particle_observer_outputs import field_route as f


def test_legendre_polynomial_gap_and_global_controls_are_exact():
    h = s.diag(s.Rational(1, 2), -s.Rational(1, 2), 8, -8)
    q = f.sign_polynomial(h, 32)
    bound = f.polynomial_bounds(32)
    epsilon = s.Rational(bound['sign_error_on_admitted_gap'].numerator,
                         bound['sign_error_on_admitted_gap'].denominator)
    for j, sign in enumerate((1, -1, 1, -1)):
        assert abs(q[j, j]-sign) <= epsilon
    low = f.sign_polynomial(s.diag(0, s.Rational(1, 8), -s.Rational(1, 8)), 12)
    assert low[0, 0] == 0
    assert 0 < low[1, 1] < 1 and -1 < low[2, 2] < 0


def test_polynomial_covariance_and_exact_word_degree():
    x = s.Symbol('x', real=True)
    expression = f.sign_polynomial(s.Matrix([[x]]), 3)[0]
    assert s.Poly(expression, x).degree() == 7
    assert s.expand(expression+expression.subs(x, -x)) == 0
    h = s.Matrix([[1, s.I], [-s.I, -1]])
    gauge = s.Matrix([[1, 1], [-1, 1]])/s.sqrt(2)
    assert f._clean(f.sign_polynomial(gauge*h*gauge.H, 3)
                    -gauge*f.sign_polynomial(h, 3)*gauge.H).is_zero_matrix


def test_cutoff_must_compress_the_polynomial_of_the_parent():
    # Multiplication by cos(theta). The halo retains every intermediate
    # Fourier mode traversed by the degree-five polynomial.
    small, degree = 2, 2
    halo = small+2*degree+1
    def cosine(k):
        n = 2*k+1
        return s.SparseMatrix(n, n, {(i, j): s.Rational(1, 2)
                                    for i in range(n) for j in range(n) if abs(i-j) == 1})
    full = f.sign_polynomial(cosine(halo), degree)
    indices = list(range(halo-small, halo+small+1))
    correct = full.extract(indices, indices)
    naive = f.sign_polynomial(cosine(small), degree)
    assert correct != naive
    larger = f.sign_polynomial(cosine(halo+1), degree)
    assert correct == larger.extract([i+1 for i in indices], [i+1 for i in indices])


def test_rational_dyson_tail_encloses_independent_high_precision_sum():
    with mp.workdps(90):
        for x, first in ((F(1, 4), 3), (F(7, 2), 6), (F(0), 5)):
            bound = f.exponential_tail_upper(x, first)
            xm = mp.mpf(x.numerator)/x.denominator
            actual = mp.exp(xm)-sum(xm**j/mp.factorial(j) for j in range(first))
            ceiling = mp.mpf(bound.numerator)/bound.denominator
            assert actual <= ceiling+mp.mpf('1e-85')
    row = f.cutoff_bound(30, 6, 6, F(1, 4))
    assert row['first_possible_omitted_Dyson_order'] == 5
    assert F(row['complete_record_TV_error_upper']) < F(1, 50000)


def test_native_gamma_and_full_Fock_lift_retain_all_occupations():
    d = f.rotor_coefficients()
    assert d['charge'] == 6
    assert d['hop'].shape == (16, 16)
    assert sorted(d['number'].diagonal()) == [0]+[1]*4+[2]*6+[3]*4+[4]
    assert f._clean(d['number']*d['hop']-d['hop']*d['number']).is_zero_matrix
    assert d['source_provenance']


def test_same_H_state_clock_current_and_complete_Gauss_sector():
    d = f.exact_rotor(6)
    assert d['H'].shape == (208, 208)
    assert d['initial_flux_mean'] == 3
    assert d['initial_current'] == 6
    assert d['initial_energy'] == 11
    assert f._clean(d['Gauss']*d['state']).is_zero_matrix
    assert f._clean(d['Gauss']*d['H']-d['H']*d['Gauss']).is_zero_matrix
    actual = s.I*(d['H']*d['E']-d['E']*d['H'])
    assert f._clean(actual-d['current']).is_zero_matrix
    assert not f._clean(actual+d['current']).is_zero_matrix
    # The charged finite shift has boundary loss. It was never made unitary.
    assert d['shift'].H*d['shift'] != s.eye(13)


def test_complete_two_record_law_includes_every_zero_and_boundary_output():
    row = f.recorded_rotor(6)
    p = np.array(row['all_joint_flux_probabilities'])
    assert p.shape == (13, 13)
    assert row['record_count_including_zero_records'] == 169
    assert row['complete_Gauss_zero_dimension'] == 78
    assert np.count_nonzero(p == 0) > 0
    assert p.min() >= 0
    assert p.sum() == pytest.approx(1, abs=2e-12)
    assert row['first_flux_mean'] > 3
    exact = [[F(x) for x in entries] for entries in row['all_joint_flux_probabilities_exact']]
    norm = sum(map(sum, exact))
    evolution = row['finite_time_evolution']
    assert abs(norm-1) == F(row['normalization_error_exact'])
    assert abs(norm-1) <= 2*F(evolution['complete_record_half_l1_error_upper'])
    assert row['evolution_error_rationally_bounded']
    assert not evolution['floating_evolution_used']
    assert not evolution['record_probabilities_renormalized']


@pytest.mark.parametrize('order', [0, 1, 2, 7])
def test_gaussian_integer_horner_matches_independent_matrix_polynomial(order):
    doubled = s.Matrix([[2, 1+s.I], [1-s.I, -3]])
    rows = tuple(tuple((j, f._gaussian(doubled[i, j])) for j in range(2)) for i in range(2))
    initial = ((2, -1), (0, 3))
    for duration in (F(0), F(1, 3), F(5, 7)):
        numerators, denominator = f._taylor_action(rows, initial, duration, order)
        actual = s.Matrix([s.Rational(a, denominator)+s.I*s.Rational(b, denominator)
                           for a, b in numerators])
        argument = -s.I*s.Rational(duration.numerator, duration.denominator)*doubled/2
        polynomial = sum((argument**j/s.factorial(j) for j in range(order+1)), s.zeros(2))
        expected = polynomial*s.Matrix([2-s.I, 3*s.I])
        assert f._clean(actual-expected).is_zero_matrix


def test_prepared_exact_matrix_is_the_full_parent_restriction():
    parent = f.exact_rotor(6)
    rows, initial, occupations, flux = f._prepared_rotor_rows(6)
    entries = {(i, j): a+s.I*b for i, row in enumerate(rows) for j, (a, b) in row}
    actual = s.SparseMatrix(len(rows), len(rows), entries)
    indices = [mask*13+ell+6 for mask in occupations for ell in flux]
    assert actual == 2*parent['H'].extract(indices, indices)
    vector = s.Matrix([a+s.I*b for a, b in initial])/s.sqrt(2)
    assert vector == parent['state'].extract(indices, [0])
    with pytest.raises(ArithmeticError):
        f._gaussian(s.Rational(1, 2))


def test_unnormalized_records_have_a_certified_truncation_error():
    row = f.recorded_rotor(6, taylor_order=4)
    probabilities = [F(v) for entries in row['all_joint_flux_probabilities_exact'] for v in entries]
    assert sum(probabilities) != 1
    assert all(p >= 0 for p in probabilities)
    assert abs(sum(probabilities)-1) <= F(row['finite_time_evolution']['any_single_event_error_upper'])
    assert not row['finite_time_evolution']['record_probabilities_renormalized']


def test_refined_complete_record_distributions_obey_declared_cutoff_bound():
    report = f.cutoff_comparison()
    assert report['comparison_is_full_record_law']
    assert [x['record_count_including_zero_records'] for x in report['runs']] == [1369, 3721]
    allowance = F(report['sum_of_infinite_cutoff_TV_bounds'])+F(report['sum_of_time_evolution_half_l1_bounds'])
    assert F(report['complete_joint_record_half_l1_difference_exact']) <= allowance
    assert F(report['sum_of_time_evolution_half_l1_bounds']) < F(1, 10**34)


def test_cubic_bound_is_full_Fock_and_not_a_chiral_claim():
    row = f.cubic_program_bound(2, 8, 200, 1, '1/1024')
    assert row['spatial_dimension'] == 3
    assert row['Dirac_CAR_modes'] == 480
    assert row['full_Fock_dimension'] == str(2**480)
    assert row['one_body_bounds']['maximum_coordinate_word_degree'] == 17
    assert row['both_chiralities_retained']
    assert not row['physical_chiral_continuum_established']


@pytest.mark.parametrize('args', [(3, 4, 1, '1'), (3, 0, 0, '1'), (3, 0, 1, 0.2), (3, 0, 1, '-1')])
def test_cutoff_admission_rejects_incomplete_or_inexact_inputs(args):
    with pytest.raises(ValueError):
        f.cutoff_bound(*args)
