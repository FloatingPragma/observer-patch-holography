from fractions import Fraction as F

import pytest
import sympy as s

from ..cube_eta_protection import spectrum
from . import vacuum


def test_degree_exhausting_exact_root_signs():
    cert = vacuum.certificate()
    assert cert['negative_frequency_lines'] == 5
    assert cert['positive_frequency_lines'] == 3
    assert cert['zero_frequency_lines'] == 0
    for row in cert['spectral_rows']:
        left, right = map(F, row['endpoint_polynomial_values'])
        assert left*right < 0
        omega = tuple(map(F, row['omega_interval']))
        assert omega[1] < 0 if row['frequency_sign'] < 0 else omega[0] > 0


def test_exact_sturm_counts_independent_of_bisection():
    z = spectrum.Z
    for _, coefficients, brackets in vacuum.FACTORS:
        polynomial = s.Poly.from_list(coefficients, z)
        assert sum(polynomial.count_roots(left, right) for left, right in brackets) == polynomial.degree()
        assert all(polynomial.count_roots(left, right) == 1 for left, right in brackets)


def test_rational_residue_intervals_cover_independent_algebraic_roots():
    numerator, denominator, _ = spectrum.probe_resolvent()
    rows = vacuum.certificate()['spectral_rows']
    for row in rows:
        coeffs = next(c for name, c, _ in vacuum.FACTORS if name == row['factor'])
        index = next(i for name, _, brackets in vacuum.FACTORS if name == row['factor']
                     for i, bracket in enumerate(brackets) if list(bracket) == row['coarse_z_bracket'])
        root = s.CRootOf(s.Poly.from_list(coeffs, spectrum.Z), index)
        value = (numerator/s.diff(denominator, spectrum.Z)).subs(spectrum.Z, root).evalf(50)
        lo, hi = [s.Rational(x) for x in row['positive_normalized_probe_weight_interval']]
        assert lo < value < hi


def test_negative_mass_is_positive_measure_not_retarded_sign():
    cert = vacuum.certificate()
    low, high = map(F, cert['normalized_probe_negative_mass_interval'])
    assert F(1, 40) < low < high < F(13, 500)
    assert all(F(r['positive_normalized_probe_weight_interval'][0]) > 0 for r in cert['spectral_rows'])
    assert cert['measure_is_retarded_commutator'] is False
    assert 'ground-state retarded' in cert['retarded_caution']


def test_closed_quadratic_residue_and_lowering_energy():
    numerator, denominator, _ = spectrum.probe_resolvent()
    z = spectrum.Z
    root = 8-2*s.sqrt(19)
    residue = s.simplify((numerator/s.diff(denominator, z)).subs(z, root))
    assert s.simplify(residue-(1-4/s.sqrt(19))/14) == 0
    assert s.Rational(1, 174) < residue < s.Rational(9, 1526)
    energy = 2+root/16
    assert energy < 2 < 3


def test_positive_energy_scale_and_same_charge_shift_cannot_change_signs():
    alpha = s.Symbol('alpha', positive=True)
    mu, shift = s.symbols('mu shift', real=True)
    lower_energy = s.Rational(5, 2)-s.sqrt(19)/8
    converted = (alpha*lower_energy+shift-6*mu)-(3*alpha+shift-6*mu)
    assert s.simplify(converted/alpha) == -s.Rational(1, 2)-s.sqrt(19)/8


def test_witness_does_not_claim_local_F_extracts_work_or_theory_invalid():
    cert = vacuum.certificate()
    witness = cert['lowering_witness']
    assert witness['executed_local_F_mean_energy_change'] == '0'
    assert witness['native_extraction_protocol_executed'] is False
    assert cert['other_vacuum_candidates_excluded'] is False
    assert cert['excited_reference_or_thermal_response_inconsistent'] is False


@pytest.mark.parametrize('steps', [0, True, 257])
def test_bad_root_bounds_rejected(steps):
    with pytest.raises(ValueError): vacuum.root_interval((1, -16, -12), (-1, 0), steps)


def test_interval_polynomial_outward_and_ratio_zero_guard():
    assert vacuum.interval_polynomial((1, 0, -1), (F(-2), F(3))) == (F(-7), F(8))
    with pytest.raises(ValueError): vacuum.ratio_interval((F(1), F(2)), (F(-1), F(1)))
