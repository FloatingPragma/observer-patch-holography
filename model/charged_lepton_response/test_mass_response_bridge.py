from fractions import Fraction as F
import json

import numpy as np
import pytest
import sympy as s

from . import mass_response_bridge as bridge
from model.standard_model_source_map import koide_source


@pytest.mark.parametrize('value', [F(0), F(1), F(2), F(1, 3), F(10**40, 7), F(1, 10**40)])
def test_dyadic_square_roots_are_certified_by_integer_squares(value):
    lo, hi = bridge.sqrt_enclosure(value, bits=96)
    assert lo * lo <= value <= hi * hi
    assert hi - lo <= F(1, 2**96)
    if value in (0, 1):
        assert lo == hi == value


@pytest.mark.parametrize('args', [(1, 0, 0), (0, -1, 0), (0, 0, 1), (True, 0, 0),
                                  (.1, 0, 0)])
def test_relative_admission_rejects_missing_positive_kinetic_or_mass_control(args):
    with pytest.raises(ValueError):
        bridge.coefficient_bound(*args)


def test_unresolved_dyadic_denominator_is_rejected():
    with pytest.raises(ValueError, match='precision'):
        bridge.coefficient_bound(F(1) - F(1, 2**20), 0, 0, bits=2)


def test_identity_and_common_scaling_do_not_change_ratios():
    bound = bridge.coefficient_bound(0, 0, 0)
    assert bound['effective_lower_multiplier'] == bound['effective_upper_multiplier'] == 1
    assert bridge.koide_interval(F(2, 3), F(2, 3), F(7), F(7))['lower'] == F(2, 3)
    assert bridge.koide_interval(F(2, 3), F(2, 3), F(7), F(7))['upper'] == F(2, 3)
    gradient = bridge.koide_log_gradient([1, 4, 9])
    assert sum(gradient) == 0
    e = s.symbols('e', positive=True)
    q = lambda m: sum(m)/sum(s.sqrt(x) for x in m)**2
    assert s.simplify(q([(1+e)*x for x in (1, 4, 9)])-q([1, 4, 9])) == 0


@pytest.mark.parametrize('seed', [7, 19, 41, 113])
def test_both_noncommuting_kinetic_forms_bound_nonnormal_mass_singular_values(seed):
    rng = np.random.default_rng(seed)
    unitary = lambda: np.linalg.qr(rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3)))[0]
    ul, ur, um, vm = [unitary() for _ in range(4)]
    eps_l, eps_r, eta = F(1, 5), F(1, 7), F(1, 20)
    z_l = ul @ np.diag([1-float(eps_l), 1.02, 1+float(eps_l)]) @ ul.conj().T
    z_r = ur @ np.diag([1-float(eps_r), .99, 1+float(eps_r)]) @ ur.conj().T
    m0 = um @ np.diag([1, 3, 9]) @ vm.conj().T
    e = float(eta) * np.outer(ul[:, 0], ur[:, 1].conj())
    inverse_sqrt = lambda z: (np.linalg.eigh(z)[1] / np.sqrt(np.linalg.eigh(z)[0])) @ np.linalg.eigh(z)[1].conj().T
    actual = np.sort(np.linalg.svd(inverse_sqrt(z_l) @ (m0+e) @ inverse_sqrt(z_r), compute_uv=False))
    assert np.linalg.norm(z_l @ m0-m0 @ z_l) > .1
    assert np.linalg.norm(e-e.conj().T) > .01
    b = bridge.coefficient_bound(eps_l, eps_r, eta)
    assert np.all(actual >= float(b['effective_lower_multiplier'])*np.array([1, 3, 9])-1e-12)
    assert np.all(actual <= float(b['effective_upper_multiplier'])*np.array([1, 3, 9])+1e-12)


def test_exact_saturation_requires_both_kinetic_factors():
    eps, eta = F(1, 4), F(1, 8)
    b = bridge.coefficient_bound(eps, eps, eta)
    assert b['effective_lower_multiplier'] == (1-eta)/(1+eps)
    assert b['effective_upper_multiplier'] == (1+eta)/(1-eps)
    m0 = s.diag(1, 1, 1)
    for sign, key in ((1, 'effective_upper_multiplier'), (-1, 'effective_lower_multiplier')):
        actual = (1+sign*s.Rational(eta))/(1-sign*s.Rational(eps))*m0
        assert actual[0, 0] == s.Rational(b[key])


def test_absolute_norm_error_keeps_zero_and_does_not_claim_relative_precision():
    rows = bridge.absolute_singular_intervals([(0, 0), (F(1, 10), F(1, 10)), (2, 3)], 0, 0, F(1, 5))
    assert rows == [(F(0), F(1, 5)), (F(0), F(3, 10)), (F(9, 5), F(16, 5))]
    for args in (([(2, 1)], 0, 0, 0), ([], 0, 0, 0), ([(1, 2)], 1, 0, 0)):
        with pytest.raises(ValueError):
            bridge.absolute_singular_intervals(*args)


def test_separate_matching_error_is_charged_exactly_once_and_not_certified():
    b = bridge.coefficient_bound(0, 0, F(1, 10), matching_error=F(1, 5))
    assert b['matched_lower_multiplier'] == F(18, 25)
    assert b['matched_upper_multiplier'] == F(33, 25)
    assert b['actual_physical_matching_bound_established'] is False


@pytest.mark.parametrize('factors', [(F(9, 10), F(11, 10)), (F(1, 100), F(100)), (F(2), F(3))])
def test_Koide_interval_contains_all_endpoint_perturbations(factors):
    from itertools import product
    a, b = factors
    m = np.array([1., 4., 9.])
    q0 = F(14, 36)
    bounds = bridge.koide_interval(q0, q0, a, b)
    for f in product((a, b), repeat=3):
        actual = m * np.array(list(map(float, f)))
        q = actual.sum()/np.sqrt(actual).sum()**2
        assert float(bounds['lower'])-1e-14 <= q <= float(bounds['upper'])+1e-14


def test_source_control_uses_actual_C_squared_and_preserves_its_original_scale():
    c = bridge.source_nonpreservation()
    actual_source = koide_source.exact_source(1)
    assert c['coefficient'] == actual_source['action_coefficient']
    assert s.trace(c['coefficient']) == s.Rational(5, 18)
    values = c['source_singular_values']
    x = s.symbols('x')
    assert s.simplify(c['coefficient'].charpoly(x).as_expr()-s.prod(x-v for v in values)) == 0
    assert s.simplify(c['corrected_coefficient'].charpoly(x).as_expr()-s.prod(x-v for v in c['corrected_singular_values'])) == 0
    assert c['strict_Koide_increase'].is_positive
    assert c['relative_mass_correction_norm'] < 1
    assert c['perturbation_is_an_actual_computed_self_energy'] is False


def test_source_whitening_alone_can_spoil_Koide_without_mass_addition():
    data = bridge.source_nonpreservation()
    p, y = data['heavy_projector'], data['coefficient']
    scale = s.Rational(10001, 10000)
    left_sqrt = s.eye(3)+(scale**2-1)*p
    z_left = left_sqrt**2
    whitened = (left_sqrt.inv()*y).applyfunc(s.simplify)
    light, _, heavy = data['source_singular_values']
    assert (whitened-light*(s.eye(3)-p)-heavy*p/scale**2).applyfunc(s.simplify) == s.zeros(3)
    z = 4+3*s.sqrt(2)
    q = s.simplify((2+(z/scale)**2)/(2+z/scale)**2)
    assert s.simplify(s.Rational(2, 3)-q).is_positive
    assert (z_left-s.eye(3)).rank() == 1


def test_exact_Koide_gradient_matches_direct_nonuniform_variation():
    e = s.symbols('e', real=True)
    masses = [s.Integer(1), s.Integer(4), s.Integer(9)]
    rates = [s.Integer(2), -s.Integer(1), s.Integer(3)]
    changed = [m*(1+e*d) for m, d in zip(masses, rates)]
    q = sum(changed)/sum(s.sqrt(x) for x in changed)**2
    exact = s.diff(q, e).subs(e, 0)/q.subs(e, 0)
    predicted = sum(x*d for x, d in zip(bridge.koide_log_gradient(masses), rates))
    assert s.simplify(exact-predicted) == 0 and exact != 0


def test_certificate_is_JSON_safe_fresh_and_has_no_physical_inputs():
    first = bridge.certificate()
    expected = json.dumps(first, sort_keys=True, allow_nan=False)
    assert len(expected) < 20000
    first['conditional_control']['effective_lower_multiplier'] = 'forged'
    assert json.dumps(bridge.certificate(), sort_keys=True, allow_nan=False) == expected
    assert first['scope']['physical_particle_assignments'] == []
    assert first['scope']['actual_source_loop_correction_computed'] is False


def test_raw_fraction_outputs_do_not_alias_supplied_fraction():
    source = F(1, 10)
    result = bridge.coefficient_bound(source, 0, 0)
    result['kinetic_left_bound']._numerator = 7
    assert source == F(1, 10)
