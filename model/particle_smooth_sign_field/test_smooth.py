from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path

import pytest
import sympy as s

from model.particle_local_gauge_state import state
from model.particle_local_history_family.history import cutoff_plan
from . import smooth as sm


def test_exact_weight_contraction_including_every_coefficient():
    u = s.symbols('u')
    remainder = s.Poly(s.expand((1-u)*(1+u/8)**4-(1-u/2)), u)
    assert remainder.degree() == 5
    assert all(c <= 0 for c in remainder.all_coeffs())
    assert remainder.nth(0) == remainder.nth(1) == 0
    for delta in ('1/100', '1', '20'):
        data = sm.settings(delta)
        assert 0 < data['weighted_series_ratio'] <= 1-data['one_minus_ratio']/2 < 1
        assert data['spatial_weight_base'] > 1


def test_star_count_and_generating_function_majorant():
    r = s.symbols('r')
    series = 1/(1-r)
    for _ in range(6):
        series = s.diff(r*series, r)
    numerator = 1+57*r+302*r**2+302*r**3+57*r**4+r**5
    assert s.cancel(series-numerator/(1-r)**7) == 0
    assert numerator.subs(r, 1) == 720
    for n in range(30):
        d = 2*n+1
        assert (2*d+1)**6 <= 4096*(n+1)**6


def test_scalar_series_preserves_sign_and_has_uniform_spectral_bound():
    # This tests actual coefficient calculus at both endpoints, at a
    # kernel zero and inside the spectral interval, without eigen fitting.
    for delta in ('1/2', '3'):
        ds = s.Rational(delta)
        for x in (-8, -1, 0, 1, 8):
            values = [sm.binomial_polynomial([[x]], m, delta)[0] for m in range(4)]
            target = s.sqrt(1+ds**2)*x/s.sqrt(x*x+ds**2)
            for value in values:
                assert s.simplify(target**2-value**2).is_nonnegative is True
            for earlier, later in zip(values, values[1:]):
                assert s.simplify(later**2-earlier**2).is_nonnegative is True
            assert values[0] == s.sqrt(1+ds**2)*x/s.sqrt(64+ds**2)


def test_polynomial_covariance_retains_whole_matrix_and_null_columns():
    kernel = s.diag(-2, 0, 1)
    change = s.Matrix([[s.Rational(3, 5), 0, s.Rational(4, 5)],
                       [0, 1, 0], [-s.Rational(4, 5), 0, s.Rational(3, 5)]])
    original = sm.binomial_polynomial(kernel, 3, '2')
    transformed = sm.binomial_polynomial(change*kernel*change.H, 3, '2')
    assert (transformed-change*original*change.H).applyfunc(s.simplify) == s.zeros(3)
    assert original == original.H and original[:, 1] == s.zeros(3, 1)
    assert original.shape == transformed.shape == (3, 3)


def test_weighted_tail_encloses_actual_finite_path_majorant():
    cfg = sm.settings('8', plaquette='2/3')
    lam, q = cfg['spatial_weight_base'], cfg['binomial_operator_ratio']
    C = cfg['rooted_nth_term_norm_prefactor']
    for N in (0, 1, 3, 20):
        partial = sum(C*(4*n+3)**6*lam**(4*n+2)*q**n for n in range(N, N+30))
        assert partial <= sm.finite_tail_bound(cfg, N)
    assert cfg['weighted_full_interaction_norm_upper']-cfg['weighted_kinetic_interaction_norm_upper'] == 48*F(2, 3)*lam**2


def test_complete_record_plan_certifies_fixed_time_and_instrument_cost():
    cfg = sm.settings('8', plaquette='1/100')
    base = sm.local_record_plan(cfg, '1', 18)
    with_records = sm.local_record_plan(cfg, '1', 18, observed_cells=3, intermediate_instruments=4)
    assert base['duration'] == with_records['duration'] == 1
    assert base['error_binary_exponent_upper'] <= -18
    assert with_records['error_binary_exponent_upper'] <= -18
    assert with_records['binomial_polynomial_degree'] >= base['binomial_polynomial_degree']
    assert with_records['untouched_references_and_all_outcomes_retained']
    assert not with_records['depends_on_total_volume']
    assert not with_records['link_derivative_or_Wilson_gap_bound_required']
    # Independently check both binary inequalities on modest rational inputs.
    q = cfg['weighted_series_ratio']; u = cfg['one_minus_ratio']
    for n in range(1, 40):
        assert q**n <= F(1, 2**int(u*n/2))


def test_open_box_keeps_observation_geometry_and_generic_word_cutoff():
    cfg = sm.settings('8', plaquette='1/100')
    plan = sm.local_record_plan(cfg, '1/10', 10, observed_cells=2, intermediate_instruments=1)
    with pytest.raises(ValueError):
        sm.finite_volume_plan(cfg, plan, 10)
    box = sm.finite_volume_plan(cfg, plan, 10, observed_radius=1)
    d = plan['maximum_coordinate_word_degree']
    assert box['boundary_separation'] == box['box_radius']-1-2*max(d, 2)
    assert box['error_binary_exponent_upper'] <= -10
    assert box['periodic_side_sufficient'] > 2*(box['box_radius']+2*max(d, 2)+d)
    code = cutoff_plan({'Legendre_degree': plan['binomial_polynomial_degree'],
                        'single_link': {'per_link_word_degree': 2},
                        'full_Fock_kinetic_norm_per_cell_upper': cfg['full_Fock_kinetic_density_rational_upper'],
                        'declared_plaquette_coefficient': cfg['plaquette_coefficient']}, box, 10)
    assert code['word_step'] == max(4, d)
    assert code['evolution_intervals'] == 2
    assert code['whole_history_error_binary_exponent_upper'] <= -10


def test_free_corner_census_has_one_zero_and_retains_doubler_energy():
    for delta in ('1/5', '1', '9'):
        energies = [(corner, sm.free_symbol_squared(corner, delta)) for corner in product((-1, 1), repeat=3)]
        assert [corner for corner, energy in energies if energy == 0] == [(1, 1, 1)]
        assert all(s.simplify(energy-4).is_nonnegative is True
                   for corner, energy in energies if corner != (1, 1, 1))


def test_free_axis_velocity_and_wilson_design_endpoint():
    k, delta = s.symbols('k delta', real=True, positive=True)
    w = -s.cos(k)
    R2 = s.simplify(w*w+s.sin(k)**2)
    assert R2 == 1
    t = s.sqrt(1+delta**2)/s.sqrt(R2+delta**2)
    assert s.simplify(t) == 1
    energy2 = s.trigsimp((1+t*w)**2+t*t*s.sin(k)**2)
    assert s.limit(energy2/k**2, k, 0) == 1
    x = s.symbols('x', real=True)
    assert s.limit(s.sqrt(1+delta**2)*x/s.sqrt(x*x+delta**2), delta, s.oo) == x
    # Coefficients of the complete three-dimensional kernel-square identity.
    u1, u2, u3 = s.symbols('u1 u2 u3')
    Rall = (u1+u2+u3-1)**2+sum(2*u-u*u for u in (u1, u2, u3))
    assert s.expand(Rall-(1+2*(u1*u2+u1*u3+u2*u3))) == 0


def test_actual_native_sixty_column_control():
    control = sm.exact_native_control()
    assert control['full_spectral_ranks'] == [16, 28, 16]
    assert sum(control['full_spectral_ranks']) == control['native_one_body_dimension'] == 60
    assert control['smooth_target_annihilates_all16_Wilson_zero_modes']
    assert control['normalized_smooth_squared_eigenvalues'] == [0, 1, s.Rational(3, 2)]
    assert not control['point_background_claimed_a_quantum_preparation']


def test_invalid_parameters_and_measurement_domains_rejected():
    for delta in ('0', '-1', 0.2):
        with pytest.raises(ValueError):
            sm.settings(delta)
    with pytest.raises(ValueError):
        sm.settings(plaquette='-1')
    with pytest.raises(ValueError):
        sm.binomial_polynomial([[1, 2], [0, 1]], 2)
    with pytest.raises(ValueError):
        sm.local_record_plan(sm.settings(), duration=0)
    with pytest.raises(ValueError):
        sm.free_symbol_squared((2, 1, 1))
    with pytest.raises(ValueError):
        sm.finite_tail_bound(sm.settings(), 10001)


def test_certificate_preserves_same_action_and_all_scope_boundaries():
    cert = sm.certificate()
    assert cert['schema'] == 'oph.smooth_sign_field.v1'
    assert cert['total_complete_local_record_error_upper'] == '3/1048576'
    assert cert['fixed_delta_complete_local_thermodynamic_limit']
    assert not cert['same_as_original_Legendre_or_discontinuous_action']
    assert not cert['physical_GeV_mass_chiral_or_Lorentz_limit_claimed']
    assert not cert['roughness_used_as_premise_of_record_comparison']
    assert cert['smooth_action']['full_Fock_kinetic_density_rational_upper'] == '180'
    assert F(cert['smooth_action_recomputed_roughness_bound']) <= F(1, 2**20)
    root = Path(__file__).resolve().parents[2]
    assert all((root/path).is_file() for path in cert['source_inputs'])
    assert len(json.dumps(cert)) < 50000
