"""Reuse the complete reciprocal source and the regional tensor error theorem."""
from fractions import Fraction as F

import sympy as s

from model.finite_area_response import response as area
from model.matched_einstein_limit import frame


def metric_readouts(r):
    """Read the metric from the unchanged scalar action density."""
    r = s.sympify(r)
    if r.is_positive is not True:
        raise ValueError('positive source geometry required')
    density = s.diag(-1/r, 1/r, r*r, r*r)
    metric = (s.sqrt(-density.det())*density.inv()).applyfunc(s.simplify)
    return {
        'metric': metric,
        'area': s.simplify(2*s.sqrt(metric[2:, 2:].det())),
        'proper_volume': s.simplify(s.sqrt(metric[1:, 1:].det())),
        'lapse': s.simplify(s.sqrt(-metric[0, 0])),
        'radial_null_coordinate_speed_squared': s.simplify(-metric[0, 0]/metric[1, 1]),
        'transverse_null_coordinate_speed_squared': s.simplify(-metric[0, 0]/metric[2, 2]),
    }


def source_response():
    """Certified complete-law response of every two-valued metric readout."""
    matrices = area.matrices()
    spectra_weights = tuple(area.spectral_weights(m) for m in matrices)
    original = area.fixed_law_response(tuple(w for _, w in spectra_weights))
    branches = tuple(metric_readouts(s.Rational(r)) for r in area.GEOMETRIES)
    readings = {}
    for name in ('area', 'proper_volume', 'lapse',
                 'radial_null_coordinate_speed_squared',
                 'transverse_null_coordinate_speed_squared'):
        values = tuple(area.rational(b[name]) for b in branches)
        difference = values[1]-values[0]
        enclosure = original['dp1']*difference
        readings[name] = {'branch_values': list(map(str, values)),
                          'exact_branch_difference': str(difference),
                          'response_per_unit_randomization': enclosure.encoded()}
    lapse = original['dp1']*F(1, 5)
    if not lapse.lo > 0:
        raise ArithmeticError('native lapse response sign was not certified')
    ratio = F(-99, 200)
    combined = original['area']*ratio
    if not (combined.lo <= lapse.hi and lapse.lo <= combined.hi):
        raise ArithmeticError('area/lapse relation failed')
    if F(readings['area']['exact_branch_difference'])*ratio != F(1, 5):
        raise ArithmeticError('exact source branch algebra failed')
    return {
        'source': 'unchanged reciprocal scalar, r=(9/10,11/10), beta=1/sqrt(14)',
        'all_16_joint_energy_states_retained': True,
        'all_64_cross_geometry_edges_and_rejections_retained': True,
        'geometry_distribution_unchanged_at_intervention': True,
        'intervention': 'matter replacement mixture, then one original thermal update',
        'geometry_one_probability_slope': original['dp1'].encoded(),
        'readouts': readings,
        'exact_lapse_to_area_response_ratio': str(ratio),
        'fixed_branch_proper_time_ratio_for_equal_coordinate_duration': '9/11',
        'same_branch_proper_generator': 'H_g/r_g',
        'proper_time_rescaling_changes_complete_quantum_probabilities': False,
        'spatially_separated_clocks_or_light_bending_computed': False,
        'physical_gravitational_redshift_identified': False,
        'reservoir_preparation_and_pulse_controller_are_supplied': True,
    }


def newton_error_bound(eta, bulk_over_l4, area_over_l4, stationarity_over_l4,
                       support_radius, observation_radius):
    """Conditional weak-field transfer; outputs are pi times absolute errors.

    The input is the matched regional Einstein component-error bound, with
    independent physical source identification, a static spherical weak-field
    metric, zero pressure, and compact residual support supplied. The scalar
    source geometry above does not instantiate these hypotheses.
    """
    radius, r = F(support_radius), F(observation_radius)
    if radius <= 0 or r < radius:
        raise ValueError('positive support radius and exterior observation required')
    eps_pi = frame.normalized_error(eta, bulk_over_l4, area_over_l4, stationarity_over_l4)
    return {'pi_times_Einstein_component_error': str(eps_pi),
            'pi_times_Psi_Poisson_source_error': str(eps_pi/2),
            'pi_times_Phi_Poisson_source_error': str(2*eps_pi),
            'pi_times_Phi_force_error': str(2*eps_pi*radius**3/(3*r*r)),
            'pi_times_Phi_potential_error': str(2*eps_pi*radius**3/(3*r)),
            'pi_times_Psi_force_error': str(eps_pi*radius**3/(6*r*r)),
            'pi_times_Psi_potential_error': str(eps_pi*radius**3/(6*r)),
            'support_radius': str(radius), 'observation_radius': str(r),
            'Newton_strength_fitted_or_source_derived_here': False,
            'matched_regional_source_hypotheses_instantiated': False}


def linearized_curvature():
    """Calculate the static linearized Einstein tensor from metric derivatives."""
    x = s.symbols('x0:4', real=True)
    phi = s.Function('Phi')(*x[1:])
    psi = s.Function('Psi')(*x[1:])
    eta = s.diag(-1, 1, 1, 1)
    h = s.diag(-2*phi, -2*psi, -2*psi, -2*psi)
    trace = sum(eta[i, i]*h[i, i] for i in range(4))
    ricci = s.zeros(4)
    for mu in range(4):
        for nu in range(4):
            ricci[mu, nu] = s.expand((
                sum(eta[k, k]*(s.diff(h[k, nu], x[k], x[mu])
                              + s.diff(h[k, mu], x[k], x[nu])
                              - s.diff(h[mu, nu], x[k], x[k])) for k in range(4))
                - s.diff(trace, x[mu], x[nu]))/2)
    scalar = sum(eta[k, k]*ricci[k, k] for k in range(4))
    einstein = (ricci-eta*scalar/2).applyfunc(s.expand)
    lap = lambda f: sum(s.diff(f, z, z) for z in x[1:])
    checks = {
        'G00_equals_2_laplacian_Psi': s.simplify(einstein[0, 0]-2*lap(psi)) == 0,
        'spatial_trace_equals_2_laplacian_Phi_minus_Psi':
            s.simplify(sum(einstein[i, i] for i in range(1, 4))-2*lap(phi-psi)) == 0,
        'static_energy_flux_zero': all(einstein[0, i] == 0 for i in range(1, 4)),
    }
    if not all(checks.values()):
        raise ArithmeticError('weak-field curvature identities failed')
    return {'checks': checks, 'metric_convention': 'g00=-(1+2Phi), gij=(1-2Psi)deltaij',
            'equation_used': 'G=8*pi*G_N*T+E, Lambda negligible or subtracted background',
            'matter': 'static pressureless source',
            'Phi_error_source': '(E00+E11+E22+E33)/2',
            'Psi_error_source': 'E00/2',
            'no_source_Einstein_equation_is_assumed_proved': True}
