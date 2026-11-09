"""Weak-field Hessian of the actual trace-metric M15 coordinate action.

The calculation uses all twelve native gauge directions and the inherited
electric Casimir. It is a tangent quadratic comparison on the same clock;
it supplies no interacting vacuum or nonperturbative particle pole.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

from model.common_source.chiral_candidate import action
from model.particle_field_preparation import preparation
from model.particle_geometric_response import geometry
from model.particle_local_wilson_density import density
from model.particle_observer_outputs import field_route

ROOT = Path(__file__).resolve().parents[2]
SOURCE_INPUTS = (
    'model/particle_field_preparation/preparation.py',
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/particle_local_wilson_density/density.py',
    'model/particle_smooth_sign_field/smooth.py',
    'model/particle_geometric_response/geometry.py',
    'model/common_source/chiral_candidate/action.py',
    'model/observer_particles/results/representation_identification_operators.json',
    'model/observer_particles/results/representation_identification.json',
    'model/observer_particles/coefficient_selection.py',
    'Lean/Screen/KineticFamilyCancellation.lean',
    'Lean/Screen/GaugeKineticInvariantForms.lean',
)


def clean(matrix):
    return s.ImmutableMatrix(matrix).applyfunc(s.simplify)


def positive(value, name):
    if type(value) not in (int, str, F):
        raise ValueError(name+' must be an exact positive rational')
    result = F(value)
    if result <= 0:
        raise ValueError(name+' must be positive')
    return result


def unit(i, j):
    out = s.zeros(5)
    out[i, j] = 1
    return out


@lru_cache(maxsize=1)
def native_generators():
    """All12 Hermitian generators in the retained ordered M15 exterior basis.

    SU3 and SU2 fundamental trace is 1/2. Y=q6/6 and Q=T3+Y.
    The inherited twelve port currents are reconstructed in this basis.
    """
    native = field_route.source_data.load()
    packet = action._pinned('representation_identification_operators.json',
                            'representation_identification.json')
    subsets = tuple(tuple(packet['basis_subsets'][i]) for i in packet['even_indices'])
    color = (
        (unit(0, 1)+unit(1, 0))/2,
        (-s.I*unit(0, 1)+s.I*unit(1, 0))/2,
        (unit(0, 0)-unit(1, 1))/2,
        (unit(0, 2)+unit(2, 0))/2,
        (-s.I*unit(0, 2)+s.I*unit(2, 0))/2,
        (unit(1, 2)+unit(2, 1))/2,
        (-s.I*unit(1, 2)+s.I*unit(2, 1))/2,
        (unit(0, 0)+unit(1, 1)-2*unit(2, 2))/(2*s.sqrt(3)),
    )
    weak = ((unit(3, 4)+unit(4, 3))/2,
            (-s.I*unit(3, 4)+s.I*unit(4, 3))/2,
            (unit(3, 3)-unit(4, 4))/2)
    hypercharge = s.diag(-s.Rational(1, 3), -s.Rational(1, 3),
                         -s.Rational(1, 3), s.Rational(1, 2), s.Rational(1, 2))
    fundamental = tuple(map(clean, (*color, *weak, hypercharge)))
    generators = tuple(clean(action._exterior_derivative(t, subsets)) for t in fundamental)
    gram = clean(s.Matrix(12, 12, lambda i, j: s.trace(generators[i]*generators[j])))
    expected = s.diag(*([2]*11+[s.Rational(10, 3)]))
    if gram != expected or any(t != t.H for t in generators):
        raise ArithmeticError('native M15 trace form differs')
    if tuple(6*x for x in generators[-1].diagonal()) != native['charges']:
        raise ArithmeticError('the same native integer charge convention is required')
    ports = tuple(action._exterior_derivative(action._load_sparse(g), subsets)
                  for g in packet['spin_lifted_one_body_generators'])
    coefficients = []
    for current in ports:
        h = clean(s.I*current)
        row = tuple(s.simplify(s.trace(t*h)/gram[i, i]) for i, t in enumerate(generators))
        if clean(sum((c*t for c, t in zip(row, generators, strict=True)), s.zeros(15))-h) != s.zeros(15):
            raise ArithmeticError('a native port current lies outside the gauge basis')
        coefficients.append(row)
    if s.Matrix(coefficients).rank() != 12:
        raise ArithmeticError('the complete twelve-dimensional source is required')
    casimir = clean(sum((t*t/gram[i, i] for i, t in enumerate(generators)), s.zeros(15)))
    casimirs = {r['q6Y']: r['trace_metric_C2'] for r in preparation.native_casimirs()}
    if casimir != s.diag(*(casimirs[q] for q in native['charges'])):
        raise ArithmeticError('the inherited electric Casimir normalization differs')
    return {'labels': tuple([f'color_{i}' for i in range(1, 9)]
                            +[f'weak_{i}' for i in range(1, 4)]+['Y']),
            'fundamental': fundamental, 'generators': generators, 'gram': gram,
            'native_port_coefficients': tuple(coefficients), 'casimir': casimir,
            'native_source_provenance': native['provenance']}


def plaquette_hessian(signed_links):
    """Expand the ordered product of four exp(i epsilon X_j) through order2.

    Link3 and link4 include the negative orientation in X_j. No commuting
    assumption is made. The result is the coefficient of epsilon squared,
    so the second derivative is twice the returned coefficient.
    """
    xs = tuple(clean(x) for x in signed_links)
    if len(xs) != 4 or any(x.shape != (15, 15) or x != x.H for x in xs):
        raise ValueError('four Hermitian native M15 signed link matrices required')
    first = sum(xs, s.zeros(15))
    second = -sum((x*x/2 for x in xs), s.zeros(15))
    second -= sum((xs[i]*xs[j] for i in range(4) for j in range(i+1, 4)), s.zeros(15))
    actual = s.simplify(-s.re(s.trace(second))/15)
    expected = s.simplify(s.trace(first*first)/30)
    if actual != expected:
        raise ArithmeticError('ordered noncommuting plaquette Hessian differs')
    return {'linear_curvature': clean(first), 'quadratic_cost_coefficient': actual,
            'cost_second_derivative': 2*actual}


def weak_field_remainder(integrated_link_norm):
    """Bound normalized trace remainder after order2 for a four-link product.

    Input x=abs(epsilon)*sum_j ||X_j||, with exact rational upper bounds.
    This applies to the magnetic multiplier, not to a prepared-state error.
    """
    if type(integrated_link_norm) not in (int, str, F):
        raise ValueError('integrated link norm requires exact rational data')
    x = F(integrated_link_norm)
    if x < 0:
        raise ValueError('integrated link norm must be nonnegative')
    return field_route.exponential_tail_upper(x, 3)


def constant_nonabelian_loop(x, y):
    """Quartic interaction of a constant noncommuting source plaquette.

    This retains the first nonzero magnetic cost of a zero-momentum
    tangent direction. It is not a physical scattering amplitude.
    """
    x, y = clean(x), clean(y)
    if any(t.shape != (15, 15) or t != t.H for t in (x, y)):
        raise ValueError('two Hermitian native M15 matrices required')
    from math import factorial
    coefficients = [s.eye(15)]+[s.zeros(15) for _ in range(4)]
    for link in (x, y, -x, -y):
        expansion = [(s.I*link)**n/factorial(n) for n in range(5)]
        coefficients = [sum((coefficients[j]*expansion[n-j] for j in range(n+1)),
                            s.zeros(15)) for n in range(5)]
    costs = [s.simplify(-s.re(s.trace(v))/15) for v in coefficients[1:]]
    bracket = x*y-y*x
    quartic = s.simplify(-s.trace(bracket*bracket)/30)
    if costs[:3] != [0, 0, 0] or costs[3] != quartic:
        raise ArithmeticError('the native nonabelian quartic interaction differs')
    return {'linear_through_cubic_cost_coefficients': costs[:3],
            'quartic_cost_coefficient': quartic,
            'nonzero_nonlinear_magnetic_interaction': quartic != 0,
            'interacting_particle_scattering_derived': False}


def curl_symbol(difference):
    """Forward lattice curl, retaining its longitudinal kernel."""
    d = clean(s.Matrix(difference))
    if d.shape != (3, 1):
        raise ValueError('three forward-difference components required')
    return clean(s.Matrix([[0, -d[2], d[1]], [d[2], 0, -d[0]], [-d[1], d[0], 0]]))


def lattice_symbol(phases):
    """Exact symbol on unit-circle phases z_i, including the zero mode."""
    z = tuple(map(s.sympify, phases))
    if len(z) != 3 or any(s.simplify(v*s.conjugate(v)-1) != 0 for v in z):
        raise ValueError('three exact unit-circle phases required')
    d = clean(s.Matrix([v-1 for v in z]))
    curl = curl_symbol(d)
    norm2 = s.simplify((d.H*d)[0])
    hessian = clean(curl.H*curl)
    if hessian != clean(norm2*s.eye(3)-d*d.H) or clean(hessian*d) != s.zeros(3, 1):
        raise ArithmeticError('the forward lattice curl or gauge kernel differs')
    rank = 0 if norm2 == 0 else 2
    if hessian.rank() != rank or clean(hessian*hessian-norm2*hessian) != s.zeros(3):
        raise ArithmeticError('the transverse lattice eigenvalues differ')
    return {'difference': d, 'curl': curl, 'hessian': hessian,
            'transverse_eigenvalue': norm2, 'rank': rank,
            'longitudinal_gauge_null_retained': True,
            'zero_momentum_flat_coordinates': 3 if norm2 == 0 else 0}


def clock_family(magnetic_scale='1/7', electric_scale='1'):
    """HE=e0*C_trace/(2a), HB=−kappa0*ReTr(U_p)/(15a).

    The inherited electric form is inverse to Tr_M15(t_a*t_b). Equal
    coefficients on canonical factor generators specify a different form.
    """
    k = positive(magnetic_scale, 'magnetic scale')
    e = positive(electric_scale, 'electric scale')
    data = native_generators()
    indices = [F(str(data['gram'][i, i])) for i in range(12)]
    velocities = [e/index*k*index/15 for index in indices]
    if len(set(velocities)) != 1:
        raise ArithmeticError('the trace-dual source lost its common gauge cone')
    return {'magnetic_scale': k, 'electric_scale': e,
            'electric_coefficients_canonical_basis': [e/i for i in indices],
            'magnetic_quadratic_coefficients_canonical_basis': [k*i/15 for i in indices],
            'all_twelve_speed_squared': velocities,
            'Dirac_principal_speed_squared': F(1),
            'same_clock_gauge_Dirac_characteristic_match': e*k == 15,
            'matching_electric_scale_at_same_magnetic_scale': 15/k,
            'is_an_architecture_family_without_mass_calibration': True,
            'physical_interacting_limiting_speed_proved': False}


def charge_normalization(scale='6'):
    """t'=s*t, theta'=theta/s, pi'=s*pi preserves energy and cone."""
    q = positive(scale, 'charge scale')
    index = F(10, 3)
    return {'generator_scale': q, 'trace_index_Y': index,
            'trace_index_rescaled': q*q*index,
            'electric_coefficient_Y_at_unit_electric_scale': 1/index,
            'electric_coefficient_rescaled': 1/(q*q*index),
            'magnetic_coefficient_rescaled_over_original': q*q,
            'unchanged_electric_times_magnetic_index': F(1),
            'canonical_pair': 'theta_prime=theta/s, pi_prime=s*pi',
            'integer_q6_equals_6Y': q == 6}


def free_symbol(momentum, spacing='1'):
    """The inherited ordinary Wilson action on the exact native Dirac rows."""
    a = s.Rational(str(positive(spacing, 'spacing')))
    k = tuple(map(s.sympify, momentum))
    if len(k) != 3 or any(x.is_real is not True for x in k):
        raise ValueError('three real lattice momenta required')
    gamma = density.gamma_data()
    beta = gamma[3]
    alpha = tuple(clean(s.I*beta*g) for g in gamma[:3])
    wilson = sum((1-s.cos(x) for x in k), s.S.Zero)
    h = clean((sum((alpha[i]*s.sin(k[i]) for i in range(3)), s.zeros(4))+beta*wilson)/a)
    energy2 = s.simplify((sum((s.sin(x)**2 for x in k), s.S.Zero)+wilson**2)/a**2)
    if clean(h*h-energy2*s.eye(4)) != s.zeros(4):
        raise ArithmeticError('native Clifford free dispersion differs')
    return {'Hamiltonian': h, 'energy_squared': energy2,
            'all_Dirac_rows': 4, 'both_chiralities_retained': True}


def dispersion_error(spacing='1/16', momentum_norm_upper='1'):
    """Uniform low-momentum error envelopes; momentum is in the source units.

    0 <= |p|²−lambda_lattice/a² <= a²|p|⁴/12.
    |E_W²−|p|²| <= 7a²|p|⁴/12.
    The Wilson operator error is <= a|p|²/2+a²|p|³/6.
    """
    a = positive(spacing, 'spacing')
    p = positive(momentum_norm_upper, 'momentum norm upper')
    return {'spacing': a, 'momentum_norm_upper': p,
            'gauge_unit_speed_squared_error_upper': a*a*p**4/12,
            'Wilson_energy_squared_error_upper': 7*a*a*p**4/12,
            'Wilson_operator_error_upper': a*p*p/2+a*a*p**3/6,
            'comparison_is_uniform_on_fixed_bounded_momentum_sets': True,
            'many_body_interacting_pole_error_included': False}


def coframe_principal(coframe, covector):
    """Both tangent fields use the inherited coframe, without a time rechart.

    This is a continuum principal symbol. The exact finite lattice curl
    statement above is for the flat lattice and retains its complex phases.
    """
    data = geometry.coframe_data(coframe)
    k = clean(s.Matrix(covector))
    if k.shape != (3, 1) or any(x.is_real is not True for x in k):
        raise ValueError('a real three-component covector is required')
    metric = data['electric_magnetic_coefficient']
    curl = curl_symbol(k)
    gauge = clean(-metric*curl*metric*curl)
    q = s.simplify((k.T*data['inverse_metric']*k)[0])
    identity = clean(q*s.eye(3)-k*k.T*data['inverse_metric'])
    if gauge != identity or clean(gauge* k) != s.zeros(3, 1):
        raise ArithmeticError('constant-coframe gauge principal identity differs')
    gamma = density.gamma_data()
    physical_k = data['inverse_frame']*k
    dirac = clean(sum((physical_k[i]*s.I*gamma[3]*gamma[i] for i in range(3)), s.zeros(4)))
    if clean(dirac*dirac-q*s.eye(4)) != s.zeros(4):
        raise ArithmeticError('the same coframe Dirac cone differs')
    return {'coframe': clean(coframe), 'inverse_metric': data['inverse_metric'],
            'covector': k, 'gauge_principal_matrix_at_matched_clock': gauge,
            'characteristic_squared': q, 'Dirac_principal_matrix': dirac,
            'gauge_null_retained': True,
            'two_transverse_responses_for_nonzero_covector': q != 0,
            'interacting_Lorentz_symmetry_or_Einstein_equation_inferred': False}


def serialize(value):
    if isinstance(value, s.MatrixBase):
        return [[str(v) for v in value.row(i)] for i in range(value.rows)]
    if isinstance(value, dict):
        return {str(k): serialize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(v) for v in value]
    if isinstance(value, (F, s.Basic)):
        return str(value)
    return value


def certificate():
    d = native_generators()
    gs = d['generators']
    xs = (gs[0]+gs[8]/3, gs[1]/2+gs[11], -gs[2]/5, -gs[9]/7)
    plaquette = plaquette_hessian(xs)
    flat = lattice_symbol((s.I, -1, 1))
    shear = s.Matrix([[1, s.Rational(1, 3), 0], [0, 1, s.Rational(1, 5)], [0, 0, 1]])
    data = {
        'schema': 'oph.native_gauge_characteristics.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_INPUTS},
        'native_source_provenance': d['native_source_provenance'],
        'native_basis': {'labels': d['labels'], 'trace_Gram': d['gram'],
                        'native_port_coefficients': d['native_port_coefficients'],
                        'all_twelve_generators_used': True,
                        'charge_convention': 'Y=q6/6, Q=T3+Y, all-left internal M15',
                        'trace_metric_Casimir_diagonal': list(d['casimir'].diagonal())},
        'plaquette_control': {'ordered_noncommuting_quadratic_cost': plaquette['quadratic_cost_coefficient'],
                             'remainder_for_integrated_link_norm_1_over_8': weak_field_remainder('1/8')},
        'constant_nonabelian_plaquette': constant_nonabelian_loop(gs[8], gs[9]),
        'lattice_curl_control': flat,
        'zero_momentum_control': lattice_symbol((1, 1, 1)),
        'inherited_unit_electric_example': clock_family(),
        'declared_common_clock_example': clock_family(electric_scale='105'),
        'primitive_charge_covariance': charge_normalization(),
        'free_Wilson_example': free_symbol((s.pi/2, 0, 0)),
        'low_momentum_error': dispersion_error(),
        'constant_coframe_control': coframe_principal(shear, (1, 2, -1)),
        'scope': {'same_native_magnetic_Hessian_and_electric_metric': True,
                  'actual_configured_quantum_field_state_changed': False,
                  'common_clock_normalization_selected_from_measured_data': False,
                  'quadratic_tangent_comparison_requires_weak_field_domain': True,
                  'weak_field_localization_of_ground_state_proved': False,
                  'nonperturbative_interacting_relativistic_limit_proved': False,
                  'named_particle_mass_or_running_coupling_prediction': False},
    }
    return serialize(data)
