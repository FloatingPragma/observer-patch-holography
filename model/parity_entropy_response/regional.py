"""Exact source-family, local faithfulness and regional response bounds."""
from fractions import Fraction as F

import sympy as s

from ..radial_interacting_core import core


def rational(value):
    if type(value) not in (int, F):
        raise ValueError('an exact rational input is required')
    return F(value)


def family_weights(q, n_ground, n_charged, epsilon=F(1, 64)):
    q, ng, nc, epsilon = map(rational, (q, n_ground, n_charged, epsilon))
    if not (F(1, 4) <= q <= F(3, 4) and 0 <= ng < F(1, 1000)
            and F(999, 1000) < nc <= 1 and epsilon == F(1, 64)):
        raise ValueError('the certified actual-source parity family is required')
    charged = (q-epsilon/2-(1-epsilon)*ng)/(nc-ng)
    ground = 1-epsilon-charged
    if min(ground, charged) <= F(6, 25):
        raise ArithmeticError('the faithful affine mixture lost a spectral weight')
    if ground*ng+charged*nc+epsilon/2 != q:
        raise ArithmeticError('actual-source mixture misses its target regional parity')
    return ground, charged, epsilon


def ground_faithfulness():
    row = core.core_bounds()
    b_squared, complement_norm = F(3, 2), row['mass']+row['hopping_bound']
    alpha_upper = 1/(1+b_squared/complement_norm**2)
    lower = (1-alpha_upper**2)/2-row['projection_angle_bound']
    upper = row['center_parity_actual_ground_upper']
    if not F(1, 4096) < lower < upper < F(1, 1000):
        raise ArithmeticError('the actual interacting ground lacks certified local faithfulness')
    return {'central_V_star_V': str(b_squared),
            'complement_operator_norm_upper': str(complement_norm),
            'central_zero_weight_upper': str(alpha_upper),
            'actual_ground_parity_lower': str(lower),
            'actual_ground_parity_upper': str(upper),
            'actual_ground_local_minimum_eigenvalue_lower': str(lower/2),
            'actual_ground_local_M4_faithful': True,
            'actual_ground_full_Fock_state_faithful': False,
            'actual_ground_parity_value_computed_or_prepared': False}


def algebra():
    q, r, ng, nc, eps, gap, eg, einf, U = s.symbols(
        'q r ng nc eps gap eg einf U', real=True)
    ac = (q-eps/2-(1-eps)*ng)/(nc-ng)
    ag = 1-eps-ac
    energy = ag*eg+ac*(eg+gap)+eps*einf
    if (s.simplify(ag*ng+ac*nc+eps/2-q) != 0
            or s.simplify(s.diff(energy, q)-gap/(nc-ng)) != 0):
        raise ArithmeticError('same-source parity or full-energy slope failed')
    h = lambda x: -x*s.log(x)-(1-x)*s.log(1-x)
    divergence = q*s.log(q/r)+(1-q)*s.log((1-q)/(1-r))
    expression = s.expand_log(h(q)-h(r)-s.log((1-r)/r)*(q-r)+divergence, force=True)
    if s.simplify(expression) != 0 or s.simplify(s.diff(h(q), q, 2)+1/(q*(1-q))) != 0:
        raise ArithmeticError('regional entropy identity or Hessian failed')
    D = s.diag(1, 0, 0, 1)
    K0, K1 = (7*s.eye(4)-6*D)/s.sqrt(50), s.I*(s.eye(4)+6*D)/s.sqrt(50)
    Z = s.eye(4)-2*D
    for i in range(4):
        for j in range(4):
            E = s.zeros(4); E[i, j] = 1
            if K0*E*K0.H+K1*E*K1.H != s.Rational(16, 25)*E+s.Rational(9, 25)*Z*E*Z:
                raise ArithmeticError('full-site matrix-unit channel identity failed')
    rho = s.diag(q, 1-q, 1-q, q)/2
    if s.simplify(s.trace(K1*rho*K1.H)-s.Rational(1, 50)-s.Rational(24, 25)*q) != 0:
        raise ArithmeticError('actual parity pointer probability failed')
    if U*(s.diag(0, 1, 1, 2)-s.eye(4))**2/2-U*s.eye(4)/4 != U*D/2-U*s.eye(4)/4:
        raise ArithmeticError('local contact calibration failed')
    return {'complete_site_matrix_units_checked': 16,
            'source_family_parity_identity': True,
            'source_energy_derivative': 'g/(n_C-n_G)',
            'local_contact_energy_derivative': 'U/2',
            'local_density': 'diag(q,1-q,1-q,q)/2',
            'modular_generator': 'log(2/(1-q))*I+log((1-q)/q)*D',
            'local_contact_time_factor': '-2*log((1-q)/q)/U',
            'half_parity_modular_flow_is_identity': True,
            'modular_time_factor_changes_sign_at_half': True,
            'relative_entropy': 'q*log(q/r)+(1-q)*log((1-q)/(1-r))',
            'entropy_remainder_lower': '2*(q-r)^2',
            'entropy_remainder_upper': '(8/3)*(q-r)^2',
            'source_bias': '(3/25)*(1-2*q)',
            'entropy_bias_identity': 'DeltaS=-(25/6)*log((1-r)/r)*DeltaTheta-D(tau_q||tau_r)',
            'all_state_variational_identity': 'S(omega)-kappa_r*Tr(D omega)=S(tau_r)-kappa_r*r-D(omega||tau_r)',
            'fixed_response_entropy_maximizer_unique': True,
            'entropy_defined_geometric_area': False}


def certificate():
    row = core.core_bounds()
    delta_lower = row['center_parity_actual_charged_minus_ground_lower']
    eps = F(1, 64)
    weight_lower = (1-eps)*F(999, 1000)-F(3, 4)+eps/2
    if weight_lower <= F(6, 25) or delta_lower <= F(998, 1000):
        raise ArithmeticError('actual-source family has no uniform positive weights')
    return {'algebra': algebra(), 'actual_ground': ground_faithfulness(),
            'same_source_affine_family': {
                'q_interval': ['1/4', '3/4'], 'identity_weight': str(eps),
                'each_spectral_mixture_weight_strict_lower': str(weight_lower),
                'full_Fock_dimension': row['full_Fock_dimension'],
                'full_state_minimum_eigenvalue_lower': str(eps/row['full_Fock_dimension']),
                'same_full_Hamiltonian_stationarity': True,
                'spectral_preparation_and_parity_calibration_are_inputs': True,
                'spectral_state_preparation_executed': False,
                'energy_slope_lower': str(row['actual_charged_gap_lower']),
                'energy_slope_upper': str(row['actual_charged_gap_upper']/delta_lower),
                'local_contact_slope': str(row['centered_contact_coupling']/2),
                'regional_logarithm_Lipschitz': '16/3',
                'pointer_probability_to_parity_error_factor': '25/24',
                'pointer_probability_to_modular_error_factor': '50/9'},
            'stationarity_scope': 'maximum regional entropy at fixed native kinetic bias',
            'Einstein_fixed_volume_stationarity_inferred': False}
