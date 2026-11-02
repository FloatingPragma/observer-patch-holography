"""Exact rational bounds; no floating fit determines a theorem threshold."""
from fractions import Fraction as F


def rational(value, *, positive=False):
    if type(value) not in (int, F):
        raise ValueError('an exact rational is required')
    out = F(value)
    if out < 0 or (positive and not out):
        raise ValueError('a positive rational is required' if positive else 'negative input')
    return out


def finite_graph(c, U, edges, time=0):
    c, U, time = rational(c, positive=True), rational(U, positive=True), rational(time)
    if type(edges) is not int or edges < 1:
        raise ValueError('a positive integer edge count is required')
    tau = c * edges
    if U < 4*tau:
        raise ValueError('this canonical-isometry bound requires U >= 4 c |E|')
    return {'tau': tau, 'charge_gap_lower': U/2-tau,
            'spin_exchange': c*c/U, 'isometry_error_upper': 2*tau/U,
            'generator_error_upper': 3*tau**3/U**2,
            'ordered_energy_error_upper': tau**3/(U-tau)**2,
            'complete_dynamics_error_upper': min(F(2), 4*tau/U+3*time*tau**3/U**2)}


def refinement_exponents(power):
    """Open physical boxes: tau=O(L^4), U=L^power, c=O(L)."""
    p = rational(power)
    return {'angle_power': 4-p, 'dynamics_remainder_power': 12-2*p,
            'exchange_power': 2-p, 'physical_spin_speed_power': 1-p,
            'global_error_proved_to_vanish': p > 6,
            'nonzero_fixed_exchange': p == 2}


def growing_schedule(k):
    """Finite exact envelope, without constructing 2**(256*k*k).

    c=2**(256 k^2), U=k^2, delta=1/k, R=k^3.
    Gaussian + exponential majorants give U^2 M(R) <= 3073/k^2.
    """
    if type(k) is not int or not 1 <= k <= 10**9:
        raise ValueError('integer 1 <= k <= 10^9 required')
    # k <= 2^k: prove each exponentially small remainder <= 1/(2 k^2).
    first_exponent = 3+9*k-256*k*k
    second_exponent = 11+8*k-195*k*k
    if max(first_exponent, second_exponent) > -1:
        raise ArithmeticError('exponential remainder proof failed')
    # k^8/(2 c^2) <= 1/2 < 2/3, so exp(prefactor) <= 2.
    if 8*k-512*k*k > 0:
        raise ArithmeticError('outer exponential bound failed')
    return {'k': k, 'U': k*k, 'log2_c': 256*k*k, 'R': k**3,
            'delta': F(1, k), 'U_squared_mass_upper': F(3073, k*k),
            'U_sqrt_mass_upper': F(56, k), 'U_over_R': F(1, k),
            'unit_packet_unit_time_vacuum_error_upper': min(F(2), F(30, k)),
            'exponential_remainder_powers': [first_exponent, second_exponent]}


def _json(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {k: _json(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_json(v) for v in value]
    return value


def certificate():
    return _json({
        'finite_graph_domain': 'U >= 4 tau, tau=c |E|, c>0; N=|V|',
        'canonical_isometry_error': '2 tau/U',
        'complete_dynamics_error': 'min(2,4 tau/U+3 |t| tau^3/U^2)',
        'ordered_low_energy_error': 'tau^3/(U-tau)^2',
        'charge_gap_lower': 'U/2-tau',
        'fixed_graph_limit': 'U=c^2/J, fixed J>0, c->infinity',
        'dimer_interaction_control': finite_graph(128, 128**2, 1, 2),
        'cube_fixed_graph_control': finite_graph(65536, 65536**2, 12, 1),
        'physical_box_power_examples': {str(p): refinement_exponents(p) for p in (2, 4, 6, 7)},
        'growing_contact_free_domain': 'U(c)=o((log c)^(3/2)), c=v/a -> infinity',
        'growing_free_schedules': [growing_schedule(k) for k in (64, 256, 1024)],
        'finite_spin_limit_is_interacting': True,
        'fixed_physical_box_interacting_continuum_proved': False,
        'global_norm_failure_is_a_model_no_go': False,
        'growing_contact_free_limit_requires_new_compactness_argument': True,
    })
