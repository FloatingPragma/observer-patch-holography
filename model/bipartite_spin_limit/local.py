"""Volume-independent local bounds for dressed fixed-spacing spin dynamics.

The proof is LOCAL_LIMIT.md. These majorants use the integer onsite defect
penalty, not a gap of the extensive interacting Hamiltonian.
"""
from fractions import Fraction as F

from .bounds import encode
from ..strong_coupling_refinement.bounds import rational


STEPS = 5


def local_bounds(c, U, time=0, support=1):
    c, U = rational(c, positive=True), rational(U, positive=True)
    time = rational(time)
    if type(support) is not int or support < 1:
        raise ValueError('positive integer observable support required')
    B, nu = 144*c, U/2
    if nu < 2**20*B:
        raise ValueError('local construction requires U/2 >= 2^20 (144 c)')
    rho = 4096*B/nu
    remainder = B*rho**STEPS
    spin_error = 2**25*B**3/nu**2
    exchange = c*c/U
    dressed_error = 2**16*support*time*(
        remainder*(1+128*B*time)**3
        + spin_error*(1+3072*exchange*time)**3)
    dressing = 8*support*B/nu
    return {
        'local_interaction_majorant': B,
        'onsite_penalty': nu,
        'iteration_contraction': rho,
        'off_diagonal_local_error': remainder,
        'spin_generator_local_error': spin_error,
        'exchange': exchange,
        'dressed_observable_error': min(F(2), dressed_error),
        'bare_observable_on_dressed_input_error': min(F(2), dressed_error+dressing),
        'observable_dressing_error': min(F(2), dressing),
        'depends_on_total_volume': False,
        'undressed_input_comparison': False,
        'interleaved_repair_comparison': False,
    }


def iteration_majorants():
    """Worst admitted scale; every operation is rational and volume free."""
    B, nu, delta = F(1), F(2**20), F(1, 12)
    diagonal, off = B/3, 2*B/3
    rho = 4096*B/nu
    rows = []
    for step in range(STEPS):
        q = 24*off/(nu*delta)
        change = q/(1-q)*(diagonal+off)
        updated_diagonal, updated_off = diagonal+change, 2*change
        if (q >= F(1, 2) or updated_off > rho*off
                or updated_diagonal+updated_off > 2*B):
            raise ArithmeticError('local iteration majorant failed')
        rows.append({'step': step+1, 'bch_ratio': q,
                     'diagonal': updated_diagonal, 'off_diagonal': updated_off})
        diagonal, off = updated_diagonal, updated_off
    first_cubic = 2*(288*B/nu)**2*B
    later_diagonal = B*rho**2
    if first_cubic+later_diagonal > 2**25*B**3/nu**2:
        raise ArithmeticError('effective-spin remainder majorant failed')
    return encode(rows)


def certificate():
    return {
        'steps': STEPS,
        'interaction_decay_initial': '1',
        'interaction_decay_final_lower': '1/2',
        'maximum_source_degree': 6,
        'integer_penalty': 'sum_x (n_x-1)^2',
        'majorant_controls': iteration_majorants(),
        'fixed_exchange_example': encode(local_bounds(2**256, 2**512, 1, 8)),
        'fixed_exchange_asymptotic_powers_in_c': {
            'off_diagonal_local_error': -4,
            'spin_generator_local_error': -1,
            'full_local_dynamics_error_at_fixed_time': -1,
            'observable_dressing_error': -1,
        },
        'thermodynamic_local_dynamics_constructed': True,
        'source_ground_state_identified': False,
        'continuum_relativistic_field_constructed': False,
    }
