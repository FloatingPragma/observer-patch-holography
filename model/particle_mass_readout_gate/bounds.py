"""Exact common-clock restriction on the specified protected spin family."""
from fractions import Fraction as F

from model.protected_spin_limit.bounds import trajectory


def charged_clock_bound(edges, energy_cap=F(1)):
    if type(edges) is not int or edges < 2:
        raise ValueError('integer edge count at least two required')
    if type(energy_cap) not in (int, F) or energy_cap <= 0:
        raise ValueError('positive exact charged excitation-energy cap required')
    native = trajectory(edges, F(1), F(0), F(0))
    coefficient = F(edges**8, 2)-edges**5
    if native['charged_gap_lower'] != coefficient or coefficient < F(edges**8, 4):
        raise ArithmeticError('owning protected-spin trajectory changed')
    return {
        'edge_count': edges,
        'charged_energy_lower_per_exchange': str(coefficient),
        'charged_energy_cap': str(F(energy_cap)),
        'rescaled_exchange_upper_if_a_charged_excitation_meets_cap':
            str(F(energy_cap)/coefficient),
        'simpler_rescaled_exchange_upper': str(4*F(energy_cap)/edges**8),
        'scope': 'fixed-spacing local spin dynamics on the same rescaled clock',
    }


def ratio(numerator, denominator):
    if type(numerator) not in (int, F) or type(denominator) not in (int, F):
        raise ValueError('exact rational inputs required')
    if denominator == 0:
        raise ValueError('a zero invariant parameter cannot normalize a mass ratio')
    return F(numerator)/F(denominator)


def charged_weight_bound(edges, exchange=F(1), mean_excitation_cap=F(1)):
    if type(edges) is not int or edges < 2:
        raise ValueError('integer edge count at least two required')
    if type(exchange) not in (int, F) or exchange <= 0:
        raise ValueError('positive exact exchange required')
    if type(mean_excitation_cap) not in (int, F) or mean_excitation_cap < 0:
        raise ValueError('nonnegative exact mean excitation-energy cap required')
    gap = trajectory(edges, F(exchange), F(0), F(0))['charged_gap_lower']
    return {
        'edge_count': edges, 'exchange': str(F(exchange)),
        'positive_ground_subtracted_clock': 'L=K-E_ground >= 0',
        'state': 'any normalized density operator with Tr(rho L)<=C',
        'charge_projection': 'Pi_(N!=sites)',
        'mean_excitation_cap': str(F(mean_excitation_cap)),
        'charged_weight_upper': str(min(F(1), F(mean_excitation_cap)/gap)),
        'simpler_charged_weight_upper': str(min(F(1), 4*F(mean_excitation_cap)/(F(exchange)*edges**8))),
        'relativistic_or_electromagnetic_identification': False,
    }
