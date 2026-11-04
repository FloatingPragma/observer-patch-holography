"""Exact finite bounds for bipartite source graphs."""
from fractions import Fraction as F

from ..strong_coupling_refinement.bounds import rational


def finite_graph(c, U, edges, time=0):
    c, U = rational(c, positive=True), rational(U, positive=True)
    time = rational(time)
    if type(edges) is not int or edges < 1:
        raise ValueError('positive integer edge count required')
    tau = c * edges
    if U < 4 * tau:
        raise ValueError('canonical dynamics requires U >= 4 tau')
    return {
        'tau': tau,
        'exchange': c*c/U,
        'charge_gap_lower': U/2-tau,
        'ordered_energy_error_upper': 2*tau**4/(U**2*(U-tau)),
        'canonical_generator_error_upper': 7*tau**4/U**3,
        'complete_dynamics_error_upper': min(F(2), 4*tau/U+7*time*tau**4/U**3),
        'bare_embedding_error_upper': 2*tau/U,
    }


def refinement_exponents(power):
    p = rational(power)
    return {
        'angle_power': 4-p,
        'dynamics_remainder_power': 16-3*p,
        'exchange_power': 2-p,
        'global_error_proved_to_vanish': p > F(16, 3),
        'fixed_nonzero_exchange': p == 2,
    }


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def certificate():
    graph_coefficient = F(2)+F(3, 2)+F(1, 2)
    canonical_coefficient = F(9, 8)*graph_coefficient+F(17, 8)
    if graph_coefficient != 4 or canonical_coefficient != F(53, 8) or canonical_coefficient > 7:
        raise ArithmeticError('canonical generator majorant failed')
    return encode({
        'graph_coefficient': graph_coefficient,
        'canonical_coefficient_before_rounding': canonical_coefficient,
        'dimer': finite_graph(128, 128**2, 1, 2),
        'cube': finite_graph(65536, 65536**2, 12, 1),
        'physical_box': {str(p): refinement_exponents(p) for p in (2, 4, F(16, 3), F(11, 2), 6)},
        'interacting_relativistic_continuum_proved': False,
        'physical_species_identified': False,
    })
