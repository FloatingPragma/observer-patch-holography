"""Exact bounds for the declared joint volume/coupling trajectory."""
from fractions import Fraction as F

from model.strong_coupling_refinement.bounds import finite_graph, rational
from model.full_fock_field.bounds import local_boundary_bound


def history(c, interaction, edges, duration, segments):
    if type(segments) is not int or segments < 1:
        raise ValueError('a positive integer number of free intervals is required')
    b = finite_graph(c, interaction, edges, duration)
    r = b['tau']/F(interaction)
    return min(F(2), 8*segments*r+2*F(duration)*b['generator_error_upper'])


def trajectory(edges, exchange=F(1), duration=F(1), repair_rate_per_edge=F(1)):
    if type(edges) is not int or edges < 2:
        raise ValueError('at least two edges are required for the stated trajectory')
    J = rational(exchange, positive=True)
    T, rate = rational(duration), rational(repair_rate_per_edge)
    c, U = J*edges**4, J*edges**8
    b = finite_graph(c, U, edges, T)
    return {'edges': edges, 'exchange': J, 'hopping': c, 'interaction': U,
            'isometry_error': b['isometry_error_upper'],
            'generator_error': b['generator_error_upper'],
            'charged_gap_lower': b['charge_gap_lower'],
            'poisson_mean_repairs': rate*edges*T,
            'complete_record_trace_norm_error': min(F(2), F(8, edges**3)*(1+rate*edges*T)+6*J*T/edges),
            'dressed_ground_initial_error': F(4, edges**3)}


def cube(radius, exchange=F(1), duration=F(1), repair_rate_per_edge=F(1)):
    if type(radius) is not int or radius < 1:
        raise ValueError('positive integer cube radius required')
    side = 2*radius+1
    result = trajectory(3*(side-1)*side**2, exchange, duration, repair_rate_per_edge)
    result.update({'radius': radius, 'vertices': side**3, 'spacing': 'fixed'})
    return result


def local_comparison(radius, support_radius, support_sites, exchange, duration, rate):
    b = cube(radius, exchange, duration, rate)
    # The bond supports/norms are exactly those used in the spin commutator-chain proof.
    tail = local_boundary_bound(radius, support_radius, support_sites, duration, exchange)
    return {'global_repaired_history_error': b['complete_record_trace_norm_error'],
            'spin_locality_error': tail,
            'local_expectation_error': min(F(2), b['complete_record_trace_norm_error']+tail)}


def torus_gap(width, exchange=F(1)):
    if type(width) is not int or width < 3 or width % 2 != 1:
        raise ValueError('an odd transverse width of at least three is required')
    J = rational(exchange, positive=True)
    length, area = 2*width**4, width**2
    vertices, edges = length*area, 3*length*area
    result = trajectory(edges, J, 0, 0)
    # pi < 22/7 and (22/7)^2 < 10. This intentionally retains a rational upper bound.
    result.update({'width': width, 'length': length, 'cross_section': area, 'vertices': vertices,
                   'spin_gap_upper': 10*J*area/length,
                   'source_gap_upper': 10*J*area/length+6*J/edges,
                   'gap_counts_ground_multiplicity': True,
                   'infinite_GNS_massless_particle_inferred': False})
    return result


def encode(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


def certificate():
    return encode({'complete_record_history_bound': 'min(2,8*m*tau/U+6*T*tau^3/U^2)',
                   'trajectory': 'c=J*E^4, U=J*E^8, J>0, E>=2',
                   'poisson_history_bound': 'min(2,8*(1+Gamma*E*T)/E^3+6*J*T/E)',
                   'cube_examples': [cube(r) for r in (1, 2, 4, 8)],
                   'local_example': local_comparison(160, 0, 1, F(1), F(1), F(1)),
                   'torus_gap_examples': [torus_gap(w) for w in (3, 5, 9, 17)]})
