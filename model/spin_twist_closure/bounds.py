"""Rational source bounds for single-axis twists and averaged bonds."""
from fractions import Fraction as F


def rational(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, F)):
        raise ValueError(name+' must be an exact rational')
    return F(value)


def twist_remainder(p, differences):
    p = rational(p, 'p')
    differences = tuple(rational(x, 'angle difference') for x in differences)
    if p <= 0:
        raise ValueError('positive scale required')
    return sum((abs(x)**3 for x in differences), F())/(12*p)


def compact_twist(p, radius, bound, axis):
    """A smooth test supported in [-R,R]^3, p<=1 and R>=1.

    For P, bound controls each coordinate derivative. For Q it controls
    the absolute test value. These are analytic input bounds, not samples.
    """
    p, radius, bound = [rational(x, n) for x, n in
                       ((p, 'p'), (radius, 'radius'), (bound, 'bound'))]
    if not 0 < p <= 1 or radius < 1 or bound < 0 or axis not in ('P', 'Q'):
        raise ValueError('valid compact test and scale required')
    cells = (2*radius+1)**3
    difference = p*p*bound*(2 if axis == 'Q' else 1)
    edges = 6*cells/(p**3)
    return {'field': axis, 'incident_edge_upper': edges,
            'relative_angle_upper': difference,
            'quadratic_coefficient_norm_upper': edges*difference**2/(2*p),
            'scaled_BCH_remainder_norm_upper': edges*difference**3/(12*p),
            'bound_is_uniform_over_all_source_states_and_references': True}


def bond_kick(degree, coefficient_max, coefficient_l1):
    if type(degree) is not int or degree < 1:
        raise ValueError('positive maximum graph degree required')
    maximum, total = [rational(x, n) for x, n in
                      ((coefficient_max, 'maximum'), (coefficient_l1, 'l1'))]
    if not 0 <= maximum <= total:
        raise ValueError('ordered nonnegative coefficient norms required')
    count = 3*degree*(2*degree-1)
    return {'local_overlap_count_upper': count,
            'double_commutator_norm_over_J_upper': count*maximum*total,
            'energy_kick_squared_over_J_upper': F(count, 2)*maximum*total}


def phase_transport(p, amplitude, l1_bound, energy_bound, time, phase=F(1)):
    """An exact squared-norm bound using rational square-root majorants.

    The caller supplies a rational amplitude a with |z_b|<=a p^3 and
    sum|z_b|<=a*l1_bound. The returned expression bounds sqrt(99*l1_bound
    *energy_bound) by an outward integer square root of its rational value.
    """
    from math import isqrt
    p, amplitude, l1_bound, energy_bound, time, phase = [rational(x, n) for x, n in
        ((p,'p'), (amplitude,'amplitude'), (l1_bound,'l1_bound'),
         (energy_bound,'energy_bound'), (time,'time'), (phase,'phase'))]
    if p <= 0 or min(amplitude, l1_bound, energy_bound) < 0:
        raise ValueError('positive scale and nonnegative bounds required')
    radicand = 99*l1_bound*energy_bound
    upper = isqrt(radicand.numerator//radicand.denominator)
    if upper*upper < radicand:
        upper += 1
    estimate = (4*abs(time*phase)*amplitude*upper*p
                +990*abs(time)*phase**2*amplitude**2*l1_bound*p*p)
    return {'scaled_time': time,
            'initial_energy_over_Jp_upper': energy_bound,
            'square_root_majorant': F(upper),
            'phase_transport_squared_norm_upper': min(F(4), estimate),
            'uncapped_squared_norm_upper': estimate}


def ground_energy_cube(side):
    """Actual shifted negative-projector bond normalization, all bounds / J."""
    if type(side) is not int or side < 1:
        raise ValueError('positive integer cube side required')
    volume = side**3
    internal = 3*side*side*(side-1)
    crossing = 6*side*side
    return {'side':side,'sites':volume,'internal_edges':internal,'crossing_edges':crossing,
            'local_ground_energy_excess_over_J_upper':F(crossing),
            'normalized_ground_energy_fluctuation_squared_over_J2_upper':F(internal*crossing,volume**2),
            'finite_ground_density_above_thermodynamic_density_over_J_upper':F(3,side),
            'forward_boundary_density_norm_over_J_upper':F(3,side),
            'ergodicity_assumed':False,
            'individual_transverse_bond_limits_proved_scalar':False}


def serialize(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(item) for item in value]
    return value
