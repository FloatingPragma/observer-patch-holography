"""Exact rational windows and independent spectral countercontrols."""
from fractions import Fraction as F


def rational(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, F)):
        raise ValueError(f'{name} must be an exact integer or Fraction')
    return F(value)


def parameters(m0, radius):
    m = rational(m0, 'm0')
    r = rational(radius, 'radius')
    if not 0 < m <= F(1, 2) or r <= 0:
        raise ValueError('require 0<m0<=1/2 and radius>0')
    return m, r


def stiffness_interval(m0):
    m, _ = parameters(m0, 1)
    return m*m/2, F(1, 2)


def window(channel, m0, radius):
    m, r = parameters(m0, radius)
    if channel == 'P':
        total, inverse, first = m*m*r/24, F(1, 12), r*r/4
        low, high = m*m*r/8, 24*r/(m*m)
        promised = m*m*r/48
    elif channel == 'Q':
        total, inverse, first = m*m/r, 4/(r*r), F(3)
        low, high = m*m*r/16, 12*r/(m*m)
        promised = m*m/(2*r)
    else:
        raise ValueError('channel must be P or Q')
    lower = total - low*inverse - first/high
    if not 0 < low < high or lower != promised or lower <= 0:
        raise ArithmeticError('positive spectral window failed')
    return {'channel': channel, 'm0': m, 'momentum_inner': r/2,
            'momentum_outer': r, 'total_weight_lower': total,
            'inverse_moment_upper': inverse, 'first_moment_upper': first,
            'energy_lower': low, 'energy_upper': high,
            'retained_weight_lower': lower}


def tail_bound(first, inverse, low, high):
    first, inverse, low, high = [rational(x, n) for x, n in
                                zip((first, inverse, low, high), ('first', 'inverse', 'low', 'high'))]
    if first < 0 or inverse < 0 or not 0 < low < high:
        raise ValueError('invalid positive-measure bounds')
    return low*inverse + first/high


def escaping_energy(n):
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError('n must be an integer at least two')
    weights = (F(n-1, n), F(1, n))
    energies = (F(1), F(n))
    return {'mass': sum(weights),
            'first': sum(w*e for w, e in zip(weights, energies)),
            'inverse': sum(w/e for w, e in zip(weights, energies)),
            'mass_at_escaping_energy': weights[1],
            'first_moment_at_escaping_energy': weights[1]*energies[1],
            'weak_limit_first': F(1)}


def continuous_countercontrol(m0, radius):
    """Uniform energy density on [r,2r]; no atoms, not a source-state claim.

    Its inverse moment is weight*log(2)/r < weight/r.  We retain the
    rational upper bound, rather than serializing an approximate logarithm.
    """
    m, r = parameters(m0, radius)
    weight = m*m*r/24
    first = F(3, 2)*weight*r
    inverse_upper = weight/r
    b = window('P', m, r)
    return {'weight': weight, 'energy_lower': r, 'energy_upper': 2*r,
            'first_moment': first, 'inverse_moment_upper': inverse_upper,
            'atom_free': True,
            'obeys_coarse_P_bounds': (first <= b['first_moment_upper']
                                     and inverse_upper <= b['inverse_moment_upper']
                                     and weight >= b['total_weight_lower']),
            'is_a_Heisenberg_state_construction': False}


def rotation_remainder_coefficient():
    # ||C||_1<=3/4; 2t+t^2<=3t; (|da|+|db|)^2<=2(da^2+db^2).
    coefficient = F(1, 2)*F(3, 4)*3*2
    assert coefficient <= 3
    return coefficient


def serialize(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k: serialize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(v) for v in value]
    return value
