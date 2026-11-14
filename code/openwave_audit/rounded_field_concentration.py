"""Finite-time concentration for the declared native rounded shear field."""
from fractions import Fraction as F


def integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError('literal integer in the admitted domain required')
    return value


def coefficients(a=2, b=3):
    return F(1, 2**integer(a)), F(1, 2**integer(b))


def heating_bound(trace, a=2, b=3):
    """One-step conditional shadow-energy variance bound, exact rational."""
    integer(trace)
    alpha, beta = coefficients(a, b)
    return F(trace, 4)*(beta*int(a > 0)+alpha*int(b > 0))


def deviation_bound(steps, trace, amplitude=1, replicas=1, a=2, b=3):
    """E Q(z/A - deterministic mean/A); independent replicas are averaged."""
    return integer(steps)*heating_bound(trace, a, b)/(integer(amplitude, 1)**2*integer(replicas, 1))


def uniform_tail(steps, trace, threshold_squared, amplitude=1, replicas=1,
                 a=2, b=3):
    """P(max_{0..steps} Q(error)>threshold_squared), under positive Q."""
    if type(threshold_squared) is not F or threshold_squared <= 0:
        raise ValueError('positive exact rational threshold squared required')
    return min(F(1), deviation_bound(steps, trace, amplitude, replicas, a, b)/threshold_squared)


def scalar_stability():
    """The frozen actual-face uniform bound gives diagonal norm <=6Q."""
    alpha, beta = coefficients()
    curvature = F(1504, 17)
    assert alpha*beta*curvature == F(47, 17)
    # sqrt(alpha beta Lambda)/2 < 5/6, checked without floating roots.
    assert alpha*beta*curvature/4 < F(25, 36)
    return {'curvature_upper': curvature, 'Q_to_diagonal_norm': F(6),
            'transverse_speed_squared_cells_per_macro': F(2, 3)*alpha*beta}
