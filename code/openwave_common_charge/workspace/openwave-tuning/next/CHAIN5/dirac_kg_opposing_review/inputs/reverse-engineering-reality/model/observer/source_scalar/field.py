"""The declared two-site scalar action and its independently specified targets.

The mass, spring and model coordinate are construction inputs. These helpers
define the target; the execution compiles it into the existing source pulses.
The occupation projection uses physical matrix elements of Q^2 and P^2,
which differ from the squares of the finite Q and P matrices at the boundary.
"""
from fractions import Fraction as F
import math
from numbers import Integral

import numpy as np

KAPPA = F(9, 32)


def dimension(d):
    if isinstance(d, (bool, np.bool_)) or not isinstance(d, Integral) or d not in (2, 4, 8, 16):
        raise ValueError("onsite occupation dimension must be 2, 4, 8 or 16")
    return int(d)


def quadratures(d):
    d = dimension(d)
    annihilation = np.diag(np.sqrt(np.arange(1, d)), 1)
    return (annihilation+annihilation.T)/np.sqrt(2), (annihilation-annihilation.T)/(1j*np.sqrt(2))


def projected_squares(d):
    q, p = quadratures(d)
    boundary = np.zeros((d, d))
    boundary[-1, -1] = d/2
    return q@q+boundary, (p@p).real+boundary


def onsite_hamiltonian(d):
    qq, pp = projected_squares(d)
    return .5*(pp+(1+float(KAPPA))*qq)


def physical_hamiltonian(d):
    """Pi H Pi in little-endian site order: |n1,n0>."""
    q, _ = quadratures(d)
    onsite = onsite_hamiltonian(d)
    return np.kron(np.eye(d), onsite)+np.kron(onsite, np.eye(d))-float(KAPPA)*np.kron(q, q)


def _root_upper(value, scale=10**15):
    integer = math.isqrt(value.numerator*scale*scale//value.denominator)
    if F(integer, scale)**2 < value:
        integer += 1
    return F(integer, scale)


def vacuum_vector(d):
    """Finite normalized projection of the analytic infinite coupled vacuum.

The common normal mode is blank; the difference mode has omega=5/4 and
squeezing parameter 1/9. Retain total occupation through d-2. Rotating that
sector by pi/4 is an exact embedding in the two d-level registers. The
returned vector is a compiler target, never injected as the executed state.
"""
    d = dimension(d)
    target = np.zeros(d*d)
    probability0 = 4*math.sqrt(5)/9
    recurrence, total = F(1), F(0)
    for pair in range(d//2):
        number = 2*pair
        coefficient = (-1)**pair*math.sqrt(probability0*float(recurrence))
        for n0 in range(number+1):
            n1 = number-n0
            target[n0+d*n1] += coefficient*(-1)**n1*math.sqrt(math.comb(number, n0)/2**number)
        total += recurrence
        recurrence *= F(2*pair+1, 81*(2*pair+2))
    retained = float(np.dot(target, target))
    target /= math.sqrt(retained)
    scale = 10**30
    sqrt5_lower = F(math.isqrt(5*scale*scale), scale)
    tail_upper = 1-F(4, 9)*sqrt5_lower*total
    return target, {"dimension": d, "maximum_total_occupation": d-2,
                    "retained_probability": retained,
                    "tail_probability_upper": str(tail_upper),
                    "infinite_vacuum_trace_distance_upper": str(_root_upper(tail_upper)),
                    "bound_scope": "initial preparation only; no propagation, truncation or continuum error transfer"}


def ideal_prediction(time=1., current=1):
    """Infinite two-oscillator ground-state reference, not an outcome oracle."""
    if type(current) is not int or current not in (-1, 0, 1) or not math.isfinite(time):
        raise ValueError("finite model time and native current required")
    mean = current*(math.sin(time)-.8*math.sin(1.25*time))/2
    variance = 9/20
    return {"mean": mean, "variance": variance,
            "plus_probability": (1+math.exp(-variance/2)*math.sin(mean))/2}


def canonical_eigenbasis(matrix):
    values, vectors = np.linalg.eigh(matrix)
    # Remove numerical parity dust and use the first supported coordinate,
    # not an argmax whose tied components can flip across LAPACK versions.
    vectors[abs(vectors) < 1e-14] = 0
    for column in vectors.T:
        pivot = int(np.flatnonzero(abs(column) > 1e-12)[0])
        if column[pivot] < 0:
            column *= -1
    return values, vectors
