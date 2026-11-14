"""Exact all-input Strang error bound for the declared finite scalar code.

Unnormalized number states make every quadratic matrix entry rational.
Factorial metric weights restore the normalized Frobenius norm exactly.
No floating propagator, fitted convergence rate or physical service is used.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import math


def _multiply(left, right):
    right_rows = defaultdict(list)
    for (row, column), value in right.items():
        right_rows[row].append((column, value))
    result = defaultdict(F)
    for (row, middle), value in left.items():
        for column, other in right_rows[middle]:
            result[row, column] += value*other
    return {index: value for index, value in result.items() if value}


def _commutator(left, right):
    result = defaultdict(F, _multiply(left, right))
    for index, value in _multiply(right, left).items():
        result[index] -= value
    return {index: value for index, value in result.items() if value}


def _raw_split(d):
    # Coordinates of (a0^dag)^n0 (a1^dag)^n1 |00>; not orthonormal.
    # The occupation boundary is the compression of the physical operator.
    onsite, bond = defaultdict(F), defaultdict(F)
    kappa = F(9, 32)
    for n1 in range(d):
        for n0 in range(d):
            column = n0+d*n1
            for occupied, stride in ((n0, 1), (n1, d)):
                onsite[column, column] += (1+kappa/2)*(occupied+F(1, 2))
                if occupied+2 < d:
                    onsite[column+2*stride, column] += kappa/4
                if occupied >= 2:
                    onsite[column-2*stride, column] += kappa*occupied*(occupied-1)/4
            for shift0, factor0 in ((-1, n0), (1, 1)):
                for shift1, factor1 in ((-1, n1), (1, 1)):
                    if 0 <= n0+shift0 < d and 0 <= n1+shift1 < d:
                        row = column+shift0+d*shift1
                        bond[row, column] -= kappa*factor0*factor1/2
    return dict(onsite), {index: value for index, value in bond.items() if value}


def _frobenius_squared(matrix, d):
    weights = [math.factorial(index % d)*math.factorial(index//d) for index in range(d*d)]
    return sum((value*value*F(weights[row], weights[column])
                for (row, column), value in matrix.items()), F())


def root_upper(value, scale=10**15):
    """Smallest multiple of 1/scale whose square dominates a rational."""
    if not isinstance(value, F) or value < 0 or type(scale) is not int or scale <= 0:
        raise ValueError("nonnegative exact rational and positive integer scale required")
    integer = math.isqrt(value.numerator*scale*scale//value.denominator)
    if F(integer, scale)**2 < value:
        integer += 1
    return F(integer, scale)


@lru_cache(maxsize=4)
def nested_commutator_bounds(d):
    if type(d) is not int or d not in (2, 4, 8, 16):
        raise ValueError("dimension must be 2, 4, 8 or 16")
    onsite, bond = _raw_split(d)
    commutator = _commutator(onsite, bond)
    first = _frobenius_squared(_commutator(onsite, commutator), d)
    second = _frobenius_squared(_commutator(bond, commutator), d)
    return first, second


def split_error_bound(d, steps, time=F(1)):
    """Uniform operator-norm, pure-state trace-distance and effect bound.

The target is exp(-i*time*H_d), not the infinite oscillator or continuum.
The same bound remains valid on arbitrary reference-entangled inputs.
"""
    if type(steps) is not int or steps <= 0 or not isinstance(time, F) or time < 0:
        raise ValueError("positive integer steps and nonnegative exact rational time required")
    first, second = nested_commutator_bounds(d)
    first_upper, second_upper = root_upper(first), root_upper(second)
    epsilon = time**3*(first_upper/24+second_upper/12)/steps**2
    return {
        "dimension": d,
        "steps": steps,
        "model_time": str(time),
        "onsite_nested_commutator_frobenius_squared": str(first),
        "bond_nested_commutator_frobenius_squared": str(second),
        "onsite_nested_commutator_norm_upper": str(first_upper),
        "bond_nested_commutator_norm_upper": str(second_upper),
        "unitary_operator_error_upper": str(epsilon),
        "state_trace_distance_upper": str(min(F(1), epsilon)),
        "detector_probability_error_upper": str(min(F(1), epsilon)),
        "scope": "ideal split versus continuous evolution of the same compressed finite Hamiltonian; no infinite-cutoff or physical-service bound",
    }
