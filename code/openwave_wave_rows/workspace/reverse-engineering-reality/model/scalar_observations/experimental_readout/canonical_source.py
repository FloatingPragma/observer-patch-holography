"""Fixed-arithmetic word generation for the executable oscillator readout.

This module chooses the requested angles, not their executed matrices. A private
80-decimal-digit mpmath context evaluates the original binary64 parameters,
oscillator matrices and complete Gray elimination. Only coefficients passed to
the existing source compiler become floats. The independently executed source
words still have to pass the SciPy finite-matrix targets in executable_source.
No precision or physical-error certificate follows from this arithmetic policy.
"""
from functools import lru_cache
import math
from numbers import Integral, Real

import mpmath
import numpy as np

from ...observer.source_scalar import compiler as c


DECIMAL_DIGITS = 80


def policy():
    return {
        'arithmetic': 'private fixed-precision mpmath context',
        'decimal_digits': DECIMAL_DIGITS,
        'input_real_parameters': 'exact ratios of the supplied binary64 values',
        'angle_conversion': 'binary64 at the existing source compiler boundary',
        'basis_sign': 'first entry with magnitude greater than 1e-12 is positive',
        'gray_elimination': 'complete fixed schedule including zero rotations',
        'source_primitive_or_execution_changed': False,
        'outward_roundoff_or_physical_error_certificate': False,
    }


def _context():
    context = mpmath.mp.clone()
    context.dps = DECIMAL_DIGITS
    return context


def _dimension(value):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
            or value < 2 or value > 64 or value & (value - 1)):
        raise ValueError('a power-of-two dimension from 2 through 64 is required')
    return int(value)


def _parameter(context, value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f'{name} must be a finite real binary64 parameter')
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f'{name} must be a finite real binary64 parameter')
    numerator, denominator = value.as_integer_ratio()
    return context.mpf(numerator) / denominator


def _ladder(context, dimension):
    result = context.matrix(dimension)
    for index in range(1, dimension):
        result[index - 1, index] = context.sqrt(index)
    return result


def _diagonal_word(context, phases):
    """Complete Walsh expansion, rounded only at each Pauli request."""
    bits = len(phases).bit_length() - 1
    word = []
    for mask in range(len(phases)):
        coefficient = context.fsum(value * (-1) ** ((index & mask).bit_count())
                                   for index, value in enumerate(phases)) / len(phases)
        label = ['Z' if mask >> bit & 1 else 'I' for bit in range(bits)]
        word.extend(c.pauli_rotation(label, float(coefficient)))
    return word


def _orthogonal_word(context, value):
    """The existing complete Gray-Givens schedule in fixed arithmetic."""
    dimension = value.rows
    bits = dimension.bit_length() - 1
    gray = [index ^ (index >> 1) for index in range(dimension)]
    residual = context.matrix([[value[row, column] for column in gray] for row in gray])
    rotations = []
    for column in range(dimension - 1):
        for row in range(dimension - 1, column, -1):
            x, y = residual[row - 1, column], residual[row, column]
            angle = context.zero if x == 0 and y == 0 else context.atan2(y, x)
            cosine, sine = context.cos(angle), context.sin(angle)
            for index in range(dimension):
                first, second = residual[row - 1, index], residual[row, index]
                residual[row - 1, index] = cosine * first + sine * second
                residual[row, index] = -sine * first + cosine * second
            rotations.append((gray[row - 1], gray[row], angle))
    signs = [1 if residual[index, index] >= 0 else -1 for index in range(dimension)]
    if max(abs(residual[row, column] - (signs[row] if row == column else 0))
           for row in range(dimension) for column in range(dimension)) > context.mpf('1e-60'):
        raise ArithmeticError('canonical Gray elimination did not retain the complete basis')
    phases = [context.zero] * dimension
    for index, sign in enumerate(signs):
        phases[gray[index]] = context.zero if sign == 1 else context.pi
    word = _diagonal_word(context, phases)
    for first, second, angle in reversed(rotations):
        difference = first ^ second
        target = difference.bit_length() - 1
        orientation = (-1) ** ((first >> target) & 1)
        for mask in range(1 << bits):
            if mask & difference:
                continue
            label = ['Y' if bit == target else 'Z' if mask >> bit & 1 else 'I'
                     for bit in range(bits)]
            coefficient = (-orientation * angle * (-1) ** ((first & mask).bit_count())
                           / (1 << (bits - 1)))
            word.extend(c.pauli_rotation(label, float(coefficient)))
    return word


def displacement_word(dimension, alpha_real):
    dimension = _dimension(dimension)
    context = _context()
    alpha = _parameter(context, alpha_real, 'real displacement')
    ladder = _ladder(context, dimension)
    matrix = context.expm(alpha * (ladder.T - ladder))
    return _orthogonal_word(context, matrix)


@lru_cache(maxsize=8)
def _coordinate_data(dimension):
    context = _context()
    ladder = _ladder(context, dimension)
    values, vectors = context.eigsy((ladder + ladder.T) / 2)
    for column in range(dimension):
        pivot = next(row for row in range(dimension)
                     if abs(vectors[row, column]) > context.mpf('1e-12'))
        if vectors[pivot, column] < 0:
            for row in range(dimension):
                vectors[row, column] = -vectors[row, column]
    return tuple(values), tuple(tuple(vectors[row, column] for column in range(dimension))
                                for row in range(dimension))


def coordinate_word(dimension):
    dimension = _dimension(dimension)
    context = _context()
    _, vectors = _coordinate_data(dimension)
    return _orthogonal_word(context, context.matrix(vectors))


def quadrature_words(dimension, u, v, theta):
    """Return canonical phi/radius and the rotation and detector source words."""
    dimension = _dimension(dimension)
    context = _context()
    u, v, theta = (_parameter(context, value, name)
                   for value, name in ((u, 'u'), (v, 'v'), (theta, 'theta')))
    radius = context.sqrt(u*u + v*v)
    # The negative real axis uses +pi, independent of binary64 signed zero.
    phi = context.zero if radius == 0 else context.atan2(-u, v)
    rotation = _diagonal_word(context, [-phi * index for index in range(dimension)])
    values, _ = _coordinate_data(dimension)
    angles = [radius * context.mpf(value) + theta / 2 for value in values]
    detector = _diagonal_word(context, [-value for value in angles] + angles)
    return float(phi), float(radius), rotation, detector
