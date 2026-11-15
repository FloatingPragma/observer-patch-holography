"""Directed comparisons of reported binary arrays with an exact proof object.

Every original array component is converted to its exact rational value.  No
source execution, numerical eigensolver, normalization or probability cap is
used.  The comparison includes the polynomial's entire occupation support.

An Enclosure is a trusted, in-process result of smooth_history.enclose.  Its
existing integrity check detects mutation; its public seal does not authenticate
a theorem or a caller-supplied error radius.  Receipt producers must recompute
the enclosure and establish the provenance of the reported arrays themselves.
These functions certify the arrays supplied, not that they came from a source
backend.  Reference systems are outside this comparison API.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import math

import numpy as np

from . import protocol, smooth_history as smooth


ROOT_BITS = 192


def _sqrt_interval(value):
    """Exact rational endpoints around sqrt(value), including tiny values."""
    if type(value) is not F or value < 0:
        raise ValueError('a nonnegative exact rational radicand is required')
    return _cached_sqrt_interval(value)


@lru_cache(maxsize=256)
def _cached_sqrt_interval(value):
    """Memoize only values already checked at the caller boundary."""
    numerator, denominator = math.isqrt(value.numerator), math.isqrt(value.denominator)
    if numerator*numerator == value.numerator and denominator*denominator == value.denominator:
        root = F(numerator, denominator)
        return root, root
    # Scale with the original rational's magnitude, rather than first rounding
    # a small radicand to a binary float or imposing an absolute error floor.
    exponent = (value.numerator.bit_length()-value.denominator.bit_length())//2
    scale = 1 << max(0, ROOT_BITS-exponent)
    lower = math.isqrt(value.numerator*scale*scale//value.denominator)
    return F(lower, scale), F(lower+1, scale)


def _squared(interval):
    lower, upper = interval
    return (F(0) if lower <= 0 <= upper else min(lower*lower, upper*upper),
            max(lower*lower, upper*upper))


def _radical_sum(coefficients):
    """Enclose sum_w coefficients[w]*sqrt(w), after exact cancellations."""
    lower = upper = F(0)
    for weight, coefficient in coefficients.items():
        lo, hi = _sqrt_interval(F(weight))
        if coefficient < 0:
            lo, hi = hi, lo
        lower += coefficient*lo
        upper += coefficient*hi
    return lower, upper


def _component(value):
    # Widening binary16/32 components to binary64 is exact; extended-precision
    # dtypes are refused before this conversion.  Integer entries stay exact.
    if isinstance(value, (int, np.integer)):
        return F(int(value)), F(0)
    return F(float(value.real)), F(float(value.imag))


def _array(value, shape):
    if np.ma.isMaskedArray(value) or not isinstance(value, np.ndarray):
        raise ValueError('an unmasked original numeric ndarray is required')
    kind, width = value.dtype.kind, value.dtype.itemsize
    if (value.shape != shape or kind not in 'iufc'
            or (kind == 'f' and width > 8) or (kind == 'c' and width > 16)
            or not np.isfinite(value).all()):
        raise ValueError('finite original numeric arrays with the declared shape and no reference axis are required')
    return value.copy()


def _inputs(program, arrays, enclosure, *, density):
    program = protocol.validate_program(program)
    smooth._validate(enclosure)
    target_program = protocol.validate_program(enclosure.program)
    target = lambda p: (p.operations, p.coupling, p.native_current)
    if target(program) != target(target_program):
        raise ValueError('the exact operations, coupling and native current must match the enclosure')
    if (enclosure.dimension is not None and
            (type(enclosure.dimension) is not int
             or enclosure.dimension != program.dimension
             or target_program.dimension != program.dimension)):
        raise ValueError('a finite enclosure must use the reported work cutoff')
    if (type(enclosure.normalization) is not F or enclosure.normalization <= 0
            or any(type(error) is not F or error < 0
                   for error in (enclosure.vector_error, enclosure.preparation_error))):
        raise ValueError('positive exact normalization and nonnegative proof errors required')
    records = tuple(product((0, 1), repeat=sum(op.kind in protocol.READS for op in program.operations)))
    if (type(arrays) is not dict or any(type(record) is not tuple
            or len(record) != len(records[0])
            or any(type(bit) is not int or bit not in (0, 1) for bit in record)
            for record in arrays) or set(arrays) != set(records)):
        raise ValueError('every complete binary record, including zero outcomes, is required')
    shape = (program.work_dimension,)* (2 if density else 1)
    copied = {record: _array(arrays[record], shape) for record in records}
    for polynomial in enclosure.branches.values():
        if type(polynomial.denominator) is not int or polynomial.denominator <= 0:
            raise ValueError('a positive exact polynomial denominator is required')
        for key, coefficient in polynomial.coefficients.items():
            if (type(key) is not tuple or len(key) != 3
                    or any(type(index) is not int or index < 0 for index in key)
                    or key[2] > 3 or type(coefficient) is not tuple or len(coefficient) != 2
                    or any(type(part) is not int for part in coefficient)
                    or (enclosure.dimension is not None
                        and max(key[:2]) >= enclosure.dimension)):
                raise ValueError('valid exact Fock polynomial coefficients are required')
    return program, copied, records


def _polynomial(polynomial, normalization, dimension):
    """Full norm and the part of its factorial-basis vector in the work code."""
    divisor = polynomial.denominator**2*normalization
    norm = F(0)
    retained = {}
    outside = False
    for (n0, n1, gauge), (real, imaginary) in polynomial.coefficients.items():
        weight = math.factorial(n0)*math.factorial(n1)
        norm += F(weight*(real*real+imaginary*imaginary), 1)/divisor
        if n0 < dimension and n1 < dimension:
            retained[gauge+4*(n0+dimension*n1)] = (real, imaginary, weight)
        elif real or imaginary:
            outside = True
    return norm, retained, divisor, outside


def _nonnegative(lower, upper):
    if upper < 0:
        raise ArithmeticError('an exact squared-norm enclosure is inconsistent')
    return max(F(0), lower), upper


def _finish(records, rows, enclosure):
    lower = sum((row['distance_lower'] for row in rows.values()), F(0))
    upper = sum((row['distance_upper'] for row in rows.values()), F(0))
    error = enclosure.output_error
    return {
        'distance_lower': max(F(0), lower-error),
        'distance_upper': upper+error,
        'polynomial_distance_lower': lower,
        'polynomial_distance_upper': upper,
        'enclosure_error': error,
        'source_mass_real': sum((row['source_mass_real'] for row in rows.values()), F(0)),
        'source_mass_imag': sum((row['source_mass_imag'] for row in rows.values()), F(0)),
        'records': records,
        'per_record': rows,
    }


def compare_vectors(program, vectors, enclosure):
    """Bound half trace norm of complete reported rank-one CQ outputs.

    Each vector must have shape (4*d*d,), with index g+4*(n0+d*n1).
    It need not be normalized or arise from the declared preparation.  The
    target is the enclosure's declared prepared history, so a different input
    is measured as an output discrepancy, not silently granted that theorem.
    """
    program, vectors, records = _inputs(program, vectors, enclosure, density=False)
    rows = {}
    for record in records:
        values = tuple(_component(value) for value in vectors[record])
        mass = sum((real*real+imag*imag for real, imag in values), F(0))
        norm, retained, divisor, _ = _polynomial(
            enclosure.branches[record], enclosure.normalization, program.dimension)
        real_terms, imag_terms = {}, {}
        for index, (real, imag, weight) in retained.items():
            vr, vi = values[index]
            real_terms[weight] = real_terms.get(weight, F(0))+vr*real+vi*imag
            imag_terms[weight] = imag_terms.get(weight, F(0))+vr*imag-vi*real
        real = _squared(_radical_sum(real_terms))
        imag = _squared(_radical_sum(imag_terms))
        overlap = (real[0]+imag[0])/divisor, (real[1]+imag[1])/divisor
        squared = _nonnegative((mass+norm)**2-4*overlap[1],
                               (mass+norm)**2-4*overlap[0])
        rows[record] = {
            'distance_lower': _sqrt_interval(squared[0])[0]/2,
            'distance_upper': _sqrt_interval(squared[1])[1]/2,
            'source_mass_real': mass, 'source_mass_imag': F(0),
            'polynomial_mass': norm,
        }
    return _finish(records, rows, enclosure)


def compare_densities(program, densities, enclosure):
    """Bound complete matrix-output half nuclear norm against pure polynomials.

    Reported matrices need not be Hermitian, positive or trace one.  For each
    block A=rho-|p><p|, compute ||A||_F squared from exact components and
    ||p|| to include all outside-code support.  Then ||A||_F <= ||A||_1 <=
    sqrt(D+1)*||A||_F; D suffices when p is wholly in the D-dimensional code.
    These are matrix-error bounds, not assertions that rho is a quantum state.
    """
    program, densities, records = _inputs(program, densities, enclosure, density=True)
    rows = {}
    for record in records:
        values = tuple(tuple(_component(value) for value in row) for row in densities[record])
        real_mass = sum((values[i][i][0] for i in range(program.work_dimension)), F(0))
        imag_mass = sum((values[i][i][1] for i in range(program.work_dimension)), F(0))
        norm, retained, divisor, outside = _polynomial(
            enclosure.branches[record], enclosure.normalization, program.dimension)
        frobenius = sum((real*real+imag*imag for row in values for real, imag in row), F(0))
        terms = {}
        for i, (ar, ai, wi) in retained.items():
            for j, (br, bi, wj) in retained.items():
                rr, ri = values[i][j]
                # Re(conj(c_i) * rho_ij * c_j), before either radical is rounded.
                coefficient = (ar*rr+ai*ri)*br-(ar*ri-ai*rr)*bi
                if coefficient:
                    weight = wi*wj
                    terms[weight] = terms.get(weight, F(0))+coefficient
        cross = _radical_sum(terms)
        squared = _nonnegative(frobenius+norm*norm-2*cross[1]/divisor,
                               frobenius+norm*norm-2*cross[0]/divisor)
        rank = program.work_dimension+int(outside)
        trace_lower = _sqrt_interval((real_mass-norm)**2+imag_mass**2)[0]/2
        rows[record] = {
            'distance_lower': max(_sqrt_interval(squared[0])[0]/2, trace_lower),
            'distance_upper': _sqrt_interval(rank*squared[1])[1]/2,
            'source_mass_real': real_mass, 'source_mass_imag': imag_mass,
            'polynomial_mass': norm, 'rank_upper': rank,
        }
    return _finish(records, rows, enclosure)
