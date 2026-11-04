"""Exact infinite-Jacobi eigenvalue enclosures from untrusted rational vectors.

The fixed witness uses v=kappa=a=1. Spectral counts and residuals are exact
rational calculations; the eigensolver that proposed the vectors is not
part of the certificate.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path

from tools import strict_json


WITNESS = Path(__file__).with_name('spectral_witness.json')


def sqrt_upper(value):
    if value < 0:
        raise ValueError('nonnegative square required')
    root = isqrt(value.numerator*value.denominator)
    if root*root != value.numerator*value.denominator:
        root += 1
    return F(root, value.denominator)


def tridiagonal(radius):
    if type(radius) is not int or radius < 2:
        raise ValueError('integer radius at least two required')
    positions = range(-radius, radius+1)
    diagonal = [F(n*(-1 if n % 2 else 1)) for n in positions]
    off = [F(-1 if n % 2 else 1, 2) for n in range(-radius, radius)]
    return diagonal, off


def sturm_count(diagonal, off, point):
    """Negative inertia of B_finite-point, by exact tridiagonal LDL."""
    pivot, count = diagonal[0]-point, 0
    for index in range(len(diagonal)):
        if index:
            pivot = diagonal[index]-point-off[index-1]**2/pivot
        if pivot == 0:
            raise ValueError('unresolved exact zero LDL pivot; choose another endpoint')
        count += int(pivot < 0)
    return count


def protected_endpoint(diagonal, off, radius, point):
    """No eigenvalue crosses this endpoint when the two tails are attached."""
    margin = F(1, 50)
    if abs(point)+margin >= radius:
        raise ValueError('endpoint outside the tail-resolvent gap')
    lower = sturm_count(diagonal, off, point-margin)
    upper = sturm_count(diagonal, off, point+margin)
    if lower != upper:
        raise ValueError('finite endpoint lacks its exact spectral separation')
    self_energy = F(1, 4)/(radius-abs(point))
    if self_energy >= margin:
        raise ValueError('tail self-energy can cross the protected endpoint')
    return lower


def residual_bound(diagonal, off, center, vector):
    if len(vector) != len(diagonal):
        raise ValueError('complete finite witness vector required')
    norm_squared = sum((x*x for x in vector), F())
    if norm_squared == 0:
        raise ValueError('nonzero witness vector required')
    residual = []
    for index, value in enumerate(vector):
        item = (diagonal[index]-center)*value
        if index:
            item += off[index-1]*vector[index-1]
        if index+1 < len(vector):
            item += off[index]*vector[index+1]
        residual.append(item)
    # The zero extension has two additional, retained boundary residuals.
    squared = sum((x*x for x in residual), F())+(vector[0]**2+vector[-1]**2)/4
    return sqrt_upper(squared/norm_squared)


def outward_decimal(value, digits=12, upper=False):
    scaled = value*10**digits
    integer = -(-scaled.numerator//scaled.denominator) if upper else scaled.numerator//scaled.denominator
    return F(integer, 10**digits)


def certify(value):
    if (value.get('schema') != 'oph.defect_jacobi_rational_witness.v1'
            or value.get('radius') != 32
            or [value.get(key) for key in ('velocity', 'gradient', 'spacing')] != ['1']*3):
        raise ValueError('the stated source Jacobi witness is required')
    radius = value['radius']
    diagonal, off = tridiagonal(radius)
    low = protected_endpoint(diagonal, off, radius, F(-1, 2))
    base = protected_endpoint(diagonal, off, radius, F(1, 2))
    if base-low != 1:
        raise ValueError('the central spectral window must contain only the exact zero')
    rows = []
    candidates = value.get('candidates', [])
    if len(candidates) != 6:
        raise ValueError('all six advertised positive levels are required')
    for index, candidate in enumerate(candidates, 1):
        if candidate['positive_level'] != index:
            raise ValueError('positive-level ordering changed')
        center = F(candidate['center'])
        left, right = center-F(1, 20), center+F(1, 20)
        first = protected_endpoint(diagonal, off, radius, left)
        last = protected_endpoint(diagonal, off, radius, right)
        if last-first != 1 or first-base != index-1:
            raise ValueError('the complete spectral count does not isolate this level')
        error = residual_bound(diagonal, off, center, [F(x) for x in candidate['vector']])
        if error >= F(1, 10**11):
            raise ValueError('infinite-source residual is too large')
        lower = outward_decimal(center-error)
        upper = outward_decimal(center+error, upper=True)
        if not left < lower < upper < right:
            raise ValueError('residual enclosure escaped the isolated source level')
        rows.append({'positive_level': index,
                     'energy_interval': [str(lower), str(upper)],
                     'infinite_source_residual_upper': str(error),
                     'spectral_count_inside_coarse_window': 1,
                     'positive_levels_below_coarse_window': index-1})
    return {'velocity': '1', 'gradient': '1', 'spacing': '1',
            'finite_witness_radius': radius,
            'infinite_tail_retained_by_resolvent': True,
            'zero_eigenvalue_multiplicity_one_dimension': 1,
            'negative_levels_are_exact_reflections': True,
            'positive_levels': rows,
            'three_dimensional_levels': '+/-sqrt(b_i^2+b_j^2+b_k^2), b_0=0; retain every sign and permutation',
            'three_dimensional_zero_multiplicity': 2,
            'physical_units_calibrated': False}


def certificate():
    return certify(strict_json.load(WITNESS))
