"""Distribution-valued limits of the actual ancestral native field.

The comparison sphere supplies a test-function topology, not a physical
metric or a target Green operator. The proof is in OCCUPATION_FIELD.md.
Observation count tends to infinity before spatial refinement.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.sparse.linalg import splu

from . import uniform_green as uniform
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*uniform.SOURCE_PATHS,
    'model/physical_space_time/continuum_field.py',
    'model/physical_space_time/__init__.py',
    # The package initializers execute these modules during a fresh CLI run.
    'model/__init__.py', 'model/common_source/__init__.py',
    'model/common_source/experiment.py', 'model/common_source/gauge_matter.py',
    'model/common_source/native_join.py', 'model/common_source/records.py',
    'model/common_source/scalar_coupling.py', 'model/observer/__init__.py',
    'model/observer/p_response/__init__.py', 'model/observer/p_response/experiment.py',
    'model/observer/scalar_response/__init__.py', 'model/observer/source_scalar/__init__.py',
    'model/observer/source_scalar/bounds.py', 'model/observer/source_scalar/compiler.py',
    'model/observer/source_scalar/cutoff.py', 'model/observer/source_scalar/experiment.py',
    'model/observer/source_scalar/field.py', 'model/observer_reconstruction/__init__.py',
    'model/relative_observer_geometry/__init__.py'))))
MEASURE_DENSITY_UPPER = F(2)
ALL_CLOCK_COVARIANCE_BOUND = uniform.GREEN_BOUND + F(5, 12)


def _integer(value):
    if type(value) is not int or value < 0:
        raise ValueError('an original nonnegative integer is required')
    return value


def harmonic_tail(max_degree):
    """Uniform expected squared H^-2 tail, for all three normalized fields.

    Real spherical harmonics are orthonormal for probability sphere area.
    The degree-zero coefficient vanishes by native load conservation.
    """
    degree = _integer(max_degree)
    return MEASURE_DENSITY_UPPER*ALL_CLOCK_COVARIANCE_BOUND/(degree+1)**2


def tail_probability(max_degree, radius):
    """Markov bound for the H^-2 tail norm exceeding an exact radius."""
    if (type(radius) is not F or type(radius.numerator) is not int
            or type(radius.denominator) is not int or radius.numerator <= 0
            or radius.denominator <= 0
            or math.gcd(radius.numerator, radius.denominator) != 1):
        raise ValueError('a positive exact Fraction radius is required')
    return min(F(1), harmonic_tail(max_degree)/radius**2)


def projection_error(level):
    """Sum of three mean-square field errors for affine coordinate readouts.

    Coordinates live on the auxiliary flat polyhedron and are transported
    radially to the sphere. They do not assert physical spatial coordinates.
    """
    level = _integer(level)
    return ALL_CLOCK_COVARIANCE_BOUND*F(2, 21*4**level)


def geometric_certificate():
    geometry = uniform.exact_base_geometry()
    a2 = sp.sympify(geometry['squared_unit_side'])
    r2 = sp.sympify(geometry['squared_inradius'])
    if not (a2 > sp.Rational(16, 15) and a2 < sp.Rational(8, 7)
            and r2 > sp.Rational(9, 16)
            and sp.Rational(5, 3)**2 < 3 < sp.Rational(7, 4)**2):
        raise ArithmeticError('comparison-area rational margins failed')
    # a0=sqrt(3)*a²/4 lies in (4/9,1/2); r>3/4.
    density_upper = F(22, 7)/(5*F(4, 9)*F(3, 4))
    density_lower = 3*F(9, 16)/(5*F(1, 2))
    if not F(1, 2) < density_lower < density_upper < MEASURE_DENSITY_UPPER:
        raise ArithmeticError('radial probability-area comparison failed')
    return {'sphere_area': 'dOmega/(4*pi)',
        'source_measure': 'radial pushforward of normalized flat-polyhedron area',
        'density_formula': 'pi*|x|^3/(5*a0*r)',
        'density_lower_strict': str(density_lower),
        'density_upper_strict': str(density_upper),
        'density_enclosure': ['1/2', '2'],
        'affine_vector_projection_squared_L2_error': 'a^2/(12*4^level)',
        'squared_unit_side': str(a2)}


def affine_moments(level):
    """Exact triangle integration, independent of a sampled covariance fit."""
    level = _integer(level)
    if level > 3:
        raise ValueError('this finite geometry control is limited to levels zero through three')
    vertices, faces = uniform.affine_prefix(level)
    root = sp.sqrt(5)
    radius2 = (5+root)/2
    vertices = [sp.Matrix([(v[i]+v[i+1]*root)/2**(level+1)
                           for i in (0, 2, 4)]) for v in vertices]
    projected = sp.zeros(3)
    integrated = sp.zeros(3)
    first = sp.zeros(3, 1)
    for face in faces:
        corners = [vertices[i] for i in face]
        total = sum(corners, sp.zeros(3, 1))
        center = total/3
        first += center
        projected += center*center.T
        # Uniform barycentric second moments: E(t_i²)=1/6, E(t_i t_j)=1/12.
        integrated += (sum((v*v.T for v in corners), sp.zeros(3))+total*total.T)/12
    projected = sp.simplify(projected/(len(faces)*radius2))
    integrated = sp.simplify(integrated/(len(faces)*radius2))
    residual = sp.simplify(integrated-projected)
    a2 = 4/radius2
    expected = sp.eye(3)*a2/(36*4**level)
    if sp.simplify(first) != sp.zeros(3, 1) or sp.simplify(residual-expected) != sp.zeros(3):
        raise ArithmeticError('exact affine projection or conserved mean failed')
    if sp.simplify(integrated-sp.eye(3)*(1+3*(5+2*root)/15)/12) != sp.zeros(3):
        raise ArithmeticError('triangle second moments changed with refinement')
    return {'level': level, 'cells': len(faces),
        'projected_L2_Gram_diagonal': str(projected[0, 0]),
        'exact_vector_projection_error': str(sp.simplify(sp.trace(residual))),
        'mean_zero_and_exact_triangle_integration_checked': True}


def native_affine_covariance(level):
    """Three continuous affine readouts on the full native port Laplacian."""
    level = _integer(level)
    if level > 2:
        raise ValueError('the retained native covariance control uses levels zero through two')
    vertices, faces = uniform.affine_prefix(level)
    root = np.sqrt(5.)
    radius = np.sqrt((5+root)/2)
    vertices = np.asarray([[(v[i]+v[i+1]*root)/(2**(level+1)*radius)
                            for i in (0, 2, 4)] for v in vertices])
    means = np.asarray([vertices[list(face)].mean(axis=0) for face in faces])
    z = (means-means.mean(axis=0))/np.sqrt(len(faces))
    _, laplacian, _, _, _, _ = uniform.ancestry.problem(level)
    w = np.repeat(z/np.sqrt(12), 12, axis=0)
    solved = np.zeros_like(w)
    solved[1:] = splu(laplacian[1:, 1:].astype(float)).solve(w[1:])
    solved -= solved.mean(axis=0)
    residual = float(np.linalg.norm(laplacian@solved-w))
    covariance = w.T@solved/4**level
    if residual > 2e-10 or np.linalg.eigvalsh(covariance)[0] <= 0:
        raise ArithmeticError('native affine-field Poisson control failed')
    if np.linalg.norm(covariance, 2) > float(uniform.GREEN_BOUND)*np.linalg.norm(z, 2)**2:
        raise ArithmeticError('native covariance exceeds the all-level operator bound')
    return {'level': level, 'covariance_diagnostic': covariance.tolist(),
        'Poisson_residual_diagnostic': residual,
        'native_ports': len(w), 'continuous_affine_smearings': 3}


def run():
    return {'schema': 'oph.native_continuum_distribution.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'geometry': geometric_certificate(),
        'exact_projection_controls': [affine_moments(level) for level in range(3)],
        'native_covariance_controls': [native_affine_covariance(level) for level in range(3)],
        'all_level_enclosures': {
            'all_clock_L2_covariance_norm_upper': str(ALL_CLOCK_COVARIANCE_BOUND),
            'expected_H_minus_2_norm_squared_upper': str(harmonic_tail(0)),
            'harmonic_tail_examples': {str(n): str(harmonic_tail(n)) for n in (0, 10, 100, 1000)},
            'affine_projection_examples': {str(n): str(projection_error(n)) for n in (0, 3, 20)},
            'full_dictionary_clock_operator_difference': {
                str(n): str(uniform.bounds(n, n)['incident_Green_difference_upper'])
                for n in (0, 3, 20)}},
        'theorem': {
            'subsequential_weak_limit_in_H_minus_2_of_S2': True,
            'Green_and_two_normalized_clock_fields_share_limit': True,
            'bounded_L2_test_function_covariance_extension': True,
            'nonzero_three_coarse_smearings_retained_in_mean_square': True,
            'limit_order': 'finite-graph observation-count CLT, then common spatial subsequence',
            'physical_spherical_metric_selected': False,
            'pointwise_or_L2_random_field_claimed': False,
            'unique_covariance_or_CFT_identified': False,
            'Lorentzian_dynamics_or_modular_time_identified': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native continuum-distribution bounds and source controls verified.')
    else:
        print(canonical(result))


if __name__ == '__main__':
    main()
