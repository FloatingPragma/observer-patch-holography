"""Exact scalar-potential constraints on the native neutral family triplet.

This concerns one real vector with an invariant scalar potential and no
additional spurions. It does not identify scalar curvatures with fermion
Yukawas, select a physical vacuum, or supply a physical clock or mass map.
"""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
import hashlib

import sympy as s

from . import native_gauge_object as objects, native_family_selection as selector
from . import slow_geometry

ROOT = objects.ROOT
RADICAL = s.sqrt(5)


def _scalar(value):
    return s.simplify(s.expand(value))


def _zero(matrix):
    return all(s.expand(x) == 0 or _scalar(x) == 0 for x in matrix)


def _poly(expression, *variables):
    return s.Poly(expression, *variables, extension=RADICAL)


@lru_cache(maxsize=1)
def architecture():
    source = objects.geometry(); B = source['frame']
    metric = s.expand((B.T*B)[0, 0]); norm2 = s.expand((B[0, :]*B[0, :].T)[0])
    rotations, traces = [], Counter()
    for g in source['rotations']:
        R = objects._clean(B.T*objects.permutation_matrix(g)*B/metric)
        if not _zero(R.T*R-s.eye(3)) or _scalar(R.det()-1) != 0:
            raise ArithmeticError('actual family action is not a proper real rotation')
        rotations.append(R); traces[_scalar(s.trace(R))] += 1
    # Symmetric-power characters follow from det(I-tR)=1-tr(R)t+tr(R)t²-t³.
    dimensions = []
    for degree in range(7):
        total = 0
        for trace, count in traces.items():
            h = [s.Integer(1)]
            for k in range(1, degree+1):
                h.append(s.expand(trace*h[k-1]-(trace*h[k-2] if k >= 2 else 0)
                                  +(h[k-3] if k >= 3 else 0)))
            total += count*h[degree]
        dimensions.append(_scalar(total/len(rotations)))
    if dimensions != [1, 0, 1, 0, 1, 0, 2]:
        raise ArithmeticError('native low-degree invariant dimensions changed')
    x = s.Matrix(s.symbols('x y z', real=True)); r2 = (x.T*x)[0]
    moments = {d: _poly(sum((B[i, :]*x)[0]**d for i in range(12))/norm2**(d//2), *x).as_expr()
               for d in (2, 4, 6)}
    if _scalar(moments[2]-4*r2) != 0 or _scalar(moments[4]-s.Rational(12, 5)*r2**2) != 0:
        raise ArithmeticError('actual port second/fourth moments are not radial')
    f6 = moments[6]; gradient = s.Matrix([s.diff(f6, q) for q in x]); hessian = s.hessian(f6, x)
    gram = sum(((B[i, :]*x)[0]**4*B[i, :].T*B[i, :]/norm2**3
                for i in range(12)), s.zeros(3))
    gram = gram.applyfunc(lambda q: _poly(q, *x).as_expr())
    if not _zero(hessian-30*gram):
        raise ArithmeticError('the native positive Gram is not the sixth-moment derivative')
    u, v, w = s.symbols('u v w', real=True)
    cubic = 0
    for powers, coefficient in _poly(f6, *x).terms():
        if any(k % 2 for k in powers): raise ArithmeticError('native moment lost coordinate sign symmetry')
        cubic += coefficient*u**(powers[0]//2)*v**(powers[1]//2)*w**(powers[2]//2)
    if _scalar(cubic-cubic.xreplace({u: v, v: w, w: u})) != 0:
        raise ArithmeticError('native squared-coordinate moment is not cyclic')
    restricted = s.expand(cubic.subs(w, 1-u-v))
    du, dv = s.diff(restricted, u), s.diff(restricted, v)
    resultant = _poly(s.resultant(du, dv, v), u).monic()
    candidates = (s.Rational(1, 3), s.Rational(1, 4), (3+RADICAL)/8, (3-RADICAL)/8)
    if resultant != _poly(s.prod(u-a for a in candidates), u):
        raise ArithmeticError('complete interior stationary resultant changed')
    interior = []
    for a in candidates:
        gcd = s.gcd(_poly(du.subs(u, a), v), _poly(dv.subs(u, a), v)).monic()
        if gcd.degree() != 1: raise ArithmeticError('a stationary coordinate has a missing or multiple lift')
        b = _scalar(-gcd.nth(0)); point = (a, b, _scalar(1-a-b))
        if any(q.is_positive is not True for q in point): raise ArithmeticError('an interior point left the sphere octant')
        interior.append(point)
    edge_derivative = _poly(s.diff(cubic.subs({v: 1-u, w: 0}), u), u).monic()
    boundary_roots = ((5-RADICAL)/10, (3+RADICAL)/6)
    if edge_derivative != _poly(s.prod(u-a for a in boundary_roots), u):
        raise ArithmeticError('complete boundary stationary polynomial changed')
    if any(q.is_positive is not True or (1-q).is_positive is not True for q in boundary_roots):
        raise ArithmeticError('a boundary stationary point is outside its edge')
    critical = Counter()
    for point in interior:
        critical[_scalar(cubic.subs(dict(zip((u, v, w), point))))] += 8
    for a in boundary_roots:
        critical[_scalar(cubic.subs({u: a, v: 1-a, w: 0}))] += 12
    critical[_scalar(cubic.subs({u: 1, v: 0, w: 0}))] += 6
    # Check the representatives and all their actual native orbits. Completeness
    # follows from the resultant, its unique lifts and all boundary strata.
    face = slow_geometry.source().oriented_faces()[0]; edge = source['seams'][0]
    seeds = {'vertex': B[0, :].T,
             'face': sum((B[i, :].T for i in face), s.zeros(3, 1)),
             'edge': sum((B[i, :].T for i in edge), s.zeros(3, 1))}
    classes = {}
    for name, seed in seeds.items():
        squared = _scalar((seed.T*seed)[0]); substitution = dict(zip(x, seed))
        value = _scalar(f6.subs(substitution)/squared**3)
        if not _zero(gradient.subs(substitution)-6*f6.subs(substitution)*seed/squared):
            raise ArithmeticError('a claimed native stationary orbit is not stationary')
        Q = s.eye(3)-seed*seed.T/squared
        tangent = (Q*(hessian.subs(substitution)/squared**2-6*value*s.eye(3))*Q).applyfunc(_scalar)
        lam = s.Symbol('lambda'); polynomial = _poly(tangent.charpoly(lam).as_expr(), lam).as_expr()
        orbit = {tuple(_scalar(q) for q in R*seed) for R in rotations}
        classes[name] = {'size': len(orbit), 'value': value, 'tangent': tangent,
                         'tangent_characteristic_polynomial': s.factor(polynomial, extension=RADICAL)}
    expected = {'vertex': (12, s.Rational(52, 25), (-s.Rational(192, 25),)*2),
                'face': (20, s.Rational(68, 45), (s.Rational(64, 15),)*2),
                'edge': (30, s.Rational(8, 5), (12*(1-RADICAL)/5, 12*(1+RADICAL)/5))}
    for name, (count, value, curvatures) in expected.items():
        row = classes[name]; lam = s.Symbol('lambda')
        if (row['size'] != count or row['value'] != value or critical[value] != count
                or _poly(row['tangent_characteristic_polynomial'], lam)
                != _poly(lam*s.prod(lam-c for c in curvatures), lam)):
            raise ArithmeticError('complete native orbit or tangent-curvature classification failed')
        row['tangent_eigenvalues'] = curvatures
        row['source_Gram_eigenvalues'] = tuple(_scalar((q+6*value)/30) for q in curvatures)+(value,)
    if sum(critical.values()) != 62 or set(critical) != {row[1] for row in expected.values()}:
        raise ArithmeticError('a stationary sphere stratum is missing')
    # Literal attachment to the selected neutral native object, not merely an
    # abstract named three-dimensional representation.
    selected = selector.geometry(); F = (selected['frame']*B/s.sqrt(metric)).applyfunc(_scalar)
    if not _zero(F.T*F-s.eye(3)) or not _zero(F*F.T-selected['P']):
        raise ArithmeticError('scalar family coordinates do not span the selected persistent object')
    probe = s.Matrix([1, 2, 3]); test_gram = gram.subs(dict(zip(x, probe)))
    discriminant = _scalar(s.discriminant(test_gram.charpoly().as_expr(), test_gram.charpoly().gen))
    if (_zero(test_gram*probe*probe.T-probe*probe.T*test_gram)
            or discriminant.is_positive is not True or _scalar(test_gram.det()).is_positive is not True):
        raise ArithmeticError('the nonstationary diagnostic failed to leave the degenerate Gram class')
    return {'frame': B, 'rotations': tuple(rotations), 'trace_counts': traces,
            'dimensions': tuple(dimensions), 'variables': x, 'moments': moments,
            'source_Gram': gram, 'nonstationary_Gram_discriminant': discriminant,
            'squared_coordinate_variables': (u, v, w), 'squared_coordinate_cubic': cubic,
            'interior_resultant_monic': resultant.as_expr(), 'interior_points': tuple(interior),
            'boundary_polynomial_monic': edge_derivative.as_expr(), 'boundary_roots': boundary_roots,
            'classes': classes, 'native_family_isometry': F}


def report():
    a = architecture()
    return {'subject': 'one real vector in the selected neutral native family triplet',
        'hypotheses': ['real polynomial or Taylor potential invariant under all60 native proper rotations',
                       'no additional spurion, oriented history or matrix-valued order parameter',
                       'uniform twelve-port moment defines the tested degree-six anisotropy'],
        'same_selected_native_family_object_isometry_exact': True,
        'actual_rotation_trace_counts': [[str(k), v] for k, v in sorted(a['trace_counts'].items(), key=lambda row: float(row[0]))],
        'homogeneous_invariant_dimensions_degrees_zero_through_six': list(map(int, a['dimensions'])),
        'complete_scalar_potential_through_degree_four': 'c0+c2*|phi|^2+c4*|phi|^4',
        'native_port_second_moment': '4*|phi|^2', 'native_port_fourth_moment': '(12/5)*|phi|^4',
        'first_anisotropic_degree': 6, 'sixth_moment_exact': str(a['moments'][6]),
        'complete_sphere_stationary_certificate': {
            'variables': 'u=x^2,v=y^2,w=z^2; u+v+w=1',
            'interior_resultant_monic': str(a['interior_resultant_monic']),
            'each_interior_resultant_root_has_one_positive_lift': True,
            'interior_squared_coordinates': [[str(q) for q in point] for point in a['interior_points']],
            'boundary_stationary_polynomial_monic': str(a['boundary_polynomial_monic']),
            'boundary_roots': list(map(str, a['boundary_roots'])),
            'coordinate_axis_endpoints_included': True, 'total_stationary_points': 62,
            'orbits': [{'native_orbit': name, 'point_count': row['size'],
                       'I6_value_exact': str(row['value']),
                       'tangent_curvatures_exact': list(map(str, row['tangent_eigenvalues'])),
                       'source_Gram_eigenvalues_exact': list(map(str, row['source_Gram_eigenvalues']))}
                      for name, row in a['classes'].items()]},
        'fixed_radius_selection': {'positive_I6_coefficient': 'twenty face directions; C3 stabilizer',
            'negative_I6_coefficient': 'twelve vertex directions; C5 stabilizer',
            'radial_stability_or_nonzero_vacuum_selected': False},
        'coupling_boundary': {
            'positive_native_family_Gram': 'K6(phi)=sum_i(n_i.dot(phi))^4*n_i*n_i^T=Hess(I6)/30',
            'real_symmetric_Gram_at_both_selected_orbits_has_a_doublet': True,
            'nonstationary_probe_1_2_3_has_three_distinct_Gram_eigenvalues': True,
            'nonstationary_probe_is_a_selected_state_or_mass_candidate': False,
            'tangent_scalar_curvature_is_fermion_Yukawa_spectrum': False,
            'real_symmetric_axis_stabilizer_commutant': 'a*I+b*u*u^T; transverse pair degenerate',
            'general_complex_Hermitian_counterexample': '3*I+2*u*u^T+i*[u]_cross has eigenvalues2,4,5 for unit u',
            'C3_or_C5_alone_forces_general_complex_Yukawa_degeneracy': False,
            'physical_kinetic_normalization_or_mass_map_derived': False},
        'potential_evolution_or_vacuum_preparation_executed': False,
        'single_vector_order_parameter_is_a_selected_physical_field': False,
        'empirical_inputs': []}


def source_bindings():
    names = set(objects.source_bindings()) | set(selector.source_bindings()) | {
        'model/common_source/native_family_potential.py',
        'model/common_source/test_native_family_potential.py'}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}
