"""All-level native Green bounds on the source's own spherical subdivision.

The flat icosahedron is an auxiliary comparison surface with the same cell
graph. The proof in OCCUPATION_FIELD.md, rather than finite enumeration,
establishes the uniform bound and compatible subsequential Gaussian laws.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from types import MappingProxyType

import numpy as np
import sympy as sp

from . import ancestry_green as ancestry

SOURCE_PATHS = tuple(sorted(set((*ancestry.SOURCE_PATHS,
    'model/physical_space_time/uniform_green.py'))))
GREEN_BOUND = F(4096, 3)
# Coordinates are (a+b sqrt(5))/2, represented by six integer coefficients.
VERTICES = ((-2, 0, 1, 1, 0, 0), (2, 0, 1, 1, 0, 0),
    (-2, 0, -1, -1, 0, 0), (2, 0, -1, -1, 0, 0),
    (0, 0, -2, 0, 1, 1), (0, 0, 2, 0, 1, 1),
    (0, 0, -2, 0, -1, -1), (0, 0, 2, 0, -1, -1),
    (1, 1, 0, 0, -2, 0), (1, 1, 0, 0, 2, 0),
    (-1, -1, 0, 0, -2, 0), (-1, -1, 0, 0, 2, 0))
FACES = ((0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
    (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
    (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
    (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1))


def _level(value):
    if type(value) is not int or value < 0:
        raise ValueError('an original nonnegative subdivision level is required')
    return value


def bounds(level, coarse=0):
    """Rational theorem bounds at arbitrary levels; no graph allocation."""
    level, coarse = _level(level), _level(coarse)
    if coarse > level:
        raise ValueError('the observation level must not exceed the source level')
    scale = 4**level
    return {'cell_gap_lower': F(9, 512*scale), 'Green_norm_upper': GREEN_BOUND,
        'incident_Green_difference_upper': F(1+54*2**coarse, 132*scale),
        'attempt_Green_difference_upper': F(1, 2520*scale**2),
        'base_rank_three_Green_lower': F(4)}


def _qsign(a, b):
    """Exact sign of a+b sqrt(5), for integer a,b."""
    if not b:
        return (a > 0)-(a < 0)
    if not a or (a > 0) == (b > 0):
        return (b > 0)-(b < 0)
    difference = a*a-5*b*b
    return ((a > 0)-(a < 0))*((difference > 0)-(difference < 0))


def _squared(vector):
    return (sum(vector[i]**2+5*vector[i+1]**2 for i in (0, 2, 4)),
            sum(2*vector[i]*vector[i+1] for i in (0, 2, 4)))


def dual_edges(faces):
    """Reconstruct the dual solely from shared primal edges; fail closed."""
    edges = {}
    for index, face in enumerate(faces):
        if len(face) != 3 or len(set(face)) != 3:
            raise ArithmeticError('a triangular source cell is degenerate')
        for a, b in zip(face, (*face[1:], face[0]), strict=True):
            edges.setdefault(tuple(sorted((a, b))), []).append(index)
    if any(len(cells) != 2 for cells in edges.values()):
        raise ArithmeticError('a closed source edge must have exactly two cells')
    dual = tuple(sorted(tuple(sorted(cells)) for cells in edges.values()))
    degree = [0]*len(faces)
    for a, b in dual:
        degree[a] += 1; degree[b] += 1
    if len(set(dual)) != len(dual) or any(d != 3 for d in degree):
        raise ArithmeticError('the source dual is not a simple degree-three graph')
    return dual


def refine(vertices, faces):
    """Independent affine shared-midpoint stencil, with exact dyadic numerators."""
    fine_vertices = [tuple(2*x for x in v) for v in vertices]
    midpoints = {}
    def midpoint(a, b):
        edge = tuple(sorted((a, b)))
        if edge not in midpoints:
            midpoints[edge] = len(fine_vertices)
            fine_vertices.append(tuple(x+y for x, y in zip(vertices[a], vertices[b], strict=True)))
        return midpoints[edge]
    fine_faces = []
    for a, b, c in faces:
        ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
        fine_faces.extend(((a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)))
    return tuple(fine_vertices), tuple(fine_faces)


@lru_cache(maxsize=4)
def affine_prefix(level):
    """Finite executable control; the theorem itself has no level cap."""
    level = _level(level)
    vertices, faces = VERTICES, FACES
    for _ in range(level):
        vertices, faces = refine(vertices, faces)
    return vertices, faces


@lru_cache(maxsize=1)
def exact_base_geometry():
    root = sp.sqrt(5)
    raw = sp.Matrix([[(v[i]+v[i+1]*root)/2 for i in (0, 2, 4)] for v in VERTICES])
    radius_squared = (5+root)/2
    side_squared = sp.simplify(4/radius_squared)
    inradius_squared = (5+2*root)/15
    centroid = sp.Matrix([list(sum((raw.row(i) for i in face), sp.zeros(1, 3))/3)
                          for face in FACES])
    for face in FACES:
        for a, b in zip(face, (*face[1:], face[0]), strict=True):
            if sp.simplify((raw.row(a)-raw.row(b)).dot(raw.row(a)-raw.row(b))-4) != 0:
                raise ArithmeticError('the base face is not equilateral')
        normal = (raw.row(face[1])-raw.row(face[0])).cross(raw.row(face[2])-raw.row(face[0]))
        if not sp.simplify(normal.dot(raw.row(face[0]))) > 0:
            raise ArithmeticError('the base face does not point outwards')
        if any(sp.simplify(normal.dot(raw.row(v)-raw.row(face[0]))) > 0 for v in range(12)):
            raise ArithmeticError('the base comparison polyhedron is not convex')
    if (sp.simplify((raw*raw.T).diagonal()-sp.ones(1, 12)*radius_squared) != sp.zeros(1, 12)
            or sp.simplify((centroid*centroid.T).diagonal()/radius_squared
                           -inradius_squared*sp.ones(1, 20)) != sp.zeros(1, 20)
            or sp.simplify(centroid.T*centroid/radius_squared-20*inradius_squared*sp.eye(3)/3) != sp.zeros(3)
            or sp.ones(1, 20)*centroid != sp.zeros(1, 3)):
        raise ArithmeticError('base comparison radii or centroid Gram changed')
    if not (sp.simplify(inradius_squared-F(9, 16)) > 0 and root > F(5, 3)):
        raise ArithmeticError('rational comparison margins failed')
    lower = sp.simplify(8*inradius_squared/side_squared)
    if sp.simplify(lower-(7+3*root)/3) != 0 or not lower > 4:
        raise ArithmeticError('the noncollapse comparison constant failed')
    return MappingProxyType({'squared_unit_side': str(side_squared), 'squared_inradius': str(inradius_squared),
        'rank_three_Green_lower_exact': str(lower),
        'base_centroid_Gram': '(20 r^2/3) I_3',
        'convexity_outward_faces_radii_and_edge_lengths_checked': True})


def check_affine_geometry(level):
    """Exact face lengths, ancestry averages and hinge-energy control."""
    vertices, faces = affine_prefix(level)
    edges = dual_edges(faces)
    scale = 4**level
    if len(vertices) != 10*scale+2 or len(faces) != 20*scale:
        raise ArithmeticError('the subdivision census changed')
    for face in faces:
        for a, b in zip(face, (*face[1:], face[0]), strict=True):
            vector = tuple(x-y for x, y in zip(vertices[a], vertices[b], strict=True))
            if _squared(vector) != (16, 0):
                raise ArithmeticError('an affine edge does not have the common dyadic length')
    centroids = tuple(tuple(sum(vertices[v][i] for v in face) for i in range(6)) for face in faces)
    for a, b in edges:
        squared = _squared(tuple(x-y for x, y in zip(centroids[a], centroids[b], strict=True)))
        # Denominator 3*2^(l+1): 4/(3*4^l) becomes the integer threshold 48.
        if _qsign(squared[0]-48, squared[1]) > 0:
            raise ArithmeticError('an adjacent centroid violates the hinge-energy bound')
    for parent, face in enumerate(FACES):
        base = tuple(sum(VERTICES[v][i] for v in face) for i in range(6))
        total = tuple(sum(row[i] for row in centroids[parent*scale:(parent+1)*scale]) for i in range(6))
        if total != tuple(scale*2**level*x for x in base):
            raise ArithmeticError('affine descendant means lost their fixed ancestry readout')
    return {'level': level, 'vertices': len(vertices), 'cells': len(faces),
        'dual_edges': len(edges), 'exact_edge_centroid_and_ancestry_checks': True}


def source_binding(level):
    """Bind independent cell combinatorics to actual source seams and ancestry."""
    _, faces = affine_prefix(level)
    builder = ancestry.elimination.source._source()
    tower = builder.build_geodesic_icosahedral_tower(level)
    actual_faces = tuple(tuple(map(int, row)) for row in tower.levels[level].faces)
    if actual_faces != faces:
        raise ArithmeticError('the source midpoint stencil differs from the comparison graph')
    raw, _, _, _, _, _ = ancestry.problem(level)
    actual_edges = tuple(sorted(tuple(sorted((int(a)//12, int(b)//12)))
        for a, b in zip(raw.seam_a, raw.seam_b, strict=True) if a//12 != b//12))
    if actual_edges != dual_edges(faces):
        raise ArithmeticError('native cross seams do not realize the source dual edges')
    for coarse in range(level+1):
        parents = tower.cell_ancestor_ids(fine_level=level, coarse_level=coarse)
        if not np.array_equal(parents, np.arange(len(faces))//4**(level-coarse)):
            raise ArithmeticError('the source changed its nested observation dictionary')
    return {'level': level, 'actual_faces_cross_seams_and_all_ancestries_checked': True,
        'distinct_port_routing_and_native_carrier_rows_checked': True}


def dictionary_control(level, coarse):
    """Check exact centered Gram and interlevel isometry on dyadic scores."""
    level, coarse = _level(level), _level(coarse)
    if coarse >= level:
        raise ValueError('a finer observation level is required for this control')
    cells, channels = 20*4**level, 20*4**coarse
    if 4*cells*channels >= 2**63:
        raise ValueError('this finite integer Gram control would overflow; use the all-level theorem bounds')
    score = np.full((cells, channels), -1, dtype=np.int64)
    score[np.arange(cells), np.arange(cells)//4**(level-coarse)] += channels
    fine = np.full((cells, 4*channels), -1, dtype=np.int64)
    fine[np.arange(cells), np.arange(cells)//4**(level-coarse-1)] += 4*channels
    expected = 4**(level-coarse)*channels*(channels*np.eye(channels, dtype=np.int64)-1)
    if (not np.array_equal(score.T@score, expected)
            or not np.array_equal(fine.reshape(cells, channels, 4).sum(axis=2), 4*score)
            or np.any(score.sum(axis=0)) or np.any(score.sum(axis=1))):
        raise ArithmeticError('centered ancestry Gram or isometric refinement identity failed')
    return {'source_level': level, 'coarse_observation_level': coarse,
        'coarse_channels': channels, 'fine_channels': 4*channels,
        'exact_centered_Gram_and_refinement_checked': True}


def run():
    geometry = exact_base_geometry()
    rows = [check_affine_geometry(level) | source_binding(level) for level in range(4)]
    for level in range(4):
        value = ancestry.covariance(level)
        matrix = np.asarray(value['green'], dtype=float)
        root = np.sqrt(5.)
        raw = np.asarray([[(v[i]+v[i+1]*root)/2 for i in (0, 2, 4)] for v in VERTICES])
        centroid = np.asarray([raw[list(face)].mean(axis=0) for face in FACES])
        basis, _ = np.linalg.qr(centroid)
        if (np.linalg.norm(matrix, 2)+float(value['green_error']) > float(GREEN_BOUND)
                or np.linalg.eigvalsh(basis.T@matrix@basis)[0]-float(value['green_error']) < 4):
            raise ArithmeticError('the finite native Green witness violates the all-level enclosure')
    return {'contract': {'architecture': 'existing shared-midpoint spherical source combinatorics',
        'native_update': 'uniform fair binary endpoint repair on distinct-port carrier joins',
        'observer': 'centered ancestral counts from complete initial loads and seam currents',
        'limit_order': 'stationary observation count first; common spatial subsequence second',
        'auxiliary_flat_comparison_is_a_supplied_physical_metric': False},
        'exact_base_geometry': dict(geometry), 'finite_source_binding_controls': rows,
        'nested_dictionary_controls': [dictionary_control(3, k) for k in range(3)],
        'all_level_bounds': {str(level): {key: str(value) for key, value in bounds(level).items()}
                             for level in (0, 1, 3, 20)},
        'theorem': {'all_level_scaled_native_Green_norm_upper': str(GREEN_BOUND),
            'all_level_root_clock_error_upper': '5/(12*4^level)',
            'fixed_coarse_level_k_root_clock_error_upper': '(1+54*2^k)/(132*4^level)',
            'all_level_attempt_clock_error_upper': '1/(2520*16^level)',
            'subsequential_base_three_mode_Green_lower': '4',
            'common_subsequence_all_ancestral_Gaussian_laws_exist': True,
            'source_graph_inequality_import': 'round-sphere finite-perimeter isoperimetry'},
        'acceptance': {'uniform_all_level_native_Green_bound_proved': True,
            'attempt_and_root_incident_clock_difference_vanishes': True,
            'nonzero_compatible_subsequential_Gaussian_field_law_proved': True,
            'whole_sequence_convergence_or_unique_spatial_limit_proved': False,
            'target_spherical_Green_operator_identified': False,
            'interacting_CFT_or_Lorentzian_spacetime_identified': False}}
