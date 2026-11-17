"""Compact, strongly local native spatial limits with spectral dimension two.

The source topology and its ancestry dictionary are the existing ones.  The
result upgrades their common weak Green subsequence to operator-norm
convergence, with an injective compact limit and a positive inverse generator.
Full-port energy lifts supply its conservative strongly local Dirichlet form;
an explicit source interpolation identifies its regular H1 form core, and
triangle tests bound both sides of its spectral growth. Exact terminal-star
elimination gives a local conductance graph with the same continuum Green
limits, an all-level Loewner sandwich against its harmonic-mean isotropic
reference, and a radial Poisson-Steklov half-cylinder bulk extension whose
boundary trace covariance is related to the native Green limit by an A^(1/4)
spectral filter. A selected intrinsic metric and physical spatial
identification are not asserted.
The all-level proof is in OCCUPATION_FIELD.md.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import shortest_path
import sympy as sp

from . import continuum_field as continuum
from . import source_green as source
from . import uniform_green as uniform
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*continuum.SOURCE_PATHS,
    'model/physical_space_time/local_continuum.py'))))
GREEN_BOUND = uniform.GREEN_BOUND
PATCH_POINCARE_CONSTANT = F(26)
COARSE_ENERGY_TAIL_CONSTANT = F(624)


def _int(value, minimum=0, maximum=None):
    if (type(value) is not int or value < minimum
            or (maximum is not None and value > maximum)):
        raise ValueError('an original integer in the certified range is required')
    return value


def _level(value):
    return _int(value, minimum=0)


def compactness_tail_squared(coarse_level):
    """Squared norm bound for (I-P_k)G_l, at every l>=k."""
    return COARSE_ENERGY_TAIL_CONSTANT*GREEN_BOUND/4**_level(coarse_level)


def covariance_eigenvalue_tail(coarse_level):
    """Bound on positive eigenvalue number 20*4^k and every later one."""
    return COARSE_ENERGY_TAIL_CONSTANT/4**_level(coarse_level)


def continuum_eigenvalue_enclosure(index):
    """All positive eigenvalues, counted from one, for every selected limit."""
    if type(index) is not int or index < 1:
        raise ValueError('an original positive integer eigenvalue index is required')
    return {'lower': F(index, 49920), 'upper': F(3*index, 2)}


def tent_eigenvalue_bound(coarse_level):
    """Min-max bound from mean-zero combinations of ancestral triangle tents."""
    scale = 4**_level(coarse_level)
    return {'index': 20*scale-1, 'upper': F(3*scale, 2)}


def harmonic_lift_tail(level):
    """Squared internal fluctuation per unit scaled full-port energy."""
    return F(1, 2*4**_level(level))


@lru_cache(maxsize=1)
def _carrier_resistances():
    """Exact resistances of the record-recovered twelve-port Laplacian."""
    _, laplacian = source.passage.recovered_carrier()
    green = (laplacian+sp.ones(12)/12).inv()-sp.ones(12)/12
    return tuple(tuple(F(green[a, a]+green[b, b]-2*green[a, b])
                       for b in range(12)) for a in range(12))


def terminal_star_arms(ports):
    """Three positive star resistances, ordered like the actual terminal ports."""
    if (type(ports) is not tuple or len(ports) != 3 or len(set(ports)) != 3
            or any(type(p) is not int or not 0 <= p < 12 for p in ports)):
        raise ValueError('three distinct original integer carrier ports are required')
    resistance = _carrier_resistances()
    arms = tuple((resistance[ports[i]][ports[(i+1) % 3]]
                  +resistance[ports[i]][ports[(i+2) % 3]]
                  -resistance[ports[(i+1) % 3]][ports[(i+2) % 3]])/2
                 for i in range(3))
    if not all(F(2, 15) <= arm <= F(3, 10) for arm in arms):
        raise ArithmeticError('the exact native terminal-star bounds failed')
    return arms


@lru_cache(maxsize=220)
def _terminal_harmonic_lift(ports):
    """Exact minimum-energy twelve-port lift of three terminal values."""
    terminal_star_arms(ports)
    _, laplacian = source.passage.recovered_carrier()
    interior = tuple(p for p in range(12) if p not in ports)
    interior_lift = -laplacian.extract(interior, interior).inv()*laplacian.extract(interior, ports)
    lift = sp.zeros(12, 3)
    for i, p in enumerate(ports):
        lift[p, i] = 1
    for i, p in enumerate(interior):
        for j in range(3):
            lift[p, j] = interior_lift[i, j]
    if any(value < 0 for value in lift) or lift*sp.ones(3, 1) != sp.ones(12, 1):
        raise ArithmeticError('native terminal harmonic lift lost its maximum principle')
    return tuple(tuple(F(lift[p, i]) for i in range(3)) for p in range(12))


def terminal_star_control():
    """Every carrier triple, with independent exact terminal-energy elimination."""
    _, laplacian = source.passage.recovered_carrier()
    distance = shortest_path((np.asarray(laplacian, dtype=int) < 0).astype(int),
                             directed=False, unweighted=True)
    resistance = _carrier_resistances()
    by_distance = {str(d): sorted({resistance[a][b] for a in range(12)
                                  for b in range(a+1, 12) if distance[a, b] == d})
                   for d in (1, 2, 3)}
    classes = Counter()
    for ports in combinations(range(12), 3):
        arms = terminal_star_arms(ports)
        lift = sp.Matrix(_terminal_harmonic_lift(ports))
        conductance = sp.Matrix([1/r for r in arms])
        star_trace = sp.diag(*conductance)-conductance*conductance.T/sum(conductance)
        if lift.T*laplacian*lift != star_trace:
            raise ArithmeticError('star and native terminal Dirichlet energies differ')
        classes[tuple(sorted(arms))] += 1
    return {'terminal_triples_checked': sum(classes.values()),
        'resistance_by_native_hop_distance': {d: [str(r) for r in values]
                                             for d, values in by_distance.items()},
        'star_arm_classes': [{'resistances': [str(r) for r in arms], 'triples': count}
                             for arms, count in sorted(classes.items())],
        'all_exact_terminal_Dirichlet_forms_and_harmonic_lifts_checked': True,
        'star_arm_lower': '2/15', 'star_arm_upper': '3/10',
        'weighted_dual_conductance_lower': '5/8', 'weighted_dual_conductance_upper': '15/19'}


def weighted_dual_network(graph):
    """Positive local conductances from the actual three distinct routed ports."""
    ports = [[] for _ in range(graph.cells)]
    for a, b in graph.inter_seams:
        ports[a//12].append(a % 12); ports[b//12].append(b % 12)
    ports = tuple(tuple(sorted(row)) for row in ports)
    arms = tuple(dict(zip(row, terminal_star_arms(row), strict=True)) for row in ports)
    edges = tuple((a//12, b//12, 1/(1+arms[a//12][a % 12]+arms[b//12][b % 12]))
                  for a, b in graph.inter_seams)
    if (len({(a, b) for a, b, _ in edges}) != 3*graph.cells//2
            or any(a == b or not F(5, 8) <= w <= F(15, 19) for a, b, w in edges)):
        raise ArithmeticError('the native weighted dual graph or its conductances changed')
    return ports, arms, edges


def weighted_green_difference_bound(level):
    """Operator norm on the common normalized, centered cell observation space."""
    return F(182, 2**_level(level))


def weighted_dual_control(level):
    """Exact energy-map controls and finite Green diagnostics on actual joins."""
    if _level(level) > 2:
        raise ValueError('weighted dual source controls use levels zero through two')
    graph, _, incidence, integer_f, denominator = source.matrices(level)
    ports, arms, edges = weighted_dual_network(graph)
    values = [F(17*p % 23-11, 5) for p in range(graph.ports)]
    centers, means = [], []
    for c, row in enumerate(ports):
        weights = [1/arms[c][p] for p in row]
        centers.append(sum(w*values[12*c+p] for w, p in zip(weights, row, strict=True))/sum(weights))
        means.append(sum(values[12*c:12*c+12])/12)
    native = sum((values[a]-values[b])**2 for a, b in graph.seams)
    reduced = sum(w*(centers[a]-centers[b])**2 for a, b, w in edges)
    if reduced > native or sum((a-b)**2 for a, b in zip(centers, means, strict=True)) > native/2:
        raise ArithmeticError('native-to-star energy or readout bound failed')
    target = [F(7*c % 19-9, 3) for c in range(graph.cells)]
    boundary = [{} for _ in range(graph.cells)]
    for (a, b), (c, d, w) in zip(graph.inter_seams, edges, strict=True):
        current = w*(target[c]-target[d])
        boundary[c][a % 12] = target[c]-arms[c][a % 12]*current
        boundary[d][b % 12] = target[d]+arms[d][b % 12]*current
    lifted = []
    for c, row in enumerate(ports):
        matrix = _terminal_harmonic_lift(row)
        lifted.extend(sum(coeff*boundary[c][p] for coeff, p in zip(line, row, strict=True))
                      for line in matrix)
    target_energy = sum(w*(target[a]-target[b])**2 for a, b, w in edges)
    lifted_energy = sum((lifted[a]-lifted[b])**2 for a, b in graph.seams)
    error = sum((sum(lifted[12*c:12*c+12])/12-target[c])**2 for c in range(graph.cells))
    if lifted_energy > target_energy or error > F(3, 10)*target_energy:
        raise ArithmeticError('star-to-native energy or readout bound failed')
    weighted = np.zeros((graph.cells, graph.cells))
    for a, b, w in edges:
        weighted[a, a] += float(w); weighted[b, b] += float(w)
        weighted[a, b] -= float(w); weighted[b, a] -= float(w)
    schur = incidence.T@np.linalg.solve(integer_f/denominator, incidence)/12
    difference = np.linalg.norm((source._green(schur)-source._green(weighted/12))/4**level, 2)
    bound = weighted_green_difference_bound(level)
    if difference > float(bound):
        raise ArithmeticError('weighted Green diagnostic violates the analytic bound')
    return {'level': level, 'cells': graph.cells, 'ports': graph.ports,
        'actual_terminal_star_classes': {','.join(map(str, key)): count for key, count in
            sorted(Counter(tuple(sorted(row.values())) for row in arms).items())},
        'exact_weighted_dual_conductances': {str(w): count for w, count in
            sorted(Counter(w for _, _, w in edges).items())},
        'both_exact_energy_decreases_and_readout_bounds_checked': True,
        'scaled_Green_difference_norm_diagnostic': float(difference),
        'scaled_Green_difference_norm_upper': str(bound)}


def cycle_variance_control(valence):
    """Exact rational PSD certificate for variance <= edge energy on C5/C6."""
    if type(valence) is not int or valence not in (5, 6):
        raise ValueError('a source vertex has original integer valence five or six')
    laplacian = sp.zeros(valence)
    for i in range(valence):
        j = (i+1) % valence
        laplacian[i, i] += 1; laplacian[j, j] += 1
        laplacian[i, j] -= 1; laplacian[j, i] -= 1
    residual = laplacian-sp.eye(valence)+sp.ones(valence)/valence
    pivots = []
    for _ in range(valence):
        pivot = residual[0, 0]
        if pivot < 0 or (pivot == 0 and any(residual[0, j] != 0
                                            for j in range(residual.cols))):
            raise ArithmeticError('the exact vertex-cycle variance certificate failed')
        pivots.append(pivot)
        tail = residual[1:, 1:]
        if pivot:
            tail -= residual[1:, :1]*residual[:1, 1:]/pivot
        residual = tail
    return {'valence': valence, 'rational_nonnegative_pivots': [str(x) for x in pivots],
        'variance_is_bounded_by_cycle_edge_energy': True}


def interpolation_geometry_control():
    """Exact six-triangle energy and mass identities for the source H1 lift."""
    vertex, midpoint, centroid = sp.symbols('vertex midpoint centroid', real=True)
    root = sp.sqrt(3)
    # One subtriangle has coordinates (0,0), (1/2,0), (1/2,sqrt(3)/6).
    gradient = sp.Matrix([2*(midpoint-vertex), 2*root*(centroid-midpoint)])
    integrated = sp.expand(root/24*(gradient.T*gradient)[0])
    target = root/6*(vertex-midpoint)**2+root/2*(centroid-midpoint)**2
    if sp.simplify(integrated-target) != 0:
        raise ArithmeticError('the exact six-triangle gradient identity changed')
    x, y = sp.symbols('x y', nonnegative=True)
    error = (1-x-y)*(vertex-centroid)+x*(midpoint-centroid)
    # The reference simplex has area1/2; division gives its normalized mean.
    mean_square = 2*sp.integrate(sp.integrate(error**2, (y, 0, 1-x)), (x, 0, 1))
    expected = ((vertex-centroid)**2+(midpoint-centroid)**2
                +(vertex-centroid)*(midpoint-centroid))/6
    if sp.simplify(mean_square-expected) != 0:
        raise ArithmeticError('the exact affine interpolation mass identity changed')
    return {'unit_side_subtriangle_energy': 'sqrt(3)/6*(v-m)^2+sqrt(3)/2*(c-m)^2',
        'subtriangle_mean_squared_error': '((v-c)^2+(m-c)^2+(v-c)*(m-c))/6',
        'source_vertex_cycle_certificates': [cycle_variance_control(n) for n in (5, 6)],
        'normalized_gradient_times_side_squared_per_dual_energy_upper': '7/30',
        'normalized_mass_error_times_cell_count_per_dual_energy_upper': '5/18'}


def source_interpolation_control(level):
    """Bind the interpolation, vertex cycles and all edge sums to native joins."""
    if _level(level) > 3:
        raise ValueError('source interpolation controls use levels zero through three')
    vertices, faces = uniform.affine_prefix(level)
    edge_cells, vertex_cells = {}, {v: set() for v in range(len(vertices))}
    for cell, face in enumerate(faces):
        for v in face:
            vertex_cells[v].add(cell)
        for a, b in zip(face, (*face[1:], face[0]), strict=True):
            edge_cells.setdefault(tuple(sorted((a, b))), []).append(cell)
    if any(len(cells) != 2 for cells in edge_cells.values()):
        raise ArithmeticError('source interpolation requires two cells at each primal edge')
    edges = set(uniform.dual_edges(faces))
    if edges != {tuple(sorted(cells)) for cells in edge_cells.values()}:
        raise ArithmeticError('interpolation and native dual-edge dictionaries differ')
    valences = Counter(map(len, vertex_cells.values()))
    if set(valences)-{5, 6} or valences[5] != 12:
        raise ArithmeticError('source interpolation vertex valences changed')
    adjacency = {cell: set() for cell in range(len(faces))}
    for a, b in edges:
        adjacency[a].add(b); adjacency[b].add(a)
    for star in vertex_cells.values():
        if any(len(adjacency[cell] & star) != 2 for cell in star):
            raise ArithmeticError('a source vertex star is not a cycle')
        reached, pending = set(), {next(iter(star))}
        while pending:
            cell = pending.pop()
            reached.add(cell)
            pending.update((adjacency[cell] & star)-reached)
        if reached != star:
            raise ArithmeticError('a source vertex cycle is disconnected')
    values = [F(17*cell % 23-11, 5) for cell in range(len(faces))]
    vertex_means = {v: sum(values[c] for c in cells)/len(cells)
                    for v, cells in vertex_cells.items()}
    dual_energy = sum((values[a]-values[b])**2 for a, b in edges)
    vm, cm, vc, mass = F(0), F(0), F(0), F(0)
    for edge, cells in edge_cells.items():
        midpoint = sum(values[c] for c in cells)/2
        for cell in cells:
            for vertex in edge:
                dv, dm = vertex_means[vertex]-values[cell], midpoint-values[cell]
                vm += (vertex_means[vertex]-midpoint)**2
                cm += dm*dm
                vc += dv*dv
                mass += (dv*dv+dm*dm+dv*dm)/(36*len(faces))
    gradient = vm/30+cm/10  # Normalized integral times original side squared.
    if (cm != dual_energy or vm > 4*dual_energy or vc > 4*dual_energy
            or gradient > F(7, 30)*dual_energy
            or mass > F(5, 18*len(faces))*dual_energy):
        raise ArithmeticError('the native interpolation energy or mass bound failed')
    return {'level': level, 'cells': len(faces),
        'vertex_valence_counts': {str(n): count for n, count in sorted(valences.items())},
        'all_native_vertex_stars_are_the_certified_cycles': True,
        'exact_dual_energy': str(dual_energy),
        'normalized_gradient_times_side_squared': str(gradient),
        'normalized_gradient_times_side_squared_upper': str(F(7, 30)*dual_energy),
        'exact_normalized_interpolation_mass_error': str(mass),
        'normalized_interpolation_mass_error_upper': str(F(5, 18*len(faces))*dual_energy)}


def tent_geometry_control():
    """Exact local honeycomb tensor and continuous barycentric tent moments."""
    root = sp.sqrt(3)
    basis = sp.Matrix([[1, sp.Rational(1, 2)], [0, root/2]])
    offsets = [basis*sp.Matrix(v)/3 for v in ((1, 1), (-2, 1), (1, -2))]
    tensor = sum((v*v.T for v in offsets), sp.zeros(2))
    if sp.simplify(tensor-sp.eye(2)/2) != sp.zeros(2):
        raise ArithmeticError('the exact honeycomb centroid energy tensor changed')
    # One of the three tent regions has corners (0,0),(1,0),(1/3,1/3).
    # The affine simplex coordinates have determinant1/3 and tent value v.
    u, v = sp.symbols('u v', nonnegative=True)
    moments = [sp.integrate(sp.integrate(v**power/3, (v, 0, 1-u)), (u, 0, 1))*6
               for power in (1, 2)]
    if moments != [sp.Rational(1, 3), sp.Rational(1, 6)]:
        raise ArithmeticError('the integrated barycentric tent moments changed')
    gradient_squared = sp.simplify((3/(root/2))**2)
    if gradient_squared != 12:
        raise ArithmeticError('the equilateral tent gradient changed')
    quotient = gradient_squared/(48*moments[1])
    if quotient != sp.Rational(3, 2):
        raise ArithmeticError('the continuum tent energy quotient changed')
    return {'unit_side_centroid_tensor': 'I_2/2',
        'normalized_tent_mean': str(moments[0]),
        'normalized_tent_second_moment': str(moments[1]),
        'unit_side_tent_gradient_squared': str(gradient_squared),
        'native_limit_tent_energy_quotient_at_coarse_level_zero': str(quotient)}


def tent_refinement_control(level):
    """Exact source-mesh energies for all twenty continuous ancestral tents."""
    level = _level(level)
    if level > 4:
        raise ValueError('finite tent refinement controls use levels zero through four')
    m = 2**level
    vertices = ((0, 0, 0, 0, 0, 0), (2, 0, 0, 0, 0, 0), (0, 2, 0, 0, 0, 0))
    faces = ((0, 1, 2),)
    for _ in range(level):
        vertices, faces = uniform.refine(vertices, faces)
    numerators = []
    for face in faces:
        x = sum(vertices[p][0]//2 for p in face)
        y = sum(vertices[p][1]//2 for p in face)
        numerators.append(min(x, y, 3*m-x-y))
    if any(not 0 < value <= m for value in numerators):
        raise ArithmeticError('a centroid is outside its ancestral tent')
    _, actual_faces = uniform.affine_prefix(level)
    edges = uniform.dual_edges(actual_faces)
    count = 20*m*m
    # Every score is at most m, and there are30*m² edges. This bound
    # certifies all integer products and partial sums before they execute.
    if 30*m**4 >= 2**40:
        raise ArithmeticError('the finite exact tent control exceeds its arithmetic budget')
    scores = np.zeros((count, 20), dtype=np.int64)
    for base in range(20):
        scores[base*m*m:(base+1)*m*m, base] = numerators
    energy = np.zeros((20, 20), dtype=np.int64)
    for a, b in edges:
        difference = scores[a]-scores[b]
        energy += np.outer(difference, difference)
    gram = scores.T@scores
    if (not np.array_equal(gram, np.eye(20, dtype=np.int64)*sum(x*x for x in numerators))
            or np.any(energy.sum(axis=1) < 0)):
        raise ArithmeticError('continuous tent norm or source energy control failed')
    diagonal = tuple(F(int(energy[i, i]), 20*m*m) for i in range(20))
    off_diagonal = max(abs(F(int(energy[i, j]), 20*m*m))
                       for i in range(20) for j in range(20) if i != j)
    return {'level': level, 'source_cells': count,
        'exact_normalized_patch_tent_mean': str(F(sum(numerators), m**3)),
        'exact_normalized_patch_tent_second_moment': str(F(sum(x*x for x in numerators), m**4)),
        'max_dual_energy_diagonal_error_from_three_twentieths': str(max(abs(x-F(3, 20)) for x in diagonal)),
        'max_dual_energy_off_diagonal': str(off_diagonal),
        'integer_source_adjacency_and_all_twenty_tent_energies_checked': True}


def full_port_form_control():
    """Exact harmonic lift, normal contraction and cutoff on the native source."""
    graph, _, r, integer_f, denominator = source.matrices(0)
    _, carrier = source.passage.recovered_carrier()
    local_green = (carrier+sp.ones(12)/12).inv()-sp.ones(12)/12
    ri, fi = sp.Matrix(r), sp.Matrix(integer_f)
    values = sp.Matrix([sp.Rational(i % 7-3, 3) for i in range(graph.cells)])
    edge_response = fi.inv()*ri*values*denominator
    drive = sp.zeros(graph.ports, 1)
    for j, (a, b) in enumerate(graph.inter_seams):
        drive[a] += edge_response[j]; drive[b] -= edge_response[j]
    internal = sp.zeros(graph.ports, 1)
    for cell in range(graph.cells):
        internal[12*cell:12*(cell+1), 0] = -local_green*drive[12*cell:12*(cell+1), 0]
    lifted = sp.Matrix([values[p//12]+internal[p] for p in range(graph.ports)])
    def energy(vector):
        return sum((vector[a]-vector[b])**2 for a, b in graph.seams)
    raw_energy = energy(lifted)
    raw_internal = sum(x*x for x in internal)
    if (any(sum(internal[12*c:12*(c+1), 0]) != 0 for c in range(graph.cells))
            or raw_energy != (values.T*ri.T*edge_response)[0]
            or 2*raw_internal > raw_energy):
        raise ArithmeticError('the full-port harmonic energy lift failed')
    contracted = sp.Matrix([max(sp.S.Zero, min(sp.S.One, x)) for x in lifted])
    coarse_contracted = sp.Matrix([max(sp.S.Zero, min(sp.S.One, values[p//12]))
                                   for p in range(graph.ports)])
    contracted_energy = energy(contracted)
    contraction_error = sum((x-y)**2 for x, y in zip(contracted, coarse_contracted, strict=True))
    if contracted_energy > raw_energy or contraction_error > raw_internal:
        raise ArithmeticError('the native full-port normal contraction bound failed')
    # A separate exact field checks the product inequality without assuming
    # that a cell-mean Schur matrix is a Markov generator.
    test = sp.Matrix([sp.Rational(7*p % 19-9, 5) for p in range(graph.ports)])
    cutoff = [sp.Rational((p//12) % 5, 4) for p in range(graph.ports)]
    altered = sp.Matrix([x*y for x, y in zip(test, cutoff, strict=True)])
    transition = {p for a, b in graph.inter_seams if cutoff[a] != cutoff[b] for p in (a, b)}
    omega = max((cutoff[a]-cutoff[b])**2 for a, b in graph.inter_seams)
    epsilon = sp.Rational(1, 3)
    allowance = ((1+epsilon)*energy(test)
                 +(1+1/epsilon)*omega*sum(test[p]**2 for p in transition))
    if energy(altered) > allowance:
        raise ArithmeticError('the exact native cutoff product bound failed')
    return {'level': 0, 'ports': graph.ports,
        'full_port_harmonic_energy_identity_checked': True,
        'internal_norm_squared_over_raw_energy': str(raw_internal/raw_energy),
        'contracted_energy_over_raw_energy': str(contracted_energy/raw_energy),
        'contraction_error_squared_over_internal_norm_squared': str(contraction_error/raw_internal),
        'cutoff_energy': str(energy(altered)), 'cutoff_energy_upper': str(allowance),
        'constant_zero_mode_adjoined_before_normal_contraction': True}


def separated_action_bound(level, cell_distance):
    """Rescaled Schur bilinear bound for unit vectors at cell distance d>=1."""
    level = _level(level)
    if type(cell_distance) is not int or cell_distance < 1:
        raise ValueError('a positive integer cell distance is required')
    return F(4**level, 2**cell_distance)


def macroscopic_action_exponent(level, separation):
    """Base-two exponent of the separated-support action upper bound.

    Separation is intrinsic distance on the auxiliary unit-circumradius flat
    polyhedron.  Adjacent centroids are less than 2^-l apart; each fine cell
    has diameter less than 2*2^-l.  This is a comparison length, not metres.
    An exponent avoids allocating a doubly exponential denominator.
    """
    level = _level(level)
    if type(separation) is not F or separation <= 0:
        raise ValueError('a positive exact Fraction separation is required')
    lower = separation*2**level-4
    distance = max(1, -((-lower.numerator)//lower.denominator))
    if lower <= 0:
        raise ValueError('the cells are not yet certified disjoint at this level')
    return 2*level-distance


def triangle_patch(depth):
    """Exact upright/downward cells and dual edges in one descendant triangle."""
    depth = _level(depth)
    if depth > 5:
        raise ValueError('the finite stencil control is limited to depth five')
    m = 2**depth
    up = tuple((i, j) for i in range(m) for j in range(m-i))
    down = tuple((i, j) for i in range(m-1) for j in range(m-1-i))
    vertices = tuple(('u', *v) for v in up)+tuple(('d', *v) for v in down)
    index = {v: i for i, v in enumerate(vertices)}
    edges = tuple(sorted(tuple(sorted((index[('d', i, j)], index[('u', x, y)])))
        for i, j in down for x, y in ((i, j), (i+1, j), (i, j+1))))
    return m, vertices, edges


def reflected_square_control(depth):
    """The square has exactly twice each upright-triangle edge energy."""
    m, vertices, edges = triangle_patch(depth)
    def reflect(i, j):
        return (i, j) if i+j < m else (m-1-j, m-1-i)
    upright_edges = Counter()
    for i in range(m):
        for j in range(m-i):
            for x, y in ((i+1, j), (i, j+1)):
                if x+y < m:
                    upright_edges[tuple(sorted(((i, j), (x, y))))] = 2
    reflected_edges = Counter()
    for i in range(m):
        for j in range(m):
            for x, y in ((i+1, j), (i, j+1)):
                if x < m and y < m:
                    reflected_edges[tuple(sorted((reflect(i, j), reflect(x, y))))] += 1
    if reflected_edges != upright_edges:
        raise ArithmeticError('square reflection changed the exact edge multiplicities')
    multiplicity = Counter(reflect(i, j) for i in range(m) for j in range(m))
    if (set(multiplicity) != {(i, j) for kind, i, j in vertices if kind == 'u'}
            or not all(n in (1, 2) for n in multiplicity.values())):
        raise ArithmeticError('upright cells are not covered once or twice')
    return {'depth': depth, 'side_subdivisions': m, 'descendant_cells': len(vertices),
        'internal_dual_edges': len(edges), 'square_edges': sum(reflected_edges.values()),
        'square_edge_energy_is_twice_upright_energy': True,
        'all_upright_cells_have_multiplicity_one_or_two': True}


def source_patch_control(depth):
    """Bind every base-face descendant adjacency to the independent stencil."""
    if _level(depth) > 3:
        raise ValueError('native ancestry binding controls use depths zero through three')
    m, vertices, edges = triangle_patch(depth)
    # The first two slots here are independent combinatorial triangle
    # coordinates; they are not the physical x/y slots of the carrier chart.
    simple = ((0, 0, 0, 0, 0, 0), (2, 0, 0, 0, 0, 0), (0, 2, 0, 0, 0, 0))
    faces = ((0, 1, 2),)
    for _ in range(depth):
        simple, faces = uniform.refine(simple, faces)
    lookup = {v: i for i, v in enumerate(vertices)}
    labels = []
    for face in faces:
        corners = tuple((simple[v][0]//2, simple[v][1]//2) for v in face)
        i, j = min(x for x, _ in corners), min(y for _, y in corners)
        kind = 'u' if (i, j) in corners else 'd'
        labels.append(lookup[(kind, i, j)])
    _, actual_faces = uniform.affine_prefix(depth)
    actual_edges = uniform.dual_edges(actual_faces)
    expected = set(edges)
    for base in range(20):
        lower, upper = base*m*m, (base+1)*m*m
        internal = {tuple(sorted((labels[a-lower], labels[b-lower])))
                    for a, b in actual_edges if lower <= a < upper and lower <= b < upper}
        if internal != expected:
            raise ArithmeticError('native descendants differ from the triangular patch stencil')
    return {'depth': depth, 'base_patches_checked': 20,
        'cell_ancestry_order_and_all_internal_edges_checked': True}


def exact_patch_energy_control(depth):
    """Independent rational LDL certificate for the full patch variance bound."""
    m, vertices, edges = triangle_patch(depth)
    if depth > 3:
        raise ValueError('exact matrix energy controls are limited to depth three')
    n = len(vertices)
    laplacian = sp.zeros(n)
    for a, b in edges:
        laplacian[a, a] += 1; laplacian[b, b] += 1
        laplacian[a, b] -= 1; laplacian[b, a] -= 1
    constant = sp.ones(n)/n
    certificate = 26*m*m*laplacian-sp.eye(n)+2*constant
    # Constant direction has eigenvalue one; the other directions certify
    # 26*m² L_patch >= I-P_constants without a floating eigenvalue test.
    _, pivots = certificate.LDLdecomposition(hermitian=False)
    if any(pivots[i, i] <= 0 for i in range(n)):
        raise ArithmeticError('the exact triangular patch Poincare bound failed')
    return {'depth': depth, 'cells': n,
        'variance_energy_multiplier': 26*m*m,
        'rational_LDL_positive_pivots': n}


def native_operator_control(level):
    """Native support, compactness and separated-action finite controls."""
    if _level(level) > 2:
        raise ValueError('native operator controls use levels zero through two')
    graph, _, r, integer_f, denominator = source.matrices(level)
    cell_l = r.T@r
    adjacency = -(cell_l-np.diag(np.diag(cell_l)))
    distances = shortest_path(adjacency, directed=False, unweighted=True)
    numerator_t = 2*denominator*np.eye(len(integer_f), dtype=np.int64)-integer_f
    term = r.copy()
    support_rows = []
    row_norm = max(sum(abs(int(value)) for value in row) for row in numerator_t)
    incidence_column_norm = max(sum(abs(int(value)) for value in column) for column in r.T)
    integer_limit = int(np.iinfo(np.int64).max)
    for power in range(4):
        term_norm = max(abs(int(value)) for value in term.flat)
        if incidence_column_norm*term_norm > integer_limit:
            raise ArithmeticError('exact polynomial coefficient could overflow int64')
        coefficient = r.T@term
        if np.any(coefficient[distances > power+1]):
            raise ArithmeticError('Schur polynomial has a forbidden long-range coefficient')
        support_rows.append({'power': power, 'maximum_cell_range': power+1,
            'integer_polynomial_support_checked': True})
        if power < 3:
            if row_norm*term_norm > integer_limit:
                raise ArithmeticError('exact polynomial power could overflow int64')
            term = numerator_t@term
    schur = r.T@np.linalg.solve(integer_f/denominator, r)/12
    green = source._green(schur)/4**level
    compact_rows = []
    for coarse in range(level+1):
        descendants = 4**(level-coarse)
        projected = np.repeat(green.reshape(20*4**coarse, descendants, graph.cells).mean(axis=1),
                              descendants, axis=0)
        squared = float(np.linalg.norm(green-projected, 2)**2)
        if squared > float(compactness_tail_squared(coarse))+1e-9:
            raise ArithmeticError('native Green violates the coarse compactness bound')
        compact_rows.append({'coarse_level': coarse,
            'tail_norm_squared_diagnostic': squared,
            'tail_norm_squared_upper': str(compactness_tail_squared(coarse))})
    distance = int(distances.max())
    a, b = map(int, np.argwhere(distances == distance)[0])
    coupling = float(abs(4**level*schur[a, b]))
    bound = separated_action_bound(level, distance)
    if coupling > float(bound)+1e-12:
        raise ArithmeticError('native distant action exceeded its exact bound')
    return {'level': level, 'cells': graph.cells,
        'exact_polynomial_support': support_rows, 'compactness_controls': compact_rows,
        'distant_native_pair': [a, b], 'cell_distance': distance,
        'absolute_rescaled_coupling_diagnostic': coupling,
        'absolute_rescaled_coupling_upper': str(bound)}


def non_markov_control():
    """Averages are not a vertex Kron reduction: a positive off-diagonal survives."""
    _, _, r, integer_f, denominator = source.matrices(0)
    ri, fi = sp.Matrix(r), sp.Matrix(integer_f)
    entry = (ri[:, 5].T*fi.inv()*ri[:, 9])[0]*denominator/12
    if entry <= 0:
        raise ArithmeticError('the actual source mean-elimination countercontrol changed')
    return {'level': 0, 'cells': [5, 9], 'exact_positive_Schur_off_diagonal': str(entry),
        'finite_cell_mean_generator_is_a_Markov_generator': False,
        'continuum_locality_or_nonlocality_decided_by_this_control': False}


@lru_cache(maxsize=3)
def _green_orthonormal_bundle(level):
    graph, _, r, integer_f, denominator = source.matrices(level)
    schur = r.T@np.linalg.solve(integer_f/denominator, r)/12.0
    cells = graph.cells
    projector = np.eye(cells, dtype=float)-np.ones((cells, cells), dtype=float)/cells
    generator = (4**level)*(projector@schur@projector)
    generator = 0.5*(generator+generator.T)
    green = source._green(schur)/4**level
    green = 0.5*(green+green.T)
    return green, generator, projector


def green_orthonormal_matrix(level=0):
    """Normalized cell-mean Green operator on the orthonormal cell basis."""
    level = _int(level, minimum=0, maximum=2)
    green, _, _ = _green_orthonormal_bundle(level)
    return green.copy()


def _denman_beavers_sqrt(matrix):
    """Independent Denman-Beavers square root and inverse square root on SPD input."""
    n = len(matrix)
    y = matrix.astype(float, copy=True)
    z = np.eye(n, dtype=float)
    identity = np.eye(n, dtype=float)
    for _ in range(32):
        y_inv = np.linalg.solve(y, identity)
        z_inv = np.linalg.solve(z, identity)
        y_next = 0.5*(y+z_inv)
        z_next = 0.5*(z+y_inv)
        tolerance = 8*np.finfo(float).eps*max(1.0, np.linalg.norm(y_next, 'fro'))
        if np.linalg.norm(y_next-y, 'fro') <= tolerance:
            y = 0.5*(y_next+y_next.T)
            z = 0.5*(z_next+z_next.T)
            return y, z
        y, z = y_next, z_next
    raise ArithmeticError('Denman-Beavers matrix square-root iteration did not converge')


def _matrix_exp_scaling_squaring(matrix):
    """Independent pure-NumPy Taylor scaling-and-squaring matrix exponential."""
    n = len(matrix)
    norm = float(np.linalg.norm(matrix, 1))
    scale = 0
    while norm > 0.5*(2**scale):
        scale += 1
    scaled = matrix/(2.0**scale)
    result = np.eye(n, dtype=float)
    term = np.eye(n, dtype=float)
    for k in range(1, 20):
        term = (term@scaled)/k
        result = result+term
        if np.linalg.norm(term, 1) <= 1e-16:
            break
    for _ in range(scale):
        result = result@result
    return result


def steklov_bulk_extension_bounds(level=0):
    """Radial Poisson-Steklov half-cylinder extension, trace filter and symbol bounds."""
    level = _int(level, minimum=0, maximum=2)
    green, generator, projector = _green_orthonormal_bundle(level)
    cells = len(green)
    constant_projector = np.ones((cells, cells), dtype=float)/cells
    sqrt_shifted, inv_sqrt_shifted = _denman_beavers_sqrt(generator+constant_projector)
    a_half_db = projector@sqrt_shifted@projector
    g_half_db = projector@inv_sqrt_shifted@projector
    a_quarter_shifted, _ = _denman_beavers_sqrt(a_half_db+constant_projector)
    a_quarter_db = projector@a_quarter_shifted@projector

    eigenvalues, eigenvectors = np.linalg.eigh(generator)
    positive = eigenvalues > 1e-10
    if int(np.sum(positive)) != cells-1:
        raise ArithmeticError('scaled cell-mean generator lost its mean-zero rank')
    pos_evals = eigenvalues[positive]
    pos_evecs = eigenvectors[:, positive]
    a_half_spec = (pos_evecs*np.sqrt(pos_evals))@pos_evecs.T
    g_half_spec = (pos_evecs*(1.0/np.sqrt(pos_evals)))@pos_evecs.T
    a_quarter_spec = (pos_evecs*(pos_evals**0.25))@pos_evecs.T

    sqrt_generator_error = float(np.linalg.norm(a_half_db@a_half_db-generator, 2))
    sqrt_green_error = float(np.linalg.norm(g_half_db@g_half_db-green, 2))
    sqrt_duality_error = float(np.linalg.norm(a_half_db@g_half_db-projector, 2))
    spectral_vs_db_error = float(max(
        np.linalg.norm(a_half_db-a_half_spec, 2),
        np.linalg.norm(g_half_db-g_half_spec, 2),
        np.linalg.norm(a_quarter_db-a_quarter_spec, 2),
    ))
    filter_identity_error = float(np.linalg.norm(
        a_quarter_db@green@a_quarter_db-g_half_db, 2))
    if (sqrt_generator_error > 1e-11 or sqrt_green_error > 1e-11
            or sqrt_duality_error > 1e-11 or spectral_vs_db_error > 1e-11
            or filter_identity_error > 1e-11):
        raise ArithmeticError('Steklov square-root or A^(1/4) filter identity failed')

    raw_f = np.array([float((11*i*i+5*i+3) % 29-14) for i in range(cells)], dtype=float)
    raw_g = np.array([float((7*i*i-13*i+9) % 31-15) for i in range(cells)], dtype=float)
    boundary_f = projector@raw_f
    boundary_f /= float(np.linalg.norm(boundary_f, 2))
    boundary_g = projector@raw_g
    boundary_g /= float(np.linalg.norm(boundary_g, 2))

    nodes, weights = np.polynomial.legendre.leggauss(64)
    radial_energy_f = 0.0
    tangential_energy_f = 0.0
    cross_energy_fg = 0.0
    for t, weight in zip(nodes, weights, strict=True):
        radius = (1.0+float(t))/(1.0-float(t))
        jacobian = 2.0/((1.0-float(t))**2)
        propagator = _matrix_exp_scaling_squaring(-radius*a_half_db)
        u_r = propagator@boundary_f
        du_r = -(a_half_db@u_r)
        v_r = propagator@boundary_g
        dv_r = -(a_half_db@v_r)
        radial_energy_f += float(weight)*jacobian*float(du_r@du_r)
        tangential_energy_f += float(weight)*jacobian*float(u_r@(generator@u_r))
        cross_energy_fg += float(weight)*jacobian*float(
            du_r@dv_r+u_r@(generator@v_r))

    quadrature_dirichlet_energy = radial_energy_f+tangential_energy_f
    boundary_steklov_energy_db = float(boundary_f@(a_half_db@boundary_f))
    boundary_steklov_energy_spec = float(boundary_f@(a_half_spec@boundary_f))
    boundary_cross_energy_db = float(boundary_f@(a_half_db@boundary_g))
    dirichlet_energy_error = abs(quadrature_dirichlet_energy-boundary_steklov_energy_db)
    equipartition_error = abs(radial_energy_f-tangential_energy_f)
    polarized_cross_error = abs(cross_energy_fg-boundary_cross_energy_db)
    if (dirichlet_energy_error > 1e-10
            or abs(boundary_steklov_energy_db-boundary_steklov_energy_spec) > 1e-12
            or equipartition_error > 1e-10
            or polarized_cross_error > 1e-10):
        raise ArithmeticError('radial Poisson-Steklov Dirichlet energy identity failed')

    compactness_rows = []
    for coarse in range(level+1):
        descendants = 4**(level-coarse)
        blocks = 20*4**coarse
        coarse_proj = np.zeros((cells, cells), dtype=float)
        for block in range(blocks):
             sl = slice(block*descendants, (block+1)*descendants)
             coarse_proj[sl, sl] = 1.0/descendants
        residual_proj = np.eye(cells, dtype=float)-coarse_proj
        sqrt_tail_sq = float(np.linalg.norm(residual_proj@g_half_db, 2)**2)
        sandwiched_green = float(np.linalg.norm(residual_proj@green@residual_proj, 2))
        green_tail_norm = float(np.linalg.norm(residual_proj@green, 2))
        energy_bound = covariance_eigenvalue_tail(coarse)
        green_bound = compactness_tail_squared(coarse)
        if (abs(sqrt_tail_sq-sandwiched_green) > 1e-11
                or sqrt_tail_sq > green_tail_norm+1e-11
                or sqrt_tail_sq > float(energy_bound)+1e-9
                or sqrt_tail_sq > float(green_bound)+1e-9):
            raise ArithmeticError('Steklov boundary covariance tail violates compactness bound')
        compactness_rows.append({
            'coarse_level': coarse,
            'sqrt_green_tail_norm_squared_diagnostic': sqrt_tail_sq,
            'sandwiched_green_tail_norm_diagnostic': sandwiched_green,
            'energy_tail_squared_upper': str(energy_bound),
            'green_tail_squared_upper': str(green_bound),
        })

    weighted_bundle = weighted_conductance_bundle(level)
    g_weighted = weighted_bundle['weighted_green']
    _, g_weighted_half_shifted = _denman_beavers_sqrt(
        (4**level)*(projector@(weighted_bundle['weighted_laplacian']/12.0)@projector)
        +constant_projector)
    g_weighted_half = projector@g_weighted_half_shifted@projector
    sqrt_diff_norm = float(np.linalg.norm(g_half_db-g_weighted_half, 2))
    green_diff_norm = float(np.linalg.norm(green-g_weighted, 2))
    square_root_upper = float(np.sqrt(green_diff_norm))
    analytic_square_root_upper = float(np.sqrt(float(weighted_green_difference_bound(level))))
    if (sqrt_diff_norm > square_root_upper+1e-12
            or sqrt_diff_norm > analytic_square_root_upper+1e-12):
        raise ArithmeticError('operator-monotone square-root bound failed')

    ell = sp.Symbol('ell', integer=True, positive=True)
    ratio_identity = sp.simplify(ell**2/(ell*(ell+1))-(1-1/(ell+1)))
    gap_rationalized = sp.simplify(
        (sp.sqrt(ell*(ell+1))-ell)-ell/(sp.sqrt(ell*(ell+1))+ell))
    if ratio_identity != 0 or gap_rationalized != 0:
        raise ArithmeticError('round-sphere vs 3-ball Steklov symbol identity changed')
    symbol_rows = []
    for degree in range(1, 33):
        ratio = F(degree, degree+1)
        gap = float(np.sqrt(degree*(degree+1))-degree)
        if not (F(1, 2) <= ratio < 1 and 0.0 < gap < 0.5):
            raise ArithmeticError('Steklov symbol comparison left its certified enclosure')
        if degree in (1, 2, 4, 8, 16, 32):
            symbol_rows.append({
                'degree': degree,
                'squared_symbol_ratio_exact': str(ratio),
                'symbol_difference_diagnostic': gap,
            })

    return {
        'level': level,
        'cells': cells,
        'quadrature_nodes': 64,
        'quadrature_dirichlet_energy_diagnostic': quadrature_dirichlet_energy,
        'radial_bulk_energy_diagnostic': radial_energy_f,
        'tangential_bulk_energy_diagnostic': tangential_energy_f,
        'boundary_steklov_energy_denman_beavers': boundary_steklov_energy_db,
        'boundary_steklov_energy_spectral': boundary_steklov_energy_spec,
        'dirichlet_energy_quadrature_error': dirichlet_energy_error,
        'radial_tangential_equipartition_error': equipartition_error,
        'polarized_bilinear_form_error': polarized_cross_error,
        'denman_beavers_vs_spectral_error': spectral_vs_db_error,
        'A_quarter_spectral_filter_identity_error': filter_identity_error,
        'compactness_tail_controls': compactness_rows,
        'square_root_sqrt_difference_diagnostic': sqrt_diff_norm,
        'square_root_sqrt_difference_upper_from_green': square_root_upper,
        'square_root_sqrt_difference_analytic_upper': analytic_square_root_upper,
        'round_sphere_vs_ball_symbol_bounds': {
            'squared_symbol_ratio_lower': '1/2',
            'squared_symbol_ratio_strict_upper': '1',
            'symbol_difference_strict_lower': '0',
            'symbol_difference_degree_one_exact': 'sqrt(2)-1',
            'symbol_difference_strict_upper': '1/2',
            'sampled_degrees': symbol_rows,
        },
        'boundary_trace_achieved_by_A_quarter_spectral_filter': True,
        'unfiltered_boundary_trace_of_smooth_3d_local_ball_action': False,
        'identifies_physical_radial_bulk_coordinate': False,
        'identifies_limit_with_round_sphere': False,
    }


@lru_cache(maxsize=3)
def _weighted_conductance_cache(level):
    graph, _, r, integer_f, denominator = source.matrices(level)
    _, _, edges = weighted_dual_network(graph)
    w_min = F(5, 8)
    w_max = F(15, 19)
    w_star = 2/(1/w_min+1/w_max)
    if w_star != F(30, 43):
        raise ArithmeticError('harmonic-mean isotropic reference conductance changed')
    cells = graph.cells
    edge_count = len(edges)
    arithmetic_mean = sum((w for _, _, w in edges), F(0))/edge_count
    histogram = {str(w): count for w, count in sorted(Counter(w for _, _, w in edges).items())}
    weighted_laplacian = np.zeros((cells, cells), dtype=float)
    isotropic_laplacian = np.zeros((cells, cells), dtype=float)
    w_star_float = float(w_star)
    for a, b, w in edges:
        wf = float(w)
        weighted_laplacian[a, a] += wf
        weighted_laplacian[b, b] += wf
        weighted_laplacian[a, b] -= wf
        weighted_laplacian[b, a] -= wf
        isotropic_laplacian[a, a] += w_star_float
        isotropic_laplacian[b, b] += w_star_float
        isotropic_laplacian[a, b] -= w_star_float
        isotropic_laplacian[b, a] -= w_star_float
    schur = r.T@np.linalg.solve(integer_f/denominator, r)/12.0
    native_green = source._green(schur)/4**level
    weighted_green = source._green(weighted_laplacian/12.0)/4**level
    isotropic_green = source._green(isotropic_laplacian/12.0)/4**level
    return {
        'level': level,
        'cells': cells,
        'edge_count': edge_count,
        'edges': edges,
        'conductance_min': w_min,
        'conductance_max': w_max,
        'harmonic_mean_conductance': w_star,
        'arithmetic_mean_conductance': arithmetic_mean,
        'conductance_histogram': histogram,
        'weighted_laplacian': weighted_laplacian,
        'isotropic_laplacian': isotropic_laplacian,
        'weighted_green': weighted_green,
        'isotropic_green': isotropic_green,
        'native_green': native_green,
    }


def weighted_conductance_bundle(level=0):
    """Exact dual conductances and normalized weighted/isotropic operators."""
    level = _int(level, minimum=0, maximum=2)
    cached = _weighted_conductance_cache(level)
    return {
        'level': cached['level'],
        'cells': cached['cells'],
        'edge_count': cached['edge_count'],
        'edges': cached['edges'],
        'conductance_min': cached['conductance_min'],
        'conductance_max': cached['conductance_max'],
        'harmonic_mean_conductance': cached['harmonic_mean_conductance'],
        'arithmetic_mean_conductance': cached['arithmetic_mean_conductance'],
        'conductance_histogram': dict(cached['conductance_histogram']),
        'weighted_laplacian': cached['weighted_laplacian'].copy(),
        'isotropic_laplacian': cached['isotropic_laplacian'].copy(),
        'weighted_green': cached['weighted_green'].copy(),
        'isotropic_green': cached['isotropic_green'].copy(),
        'native_green': cached['native_green'].copy(),
    }


@lru_cache(maxsize=1)
def _exact_level_zero_loewner_certificate():
    """Exact rational SymPy LDL certificate for the level-zero Loewner sandwiches."""
    bundle = _weighted_conductance_cache(0)
    cells = bundle['cells']
    w_star = sp.Rational(30, 43)
    l_w = sp.zeros(cells)
    l_iso = sp.zeros(cells)
    for a, b, w in bundle['edges']:
        wr = sp.Rational(w.numerator, w.denominator)
        l_w[a, a] += wr; l_w[b, b] += wr
        l_w[a, b] -= wr; l_w[b, a] -= wr
        l_iso[a, a] += w_star; l_iso[b, b] += w_star
        l_iso[a, b] -= w_star; l_iso[b, a] -= w_star
    basis = sp.zeros(cells, cells-1)
    for i in range(cells-1):
        basis[i, i] = 1
        basis[i+1, i] = -1
    lap_lower = basis.T*(l_w-sp.Rational(43, 48)*l_iso)*basis
    lap_upper = basis.T*(sp.Rational(43, 38)*l_iso-l_w)*basis
    ones = sp.ones(cells)/cells
    g_w = (l_w/12+ones).inv()-ones
    g_iso = (l_iso/12+ones).inv()-ones
    green_lower = basis.T*(g_w-sp.Rational(38, 43)*g_iso)*basis
    green_upper = basis.T*(sp.Rational(48, 43)*g_iso-g_w)*basis
    for matrix in (lap_lower, lap_upper, green_lower, green_upper):
        _, pivots = matrix.LDLdecomposition(hermitian=False)
        if any(pivots[i, i] <= 0 for i in range(cells-1)):
            raise ArithmeticError('exact rational Loewner LDL certificate failed at level zero')
    return cells-1


def isotropic_comparison_bounds(level=0):
    """Exact Loewner sandwiches and uniform error bounds against the isotropic dual."""
    level = _int(level, minimum=0, maximum=2)
    bundle = weighted_conductance_bundle(level)
    cells = bundle['cells']
    edges = bundle['edges']
    w_min = bundle['conductance_min']
    w_max = bundle['conductance_max']
    w_star = bundle['harmonic_mean_conductance']
    lap_lower_mult = w_min/w_star
    lap_upper_mult = w_max/w_star
    green_lower_mult = w_star/w_max
    green_upper_mult = w_star/w_min
    rel_bound = (green_upper_mult-1)
    if (lap_lower_mult != F(43, 48) or lap_upper_mult != F(43, 38)
            or green_lower_mult != F(38, 43) or green_upper_mult != F(48, 43)
            or 1-green_lower_mult != F(5, 43) or rel_bound != F(5, 43)):
        raise ArithmeticError('exact Loewner multiplier arithmetic changed')

    lower_slacks = [w-lap_lower_mult*w_star for _, _, w in edges]
    upper_slacks = [lap_upper_mult*w_star-w for _, _, w in edges]
    if any(slack < 0 for slack in lower_slacks) or any(slack < 0 for slack in upper_slacks):
        raise ArithmeticError('a dual edge conductance violated the Loewner sandwich')

    ldl_pivots = _exact_level_zero_loewner_certificate()

    l_w = bundle['weighted_laplacian']/12.0
    l_iso = bundle['isotropic_laplacian']/12.0
    g_w = bundle['weighted_green']
    g_iso = bundle['isotropic_green']
    g_native = bundle['native_green']
    q_basis, _ = np.linalg.qr(
        np.eye(cells, dtype=float)[:, :cells-1]
        -np.ones((cells, cells-1), dtype=float)/cells)

    def min_mean_zero_eig(sym_matrix):
        reduced = q_basis.T@sym_matrix@q_basis
        reduced = 0.5*(reduced+reduced.T)
        return float(np.linalg.eigvalsh(reduced)[0])

    lap_low_eig = min_mean_zero_eig(l_w-float(lap_lower_mult)*l_iso)
    lap_high_eig = min_mean_zero_eig(float(lap_upper_mult)*l_iso-l_w)
    green_low_eig = min_mean_zero_eig(g_w-float(green_lower_mult)*g_iso)
    green_high_eig = min_mean_zero_eig(float(green_upper_mult)*g_iso-g_w)
    if min(lap_low_eig, lap_high_eig, green_low_eig, green_high_eig) < -1e-12:
        raise ArithmeticError('numerical Loewner PSD check failed on mean-zero subspace')

    iso_norm = float(np.linalg.norm(g_iso, 2))
    weighted_iso_diff = float(np.linalg.norm(g_w-g_iso, 2))
    relative_diff = weighted_iso_diff/iso_norm
    native_iso_diff = float(np.linalg.norm(g_native-g_iso, 2))
    native_weighted_bound = weighted_green_difference_bound(level)
    combined_upper = float(native_weighted_bound)+float(rel_bound)*iso_norm
    if relative_diff > float(rel_bound)+1e-12 or native_iso_diff > combined_upper+1e-12:
        raise ArithmeticError('Green comparison norm difference exceeded its bound')

    return {
        'level': level,
        'cells': cells,
        'dual_edges': bundle['edge_count'],
        'conductance_lower': str(w_min),
        'conductance_upper': str(w_max),
        'harmonic_mean_reference_conductance': str(w_star),
        'arithmetic_mean_conductance': str(bundle['arithmetic_mean_conductance']),
        'exact_conductance_histogram': bundle['conductance_histogram'],
        'laplacian_loewner_lower_multiplier': str(lap_lower_mult),
        'laplacian_loewner_upper_multiplier': str(lap_upper_mult),
        'green_loewner_lower_multiplier': str(green_lower_mult),
        'green_loewner_upper_multiplier': str(green_upper_mult),
        'uniform_relative_green_error_upper': str(rel_bound),
        'min_lower_edge_slack': str(min(lower_slacks)),
        'min_upper_edge_slack': str(min(upper_slacks)),
        'all_edge_conductance_slacks_nonnegative': True,
        'exact_rational_level_zero_loewner_LDL_pivots': ldl_pivots,
        'laplacian_lower_slack_min_eigenvalue': lap_low_eig,
        'laplacian_upper_slack_min_eigenvalue': lap_high_eig,
        'green_lower_slack_min_eigenvalue': green_low_eig,
        'green_upper_slack_min_eigenvalue': green_high_eig,
        'relative_green_norm_difference_diagnostic': relative_diff,
        'native_to_isotropic_green_norm_difference_diagnostic': native_iso_diff,
        'native_to_isotropic_green_norm_difference_upper': combined_upper,
        'all_level_loewner_sandwich_certified': True,
        'uniform_relative_green_error_bound_certified': True,
        'continuum_green_loewner_sandwich_certified': True,
        'exact_continuum_metric_isotropy_proved': False,
    }


def run():
    return {'schema': 'oph.native_compact_spatial_generator.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'contract': {'source': 'existing spherical shared-midpoint carrier joins with distinct ports',
            'space': 'mean-zero L2 of normalized auxiliary flat-polyhedron area',
            'observation_dictionary': 'normalized ancestral cell means',
            'limit': 'the same common spatial subsequence as the native Green field',
            'physical_spatial_metric_or_radial_coordinate_supplied_by_this_result': False},
        'exact_patch_controls': [reflected_square_control(k) | source_patch_control(k)
                                 | exact_patch_energy_control(k) for k in range(3)],
        'native_operator_controls': [native_operator_control(k) for k in range(3)],
        'mean_elimination_control': non_markov_control(),
        'full_port_form_control': full_port_form_control(),
        'terminal_star_control': terminal_star_control(),
        'weighted_dual_controls': [weighted_dual_control(k) for k in range(3)],
        'isotropic_comparison_controls': [isotropic_comparison_bounds(k) for k in range(3)],
        'steklov_bulk_extension_controls': [steklov_bulk_extension_bounds(k) for k in range(2)],
        'exact_interpolation_geometry': interpolation_geometry_control(),
        'exact_source_interpolation_controls': [source_interpolation_control(k) for k in range(4)],
        'exact_tent_geometry': tent_geometry_control(),
        'exact_tent_refinement_controls': [tent_refinement_control(k) for k in range(5)],
        'all_level_bounds': {
            'patch_variance_energy_multiplier': '26*4^(l-k)',
            'coarse_projection_energy_multiplier': '624*4^(-k)',
            'Green_norm_upper': str(GREEN_BOUND),
            'Green_tail_squared': {str(k): str(compactness_tail_squared(k)) for k in (0, 4, 20)},
            'Steklov_boundary_covariance_tail_squared_upper': {
                str(k): str(covariance_eigenvalue_tail(k)) for k in (0, 4, 20)},
            'covariance_eigenvalue_number_20_times_4_power_k_upper': {
                str(k): str(covariance_eigenvalue_tail(k)) for k in (0, 4, 20)},
            'separated_cell_support_bilinear_action': '4^l*2^(-d)*||u||*||v||, d>=1',
            'macroscopic_separation_one_quarter_action_base_two_exponents': {
                str(k): macroscopic_action_exponent(k, F(1, 4)) for k in (5, 8, 12)},
            'full_port_internal_tail_per_unit_energy': {
                str(k): str(harmonic_lift_tail(k)) for k in (0, 4, 20)},
            'weighted_dual_Green_difference_upper': {
                str(k): str(weighted_green_difference_bound(k)) for k in (0, 4, 20)},
            'harmonic_mean_reference_conductance': '30/43',
            'weighted_to_isotropic_laplacian_loewner_enclosure': '43/48*L_l^iso <= L_l^w <= 43/38*L_l^iso',
            'weighted_to_isotropic_green_loewner_enclosure': '38/43*G_l^iso <= G_l^w <= 48/43*G_l^iso',
            'uniform_weighted_to_isotropic_relative_green_error_upper': '5/43',
            'continuum_green_loewner_enclosure': '38/43*G_iso <= G <= 48/43*G_iso',
            'interpolation_gradient_squared_upper': '112/a^2*E_l(f_l)',
            'interpolation_mass_error_upper': '20/(3*4^l)*E_l(f_l)',
            'continuum_form_comparison_on_H1': 'a^2/112*integral|grad f|^2 <= E(f) <= a^2/48*integral|grad f|^2',
            'continuum_positive_eigenvalue_enclosures': {
                str(n): {key: str(value) for key, value in continuum_eigenvalue_enclosure(n).items()}
                for n in (1, 19, 20, 79, 80, 1000000)},
            'mean_zero_heat_trace_lower': '1/(exp(3*t/2)-1), t>0',
            'mean_zero_heat_trace_upper': '1/(exp(t/49920)-1), t>0'},
        'theorem': {'common_Green_subsequence_converges_in_operator_norm': True,
            'limit_is_compact_positive_and_injective_on_mean_zero_L2': True,
            'inverse_is_positive_self_adjoint_spatial_generator': True,
            'finite_reduced_forms_Mosco_converge_on_full_L2_with_constant_zero_mode': True,
            'limit_is_conservative_strongly_local_Dirichlet_form': True,
            'form_domain_is_H1_of_the_auxiliary_flat_surface': True,
            'continuous_piecewise_affine_functions_are_a_regular_form_core': True,
            'heat_trace_spectral_dimension': 2,
            'same_dictionary_scaled_compact_smooth_second_order_three_dimensional_elliptic_identification_excluded': True,
            'radial_Poisson_Steklov_half_cylinder_extension_has_boundary_form_A_half': True,
            'steklov_boundary_trace_covariance_is_G_half_via_A_quarter_spectral_filter': True,
            'all_level_isotropic_conductance_loewner_sandwich_and_five_forty_thirds_relative_bound': True,
            'finite_rescaled_action_couplings_at_separated_continuum_supports_vanish': True,
            'explicit_positive_weighted_dual_graph_has_same_continuum_Green_limits': True,
            'weighted_dual_graph_reproduces_finite_native_repair_histories': False,
            'unique_whole_sequence_generator_selected': False,
            'exact_intrinsic_metric_identified': False,
            'limit_identified_as_round_sphere_Laplacian': False,
            'three_dimensional_space_or_Lorentzian_metric_derived': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native compact spatial-generator bounds and source controls verified.')
    else:
        print(canonical(result))


if __name__ == '__main__':
    main()
