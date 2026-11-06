"""Exact certificate for the conditional selector in PROOF.md.

The response dimension, quadratic readout and optimization criterion are
supplied premises. No carrier template enters the sign-matrix enumeration.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path

import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def clean(value):
    if isinstance(value, s.MatrixBase):
        return value.applyfunc(lambda x: s.expand(s.radsimp(x)))
    return s.expand(s.radsimp(value))


def projector(vector):
    vector = s.Matrix(vector)
    norm2 = clean(vector.dot(vector))
    if vector.cols != 1 or norm2.is_positive is not True:
        raise ValueError('a nonzero real vector is required')
    return clean(vector * vector.T / norm2)


def coordinates(matrix):
    matrix = s.Matrix(matrix)
    n = matrix.rows
    if matrix.cols != n or clean(matrix - matrix.T) != s.zeros(n):
        raise ValueError('a real symmetric matrix is required')
    return s.Matrix([matrix[i, i] for i in range(n)] +
                    [s.sqrt(2)*matrix[i, j] for i in range(n) for j in range(i+1, n)])


def information(projectors, weights):
    projectors = [s.Matrix(p) for p in projectors]
    weights = [s.sympify(w) for w in weights]
    if not projectors or len(projectors) != len(weights):
        raise ValueError('one positive weight is required per channel')
    if any(w.is_positive is not True for w in weights) or clean(sum(weights)) != 1:
        raise ValueError('positive weights summing to one are required')
    n = projectors[0].rows
    size = n*(n+1)//2
    total = s.zeros(size)
    for p, w in zip(projectors, weights):
        if p.shape != (n, n) or clean(p-p.T) != s.zeros(n):
            raise ValueError('same-size real orthogonal rank-one projectors required')
        if clean(p*p-p) != s.zeros(n) or clean(s.trace(p)) != 1:
            raise ValueError('orthogonal rank-one projectors required')
        if any(entry.is_real is not True for entry in p):
            raise ValueError('real projectors required')
        v = coordinates(p)
        total += w*v*v.T
    return clean(total)


def isotropic_information(dimension):
    size = dimension*(dimension+1)//2
    tr = s.Matrix([1]*dimension+[0]*(size-dimension))
    return (2*s.eye(size)+tr*tr.T)/(dimension*(dimension+2))


@lru_cache(maxsize=1)
def icosahedral_axes():
    phi = (1+s.sqrt(5))/2
    vectors = [(0, 1, phi), (0, -1, phi), (1, phi, 0),
               (-1, phi, 0), (phi, 0, 1), (phi, 0, -1)]
    return tuple(s.ImmutableMatrix(projector(v)) for v in vectors)


def gram_signs():
    phi = (1+s.sqrt(5))/2
    vectors = [s.Matrix(v) for v in [(0, 1, phi), (0, -1, phi),
               (1, phi, 0), (-1, phi, 0), (phi, 0, 1), (phi, 0, -1)]]
    g = clean(s.Matrix(6, 6, lambda i, j: vectors[i].dot(vectors[j])/(phi+2)))
    seidel = clean(s.sqrt(5)*(g-s.eye(6)))
    assert all(x in (-1, 0, 1) for x in seidel)
    return [[int(x) for x in row] for row in seidel.tolist()]


def square(matrix):
    n = len(matrix)
    return [[sum(matrix[i][k]*matrix[k][j] for k in range(n))
             for j in range(n)] for i in range(n)]


def is_conference(matrix):
    if len(matrix) != 6 or any(len(row) != 6 for row in matrix):
        return False
    if any(matrix[i][j] != matrix[j][i] or matrix[i][j] not in ((0,) if i == j else (-1, 1))
           for i in range(6) for j in range(6)):
        return False
    return square(matrix) == [[5 if i == j else 0 for j in range(6)] for i in range(6)]


def connected_cycle(vertices, edges):
    vertices = set(vertices)
    neighbors = {v: set() for v in vertices}
    for a, b in edges:
        neighbors[a].add(b)
        neighbors[b].add(a)
    if any(len(n) != 2 for n in neighbors.values()):
        return False
    reached, todo = set(), [next(iter(vertices))]
    while todo:
        v = todo.pop()
        if v not in reached:
            reached.add(v)
            todo.extend(neighbors[v]-reached)
    return reached == vertices


def paired_incidence(seidel):
    if not is_conference(seidel):
        raise ValueError('a six-axis symmetric conference matrix is required')
    ports = list(itertools.product(range(6), (-1, 1)))
    edges = {tuple(sorted((a, b))) for a, b in itertools.combinations(range(12), 2)
             if ports[a][0] != ports[b][0]
             and ports[a][1]*ports[b][1]*seidel[ports[a][0]][ports[b][0]] > 0}
    neighbors = {a: {b for b in range(12) if tuple(sorted((a, b))) in edges} for a in range(12)}
    faces = [f for f in itertools.combinations(range(12), 3)
             if all(tuple(sorted(e)) in edges for e in itertools.combinations(f, 2))]
    face_counts = [sum(set(e) <= set(f) for f in faces) for e in sorted(edges)]
    link_cycles = [connected_cycle(neighbors[a],
                   [e for e in edges if set(e) <= neighbors[a]]) for a in range(12)]
    return {'ports': len(ports), 'edges': len(edges), 'triangle_faces': len(faces),
            'degrees': [len(neighbors[a]) for a in range(12)],
            'faces_per_edge': sorted(set(face_counts)),
            'all_vertex_links_are_five_cycles': all(link_cycles),
            'euler_value': len(ports)-len(edges)+len(faces)}


@lru_cache(maxsize=1)
def sign_classification():
    solutions = []
    pairs = list(itertools.combinations(range(1, 6), 2))
    for values in itertools.product((-1, 1), repeat=len(pairs)):
        matrix = [[0 if i == j else 1 for j in range(6)] for i in range(6)]
        for (i, j), value in zip(pairs, values):
            matrix[i][j] = matrix[j][i] = value
        if is_conference(matrix):
            negative_edges = [e for e in pairs if matrix[e[0]][e[1]] == -1]
            assert connected_cycle(range(1, 6), negative_edges)
            incidence = paired_incidence(matrix)
            assert incidence == {'ports': 12, 'edges': 30, 'triangle_faces': 20,
                                 'degrees': [5]*12, 'faces_per_edge': [2],
                                 'all_vertex_links_are_five_cycles': True, 'euler_value': 2}
            solutions.append(matrix)
    assert len(solutions) == 12
    return {'normalized_sign_matrices_checked': 1024,
            'conference_solutions': len(solutions),
            'labelled_five_cycles': 12,
            'every_solution_has_icosahedral_incidence': True}


def stability_arithmetic():
    epsilon, gamma = s.Rational(1, 10**7), s.Rational(2, 15)
    b_squared = epsilon+5*epsilon**2
    b_upper = s.Rational(1, 3000)
    assert b_squared < b_upper**2
    delta = b_upper+4*epsilon
    inverse_error = delta/(gamma*(gamma-epsilon))
    assert inverse_error < s.Rational(1, 50)
    tau, sigma = s.Rational(3, 125), s.Rational(3, 100)
    assert s.Rational(7, 2000)+s.Rational(6, 299) < tau
    assert s.Rational(7, 4)**2 > 3
    assert (s.Rational(3, 500)*375)**2 > 5
    integer_gap_bound = 5*((2+tau)*tau+sigma*(6+2*tau+sigma))
    assert integer_gap_bound == s.Rational(57729, 50000) < 2
    return {key: str(value) for key, value in {
        'epsilon': epsilon, 'optimal_minimum_eigenvalue': gamma,
        'cross_block_squared_bound': b_squared, 'cross_block_norm_upper': b_upper,
        'information_norm_difference_upper': delta,
        'inverse_norm_difference_upper': inverse_error,
        'inverse_norm_simplified_upper': s.Rational(1, 50),
        'weight_lower': s.Rational(50, 301), 'weight_upper': s.Rational(50, 299),
        'squared_overlap_error_upper': s.Rational(1, 375),
        'unweighted_frame_error_upper': tau, 'gram_sign_error_upper': sigma,
        'conference_integer_gap_bound': integer_gap_bound,
        'nonzero_even_integer_gap': 2}.items()}


def witness_certificate():
    projectors = icosahedral_axes()
    info = information(projectors, [s.Rational(1, 6)]*6)
    ideal = isotropic_information(3)
    assert info == ideal and info.rank() == 6
    assert clean(sum(projectors, s.zeros(3))) == 2*s.eye(3)
    overlaps = [clean(s.trace(p*q)) for p, q in itertools.combinations(projectors, 2)]
    assert overlaps == [s.Rational(1, 5)]*15
    inverse = info.inv()
    contractions = clean(s.Matrix(6, 6, lambda i, j:
                         (coordinates(projectors[i]).T*inverse*coordinates(projectors[j]))[0]))
    assert contractions == 6*s.eye(6)
    seidel = gram_signs()
    assert is_conference(seidel)
    return {'response_dimension': 3, 'symmetric_response_dimension': 6,
            'channels': 6, 'weights': ['1/6']*6,
            'information_eigenvalues': {'1/3': 1, '2/15': 5},
            'all_distinct_squared_overlaps': '1/5',
            'inverse_gram_contractions': '6 I_6',
            'paired_positive_overlap_incidence': paired_incidence(seidel)}


def control_certificates():
    axes = icosahedral_axes()
    m0 = isotropic_information(3)
    tetra = [projector(v) for v in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]
    tetra_rank = information(tetra, [s.Rational(1, 4)]*4).rank()
    octa = [projector(v) for v in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]]
    octa_rank = information(octa, [s.Rational(1, 3)]*3).rank()
    assert (tetra_rank, octa_rank) == (4, 3)

    generic = octa+[projector(v) for v in [(1, 1, 0), (1, 0, 1), (0, 1, 1)]]
    generic_info = information(generic, [s.Rational(1, 6)]*6)
    assert generic_info.rank() == 6
    # The unit Frobenius xy response is the fourth coordinate.
    assert generic_info[3, 3] == s.Rational(1, 12) < s.Rational(2, 15)

    unequal = information(axes, [s.Rational(1, 3)]+[s.Rational(2, 15)]*5)
    q0 = coordinates(axes[0])
    assert unequal == clean(s.Rational(4, 5)*m0+s.Rational(1, 5)*q0*q0.T)
    # A positive rank-one perturbation leaves four trace-free gamma modes.
    assert (unequal-s.Rational(8, 75)*s.eye(6)).rank() == 2
    assert unequal.is_positive_definite is True

    rotation = s.Matrix([[s.Rational(3, 5), -s.Rational(4, 5), 0],
                         [s.Rational(4, 5), s.Rational(3, 5), 0], [0, 0, 1]])
    union = list(axes)+[clean(rotation*p*rotation.T) for p in axes]
    assert all(clean(s.trace(p*q)) != 1 for p, q in itertools.combinations(union, 2))
    assert information(union, [s.Rational(1, 12)]*12) == m0

    planar = [projector(v) for v in [(1, 0), (-1, s.sqrt(3)), (-1, -s.sqrt(3))]]
    plane_info = information(planar, [s.Rational(1, 3)]*3)
    assert plane_info == isotropic_information(2)

    t = s.Rational(1, 10**8)
    cayley = s.Matrix([[(1-t*t)/(1+t*t), -2*t/(1+t*t), 0],
                       [2*t/(1+t*t), (1-t*t)/(1+t*t), 0], [0, 0, 1]])
    assert cayley.T*cayley == s.eye(3)
    perturbed = clean(cayley*axes[0]*cayley.T)
    changed_overlap = clean(s.trace(perturbed*axes[1]))
    assert changed_overlap != s.Rational(1, 5)
    perturbation_bound = 4*t/3
    assert perturbation_bound < s.Rational(1, 10**7)
    # Sign preservation uses a vector perturbation bound, independently of
    # the general information-to-incidence theorem.
    assert (2*t)**2 < s.Rational(1, 5)

    return {'tetrahedral_four_axis_rank': tetra_rank,
            'orthogonal_three_axis_rank': octa_rank,
            'nonoptimal_minimal_six_axis_rank': int(generic_info.rank()),
            'nonoptimal_minimal_unit_response_information': '1/12',
            'unequal_icosahedral_weights': ['1/3']+['2/15']*5,
            'unequal_icosahedral_minimum_information': '8/75',
            'nonminimal_optimum_distinct_axes': 12,
            'nonminimal_optimum_paired_ports': 24,
            'without_inverse_port_convention_representatives': 6,
            'two_dimensional_optimum_channels': 3,
            'two_dimensional_optimum_paired_ports': 6,
            'two_dimensional_information_eigenvalues': {'1/2': 1, '1/4': 2},
            'perturbed_optimum': {
                'cayley_parameter': str(t),
                'information_norm_difference_bound': str(perturbation_bound),
                'meets_stability_hypothesis_by_weyl_bound': True,
                'metric_is_equiangular': False,
                'exact_icosahedral_sign_incidence': True}}


def build_certificate():
    return {'schema': 'carrier-selection-reduction-v1',
            'result': 'conditional weighted minimal quadratic readout and stable incidence',
            'proof_domain': 'real dimension three; positive normalized channel weights',
            'selection_premises': ['full symmetric-operator quadratic readout',
                'minimum number of distinct unoriented channels',
                'maximum minimum information eigenvalue at that count'],
            'paired_ports_premise': 'both signs exposed as distinct inverse ports',
            'source_independent_rank_three_derivation': False,
            'physical_native_carrier_selection_theorem': False,
            'lean_kernel_proof': False,
            'weighted_optimum': witness_certificate(),
            'complete_sign_classification': sign_classification(),
            'stable_incidence_exact_bounds': stability_arithmetic(),
            'premise_controls': control_certificates(),
            'source_sha256': {
                str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in (HERE/'certificate.py', HERE/'PROOF.md')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    payload = json.dumps(build_certificate(), indent=2, sort_keys=True)+'\n'
    target = HERE/'results.json'
    if args.check:
        if not target.exists() or target.read_text() != payload:
            raise SystemExit('carrier selection receipt differs from its exact producer')
        print('carrier selection receipt verified')
    else:
        target.write_text(payload)
        print('wrote model/carrier_selection_reduction/results.json')


if __name__ == '__main__':
    main()
