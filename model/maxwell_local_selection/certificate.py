"""Classify local gauge quadratics and reconstruct conservative evolution.

The carrier is read from the committed Lean tables. Quadratic response,
face-sharing locality, symmetry covariance, energy and time reversal are
explicit hypotheses, not consequences of the native record law.
"""
from __future__ import annotations

import argparse
import ast
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re

import sympy as s

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SOURCE_PATHS = (
    'Lean/ObserverPatchHolography/CoreAxioms.lean',
    'Lean/Screen/SeamCurrentCarrierQuotient.lean',
    'Lean/Screen/A5PortAction.lean',
    'Lean/Screen/LocalFaceMaxwellAction.lean',
    'Lean/Screen/TemporalMaxwellEvolution.lean',
    'Lean/Screen/WhitneyMaxwellDynamics.lean',
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exact_matrix(value):
    matrix = s.Matrix(value)
    require(not any(x.has(s.Float) for x in matrix), 'exact matrix entries required')
    return matrix


@lru_cache(maxsize=1)
def carrier():
    text = (ROOT/SOURCE_PATHS[1]).read_text()
    ends = []
    for name in ('seamLeft', 'seamRight'):
        match = re.search(rf'def {name} : Fin 30 → Fin 12 :=\s*!\[([^\]]*)\]', text, re.S)
        require(match is not None, f'missing committed {name}')
        ends.append(tuple(int(x) for x in match[1].replace('\n', ' ').split(',')))
    edges = tuple(zip(*ends, strict=True))
    require(len(edges) == len(set(edges)) == 30 and all(i < j for i, j in edges),
            'committed oriented seam table changed')
    core = (ROOT/SOURCE_PATHS[0]).read_text()
    block = core[core.index('def orientedFaces'):core.index('def faceEdges')]
    faces = tuple(tuple(map(int, row)) for row in re.findall(r'\((\d+),\s*(\d+),\s*(\d+)\)', block))
    require(len(faces) == 20, 'twenty committed faces required')
    group_text = (ROOT/SOURCE_PATHS[2]).read_text()
    group_block = re.search(r'def perms : List \(List Nat\) := (\[.*?\]\])', group_text, re.S)
    require(group_block is not None, 'missing committed rotation table')
    permutations = tuple(tuple(row) for row in ast.literal_eval(group_block[1]))
    require(len(set(permutations)) == 60, 'sixty distinct committed rotations required')
    D, C = s.zeros(30, 12), s.zeros(20, 30)
    edge_index = {edge: i for i, edge in enumerate(edges)}
    for e, (i, j) in enumerate(edges):
        D[e, i], D[e, j] = -1, 1
    for f, (i, j, k) in enumerate(faces):
        for a, b in ((i, j), (j, k), (k, i)):
            C[f, edge_index[tuple(sorted((a, b)))]] = 1 if a < b else -1
    require(C*D == s.zeros(20, 12), 'boundary of boundary is nonzero')
    require(D.rank() == 11 and C.rank() == 19, 'carrier rank drift')
    face_index = {frozenset(face): f for f, face in enumerate(faces)}
    actions = []
    for row in permutations:
        require(sorted(row) == list(range(12)), 'invalid port permutation')
        U, V = s.zeros(30), s.zeros(20)
        face_map = []
        for e, (i, j) in enumerate(edges):
            a, b = row[i], row[j]
            U[edge_index[tuple(sorted((a, b)))], e] = 1 if a < b else -1
        for f, face in enumerate(faces):
            image = tuple(row[i] for i in face)
            target = face_index[frozenset(image)]
            oriented = faces[target]
            require(image in (oriented, oriented[1:]+oriented[:1], oriented[2:]+oriented[:2]),
                    'rotation reversed committed face orientation')
            V[target, f] = 1
            face_map.append(target)
        require(C*U == V*C and U.T*U == s.eye(30), 'source incidence covariance failed')
        actions.append((U, V, tuple(face_map)))
    return D, C, tuple(actions)


def local_pairs(C):
    """Diagonal and pairs contained in one committed triangular face."""
    return tuple((i, j) for i in range(C.cols) for j in range(i, C.cols)
                 if i == j or any(C[f, i]*C[f, j] != 0 for f in range(C.rows)))


def gauge_constraints(D, pairs):
    """Coefficients of H D=0 for every symmetric allowed Hessian entry."""
    result = s.zeros(D.rows*D.cols, len(pairs))
    for variable, (i, j) in enumerate(pairs):
        for p in range(D.cols):
            result[D.cols*i+p, variable] += D[j, p]
            if i != j:
                result[D.cols*j+p, variable] += D[i, p]
    return result


def classify():
    D, C, actions = carrier()
    pairs = local_pairs(C)
    constraints = gauge_constraints(D, pairs)
    face_embedding = s.Matrix([[C[f, i]*C[f, j] for f in range(20)] for i, j in pairs])
    # The analytic proof uses only these exact local triangle incidences.
    for i, j in pairs:
        if i != j:
            require(sum(C[f, i]*C[f, j] != 0 for f in range(20)) == 1,
                    'an off-diagonal pair belongs to several faces')
    triangle_equations = 0
    for f in range(20):
        edges = [e for e in range(30) if C[f, e]]
        vertices = {p for e in edges for p in range(12) if D[e, p]}
        require(len(edges) == len(vertices) == 3, 'nontriangular face')
        for e in edges:
            third = next(p for p in vertices if D[e, p] == 0)
            contributors = [j for j in range(30)
                            if (min(e, j), max(e, j)) in pairs and D[j, third]]
            require(set(contributors) == set(edges)-{e}, 'third-vertex locality failed')
            require(sum(C[f, j]*D[j, third] for j in contributors) == 0,
                    'triangle orientation failed')
            triangle_equations += 1
    gauge_rank, face_rank = constraints.rank(), face_embedding.rank()
    require(len(pairs) == 90 and gauge_rank == 70 and face_rank == 20,
            'local gauge-kernel classification failed')
    require(constraints*face_embedding == s.zeros(360, 20), 'face forms violate gauge invariance')
    # Injectivity of face_embedding plus dimension equality proves that these
    # twenty weighted face forms exhaust the entire local gauge kernel.
    weight_equations = []
    for _, _, face_map in actions:
        for f, target in enumerate(face_map):
            row = [0]*20
            row[f] += 1
            row[target] -= 1
            weight_equations.append(row)
    invariant_rank = s.Matrix(weight_equations).rank()
    require(invariant_rank == 19, 'invariant face weights are not one-dimensional')
    require({mapping[0] for _, _, mapping in actions} == set(range(20)), 'faces not transitive')
    # Averaging the representation is the projector onto invariant vectors.
    # Its exact zero proves there is no invariant linear term in A.
    require(sum((U for U, _, _ in actions), s.zeros(30)) == s.zeros(30),
            'unexpected invariant seam linear term')
    H = C.T*C
    require(all(U.T*H*U == H for U, _, _ in actions), 'selected form lost covariance')
    return {'symmetric_face_local_unknowns': len(pairs),
            'analytic_triangle_equations_checked': triangle_equations,
            'gauge_constraint_rank_exact': gauge_rank,
            'local_gauge_quadratic_dimension_exact': len(pairs)-gauge_rank,
            'weighted_face_embedding_rank_exact': face_rank,
            'rotation_weight_constraint_rank_exact': invariant_rank,
            'rotation_covariant_local_gauge_dimension_exact': 20-invariant_rank,
            'invariant_linear_seam_dimension_exact': 0,
            'complete_form': 'S(A)=constant+(kappa/2)||C A||^2',
            'positive_nonzero_form_iff': 'kappa>0',
            'normalization_selected': False}


def is_local_gauge_hessian(matrix):
    D, C, _ = carrier()
    H = exact_matrix(matrix)
    if H.shape != (30, 30) or H != H.T or H*D != s.zeros(30, 12):
        return False
    pairs = set(local_pairs(C))
    return all(H[i, j] == 0 or (min(i, j), max(i, j)) in pairs
               for i in range(30) for j in range(30))


def evolution(electric=1, magnetic=1):
    """Exact generator for supplied positive electric/magnetic energy weights."""
    electric, magnetic = map(s.sympify, (electric, magnetic))
    require(all(not x.has(s.Float) and x.is_positive is True for x in (electric, magnetic)),
            'exact positive energy coefficients required')
    D, C, _ = carrier()
    G = s.BlockMatrix([[s.zeros(30), magnetic/electric*C.T], [-C, s.zeros(20)]]).as_explicit()
    W = s.diag(*([electric]*30+[magnetic]*20))
    time_reversal = s.diag(*([1]*30+[-1]*20))
    require(G.T*W+W*G == s.zeros(50), 'energy conservation failed')
    require(time_reversal*G*time_reversal == -G, 'time reversal failed')
    require(D.T*G[:30, :] == s.zeros(12, 50), 'Gauss constraint failed')
    require(s.ones(1, 20)*G[30:, :] == s.zeros(1, 50), 'closed flux failed')
    return G, W, time_reversal


def controls():
    D, C, actions = carrier()
    H = C.T*C
    unequal = C.T*s.diag(2, *([1]*19))*C
    require(is_local_gauge_hessian(unequal), 'unequal-weight control not local/gauge')
    require(any(U.T*unequal*U != unequal for U, _, _ in actions), 'covariance mutation undetected')
    nonlocal_form = H*H
    require(nonlocal_form*D == s.zeros(30, 12), 'nonlocal control not gauge invariant')
    require(all(U.T*nonlocal_form*U == nonlocal_form for U, _, _ in actions),
            'nonlocal control not covariant')
    require(not is_local_gauge_hessian(nonlocal_form), 'nonlocal control unexpectedly local')
    G, W, reversal = evolution()
    # A skew electric term preserves energy and Faraday, but violates the
    # declared reversal. Use two actual curl-sector vectors to retain Gauss.
    u, v = C.row(0).T, C.row(1).T
    skew = u*v.T-v*u.T
    gyro = G+s.diag(skew, s.zeros(20))
    require(skew != s.zeros(30) and D.T*skew == s.zeros(12, 30), 'gyro control failed')
    require(gyro.T*W+W*gyro == s.zeros(50), 'gyro energy control failed')
    require(gyro[30:, :] == G[30:, :], 'gyro changed Faraday')
    require(reversal*gyro*reversal != -gyro, 'time reversal mutation undetected')
    G2, _, _ = evolution(1, 2)
    require(G2 != G, 'coefficient freedom disappeared')
    # The carrier action/ODE does not select a finite-step integrator.
    h, lam = s.Rational(1, 2), s.Integer(2)
    oscillator = s.Matrix([[0, 1], [-lam, 0]])
    midpoint = (s.eye(2)-h*oscillator/2).inv()*(s.eye(2)+h*oscillator/2)
    verlet = s.Matrix([[1-h*h*lam/2, h], [-h*lam*(1-h*h*lam/4), 1-h*h*lam/2]])
    R = s.diag(1, -1)
    require(midpoint.det() == verlet.det() == 1 and midpoint != verlet, 'distinct integrators missing')
    require(R*midpoint*R == midpoint.inv() and R*verlet*R == verlet.inv(), 'integrator reversal failed')
    return {'unequal_face_weights_keep_gauge_locality_but_break_covariance': True,
            'H_squared_keeps_gauge_covariance_but_breaks_face_locality': True,
            'gyroscopic_term_keeps_energy_Faraday_Gauss_but_breaks_time_reversal': True,
            'electric_magnetic_ratio_changes_generator': True,
            'distinct_reversible_symplectic_steps_at_h_half_lambda_two': {
                'implicit_midpoint': [[str(x) for x in row] for row in midpoint.tolist()],
                'velocity_Verlet': [[str(x) for x in row] for row in verlet.tolist()]}}


def build():
    D, C, _ = carrier()
    G, _, _ = evolution()
    paths = SOURCE_PATHS+('model/maxwell_local_selection/certificate.py',
                         'model/maxwell_local_selection/PROOF.md')
    return {'schema': 'oph.maxwell_local_selection.v1',
            'carrier': {'ports': 12, 'seams': 30, 'faces': 20,
                        'gradient_rank': D.rank(), 'curl_rank': C.rank(), 'rotations': 60},
            'quadratic_classification': classify(),
            'evolution': {'generator': 'dE/dt=(beta/epsilon) C^T B; dB/dt=-C E',
                          'energy': '(epsilon ||E||^2+beta ||B||^2)/2',
                          'electric_time_reversal_parity': 1, 'magnetic_time_reversal_parity': -1,
                          'generator_rank_exact': G.rank(),
                          'transverse_oscillator_count': 19,
                          'canonical_ambient_zero_modes': 12,
                          'unique_on_closed_flux_physical_carrier': True,
                          'uniqueness_inputs': ['linear autonomous real evolution', 'Faraday identity',
                                               'positive diagonal energy with supplied epsilon,beta',
                                               'energy conservation for every state',
                                               'specified electric/magnetic time reversal'],
                          'free_coefficients': ['electric energy epsilon>0', 'magnetic energy beta>0']},
            'controls': controls(),
            'scope': {'native_incidence_and_rotation_tables_used': True,
                      'local_quadratic_form_classified_exhaustively': True,
                      'Ampere_is_derived_from_stated_energy_and_reversal': True,
                      'native_record_law_selects_field_action_hypotheses': False,
                      'relative_energy_normalization_or_clock_selected': False,
                      'native_repair_implements_this_conservative_flow': False,
                      'continuum_or_photon_pole_identification': False,
                      'new_result_is_Lean_kernel_checked': False},
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    text = json.dumps(build(), sort_keys=True, indent=2)+'\n'
    out = HERE/'results.json'
    if args.check:
        require(out.read_text() == text, 'Maxwell selection receipt is stale')
        print('Exact local Maxwell classification and conditional evolution replay passed.')
    else:
        out.write_text(text)


if __name__ == '__main__':
    main()
