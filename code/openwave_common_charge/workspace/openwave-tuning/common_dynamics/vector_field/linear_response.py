"""Transverse harmonic response of the declared local face complex.

This is the harmonic, unwrapped field response of the same plaquette kick.
Finite-register boundary and small-angle state errors remain explicit.
"""
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import hashlib
import json

import numpy as np
import sympy as s

import cell_complex as cc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def bloch_data():
    g, B, A, _, _, labels = cc.build(3)
    ico = tuple(cc.geometry.slow_geometry.source().seams())
    ico_index = {edge: i for i, edge in enumerate(ico)}
    edge_index = {edge: e for e, edge in enumerate(g.seams)}
    coords = g.cell_coordinates
    cells = {c: i for i, c in enumerate(coords)}
    metadata = {}
    for cell, c in enumerate(coords):
        for (a, b), local in ico_index.items():
            metadata[edge_index[(12 * cell + a, 12 * cell + b)]] = (c, local, 1)
        for axis, row in enumerate(cc.geometry.links()):
            nxt = list(c); nxt[axis] = (nxt[axis] + 1) % 3
            other = cells[tuple(nxt)]
            for channel, (a, b) in enumerate(row):
                x, y = 12 * cell + a, 12 * other + b
                metadata[edge_index[tuple(sorted((x, y)))]] = (c, 30 + 2 * axis + channel, 1 if x < y else -1)
    columns = [i for i, label in enumerate(labels) if label[1] == 0]
    assert len(columns) == 26 and len(metadata) == len(g.seams)
    terms = []
    for local_face, col in enumerate(columns):
        for edge in range(A.rows):
            if A[edge, col]:
                c, local_edge, sign = metadata[edge]
                shift = tuple(x if x <= 1 else x - 3 for x in c)
                terms.append((local_edge, local_face, shift, sign * int(A[edge, col])))
    A0 = s.zeros(36, 26)
    D = [s.zeros(36, 26) for _ in range(3)]
    for e, f, shift, sign in terms:
        A0[e, f] += sign
        for axis in range(3):
            D[axis][e, f] -= shift[axis] * sign
    B0 = s.zeros(12, 36)
    for (a, b), edge in ico_index.items():
        B0[a, edge], B0[b, edge] = 1, -1
    for axis, row in enumerate(cc.geometry.links()):
        for channel, (a, b) in enumerate(row):
            B0[a, 30 + 2 * axis + channel] = 1
            B0[b, 30 + 2 * axis + channel] = -1
    return A0, B0, D, terms, [labels[i] for i in columns]


def run():
    paths = [Path(__file__), Path(cc.__file__), Path(cc.geometry.__file__), HERE / 'CELL_COMPLEX.json']
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    A0, B0, D, terms, labels = bloch_data()
    assert A0.rank() == 22 and B0.rank() == 11 and B0 * A0 == s.zeros(12, 26)
    basis = s.Matrix.hstack(*(B0.col_join(A0.T)).nullspace())
    assert basis.shape == (36, 3)
    P = basis * (basis.T * basis).inv() * basis.T
    kernel = s.Matrix.hstack(*A0.nullspace())
    Q = kernel * (kernel.T * kernel).inv() * kernel.T
    assert P * P == P and Q * Q == Q
    macro = {tuple(label[2:]): i for i, label in enumerate(labels) if label[0] == 'macro_plaquette'}
    # D_b of the (a,b) plaquette is + the a-directed coarse chain;
    # D_a is its negative b-directed chain, with the chosen Bloch convention.
    w = [None] * 3
    for (a, b), f in macro.items():
        wa, wb = P * D[b][:, f], -P * D[a][:, f]
        if w[a] is None: w[a] = wa
        else: assert w[a] == wa
        if w[b] is None: w[b] = wb
        else: assert w[b] == wb
    W = s.Matrix.hstack(*w)
    gram = W.T * W
    assert gram == s.eye(3) * s.Rational(2, 3)
    assert W * gram.inv() * W.T == P
    momenta = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 2, 3)]
    for p in momenta:
        derivative = sum((p[a] * D[a] for a in range(3)), s.zeros(36, 26))
        effective = P * derivative * Q * derivative.T * P
        p = s.Matrix(p)
        expected = W * (s.eye(3) * p.dot(p) - p * p.T) * W.T
        assert effective == expected
    # Exact convolution row sums give a uniform bound before periodic
    # folding, so this bound covers every refinement rather than only L2.
    columns = defaultdict(list)
    for e, f, shift, coefficient in terms:
        columns[f].append((e, shift, coefficient))
    convolution = defaultdict(int)
    for column in columns.values():
        for e, shift, coefficient in column:
            for v, other_shift, other_coefficient in column:
                delta = tuple(a-b for a, b in zip(shift, other_shift))
                convolution[e, v, delta] += coefficient * other_coefficient
    row_sums = defaultdict(int)
    for (e, v, delta), coefficient in convolution.items():
        row_sums[e] += abs(coefficient)
    uniform_bound = max(row_sums.values())
    assert uniform_bound == 158
    # Continuous Bloch exactness checks reconstruct the incidence and face
    # symbols from actual size3 edges, without a continuum input matrix.
    numeric = []
    for k in ((.05, 0., 0.), (.025, 0., 0.), (.01, .02, .03)):
        A = np.zeros((36, 26), dtype=complex)
        for e, f, shift, sign in terms:
            A[e, f] += sign * np.exp(-1j * np.dot(k, shift))
        B = np.asarray(B0, dtype=complex)
        for axis, row in enumerate(cc.geometry.links()):
            for channel, (a, b) in enumerate(row):
                B[b, 30 + 2 * axis + channel] = -np.exp(-1j * k[axis])
        assert np.max(abs(B @ A)) < 2e-13
        eig = np.linalg.eigvalsh(A @ A.conj().T)
        positive = eig[eig > 1e-10]
        assert len(positive) == 24
        leading = 2 * sum(x*x for x in k) / 3
        numeric.append({'k': k, 'transverse_eigenvalues': positive[:2].tolist(),
                        'quadratic_prediction': leading,
                        'relative_errors': (abs(positive[:2] - leading) / leading).tolist()})
    # The one-step cosine-to-quadratic bound is a genuine interval bound,
    # conditional on the retained state's plaquette-angle support/leakage.
    bounds = [{'angle_cap': str(r), 'one_face_beta_quartic_error_for_beta_1_4': str(s.Rational(1, 4) * r**4 / 24)}
              for r in (s.Rational(1, 2), s.Rational(1, 4), s.Rational(1, 8))]
    assert pins == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result = {'schema': 'oph.local12.transverse_field_harmonic.v1',
              'checks': {'actual_size3_local_face_symbol': True, 'integer_BA_zero': True,
                         'exact_three_harmonic_edge_coordinates': True,
                         'exact_isotropic_transverse_quadratic_coefficient': True,
                         'two_nonzero_transverse_branches': True,
                         'finite_Bloch_response_controls': True,
                         'uniform_all_refinement_curvature_bound': True,
                         'cosine_remainder_bound_explicit': True, 'source_bytes_unchanged': True},
              'zero_symbol_face_rank': 22, 'zero_symbol_incidence_rank': 11,
              'harmonic_edge_rank': 3, 'coarse_cycle_projection_Gram': [['2/3', '0', '0'], ['0', '2/3', '0'], ['0', '0', '2/3']],
              'transverse_quadratic_symbol': '(2/3)(|k|^2 I3-k k^T) in an orthonormal harmonic frame',
              'linear_kick_equations': 'For V_field=U_E U_B and canonical angle defined by U_e=exp(i theta_e): E_next=E-beta A A* theta; theta_next=theta+alpha E_next. On each curvature mode lambda: cos(omega)=1-alpha*beta*lambda/2.',
              'linearized_speed_per_cell_step': 'sqrt(2 alpha beta/3)',
              'same_free_matter_speed_choice': {'alpha': '1/8', 'beta': '1/4', 'product': '1/32', 'speed': '1/sqrt48', 'ancestry': 'explicit dimensionless architecture choice, not a derived coupling or fine-structure constant'},
              'uniform_curvature_spectral_upper': uniform_bound,
              'Dirac_speed_subcycle_adapter': {'field_substeps_per_matter_macro': 8,
                  'alpha': '1/16', 'beta': '3/8', 'product': '3/128',
                  'harmonic_UV_stability_product_upper': '237/64', 'distance_below_four': '19/64',
                  'leading_speed_per_matter_macro': '1',
                  'scope': 'A declared speed/stability adapter for a separate c=1 matter amendment. It is not a transfer of the base scalar/fermion dynamics or its finite records, and not a quantum field continuum theorem.'},
              'numeric_Bloch_controls': numeric, 'one_step_quartic_bounds': bounds,
              'state_error_formula': 'For plaquette-angle support |theta_f|<=r_f with leakage probability p, the cosine and quadratic field kicks differ on the state by at most sum_f |beta_f| r_f^4/24+2 sqrt(p). Sum along actual intermediate states for multiple steps. Cyclic electric-wrap errors are additional.',
              'scope': 'The homogeneous local field extension has two transverse harmonic-response branches derived from its actual face complex. A quantum Maxwell continuum theorem additionally requires field-state preparation, retained leakage/flux-cutoff bounds, clock/precision refinement and interacting matter control. Finite modular shifts do not satisfy the exact unbounded U1 electric commutator.',
              'source_sha256': pins, 'official_rows_changed': []}
    (HERE / 'LINEAR_RESPONSE.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'checks': result['checks'], 'harmonic_Gram': result['coarse_cycle_projection_Gram'], 'controls': numeric}))


if __name__ == '__main__':
    run()
