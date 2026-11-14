"""Check the graph-independent native pair record second-moment identity.

The proof uses a product-graph Poisson equation and a cycle martingale.
Small graph traversal identities are verified in rational arithmetic. The
carrier projector is compared with Fable's independently saved moment solve.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import sympy as sp


def exact_case(vertices, edges):
    count = len(edges)
    incidence = sp.zeros(vertices, count)
    for j, (a, b) in enumerate(edges):
        incidence[a, j], incidence[b, j] = 1, -1
    laplacian = incidence * incidence.T
    pairs = [(x, y) for x in range(vertices) for y in range(vertices) if x != y]
    index = {state: i for i, state in enumerate(pairs)}
    ids = [x * vertices + y for x, y in pairs]
    product = sp.kronecker_product(laplacian, sp.eye(vertices)) + sp.kronecker_product(sp.eye(vertices), laplacian)
    killed = product.extract(ids, ids)
    # Multiplication by 2M converts each event expectation to this generator.
    rhs = sp.Matrix([[int(x in edge) + int(y in edge) for edge in edges] for x, y in pairs])
    inverse = killed.inv()
    expected = inverse * rhs
    lifetime = inverse * (2 * count * sp.ones(len(pairs), 1))
    initial = sp.zeros(1, len(pairs))
    for a, b in edges:
        initial[0, index[a, b]] += sp.Rational(1, 2 * count)
        initial[0, index[b, a]] += sp.Rational(1, 2 * count)
    per_edge = list(initial * expected)
    mean_lifetime = (initial * lifetime)[0]
    reduced = incidence[:-1, :]
    projector = sp.eye(count) - reduced.T * (reduced * reduced.T).inv() * reduced
    predicted = sp.Rational(vertices, count) * projector
    alternate = list(expected[index[edges[0]], :])
    return {'vertices': vertices, 'edges': count,
            'per_edge_transfer_means': list(map(str, per_edge)),
            'adjacent_mean_lifetime': str(mean_lifetime),
            'predicted_record_squared_norm': str(sp.trace(predicted)),
            'expected_record_second_moment': [[str(v) for v in predicted.row(i)] for i in range(count)],
            'every_edge_mean_exact': all(v == sp.Rational(vertices - 1, count) for v in per_edge),
            'lifetime_exact': mean_lifetime == sp.Rational(vertices * (vertices - 1), 2),
            'projector_exact': projector * projector == projector and incidence * projector == sp.zeros(vertices, count),
            'fixed_first_edge_changes_traversal_vector': alternate != per_edge}


def saved_carrier_comparison(fable):
    # The saved covariance uses Fable's seam ordering, so parse its literal
    # carrier constants and verify the graph rather than guess an isomorphism.
    import ast
    geometry = fable.parent / 'reverse-engineering-reality/evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'
    matrix_file = fable / 'issues/LA_L02/carrier_record_second_moment.npy'
    sources = [geometry, matrix_file, fable / 'kernel/ophk.py', fable / 'issues/LA_L02/audit.py']
    before = {os.path.relpath(p, fable): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    tree = ast.parse(geometry.read_text())
    constants = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    try:
                        constants[target.id] = ast.literal_eval(node.value)
                    except (ValueError, TypeError):
                        pass
    # Exported geometry is compared as data; no Fable code is executed.
    faces = constants['ORIENTED_BASE_FACES']
    edges = sorted({tuple(sorted((face[k], face[(k + 1) % 3]))) for face in faces for k in range(3)})
    incidence = np.zeros((12, len(edges)))
    for j, (a, b) in enumerate(edges):
        incidence[a, j], incidence[b, j] = 1, -1
    assert len(edges) == 30 and np.all(np.sum(abs(incidence), axis=1) == 5)
    reduced = incidence[:-1]
    projector = np.eye(30) - reduced.T @ np.linalg.solve(reduced @ reduced.T, reduced)
    observed = np.load(matrix_file, allow_pickle=False)
    expected = 12 / 30 * projector
    assert observed.shape == expected.shape
    return {'matrix_max_error': float(np.max(abs(observed - expected))),
            'trace_observed': float(np.trace(observed)), 'trace_predicted': float(np.trace(expected)),
            'expected_trace_exact': '38/5',
            'source_sha256': before,
            'sources_unchanged': before == {os.path.relpath(p, fable): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}


def run(fable):
    examples = {'path4': (4, [(0, 1), (1, 2), (2, 3)]),
                'star5': (5, [(0, i) for i in range(1, 5)]),
                'triangle_tail': (4, [(0, 1), (1, 2), (0, 2), (2, 3)]),
                'diamond': (4, [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]),
                'house5': (5, [(0, 1), (1, 2), (2, 3), (0, 3), (1, 4), (2, 4)]),
                'K4': (4, [(i, j) for i in range(4) for j in range(i + 1, 4)])}
    cases = {name: exact_case(n, edges) for name, (n, edges) in examples.items()}
    comparison = saved_carrier_comparison(fable)
    checks = {'exact_per_edge_counts_all_six_graphs': all(x['every_edge_mean_exact'] for x in cases.values()),
              'exact_lifetimes_all_six_graphs': all(x['lifetime_exact'] for x in cases.values()),
              'exact_cycle_projectors_all_six_graphs': all(x['projector_exact'] for x in cases.values()),
              'actual_carrier_moment_matrix_matches': comparison['matrix_max_error'] < 1e-11,
              'comparison_sources_unchanged': comparison['sources_unchanged'],
              'arbitrary_preparation_transfer_rejected': cases['triangle_tail']['fixed_first_edge_changes_traversal_vector']}
    return {'schema': 'oph.audit.native_pair_covariance.v1', 'checks': checks, 'exact_cases': cases,
            'saved_carrier_comparison': comparison,
            'theorem': 'For a uniformly chosen oriented initial edge with one excess, one deficit and its compensating integer C record, E[F_tau F_tau^T]=(N/M)P_cycle on every finite connected simple graph. Mean nonzero transfers on each edge are (N-1)/M.',
            'scope': 'Uniform independent seam attempts, fair coins, classical integer current records and absorption of one signed pair only. No coherent Hamiltonian, modular representative energy, physical photon or in-platform row award.',
            'reviewer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fable', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = run(args.fable)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'checks': result['checks'], 'saved_carrier_comparison': result['saved_carrier_comparison']}, indent=2))
    raise SystemExit(0 if all(result['checks'].values()) else 1)
