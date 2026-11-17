"""Exact local check: an equivalent Lindblad unravelling changes jump records.

The three coordinates are seam endpoint a, endpoint b, and an unaffected
position. All arithmetic below is in exactly represented dyadic rationals.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path

import numpy as np


def run(producer):
    eye = np.eye(3)
    a, b, other = [np.diag(eye[i]) for i in range(3)]
    ba, ab = np.outer(eye[1], eye[0]), np.outer(eye[0], eye[1])
    native = [
        (('C', 0, 1), ba), (('C', 0, 0), b+other),
        (('C', 1, -1), ab), (('C', 1, 0), a+other),
        (('E', 0, 1), ba), (('E', 0, -1), b), (('E', 0, 0), other),
        (('E', 1, 1), a), (('E', 1, -1), ab), (('E', 1, 0), other),
    ]  # actual K is one half times each displayed partial permutation
    channel = .25*sum(np.kron(k, k) for _, k in native)
    completeness = .25*sum(k.T@k for _, k in native)
    local = a+b
    recycling = [local, a, b, ba, ab]  # actual J is 1/sqrt(2) times each
    recycled = .5*sum(np.kron(j, j) for j in recycling)
    loss = .75*(np.kron(local, eye)+np.kron(eye, local))
    gamma = .5*sum(j.T@j for j in recycling)
    assert np.array_equal(completeness, eye)
    assert np.array_equal(channel-np.eye(9), recycled-loss)
    assert np.array_equal(gamma, 1.5*local)

    # The target helper is extracted solely to compare actual output matrices.
    # No graph construction, producer runner or acceptance case is executed.
    before = hashlib.sha256(producer.read_bytes()).hexdigest()
    tree = ast.parse(producer.read_text())
    selected = ast.Module(body=[node for node in tree.body if isinstance(node, ast.FunctionDef)
                               and node.name == 'kraus_bare'], type_ignores=[])
    namespace = {'np': np}
    exec(compile(selected, str(producer), 'exec'), namespace)
    actual = {}
    for label, source, target in namespace['kraus_bare'](2, [(0, 1)], 0):
        matrix = np.zeros((3, 3))
        matrix[target, source] = 1
        actual[label] = matrix
    assert set(actual) == {label for label, _ in native}
    assert all(np.array_equal(actual[label], k) for label, k in native)
    nonzero_record = .25*sum(np.linalg.norm(k@eye[0])**2 for label, k in native if label[2])
    moved_load = .25*sum(abs(k[1, 0])**2 for _, k in native)
    assert nonzero_record == .75 and moved_load == .5
    assert before == hashlib.sha256(producer.read_bytes()).hexdigest()
    return {'schema': 'oph.openwave.native_jump_record_scope.v1',
            'checks': {'native_completeness_exact': True,
                       'same_unconditional_generator_exact': True,
                       'actual_producer_Kraus_matrices_match': True,
                       'different_conditional_count_rates': True,
                       'source_unchanged': True},
            'degree_six_port_rates_per_unit_seam_rate': {
                'literal_nonzero_C_or_E_record': 6*nonzero_record,
                'literal_load_transfer': 6*moved_load,
                'recycling_unravelling_jump': 6*float(gamma[0, 0])},
            'all_scheduled_native_events_rate': 'M times the per-seam Poisson rate',
            'interpretation': 'H_eff=H-i(9r/2)P_port is the no-jump term of the recycling unravelling on a degree-six graph. Its norm is not the probability of no scheduled native event, no nonzero recorded value, or no load transfer under the declared C/E/MIX records.',
            'producer_sha256': before, 'producer': 'issues/N3D_cubic_routing/n3dlib.py',
            'auditor_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope': 'Exact one-seam superoperator identity extended by summing seams. No native full-state evolution is disputed. Conditional jump records require an instrument specification beyond the unconditional generator.',
            'official_rows_changed': []}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--producer', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.producer)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'checks': result['checks'], 'rates': result['degree_six_port_rates_per_unit_seam_rate']}))
