"""Bounded attachment of native source semantics to the L=2 regional input."""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import re

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SCHEMA = 'native-binary-source-refinement-v1'

from .python_epoch import PythonSourceEpoch

_PYTHON_EPOCH = PythonSourceEpoch(ROOT, __name__)
with _PYTHON_EPOCH.imports():
    from tools import strict_json
    from model.observer import native
    from model.native_quantum_reconstruction import certificate as reconstruction
    from model.native_gravity_response import certificate as gravity
    from . import native_validation as _validation_module
    from . import source_profile as _profile_module
    from .native_validation import NativeValidationSession, PREFIX, current_session
    from .source_profile import compile_native, grammar_cases, native_text, require, source_bindings


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def exact(matrix):
    return [[str(Q(x)) for x in row] for row in matrix]


def blocks(value):
    return {str(q): exact(matrix) for q, matrix in value.items()}


def carrier():
    states = sorted(tuple(int(k in (a, b)) for k in range(12))
                    for a in range(12) for b in range(a+1, 12))
    require(tuple(states) == tuple(native.states(binary=True)) == reconstruction.STATES,
            'source:carrier-order')
    seams = reconstruction.quantum_mixing.seams()
    require(len(states) == 66 and len(seams) == 30 and len(set(seams)) == 30
            and all(type(a) is int and type(b) is int and 0 <= a < b < 12 for a, b in seams),
            'source:seam-domain')
    return {'states': states, 'seams': seams}


def check_source_case(text, data, session):
    """Compare the live Lean evaluator with Python and the endpoint/record law."""
    function = compile_native(text)
    case = hashlib.sha256((text+'\n'+canonical(data)).encode('utf-8')).hexdigest()
    folder = session.output/'cases'/case
    folder.mkdir(parents=True, exist_ok=True)
    source, inputs = folder/'source.py', folder/'carrier.json'
    for path, body in ((source, (text+'\n').encode('utf-8')),
                       (inputs, (canonical(data)+'\n').encode('utf-8'))):
        if path.exists():
            require(path.read_bytes() == body, 'source:enrolled-case-input-drift')
        else:
            with path.open('xb') as stream:
                stream.write(body)
    raw = session.execute(PREFIX+'Replay', [source, inputs])
    replay = strict_json.loads(raw.decode('utf-8'))
    require(set(replay) == {'native_text', 'states', 'seams', 'counts', 'branches'}
            and replay['native_text'] == text
            and canonical(replay['states']) == canonical(data['states'])
            and canonical(replay['seams']) == canonical(data['seams']),
            'source:replay-input-binding')
    states, seams = data['states'], data['seams']
    index = {tuple(state): i for i, state in enumerate(states)}
    expected, counts = [], [[0]*66 for _ in states]
    for i, state in enumerate(states):
        for e, edge in enumerate(seams):
            for coin in (0, 1):
                final, record = function(tuple(state), edge, coin)
                require(final in index, 'source:binary-carrier-closure')
                a, b = edge
                endpoint_a = state[a] if state[a] == state[b] else coin
                endpoint_b = state[b] if state[a] == state[b] else 1-coin
                target = list(state)
                target[a], target[b] = endpoint_a, endpoint_b
                require(final == tuple(target), 'source:binary-endpoint-law')
                require(record == (coin, state[a]-final[a]), 'source:signed-record-law')
                j = index[final]
                expected.append([i, e, coin, j, *record])
                counts[i][j] += 1
    require(canonical(replay['branches']) == canonical(expected), 'source:full-branch-record-equality')
    require(canonical(replay['counts']) == canonical(counts), 'source:count-table-equality')
    return replay


def check_grammar_cases(text, session):
    cases = grammar_cases(text)
    body = canonical([[name, body] for name, body in cases.items()]).encode('utf-8')
    path = session.output/('grammar-'+hashlib.sha256(body).hexdigest()+'.json')
    if path.exists():
        require(path.read_bytes() == body, 'grammar:enrolled-input-drift')
    else:
        with path.open('xb') as stream:
            stream.write(body)
    actual = strict_json.loads(session.execute(PREFIX+'GrammarCheck', [path]).decode('utf-8'))
    expected = [{'name': name, 'source_text': body, 'grammar_admitted': False}
                for name, body in cases.items()]
    require(canonical(actual) == canonical(expected), 'grammar:case-input-output-binding')
    return actual


def _receipt_math(replay, data):
    counts = replay['counts']
    require(all(sum(row) == 60 for row in counts)
            and all(counts[a][b] == counts[b][a] for a in range(66) for b in range(66)),
            'count:normalization-or-symmetry')
    square = [[sum(counts[a][k]*counts[k][b] for k in range(66))
               for b in range(66)] for a in range(66)]
    total = sum(square[a][a] for a in range(66))
    require(total == 171720, 'count:return-normalizer')
    sigma = [[Q(value, total) for value in row] for row in square]
    actual_counts, actual_denominator = reconstruction.source_counts()
    require(counts == list(map(list, actual_counts)) and actual_denominator == 60,
            'consumer:actual-source-counts')
    C, Z, actual, floor = gravity.source_data()
    require(C.tolist() == counts and Z == total and actual.tolist() == sigma
            and list(map(list, reconstruction.reconstructed_density(2))) == sigma,
            'consumer:full-state-equality')
    lower = min(2*counts[a][a]-sum(counts[a]) for a in range(66))
    require(lower == 40 and floor == Q(lower**2, total), 'consumer:faithful-floor-value')
    paths = reconstruction.matrix_unit_paths()
    moments = []
    for a in range(66):
        for b in range(66):
            path, supplied_coefficient = paths[a, b]
            require(path[0] == a and path[-1] == b, 'moment:path-endpoints')
            coefficient = Q(1)
            for left, right in zip(path, path[1:]):
                coefficient *= Q(counts[left][right], 60)
            require(coefficient == supplied_coefficient and coefficient > 0,
                    'moment:source-word-coefficient')
            closed = coefficient*Q(square[b][a], 66*60**2)
            returned = Q(total, 66*60**2)
            moment = closed/(coefficient*returned)
            require(closed == reconstruction.closed_path_probability(2, path)
                    and returned == reconstruction.return_probability(2)
                    and moment == reconstruction.matrix_unit_moment(2, a, b) == sigma[b][a],
                    'moment:count-first-orientation')
            moments.append([a, b, list(path), str(coefficient), str(moment)])
    a, b = next((a, b) for (a, b), (path, _) in paths.items() if a != b and len(path) == 2)
    path, coefficient = paths[a, b]
    diagonal = np.diag(np.diag(actual))
    require(all(diagonal[i, i] >= 0 for i in range(66)) and gravity.trace(diagonal) == 1,
            'control:valid-diagonal-state')
    one_ports = []
    for port in range(12):
        reference = gravity.restrict(actual, (port,))
        require(blocks(reference) == blocks(gravity.restrict(diagonal, (port,))),
                'control:one-port-preservation')
        require(blocks(reference) == {'0': [['5/6']], '1': [['1/6']]},
                'control:one-port-values')
        one_ports.append(blocks(reference))
    full_plus = (sigma[a][a]+sigma[b][b]+sigma[a][b]+sigma[b][a])/2
    diagonal_plus = (diagonal[a, a]+diagonal[b, b])/2
    require(full_plus == Q(181, 11448) and diagonal_plus == Q(2611, 171720),
            'control:selected-effect-values')
    require(sigma[b][a] == Q(13, 21465) and full_plus-diagonal_plus == sigma[b][a]
            and diagonal[b, a] != sigma[b][a], 'control:off-diagonal-separation')
    omitted = reconstruction.closed_path_probability(2, path)/reconstruction.return_probability(2)
    require(omitted == Q(13, 1287900) and omitted != sigma[b][a], 'control:coefficient-omission')
    regional = gravity.restrict(actual, gravity.SEAM)
    return {'counts': counts, 'branches': replay['branches'], 'moments': moments,
            'normalizer': total, 'denominator': 60, 'density': exact(sigma),
            'faithful_floor': str(floor), 'region': list(gravity.SEAM),
            'regional_input': blocks(regional),
            'entropy': [[str(p), str(c)] for p, c in gravity.bulk_entropy(regional)],
            'one_port_states': one_ports,
            'scene': {'indices': [a, b], 'states': [list(data['states'][a]), list(data['states'][b])],
                      'path': list(path), 'coefficient': str(coefficient),
                      'tail_count': square[b][a],
                      'closed_probability': str(reconstruction.closed_path_probability(2, path)),
                      'return_probability': str(reconstruction.return_probability(2)),
                      'moment': str(sigma[b][a]), 'plus_effect': str(full_plus),
                      'diagonal_plus_effect': str(diagonal_plus),
                      'coefficient_omission': str(omitted)}}


@_PYTHON_EPOCH.protect
def produce(root=ROOT, *, session=None):
    root = Path(root).resolve()
    require(root == ROOT.resolve() and all(Path(module.__file__).resolve().is_relative_to(root)
            for module in (native, reconstruction, gravity)), 'source:loaded-module-root')
    before = source_bindings(root)
    text = native_text(root)
    endpoints = (root/'Lean/EventAlgebra/NativeBinarySourceEndpoints.lean').read_text(encoding='utf-8')
    literal = re.search(r'def sourceText : String := r#"(.*?)"#', endpoints, re.S)
    require(literal is not None and literal.group(1) == text, 'source:literal-slice-equality')
    session = session or current_session(root)
    require(session.root == root, 'kernel:session-root-binding')
    session.validate()
    data = carrier()
    replay = check_source_case(text, data, session)
    grammar = check_grammar_cases(text, session)
    result = {'schema': SCHEMA, 'source_bindings': before,
              'source_text': text, 'carrier': json.loads(canonical(data)),
              'result_class': 'bounded source/implementation evidence',
              'data': _receipt_math(replay, data)}
    result['data']['grammar_controls'] = grammar
    require(source_bindings(root) == before, 'source:changed-during-verification')
    return result


def verify(receipt, root=ROOT, *, session=None):
    require(type(receipt) is dict and set(receipt) == {
        'schema', 'source_bindings', 'source_text', 'carrier', 'result_class', 'data'}, 'receipt:schema')
    require(receipt['schema'] == SCHEMA, 'receipt:profile')
    require(receipt['source_bindings'] == source_bindings(root), 'receipt:source-custody')
    fresh = produce(root, session=session)
    compare_receipt(receipt, fresh)
    return fresh


def compare_receipt(receipt, receiver_result):
    """Compare every payload field with the receiver's complete fresh result.

    The caller supplies the result of produce(), never a receipt-provided
    authority. The public verification entry always performs that production.
    """
    require(canonical(receipt) == canonical(receiver_result), 'receipt:complete-semantic-equality')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    parser.add_argument('--evidence', type=Path, help='retain receiver-owned native evidence here')
    args = parser.parse_args()
    session = current_session(ROOT, args.evidence)
    if args.check:
        report = verify(strict_json.loads(args.check.read_text(encoding='utf-8')), session=session)
    else:
        report = produce(session=session)
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'schema': SCHEMA, 'branches': len(report['data']['branches']),
                      'moments': len(report['data']['moments']),
                      'off_diagonal_moment': report['data']['scene']['moment'],
                      'native_evidence': str(session.output)}))


_PYTHON_EPOCH.seal(source_bindings(ROOT), globals(),
    (strict_json, native, reconstruction, gravity, _validation_module, _profile_module))


if __name__ == '__main__':
    main()
