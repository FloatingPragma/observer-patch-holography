"""Independent Q(phi) and Choi-block replay; never imports the probe producer.

Source export arrays are reconstructed from the audited source function AST,
not obtained from NumPy or the compiler's transfer_kraus implementation.
"""
import argparse
import ast
import difflib
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess

import verify_marked_triplet_attachment_independent as m
from verify_galois_port_frame_independent import ZERO, ONE, add, sub, mul, div, scalar, total, matrix

ROOT = Path(__file__).resolve().parents[2]
BASELINE = 'a44330802c726037f88356f9f62044184ed429a5'
HISTORICAL_SOURCE_COMMIT = '53e365b248a4611bb911e75d7afae0f52bfd8675'
# Reviewed integration: detached Fraction returns only; see INTEGRATION_2026-10-29.md.
SOURCE_COMMIT = '83113a9cc33d82a56c84a8634445c8edd1e599a2'
REFERENCE = ROOT/'code/a5_closure/receipts/proper_code_attachment_probe.json'
COMPILER = 'model/observer/source_scalar/compiler.py'
INPUTS = sorted(set(m.INPUTS + [
    'code/a5_closure/manifests/marked_triplet_attachment_reference.json',
    'code/a5_closure/marked_triplet_attachment_certificate.py',
    'code/a5_closure/galois_port_frame_certificate.py',
    'code/a5_closure/verify_marked_triplet_attachment_independent.py',
    'code/a5_closure/verify_galois_port_frame_independent.py',
    'code/a5_closure/galois_source_response_control.py',
    'code/source_selection_model/DERIVATION.md', 'code/source_selection_model/geometry.py',
    'code/m1_source_realization/receipt.json', COMPILER,
    'model/native_clock/configuration.py', 'model/native_interface/tomography.py',
    'model/predictions/prepared_clock_source.py', 'model/predictions/prepared_clock.py',
    'model/predictions/inverse_clock.py',
]))
IMPLEMENTATIONS = ('code/a5_closure/proper_code_attachment_probe.py',
                   'code/a5_closure/verify_proper_code_attachment_independent.py')


def encoded(value):
    if isinstance(value, tuple) and len(value) == 2 and all(type(x) is Q for x in value):
        a, b = value
        return [str(a+b/2), str(b/2)]
    if isinstance(value, (list, tuple)):
        return [encoded(x) for x in value]
    if type(value) is dict:
        return {k: encoded(v) for k, v in value.items()}
    return value


def trans(a):
    return [list(x) for x in zip(*a)]


def exact_export_source(text):
    """Fail closed if the actual source service ceases to be this contract.

    This structural pin complements the Git/byte pins. Re-signing a changed
    source hash cannot turn a different service into the specified one.
    """
    expected = ast.parse('''
def transfer_kraus(d, coherent=True):
    if type(d) is not int or d not in (2, 4) or type(coherent) is not bool:
        raise ValueError("a declared d=2 or d=4 code and Boolean coherence required")
    kraus = []
    if coherent:
        success = np.zeros((d+1, 6), complex)
        success[:d, :d] = np.eye(d)
        kraus.append(success)
    else:
        for j in range(d):
            branch = np.zeros((d+1, 6), complex)
            branch[j, j] = 1
            kraus.append(branch)
    for j in range(d, 6):
        branch = np.zeros((d+1, 6), complex)
        branch[d, j] = 1
        kraus.append(branch)
    return kraus
''').body[0]
    functions = [n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef)
                 and n.name == 'transfer_kraus']
    m.need(len(functions) == 1, 'one source export service required')
    actual = functions[0]
    if isinstance(actual.body[0], ast.Expr) and isinstance(actual.body[0].value, ast.Constant):
        actual.body = actual.body[1:]
    m.need(ast.dump(actual) == ast.dump(expected), 'source export AST contract changed')


def service(d, coherent, upper, index):
    # Independent Choi construction from the source's index rules. Rows of
    # C are (input,output). Every outcome is a SUM of Kraus outer products.
    size = 6*(d+1)
    choi = {l: [[ZERO for _ in range(size)] for _ in range(size)] for l in ('success', 'failure')}
    columns = [list(range(d))] if coherent else [[j] for j in range(d)]
    for js in columns:
        for a in js:
            for b in js:
                choi['success'][a*(d+1)+a][b*(d+1)+b] = ONE
    for a in range(d, 6):
        choi['failure'][a*(d+1)+d][a*(d+1)+d] = ONE
    images = [{l: [[c[a*(d+1)+i][b*(d+1)+j] for j in range(d+1)] for i in range(d+1)]
               for l, c in choi.items()} for a in range(6) for b in range(6)]
    for a in range(6):
        for b in range(6):
            m.need(total(images[6*a+b][l][i][i] for l in choi for i in range(d+1))
                   == scalar(int(a == b)), 'Choi output partial trace must be I6')
    ident = all(images[6*a+b]['success'] ==
                [[scalar(int(i == a and j == b)) for j in range(d+1)] for i in range(d+1)]
                and all(x == ZERO for row in images[6*a+b]['failure'] for x in row)
                for a in range(3) for b in range(3))
    effect = [[total(images[6*j+i]['success'][r][r] for r in range(d+1))
               for j in range(3)] for i in range(3)]
    commutes = [m.mm(effect, u) == m.mm(u, effect) for u in upper]
    u = upper[index]
    probabilities = [total(mul(u[i][j], u[i][j]) for i in range(min(3, d))) for j in range(3)]
    column = next((j for j in range(3) if probabilities[j] != scalar(int(j < d))), None)
    witness = None
    if column is not None:
        before, after = scalar(int(column < d)), probabilities[column]
        witness = {'upper_basis_index': column,
                   'before': {'success': before, 'failure': sub(ONE, before)},
                   'after': {'success': after, 'failure': sub(ONE, after)}}
    # On full M6, the d=4 service also truncates part of the LOWER block.
    from verify_galois_port_frame_independent import sigma
    all_actions = [[[blocks[i//3][i%3][j%3] if i//3 == j//3 else ZERO
                     for j in range(6)] for i in range(6)]
                   for blocks in ((u, [[sigma(x) for x in row] for row in u]) for u in upper)]
    full_effect = [[scalar(int(i == j and i < d)) for j in range(6)] for i in range(6)]
    full_count = sum(m.mm(full_effect, u) == m.mm(u, full_effect) for u in all_actions)
    us = all_actions[index]
    full_probs = [total(mul(us[i][j], us[i][j]) for i in range(d)) for j in range(6)]
    full_column = next(j for j in range(6) if full_probs[j] != scalar(int(j < d)))
    # Restrict the Choi input to |0,0>+|1,1>, with half-weight density.
    bell = []
    for label in ('success', 'failure'):
        bell.append([[div(choi[label][a*(d+1)+i][b*(d+1)+j], scalar(2))
                      if a < 2 and b < 2 else ZERO
                      for j in range(d+1) for b in range(2)]
                     for i in range(d+1) for a in range(2)])
    return {'input_dimension': 6, 'output_dimension': d+1,
            'outcomes': ['success', 'failure'],
            'kraus_counts': {'success': 1 if coherent else d, 'failure': 6-d},
            'all_input_matrix_units_checked': 36, 'upper_matrix_units_checked': 9,
            'trace_preserving': True, 'upper_identity_including_failure': ident,
            'success_effect_upper': effect,
            'success_effect_commutes_with_all_marked_actions': all(commutes),
            'success_effect_commuting_action_count': sum(commutes),
            'full_M6_success_effect_commuting_action_count': full_count,
            'full_M6_fixed_mark_probability_witness': {'processor_basis_index': full_column,
                'before_success': scalar(int(full_column < d)), 'after_success': full_probs[full_column],
                'after_failure': sub(ONE, full_probs[full_column])},
            'fixed_mark_probability_witness': witness,
            'complete_matrix_unit_map_sha256': hashlib.sha256(json.dumps(
                encoded(images), sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
            'entangled_reference_fixture_output': bell,
            'reference_claim': 'all finite untouched references by complex matrix-unit linearity',
            'reference_fixture_dimension': 2,
            'fine_failure_kraus_index_is_extra_read_access': False}


def verify(packet, root=ROOT):
    pins = {}
    for path in INPUTS:
        data = (root/path).read_bytes()
        original = subprocess.check_output(['git', 'show', SOURCE_COMMIT+':'+path], cwd=ROOT)
        m.need(data == original, 'source custody mismatch: '+path)
        pins[path] = {'repository': 'FloatingPragma/observer-patch-holography',
                      'commit': SOURCE_COMMIT,
                      'blob': hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest(),
                      'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)}
    exact_export_source((root/COMPILER).read_text())
    frozen = m.load(root/'code/a5_closure/manifests/marked_triplet_attachment_reference.json')
    vertices = m.response_vertices()
    relabel = [vertices.index(v) for v in m.source_frame(root)]
    source_marker = m.source_marked_action(root, vertices, relabel)
    index = [tuple(a) for a in frozen['actions']].index(source_marker)
    m.need(index == 4, 'frozen source marker is not historical path 1')
    m.verify(frozen, root=root)
    t = trans([vertices[i] for i in relabel])
    ps = [matrix(p, 12, 12) for p in frozen['projectors']]
    scale = m.mm(t, trans(t))[0][0]
    m.need(m.mm(t, trans(t)) == [[mul(scale, x) for x in row] for row in m.eye(3)], 'frame surjectivity')
    m.need(m.mm(trans(t), t) == [[mul(scale, x) for x in row] for row in ps[0]], 'source slow Gram')
    m.need(m.mm(t, ps[0]) == t, 'frame factors through slow quotient')
    rows = m.load(root/'code/source_selection_model/response.json')['closed_paths']
    upper = [matrix([[cell[:2] for cell in r] for r in row['observed_response'][0]], 3, 3) for row in rows]
    lower = [matrix([[cell[:2] for cell in r] for r in row['observed_response'][1]], 3, 3) for row in rows]
    traces = lambda i: [total(u[j][j] for j in range(3)) for u in (upper[i], lower[i])]
    chars = [total(ps[0][i][p[i]] for i in range(12)) for p in frozen['actions']]
    hom = [div(total(mul(c, total(u[i][i] for i in range(3))) for c, u in zip(chars, us)), scalar(60))
           for us in (upper, lower)]
    m.need(hom == [ONE, ZERO], 'independent character Hom positive and negative controls')
    for p, u in zip(frozen['actions'], upper, strict=True):
        m.need([[r[p[j]] for j in range(12)] for r in t] == m.mm(u, t), 'source action covariance')
    tree = ast.parse((root/'code/a5_closure/galois_source_response_control.py').read_text())
    old = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
               and any(isinstance(x, ast.Name) and x.id == 'SOURCE_PINS' for x in n.targets))
    drift = {p: {'expected': digest, 'actual': hashlib.sha256(
        (root/'code/source_selection_model'/p).read_bytes()).hexdigest()} for p, digest in old.items()}
    mismatches = [p for p, row in drift.items() if row['expected'] != row['actual']]
    m.need(mismatches == ['DERIVATION.md'], 'historical drift scope changed; review required')
    old_commit = 'a84438d99a355dfd31e39aeb667e00d8ea60e445'
    old_text = subprocess.check_output(['git', 'show', old_commit+':code/source_selection_model/DERIVATION.md'], cwd=ROOT)
    m.need(hashlib.sha256(old_text).hexdigest() == old['DERIVATION.md'], 'historical derivation object')
    delta = ''.join(difflib.unified_diff(old_text.decode().splitlines(True),
        (root/'code/source_selection_model/DERIVATION.md').read_text().splitlines(True),
        fromfile=old_commit, tofile=HISTORICAL_SOURCE_COMMIT))
    expected = encoded({'schema': 'oph.proper-code-attachment-probe.v1',
        'work_baseline': BASELINE, 'source_pins': pins,
        'implementation_sha256': {p: hashlib.sha256((root/p).read_bytes()).hexdigest() for p in IMPLEMENTATIONS},
        'marking': {'selection': 'transport frozen generator before evaluating characters',
            'current_index': index, 'historical_index': 1, 'native_action': list(source_marker),
            'current_traces': traces(index), 'historical_traces': traces(1)},
        'mathematical_attachment': {'T_upper_from_spatial_loads': t, 'TT_star_scale': scale,
            'T_star_T_is_scale_times_slow_projector': True, 'actions_checked': 60,
            'inverse_on_slow': 'T_star / scale', 'hom_slow_to_upper': 1, 'hom_slow_to_lower': 0,
            'spatial_vectors_promoted_to_density_operators': False},
        'services': {f'd{d}_'+('coherent' if c else 'dephased'): service(d, c, upper, index)
                     for d in (2, 4) for c in (True, False)},
        'historical_galois_drift': {'pins': drift, 'mismatches': mismatches,
            'pinned_derivation_commit': old_commit,
            'derivation_diff_sha256': hashlib.sha256(delta.encode()).hexdigest(),
            'historical_full_replay_passed': False, 'hashes_updated': False},
        'scope': {'spatial_domain': 'real rank-three slow subspace of twelve-component load readouts',
            'internal_domain': 'upper complex three-dimensional subspace of declared M6 processor',
            'existing_protected_comparison_domain': 'declared B(C12) in a sixteen-address work register; native load-to-state identification not supplied by this export probe',
            'service_domain': 'declared M6 to C^(d+1), d in {2,4}, including failure',
            'coherent_proper_code_service_is_declared': True,
            'compiled_marked_action_words_certified': False,
            'source_spatial_quantum_preparation_encoding_derived': False,
            'complete_native_history_attachment_derived': False,
            'physical_golden_embedding_selected': False, 'target_character_used_for_selection': False,
            'universal_compiler_obstruction': False,
            'missing_premises': ['source identification of the spatial load quotient with an admitted quantum state and observable domain',
                'source operation encoding that domain into the internal upper processor',
                'frozen compiled marked-action words with complete service and error contracts']},
        'verdict': 'SOURCE_INTERFACE_UNDERDETERMINED'})
    # Serialization distinguishes integers, Booleans, strings and floats and
    # rejects additional fields. The parser rejects duplicates/nonfinite data.
    m.need(json.dumps(packet, sort_keys=True, allow_nan=False) ==
           json.dumps(expected, sort_keys=True, allow_nan=False), 'independent proper-code probe mismatch')
    return {'verified': True, 'verdict': expected['verdict'], 'actions_checked': 60,
            'service_cases': 4, 'matrix_units_per_service': 36,
            'method': 'Q(phi) frame/character reconstruction and complete Choi blocks'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', nargs='?', type=Path, default=REFERENCE)
    args = parser.parse_args()
    print(json.dumps(verify(m.load(args.receipt)), indent=2, sort_keys=True))
