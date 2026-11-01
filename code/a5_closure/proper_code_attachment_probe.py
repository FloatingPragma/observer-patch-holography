"""Finite exact probe of the existing proper-code export, not a spatial encoding.

The only quantum inputs tested belong to the source's declared M6 processor.
No real spatial vector is admitted as a quantum state by this probe.
"""
import argparse
import ast
import difflib
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import marked_triplet_attachment_certificate as marked
from galois_port_frame_certificate import F, Z, O, encode, mm, transpose, require
from model.observer.source_scalar import compiler

BASELINE = 'a44330802c726037f88356f9f62044184ed429a5'
HISTORICAL_SOURCE_COMMIT = '53e365b248a4611bb911e75d7afae0f52bfd8675'
# Reviewed integration: detached Fraction returns only; see INTEGRATION_2026-10-29.md.
SOURCE_COMMIT = '31e1506c30114e37b45f6858d0029a8ff0b89db5'
REFERENCE = ROOT/'code/a5_closure/receipts/proper_code_attachment_probe.json'
COMPILER = 'model/observer/source_scalar/compiler.py'
INPUTS = sorted(set(marked.INPUTS + [
    'code/a5_closure/manifests/marked_triplet_attachment_reference.json',
    'code/a5_closure/marked_triplet_attachment_certificate.py',
    'code/a5_closure/galois_port_frame_certificate.py',
    'code/a5_closure/verify_marked_triplet_attachment_independent.py',
    'code/a5_closure/verify_galois_port_frame_independent.py',
    'code/a5_closure/galois_source_response_control.py',
    'code/source_selection_model/DERIVATION.md',
    'code/source_selection_model/geometry.py',
    'code/m1_source_realization/receipt.json', COMPILER,
    'model/native_clock/configuration.py', 'model/native_interface/tomography.py',
    'model/predictions/prepared_clock_source.py', 'model/predictions/prepared_clock.py',
    'model/predictions/inverse_clock.py',
]))
IMPLEMENTATIONS = ('code/a5_closure/proper_code_attachment_probe.py',
                   'code/a5_closure/verify_proper_code_attachment_independent.py')


def eye(n):
    return [[O if i == j else Z for j in range(n)] for i in range(n)]


def source_pins():
    pins = {}
    for path in INPUTS:
        data = (ROOT/path).read_bytes()
        original = subprocess.check_output(['git', 'show', SOURCE_COMMIT+':'+path], cwd=ROOT)
        require(data == original, 'scientific source differs from pinned source commit: '+path)
        blob = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        pins[path] = {'repository': 'FloatingPragma/observer-patch-holography',
                      'commit': SOURCE_COMMIT, 'blob': blob,
                      'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)}
    return pins


def export_maps(d, coherent):
    """Read actual source arrays; admit only their exact integral entries."""
    matrices = []
    for array in compiler.transfer_kraus(d, coherent=coherent):
        require(array.shape == (d+1, 6), 'source export dimensions')
        matrix = []
        for row in array:
            line = []
            for value in row:
                require(value.imag == 0 and value.real in (0, 1), 'exact source export entry')
                line.append(F(int(value.real)))
            matrix.append(line)
        matrices.append(matrix)
    return {'success': matrices[:1 if coherent else d],
            'failure': matrices[1 if coherent else d:]}


def image(maps, a, b):
    """Both complete CP branches on a complex-linear matrix-unit basis."""
    return {label: [[sum((k[i][a]*k[j][b] for k in ks), Z)
                     for j in range(len(ks[0]))] for i in range(len(ks[0]))]
            for label, ks in maps.items()}


def service_report(d, coherent, upper, marked_index):
    maps = export_maps(d, coherent)
    images = [image(maps, a, b) for a in range(6) for b in range(6)]
    for a in range(6):
        for b in range(6):
            out = images[6*a+b]
            require(sum((out[l][i][i] for l in out for i in range(d+1)), Z)
                    == F(int(a == b)), 'all-input trace preservation')
    upper_identity = True
    for a in range(3):
        for b in range(3):
            out = images[6*a+b]
            target = [[F(int(i == a and j == b)) for j in range(d+1)] for i in range(d+1)]
            upper_identity &= out['success'] == target and all(x == Z for r in out['failure'] for x in r)
    success_effect = [[sum((sum((k[r][i]*k[r][j] for r in range(d+1)), Z)
                           for k in maps['success']), Z) for j in range(3)] for i in range(3)]
    commuting = [mm(success_effect, u) == mm(u, success_effect) for u in upper]
    u = upper[marked_index]
    witness = next((j for j in range(3)
                    if sum((u[i][j]*u[i][j] for i in range(min(d, 3))), Z)
                    != F(int(j < d))), None)
    before = F(int(witness < d)) if witness is not None else None
    after = sum((u[i][witness]*u[i][witness] for i in range(min(d, 3))), Z) if witness is not None else None
    full_actions = [[[block[i//3][i%3][j%3] if i//3 == j//3 else Z
                      for j in range(6)] for i in range(6)]
                    for block in ((u, [[x.sigma() for x in r] for r in u]) for u in upper)]
    full_effect = [[F(int(i == j and i < d)) for j in range(6)] for i in range(6)]
    full_count = sum(mm(full_effect, u) == mm(u, full_effect) for u in full_actions)
    full_u = full_actions[marked_index]
    full_column = next(j for j in range(6) if sum((full_u[i][j]*full_u[i][j] for i in range(d)), Z)
                       != F(int(j < d)))
    full_after = sum((full_u[i][full_column]*full_u[i][full_column] for i in range(d)), Z)
    # A two-level untouched reference is a test input, not a new source register.
    bell = []
    for label in ('success', 'failure'):
        out = [[Z for _ in range(2*(d+1))] for _ in range(2*(d+1))]
        for a in (0, 1):
            for b in (0, 1):
                block = images[6*a+b][label]
                for i in range(d+1):
                    for j in range(d+1):
                        out[2*i+a][2*j+b] += block[i][j]/2
        bell.append(out)
    return {'input_dimension': 6, 'output_dimension': d+1,
            'outcomes': ['success', 'failure'],
            'kraus_counts': {k: len(v) for k, v in maps.items()},
            'all_input_matrix_units_checked': 36, 'upper_matrix_units_checked': 9,
            'trace_preserving': True, 'upper_identity_including_failure': bool(upper_identity),
            'success_effect_upper': success_effect,
            'success_effect_commutes_with_all_marked_actions': all(commuting),
            'success_effect_commuting_action_count': sum(commuting),
            'full_M6_success_effect_commuting_action_count': full_count,
            'full_M6_fixed_mark_probability_witness': {'processor_basis_index': full_column,
                'before_success': F(int(full_column < d)), 'after_success': full_after,
                'after_failure': O-full_after},
            'fixed_mark_probability_witness': None if witness is None else {
                'upper_basis_index': witness, 'before': {'success': before, 'failure': O-before},
                'after': {'success': after, 'failure': O-after}},
            'complete_matrix_unit_map_sha256': hashlib.sha256(
                json.dumps(encode(images), sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
            'entangled_reference_fixture_output': bell,
            'reference_claim': 'all finite untouched references by complex matrix-unit linearity',
            'reference_fixture_dimension': 2,
            'fine_failure_kraus_index_is_extra_read_access': False}


def drift_report():
    tree = ast.parse((ROOT/'code/a5_closure/galois_source_response_control.py').read_text())
    pins = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == 'SOURCE_PINS' for t in n.targets))
    rows = {p: {'expected': h, 'actual': hashlib.sha256(
        (ROOT/'code/source_selection_model'/p).read_bytes()).hexdigest()} for p, h in pins.items()}
    old_commit = 'a84438d99a355dfd31e39aeb667e00d8ea60e445'
    old = subprocess.check_output(['git', 'show', old_commit+':code/source_selection_model/DERIVATION.md'], cwd=ROOT)
    require(hashlib.sha256(old).hexdigest() == pins['DERIVATION.md'], 'historical derivation custody')
    delta = ''.join(difflib.unified_diff(old.decode().splitlines(True),
        (ROOT/'code/source_selection_model/DERIVATION.md').read_text().splitlines(True),
        fromfile=old_commit, tofile=HISTORICAL_SOURCE_COMMIT))
    return {'pins': rows, 'mismatches': [p for p, r in rows.items() if r['actual'] != r['expected']],
            'pinned_derivation_commit': old_commit,
            'derivation_diff_sha256': hashlib.sha256(delta.encode()).hexdigest(),
            'historical_full_replay_passed': False, 'hashes_updated': False}


def build():
    pins = source_pins()
    frozen = json.loads((ROOT/'code/a5_closure/manifests/marked_triplet_attachment_reference.json').read_text())
    frame, vertices = marked.raw_frame(), marked.response_vertices()
    relabel = [vertices.index(v) for v in frame]
    mark = marked.source_marked_action(frame, vertices, relabel)
    index = frozen['actions'].index(list(mark))
    require(index == 4 and index != 1, 'source marker precedes character evaluation')
    live = marked.build()
    require(live == frozen, 'fresh fixed marked action must match its frozen receipt')
    t = transpose(frame)
    scale = mm(t, transpose(t))[0][0]
    p = [[F(Q(a), Q(b)) for a, b in row] for row in frozen['projectors'][0]]
    require(mm(t, transpose(t)) == [[scale*x for x in r] for r in eye(3)], 'frame full rank')
    require(mm(transpose(t), t) == [[scale*x for x in r] for r in p], 'frame slow Gram')
    response = json.loads((ROOT/marked.RESPONSE).read_text())
    upper = [[[F(Q(x[0]), Q(x[1])) for x in row] for row in path['observed_response'][0]]
             for path in response['closed_paths']]
    for action, u in zip(frozen['actions'], upper, strict=True):
        require([[r[action[j]] for j in range(12)] for r in t] == mm(u, t), 'full fixed-action covariance')
    services = {f'd{d}_'+('coherent' if c else 'dephased'): service_report(d, c, upper, index)
                for d in (2, 4) for c in (True, False)}
    require(services['d4_coherent']['upper_identity_including_failure'], 'positive internal export control')
    require(not services['d2_coherent']['success_effect_commutes_with_all_marked_actions'], 'd2 truncation control')
    require(not services['d4_dephased']['upper_identity_including_failure'], 'dephasing is not identity transfer')
    traces = lambda i: [sum((u[j][j] for j in range(3)), Z) for u in (
        upper[i], [[x.sigma() for x in row] for row in upper[i]])]
    return encode({'schema': 'oph.proper-code-attachment-probe.v1',
        'work_baseline': BASELINE, 'source_pins': pins,
        'implementation_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in IMPLEMENTATIONS},
        'marking': {'selection': 'transport frozen generator before evaluating characters',
            'current_index': index, 'historical_index': 1, 'native_action': list(mark),
            'current_traces': traces(index), 'historical_traces': traces(1)},
        'mathematical_attachment': {'T_upper_from_spatial_loads': t, 'TT_star_scale': scale,
            'T_star_T_is_scale_times_slow_projector': True, 'actions_checked': len(upper),
            'inverse_on_slow': 'T_star / scale', 'hom_slow_to_upper': 1, 'hom_slow_to_lower': 0,
            'spatial_vectors_promoted_to_density_operators': False},
        'services': services, 'historical_galois_drift': drift_report(),
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true')
    group.add_argument('--write', type=Path)
    args = parser.parse_args()
    payload = json.dumps(build(), indent=2, sort_keys=True)+'\n'
    if args.check:
        require(REFERENCE.read_text() == payload, 'proper-code probe receipt differs')
        print('PROPER_CODE_ATTACHMENT_PROBE_VALID: SOURCE_INTERFACE_UNDERDETERMINED')
    elif args.write:
        with args.write.open('x') as stream:
            stream.write(payload)
    else:
        print(payload, end='')
