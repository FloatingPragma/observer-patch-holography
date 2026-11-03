"""Finite source-routed phase access on the fixed GV code; no decoder added."""
import argparse
from functools import lru_cache
from fractions import Fraction as Q
import hashlib
from itertools import combinations
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import c12_current_attachment_probe as previous
from galois_port_frame_certificate import F, Z, O, encode, mm, transpose, require
from model.native_interface import tomography, pointer

BASELINE = '3800f801283078b6b86c23bbb9b85b545a70b69c'
# Reviewed current prerequisites, separate from the historical work baseline.
PRIOR_COMMIT = 'f5486275bea27c31ef7b7a02ae8212bec6a55db9'
SOURCE_COMMIT = previous.SOURCE_COMMIT
INPUTS = sorted(set(previous.INPUTS + [
    'model/compiled_native_interface/construction.py',
    'model/native_interface/construction.py',
    'model/native_quantum_history/construction.py',
    'model/observer/source_discrimination/circuit.py',
    'model/common_source/quantum_native.py', 'model/common_source/native_routing.py',
    'model/common_source/exact_output.py', 'model/common_source/QUANTUM_NATIVE.md',
    'model/native_interface/pointer.py', 'model/native_interface/README.md',
    'model/observer/source_discrimination/experiment.py',
]))
PRIOR_FILES = previous.PRIOR_FILES + previous.IMPLEMENTATIONS + (
    'code/a5_closure/receipts/c12_current_attachment_probe.json',)
IMPLEMENTATIONS = ('code/a5_closure/coherent_native_access_probe.py',
                  'code/a5_closure/verify_coherent_native_access_independent.py')
REFERENCE = ROOT/'code/a5_closure/receipts/coherent_native_access_probe.json'


def digest(value):
    return previous.digest(value)


def fixed_code():
    states = previous.native.states(binary=True)
    frame = previous.prior.marked.raw_frame()
    vertices = previous.prior.marked.response_vertices()
    mark = previous.prior.marked.source_marked_action(frame, vertices, [vertices.index(v) for v in frame])
    frozen = previous.prior.marked.build()
    require(frozen['actions'].index(list(mark)) == 4, 'marking before quadrature evaluation')
    b = [[sum((F(6*x[p]-1)*frame[p][j] for p in range(12)), Z) for j in range(3)] for x in states]
    k = F(3600, 720)
    require(mm(transpose(b), b) == [[k*int(i == j) for j in range(3)] for i in range(3)], 'fixed GV Gram')
    a, c = [r[0] for r in b], [r[1] for r in b]
    first = next((i, j) for i in range(66) for j in range(i+1, 66) if c[i]*a[j] != a[i]*c[j])
    require(first == (0, 6), 'source-frozen first nonzero imaginary coordinate')
    diagonal = [(x*x+y*y)/(2*k) for x, y in zip(a, c)]
    real = [[(a[i]*a[j]+c[i]*c[j])/(2*k) for j in range(66)] for i in range(66)]
    imag = [[(c[i]*a[j]-a[i]*c[j])/(2*k) for j in range(66)] for i in range(66)]
    return b, k, diagonal, real, imag, frozen


def active_pairs(operations):
    permutation = tomography.permutation(operations)
    require(len(set(map(int, permutation))) == 512, 'complete routing permutation')
    lookup = {int(permutation[c]): i for i, c in enumerate(tomography.codes())}
    return [(lookup.get(8*s+2), lookup.get(8*s+1)) for s in range(64)]


def probabilities(kind, pairs, diagonal, real, imag):
    """Ideal full 512-address target, then the source's record/R marginal.

    Pointer labels 0..65, timeout, export_abort; ideal detector has timeout 0.
    Padding is summed only into its declared classical abort event.
    """
    plus, minus = list(diagonal)+[Z, Z], list(diagonal)+[Z, Z]
    if kind == 'basis':
        return plus, minus
    for left, right in pairs:
        if left is None and right is None:
            continue
        p, q = diagonal[left] if left is not None else Z, diagonal[right] if right is not None else Z
        r = real[left][right] if left is not None and right is not None else Z
        c = imag[left][right] if left is not None and right is not None else Z
        for index, orientation in ((left, 1), (right, -1)):
            target = 67 if index is None else index
            base = (p+q)/2
            if index is not None:
                plus[target] = minus[target] = Z
            if kind == 'mixer':
                plus[target] += base+orientation*r
                minus[target] += base+orientation*r
            else:
                plus[target] += base-orientation*c
                minus[target] += base+orientation*c
    require(sum(plus, Z) == sum(minus, Z) == O, 'all ideal outcomes normalized')
    require(all(x.sign() >= 0 for x in plus+minus), 'complete ideal record positivity')
    return plus, minus


@lru_cache(maxsize=1)
def exact_probe():
    b, k, diagonal, real, imag, frozen = fixed_code()
    census = {'basis': 0, 'mixer': 0, 'y_analyzer': 0}
    records, route_rows = [], []
    selected = None
    minimum_gap = None
    # This is the existing tomography menu, fixed by pair labels and kind,
    # not a menu chosen by measured probabilities or a golden character.
    for left, right in combinations(range(66), 2):
        route = tomography.route(left, right)
        pairs = active_pairs(route)
        require((left, right) in pairs, 'source pair isolation')
        route_rows.append([left, right, [list(x) for x in route]])
        for kind in ('mixer', 'y_analyzer'):
            # Only imaginary cross terms contribute to this difference.
            delta = [Z]*68
            if kind == 'y_analyzer':
                for a, c in pairs:
                    if a is not None and c is not None:
                        delta[a] -= 2*imag[a][c]
                        delta[c] += 2*imag[a][c]
            census[kind] += int(any(x != Z for x in delta))
            if any(x != Z for x in delta):
                first = next(x for x in delta if x != Z)
                gap = first if first.sign() > 0 else -first
                if minimum_gap is None or (gap-minimum_gap).sign() < 0:
                    minimum_gap = gap
            require(sum(delta, Z) == Z, 'all record differences include abort and timeout')
            records.append([kind, left, right, delta])
            if (left, right, kind) == (0, 6, 'y_analyzer'):
                plus, minus = probabilities(kind, pairs, diagonal, real, imag)
                require([p-q for p, q in zip(plus, minus)] == delta, 'complete selected law')
                padding = next((a, c) for a, c in pairs if (a is None) != (c is None)
                               and any(x != Z for x in b[a if a is not None else c]))
                row = padding[0] if padding[0] is not None else padding[1]
                column = next(j for j, x in enumerate(b[row]) if x != Z)
                leakage = b[row][column]*b[row][column]/(2*k)
                selected = {'pair': [left, right], 'kind': kind,
                    'selection_rule': 'first nonzero imaginary entry in the previously frozen coordinate order',
                    'route': [list(x) for x in route], 'source_pulses': len(tomography.word(left, right, kind)),
                    'ideal_plus': plus, 'ideal_minus': minus,
                    'ideal_label_0_difference': delta[0],
                    'ideal_record_TV': sum((x if x.sign() >= 0 else -x for x in delta), Z)/2,
                    'padding_leakage_control': {'native_index': row, 'code_column': column,
                        'input_pair': list(padding), 'probability_lower': leakage},
                    'word_preserves_GV_code': False}
    require(selected is not None and selected['padding_leakage_control']['probability_lower'].sign() > 0,
            'selected unitary escapes the native GV code before destructive readout')
    # All 60 input actions; the readout is not assigned an invented output
    # triplet action. Its fixed score rotates nontrivially on the real code.
    score = [[(b[6][i]*b[0][j]-b[0][i]*b[6][j])/k for j in range(3)] for i in range(3)]
    rows = json.loads((ROOT/'code/source_selection_model/response.json').read_text())['closed_paths']
    states = previous.native.states(binary=True)
    lookup = {x: i for i, x in enumerate(states)}
    orbit = []
    for action, row in zip(frozen['actions'], rows, strict=True):
        u = [[F(Q(x[0]), Q(x[1])) for x in r] for r in row['observed_response'][0]]
        for i, x in enumerate(states):
            y = [0]*12
            for p, value in enumerate(x):
                y[action[p]] = value
            require(b[lookup[tuple(y)]] == [sum((u[p][q]*b[i][q] for q in range(3)), Z) for p in range(3)],
                    'all sixty fixed marked input actions')
        orbit.append(mm(mm(transpose(u), score), u))
    return {'menu': {'settings': 4291, 'pairs': 2145,
        'scope': 'basis record readout and one existing routed mixer or y_analyzer followed by tomography.records',
        'ideal_record_distinguishing_settings': census,
        'minimum_first_distinguishing_label_gap': minimum_gap, 'routes_sha256': digest(route_rows),
        'complete_ideal_record_differences_sha256': digest(records)}, 'selected': selected,
        'marked_action': {'current_index': 4, 'historical_index': 1, 'input_actions_checked': 60,
            'fixed_score_stabilizer_count': sum(x == score for x in orbit),
            'fixed_score_orbit_size': len({digest(x) for x in orbit}), 'score_orbit_sha256': digest(orbit),
            'output_triplet_action_supplied': False, 'action_twist_selected': False,
            'score_invariance_is_not_an_attachment_criterion': True}}


def service_bounds():
    certificate = tomography.certificate()  # Fresh original-binary source gates.
    require(certificate['maximum_certified_blocks'] == 31 and certificate['all_failures_retained'], 'native source certificate')
    bound = Q(certificate['complete_record_half_l1_upper'])
    require(bound < Q(1, 10**7), 'reference-safe complete record bound')
    gap = -exact_probe()['selected']['ideal_label_0_difference']
    require((gap-F(2*bound)).sign() > 0, 'finite record difference survives both model error bounds')
    require((exact_probe()['menu']['minimum_first_distinguishing_label_gap']-F(2*bound)).sign() > 0,
            'every distinguishing ideal setting survives both model errors')
    return {'primitive_operator_gate': '1/100000000000', 'maximum_certified_blocks': 31,
        'complete_record_half_l1_upper': str(bound), 'portable_upper': '1/10000000',
        'source_detector_timeout_upper': certificate['detector_timeout_upper'],
        'raw_state_propagation_used_for_witness': False,
        'physical_service_error_included': False, 'arbitrary_untouched_reference_bound': True,
        'unknown_input_code_preparation_is_supplied_not_derived': True,
        'all_original_binary_source_gates_replayed': True, 'gap_exceeds_two_model_error_bounds': True,
        'every_ideal_distinguishing_setting_is_model_robust': True,
        'exact_raw_equality_claimed_for_ideal_null_settings': False}


def pointer_probe():
    """The EXISTING auxiliary/native bit read, with surviving quantum work."""
    audit = pointer.certificate()
    require(audit['exact_instrument']['logical_superoperator_entries_checked'] == 64,
            'complete source pointer matrix-unit instrument')
    b, k, _, _, imag, _ = fixed_code()
    codes = tomography.codes()
    reports = []
    for bit in range(9):
        effect = [[sum((b[n][i]*b[n][j] for n in range(66) if codes[n] >> bit & 1), Z)/k
                   for j in range(3)] for i in range(3)]
        witness = next((i, j) for i in range(66) for j in range(i+1, 66)
                       if ((codes[i] ^ codes[j]) >> bit & 1) == 0 and imag[i][j] != Z)
        reports.append({'bit': bit, 'ideal_record_marginals_equal': True,
            'ideal_surviving_quantum_work_distinguishes': True,
            'retained_entry': list(witness), 'branch': (codes[witness[0]] >> bit) & 1,
            'plus_imaginary_coefficient': imag[witness[0]][witness[1]],
            'compressed_effect_sha256': digest(effect),
            'GV_code_invariant_under_bit_projection': mm(effect, effect) == effect})
    reset = tomography.configuration.budget(mode='reset', target=pointer.ready_state())
    error = reset['failure_upper']+audit['half_diamond_upper']
    chosen = reports[7]
    require(chosen['retained_entry'] == [0, 6] and (2*imag[0][6]-F(2*error)).sign() > 0,
            'complete surviving-work phase witness exceeds both pointer errors')
    return {'source_word_pulses': audit['source_pulses'], 'work_qubits': 9,
        'existing_native_auxiliary_ready_code': 2, 'logical_superoperator_entries_checked': 64,
        'input_auxiliary_independent_of_unknown_work_and_reference': True,
        'reset_timeout_included_in_complete_reduced_work_comparison': True,
        'finer_quantum_auxiliary_outputs_reconstructed_by_probe': False,
        'finer_quantum_auxiliary_outputs_claimed_equal': False,
        'reduced_work_witness_distinguishes_any_finer_retained_output': True,
        'ideal_success_work_maps': 'P_bit0 rho P_bit0, P_bit1 rho P_bit1',
        'records_only_equal_does_not_mean_quantum_outputs_equal': True,
        'bit_reports': reports, 'selected_bit': 7,
        'reset_timeout_upper': str(reset['failure_upper']),
        'reset_routing_diameter': tomography.configuration.routing_diameter(),
        'complete_reduced_work_error_upper': str(error),
        'raw_word_gate_passed': True, 'complete_quantum_output_gap_certified': True,
        'physical_service_error_included': False,
        'coherent_GV_decoder_or_internal_triplet_export_supplied': False}


def build():
    pins = previous.custody(INPUTS, SOURCE_COMMIT, 'FloatingPragma/observer-patch-holography')
    old = previous.custody(PRIOR_FILES, PRIOR_COMMIT, 'FloatingPragma/observer-patch-holography')
    packet = exact_probe()
    return encode({'schema': 'oph.coherent-native-access-probe.v1', 'work_baseline': BASELINE,
        'source_pins': pins, 'prior_probe_pins': old,
        'implementation_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in IMPLEMENTATIONS},
        **packet, 'service_bounds': service_bounds(), 'existing_native_work_pointer': pointer_probe(),
        'output_contract': {'labels': [str(i) for i in range(66)]+['timeout', 'export_abort'],
            'records_retains_untouched_reference': True, 'record_service_discards_final_native_work': True,
            'finer_controller_timeout_quantum_outputs_claimed_equal': False,
            'pre_readout_512_native_quantum_output_preserved': True,
            'unitary_target_preserves_orthogonal_inputs_and_references': True,
            'tomography_quantum_payload_exported': False,
            'tomography_uses_coherent_native_work_coupling': False,
            'tomography_uses_uncommitted_native_records': False},
        'verdict': 'COHERENT_ACCESS_WITHOUT_ATTACHMENT',
        'phase_0_stop_preserved': True, 'previous_record_only_grammar_not_overruled': True,
        'universal_coherent_grammar_obstruction_claimed': False,
        'physical_golden_embedding_selected': False,
        'conjugation': 'complex at a fixed real embedding; not Galois',
        'complete_operational_attachment_certified': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true')
    group.add_argument('--write', type=Path)
    args = parser.parse_args()
    text = json.dumps(build(), indent=2, sort_keys=True)+'\n'
    if args.check:
        require(REFERENCE.read_text() == text, 'coherent access receipt differs')
        print('COHERENT_ACCESS_WITHOUT_ATTACHMENT_VALID')
    elif args.write:
        with args.write.open('x') as stream:
            stream.write(text)
    else:
        print(text, end='')
