"""Independent GF(2) routes, Q(phi) effects and raw-block interval bounds.

Imports no new producer, tomography executor or compiled-interface checker.
The existing source compiler is used only to replay small original-binary
blocks; all targets, norm enclosures, routes and code effects are rebuilt.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
import sys

import verify_c12_current_attachment_independent as previous
from verify_galois_port_frame_independent import ZERO as Z, ONE as O, add, sub, mul, div, scalar, total, sign

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from model.observer.source_scalar import compiler

BASELINE = '3800f801283078b6b86c23bbb9b85b545a70b69c'
# Reviewed current prerequisites, separate from the historical work baseline.
PRIOR_COMMIT = '2da8c4b09dca49b622b9b5454ebc1a4e91e7cc3c'
SOURCE_COMMIT = previous.SOURCE_COMMIT
INPUTS = sorted(set(previous.INPUTS + [
    'model/compiled_native_interface/construction.py',
    'model/native_interface/construction.py', 'model/native_quantum_history/construction.py',
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
m = previous.m


def digest(value):
    return previous.digest(value)


def source_codes():
    # Reconstruct the source's occupation-shell address convention, including
    # the preferred shells; no native encoder is called.
    states = previous.configurations(full=True)
    shells = [sorted({x[2:] for x in states if sum(x[2:]) == n}) for n in range(3)]
    for n, preferred in ((1, (1,)+(0,)*9), (2, (1, 1)+(0,)*8)):
        shells[n].remove(preferred); shells[n].insert(0, preferred)
    local = ((0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0))
    return [local.index(x[:2])+8*shells[sum(x[2:])].index(x[2:]) for x in previous.configurations()]


def affine(value, operations):
    for op in operations:
        if op[0] == 'X':
            value ^= 1 << op[1]
        else:
            value ^= ((value >> op[1]) & 1) << op[2]
    return value


def route(codes, a, b):
    difference = codes[a] ^ codes[b]
    pivot = min(j for j in range(9) if difference >> j & 1)
    operations = [('X', j) for j in range(9) if codes[a] >> j & 1]
    operations += [('C', pivot, j) for j in range(9) if j != pivot and difference >> j & 1]
    if pivot != 0:
        operations += [('C', pivot, 0), ('C', 0, pivot), ('C', pivot, 0)]
    operations += [('C', 0, 1), ('X', 1)]
    m.need((affine(codes[a], operations), affine(codes[b], operations)) == (2, 1), 'GF(2) native route')
    return operations


def interval_block_gate(array, target):
    """Bound ||A-U||F² with rational component intervals for sqrt(2).

    Target entries are (r,s,t,u) for r+s*sqrt(2)+i(t+u*sqrt(2)).
    This differs from the source checker's radical-sum cross-term norm.
    """
    scale = 10**80
    lo = Q(math.isqrt(2*scale*scale), scale)
    hi = lo+Q(1, scale)
    m.need(lo*lo <= 2 <= hi*hi, 'sqrt(2) interval')
    squared = Q()
    for i in range(len(array)):
        for j in range(len(array)):
            r, s, t, u = target.get((i, j), (Q(),)*4)
            for actual, p, q in ((array[i, j].real, r, s), (array[i, j].imag, t, u)):
                x = Q.from_float(float(actual))
                endpoints = (x-p-q*lo, x-p-q*hi)
                squared += max(abs(e) for e in endpoints)**2
    m.need(squared <= Q(1, 10**22), 'independent original-binary source block gate')
    return squared


@lru_cache(maxsize=1)
def raw_source_gates():
    # Rebuild the seven actual block schedules and full targets from source
    # Walsh/parity rules, not by importing tomography or compiled.target.
    pulse = compiler.pulse
    xword = [pulse('I', -math.pi/2, [0]), pulse('X', math.pi/2, [0])]
    cxword = compiler.cnot(0, 1)
    reflection = [pulse('I', -math.pi/8, [0])]
    for pivot in range(3):
        for j in range(1 << pivot):
            mask = (1 << pivot) | (j ^ (j >> 1))
            reflection.append(pulse('Z', -(-1)**((3 & mask).bit_count())*math.pi/8, [pivot]))
            following = (j+1) % (1 << pivot)
            change = (j ^ (j >> 1)) ^ (following ^ (following >> 1))
            if change:
                reflection += compiler.cnot(change.bit_length()-1, pivot)
    balanced = (compiler.cnot(0, 1)+[pulse('Y', math.pi/8, [0])]+reflection
                +[pulse('Y', -math.pi/8, [0])]+compiler.cnot(0, 1))
    def inverse(word):
        return [pulse(x['kind'], -x['angle'], x['targets']) for x in reversed(word)]
    integer = lambda n, p: {(p[j], j): (Q(1), Q(), Q(), Q()) for j in range(n)}
    mixer = {(j, j): (Q(1), Q(), Q(), Q()) for j in range(8) if j not in (1, 2)}
    y = {(j, j): (Q(), Q(1, 2), Q(), Q((-1)**(j & 1), 2)) for j in range(8) if j not in (1, 2)}
    for i, j, s in ((2, 2, 1), (2, 1, 1), (1, 2, 1), (1, 1, -1)):
        mixer[i, j] = Q(), Q(s, 2), Q(), Q()
        y[i, j] = Q(s, 2), Q(), Q(s*(-1)**(j & 1), 2), Q()
    blocks = [(xword, 1, integer(2, [1, 0])), (cxword, 2, integer(4, [0, 3, 2, 1])),
        (inverse(xword), 1, integer(2, [1, 0])), (inverse(cxword), 2, integer(4, [0, 3, 2, 1])),
        (balanced, 3, mixer), ([pulse('Z', math.pi/4, [0])]+balanced, 3, y)]
    permutation = [affine(x, [('C', 2, 0), ('C', 0, 1)]) for x in range(8)]
    blocks.append((compiler.cnot(2, 0)+compiler.cnot(0, 1), 3, integer(8, permutation)))
    for word, bits, target in blocks:
        interval_block_gate(compiler.execute(word, bits), target)
    return True


@lru_cache(maxsize=1)
def mathematical_packet():
    frame = m.source_frame(ROOT)
    vertices = m.response_vertices()
    mark = m.source_marked_action(ROOT, vertices, [vertices.index(v) for v in frame])
    frozen = m.load(ROOT/'code/a5_closure/manifests/marked_triplet_attachment_reference.json')
    m.need(frozen['actions'].index(list(mark)) == 4, 'mark fixed before phase comparison')
    states = previous.configurations()
    b = [[total(mul(scalar(6*x[p]-1), frame[p][j]) for p in range(12)) for j in range(3)] for x in states]
    k = (Q(2880), Q(1440))
    m.need(m.mm(previous.prior.trans(b), b) == [[k if i == j else Z for j in range(3)] for i in range(3)], 'GV Gram')
    a, c = [r[0] for r in b], [r[1] for r in b]
    real = [[div(add(mul(a[i], a[j]), mul(c[i], c[j])), mul(scalar(2), k)) for j in range(66)] for i in range(66)]
    imag = [[div(sub(mul(c[i], a[j]), mul(a[i], c[j])), mul(scalar(2), k)) for j in range(66)] for i in range(66)]
    m.need(next((i, j) for i in range(66) for j in range(i+1, 66) if imag[i][j] != Z) == (0, 6), 'previous first coordinate')
    diagonal = [real[i][i] for i in range(66)]
    codes = source_codes()
    census = {'basis': 0, 'mixer': 0, 'y_analyzer': 0}
    records, route_rows, maximum = [], [], 0
    selected = None
    minimum_gap = None
    for left, right in combinations(range(66), 2):
        operations = route(codes, left, right)
        maximum = max(maximum, 2*len(operations)+1)
        route_rows.append([left, right, [list(x) for x in operations]])
        lookup = {affine(code, operations): i for i, code in enumerate(codes)}
        pairs = [(lookup.get(8*s+2), lookup.get(8*s+1)) for s in range(64)]
        delta = [Z]*68
        for i, j in pairs:
            if i is not None and j is not None:
                delta[i] = sub(delta[i], mul(scalar(2), imag[i][j]))
                delta[j] = add(delta[j], mul(scalar(2), imag[i][j]))
        records += [['mixer', left, right, [Z]*68], ['y_analyzer', left, right, delta]]
        census['y_analyzer'] += int(any(x != Z for x in delta))
        if any(x != Z for x in delta):
            first = next(x for x in delta if x != Z)
            gap = first if sign(first) > 0 else tuple(-v for v in first)
            if minimum_gap is None or sign(sub(gap, minimum_gap)) < 0:
                minimum_gap = gap
        if (left, right) == (0, 6):
            plus, minus = diagonal+[Z, Z], diagonal+[Z, Z]
            for i, j in pairs:
                if i is None and j is None:
                    continue
                mean = div(add(diagonal[i] if i is not None else Z, diagonal[j] if j is not None else Z), scalar(2))
                ci = imag[i][j] if i is not None and j is not None else Z
                for index, s in ((i, 1), (j, -1)):
                    at = index if index is not None else 67
                    if index is not None:
                        plus[at] = minus[at] = Z
                    plus[at] = add(plus[at], sub(mean, mul(scalar(s), ci)))
                    minus[at] = add(minus[at], add(mean, mul(scalar(s), ci)))
            m.need(total(plus) == total(minus) == O and all(sign(x) >= 0 for x in plus+minus), 'complete ideal positive law')
            padding = next((i, j) for i, j in pairs if (i is None) != (j is None)
                           and any(x != Z for x in b[i if i is not None else j]))
            row = padding[0] if padding[0] is not None else padding[1]
            col = next(j for j, x in enumerate(b[row]) if x != Z)
            leakage = div(mul(b[row][col], b[row][col]), mul(scalar(2), k))
            # Each routed X uses two pulses, each CNOT three.
            pulses = 2*sum(2 if op[0] == 'X' else 3 for op in operations)+35
            selected = {'pair': [0, 6], 'kind': 'y_analyzer',
                'selection_rule': 'first nonzero imaginary entry in the previously frozen coordinate order',
                'route': [list(x) for x in operations], 'source_pulses': pulses,
                'ideal_plus': plus, 'ideal_minus': minus, 'ideal_label_0_difference': delta[0],
                'ideal_record_TV': div(total(x if sign(x) >= 0 else tuple(-v for v in x) for x in delta), scalar(2)),
                'padding_leakage_control': {'native_index': row, 'code_column': col,
                    'input_pair': list(padding), 'probability_lower': leakage}, 'word_preserves_GV_code': False}
    m.need(maximum == 31, 'all-route block envelope')
    score = [[div(sub(mul(b[6][i], b[0][j]), mul(b[0][i], b[6][j])), k) for j in range(3)] for i in range(3)]
    source = m.load(ROOT/'code/source_selection_model/response.json')['closed_paths']
    orbit = []
    lookup = {x: i for i, x in enumerate(states)}
    for action, row in zip(frozen['actions'], source, strict=True):
        u = [[(Q(x[0])-Q(x[1]), 2*Q(x[1])) for x in line] for line in row['observed_response'][0]]
        inverse = [action.index(i) for i in range(12)]
        for i, x in enumerate(states):
            image = lookup[tuple(x[inverse[p]] for p in range(12))]
            m.need(b[image] == [total(mul(u[p][q], b[i][q]) for q in range(3)) for p in range(3)], 'sixty marked code actions')
        orbit.append(m.mm(m.mm(previous.prior.trans(u), score), u))
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


@lru_cache(maxsize=1)
def pointer_packet():
    states = previous.configurations()
    frame = m.source_frame(ROOT)
    b = [[total(mul(scalar(6*x[p]-1), frame[p][j]) for p in range(12)) for j in range(3)] for x in states]
    k, codes = (Q(2880), Q(1440)), source_codes()
    imag = lambda i, j: div(sub(mul(b[i][1], b[j][0]), mul(b[i][0], b[j][1])), mul(scalar(2), k))
    reports = []
    for bit in range(9):
        effect = [[div(total(mul(b[n][i], b[n][j]) for n in range(66) if codes[n] >> bit & 1), k)
                   for j in range(3)] for i in range(3)]
        pair = next((i, j) for i in range(66) for j in range(i+1, 66)
                    if ((codes[i] ^ codes[j]) >> bit & 1) == 0 and imag(i, j) != Z)
        reports.append({'bit': bit, 'ideal_record_marginals_equal': True,
            'ideal_surviving_quantum_work_distinguishes': True, 'retained_entry': list(pair),
            'branch': (codes[pair[0]] >> bit) & 1, 'plus_imaginary_coefficient': imag(*pair),
            'compressed_effect_sha256': digest(effect),
            'GV_code_invariant_under_bit_projection': m.mm(effect, effect) == effect})
    # Independently reconstruct both native effects on the four actual aux
    # codes, then all 64 entries of the reduced CP instrument.
    seam = ((0, 0), (0, 1), (1, 0), (0, 2))
    for law in ('C', 'E'):
        effects = [[0]*4 for _ in range(2)]
        groups = {}
        for code, (a, c) in enumerate(seam):
            for coin in (0, 1):
                out = (a+c)//2+coin*((a+c) % 2)
                current = a-out
                label = current if law == 'C' else a-c
                groups.setdefault((coin, label), []).append((out, a+c-out))
                effects[int(2*current+2*coin-1 < 0)][code] += 1
        m.need(all(len(outputs) == len(set(outputs)) for outputs in groups.values()),
               'native pointer fine branches have injective support')
        m.need(effects == [[1, 0, 2, 0], [1, 2, 0, 2]], 'native pointer effect counts')
        perm = [affine(x, [('C', 2, 0), ('C', 0, 1)]) for x in range(8)]
        for bit in (0, 1):
            for a in (0, 1):
                for b0 in (0, 1):
                    for c in (0, 1):
                        for d in (0, 1):
                            observed = sum(effects[bit][z]*int(perm[2+4*c] == z+4*a)
                                           *int(perm[2+4*d] == z+4*b0) for z in range(4))
                            m.need(observed == 2*int(a == b0 == c == d == bit), 'complete reduced pointer CP map')
    # Diameter by all-pairs distances on binary mixed-edge moves, rather
    # than the source's endpoint-based reverse BFS routes.
    _, ds = m.incidence(ROOT)
    edges = [(i, j) for i in range(12) for j in range(i+1, 12) if ds[i][j] == 1]
    lookup = {x: i for i, x in enumerate(states)}
    distance = [[0 if i == j else 100 for j in range(66)] for i in range(66)]
    for i, x in enumerate(states):
        for a, c in edges:
            y = list(x); y[a], y[c] = y[c], y[a]
            distance[i][lookup[tuple(y)]] = min(distance[i][lookup[tuple(y)]], 1)
    for pivot in range(66):
        for i in range(66):
            for j in range(66):
                distance[i][j] = min(distance[i][j], distance[i][pivot]+distance[pivot][j])
    diameter = max(map(max, distance))
    m.need(diameter == 4, 'source binary routing diameter')
    timeout = Q(66+diameter, 2**32)
    gate = Q(1, 10**11)
    error = timeout+gate*(2+gate)/2
    m.need(reports[7]['retained_entry'] == [0, 6] and sign(sub(mul(scalar(2), imag(0, 6)), scalar(2*error))) > 0,
           'complete pointer gap survives reset and raw word errors')
    return {'source_word_pulses': 6, 'work_qubits': 9,
        'existing_native_auxiliary_ready_code': 2, 'logical_superoperator_entries_checked': 64,
        'input_auxiliary_independent_of_unknown_work_and_reference': True,
        'reset_timeout_included_in_complete_reduced_work_comparison': True,
        'finer_quantum_auxiliary_outputs_reconstructed_by_probe': False,
        'finer_quantum_auxiliary_outputs_claimed_equal': False,
        'reduced_work_witness_distinguishes_any_finer_retained_output': True,
        'ideal_success_work_maps': 'P_bit0 rho P_bit0, P_bit1 rho P_bit1',
        'records_only_equal_does_not_mean_quantum_outputs_equal': True,
        'bit_reports': reports, 'selected_bit': 7, 'reset_timeout_upper': str(timeout),
        'reset_routing_diameter': diameter, 'complete_reduced_work_error_upper': str(error),
        'raw_word_gate_passed': True, 'complete_quantum_output_gap_certified': True,
        'physical_service_error_included': False,
        'coherent_GV_decoder_or_internal_triplet_export_supplied': False}


def verify(packet, root=ROOT):
    pins = previous.custody(INPUTS, SOURCE_COMMIT, 'FloatingPragma/observer-patch-holography', root)
    old = previous.custody(PRIOR_FILES, PRIOR_COMMIT, 'FloatingPragma/observer-patch-holography', root)
    # The prior receipt remains frozen. Its separate independent verifier is
    # run in CI; its custody here prevents a changed premise silently passing.
    raw_source_gates()
    core = mathematical_packet()
    eps = (1+Q(1, 10**11))**31-1
    bound = eps*(2+eps)/2+Q(33, 2**31)
    m.need(bound < Q(1, 10**7) and sign(sub(tuple(-x for x in core['selected']['ideal_label_0_difference']), scalar(2*bound))) > 0,
           'complete record gap with both errors')
    m.need(sign(sub(core['menu']['minimum_first_distinguishing_label_gap'], scalar(2*bound))) > 0,
           'all distinguishing settings are robust in the bounded model')
    expected = previous.prior.encoded({'schema': 'oph.coherent-native-access-probe.v1', 'work_baseline': BASELINE,
        'source_pins': pins, 'prior_probe_pins': old,
        'implementation_sha256': {p: hashlib.sha256((root/p).read_bytes()).hexdigest() for p in IMPLEMENTATIONS},
        **core, 'existing_native_work_pointer': pointer_packet(),
        'service_bounds': {'primitive_operator_gate': '1/100000000000', 'maximum_certified_blocks': 31,
            'complete_record_half_l1_upper': str(bound), 'portable_upper': '1/10000000',
            'source_detector_timeout_upper': '33/2147483648',
            'raw_state_propagation_used_for_witness': False, 'physical_service_error_included': False,
            'arbitrary_untouched_reference_bound': True,
            'unknown_input_code_preparation_is_supplied_not_derived': True,
            'all_original_binary_source_gates_replayed': True, 'gap_exceeds_two_model_error_bounds': True,
            'every_ideal_distinguishing_setting_is_model_robust': True,
            'exact_raw_equality_claimed_for_ideal_null_settings': False},
        'output_contract': {'labels': [str(i) for i in range(66)]+['timeout', 'export_abort'],
            'records_retains_untouched_reference': True, 'record_service_discards_final_native_work': True,
            'finer_controller_timeout_quantum_outputs_claimed_equal': False,
            'pre_readout_512_native_quantum_output_preserved': True,
            'unitary_target_preserves_orthogonal_inputs_and_references': True,
            'tomography_quantum_payload_exported': False,
            'tomography_uses_coherent_native_work_coupling': False,
            'tomography_uses_uncommitted_native_records': False},
        'verdict': 'COHERENT_ACCESS_WITHOUT_ATTACHMENT', 'phase_0_stop_preserved': True,
        'previous_record_only_grammar_not_overruled': True,
        'universal_coherent_grammar_obstruction_claimed': False,
        'physical_golden_embedding_selected': False,
        'conjugation': 'complex at a fixed real embedding; not Galois',
        'complete_operational_attachment_certified': False})
    m.need(json.dumps(packet, sort_keys=True, allow_nan=False) == json.dumps(expected, sort_keys=True, allow_nan=False),
           'independent coherent-native probe mismatch')
    return {'verified': True, 'verdict': expected['verdict'], 'settings': 4291,
            'marked_input_actions': 60, 'raw_blocks_replayed': 7,
            'method': 'GF(2) routing, Q(phi) effects, independent component intervals on original-binary blocks'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', nargs='?', type=Path, default=REFERENCE)
    args = parser.parse_args()
    print(json.dumps(verify(m.load(args.receipt)), indent=2, sort_keys=True))
