"""M6 source-pulse and reversible-word compiler for the shared LOCAL12 walk.

No target Hamiltonian, propagator or spectral eigenbasis is an execution input.
Finite logical contractions are proved from complete native primitive words.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
import json
import math
from pathlib import Path

import numpy as np

import common_walk as common
from model.observer.source_scalar import compiler
from model.openwave_coherent_native_architecture import source as architecture
from model.openwave_native_process import source as original_native
from model.native_interface import pointer
from model.native_quantum_service_bridge import construction as service

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RER = ROOT / 'reverse-engineering-reality'
DATA_BITS = 11
WORK_SOURCE = (11, 12, 13)
RAW_RECORD = (14, 15, 16)
LAW_BIT = 17
SCRATCH = tuple(range(18, 23))
MEMORY_BITS = 23


def remap_gates(gates, mapping):
    return [(kind, tuple(mapping[q] for q in targets)) for kind, targets in gates]


def multi_control(values, target, scratch=SCRATCH):
    """Complete reversible Toffoli ladder, including both zero-control flips."""
    controls = sorted(values)
    if len(controls) < 2 or target in controls or set(controls+[target]) & set(scratch):
        raise ValueError('distinct controls, target and scratch required')
    if len(scratch) < len(controls)-2:
        raise ValueError('insufficient clean scratch for ladder')
    zeros = [q for q in controls if values[q] == 0]
    out = [('flip', (q,)) for q in zeros]
    if len(controls) == 2:
        out += [('ccx', (*controls, target))]
    else:
        chain = [(controls[0], controls[1], scratch[0])]
        for j in range(2, len(controls)-1):
            chain.append((scratch[j-2], controls[j], scratch[j-1]))
        out += [('ccx', row) for row in chain]
        out += [('ccx', (scratch[len(controls)-3], controls[-1], target))]
        out += [('ccx', row) for row in reversed(chain)]
    out += [('flip', (q,)) for q in reversed(zeros)]
    return out


def execute_reversible(gates, inputs):
    out = np.asarray(inputs, dtype=np.int64).copy()
    for kind, targets in gates:
        if kind == 'flip':
            out ^= 1 << targets[0]
        elif kind == 'ccx':
            a, b, c = targets
            out ^= (((out >> a)&1) & ((out >> b)&1)) << c
        elif kind == 'cx':
            a, b = targets
            out ^= ((out >> a)&1) << b
        else:
            raise ValueError(kind)
    return out


def endpoint_extraction():
    out = []
    for direction, matching in enumerate(common.INTERNAL_MATCHINGS):
        for a, b in matching:
            for port, work_bit in ((a, WORK_SOURCE[1]), (b, WORK_SOURCE[0])):
                values = {j+1: (direction >> j)&1 for j in range(3)}
                values.update({j+4: (port >> j)&1 for j in range(4)})
                out += multi_control(values, work_bit)
    # Canonical external seam ordering at size2 is determined by the
    # toggled cell coordinate. Port groups0,1,2 route z,x,y respectively.
    for group, cell_bit in ((0, 8), (1, 10), (2, 9)):
        for coordinate in (0, 1):
            values = {j+1: (5 >> j)&1 for j in range(3)}
            values.update({6: group&1, 7: group >> 1, cell_bit: coordinate})
            out += multi_control(values, WORK_SOURCE[1 if coordinate == 0 else 0])
    return out


def internal_swap():
    permutation = list(range(256))
    for direction, matching in enumerate(common.INTERNAL_MATCHINGS):
        for a, b in matching:
            left, right = 16*a+2*direction+1, 16*b+2*direction+1
            permutation[left], permutation[right] = permutation[right], permutation[left]
    program = architecture.compile_permutation(permutation)
    # Complete clean and dirty local input domain, tensor identity on cells.
    actual = architecture.execute_basis(program)
    assert tuple(map(int, actual[:256])) == tuple(permutation)
    assert len(set(map(int, actual))) == len(actual) == 8192
    mapping = {q: q for q in range(8)}
    mapping.update({8+j: SCRATCH[j] for j in range(5)})
    return remap_gates(program['gates'], mapping), program


def external_swap():
    out = []
    for group, cell_bit in ((0, 8), (1, 10), (2, 9)):
        values = {j: (11 >> j)&1 for j in range(4)}
        values.update({6: group&1, 7: group >> 1})
        out += multi_control(values, cell_bit)
        out += multi_control(values, 5)
    return out


def uniform_word():
    # 12=3*4: a source-compiled three-component block and two actual Hadamards.
    high = compiler.compile_real_state(np.array([1., 1., 1., 0.])/math.sqrt(3))
    word = compiler.remap(high, [2, 3])
    H = architecture.primitive_routines()['hadamard']
    word += compiler.remap(H, [0]) + compiler.remap(H, [1])
    return word


@lru_cache(maxsize=1)
def controls():
    g = common.graph()
    primitive_matrices, primitive_errors, primitive_resources = architecture.execute_primitives()
    cx_word = compiler.cnot(0, 1)
    cx = compiler.execute(cx_word, 2)
    cx_target = np.eye(4)[:, [0, 3, 2, 1]]
    primitive_errors['cx'] = float(np.linalg.norm(cx-cx_target, 'fro'))
    assert primitive_errors['cx'] < 1e-10
    primitive_resources['cx'] = compiler.resources(cx_word)
    extraction = endpoint_extraction()
    inside, inside_program = internal_swap()
    outside = external_swap()
    local_mapping = {0:11, 1:12, 2:13, 3:14, 4:15, 5:16, 6:0, 7:17,
                     8:18, 9:19, 10:20, 11:21, 12:22}
    local_programs = {}
    for name, permutation in (('native_event', architecture.event_permutation()),
                              ('record_recoder', architecture.record_recoder())):
        program = architecture.compile_permutation(permutation)
        truth = architecture.execute_basis(program)
        assert tuple(map(int, truth[:256])) == permutation
        assert len(np.unique(truth)) == 8192
        local_programs[name] = program
    native_gates = remap_gates(local_programs['native_event']['gates'], local_mapping)
    recoder_gates = remap_gates(local_programs['record_recoder']['gates'], local_mapping)
    source_program = (extraction + [('cx', (12, 0))] + native_gates + recoder_gates
                      + inside + outside + list(reversed(extraction)))
    original_label_program = extraction + list(reversed(recoder_gates)) + list(reversed(extraction))
    inputs = np.concatenate([np.arange(2048, dtype=np.int64), np.arange(2048, dtype=np.int64) | (1 << LAW_BIT)])
    outputs = execute_reversible(source_program, inputs)
    expected = common.expected_swap_permutation()
    assert np.array_equal(outputs[:2048], expected)
    assert np.array_equal(outputs[2048:], expected | (1 << LAW_BIT))
    assert len(np.unique(outputs)) == 4096
    # Emit and execute the inverse recoder before a raw-label intervention.
    readable = execute_reversible(original_label_program, inputs)
    from model.observer.source_discrimination import circuit
    for vertex, event, law in product(range(96), range(12), range(2)):
        data = 16*common.encode_position(vertex)+event
        value = int(readable[data+2048*law])
        record = common.raw_record(common.encode_position(vertex), event, 'E' if law else 'C')
        assert (value >> 4)&127 == common.encode_position(vertex)
        assert value&1 == record['coin']
        assert (value >> 14)&7 == circuit.CURRENT_CODES[record['current']]+4*law
        assert (value >> 11)&7 == 0 and value >> 18 == 0
    # The native committed instrument is checked independently on every local
    # binary source matrix unit, including all25 raw slots and original labels.
    slots = 0
    for x, y in product(range(4), repeat=2):
        actual = architecture.committed_column(architecture.BINARY_CODES[x], architecture.BINARY_CODES[y])
        wanted = {}
        for row in original_native.raw_local_branches():
            if not row.weight or not row.pairs:
                continue
            action = dict(row.pairs)
            if x in action and y in action:
                _, law, coin, raw, _ = row.record
                wanted[raw+8*coin+16*int(law == 'E')] = ((architecture.BINARY_CODES[action[x]], architecture.BINARY_CODES[action[y]]), row.weight)
        assert actual == wanted
        slots += 25
    preparation_word = uniform_word()
    U = compiler.execute(preparation_word, 4)
    reflection_word = (compiler.inverse(preparation_word)
                       + compiler.compile_diagonal([0.]+[math.pi]*15)
                       + preparation_word)
    reflection = compiler.execute(reflection_word, 4)
    u = np.r_[np.ones(12)/math.sqrt(12), np.zeros(4)]
    prep_error = float(np.linalg.norm(U[:, 0]-u))
    reflection_error = float(np.linalg.norm(reflection-(2*np.outer(u, u)-np.eye(16)), 'fro'))
    assert prep_error < 1e-11 and reflection_error < 1e-10
    assert np.linalg.norm(reflection[:, 12:]+np.eye(16)[:, 12:], 'fro') < 1e-10
    counts = Counter(kind for kind, _ in source_program)
    pulses_S = sum(counts[kind]*primitive_resources[kind]['native_pulses'] for kind in counts)
    primitive_words = architecture.primitive_routines()
    primitive_words['cx'] = cx_word
    return {'uniform12_preparation': U, 'reflection': reflection,
            'swap_permutation': outputs[:2048].copy(),
            'source_program': source_program, 'primitive_words': primitive_words,
            'original_label_program': original_label_program,
            'preparation_word': preparation_word, 'reflection_word': reflection_word,
            'primitive_errors': primitive_errors, 'primitive_resources': primitive_resources,
            'local_programs': local_programs, 'inside_program': inside_program,
            'local_native_matrix_unit_slots': slots,
            'preparation_error': prep_error, 'reflection_full_phase_error': reflection_error,
            'source_swap_gate_counts': dict(counts), 'source_swap_native_pulses': pulses_S,
            'source_swap_complete_clean_columns': len(inputs)}


def pulse_word():
    """One complete W word as pulse-address records; no target matrix input."""
    built = controls()
    word = list(built['reflection_word'])
    for kind, targets in built['source_program']:
        word += compiler.remap(built['primitive_words'][kind], list(targets))
    return word


def preparation_word():
    """Full source/event/law preparation after retained native blank resets."""
    built = controls()
    return (built['preparation_word']
            + compiler.remap(built['preparation_word'], [4, 5, 6, 7])
            + compiler.remap(built['primitive_words']['hadamard'], [LAW_BIT]))


def classical_thinning_contract():
    """Complete global-to-local population comparison, with finite wait failure."""
    g = common.graph()
    M = len(g.seams)
    contact = F(6, M)
    no_contact = 1-contact
    cap = 16
    # One of six seams, two coins and two law values is the terminal local
    # letter. Every omitted global letter still has its seam/coin/law record.
    each_local_letter_success = sum(no_contact**k/F(4*M) for k in range(cap))
    timeout = no_contact**cap
    assert 24*each_local_letter_success+timeout == 1
    assert each_local_letter_success == (1-timeout)/24
    from model.observer.native import native
    examples = []
    for vertex in (0, 47, 95):
        n = tuple(int(x == vertex) for x in range(96))
        incident = next(edge for edge in g.seams if vertex in edge)
        idle = [edge for edge in g.seams if vertex not in edge][:2]
        records = []
        for event, edge in enumerate(idle+[incident]):
            coin = event % 2
            after, (_, current) = native(n, edge, coin)
            records.append({'seam': list(edge), 'coin': coin, 'law': 'C' if event%2 == 0 else 'E',
                            'current': current, 'nonincident': vertex not in edge})
            if vertex not in edge:
                assert after == n and current == 0
            n = after
        local, (_, current) = native(tuple(int(x == vertex) for x in range(96)), incident, 0)
        assert n == local and records[-1]['current'] == current
        examples.append({'initial_vertex': vertex, 'complete_global_tape': records,
                         'global_fine_probability': str(F(1, (4*M)**len(records))),
                         'final_source': list(n), 'selected_local_seam': list(incident)})
    return {'contact_probability': str(contact), 'nonincident_probability': str(no_contact),
            'uncapped_wait_law': 'P(k nonincident letters then contact)=(47/48)^k/48',
            'mean_global_attempts_per_local_attempt': 48, 'finite_wait_cap': cap,
            'each_of24_local_letter_success_probabilities': str(each_local_letter_success),
            'finite_wait_timeout_probability': str(timeout), 'complete_probability': '1',
            'complete_examples': examples,
            'Poisson_population_identity': 'r_global*(T_global-I)=(r_global/48)*(T_local-I)',
            'scope': 'Diagonal classical position process and complete classical attempted-event tapes only; no equality of raw quantum channels or coherent walks follows.'}


def run():
    paths = [Path(__file__), HERE/'common_walk.py',
             RER/'model/native_spatial_refinement/source.py',
             RER/'model/openwave_coherent_native_architecture/source.py',
             RER/'model/openwave_native_process/source.py',
             RER/'model/observer/source_scalar/compiler.py',
             RER/'model/observer/source_discrimination/circuit.py',
             RER/'model/native_interface/pointer.py',
             RER/'model/native_clock/configuration.py',
             RER/'model/native_quantum_service_bridge/construction.py']
    summary = compiler.source_summary()
    paths += [RER/name for name in summary['source_sha256']]
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    built = controls()
    g = common.graph()
    D, N = 12, 96
    # Compression is reconstructed from the executed reversible circuit.
    compressed = np.zeros((N, N), dtype=np.int64)
    for x in range(N):
        for event in range(D):
            image = int(built['swap_permutation'][16*common.encode_position(x)+event])
            y = common.decode_position(image//16)
            assert y is not None
            compressed[y, x] += 1
    assert np.array_equal(compressed, common.transfer_numerators())
    original_columns = common.spatial.native_columns(g)
    L = D*np.eye(N)-compressed
    global_numerators = (2*len(g.seams))*np.eye(N)-L
    transfer = compressed/D
    assert np.max(abs(global_numerators/(2*len(g.seams)) - ((1-F(6, len(g.seams)))*np.eye(N)+F(6, len(g.seams))*transfer))) < 1e-14
    eigenvalues, eigenvectors = np.linalg.eigh(transfer)
    assert min(eigenvalues) > 0 and abs(max(eigenvalues)-1) < 1e-13
    spectral_residual = 0.
    u = np.r_[np.ones(12)/math.sqrt(12), np.zeros(4)]
    for i, value in enumerate(eigenvalues):
        A = np.zeros((128, 16), dtype=complex)
        for vertex in range(96):
            A[common.encode_position(vertex)] = eigenvectors[vertex, i]*u
        SA = np.empty(2048, dtype=complex)
        SA[built['swap_permutation']] = A.reshape(-1)
        SA = SA.reshape(128, 16)
        if value < 1-1e-12:
            phase = np.exp(1j*math.acos(float(value)))
            state = A-phase*SA
            residual = np.linalg.norm(common.apply_walk(state)-phase*state)/np.linalg.norm(state)
        else:
            residual = np.linalg.norm(common.apply_walk(A)-A)
        spectral_residual = max(spectral_residual, float(residual))
    assert spectral_residual < 2e-10
    source_packet = common.initial_cell_state()
    wave_checks = []
    for step in range(5):
        record = common.complete_position_records(source_packet)
        assert abs(record['original_norm_squared']-1) < 2e-10
        wave_checks.append({'coherent_step': step, 'complete_probability': record['original_norm_squared'],
                           'source_padding': record['source_padding_probability'],
                           'event_padding': record['event_padding_probability_on_valid_source'],
                           'position_probabilities': record['valid_position_probabilities'].tolist()})
        source_packet = common.apply_walk(source_packet)
    # Every terminal physical source/event/law letter replays the literal
    # endpoint rule. Original-label measurement is a declared intervention.
    raw_checks = 0
    from model.observer.native import native
    for y, event, law in product(range(96), range(12), ('C', 'E')):
        record = common.raw_record(common.encode_position(y), event, law)
        before = tuple(int(x == record['source_before']) for x in range(96))
        after, (_, current) = native(before, tuple(record['seam']), record['coin'])
        assert current == record['current'] and after[y] == 1
        raw_checks += 1
    all_words = [*built['primitive_words'].values(), built['preparation_word'], built['reflection_word']]
    export_checks = 0
    maximum_failure = 0.
    for kind, angle in {(p['kind'], p['angle']) for word in all_words for p in word}:
        d = 4 if kind == 'J' else 2
        full = compiler.native_pulse(kind, angle)
        channels = compiler.transfer_kraus(d)
        assert np.linalg.norm(sum(K.conj().T@K for K in channels)-np.eye(6), 'fro') < 1e-14
        gram = sum((K@full[:, :d]).conj().T@(K@full[:, :d]) for K in channels)
        assert np.linalg.norm(gram-np.eye(d), 'fro') < 1e-12
        maximum_failure = max(maximum_failure, float(sum(np.linalg.norm(K@full[:, :d], 'fro')**2 for K in channels[1:])))
        export_checks += 1
    blank_reset = service.reset_receipt(0)
    auxiliary_reset = service.reset_receipt(pointer.READY_CODE)
    pointer_certificate = pointer.certificate()
    word = pulse_word()
    resources = compiler.resources(word)
    prep_word = preparation_word()
    thinning = classical_thinning_contract()
    assert resources['native_pulses'] == built['source_swap_native_pulses']+len(built['reflection_word'])
    assert pins == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    program = {'model_id': common.MODEL_ID, 'descriptor': common.source_descriptor(),
               'memory_bits': MEMORY_BITS, 'data_bits': DATA_BITS,
               'source_endpoint_work_bits': WORK_SOURCE, 'original_raw_record_bits': RAW_RECORD,
               'law_bit': LAW_BIT, 'scratch_bits': SCRATCH,
               'preparation_word_event_or_local_port': built['preparation_word'],
               'full_preparation_word': prep_word,
               'reflection_word': built['reflection_word'],
               'primitive_words': built['primitive_words'], 'swap_routine_calls': built['source_program'],
               'original_label_readout_routine_calls': built['original_label_program'],
               'terminal_native_bit_readout': {'work_bits': list(range(MEMORY_BITS)),
                   'local_source_word': pointer._local_word(),
                   'local_wire_order': 'auxiliary_bit0,auxiliary_bit1,selected_work_bit',
                   'global_auxiliary_wires': [MEMORY_BITS, MEMORY_BITS+1],
                   'reset_before_each_bit': 'native configuration ready_code2; timeout retained'},
               'full_pulse_expander': 'compile_walk.pulse_word',
               'source_sha256': pins}
    (HERE/'COMMON_WALK_PROGRAM.json').write_text(json.dumps(program, indent=2)+'\n')
    return {'schema': 'oph.common.local12.native_walk.v1', 'model_id': common.MODEL_ID,
            'checks': {'actual_geometry_six_perfect_matchings': True,
                       'full_native_event_and_inverse_record_recoder_reused': True,
                       'all4096_clean_data_and_law_columns_compiled': True,
                       'complete_local_dirty_scratch_circuits_retained': True,
                       'uniform12_and_reflection_full_phase_source_executed': True,
                       'source_compression_equals_actual_local_native_transfer': True,
                       'exact_global_population_thinning_identity': True,
                       'same_walk_spectral_and_packet_checks': True,
                       'all2304_terminal_original_native_letters_replayed': True,
                       'complete_M6_export_and_native_reset_failures_retained': True,
                       'source_pins_unchanged': True},
            'descriptor': common.source_descriptor(), 'native_global_columns': original_columns,
            'local_binary_matrix_unit_slots': built['local_native_matrix_unit_slots'],
            'complete_clean_data_law_columns': 4096, 'local_dirty_basis_columns_per_event_or_recoder': 8192,
            'primitive_source_errors': built['primitive_errors'],
            'uniform12_preparation_error': built['preparation_error'],
            'reflection_full_phase_error': built['reflection_full_phase_error'],
            'all96_spectral_modes_maximum_residual': spectral_residual,
            'spectrum_contract': {'cos_walk_phase': '1-laplacian_eigenvalue/12',
                'valid_source_event_dimension': 1152, 'dispersive_eigenvalues': 190,
                'active_plus_one_multiplicity': 194, 'active_minus_one_multiplicity': 768,
                'event_padding_minus_one_modes': 384,
                'source_padding_plus_one_modes': 32, 'source_padding_minus_one_modes': 480,
                'physical_speed_or_Lorentz_limit_inferred': False},
            'packet_records': wave_checks, 'terminal_native_record_cases': raw_checks,
            'primitive_export_completeness_checks': export_checks,
            'maximum_numerical_export_failure_effect_trace': maximum_failure,
            'one_W_source_resources': resources, 'swap_gate_counts': built['source_swap_gate_counts'],
            'source_preparation_resources': compiler.resources(prep_word),
            'classical_global_attempt_reconstruction': thinning,
            'preparation_reset': blank_reset, 'native_readout_auxiliary_reset': auxiliary_reset,
            'native_bit_readout_certificate': {k: str(v) if isinstance(v, F) else v for k, v in pointer_certificate.items()},
            'complete_reset_readout_budget': {
                'blank_work_resets': MEMORY_BITS, 'terminal_auxiliary_resets': MEMORY_BITS,
                'reset_failure_upper_union_bound': str(MEMORY_BITS*(F(blank_reset['certificate']['failure_upper'])+F(auxiliary_reset['certificate']['failure_upper']))),
                'native_reset_and_read_attempt_cap': MEMORY_BITS*(blank_reset['certificate']['native_attempts_max']+auxiliary_reset['certificate']['native_attempts_max']+1),
                'terminal_read_source_pulses': MEMORY_BITS*pointer_certificate['source_pulses'],
                'bound_applies_without_independent_input_preparation': True,
                'private_reset_records_retained': True},
            'accuracy_scope': {'ideal_source_circuit_identity': 'Exact Boolean/phase identities under the declared ideal native primitive and code services.',
                'floating_diagnostics': 'Actual fullM6 primitive matrices and complete4bit reflection are executed; the23bit reversible circuit is executed by its exact certified Boolean factorization, not by a dense raw floating432425pulse product.',
                'accumulated_physical_pulse_error_certified': False,
                'whole_word_channel_accuracy_inferred_from_local_float_checks': False,
                'conditional_error_composition': 'A separately certified per-routine operator error can be composed through the emitted routine count; the displayed floating residuals are not such a physical certificate.'},
            'complete_instrument_recipe': 'Run conclusive native blank resets for23work cells; on any failure retain timeout and final work. Prepare localport/event uniform12 and law plus with emitted pulses. Run requested W words with every M6 export complement. Terminal readout uses the complete native pointer perworkbit, with fresh auxiliary reset or retained timeout. Classify source/event/scratch padding explicitly; never divide original outcomes by success.',
            'clock_contract': 'Native globalattempt, native localattempt, pulsewordcount, proper-code-service count and coherentWstep are separate. Global-to-local thinning is a population/process comparison with original nonincident coin/current letters retained. There is no identification of coherent W steps with physical time.',
            'thinning_probability_exact': '1/48', 'mean_global_attempts_per_local_attempt': 48,
            'native_global_quantum_channel_equals_local_or_W': False,
            'program_sha256': hashlib.sha256((HERE/'COMMON_WALK_PROGRAM.json').read_bytes()).hexdigest(),
            'source_sha256': pins, 'official_rows_changed': []}


if __name__ == '__main__':
    result = run()
    (HERE/'COMMON_WALK_RECEIPT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('checks', 'uniform12_preparation_error',
        'reflection_full_phase_error', 'all96_spectral_modes_maximum_residual', 'one_W_source_resources')}, indent=2))
