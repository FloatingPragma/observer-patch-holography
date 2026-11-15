"""Finite modular current attachment to the one shared LOCAL12 source.

The full construction has 288 Z6 edge registers. Executed arrays retain one
cycle coordinate, with the other 192 coordinates in their uniform states.
That subspace is invariant, so this compression is exact, not a field cutoff.
"""
from collections import Counter, deque
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
WALK = HERE.parent / 'native_walk_bridge'
sys.path.insert(0, str(WALK))
import common_walk as common
import compile_walk as compiled
from model.openwave_coherent_native_architecture import source as architecture
from model.openwave_native_process import source as original_native
from model.observer.source_scalar import compiler

MODULUS = 6
BASELINE_VERTEX = 0


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@lru_cache(maxsize=1)
def geometry():
    edges = np.array(common.graph().seams, dtype=np.int64)
    B = np.zeros((96, 288), dtype=np.int64)
    adjacency = [[] for _ in range(96)]
    for e, (a, b) in enumerate(edges):
        B[a, e] = 1
        B[b, e] = -1
        adjacency[a].append((int(b), e, 1))
        adjacency[b].append((int(a), e, -1))
    previous = {0: None}
    todo = deque([0])
    while todo:
        a = todo.popleft()
        for b, e, sign in sorted(adjacency[a]):
            if b not in previous:
                previous[b] = (a, e, sign)
                todo.append(b)
    assert len(previous) == 96
    P = np.zeros((288, 96), dtype=np.int64)
    tree = set()
    for x in range(1, 96):
        cursor = x
        while cursor:
            parent, e, sign = previous[cursor]
            P[e, x] = sign
            tree.add(e)
            cursor = parent
    chords = [e for e in range(288) if e not in tree]
    C = np.zeros((288, len(chords)), dtype=np.int64)
    for k, e in enumerate(chords):
        a, b = edges[e]
        C[:, k] = P[:, a] - P[:, b]
        C[e, k] += 1
    G = np.eye(96, dtype=np.int64)[:, 0]
    assert len(tree) == 95 and len(chords) == 193
    assert np.array_equal(B @ P, G[:, None] - np.eye(96, dtype=np.int64))
    assert not np.any(B @ C)
    assert np.array_equal(C[chords], np.eye(193, dtype=np.int64))
    # Chord identity plus leaf elimination on the tree proves integral and
    # modulo-k uniqueness, even when k is composite (no field division).
    chosen = next(e for e in chords if edges[e, 1] < 12)
    return {'edges': edges, 'B': B, 'P': P, 'C': C, 'chords': chords,
            'tree': sorted(tree), 'chosen_edge': chosen, 'G': G}


@lru_cache(maxsize=1)
def current_atoms():
    """Native current for every padded pre-S source/event basis column."""
    g = geometry()
    neighbours, edge_ids = common.geometry_tables()
    after = compiled.controls()['swap_permutation'].copy()
    edge = np.full(2048, -1, dtype=np.int64)
    current = np.zeros(2048, dtype=np.int64)
    for x in range(96):
        for event in range(12):
            d, sigma = divmod(event, 2)
            index = 16 * common.encode_position(x) + event
            e = int(edge_ids[x, d])
            y = int(neighbours[x, d]) if sigma else x
            a, b = g['edges'][e]
            edge[index] = e
            current[index] = int(x == a) - int(y == a)
            assert after[index] == 16 * common.encode_position(y) + event
    assert np.array_equal(after[after], np.arange(2048))
    assert np.array_equal(current[after], -current)
    return after, edge, current


def apply_current_walk(state, selected_edge=None):
    """Actual compiled R/S with one cycle coordinate and its padding retained.

    z=0..5 is the selected chord coordinate; labels6/7 are explicit unitary
    complement labels. Other193-1 cycle coordinates are exactly uniform.
    The source is not reset, measured, or renormalized between coherent steps.
    """
    value = np.asarray(state, dtype=complex)
    if value.shape != (128, 16, 8) or not np.isfinite(value).all():
        raise ValueError('finite padded128x16x8 source/field state required')
    selected_edge = geometry()['chosen_edge'] if selected_edge is None else selected_edge
    if selected_edge not in geometry()['chords']:
        raise ValueError('retained coordinate must be an actual chord')
    reflected = np.einsum('ab,pbz->paz', compiled.controls()['reflection'], value)
    swap, edges, currents = current_atoms()
    result = np.empty((2048, 8), dtype=complex)
    flat = reflected.reshape(2048, 8)
    increment = np.where(edges == selected_edge, currents, 0)
    for z in range(8):
        target_z = (z + increment) % 6 if z < 6 else np.full(2048, z)
        result[swap, target_z] = flat[:, z]
    return result.reshape(128, 16, 8)


def apply_character_walk(state, character, selected_edge=None):
    selected_edge = geometry()['chosen_edge'] if selected_edge is None else selected_edge
    value = np.asarray(state, dtype=complex)
    if value.shape != (128, 16):
        raise ValueError('source shape128x16 required')
    reflected = (value @ compiled.controls()['reflection'].T).reshape(-1)
    swap, edges, currents = current_atoms()
    phase = np.exp(-2j * math.pi * character * np.where(edges == selected_edge, currents, 0) / 6)
    result = np.empty(2048, dtype=complex)
    result[swap] = reflected * phase
    return result.reshape(128, 16)


def ideal_character(character):
    return np.r_[np.exp(2j * math.pi * character * np.arange(6) / 6) / math.sqrt(6), 0j, 0j]


@lru_cache(maxsize=1)
def field_preparation():
    word = compiler.compile_real_state(np.r_[np.ones(6) / math.sqrt(6), 0., 0.])
    U = compiler.execute(word, 3)
    return word, U


@lru_cache(maxsize=1)
def compiled_retained_increment():
    """Actual reversible gates for selected cycle's two oriented increments.

    Source11 bits match DYN01, existing work11..22 is untouched; field23..25,
    extra clean ladder26..36. Complete dirty action is the same reversible
    word, though only all clean columns are exhaustively replayed here.
    """
    return compile_edge_increment(geometry()['chosen_edge'])


def compile_edge_increment(e):
    _, edges, currents = current_atoms()
    gates = []
    local = {}
    for sign in (-1, 1):
        permutation = [(z + sign) % 6 if z < 6 else z for z in range(8)]
        program = architecture.compile_permutation(permutation)
        assert list(map(int, architecture.execute_basis(program))) == permutation
        local[sign] = program
        source_index, = np.flatnonzero((edges == e) & (currents == sign))
        condition = {q: (int(source_index) >> q) & 1 for q in range(11)}
        for kind, targets in program['gates']:
            values = dict(condition)
            if kind == 'ccx':
                values.update({23 + q: 1 for q in targets[:2]})
            gates += compiled.multi_control(values, 23 + targets[-1], tuple(range(26, 37)))
    return gates, local


def positions(state):
    return (abs(state) ** 2).sum(axis=tuple(range(1, state.ndim)))


def run():
    files = [Path(__file__), WALK/'common_walk.py', WALK/'compile_walk.py',
             WALK/'COMMON_WALK_RECEIPT.json']
    dependency_receipt = json.loads((WALK/'COMMON_WALK_RECEIPT.json').read_text())
    for name in dependency_receipt['source_sha256']:
        files.append(ROOT/name)
    files = list(dict.fromkeys(files))
    pins = {str(p.relative_to(ROOT)): sha(p) for p in files}
    assert all(pins[name] == expected for name, expected in dependency_receipt['source_sha256'].items())
    g = geometry()
    swap, edges, currents = current_atoms()
    all_atoms = 0
    all_modular_cases = 0
    native_letters = 0
    local_native = architecture.event_permutation()
    raw_branches = {(row.record[1],row.record[2],row.record[3]): row
                    for row in original_native.raw_local_branches() if row.record[0] == 'repair'}
    for x in range(96):
        for event in range(12):
            i = common.encode_position(x)*16 + event
            y = common.decode_position(int(swap[i])//16)
            e, j = int(edges[i]), int(currents[i])
            difference = g['P'][:, x] - g['P'][:, y]
            difference[e] += j
            h = difference[g['chords']]
            assert np.array_equal(difference, g['C'] @ h)
            assert np.count_nonzero(h) <= 1
            assert np.array_equal(g['B'] @ difference, np.zeros(96, dtype=np.int64))
            all_atoms += 1
            for modulus in (3, 5, 6):
                for field in range(modulus):
                    target = (field+j) % modulus
                    assert (target-j) % modulus == field
                    assert (field+j-target) % modulus == 0
                    all_modular_cases += 1
            for law in ('C', 'E'):
                record = common.raw_record(int(swap[i])//16, event, law)
                a, b = common.graph().seams[e]
                code = 2 if x == a else 1
                native_after = local_native[code | record['coin'] << 6 | int(law == 'E') << 7]
                expected_code = 2 if y == a else 1
                assert native_after & 7 == expected_code
                assert record['current'] == j and record['source_before'] == x
                raw = (native_after >> 3) & (3 if law == 'C' else 7)
                branch = raw_branches[law,record['coin'],raw]
                assert branch.current == j
                assert branch.record[4] == (j if law == 'C' else record['original_E_contrast'])
                assert (1 if x == a else 2,1 if y == a else 2) in branch.pairs
                native_letters += 1
    # Exact finite encoding of the full edge configuration and Gauss sector.
    rng = np.random.default_rng(1979)
    random_fibres = 0
    for modulus in (3, 5, 6):
        for x in range(96):
            z = rng.integers(modulus, size=193)
            F = (g['P'][:, x] + g['C'] @ z) % modulus
            assert np.array_equal(F[g['chords']], z)
            charge = g['G'].copy()
            charge[x] -= 1
            assert not np.any((g['B'] @ F-charge) % modulus)
            random_fibres += 1
    # Exhaustive complete clean basis of the source-controlled extra register.
    gates, local = compiled_retained_increment()
    source = np.tile(np.arange(2048, dtype=np.int64), 8)
    field = np.repeat(np.arange(8, dtype=np.int64), 2048)
    inputs = source | (field << 23)
    actual = compiled.execute_reversible(gates, inputs)
    increment = np.where(edges[source] == g['chosen_edge'], currents[source], 0)
    expected_field = np.where(field < 6, (field+increment) % 6, field)
    expected = source | (expected_field << 23)
    assert np.array_equal(actual, expected)
    assert np.array_equal(compiled.execute_reversible(list(reversed(gates)), actual), inputs)
    gate_counts = Counter(kind for kind, _ in gates)
    primitive_resources = compiled.controls()['primitive_resources']
    increment_pulses = sum(gate_counts[k]*primitive_resources[k]['native_pulses'] for k in gate_counts)
    full_gate_counts = Counter()
    full_blocks = []
    for e in range(288):
        edge_gates, _ = compile_edge_increment(e)
        counts = Counter(kind for kind, _ in edge_gates)
        full_gate_counts.update(counts)
        full_blocks.append({'edge_index': e, 'field_bits': list(range(23+3*e,26+3*e)),
                            'clean_scratch': list(range(887,898)),
                            'gate_counts': dict(counts)})
    full_increment_pulses = sum(full_gate_counts[k]*primitive_resources[k]['native_pulses'] for k in full_gate_counts)
    # Actual source-word field preparation. Unitary complement labels stay in
    # every array; ideal-character calculations below are exact theorem tests.
    field_word, field_U = field_preparation()
    field_error = float(np.linalg.norm(field_U[:, 0]-ideal_character(0)))
    assert field_error < 1e-11
    initial = common.initial_cell_state()
    source_states = [initial.copy() for _ in range(6)]
    joint_states = [initial[:, :, None]*ideal_character(a) for a in range(6)]
    definite = np.zeros((128, 16, 8), complex)
    definite[:, :, 0] = initial
    character_words = []
    actual_prepared = []
    actual_character_preparation_errors = []
    for a in range(6):
        phase_word = compiler.compile_diagonal([2*math.pi*a*z/6 for z in range(6)]+[0.,0.])
        word = field_word + phase_word
        prepared = compiler.execute(word,3)[:,0]
        error = float(np.linalg.norm(prepared-ideal_character(a)))
        assert error < 1e-11
        character_words.append(word)
        actual_character_preparation_errors.append(error)
        actual_prepared.append(initial[:,:,None]*prepared)
    export_checks = 0
    maximum_field_export_complement_weight = 0.
    for kind, angle in {(p['kind'],p['angle']) for word in character_words for p in word}:
        d = 4 if kind == 'J' else 2
        full = compiler.native_pulse(kind,angle)
        channels = compiler.transfer_kraus(d)
        assert np.linalg.norm(sum(K.conj().T@K for K in channels)-np.eye(6),'fro') < 1e-13
        assert np.linalg.norm(sum((K@full[:,:d]).conj().T@(K@full[:,:d]) for K in channels)-np.eye(d),'fro') < 1e-11
        maximum_field_export_complement_weight = max(maximum_field_export_complement_weight,
            float(sum(np.linalg.norm(K@full[:,:d],'fro')**2 for K in channels[1:])))
        export_checks += 1
    rows = []
    maximum_intertwining = maximum_norm_error = maximum_mixture_error = 0.
    for step in range(9):
        for a in range(6):
            error = float(np.linalg.norm(joint_states[a]-source_states[a][:, :, None]*ideal_character(a)))
            maximum_intertwining = max(maximum_intertwining, error)
            maximum_norm_error = max(maximum_norm_error, abs(float(np.vdot(joint_states[a], joint_states[a]).real)-1))
        uniform = source_states[0]
        mixture = sum(positions(s) for s in source_states)/6
        mixture_error = float(np.max(abs(positions(definite)-mixture)))
        maximum_mixture_error = max(maximum_mixture_error, mixture_error)
        record = {'coherent_step': step,
                  'character1_vs_zero_position_TV': float(abs(positions(source_states[1])-positions(uniform)).sum()/2),
                  'definite_cycle_vs_zero_position_TV': float(abs(positions(definite)-positions(uniform)).sum()/2),
                  'definite_cycle_equals_six_character_mixture_max_error': mixture_error,
                  'compiled_character_prep_vs_ideal_state_norm': [float(np.linalg.norm(actual_prepared[a]-joint_states[a])) for a in range(6)],
                  'complete_position_probabilities_character0': positions(uniform).tolist(),
                  'complete_position_probabilities_character1': positions(source_states[1]).tolist(),
                  'complete_position_probabilities_definite_cycle': positions(definite).tolist(),
                  'source_event_field_total_probability': float(np.vdot(definite, definite).real),
                  'field_padding_probability': float((abs(definite[:, :, 6:])**2).sum()),
                  'renormalized': False}
        rows.append(record)
        if step < 8:
            for a in range(6):
                source_states[a] = apply_character_walk(source_states[a], a)
                joint_states[a] = apply_current_walk(joint_states[a])
                actual_prepared[a] = apply_current_walk(actual_prepared[a])
            bare = common.apply_walk(uniform)
            assert np.linalg.norm(source_states[0]-bare) < 1e-12
            definite = apply_current_walk(definite)
    assert maximum_intertwining < 1e-11 and maximum_mixture_error < 1e-11
    assert max(r['character1_vs_zero_position_TV'] for r in rows) > .001
    assert max(r['definite_cycle_vs_zero_position_TV'] for r in rows) > .001
    # Arbitrary reference is untouched; full padding is explicitly exercised.
    random_state = rng.normal(size=(128,16,8,2))+1j*rng.normal(size=(128,16,8,2))
    random_state /= np.linalg.norm(random_state)
    out = np.stack([apply_current_walk(random_state[:, :, :, r]) for r in range(2)], axis=3)
    reference_gram_error = float(np.max(abs(out.reshape(-1,2).conj().T @ out.reshape(-1,2)-random_state.reshape(-1,2).conj().T @ random_state.reshape(-1,2))))
    assert reference_gram_error < 1e-11
    # No modular field can distinguish accumulated winding k from zero.
    aliases = [{'modulus': k, 'integer_histories': [0,k], 'same_register': 0,
                'bounded_charge_labels': [k-1,0,1], 'bounded_decoder_injective': len({k-1,0,1}) == 3}
               for k in (3,5,6)]
    packet_bounds = [{'integer_packet_half_width': K, 'one_step_state_error_bound': math.sqrt(2/(2*K+1)),
                      'eight_step_state_error_bound': min(2.,8*math.sqrt(2/(2*K+1))),
                      'eight_step_hardware_half_width_without_wrap': K+8}
                     for K in (8,64,1024,65536)]
    assert pins == {str(p.relative_to(ROOT)): sha(p) for p in files}
    result = {
        'schema': 'oph.common.local12.current_attachment.v1', 'model_id': common.MODEL_ID,
        'checks': {'actual_source_R_and_reversible_S_reused': True,
                   'all1152_active_atoms_exact_integer_current_Gauss': all_atoms == 1152,
                   'all2304_C_E_terminal_labels_replayed': native_letters == 2304,
                   'finite_modular_current_inverse_and_padding': True,
                   'primitive_integral_tree_chord_coordinates': True,
                   'all16384_source_field_clean_columns_source_gate_executed': True,
                   'finite_normalized_affine_fibre_only': True,
                   'six_character_intertwiners_and_reference_completeness': True,
                   'definite_cycle_which_path_control_changes_interference': True,
                   'integer_alias_and_packet_controls_explicit': True,
                   'source_pins_unchanged': True},
        'geometry': {'vertices': 96, 'edges': 288, 'tree_edges': 95, 'cycle_coordinates': 193,
                     'chosen_chord_index': g['chosen_edge'], 'chosen_chord': g['edges'][g['chosen_edge']].tolist(),
                     'baseline': 'G=delta_vertex0', 'incidence': '+ at canonical lower endpoint, - at upper',
                     'finite_fibre_cardinality': '6^193', 'normalization': '6^(-193/2)',
                     'moduli_tested': [3,5,6], 'random_full_fibre_cases': random_fibres,
                     'all_local_modular_current_cases': all_modular_cases},
        'construction': {'full_field': '(Z6)^288; each edge uses 3 bits with6/7 unitary complement',
                         'update': 'S_F|x,e,F> = |S(x,e),e,F+j_e epsilon_e>; R changes only event',
                         'affine_code': 'F=P_x+Cz; BP_x=G-delta_x, BC=0, C[chords]=I193',
                         'normalizable_exact_intertwiner': 'J_G|x,e>=6^(-193/2) sum_z |x,e,P_x+Cz>',
                         'not_claimed': ['uniform normalized integer fibre', 'tracing F recovers bare source density', 'autonomous field energy or radiation', 'native attempt clock equals coherent step'],
                         'executed_compression': 'one explicit Z6 chord coordinate; other192 coordinates exactly uniform; source/event/field padding retained',
                         'integer_packet_bound': 'For integer cube[-K,K]^193, T-step state-vector error <= min(2,T sqrt(2/(2K+1))); finite implementation needs coordinate range at least[-K-T,K+T] to avoid wrap.'},
        'native_increment': {'retained_coordinate_reversible_gates': len(gates), 'gate_counts': dict(gate_counts),
                             'native_pulses': increment_pulses, 'complete_clean_columns': len(inputs),
                             'extra_field_bits': 3, 'extra_clean_scratch_bits': 11,
                             'full_physical_edge_register_synthesis': 'Apply the same controlled six-cycle on each of288 edge registers, using its two actual endpoint/event predicates; then execute unchanged source S. All blocks commute. No time step or field-energy phase is inserted.',
                             'full_edge_increment_pulse_count': full_increment_pulses,
                             'full_edge_increment_gate_counts': dict(full_gate_counts),
                             'full_edge_field_bits': 864, 'full_common_work_bits_with_increment_scratch': 898,
                             'dirty_input_scope': 'Each emitted flip/ccx gate defines a complete unitary dirty-scratch extension. Gauss/current correctness is on the declared clean scratch/code sector; failures and complements are retained.'},
        'field_preparation': {'actual_three_bit_source_pulse_resources': compiler.resources(field_word),
                              'actual_state_norm_error': field_error,
                              'six_actual_character_source_preparation_errors': actual_character_preparation_errors,
                              'distinct_field_prep_native_pulse_complete_export_checks': export_checks,
                              'maximum_field_prep_export_complement_Frobenius_weight': maximum_field_export_complement_weight,
                              'full_affine_preparation_recipe': 'Reset field and scratch with retained native failure branches. Source-word prepare uniform6 on193 chord registers; coherently add C tree rows from chord values and P_x from source-position predicates. This is a finite reversible modular-addition circuit and declared joint correlated preparation, not a product definite edge tape.',
                              'full_898_cell_reset_or_full_affine_pulse_word_executed': False,
                              'complete_failure_contract': 'Inherit every M6 proper-code load/export complement and native reset failure. Any extra register reset failure is a tagged terminal record with final carrier state; no division by success. This receipt executes the prepared code-sector map and its one-coordinate native circuit, not a full898-cell reset instrument.'},
        'numerics': {'maximum_character_intertwining_state_error': maximum_intertwining,
                     'maximum_norm_error': maximum_norm_error, 'maximum_mixture_position_error': maximum_mixture_error,
                     'arbitrary_reference_Gram_error': reference_gram_error},
        'records': rows, 'modular_alias_controls': aliases, 'integer_packet_bounds': packet_bounds,
        'physical_field_trace_control': {'input_logical_state': '(|vertex0>+|vertex1>)/sqrt2',
                                         'physical_field_fibres_disjoint': True,
                                         'reason': 'BF=G-delta_x differs for the two vertices modulo6',
                                         'traced_source_density': [[0.5,0.],[0.,0.5]],
                                         'bare_logical_density': [[0.5,0.5],[0.5,0.5]],
                                         'trace_distance': 0.5,
                                         'no_bare_source_density_transfer_claimed': True},
        'clock': {'evolution': 'same coherent LOCAL12 W step', 'extra_current_service_pulses_separate': True,
                  'physical_time_calibration': None},
        'source_pins': pins, 'official_rows_changed': [], 'reserved_acceptance_cases_run': False}
    (HERE/'CURRENT_WALK.json').write_text(json.dumps(result, indent=2)+'\n')
    program = {'schema': 'oph.common.current.native_increment_program.v1',
               'model_id': common.MODEL_ID, 'retained_cycle_increment_gates': gates,
               'full_physical_edge_blocks': full_blocks,
               'block_generator': 'compile_edge_increment(e); remap field23..25 to block field_bits, scratch26..36 to887..897; source bits0..10 unchanged',
               'primitive_native_words': compiled.controls()['primitive_words'],
               'field_uniform_preparation_word': field_word,
               'six_field_character_preparation_words': character_words,
               'source_clock_order': 'actual R, controlled current blocks on pre-S source/event, unchanged actual S',
               'source_pins': pins}
    (HERE/'CURRENT_PROGRAM.json').write_text(json.dumps(program, indent=2)+'\n')
    print(json.dumps({'checks': result['checks'], 'numerics': result['numerics'],
                      'increment_pulses': increment_pulses,
                      'max_character_position_TV': max(r['character1_vs_zero_position_TV'] for r in rows),
                      'max_definite_position_TV': max(r['definite_cycle_vs_zero_position_TV'] for r in rows)}, indent=2))


if __name__ == '__main__':
    run()
