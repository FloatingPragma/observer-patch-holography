"""Execute typed reducing factors without expanding repeated full-source words."""
from collections import Counter
from fractions import Fraction as F
import json
import math
import sys


def quantum(workspace):
    import numpy as np
    directory = workspace / 'openwave-tuning/common_dynamics/quantum_field'
    sys.path.insert(0, str(directory))
    import native_cycle as nc
    frozen = json.loads((directory / 'NATIVE_CYCLE.json').read_text())
    saved_program = json.loads((directory / 'NATIVE_CYCLE_PROGRAM.json').read_text())
    instruments, _ = nc.phase_clock.actual_native_pointer()
    primitive = json.loads((workspace / 'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json').read_text())
    primitive = {k: F(v) for k, v in primitive['primitive_operator_errors_exact'].items()}
    adder = sum((number * primitive[k] for k, number in Counter(k for k, _ in nc.nf.add_program()).items()), F())
    pivot, blocks = nc.encoder_blocks()
    assert pivot == saved_program['physical_pivot'] and blocks == saved_program['encoder_blocks']
    p0 = min(map(F, nc.configuration.success_probabilities(4, 32, full=True, mode='reset', target=nc.circuit.decode(0))))
    pr = min(map(F, nc.configuration.success_probabilities(4, 32, full=True, mode='reset', target=nc.circuit.decode(nc.pointer.READY_CODE))))
    reset = 3204 * (1-p0) + 5 * (1-pr)
    assert str(reset) == frozen['full_initialization_reset_union_bound_exact']
    rows = []
    for setting, expected, resource in zip(('X', 'Y'), frozen['experiments'], frozen['full_resources']):
        word, error, parts = nc.compressed_word(setting)
        assert word == saved_program['compressed_five_bit_words'][setting]
        actual = nc.recovery.execute_complete(word, 5, np.eye(32, dtype=complex)[:, [0]])
        state = actual.success[:, 0]
        target = parts['initial_state'].copy()
        for _ in range(8):
            target = nc.qs.source_step(target, parts['data']['alpha'], parts['data']['beta'])
        clock = nc.compiler.execute(nc.phase_clock.local_words(8, setting, F(1, 8)), 1, np.array([1, 0], complex))
        residual = float(np.linalg.norm(state - np.kron(clock, np.fft.fft(target, norm='ortho'))))
        assert residual < float(error) + 2e-12
        assert str(error) == expected['compressed_native_error_exact']
        laws = {}
        maximum = 0.
        for law, branches in instruments.items():
            effects = [np.zeros((2, 2), complex), np.zeros((2, 2), complex)]
            fine = []
            for metadata, kraus in branches:
                effect = sum((K.conj().T @ K for K in kraus), np.zeros((2, 2), complex))
                effects[metadata['bit']] += effect
                fine.append({'record': metadata, 'effect_real': effect.real.tolist(), 'effect_imag': effect.imag.tolist()})
            values = []
            for outcome in range(32):
                effect = nc.joint_effect(effects, [(outcome >> j) & 1 for j in range(5)])
                values.append(float(np.vdot(state, effect @ state).real))
            archived = expected['native_pointer_laws'][law]
            maximum = max(maximum, float(np.max(abs(np.asarray(values)-archived['joint_angle_bin_clock_probabilities']))))
            assert fine == archived['fine_one_bit_effects']
            assert abs(sum(values)-float(np.vdot(state, state).real)) < 2e-11
            laws[law] = values
        assert maximum < 2e-12
        electric, electric_error = nc.nf.electric_word(parts['data']['alpha'])
        full_step = 288 * electric_error + parts['magnetic_error'] + 2 * len(blocks) * adder
        bound = reset + parts['preparation_error'] + 2 * len(blocks) * adder + 8 * full_step
        bound += parts['qft_error'] + nc.phase_clock.exact_word_bound(8, setting, F(1, 8)) + 4 * nc.phase_clock.pointer_bound()
        assert str(bound) == expected['full_physical_program_record_TV_bound_exact']
        # The full generator consists of this exact face factor, 288 remapped
        # electric factors, and the declared preparation/encoder/clock words.
        prep = nc.state_preparation_word()[0]
        count = len(prep) + 2 * len(nc.encoder_word()) + 2 * len(nc.nf.H(0))
        count += 8 * (len(nc.nf.face_word(0, parts['data']['beta'])) + 288 * len(electric) + 1)
        count += len(nc.nf.qft_word()) + int(setting == 'Y')
        assert count == resource['full_word_pulses']
        rows.append({'setting': setting, 'native_pulses_executed': len(word),
                     'complete_32_outcome_laws': laws, 'all_fine_pointer_effects_exact': True,
                     'native_export_failure_probability': float(actual.failure_gram[0, 0].real),
                     'maximum_archived_probability_error': maximum, 'native_state_error': residual,
                     'full_word_pulses_by_exact_factor_count': count,
                     'complete_record_TV_bound_exact': str(bound)})
    return {'schema': 'oph.common_dynamics.public_quantum_factors.v1', 'experiments': rows,
            'execution_kind': 'Actual five-bit complete native instrument plus exact physical encoder and repeated-factor count.',
            'full_primitive_stream_reexecuted': False,
            'domain': 'Pure-gauge one-face reducing code, zero fermion occupations; distinct preparation and coefficients from the massive source.',
            'official_rows_changed': []}


def reference(workspace):
    import numpy as np
    directory = workspace / 'openwave-tuning/common_dynamics/vacuum'
    sys.path.insert(0, str(directory))
    import native_preparation as npref
    ref, compiler, cert, nf = npref.ref, npref.compiler, npref.cert, npref.nf
    frozen = json.loads((directory / 'NATIVE_PREPARATION.json').read_text())
    saved_program = json.loads((directory / 'PREPARATION_PROGRAM.json').read_text())
    values = ref.rational_gaussian_amplitudes(16)
    target = np.array(values, dtype=float) / math.sqrt(sum(v*v for v in values))
    word = compiler.compile_real_state(target)
    error, _ = npref.native_cycle.gray_preparation_bound(word, values)
    assert word == saved_program['gaussian_word'] and str(error) == frozen['gaussian_source_error_exact']
    actual = ref.dp.spin.recovery.execute_complete(word, 4, np.eye(16, dtype=complex))
    residual = float(np.linalg.norm(actual.success[:, 0]-target))
    assert residual < float(error) + 2e-13
    program = npref.prepare_blocks()
    for key, value in program.items():
        assert json.loads(json.dumps(value)) == saved_program[key]
    geometry = ref.dp.field.current_walk.geometry()
    encoded_C = np.zeros((288, 193), dtype=np.int64)
    for j, edge in enumerate(program['chord_registers']):
        encoded_C[edge, j] = 1
    chord_index = {e: j for j, e in enumerate(program['chord_registers'])}
    for block in program['cycle_shears']:
        encoded_C[block['target_edge'], chord_index[block['source_edge']]] += block['sign']
    assert np.array_equal(encoded_C, geometry['C'])
    encoded_P = np.zeros((288, 32), dtype=np.int64)
    for block in program['occupation_field_additions']:
        encoded_P[block['target_edge'], block['occupation']] += block['sign']
    assert np.array_equal(encoded_P, geometry['P'][:, np.repeat(np.arange(8)*12, 4)])
    controlled = npref.compiled.execute_reversible(npref.increment_program(), np.arange(256, dtype=np.int64))
    expected = np.array([(x & 1) + ((((x >> 1) & 15) + (x & 1)) % 16)*2 for x in range(32)])
    assert np.array_equal(controlled[:32], expected) and len(np.unique(controlled)) == 256
    assert np.array_equal(controlled >> 5, np.arange(256) >> 5)
    assert np.array_equal(npref.compiled.execute_reversible(list(reversed(npref.increment_program())), controlled), np.arange(256))
    primitive = json.loads((workspace / 'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json').read_text())
    primitive = {k: F(v) for k, v in primitive['primitive_operator_errors_exact'].items()}
    counts = Counter({'flip': len(program['matter_seed_occupied_modes']) + len(program['tree_offset_flip_wires'])})
    for k, n in Counter(k for k, _ in nf.add_program()).items():
        counts[k] += n * len(program['cycle_shears'])
    for k, n in Counter(k for k, _ in npref.increment_program()).items():
        counts[k] += n * len(program['occupation_field_additions'])
    total_error = 193 * error + sum((n * primitive[k] for k, n in counts.items()), F())
    slater_count = 0
    for a, b, sign in ref.pairs():
        first = ['I'] * (b+1)
        second = first.copy()
        first[a], first[b], second[a], second[b] = 'X', 'Y', 'Y', 'X'
        for k in range(a+1, b):
            first[k] = second[k] = 'Z'
        template = cert.pauli(''.join(first), cert._pi(-sign, 8)) + cert.pauli(''.join(second), cert._pi(sign, 8))
        native = ref.dp.field.mb.givens_word(a, b, sign*math.pi/4)
        total_error += nf.root_certificate.exact_template_bound(native, template)
        slater_count += len(native)
    primitives = npref.compiled.controls()['primitive_words']
    count = slater_count + 193 * len(word) + sum(n * len(primitives[k]) for k, n in counts.items())
    assert count == frozen['native_preparation_pulses'] and str(total_error) == frozen['full_preparation_source_error_exact']
    p0 = min(map(F, npref.configuration.success_probabilities(4, 32, full=True, mode='reset', target=npref.circuit.decode(0))))
    pr = min(map(F, npref.configuration.success_probabilities(4, 32, full=True, mode='reset', target=npref.circuit.decode(npref.pointer.READY_CODE))))
    reset = 1187 * (1-p0) + 1184 * (1-pr)
    pointer_error = nf.root_certificate.exact_template_bound(npref.pointer._local_word(), cert.cnot(2, 0) + cert.cnot(0, 1))
    bound = min(F(1), total_error + reset + 1184 * pointer_error)
    assert str(bound) == frozen['complete_preparation_data_record_TV_bound_exact']
    laws, denominator = ref.complete_occupation_law()
    assert len(laws) == denominator == 65536 and all(sign in (-1, 1) for _, sign in laws)
    rng = np.random.default_rng(141413)
    for index in (0, 1, 17, 511, 4095, 65535):
        for _ in range(3):
            z = rng.integers(0, 16, size=193)
            assert np.array_equal(npref.classical_encoder(laws[index][0], z, program), ref.dressed_field(laws[index][0], z, ref.reference_charge()))
    return {'schema': 'oph.common_dynamics.public_reference_factors.v1',
            'gaussian_all_16_columns_executed': True, 'gaussian_native_residual': residual,
            'complete_shared_193_cycle_and_32_occupation_encoder': True,
            'controlled_increment_all_256_columns_checked': True,
            'complete_occupation_outcomes': len(laws), 'encoder_record_checks': 18,
            'native_preparation_pulses_by_exact_factor_count': count,
            'complete_record_TV_bound_exact': str(bound),
            'full_primitive_stream_reexecuted': False,
            'domain': 'Distinct 32-mode Slater/Gauss preparation and terminal occupation/field instrument. No clock or full interacting Fock evolution.',
            'official_rows_changed': []}
