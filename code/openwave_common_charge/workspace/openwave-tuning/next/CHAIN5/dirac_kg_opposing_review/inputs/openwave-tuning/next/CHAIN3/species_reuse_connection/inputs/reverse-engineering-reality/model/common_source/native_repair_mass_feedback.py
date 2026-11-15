"""Literal native-word implementation of the fixed five-memory Higgs feedback.

The neutral family pair, its Nambu qubit, complete native records, old decoder
syndromes and fresh encoder blanks are explicit. This implements the existing
specified kernel; it neither changes its failed hierarchy nor transfers the
effective tiny-time readout intervals to a finite source history.
"""
from __future__ import annotations

from functools import lru_cache
from fractions import Fraction
import hashlib
import math

import numpy as np
import sympy as s

from . import native_family_selection as family
from . import native_repair_mass as feedback
from . import native_code_transport as transport
from . import native_gauge_object as objects
from ..observer.source_scalar import compiler

ROOT = objects.ROOT
STEP = 1/64
SOURCE_ADMISSION = 2e-8
ARRAY_ADMISSION = 2e-10
NATIVE_SCHEMA = {key: family.ARRAY_SCHEMA[key] for key in
                 ('family_C', 'family_E', 'family_C_labels', 'family_E_labels')}
ARRAY_SCHEMA = {
    'feedback_port_rays': ((5, 64), 'float64'),
    'feedback_source_unitaries': ((5, 128, 128), 'complex128'),
    'feedback_source_preparation': ((64, 3), 'complex128'),
    'feedback_source_input_images': ((5, 128, 6), 'complex128'),
    'feedback_encoder_C': ((512, 8), 'complex128'),
    'feedback_encoder_E': ((512, 8), 'complex128'),
    'feedback_complete_grams': ((5, 128, 128), 'complex128'),
    'feedback_branch_labels': ((640, 6), 'int64')}


@lru_cache(maxsize=1)
def geometry():
    d, law = family.geometry(), feedback.architecture()
    item = objects.geometric_intertwiners()['three']
    T = item['frame']/s.sqrt(item['metric'])
    # P has independently certified rank3. An isometric three-column frame
    # inside its range therefore spans it, without a redundant36x36 expansion.
    if not family._zero(T.T*T-s.eye(3)) or not family._zero(d['P']*T-T):
        raise ArithmeticError('the feedback triplet is not the selected native family band')
    ports = [law['seed']]
    for _ in range(4): ports.append(law['rotation'][ports[-1]])
    rays = []
    for port, Q in zip(ports, law['projectors'], strict=True):
        ray = 2*d['frame'][:, port]
        coordinate = T.T*ray
        if not family._zero(T*coordinate-ray) or not family._zero(coordinate*coordinate.T-Q):
            raise ArithmeticError('a literal port ray differs from the fixed feedback projector')
        full = np.zeros(64)
        full[transport.product_addresses()] = np.asarray(ray, float).reshape(36)
        rays.append(full)
    frame = np.zeros((64, 3), complex)
    frame[transport.product_addresses()] = np.asarray(T, complex)
    return {'ports': tuple(ports), 'rays': np.asarray(rays), 'frame': frame}


def _word(ray):
    # exp(-ih X tensor P) has phases -h,+h on |+ ray>,|- ray>.
    # The Y rotation supplies these two basis rays up to harmless signs.
    basis = compiler._compile_real_state_tree(ray, 6)+[compiler.pulse('Y', -math.pi/4, [6])]
    phases = np.zeros(128); phases[0] = -STEP; phases[64] = STEP
    return compiler.inverse(basis)+compiler.compile_diagonal(phases)+basis


def ideal_unitaries():
    X = np.asarray([[0, 1], [1, 0]])
    return np.asarray([np.eye(128)+(math.cos(STEP)-1)*np.kron(np.eye(2), np.outer(v, v))
                       -1j*math.sin(STEP)*np.kron(X, np.outer(v, v)) for v in geometry()['rays']])


def branch_labels():
    return np.asarray([(cc, cr, ec, er, old, (old+1) % 5 if event else 0)
                       for old in range(5) for cc, cr, ec, er, event, _ in feedback.native_records()], np.int64)


def _channel(values):
    value = values.reshape(-1, 64, 8, 8)
    return np.einsum('rhia,rhjc->ijac', value.conj(), value, optimize=True)


def complete_grams(arrays):
    """All128 input/reference directions, all640 labels and old syndromes."""
    labels = arrays['family_C_labels']
    advance = (labels[:, 0] == 0)&(labels[:, 1] == 1)
    C = [_channel(arrays['family_C'][~advance]), _channel(arrays['family_C'][advance])]
    E = _channel(arrays['family_E'])
    encoders = [arrays['feedback_encoder_'+law] for law in ('C', 'E')]
    final_gram = np.kron(np.eye(2), np.kron(*(U.conj().T@U for U in encoders)))
    results = []
    for old in range(5):
        gram = np.zeros((128, 128), complex)
        for event in (0, 1):
            new = (old+1) % 5 if event else 0
            U = arrays['feedback_source_unitaries'][new]
            effect = (U.conj().T@final_gram@U).reshape(2, 8, 8, 2, 8, 8)
            gram += np.einsum('ijac,nikmjl,klbd->nabmcd', C[event], effect, E,
                              optimize=True).reshape(128, 128)
        results.append(gram)
    return np.asarray(results)


def validate_arrays(arrays):
    schema = dict(NATIVE_SCHEMA, **ARRAY_SCHEMA)
    if set(arrays) != set(schema): raise ValueError('complete feedback source factors required')
    for name, (shape, dtype) in schema.items():
        value = arrays[name]
        if (type(value) is not np.ndarray or value.shape != shape or value.dtype != np.dtype(dtype)
                or not np.isfinite(value).all()): raise ValueError('invalid feedback array '+name)
    if not np.array_equal(arrays['feedback_branch_labels'], branch_labels()):
        raise ValueError('a native record or old/new memory label differs')
    if np.max(abs(arrays['feedback_port_rays']-geometry()['rays'])) > ARRAY_ADMISSION:
        raise ValueError('a native port ray changed')


def _allowances():
    # Exact arithmetic on the admitted binary float, with one upward rounding
    # at export. No platform libm power enters a fixed budget comparison.
    a = Fraction.from_float(SOURCE_ADMISSION); pair = a*(1+a)+a
    cycle = (1+a)**2*pair+(1+a)*a+pair
    prepared = cycle+(1+a)**3*a
    low_dilation = (1+a)*pair+a
    low_channel = ((1+a)**2+1)*low_dilation
    kappa = (1+a)**4
    geometric_sum = sum(kappa**j for j in range(64))
    repeated_low = low_channel*geometric_sum
    complete_repeated = ((1+a)**2*(repeated_low+kappa**64*(2+a)*a)+(2+a)*pair)
    exact = {'cycle_dilation_allowance': cycle,
        'prepared_six_input_dilation_allowance': prepared,
        'cycle_channel_allowance': (2+cycle)*cycle,
        'low_register_cycle_channel_allowance': low_channel,
        'retained_low_cycle_CP_norm_upper': kappa,
        'sixty_four_cycle_geometric_amplification_sum': geometric_sum,
        'sixty_four_low_cycle_telescope_channel_allowance': repeated_low,
        'prepared_sixty_four_cycle_final_native_output_channel_allowance': complete_repeated}
    return {key: math.nextafter(float(value), math.inf) for key, value in exact.items()}


def metrics(arrays):
    validate_arrays(arrays)
    ideal = ideal_unitaries(); pulses = arrays['feedback_source_unitaries']
    pulse_error = max(float(np.linalg.norm(a-b)) for a, b in zip(pulses, ideal, strict=True))
    pulse_gram = max(float(np.linalg.norm(U.conj().T@U-np.eye(128), 2)) for U in pulses)
    T = np.kron(np.eye(2), geometry()['frame'])
    A = np.kron(np.eye(2), arrays['feedback_source_preparation'])
    preparation_error = float(np.linalg.norm(A-T))
    images = np.asarray([U@A for U in pulses])
    if np.max(abs(images-arrays['feedback_source_input_images'])) > ARRAY_ADMISSION:
        raise ValueError('six-input source feedback images differ from their complete word factors')
    selected_error = max(float(np.linalg.norm(images[k]-ideal[k]@T)) for k in range(5))
    native_errors, encoder_errors, native_norms, encoder_norms, local_grams = [], [], [], [], []
    for law, seam in (('C', (2, 3)), ('E', (2, 10))):
        code = transport.PAIR_CODES[law]; values = arrays['family_'+law]
        labels = arrays['family_'+law+'_labels']
        expected_labels = np.asarray([(coin, raw) for coin in (0, 1)
                                      for raw in range(4 if law == 'C' else 8)], np.int64)
        if not np.array_equal(labels, expected_labels): raise ValueError('all native raw branches are required')
        target = np.zeros((len(labels), 512, 6), complex)
        for n, (coin, raw) in enumerate(labels.tolist()):
            target[n, :6] = math.sqrt(float(code.branch_weight(law, seam, coin, raw)))*np.eye(6)
        native_errors.append(float(np.linalg.norm(values[:, :, :6]-target)))
        gram = np.einsum('rha,rhb->ab', values.conj(), values, optimize=True)
        local_grams.append(float(np.linalg.norm(gram-np.eye(8), 2)))
        native_norms.append(math.sqrt(float(np.linalg.eigvalsh(gram).max())))
        encoder = arrays['feedback_encoder_'+law]
        encoder_errors.append(float(np.linalg.norm(encoder[:, :6]-code.ideal_code_matrix())))
        eg = encoder.conj().T@encoder
        local_grams.append(float(np.linalg.norm(eg-np.eye(8), 2)))
        encoder_norms.append(math.sqrt(float(np.linalg.eigvalsh(eg).max())))
    if max(native_errors+encoder_errors+[pulse_error, preparation_error, selected_error/2]) > SOURCE_ADMISSION:
        raise ValueError('complete feedback source exceeds the fixed component allowance')
    if max(local_grams+[pulse_gram]) > SOURCE_ADMISSION:
        raise ValueError('a complete source operation lost an input/reference direction')
    grams = complete_grams(arrays)
    if np.max(abs(grams-arrays['feedback_complete_grams'])) > ARRAY_ADMISSION:
        raise ValueError('complete native feedback Gram differs from its factors')
    gram_error = max(float(np.linalg.norm(g-np.eye(128), 2)) for g in grams)
    if gram_error > 8*SOURCE_ADMISSION: raise ValueError('complete native feedback lost probability')
    native_delta = native_errors[0]*(1+native_errors[1])+native_errors[1]
    encoder_delta = encoder_errors[0]*(1+encoder_errors[1])+encoder_errors[1]
    norm_encoder = math.prod(encoder_norms)
    norm_pulse = max(float(np.linalg.norm(U, 2)) for U in pulses)
    cycle = norm_encoder*norm_pulse*native_delta+norm_encoder*pulse_error+encoder_delta
    prepared = cycle+norm_encoder*norm_pulse*math.prod(native_norms)*preparation_error
    # Iterate at the decoded low-register interface. The next native encoder
    # is the preceding final encoder, not a second encoding of encoded data.
    # These factors are numerically admitted, so use a CP norm bound kappa
    # instead of silently assuming that their retained Gram is exactly I.
    return {'full128_source_unitarity_defect_diagnostic': pulse_gram,
        'full128_feedback_word_error_diagnostic': pulse_error,
        'six_input_preparation_error_diagnostic': preparation_error,
        'six_unknown_input_feedback_error_diagnostic': selected_error,
        'local_native_dilation_errors_diagnostic': native_errors,
        'local_reencoder_errors_diagnostic': encoder_errors,
        'complete128_native_feedback_Gram_defect_diagnostic': gram_error,
        'factor_proved_cycle_dilation_error_diagnostic': cycle,
        'factor_proved_prepared_dilation_error_diagnostic': prepared,
        'component_dilation_allowance': SOURCE_ADMISSION,
        **_allowances(),
        'fixed_budget_arithmetic': 'exact rational arithmetic on the admitted binary float; one upward-rounded float export',
        'effective_correlated_averaging_bound_at_t1_N64': feedback.averaging_bound(64),
        'input_domain_for_closeness': 'selected native family3 tensor Nambu2 with arbitrary untouched quantum reference',
        'full128_completeness_implies_identity_closeness_on_padding': False,
        'numerical_allowances_are_outward_hardware_certificates': False}


def execute(native_arrays=None):
    values = family._native_factors() if native_arrays is None else native_arrays
    arrays = {key: values[key] for key in NATIVE_SCHEMA}
    d = geometry(); resources, pulses = [], []
    for ray in d['rays']:
        word = _word(ray)
        pulses.append(compiler.execute(word, 7, np.eye(128, dtype=complex)))
        resources.append(compiler.resources(word))
    arrays.update({'feedback_port_rays': d['rays'].copy(),
        'feedback_source_unitaries': np.asarray(pulses),
        'feedback_source_preparation': transport.prepare_geometric_source()[:, :3],
        'feedback_branch_labels': branch_labels()})
    for law in ('C', 'E'):
        arrays['feedback_encoder_'+law] = transport.PAIR_CODES[law]._encode(np.eye(512, 8, dtype=complex))
    A = np.kron(np.eye(2), arrays['feedback_source_preparation'])
    arrays['feedback_source_input_images'] = np.asarray([U@A for U in pulses])
    arrays['feedback_complete_grams'] = complete_grams(arrays)
    report = {'classification': 'literal fixed five-memory feedback composed with complete native repair',
        'dimensionless_coherent_step': STEP, 'fixed_native_ports': list(d['ports']),
        'ordering': 'native encoders; all C/E branches and recovery; full decoders; record memory update; controlled feedback word; fresh-syndrome native reencoders',
        'source_word': 'S diag(exp(-ih) at0, exp(ih) at64,1 elsewhere) S_inverse; S prepares the native ray and the Nambu X eigenbasis',
        'resources_for_all_five_words': resources,
        'complete_domains': {'native_family_pair_low_register': 64, 'Nambu_dimension': 2,
            'all_feedback_word_input_columns': 128, 'selected_input_reference_columns': 6,
            'old_memory_states': 5, 'native_raw_labels': 128, 'joint_classical_labels': 640,
            'old_decoder_syndrome_quantum_dimension': 4096,
            'new_native_output_rows_per_carrier': 512, 'new_encoder_input_columns_per_carrier': 8,
            'all_padding_inputs_and_abort_branches_retained': True,
            'old_decoder_syndromes_erased_or_conditioned_away': False,
            'fresh_blank_syndrome_service_for_reencoding': True},
        'composition': {'all_joint_branch_maps_defined_by_executed_source_factors': True,
            'all_joint_output_amplitudes_materialized': False,
            'complete_Gram_contracts_every_raw_label_and_syndrome': True,
            'iterable_interface': 'decoded low128 after feedback and before reencoding; that encoder is reused as the next native input encoder',
            'additional_final_reencoder_executed_once': True,
            'retained_matrix_telescope_includes_CP_norm_amplification': True,
            'exact_word_CPTP_and_floating_Gram_admission_are_distinct': True,
            'target_projectors_are_applied_instead_of_source_words': False,
            'ideal_comparison_maps_are_used_only_for_error_certificates': True},
        'execution': metrics(arrays),
        'scope': {'implements_the_existing_fixed_failing_kernel': True,
            'coherent_compiler_memory_controller_and_blank_services_are_declared': True,
            'primitive_C_E_alone_produces_the_feedback': False,
            'family_pair_is_separate_and_neutral_under_internal_gauge_currents': True,
            'full_physical_matter_or_spacetime_field_history_executed': False,
            'effective_tiny_time_readout_intervals_transferred_to_this_source_history': False,
            'finite_cycle_source_error_is_effective_averaging_error': False,
            'physical_mass_scale_or_measured_hierarchy_derived': False},
        'empirical_inputs': []}
    return report, arrays


def report(): return execute()[0]


def source_bindings():
    names = set(family.source_bindings()) | {
        'model/common_source/native_repair_mass.py', 'model/common_source/test_native_repair_mass.py',
        'model/common_source/native_repair_mass_feedback.py', 'model/common_source/test_native_repair_mass_feedback.py'}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}
