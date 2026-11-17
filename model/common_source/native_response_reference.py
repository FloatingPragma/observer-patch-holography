"""A supplied native Bell reference tests a maximal compatible current algebra.

The actual exchange response is enlarged by the native fivefold class-sum
phase. Preserving a pure second-triplet Bell ray then leaves thirteen raw
currents, or twelve on End(M6) after the common scalar phase is quotiented.
This is an operational compatibility condition, not autonomous selection of
the reference, protection, phase control, or physically inaccessible controls.

The source experiment retains a factorization of four native carriers. Every
raw record and coherent decoder-syndrome amplitude is reconstructible; no
512**4 array, branch postselection, or replacement by ideal scalar Kraus
operators is used. Floating diagnostics are not physical error certificates.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import sympy as s

from . import native_response_comparison as comparison
from . import native_sm_current as current, native_code_transport as transport
from . import native_persistent_response as response, record_code, native_routing
from .mass_response import history
from ..observer.source_scalar import compiler
from ..observer.source_discrimination import circuit

ROOT = Path(__file__).resolve().parents[2]
SOURCE_ADMISSION = 2e-7
REPLAY_ADMISSION = 4e-9
LAWS = ('C', 'E', 'C', 'E')


def _entries(matrix):
    return [[str(x) for x in matrix.row(i)] for i in range(matrix.rows)]


def _bell(offset=3):
    value = s.zeros(36, 1)
    for j in range(offset, offset+3):
        value[6*j+j] = 1
    return s.ImmutableMatrix(value)


def _ray_kernel(basis, ray):
    norm = (ray.conjugate().T*ray)[0]
    projection = s.eye(36)-ray*ray.conjugate().T/norm
    images = [comparison._clean(projection*(s.kronecker_product(g, s.eye(6))
                +s.kronecker_product(s.eye(6), g))*ray) for g in basis]
    constraints = comparison._columns(tuple(images))
    vectors = constraints.to_DM(extension=True).nullspace().to_Matrix()
    generators = tuple(comparison._clean(sum((a*g for a, g in zip(vectors.row(j), basis, strict=True)),
                        s.zeros(6))) for j in range(vectors.rows))
    return constraints, generators


@lru_cache(maxsize=1)
def _exact_certificate():
    data = comparison.exact_data()
    union, native = data['union_basis'], data['native_basis']
    constraints, compatible = _ray_kernel(union, _bell())
    old_constraints, old = _ray_kernel(native, _bell())
    mixed_constraints, mixed = _ray_kernel(union, _bell(0)+_bell())
    quotient = tuple(comparison._clean(g-s.trace(g[3:, 3:])*s.eye(6)/3) for g in compatible)
    nonzero = tuple(g for g in quotient if any(g))
    pivots = comparison._columns(nonzero).to_DM(extension=True).rref()[1]
    section = tuple(nonzero[j] for j in pivots)
    comparison._coordinates(data['declared_basis'], section)
    comparison._coordinates(section, data['declared_basis'])
    observed = lambda basis: comparison._rank(comparison._columns(tuple(
        comparison._clean(g-s.trace(g)*s.eye(6)/6) for g in basis)))
    if ((len(compatible), observed(compatible), len(old), observed(old),
         len(mixed), observed(mixed), comparison._rank(constraints)) != (13, 12, 12, 11, 7, 6, 5)):
        raise ArithmeticError('native Bell-reference stabilizer or central quotient changed')
    # The native class sum is independently calculated from all actual port
    # rotations. Its two eigenvalues supply the missing relative phase.
    native_phase = response.rotation_endogeneity()
    class_sum = s.diag(*([2+2*s.sqrt(5)]*3+[2-2*s.sqrt(5)]*3))
    coordinates = response.geometric_coordinates()
    actual = s.Matrix(native_phase['near_fivefold_class_sum_exact'])
    if not comparison._zero(data['native_chart_inverse']*actual*data['native_chart']-class_sum):
        raise ArithmeticError('the promoted phase is not the actual native class sum')
    augmented = native+(comparison._clean(s.I*class_sum),)
    comparison._coordinates(augmented, union)
    comparison._coordinates(union, augmented)
    witnesses = data['extra_hermitian_basis']
    witness_images = tuple(comparison._clean((s.kronecker_product(s.I*h, s.eye(6))
                          +s.kronecker_product(s.eye(6), s.I*h))*_bell()/s.sqrt(3)) for h in witnesses)
    witness_gram = comparison._clean(s.Matrix.hstack(*witness_images).conjugate().T
                                     *s.Matrix.hstack(*witness_images))
    if witness_gram.det() == 0:
        raise ArithmeticError('the reference does not detect all five excluded phase directions')
    return {
        'native_exchange_dimension': len(native), 'class_sum_augmented_dimension': len(union),
        'native_class_sum_exact': _entries(class_sum),
        'reference': '(e3 tensor e3 + e4 tensor e4 + e5 tensor e5)/sqrt(3)',
        'collective_action': 'X tensor I6 + I6 tensor X',
        'ray_constraint_real_rank': comparison._rank(constraints),
        'ray_stabilizer_raw_dimension': len(compatible), 'observable_dimension_on_EndM6': observed(compatible),
        'ray_stabilizer_basis_exact': [_entries(g) for g in compatible],
        'constraint_matrix_exact': _entries(constraints),
        'canonical_section_basis_exact': [_entries(g) for g in section],
        'canonical_section_rule': 'X -> X-tr(X_second)/3*I6',
        'canonical_section_equals_existing_twelve_port_current_span': True,
        'maximality_identity': 'For X=diag(A,B), ray preservation iff B+B^T=(2 tr(B)/3)I3; anti-Hermitian B=R+i b I3, R real skew.',
        'without_class_sum': {'constraint_rank': comparison._rank(old_constraints),
            'raw_dimension': len(old), 'observable_dimension': observed(old)},
        'coherent_mixed_Bell_reference_control': {'constraint_rank': comparison._rank(mixed_constraints),
            'raw_dimension': len(mixed), 'observable_dimension': observed(mixed)},
        'five_excluded_current_reference_Gram_exact': _entries(witness_gram),
        'excluded_diagonal_current': 'i diag(0,0,0,1,-1,0)',
        'excluded_diagonal_survival': '(1+2*cos(2*t))^2/9',
        'central_action': 'The raw common iI6 is silent on EndM6; after quotient the remaining central action is the relative triplet phase.',
        'reference_and_class_sum_phase_are_added_source_choices': True,
        'native_four_carrier_ground_reference_identified_with_this_pure_ray': False,
        'incompatible_controls_physically_inaccessible': False,
        'local_gauge_redundancy_or_SM_matter_derived': False,
    }


def exact_certificate():
    from copy import deepcopy
    return deepcopy(_exact_certificate())


def _reference_vector():
    value = np.zeros(64)
    value[[8*j+j for j in range(3, 6)]] = 1/math.sqrt(3)
    return value


def _phase_word(matrix, parameter):
    """Compile a real symmetric phase with all eight padded columns retained."""
    full = np.zeros((8, 8))
    full[:6, :6] = np.asarray(matrix, float)
    values, basis = np.linalg.eigh(full)
    return (compiler.compile_real_orthogonal(basis.T)
            +compiler.compile_diagonal(float(parameter)*values)
            +compiler.compile_real_orthogonal(basis))


def _control_menu():
    data = comparison.exact_data()
    rows = [('identity', F(0), s.zeros(6), compiler.compile_diagonal(np.zeros(8)))]
    for port, generator in enumerate(data['declared_basis']):
        request = current.current_request('port', F(1, 64), port=port)
        rows.append((f'compatible_port_{port}', F(1, 64), generator,
                     current._logical_word(request, targets=(0, 1, 2))))
    class_sum = s.diag(*([2+2*s.sqrt(5)]*3+[2-2*s.sqrt(5)]*3))
    rows.append(('native_class_sum', F(1, 16), s.I*class_sum, _phase_word(class_sum, F(1, 16))))
    rows.append(('common_phase', F(1, 4), s.I*s.eye(6), _phase_word(s.eye(6), F(1, 4))))
    for i, h in enumerate(data['extra_hermitian_basis']):
        rows.append((f'excluded_phase_{i}', F(1, 4), s.I*h, _phase_word(h, F(1, 4))))
    return tuple(rows)


def _targets():
    result = []
    for _, parameter, generator, _ in _control_menu():
        value = np.eye(8, dtype=complex)
        value[:6, :6] = expm(float(parameter)*np.asarray(generator, complex))
        result.append(value)
    return np.asarray(result)


def _index(value, length, name):
    if type(value) is not int or not 0 <= value < length:
        raise ValueError(name+' must be an original in-range integer')
    return value


def contract_output(arrays, control, records, syndromes, pointer, input_columns=None, *, geometric_flag=0):
    """One block of the actual factored instrument, with 4096 output rows.

    Four native raw records are classical. Four high-syndrome indices label
    coherent output tensor legs, NOT additional measurements. Varying these
    indices reconstructs the complete native output, including aborts. The
    remaining rows retain both geometric blank/padding registers. All source
    matrices used below are actual executed factors, never target matrices.
    """
    control = _index(control, len(arrays['controls']), 'control')
    pointer = _index(pointer, 2, 'pointer')
    geometric_flag = _index(geometric_flag, 2, 'geometric return flag')
    if type(records) is not tuple or len(records) != 4 or type(syndromes) is not tuple or len(syndromes) != 4:
        raise ValueError('four original record pairs and four syndrome indices required')
    selected = []
    for law, record, high in zip(LAWS, records, syndromes, strict=True):
        if type(record) is not tuple or len(record) != 2:
            raise ValueError('each complete native record is a coin/raw pair')
        coin = _index(record[0], 2, 'coin')
        raw = _index(record[1], 1 << circuit.RECORD_BITS[law], 'raw record')
        high = _index(high, 64, 'coherent high syndrome')
        row = coin*(1 << circuit.RECORD_BITS[law])+raw
        selected.append(arrays['decoded_'+law][row, 8*high:8*(high+1)])
    state = arrays['bell_preparation'][:, 0] if input_columns is None else input_columns
    if (type(state) is not np.ndarray or state.ndim not in (1, 2) or state.shape[0] != 64
            or state.dtype.kind not in 'iufc' or not np.isfinite(state).all()):
        raise ValueError('finite64-row input/reference columns required')
    vector = state.ndim == 1
    state = np.asarray(state, complex).reshape(64, -1)
    state = np.kron(arrays['controls'][control], arrays['controls'][control])@state
    A = arrays['geometric_preparation'][:, 8*np.arange(8)]
    state = np.einsum('ai,ijr,bj->abr', A, state.reshape(8, 8, -1), A, optimize=True).reshape(4096, -1)
    for i, factor in enumerate(selected):
        state = compiler._apply(state, factor, tuple(range(11-3*i, 8-3*i, -1)), 12)
    inverse = arrays['geometric_inverse']
    state = compiler._apply(state, inverse, tuple(range(11, 5, -1)), 12)
    state = compiler._apply(state, inverse, tuple(range(5, -1, -1)), 12)
    state = compiler._apply(state, arrays['geometric_pointer'][geometric_flag], (8, 7, 6, 2, 1, 0), 12)
    state = compiler._apply(state, arrays['pointer'][pointer], (11, 10, 9, 5, 4, 3), 12)
    return state[:, 0] if vector else state


def _native_costs(native_resources):
    result, macros = {}, {}
    for law in ('C', 'E'):
        code = transport.PAIR_CODES[law]
        counts = native_resources[law]
        result[law] = sum(counts[key]['native_pulses'] for key in ('encoding', 'decoding', 'repair'))+counts['recovery_max_native_pulses']
        encoder = len(native_routing.permutation_edges(code.encoder_permutation()))
        route = len(native_routing.permutation_edges(record_code._route_permutation(transport.SEAMS[law])))
        recovery = max(len(native_routing.permutation_edges(code.branch_permutation(law, transport.SEAMS[law], coin, raw)))
            for coin in (0, 1) for raw in circuit.record_values(law, coin)
            if code.branch_weight(law, transport.SEAMS[law], coin, raw))
        macros[law] = 3*encoder+2*route+recovery
    return result, macros


def _export_diagnostic(words, pulses):
    families = list(words)+[native_routing._canonical_word(), native_routing._canonical_word(True),
        [compiler.pulse('X', sign*math.pi/2, [0]) for sign in (-1, 1)],
        [compiler.pulse('Y', -math.pi/4, [0])]]
    families += [circuit.repair_word(law, coin) for law in ('C', 'E') for coin in (0, 1)]
    families += [record_code._mixing_word(count, inverse) for count in (1, 4, 5, 6) for inverse in (False, True)]
    menu = {(event['kind'], event['angle']) for word in families for event in compiler._validate(word)}
    leak, excess = 0., 0.
    for kind, angle in sorted(menu):
        full = compiler.native_pulse(kind, angle)
        d = 4 if kind == 'J' else 2
        if full.shape != (6, 6) or not np.isfinite(full).all():
            raise ArithmeticError('finite complete native pulse required')
        leak = max(leak, float(np.linalg.norm(full[d:, :d], 2)))
        excess = max(excess, float(np.linalg.norm(full[:d, :d], 2))-1)
    if max(leak, excess) > 2e-11:
        raise ArithmeticError('proper-code pulse complement exceeds admission')
    # Fixed ceilings make the retained budget portable and prevent an absolute
    # replay tolerance from silently zeroing these tiny diagnostic quantities.
    ceiling = 2e-11
    failure = math.sqrt(pulses)*ceiling*math.exp((pulses-1)*math.log1p(ceiling))
    # A valid eigenbasis choice within a repeated eigenspace can change the
    # emitted angles without changing the requested operator or pulse budget.
    # Inspect every fresh angle above, but do not retain its backend-dependent
    # distinct-angle count as an exact cross-platform budget claim.
    return {'complement_operator_norm_admission': ceiling,
        'success_norm_excess_admission': ceiling, 'maximum_history_pulses': pulses,
        'first_failure_amplitude_bound': failure,
        'argument': 'orthogonal retained first-export-failure position; sqrt(N)*L*(1+e)^(N-1), with complete prefixes and no postselection',
        'physical_or_outward_rounded_certificate': False}


def _schema(arrays):
    shapes = {'native_C': (8, 512, 8), 'native_E': (16, 512, 8),
        'decoded_C': (8, 512, 8), 'decoded_E': (16, 512, 8),
        'labels_C': (8, 2), 'labels_E': (16, 2), 'controls': (20, 8, 8),
        'bell_preparation': (64, 64), 'geometric_preparation': (64, 64),
        'geometric_inverse': (64, 64), 'pointer': (2, 64, 64), 'geometric_pointer': (2, 64, 64)}
    if type(arrays) is not dict or set(arrays) != set(shapes):
        raise ValueError('complete reference source array menu required')
    for name, shape in shapes.items():
        value = arrays[name]
        dtype = np.dtype('int64') if name.startswith('labels') else np.dtype(complex)
        if type(value) is not np.ndarray or value.shape != shape or value.dtype != dtype or not np.isfinite(value).all():
            raise ValueError('original finite source array schema differs: '+name)


def _resource_record(value):
    """Keep charged counts, omitting redundant backend residuals and QR zeros.

    Every zero-angle pulse remains charged in native_pulses. Its diagnostic
    classification can vary across numerical eigensolvers. The sole retained
    resource floats are integrated processor norms of the emitted words;
    proper-code numerical errors belong to the fixed global export admission.
    """
    if type(value) is dict:
        return {key: _resource_record(item) for key, item in value.items()
                if key not in ('maximum_numerical_valid_code_leakage', 'zero_angle_requests')}
    if type(value) is list:
        return [_resource_record(item) for item in value]
    return value


def metrics(arrays):
    _schema(arrays)
    targets = _targets()
    control_errors = [float(np.linalg.norm(a-b, 2)) for a, b in zip(arrays['controls'], targets, strict=True)]
    extension = transport.preparation_plan()[1]
    geometry_error = float(np.linalg.norm(arrays['geometric_preparation']-extension, 2))
    inverse_error = float(np.linalg.norm(arrays['geometric_inverse']-extension.conj().T, 2))
    phi = _reference_vector().astype(complex)
    bell_error = float(np.linalg.norm(arrays['bell_preparation'][:, 0]-phi))
    bell_gram = float(np.linalg.norm(arrays['bell_preparation'].conj().T@arrays['bell_preparation']-np.eye(64), 2))
    P = np.outer(phi, phi.conj())
    pointer_target = np.asarray([P, np.eye(64)-P])
    pointer_error = float(np.linalg.norm((arrays['pointer']-pointer_target).reshape(128, 64), 2))
    zero = np.zeros((64, 64), complex); zero[0, 0] = 1
    geometric_pointer_error = float(np.linalg.norm((arrays['geometric_pointer']-np.asarray(
        [zero, np.eye(64)-zero])).reshape(128, 64), 2))
    native_errors, native_norms, native_grams = {}, {}, {}
    for law in ('C', 'E'):
        code = transport.PAIR_CODES[law]
        labels = np.asarray([(coin, raw) for coin in (0, 1) for raw in range(1 << circuit.RECORD_BITS[law])], np.int64)
        if not np.array_equal(labels, arrays['labels_'+law]):
            raise ValueError('a native raw record label is missing or reordered')
        weights = [float(code.branch_weight(law, transport.SEAMS[law], *label)) for label in labels.tolist()]
        native_target = np.asarray([math.sqrt(w)*code.ideal_code_matrix() for w in weights])
        decoded_target = np.zeros_like(native_target)
        for i, w in enumerate(weights):
            decoded_target[i, :6] = math.sqrt(w)*np.eye(6)
        native_errors[law] = float(np.linalg.norm((arrays['decoded_'+law][:, :, :6]-decoded_target).reshape(-1, 6), 2))
        raw_error = float(np.linalg.norm((arrays['native_'+law][:, :, :6]-native_target).reshape(-1, 6), 2))
        stack = arrays['decoded_'+law].reshape(-1, 8)
        native_norms[law] = float(np.linalg.norm(stack, 2))
        native_grams[law] = float(np.linalg.norm(stack.conj().T@stack-np.eye(8), 2))
        if max(raw_error, native_errors[law], native_grams[law]) > SOURCE_ADMISSION:
            raise ValueError('complete native source instrument exceeds admission')
    if max(control_errors+[geometry_error, inverse_error, bell_error, bell_gram, pointer_error, geometric_pointer_error]) > SOURCE_ADMISSION:
        raise ValueError('source reference preparation, current or pointer exceeds independent target admission')
    # Two copies of each geometric/current map, four local native factors.
    native_product_error = math.prod(1+native_errors[law] for law in LAWS)-1
    native_product_norm = math.prod(native_norms[law] for law in LAWS)
    prefix_error = (1+bell_error)*(1+max(control_errors))**2*(1+geometry_error)**2-1
    after_native = native_product_norm*prefix_error+native_product_error
    successful_error = (1+pointer_error)*(1+geometric_pointer_error)*(1+inverse_error)**2*(1+after_native)-1
    survivals, ideal_survivals = [], []
    for source, ideal in zip(arrays['controls'], targets, strict=True):
        moved = np.kron(source, source)@arrays['bell_preparation'][:, 0]
        survivals.append(float(np.linalg.norm(arrays['pointer'][0]@moved)**2))
        ideal_survivals.append(float(abs(np.vdot(phi, np.kron(ideal, ideal)@phi))**2))
    if max(abs(a-b) for a, b in zip(survivals, ideal_survivals, strict=True)) > 12*SOURCE_ADMISSION:
        raise ValueError('the actual pointer did not measure the independent reference survival')
    if min(survivals[:15]) < 1-12*SOURCE_ADMISSION or max(survivals[15:]) > .95:
        raise ValueError('compatible and excluded current reference witnesses are not resolved')
    # The same source current also acts nontrivially on a cross-triplet probe.
    marker = np.zeros(8, complex); marker[0] = marker[3] = 1/math.sqrt(2)
    class_marker = float(abs(np.vdot(marker, arrays['controls'][13]@marker))**2)
    if abs(class_marker-math.cos(math.sqrt(5)/8)**2) > SOURCE_ADMISSION:
        raise ValueError('the promoted class sum lacks its relative central action')
    # Actual complete output-amplitude checks supplement the all-input tensor
    # bound. Complement flags remain separate outcomes with zero ideal target.
    records, weight = [], 1.
    for law in LAWS:
        code = transport.PAIR_CODES[law]
        record = next((c, r) for c in (0, 1) for r in range(1 << circuit.RECORD_BITS[law])
                      if code.branch_weight(law, transport.SEAMS[law], c, r))
        records.append(record)
        weight *= float(code.branch_weight(law, transport.SEAMS[law], *record))
    output_errors = []
    for control in (0, 1, 13, 15, 19):
        moved = np.kron(targets[control], targets[control])@phi
        first = phi*np.vdot(phi, moved)
        for flag in (0, 1):
            for pointer, target_block in enumerate((first, moved-first)):
                actual = contract_output(arrays, control, tuple(records), (0,)*4, pointer, geometric_flag=flag)
                target = np.zeros(4096, complex)
                if flag == 0:
                    addresses = np.asarray([512*a+8*b for a in range(8) for b in range(8)])
                    target[addresses] = math.sqrt(weight)*target_block
                output_errors.append(float(np.linalg.norm(actual-target)))
    if max(output_errors) > SOURCE_ADMISSION:
        raise ValueError('complete native reference output amplitudes differ from the independent tensor target')
    return {'full_source_control_dimension': 8, 'source_control_count': len(targets),
        'control_operator_errors': control_errors, 'geometric_preparation_operator_error': geometry_error,
        'geometric_inverse_operator_error': inverse_error, 'Bell_preparation_error': bell_error,
        'Bell_complete_Gram_defect': bell_gram, 'pointer_operator_error': pointer_error,
        'geometric_pointer_operator_error': geometric_pointer_error,
        'native_complete_dilation_errors': native_errors, 'native_full8_input_Gram_defects': native_grams,
        'native_tensor_dilation_error': native_product_error,
        'native_tensor_full_input_norm': native_product_norm,
        'successful_factorized_history_operator_error': successful_error,
        'source_reference_survival_probabilities': survivals,
        'independent_reference_survival_probabilities': ideal_survivals,
        'class_sum_cross_triplet_marker_probability': class_marker,
        'class_sum_cross_triplet_marker_probability_exact': 'cos(sqrt(5)/8)^2',
        'selected_complete_output_amplitude_errors': output_errors,
        'selected_output_controls': [0, 1, 13, 15, 19],
        'selected_output_raw_records': [list(record) for record in records],
        'error_domain': 'fixed prepared reference, or replace its preparation error by zero for arbitrary exact inputs in geometric6 tensor geometric6 with untouched reference columns',
        'native_high_syndromes_classically_measured_or_postselected': False}


def execute():
    exact_certificate()
    arrays, native_resources = current._native_factors()
    menu = _control_menu()
    bell_word = compiler.compile_real_state(_reference_vector())
    preparation = [event for stage in transport._preparation_stages() for event in stage]
    inverse = compiler.inverse(preparation)
    arrays['controls'] = np.asarray([compiler.execute(row[3], 3) for row in menu])
    arrays['bell_preparation'] = compiler.execute(bell_word, 6)
    arrays['geometric_preparation'] = compiler.execute(preparation, 6)
    arrays['geometric_inverse'] = compiler.execute(inverse, 6)
    detector = history.projector_read(np.eye(64, dtype=complex), _reference_vector())
    if detector.failure_outputs:
        raise ArithmeticError('six-bit pointer added unregistered scratch outcomes')
    arrays['pointer'] = np.asarray(detector.outputs)
    zero = np.zeros(64); zero[0] = 1
    geometric_detector = history.diagonal_read(np.eye(64, dtype=complex), zero)
    if geometric_detector.failure_outputs:
        raise ArithmeticError('geometric-return pointer added unregistered scratch outcomes')
    arrays['geometric_pointer'] = np.asarray(geometric_detector.outputs)
    numerical = metrics(arrays)
    native_cost, macros = _native_costs(native_resources)
    control_costs = [compiler.resources(row[3])['native_pulses'] for row in menu]
    pulse_count = (compiler.resources(bell_word)['native_pulses']+2*max(control_costs)
        +2*compiler.resources(preparation)['native_pulses']+2*compiler.resources(inverse)['native_pulses']
        +sum(native_cost[law] for law in LAWS)+detector.resources['native_pulses']
        +geometric_detector.resources['native_pulses'])
    macro_count = sum(macros[law] for law in LAWS)
    projection = 5.1e-11
    if any(native_routing._executed_macro(direction)[1]['projection_change'] > projection for direction in (False, True)):
        raise ArithmeticError('raw native macro correction exceeds its fixed ceiling')
    raw_adjustment = math.exp(pulse_count*math.log1p(2e-11))*math.expm1(macro_count*math.log1p(projection))
    export = _export_diagnostic([bell_word, preparation, inverse, list(detector.word), list(geometric_detector.word)]
                               +[row[3] for row in menu], pulse_count)
    ceiling = SOURCE_ADMISSION
    native_norm = (1+ceiling)**2  # Four local Gram bounds give this norm.
    prefix = (1+ceiling)**5-1   # Bell, two currents, two geometric embeddings.
    after_native = native_norm*prefix+(1+ceiling)**4-1
    fixed_successful = (1+ceiling)**4*(1+after_native)-1  # Two inverses, two pointers.
    successful = (1+fixed_successful)*(1+raw_adjustment)-1
    complete = successful+export['first_failure_amplitude_bound']
    if not math.isfinite(complete) or complete > 1e-5:
        raise ArithmeticError('complete factorized reference history exceeds numerical admission')
    services = {'initial_blank_work_qubits': 36, 'native_attempt_coin_and_record_initializations': 8,
        'native_record_bits': 14, 'terminal_pointer_bits': 2,
        'first_export_failure_position_bits': pulse_count.bit_length(),
        'first_export_failure_position_commit': 1,
        'maximum_archive_bits_including_export_position': 16+pulse_count.bit_length(),
        'coherent_high_decoder_syndrome_qubits_retained': 24,
        'geometric_blank_qubits_retained': 6,
        'proper_code_loads': pulse_count, 'proper_code_exports': pulse_count,
        'clock': 'one declared dimensionless coherent parameter per control; source ordinals are not physical duration',
        'physical_error_sum': 'sum_j count_j*delta_j for independent supplied preparation, current, load/export, record commit, storage and pointer contracts; no values are calibrated here'}
    report = {'schema': 'oph.native_response_reference.v1', 'exact': exact_certificate(),
        'control_menu': [{'name': name, 'parameter': str(parameter), 'generator_exact': _entries(generator),
                          'source_native_pulses': count} for (name, parameter, generator, _), count in zip(menu, control_costs, strict=True)],
        'execution': numerical,
        'retained_output_contract': {'native_carrier_order': list(LAWS), 'raw_records': 8*16*8*16,
            'pointer_outcomes_per_raw_record': 4, 'native_rows_per_carrier': 512,
            'coherent_high_syndrome_dimension_per_carrier': 64,
            'output_rows_per_syndrome_block': 4096,
            'formula': 'M_pointer F_geometric_return (G_inverse tensor G_inverse) (D_C[r0,h0] tensor D_E[r1,h1] tensor D_C[r2,h2] tensor D_E[r3,h3]) (A tensor A) (U tensor U) B|0>',
            'contraction_api': 'contract_output(arrays, control, records, syndromes, pointer, input_columns=None, geometric_flag=0)',
            'survival_event': 'geometric_flag0 AND Bell_pointer0; all complementary flags/padding outcomes remain outputs',
            'terminal_policy': 'all four native attempts and terminal pointer run for every raw label; aborts select no recovery, remain outputs, and are never postselected',
            'proper_code_export_policy': 'first failure halts, retaining its pulse ordinal; its numerical complement is charged separately',
            'ideal_scalar_Kraus_replacements_used_in_source_contraction': False,
            'entire_512_power4_array_materialized': False},
        'resources': {'native_pulse_maximum_per_history': pulse_count, 'native_pulses_per_local_law': native_cost,
            'native_source_macro_invocations': macro_count, 'native_factors': native_resources,
            'Bell_preparation': compiler.resources(bell_word), 'geometric_preparation_per_object': compiler.resources(preparation),
            'geometric_inverse_per_object': compiler.resources(inverse), 'pointer': detector.resources,
            'geometric_return_pointer': geometric_detector.resources,
            'declared_services': services},
        'error_budget': {'component_source_admission': SOURCE_ADMISSION,
            'fixed_admission_successful_block_error': fixed_successful,
            'raw_macro_projection_change_ceiling_per_invocation': projection, 'raw_macro_adjustment': raw_adjustment,
            'raw_macro_adjustment_formula': '(1+2e-11)^N * ((1+5.1e-11)^M-1), retaining intervening success-block norm inflation',
            'successful_history_including_raw_macro_error': successful, 'proper_code_exports': export,
            'complete_flagged_operator_error': complete, 'complete_flagged_channel_error': (2+complete)*complete,
            'numerical_not_outward_rounded_or_physical': True,
            'references_are_retained_by_operator_norm_bound': True},
        'scope': {'reference_and_class_sum_are_supplied': True, 'reference_compatibility_is_a_maximal_testable_constraint': True,
            'native_exchange_selects_this_reference_or_constraint': False,
            'signed_native_exchange_words_implement_all_controls_in_this_experiment': False,
            'controls_use_existing_declared_coherent_apparatus': True,
            'laboratory_clock_or_SM_gauge_redundancy_derived': False,
            'same_reference_is_the_existing_mixed_native_ground_ray': False}}
    report['resources'] = _resource_record(report['resources'])
    return report, arrays


def compare_execution(report, arrays, fresh_report, fresh_arrays):
    own = metrics(arrays)
    _owned_compare(own, report['execution'])
    # Fixed admissions/count-derived budgets have no absolute replay slack.
    _fixed_compare(fresh_report['error_budget'], report['error_budget'])
    from copy import deepcopy
    normalized = deepcopy(report)
    normalized['execution'] = fresh_report['execution']
    normalized['error_budget'] = fresh_report['error_budget']
    current._compare(fresh_report, normalized)
    _schema(fresh_arrays)
    for key, actual in arrays.items():
        expected = fresh_arrays[key]
        if key.startswith('labels'):
            if not np.array_equal(actual, expected):
                raise ValueError('native reference record labels differ')
        elif np.linalg.norm(actual-expected) > REPLAY_ADMISSION:
            raise ValueError('complete native reference source array differs: '+key)


def _owned_compare(expected, actual):
    if type(expected) is not type(actual):
        raise ValueError('own-array diagnostic original type differs')
    if type(expected) is dict:
        if expected.keys() != actual.keys():
            raise ValueError('own-array diagnostic keys differ')
        for key in expected:
            _owned_compare(expected[key], actual[key])
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise ValueError('own-array diagnostic list differs')
        for a, b in zip(expected, actual, strict=True):
            _owned_compare(a, b)
    elif type(expected) is float:
        scale = max(1., abs(expected), abs(actual))
        if not math.isfinite(actual) or abs(expected-actual) > 32*np.finfo(float).eps*scale:
            raise ValueError('own-array diagnostic differs at machine scale')
    elif expected != actual:
        raise ValueError('own-array diagnostic exact value differs')


def _fixed_compare(expected, actual):
    if type(expected) is not type(actual):
        raise ValueError('fixed budget original type differs')
    if type(expected) is dict:
        if expected.keys() != actual.keys():
            raise ValueError('fixed budget keys differ')
        for key in expected:
            _fixed_compare(expected[key], actual[key])
    elif type(expected) is float:
        if not math.isfinite(actual) or not math.isclose(expected, actual, rel_tol=8e-15, abs_tol=0.):
            raise ValueError('fixed admission budget differs')
    elif expected != actual:
        raise ValueError('fixed admission budget exact value differs')


def source_bindings():
    paths = set(comparison.source_bindings())|set(current.source_bindings())|{
        'model/common_source/native_response_reference.py', 'model/common_source/test_native_response_reference.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}
