"""Operational family-band selection on a complete protected native pair.

Incidence selects a twelve-rotation class. Its complete36-state spectrum
selects the highest non-invariant band. A declared covariant cooling service
preserves this band and resets its complement with every outcome retained.
Its dimension is computed before comparison with the finite family theorem.
The neutral family pair is separate from the pair carrying gauge currents.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np
import sympy as s

from . import native_gauge_object as objects
from . import native_code_transport as transport
from .mass_response import history
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler

ROOT = Path(__file__).resolve().parents[2]
SOURCE_ADMISSION = 2e-8
ARRAY_ADMISSION = 2e-9
ARRAY_SCHEMA = {
    'family_projector': ((36, 36), 'complex128'),
    'family_port_frame': ((36, 12), 'complex128'),
    'family_coefficient_isometry': ((36, 12), 'complex128'),
    'family_cooling_labels': ((605, 3), 'int64'),
    'family_band_pointer': ((2, 64, 64), 'complex128'),
    'family_ray_preparations': ((12, 64, 64), 'complex128'),
    'family_selector_outputs': ((16, 16), 'complex128'),
    'family_coordinate_local': ((2, 2, 2), 'complex128'),
    'family_flip_local': ((2, 2), 'complex128'),
    'family_C': ((8, 512, 8), 'complex128'),
    'family_E': ((16, 512, 8), 'complex128'),
    'family_C_labels': ((8, 2), 'int64'),
    'family_E_labels': ((16, 2), 'int64')}


def _clean(matrix): return objects._clean(matrix)
def _zero(matrix): return not _clean(matrix).todok()
def _entries(matrix): return [[str(x) for x in row] for row in _clean(matrix).tolist()]


def _projector(operator, eigenvalue, roots):
    P = s.eye(operator.rows)
    for value in roots:
        if value != eigenvalue:
            P = P*(operator-value*s.eye(operator.rows))/(eigenvalue-value)
    return _clean(P)


@lru_cache(maxsize=1)
def geometry():
    data = objects.geometry(); rotations = tuple(sorted(data['rotations']))
    edges = set(data['seams']); neighbors = [set() for _ in range(12)]
    for i, j in edges: neighbors[i].add(j); neighbors[j].add(i)
    axes = []
    for i in range(12):
        distance = {i: 0}; front = {i}
        while front:
            following = set().union(*(neighbors[j] for j in front))-set(distance)
            for j in following: distance[j] = distance[next(iter(front))]+1
            front = following
        axes.extend((i, j) for j in range(i+1, 12) if distance[j] == 3)
    axes = tuple(axes)
    if axes != data['axes']:
        raise ArithmeticError('incidence antipodes differ from the native code axes')
    index = {pair: i for i, pair in enumerate(axes)}
    axis_actions = tuple(tuple(index[tuple(sorted((g[i], g[j])))] for i, j in axes) for g in rotations)
    pair_actions = tuple(tuple(6*a[j//6]+a[j % 6] for j in range(36)) for a in axis_actions)
    near = tuple(k for k, g in enumerate(rotations) if objects.permutation_order(g) == 5
                 and sum(tuple(sorted((i, g[i]))) in edges for i in range(12)) == 10)
    if len(near) != 12: raise ArithmeticError('incidence-defined near-fivefold class changed')
    matrices = tuple(objects.permutation_matrix(p) for p in pair_actions)
    C = s.ImmutableMatrix(sum((matrices[k] for k in near), s.zeros(36)))
    fixed = s.ImmutableMatrix(sum(matrices, s.zeros(36))/len(matrices))
    x = s.Symbol('x'); polynomial = C.charpoly(x).as_expr(); roots = s.roots(polynomial, x)
    if sum(roots.values()) != 36 or not _zero(C*fixed-len(near)*fixed):
        raise ArithmeticError('the complete class spectrum or invariant sector is missing')
    selected = max((z for z in roots if z != len(near)), key=lambda z: float(z))
    if not all(z == selected or (selected-z).is_positive is True for z in roots if z != len(near)):
        raise ArithmeticError('the largest non-invariant eigenvalue is not certified')
    projectors = {value: _projector(C, value, roots) for value in roots}
    if not _zero(sum(projectors.values(), s.zeros(36))-s.eye(36)):
        raise ArithmeticError('a full native sector was discarded')
    for value, P in projectors.items():
        if not _zero(P*P-P) or not _zero(C*P-value*P) or s.trace(P) != roots[value]:
            raise ArithmeticError('complete spectral projector certificate failed')
    P = projectors[selected]; rank = int(s.trace(P))

    def frame(value):
        # Reynolds average of a matrix unit; no geometric coordinate frame
        # or assumed family dimension enters the producer.
        projector = projectors[value]
        for seed in range(36):
            R = s.zeros(36, 12)
            for g, action in zip(rotations, pair_actions, strict=True): R[action[seed], g[0]] += 1
            R = _clean(projector*R)
            if any(R): break
        else: raise ArithmeticError('native port coefficients do not realize the selected band')
        scale = s.simplify(s.trace(R*R.T)/s.trace(projector))
        A = _clean(R/s.sqrt(scale))
        if not _zero(A*A.T-projector): raise ArithmeticError('Reynolds frame is not tight')
        return A, seed, scale

    A, seed, scale = frame(selected)
    C12 = sum((objects.permutation_matrix(rotations[k]) for k in near), s.zeros(12))
    adjacency = s.zeros(12)
    for i, j in edges: adjacency[i, j] = adjacency[j, i] = 1
    if C12 != 2*adjacency+2*s.eye(12):
        raise ArithmeticError('the native coefficient operator is not the incidence operator')
    coefficient_roots = s.roots(C12.charpoly(x).as_expr(), x)
    other = next(value for value in coefficient_roots if value not in (0, len(near), selected))
    A_other, _, _ = frame(other)
    diagonal = s.zeros(36, 12)
    for color, pair in enumerate(axes):
        for port in pair: diagonal[7*color, port] = 1/s.sqrt(2)
    W = _clean(diagonal+A+A_other)
    if not _zero(W.T*W-s.eye(12)) or not _zero(C*W-W*C12):
        raise ArithmeticError('the complete12-coefficient native isometry failed')
    if not _zero(W*_projector(C12, selected, coefficient_roots)*W.T-P):
        raise ArithmeticError('coefficient-band selection differs from the operational band')
    for g, R in zip(rotations, matrices, strict=True):
        port = objects.permutation_matrix(g)
        if not _zero(R*A-A*port) or not _zero(R*W-W*port):
            raise ArithmeticError('native state or recorded cooling covariance failed')
    return {'rotations': rotations, 'axis_actions': axis_actions, 'pair_actions': pair_actions,
            'near_indices': near, 'C': C, 'fixed': fixed, 'polynomial': s.factor(polynomial),
            'roots': roots, 'projectors': projectors, 'selected': selected, 'P': P,
            'Q': s.eye(36)-P, 'rank': rank, 'frame': A, 'seed': seed, 'scale': scale,
            'W': W, 'C12': C12, 'coefficient_roots': coefficient_roots}


def cooling_maps():
    """All logical and terminal-padding Kraus maps on the full64-state input."""
    d = geometry(); addresses = transport.product_addresses()
    P = np.zeros((64, 64), complex); Q = np.zeros((64, 64), complex)
    A = np.zeros((64, 12), complex)
    P[np.ix_(addresses, addresses)] = np.asarray(d['P'], complex)
    Q[np.ix_(addresses, addresses)] = np.asarray(d['Q'], complex)
    A[addresses] = np.asarray(d['frame'], complex)
    maps, labels = [P], [(0, -1, -1)]
    for port in range(12):
        for color, address in enumerate(addresses):
            maps.append(np.outer(A[:, port], Q[address])/math.sqrt(d['rank']))
            labels.append((1, port, color))
    for port in range(12, 16):
        for color in range(36):
            maps.append(np.zeros((64, 64), complex)); labels.append((3, port, color))
    for address in sorted(set(range(64))-set(addresses)):
        value = np.zeros((64, 64), complex); value[address, address] = 1
        maps.append(value); labels.append((2, address, -1))
    return np.asarray(maps), np.asarray(labels, np.int64)


def apply_cooling(rho):
    """Complete36-state unread map; no conditioning on the selected outcome."""
    if type(rho) is not np.ndarray or rho.shape != (36, 36) or not np.isfinite(rho).all():
        raise ValueError('a complete finite36-state operator is required')
    d = geometry(); P, Q = np.asarray(d['P'], complex), np.asarray(d['Q'], complex)
    return P@rho@P+np.trace(Q@rho)*P/d['rank']


def adjoint_cooling(observable):
    d = geometry(); P, Q = np.asarray(d['P'], complex), np.asarray(d['Q'], complex)
    return P@observable@P+np.trace(P@observable)*Q/d['rank']


def exact_report():
    d = geometry(); P, Q, A, rank = d['P'], d['Q'], d['frame'], d['rank']
    # Completeness of every raw map follows from this exact full-space Gram,
    # rather than checking one initialized family state.
    gram = P.T*P+s.trace(A.T*A)/rank*Q.T*Q
    if not _zero(gram-s.eye(36)) or not _zero(P*Q):
        raise ArithmeticError('the complete cooling instrument lost probability')
    B = objects.geometry()['frame']; metric = (B.T*B)[0, 0]
    if not _zero(A.T*A-B*B.T/metric):
        raise ArithmeticError('the produced coefficient band differs from the finite theorem band')
    return {'classification': 'specified recorded cooling selects a computed persistent family band',
        'architecture': {'native_twelve_port_incidence_unchanged': True,
            'family_pair_is_separate_from_gauge_current_pair': True,
            'gauge_action_on_family_pair': 'identity; neutral ancillary family object',
            'near_class_rule': 'order5 rotations carrying ten ports to adjacent ports',
            'preserved_band_rule': 'largest class-sum eigenvalue outside the invariant sector',
            'cooling_service_is_an_added_coherent_repair_choice': True,
            'unchanged_C_E_selects_the_band': False,
            'physical_family_role_requires_the_declared_matter_tensor_attachment': True},
        'complete_native_spectrum': {'dimension': 36, 'near_class_size': len(d['near_indices']),
            'characteristic_polynomial': str(d['polynomial']),
            'eigenvalues_and_dimensions': [[str(value), int(n)] for value, n in sorted(d['roots'].items(), key=lambda row: float(row[0]))],
            'invariant_sector_dimension': int(s.trace(d['fixed'])),
            'selected_eigenvalue': str(d['selected']), 'computed_persistent_dimension': rank,
            'retained_complement_dimension': 36-rank},
        'coefficient_realization': {'isometry36_by12_exact': _entries(d['W']),
            'all_twelve_coefficient_inputs_realized': True, 'remaining_native_dimensions_retained': 24,
            'native_class_operator_equals_two_adjacency_plus_two_identity': True,
            'coefficient_band_dimensions': sorted(map(int, d['coefficient_roots'].values())),
            'native_coefficient_projector_equals_tested_finite_family_band': True,
            'family_dimension_inserted_into_Reynolds_producer': False,
            'reynolds_seed_port': 0, 'reynolds_seed_color_pair_index': d['seed'],
            'reynolds_normalization_squared': str(d['scale'])},
        'complete_recorded_cooling': {'preserve_map': 'P',
            'reset_maps': '|a_port><color|Q/sqrt(tr(P)); all12ports and36colors',
            'ideal_logical_outcomes': 433, 'complete_source_logical_outcomes': 577,
            'terminal_padding_outcomes': 28,
            'complete_register_dimension': 64, 'all_register_outcomes': 605,
            'selector_abort_outcomes_zero_in_ideal_limit': 144,
            'full36_completeness_exact': True, 'all60_classical_record_covariance_actions_exact': True,
            'record_action': 'port->g(port), color_pair->axis(g) tensor axis(g)',
            'unread_map': 'P rho P+Tr(Q rho)P/tr(P)',
            'postselection_or_hidden_complement_removal': False,
            'reset_records_retain_the_removed_input_color_label': True},
        'persistence': {'Schrodinger_fixed_states': 'all density matrices supported on P',
            'effective_protected_output_algebra': 'End(range(P))',
            'protected_output_complex_operator_dimension': rank*rank,
            'Heisenberg_fixed_operator_system': 'B_P+(Tr(B_P)/tr(P))*Q',
            'full_Heisenberg_fixed_system_closed_under_ordinary_product': False,
            'coherent_superpositions_and_external_references_within_P_preserved': True,
            'complement_removed_from_full_record_accounting': False,
            'terminal_padding_flags_are_continuing_family_states': False,
            'terminal_selector_abort_flags_are_continuing_family_states': False,
            'decay_eigenvalues_identified_as_particle_masses': False},
        'empirical_inputs': [],
        'scope': {'operational_coefficient_band_on_actual_correctable_pair': True,
            'source_pulses_for_cooling_instrument_executed_on_unknown64_inputs': True,
            'coherent_compiler_pointer_and_blank_services_are_declared': True,
            'physical_family_count_forced_by_incidence_without_this_repair_choice': False,
            'particle_mass_spectrum_or_fermion_statistics_derived': False}}


def _native_factors():
    arrays = {}
    for law, seam in (('C', (2, 3)), ('E', (2, 10))):
        local = transport.PAIR_CODES[law]
        encoded = local._encode(np.eye(512, 8, dtype=complex))
        branches = local.repair_branches(law, seam, encoded)
        arrays['family_'+law] = np.asarray([local.decoding(row['amplitude']) for row in branches])
        arrays['family_'+law+'_labels'] = np.asarray([(row['coin'], row['record_code']) for row in branches], np.int64)
    return arrays


def _bit_instrument(record, local, bits=6):
    result = np.asarray([[1]], complex)
    for bit in reversed(range(bits)):
        result = np.kron(result, local[(record >> bit) & 1])
    return result


def _zero_shift(address, flip):
    return _bit_instrument(address, np.asarray([np.eye(2), flip]))


def _cooling_source():
    """Compile all factors, preserving their action on every64 input column."""
    d = geometry(); addresses = transport.product_addresses()
    preparation = [event for stage in transport._preparation_stages() for event in stage]
    effect = np.zeros(64); effect[8*np.arange(d['rank'])] = 1
    pointer = history.diagonal_read(np.eye(64, dtype=complex), effect,
        basis_word=compiler.inverse(preparation))
    arrays = {'family_band_pointer': np.asarray(pointer.outputs)}
    source_rays, resources = [], []
    for port in range(12):
        vector = np.zeros(64)
        vector[addresses] = np.asarray(d['frame'][:, port], float).reshape(36)/math.sqrt(d['rank']/12)
        word = compiler._compile_real_state_tree(vector, 6)
        source_rays.append(compiler.execute(word, 6, np.eye(64, dtype=complex)))
        resources.append(compiler.resources(word))
    arrays['family_ray_preparations'] = np.asarray(source_rays)
    coordinate = history.diagonal_read(np.eye(2, dtype=complex), [1, 0])
    arrays['family_coordinate_local'] = np.asarray(coordinate.outputs)
    flip_word = [compiler.pulse('I', -math.pi/2, [0]), compiler.pulse('X', math.pi/2, [0])]
    arrays['family_flip_local'] = compiler.execute(flip_word, 1, np.eye(2, dtype=complex))
    vector = np.zeros(16); vector[:12] = 1/math.sqrt(12)
    selector_word = compiler._compile_real_state_tree(vector, 4)
    selector = compiler.execute(selector_word, 4, np.eye(16, dtype=complex)[:, :1])
    arrays['family_selector_outputs'] = np.asarray([
        (_bit_instrument(i, arrays['family_coordinate_local'], 4)@selector).reshape(16)
        for i in range(16)])
    return arrays, {'band_pointer': pointer.resources,
        'each_full64_ray_preparation': resources,
        'one_bit_coordinate_pointer': coordinate.resources,
        'selector_preparation': compiler.resources(selector_word),
        'one_bit_reset_flip': compiler.resources(flip_word),
        'source_basis_columns_for_band_read_and_each_ray_unitary': 64,
        'coordinate_pointer_repetitions_per_Q_branch': 6,
        'selector_pointer_repetitions_per_valid_coordinate_branch': 4,
        'no_target_vector_replaces_the_finite_source_postmeasurement_state': True}


def cooling_gram(arrays):
    """Sum every actual cooling branch Gram by exact factor contraction."""
    G0, G1 = arrays['family_band_pointer']; local = arrays['family_coordinate_local']
    weights = np.sum(abs(arrays['family_selector_outputs'])**2, axis=1)
    reset_gram = sum((weights[i]*U.conj().T@U
        for i, U in enumerate(arrays['family_ray_preparations'])), np.zeros((64, 64), complex))
    total = G0.conj().T@G0; valid = set(transport.product_addresses())
    for j in range(64):
        branch = _bit_instrument(j, local)@G1
        if j in valid:
            flip = _zero_shift(j, arrays['family_flip_local'])
            middle = flip.conj().T@reset_gram@flip+float(weights[12:].sum())*np.eye(64)
            total += branch.conj().T@middle@branch
        else:
            total += branch.conj().T@branch
    return total


def cooling_metrics(arrays):
    """Full-source dilation comparison, including selector and padding failures."""
    d = geometry(); addresses = transport.product_addresses()
    ray_defect = max(float(np.linalg.norm(U.conj().T@U-np.eye(64), 2))
                     for U in arrays['family_ray_preparations'])
    if ray_defect > 4*SOURCE_ADMISSION:
        raise ValueError('full source ray preparation operator lost unitarity')
    P = np.zeros((64, 64), complex); Q = np.zeros((64, 64), complex)
    P[np.ix_(addresses, addresses)] = np.asarray(d['P'], complex)
    Q[np.ix_(addresses, addresses)] = np.asarray(d['Q'], complex)
    G0, G1 = arrays['family_band_pointer']; local = arrays['family_coordinate_local']
    squared = float(np.linalg.norm(G0-P)**2); padding = set(range(64))-set(addresses)
    selector_abort = 0.
    for j in range(64):
        branch = _bit_instrument(j, local)@G1
        if j in padding:
            target = np.zeros((64, 64), complex); target[j, j] = 1
            squared += float(np.linalg.norm(branch-target)**2)
            continue
        shifted = _zero_shift(j, arrays['family_flip_local'])@branch
        for port in range(16):
            state = arrays['family_selector_outputs'][port]
            work = arrays['family_ray_preparations'][port]@shifted if port < 12 else branch
            difference = state[:, None, None]*work
            if port < 12:
                vector = np.zeros(64, complex)
                vector[addresses] = np.asarray(d['frame'][:, port], complex).reshape(36)/math.sqrt(d['rank'])
                difference[port] -= np.outer(vector, Q[j])
            else:
                selector_abort += float(np.vdot(difference, difference).real)
            squared += float(np.vdot(difference, difference).real)
    delta = math.sqrt(squared); gram = cooling_gram(arrays)
    defect = float(np.linalg.norm(gram-np.eye(64), 2))
    if delta > SOURCE_ADMISSION or defect > 4*SOURCE_ADMISSION:
        raise ValueError('complete compiled cooling exceeds unknown-input/reference admission')
    return {'source_cooling_dilation_error_diagnostic': delta,
        'full64_ray_preparation_unitarity_defect_diagnostic': ray_defect,
        'full64_source_cooling_Gram_defect_diagnostic': defect,
        'selector_abort_full_input_Frobenius_probability_diagnostic': selector_abort,
        'cooling_dilation_allowance': SOURCE_ADMISSION,
        'cooling_channel_allowance': (2+SOURCE_ADMISSION)*SOURCE_ADMISSION,
        'all64_unknown_input_and_reference_columns_executed': True,
        'all605_raw_cooling_records_retained': True,
        'nonzero_source_coordinate_read_leakage_is_propagated': True,
        'all16_selector_outputs_and_quantum_rows_retained': True,
        'target_rays_inserted_in_place_of_source_preparations': False}


def validate_arrays(arrays):
    if set(arrays) != set(ARRAY_SCHEMA): raise ValueError('complete family-selection array menu required')
    for name, (shape, dtype) in ARRAY_SCHEMA.items():
        value = arrays[name]
        if (type(value) is not np.ndarray or value.shape != shape or value.dtype != np.dtype(dtype)
                or not np.isfinite(value).all()): raise ValueError('invalid family-selection array '+name)
    d = geometry()
    for name, value in (('family_projector', d['P']), ('family_port_frame', d['frame']), ('family_coefficient_isometry', d['W'])):
        if np.max(abs(arrays[name]-np.asarray(value, complex))) > 1e-12:
            raise ValueError('source-derived family map differs: '+name)
    K, labels = cooling_maps()
    if not np.array_equal(arrays['family_cooling_labels'], labels): raise ValueError('cooling outcome labels differ')
    gram = np.einsum('rij,rik->jk', K.conj(), K, optimize=True)
    if np.linalg.norm(gram-np.eye(64)) > 1e-11: raise ValueError('full cooling register is not trace preserving')


def native_metrics(arrays):
    validate_arrays(arrays)
    errors, grams, raw_weights = [], [], []
    for law, seam in (('C', (2, 3)), ('E', (2, 10))):
        values = arrays['family_'+law]; labels = arrays['family_'+law+'_labels']
        expected_labels = np.asarray([(c, r) for c in (0, 1) for r in range(1 << circuit.RECORD_BITS[law])], np.int64)
        if not np.array_equal(labels, expected_labels): raise ValueError('native raw labels changed')
        target = np.zeros((len(labels), 512, 6), complex)
        weights = []
        for n, (coin, raw) in enumerate(labels.tolist()):
            weight = transport.PAIR_CODES[law].branch_weight(law, seam, coin, raw)
            target[n, :6] = math.sqrt(float(weight))*np.eye(6)
            weights.append(weight)
        error = float(np.linalg.norm(values[:, :, :6]-target))
        if error > SOURCE_ADMISSION: raise ValueError('complete native family source exceeds admission')
        errors.append(error); raw_weights.append(weights)
        grams.append(np.einsum('raj,rak->jk', values.conj(), values, optimize=True))
    full_gram = np.kron(*grams); addresses = transport.product_addresses()
    full_error = float(np.linalg.norm(full_gram-np.eye(64), 2))
    if full_error > 4*SOURCE_ADMISSION: raise ValueError('native unknown-input/reference completeness failed')
    dilation = errors[0]*(1+errors[1])+errors[1]
    cooling = cooling_metrics(arrays)
    native_allowance = SOURCE_ADMISSION*(1+SOURCE_ADMISSION)+SOURCE_ADMISSION
    allowance = SOURCE_ADMISSION*(1+native_allowance)+native_allowance
    local_channels = []
    for law in ('C', 'E'):
        values = arrays['family_'+law].reshape(-1, 64, 8, 8)
        local_channels.append(np.einsum('rhia,rhjc->ijac', values.conj(), values, optimize=True))
    gram_cooling = cooling_gram(arrays).reshape(8, 8, 8, 8)
    complete_gram = np.einsum('ijac,ikjl,klbd->abcd', local_channels[0], gram_cooling,
        local_channels[1], optimize=True).reshape(64, 64)
    if np.linalg.norm(complete_gram-np.eye(64), 2) > 8*SOURCE_ADMISSION:
        raise ValueError('complete native-plus-cooling instrument lost an input/reference direction')
    return {'local_native_rows_retained': 512, 'local_decoded_input_columns_retained': 8,
        'local_high_syndrome_values_retained': 64, 'complete_native_joint_raw_records': 128,
        'decoder_syndromes_are_measured_classical_record_labels': False,
        'joint_decoder_syndrome_quantum_output_dimension': 4096,
        'unknown_native_logical_inputs': 36, 'full_encoded_register_inputs': 64,
        'ancillary_quantum_reference_dimension_sufficient': 36,
        'native_joint_weights_exact': [str(a*b) for a in raw_weights[0] for b in raw_weights[1]],
        'full_logical_native_plus_cooling_records': 128*577,
        'full_register_native_plus_cooling_records': 128*605,
        'classical_label_count_definition': '128 raw native C/E labels times605 cooling labels; both64-valued decoder syndromes remain quantum outputs',
        'joint_output_definition': 'decoded local512-row factors tensor; cooling acts on both low3-bit registers, identity on every retained high-syndrome row',
        'all_inputs_and_all_raw_records_in_composition': True,
        'coherent_cooling_factors_are_actual_compiled_source_outputs': True,
        'native_C_dilation_error_diagnostic': errors[0], 'native_E_dilation_error_diagnostic': errors[1],
        'full64_input_Gram_defect_diagnostic': full_error,
        'logical36_input_Gram_defect_diagnostic': float(np.linalg.norm(full_gram[np.ix_(addresses, addresses)]-np.eye(36), 2)),
        'tensor_dilation_error_diagnostic': dilation,
        'complete_native_plus_cooling_Gram_defect_diagnostic': float(np.linalg.norm(complete_gram-np.eye(64), 2)),
        'source_cooling': cooling,
        'fixed_native_composition_dilation_allowance': allowance,
        'fixed_native_plus_cooling_channel_allowance': (2+allowance)*allowance,
        'native_plus_cooling_error_bound_input_domain': 'all36 protected logical pair inputs with an arbitrary untouched quantum reference',
        'full64_native_Gram_check_proves_identity_extension_closeness': False,
        'finite_cooling_dilation_error_propagated_with_arbitrary_reference': True,
        'numerical_allowances_are_outward_or_physical_hardware_bounds': False,
        'native_C_E_and_decoder_words_executed': True,
        'cooling_source_words_executed': True, 'matter_field_history_executed': False}


def execute():
    d = geometry(); arrays = _native_factors()
    cooling, resources = _cooling_source(); arrays.update(cooling)
    arrays.update({'family_projector': np.asarray(d['P'], complex),
        'family_port_frame': np.asarray(d['frame'], complex),
        'family_coefficient_isometry': np.asarray(d['W'], complex),
        'family_cooling_labels': cooling_maps()[1]})
    report = exact_report(); report['actual_native_composition'] = native_metrics(arrays)
    report['cooling_source_resources'] = resources
    return report, arrays


def report(): return execute()[0]


def source_bindings():
    names = set(transport.source_bindings()) | set(objects.source_bindings()) | {
        'model/common_source/native_family_selection.py', 'model/common_source/test_native_family_selection.py',
        'Lean/Screen/A5FamilyBand.lean', 'Lean/Screen/A5PortAction.lean'}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}
