"""A source-executed persistent code with its actual native port action.

A source-derived axis coloring gives a six-dimensional matching code carrying
the actual A5 representation 1+5. Every proper port rotation preserves this
same code. Its sum-zero five-dimensional part is a persistent native object.
Canonical seam reversal complements the coin and reverses signed current;
all raw records, native rows and arbitrary logical references are retained.

The controller chooses an actual port rotation from a committed C record.
A following E repair and a transported quintet marker retain the same object
and its source-defined transformation. Two opposite-color objects support
the source-derived geometric 3+3prime sector: an executed preparation word
and complete local instruments retain its references, complement and aborts.
This does not derive autonomous sector selection or a continuous SM response.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from . import record_code as base, slow_geometry, native_gauge_object as object_geometry
from . import repair_field_attachment as attachment
from .mass_response import history
from ..observer import native
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler


ROOT = Path(__file__).resolve().parents[2]
SEAMS = {'C': (2, 3), 'E': (2, 10)}
SOURCE_ADMISSION = 2e-8
REPLAY_ADMISSION = 2e-9


class _EquivariantCode:
    """Reuse the source compiler with the independently derived matching table."""
    routing = base.routing
    SEAMS = base.SEAMS
    MATCHINGS = object_geometry.equivariant_matchings(0)
    _apply_permutation = staticmethod(base._apply_permutation)
    _permutation_resources = staticmethod(base._permutation_resources)

    @classmethod
    def encoder_permutation(cls):
        return base.routing.complete_permutation(tuple((logical+8*term,
            circuit.encode(tuple(int(i in pair) for i in range(12))))
            for logical, matching in enumerate(cls.MATCHINGS)
            for term, pair in enumerate(matching)))

    @classmethod
    def _encode(cls, value):
        return base._apply_permutation(cls.encoder_permutation(),
            compiler.execute(base._mixing_word(6), 9, value))

    @classmethod
    @lru_cache(maxsize=2)
    def code_matrix(cls):
        return base.routing._immutable(cls._encode(np.eye(512, 6, dtype=complex)))

    @classmethod
    def ideal_code_matrix(cls):
        matrix = np.zeros((512, 6), complex)
        for logical, matching in enumerate(cls.MATCHINGS):
            for pair in matching:
                matrix[circuit.encode(tuple(int(i in pair) for i in range(12))), logical] = 1/math.sqrt(6)
        return matrix

    @classmethod
    def decoding(cls, value):
        moved = base._apply_permutation(cls.encoder_permutation(), value, True)
        return compiler.execute(base._mixing_word(6, True), 9, moved)

    @classmethod
    def branch_images(cls, law, seam, coin, raw):
        base._label(law, seam, coin, raw)
        rows = []
        for matching in cls.MATCHINGS:
            row = []
            for pair in matching:
                state = tuple(int(i in pair) for i in range(12))
                out, emitted = _raw(law, state, seam, coin)
                if emitted == raw:
                    row.append(circuit.encode(out))
            rows.append(tuple(row))
        if len({len(row) for row in rows}) != 1:
            raise ArithmeticError('equivariant code lost recorded branch scalarity')
        flattened = [i for row in rows for i in row]
        if len(flattened) != len(set(flattened)):
            raise ArithmeticError('equivariant code lost orthogonal recorded branch images')
        return tuple(rows)

    @classmethod
    def branch_weight(cls, law, seam, coin, raw):
        return F(len(cls.branch_images(law, seam, coin, raw)[0]), 12)

    @classmethod
    def branch_permutation(cls, law, seam, coin, raw):
        return base.routing.complete_permutation(tuple((logical+8*term, address)
            for logical, row in enumerate(cls.branch_images(law, seam, coin, raw))
            for term, address in enumerate(row)))

    @classmethod
    def apply_recovery(cls, law, seam, coin, raw, value):
        if raw not in circuit.record_values(law, coin):
            raise ValueError('raw aborts cannot select recovery')
        count = len(cls.branch_images(law, seam, coin, raw)[0])
        if not count:
            return value.copy()
        moved = base._apply_permutation(cls.branch_permutation(law, seam, coin, raw), value, True)
        decoded = compiler.execute(base._mixing_word(count, True), 9, moved)
        return cls._encode(decoded)

    @classmethod
    def repair_branches(cls, law, seam, value, *, recover=True):
        branches = base.repair_branches(law, seam, value, recover=False)
        for row in branches:
            coin, raw = row['coin'], row['record_code']
            row['code_branch_weight'] = cls.branch_weight(law, seam, coin, raw)
            if recover and not row['abort']:
                row['amplitude'] = cls.apply_recovery(law, seam, coin, raw, row['amplitude'])
                row['recovery_applied'] = True
                row['recovery_resources'] = cls.recovery_resources(law, seam, coin, raw)
        return branches

    @classmethod
    def encoding_resources(cls):
        return base._counts(base._permutation_resources(cls.encoder_permutation()), base._mixing_resources(6))

    decoding_resources = encoding_resources
    repair_resources = staticmethod(base.repair_resources)

    @classmethod
    def recovery_resources(cls, law, seam, coin, raw):
        count = len(cls.branch_images(law, seam, coin, raw)[0])
        if not count:
            return dict.fromkeys(base.routing.COUNT_KEYS, 0)
        return base._counts(cls.encoding_resources(),
            base._permutation_resources(cls.branch_permutation(law, seam, coin, raw)),
            base._mixing_resources(count, True))


code = _EquivariantCode


class _OppositeEquivariantCode(_EquivariantCode):
    MATCHINGS = object_geometry.equivariant_matchings(1)


PAIR_CODES = {'C': code, 'E': _OppositeEquivariantCode}


def marker_vector():
    value = np.zeros(6)
    value[0], value[1] = 1/math.sqrt(2), -1/math.sqrt(2)
    return value


def logical_rotation(p):
    matrix = np.zeros((6, 6))
    matrix[list(object_geometry.axis_action(p)), np.arange(6)] = 1.
    return matrix


@lru_cache(maxsize=1)
def rotations():
    result = tuple(sorted(tuple(map(int, p)) for p in slow_geometry.source().rotations()))
    if len(result) != 60 or len(set(result)) != 60:
        raise ArithmeticError('the source proper-rotation group differs')
    return result


def _rotation(p):
    if type(p) is not tuple or any(type(i) is not int for i in p) or p not in rotations():
        raise ValueError('a declared proper native port rotation is required')
    return p


def rotate_state(state, p):
    p = _rotation(p)
    result = [0]*12
    for i, load in enumerate(state):
        result[p[i]] = load
    return tuple(result)


@lru_cache(maxsize=60)
def carrier_permutation(p):
    """Full78-state port map, extended by identity on the other434 addresses."""
    p = _rotation(p)
    result = list(range(512))
    for state in native.states():
        result[circuit.encode(state)] = circuit.encode(rotate_state(state, p))
    return code.routing._permutation(tuple(result))


def transported_record(law, seam, coin, raw, p):
    if (law not in ('C', 'E') or seam not in code.SEAMS
            or type(coin) is not int or coin not in (0, 1)
            or type(raw) is not int or not 0 <= raw < 1 << circuit.RECORD_BITS[law]):
        raise ValueError('a complete native seam, coin and raw-record address is required')
    p = _rotation(p)
    a, b = p[seam[0]], p[seam[1]]
    reverse = a > b
    if reverse:
        current = raw & 3
        raw = (raw & ~3)+{0: 0, 1: 2, 2: 1, 3: 3}[current]
        coin = 1-coin
    return tuple(sorted((a, b))), coin, raw


def _raw(law, state, seam, coin):
    out, (_, current) = native.native(state, seam, coin)
    raw = {0: 0, 1: 1, -1: 2}[current]
    if law == 'E':
        raw += 4*((state[seam[0]]+state[seam[1]]) % 2)
    return out, raw


def _matching_sets(p):
    return frozenset(frozenset(tuple(sorted((p[a], p[b]))) for a, b in matching)
                     for matching in code.MATCHINGS)


@lru_cache(maxsize=1)
def selected_rotation():
    # A deterministic geometric controller choice, independent of state or mass.
    marker = marker_vector()
    return next(p for p in rotations() if abs(marker@logical_rotation(p)@marker) < .9
                and p[SEAMS['E'][0]] > p[SEAMS['E'][1]])


def selected_chart(coin, raw):
    if raw not in circuit.record_values('C', coin):
        raise ValueError('an aborted record cannot schedule a later chart')
    current = circuit.record_values('C', coin)[raw]
    return selected_rotation() if current else tuple(range(12))


@lru_cache(maxsize=1)
def exact_covariance():
    """Basis-map equality proves the full coherent recorded instrument law."""
    states = native.states()
    count = 0
    reversal = 0
    for p in rotations():
        for seam in code.SEAMS:
            if p[seam[0]] > p[seam[1]]:
                reversal += 1
            for coin in (0, 1):
                for state in states:
                    for law in ('C', 'E'):
                        out, raw = _raw(law, state, seam, coin)
                        moved_seam, moved_coin, moved_raw = transported_record(law, seam, coin, raw, p)
                        right, actual_raw = _raw(law, rotate_state(state, p), moved_seam, moved_coin)
                        if right != rotate_state(out, p) or actual_raw != moved_raw:
                            raise ArithmeticError('full native recorded port covariance failed')
                        count += 1
            for law in ('C', 'E'):
                seen = set()
                for coin in (0, 1):
                    for raw in range(1 << circuit.RECORD_BITS[law]):
                        _, c, r = transported_record(law, seam, coin, raw, p)
                        seen.add((c, r))
                        if ((raw in circuit.record_values(law, coin)) != (r in circuit.record_values(law, c))):
                            raise ArithmeticError('port transport changes raw-abort status')
                if len(seen) != 2*(1 << circuit.RECORD_BITS[law]):
                    raise ArithmeticError('complete raw-record transport is not a bijection')
    families = {_matching_sets(p) for p in rotations()}
    stabilizer = [p for p in rotations() if _matching_sets(p) == _matching_sets(tuple(range(12)))]
    if len(families) != 1 or len(stabilizer) != 60:
        raise ArithmeticError('matching-code fiber orbit or stabilizer differs')
    return {'proper_rotations': 60, 'native_integer_carrier_states': 78,
        'native_seams': 30, 'basis_record_equalities_checked': count,
        'rotated_ordered_seams_requiring_orientation_reversal': reversal,
        'raw_record_bijections_include_all_abort_addresses': True,
        'canonical_seam_reversal': 'coin -> 1-coin; raw current code1 <-> code2; parity bit unchanged',
        'complete_coherent_instrument_identity': 'K[g.seam,g.record] U_g = U_g K[seam,record] on all78 native states',
        'same_fixed_code_preserved_by_all_rotations': True,
        'fixed_code_stabilizer_order': len(stabilizer), 'distinct_code_projectors': len(families),
        'general_port_intervention_values': 60, 'general_port_intervention_bits_upper': 6,
        'coherent_chart_register_supplied': False,
        'actual_code_rotation_representation': '1+5; the protected sum-zero part is the geometric quintet',
        'continuous_M6_or_SM_action_follows_from_port_covariance': False}


def transport(p, value, *, inverse=False):
    return code._apply_permutation(carrier_permutation(p), value, inverse)


def chart_recovery(law, seam, coin, raw, p, value):
    if raw not in circuit.record_values(law, coin):
        raise ValueError('raw aborts halt instead of selecting chart recovery')
    physical_seam, physical_coin, physical_raw = transported_record(law, seam, coin, raw, p)
    return code.apply_recovery(law, physical_seam, physical_coin, physical_raw, value)


def marker_read(p, value):
    """A transported rank-one quintet marker preserves the five-state sector."""
    decoded = code.decoding(value)
    vector = np.zeros(8)
    vector[:6] = logical_rotation(p)@marker_vector()
    reading = history.projector_read(decoded, vector, targets=(0, 1, 2))
    return tuple(code._encode(v) for v in reading.outputs), reading.resources


def record_menu():
    rows = []
    for cc in (0, 1):
        for cr in range(4):
            if cr not in circuit.record_values('C', cc):
                rows.append((1, cc, cr, -1, -1, -1, -1))
                continue
            p = selected_chart(cc, cr)
            frame = rotations().index(p)
            for ec in (0, 1):
                for er in range(8):
                    if er not in circuit.record_values('E', ec):
                        rows.append((2, cc, cr, frame, ec, er, -1))
                    else:
                        for bit in (0, 1):
                            rows.append((0, cc, cr, frame, ec, er, bit))
    return np.asarray(rows, np.int64)


def _original_record(p, ec, er):
    inverse = tuple(p.index(i) for i in range(12))
    seam = tuple(sorted((p[SEAMS['E'][0]], p[SEAMS['E'][1]])))
    original_seam, coin, raw = transported_record('E', seam, ec, er, inverse)
    if original_seam != SEAMS['E']:
        raise ArithmeticError('transported native seam did not return to the original chart')
    return coin, raw


def execute_source():
    initial = code.code_matrix()
    crows = code.repair_branches('C', SEAMS['C'], initial)
    outputs = {}
    prefix_frames = []
    pointer = None
    for p in (tuple(range(12)), selected_rotation()):
        relevant = [row for row in crows if not row['abort'] and selected_chart(row['coin'], row['record_code']) == p]
        if not relevant:
            continue
        incoming = np.concatenate([row['amplitude'] for row in relevant], axis=1)
        changed = transport(p, incoming)
        prefix_frames.append(changed[:, :6])
        seam = tuple(sorted((p[SEAMS['E'][0]], p[SEAMS['E'][1]])))
        erows = code.repair_branches('E', seam, changed, recover=False)
        recovered = []
        for erow in erows:
            ec, er = erow['coin'], erow['record_code']
            if erow['abort']:
                for i, crow in enumerate(relevant):
                    label = (2, crow['coin'], crow['record_code'], rotations().index(p), ec, er, -1)
                    outputs[label] = erow['amplitude'][:, 6*i:6*(i+1)]
                continue
            old_coin, old_raw = _original_record(p, ec, er)
            value = chart_recovery('E', SEAMS['E'], old_coin, old_raw, p, erow['amplitude'])
            recovered.append((ec, er, value))
        combined = np.concatenate([value for _, _, value in recovered], axis=1)
        marked, pointer = marker_read(p, combined)
        width = 6*len(relevant)
        for k, (ec, er, _) in enumerate(recovered):
            for i, crow in enumerate(relevant):
                for bit in (0, 1):
                    label = (0, crow['coin'], crow['record_code'], rotations().index(p), ec, er, bit)
                    outputs[label] = marked[bit][:, k*width+6*i:k*width+6*(i+1)]
    for row in crows:
        if row['abort']:
            outputs[(1, row['coin'], row['record_code'], -1, -1, -1, -1)] = row['amplitude']
    menu = record_menu()
    if set(outputs) != set(map(tuple, menu.tolist())):
        raise ArithmeticError('source history omitted a terminal native record')
    arrays = {'outputs': np.asarray([outputs[tuple(label)] for label in menu.tolist()]),
        'labels': menu, 'transported_prefixes': np.asarray(prefix_frames)}
    return arrays, pointer


def product_addresses():
    """Six physical code colors within each complete three-bit register."""
    return np.asarray([8*a+b for a in range(6) for b in range(6)], np.int64)


@lru_cache(maxsize=1)
def geometric_frame():
    import sympy as symbolic
    frames = object_geometry.geometric_intertwiners()
    exact = symbolic.Matrix.hstack(*[frames[name]['frame']/symbolic.sqrt(frames[name]['metric'])
        for name in ('three', 'three_prime')])
    return base.routing._immutable(np.asarray(exact.evalf(20), float).astype(complex))


@lru_cache(maxsize=1)
def preparation_plan():
    """Householder extension derived only from the actual port intertwiner.

    Each reflection maps the next j tensor blank input to its target column
    while fixing all preceding columns. The resulting complete 64-dimensional
    orthogonal matrix is compiled as U diag(-1,1,...) U inverse; U itself is
    a source state-preparation word. No target matrix is applied to the state.
    """
    target = np.zeros((64, 6))
    target[product_addresses()] = geometric_frame().real
    extension = np.eye(64)
    vectors = []
    for j in range(6):
        difference = extension[:, 8*j]-target[:, j]
        norm = np.linalg.norm(difference)
        if norm < 1e-12:
            raise ArithmeticError('the fixed geometric preparation reflection schedule changed')
        vector = difference/norm
        extension -= 2*np.outer(vector, vector@extension)
        vectors.append(base.routing._immutable(vector))
    if (np.linalg.norm(extension.T@extension-np.eye(64), 2) > 1e-12
            or np.linalg.norm(extension[:, 8*np.arange(6)]-target) > 1e-12):
        raise ArithmeticError('geometric source-preparation orthogonal extension failed')
    return tuple(vectors), base.routing._immutable(extension.astype(complex))


def _preparation_stages():
    phases = np.zeros(64)
    phases[0] = math.pi
    reflection = compiler.compile_diagonal(phases)
    for vector in preparation_plan()[0]:
        # The existing fixed binary-tree compiler works at six bits as well;
        # its public state helper preserves older Gray words at this size.
        # Use the tree word directly without changing those frozen callers.
        word = compiler._compile_real_state_tree(vector.real, 6)
        if len(word) != 963:
            raise ArithmeticError('the complete six-bit amplitude-tree schedule changed')
        yield compiler.inverse(word)+reflection+word


def prepare_geometric_source():
    state = np.eye(64, dtype=complex)[:, 8*np.arange(6)]
    for word in _preparation_stages():
        state = compiler.execute(word, 6, state)
    return state


@lru_cache(maxsize=1)
def _preparation_resources():
    rows = [compiler.resources(word) for word in _preparation_stages()]
    counts = {key: sum(row[key] for row in rows) for key in base.routing.COUNT_KEYS}
    return tuple(counts[key] for key in base.routing.COUNT_KEYS)


def preparation_resources():
    return dict(zip(base.routing.COUNT_KEYS, _preparation_resources(), strict=True))


def pair_source():
    """Actual preparation and two complete native instruments, kept factored.

    The 64-row preparation word acts on all six unknown j tensor blank inputs.
    Each local source encoder and subsequent instrument retains all eight
    decoded-color columns, including padded 6/7. Their tensor factors therefore
    propagate every prepared row, with all 128 joint raw outcomes retained.
    Both local attempts are scheduled regardless of the other local record.
    """
    p = selected_rotation()
    result = {'prepared_geometric_inputs': prepare_geometric_source(),
        'preparation_orthogonal_extension': preparation_plan()[1].copy()}
    for law in ('C', 'E'):
        local_code = PAIR_CODES[law]
        rotated = transport(p, local_code._encode(np.eye(512, 8, dtype=complex)))
        physical_seam = tuple(sorted((p[SEAMS[law][0]], p[SEAMS[law][1]])))
        branches = local_code.repair_branches(law, physical_seam, rotated)
        result['pair_'+law] = np.asarray([row['amplitude'] for row in branches])
        result['pair_'+law+'_labels'] = np.asarray([(row['coin'], row['record_code']) for row in branches], np.int64)
    import sympy as symbolic
    result['geometric_frame'] = geometric_frame().copy()
    result['product_complement'] = np.asarray((symbolic.eye(36)-object_geometry.geometric_projectors()['six']).evalf(20), float).astype(complex)
    return result


def pair_metrics(arrays):
    p = selected_rotation()
    logical = logical_rotation(p)
    grams, errors, failures = [], [], []
    for law in ('C', 'E'):
        local_code = PAIR_CODES[law]
        target_code = local_code.ideal_code_matrix()@logical
        values, labels = arrays['pair_'+law], arrays['pair_'+law+'_labels']
        expected_labels = np.array([(c, r) for c in (0, 1) for r in range(1 << circuit.RECORD_BITS[law])], np.int64)
        if (type(values) is not np.ndarray or values.dtype != np.dtype(complex)
                or values.shape != (len(expected_labels), 512, 8) or not np.isfinite(values).all()
                or type(labels) is not np.ndarray or labels.dtype != np.dtype('int64')
                or not np.array_equal(labels, expected_labels)):
            raise ValueError('complete local pair instrument and raw labels required')
        physical_seam = tuple(sorted((p[SEAMS[law][0]], p[SEAMS[law][1]])))
        target = np.zeros((len(values), 512, 6), complex)
        failure = np.zeros((8, 8), complex)
        for index, (coin, raw) in enumerate(labels.tolist()):
            if raw in circuit.record_values(law, coin):
                target[index] = math.sqrt(float(local_code.branch_weight(law, physical_seam, coin, raw)))*target_code
            else:
                failure += values[index].conj().T@values[index]
        errors.append(float(np.linalg.norm(values[:, :, :6]-target)))
        grams.append(np.einsum('nac,nad->cd', values.conj(), values, optimize=True))
        failures.append(failure)
    frame, complement = arrays['geometric_frame'], arrays['product_complement']
    prepared, extension = arrays['prepared_geometric_inputs'], arrays['preparation_orthogonal_extension']
    for value, shape in ((frame, (36, 6)), (complement, (36, 36)), (prepared, (64, 6)), (extension, (64, 64))):
        if type(value) is not np.ndarray or value.dtype != np.dtype(complex) or value.shape != shape or not np.isfinite(value).all():
            raise ValueError('original finite geometric preparation, input and complement arrays required')
    product_action = np.kron(logical, logical)
    projector = frame@frame.conj().T
    if (np.linalg.norm(frame-geometric_frame()) > 1e-12
            or np.linalg.norm(frame.conj().T@frame-np.eye(6), 2) > 1e-11
            or np.linalg.norm(projector+complement-np.eye(36), 2) > 1e-11
            or np.linalg.norm(complement-complement.conj().T, 2) > 1e-11
            or np.linalg.norm(complement@complement-complement, 2) > 1e-11
            or np.linalg.norm(complement@product_action@frame, 2) > 1e-11):
        raise ValueError('source-derived geometric domain/complement/action identity failed')
    addresses = product_addresses()
    target_preparation = np.zeros((64, 6), complex)
    target_preparation[addresses] = frame
    if (np.linalg.norm(extension-preparation_plan()[1]) > 1e-12
            or np.linalg.norm(extension.conj().T@extension-np.eye(64), 2) > 1e-11
            or np.linalg.norm(extension[:, 8*np.arange(6)]-target_preparation) > 1e-11):
        raise ValueError('the complete source-preparation orthogonal extension differs')
    preparation_error = float(np.linalg.norm(prepared-target_preparation))
    complete = np.kron(*grams)
    restricted = prepared.conj().T@complete@prepared
    full_defect = float(np.linalg.norm(complete-np.eye(64), 2))
    admitted_defect = float(np.linalg.norm(complete[np.ix_(addresses, addresses)]-np.eye(36), 2))
    restricted_defect = float(np.linalg.norm(restricted-np.eye(6), 2))
    source_bound = errors[0]*(1+errors[1])+errors[1]
    combined_bound = source_bound+math.sqrt(1+full_defect)*preparation_error
    if (max(errors+[preparation_error]) > SOURCE_ADMISSION
            or max(full_defect, restricted_defect) > 4*SOURCE_ADMISSION):
        raise ValueError('factored native pair exceeds complete input/reference source admission')
    local_allowance = SOURCE_ADMISSION
    allowance = local_allowance*(1+local_allowance)+local_allowance+math.sqrt(1+4*local_allowance)*local_allowance
    total_halt = np.kron(failures[0], grams[1])+np.kron(grams[0]-failures[0], failures[1])
    padding = np.setdiff1d(np.arange(64), addresses)
    return {'source_matching_choices': [0, 1],
        'pair_persistence_schedule_includes_candidate_interaction_dynamics': False,
        'geometric_preparation_uses_declared_entangling_source_controls': True,
        'geometric_preparation_is_free_or_autonomous': False,
        'local_full_native_rows': 512, 'local_complete_decoded_input_columns': 8,
        'full_product_register_input_columns': 64, 'full_product_logical_input_columns': 36,
        'retained_padded_product_register_levels': 28,
        'selected_source_geometric_input_columns': 6, 'retained_complementary_product_levels': 30,
        'complete_joint_raw_record_outcomes': 128,
        'joint_output_definition': '(K_C[r] tensor K_E[s]) A with actual 64-row source preparation A; factors retain all 64 register inputs',
        'geometric_frame_is_source_derived_not_a_copied_M6_current_table': True,
        'source_local_C_dilation_error_diagnostic': errors[0],
        'source_local_E_dilation_error_diagnostic': errors[1],
        'source_preparation_dilation_error_diagnostic': preparation_error,
        'source_preparation_padded_row_norm_diagnostic': float(np.linalg.norm(prepared[padding])),
        'complete_register_Gram_defect_diagnostic': full_defect,
        'complete_product_Gram_defect_diagnostic': admitted_defect,
        'geometric_input_Gram_defect_diagnostic': restricted_defect,
        'tensor_dilation_error_upper_diagnostic': source_bound,
        'preparation_and_tensor_dilation_error_upper_diagnostic': combined_bound,
        'full_register_raw_halt_Frobenius_probability_diagnostic': float(np.trace(total_halt).real),
        'prepared_geometric_raw_halt_Frobenius_probability_diagnostic': float(np.trace(prepared.conj().T@total_halt@prepared).real),
        'fixed_product_dilation_allowance': allowance,
        'fixed_product_channel_error_allowance': (2+allowance)*allowance,
        'actual_simultaneous_port_action_preserves_geometric6': True,
        'independent_local_rotations_preserve_geometric6': False,
        'geometric_sector_state_preparation_source_word_executed': True,
        'preparation_unknown_input_domain': 'j tensor blank 0 for j=0,...,5 on two three-bit registers',
        'preparation_uses_measured_particle_parameters': False,
        'complete_orthogonal_extension_columns_source_executed': 6,
        'native_repair_selects_the_geometric_six_subspace': False,
        'spin_half_or_fermionic_statistics_derived': False,
        'tensor_factor_antisymmetry_identified_with_physical_exchange_statistics': False,
        'native_product_encoders_and_local_repairs_executed': True,
        'other_product_sectors_discarded_or_declared_unphysical': False,
        'physical_SM_particle_or_continuous_gauge_action_identified': False}


def _ideal_rotated_code(p):
    result = np.zeros((512, 6), complex)
    original = code.ideal_code_matrix()
    result[list(carrier_permutation(p))] = original
    return result


def metrics(arrays):
    menu = record_menu()
    shapes = {'outputs': (len(menu), 512, 6), 'labels': menu.shape, 'transported_prefixes': (2, 512, 6),
        'pair_C': (8, 512, 8), 'pair_C_labels': (8, 2),
        'pair_E': (16, 512, 8), 'pair_E_labels': (16, 2),
        'geometric_frame': (36, 6), 'product_complement': (36, 36),
        'prepared_geometric_inputs': (64, 6), 'preparation_orthogonal_extension': (64, 64)}
    if type(arrays) is not dict or set(arrays) != set(shapes):
        raise ValueError('complete transported-object array menu differs')
    for key, shape in shapes.items():
        a = arrays[key]
        dtype = np.dtype('int64') if key == 'labels' or key.endswith('_labels') else np.dtype(complex)
        if type(a) is not np.ndarray or a.shape != shape or a.dtype != dtype or not np.isfinite(a).all():
            raise ValueError('original full native output array schema differs: '+key)
    if not np.array_equal(arrays['labels'], menu):
        raise ValueError('complete record/abort/chart menu differs')
    maps = arrays['outputs']
    target = np.zeros_like(maps)
    projector = np.outer(marker_vector(), marker_vector())
    marker = (projector, np.eye(6)-projector)
    for i, label in enumerate(menu.tolist()):
        stop, cc, cr, frame, ec, er, bit = label
        if stop:
            continue
        p = rotations()[frame]
        old_coin, old_raw = _original_record(p, ec, er)
        weight = code.branch_weight('C', SEAMS['C'], cc, cr)*code.branch_weight('E', SEAMS['E'], old_coin, old_raw)
        target[i] = math.sqrt(float(weight))*_ideal_rotated_code(p)@marker[bit]
    error = float(np.linalg.norm(maps-target))
    gram = np.einsum('nac,nad->cd', maps.conj(), maps, optimize=True)
    defect = float(np.linalg.norm(gram-np.eye(6), 2))
    if error > SOURCE_ADMISSION or defect > SOURCE_ADMISSION:
        raise ValueError('transported native object exceeds complete input/reference admission')
    abort = float(np.linalg.norm(maps[menu[:, 0] != 0])**2)
    wrong_probability = float(abs(marker_vector()@logical_rotation(selected_rotation())@marker_vector())**2)
    if wrong_probability >= 1-1e-12:
        raise ArithmeticError('wrong-chart marker is not a discriminating control')
    return {'complete_unknown_input_columns': 6, 'terminal_records': len(menu),
        'retained_native_output_rows': 512,
        'complete_history_dilation_error_diagnostic': error,
        'complete_input_Gram_defect_diagnostic': defect,
        'raw_halt_Frobenius_probability_diagnostic': abort,
        'logical_marker_correct_chart_probability': 1.,
        'logical_marker_wrong_fixed_chart_probability': wrong_probability,
        'five_state_sum_zero_sector_preserved_before_and_after_marker': True,
        'chart_choice_uses_only_committed_native_current': True,
        'target_history_preserves_all_logical_coherences_until_the_declared_marker_read': True,
        'source_history_equals_ideal_code_projection_or_reset': False,
        'fixed_history_channel_error_allowance': (2+SOURCE_ADMISSION)*SOURCE_ADMISSION,
        'numerical_admission_is_outward_certificate_or_physical_service_error': False}


def resource_counts(pointer):
    costs = {}
    for p in (tuple(range(12)), selected_rotation()):
        physical_e = tuple(sorted((p[SEAMS['E'][0]], p[SEAMS['E'][1]])))
        rotation = code._permutation_resources(carrier_permutation(p))['native_pulses']
        c = code.repair_resources('C', SEAMS['C'])['native_pulses']+max(
            code.recovery_resources('C', SEAMS['C'], coin, raw)['native_pulses']
            for coin in (0, 1) for raw in circuit.record_values('C', coin))
        e = code.repair_resources('E', physical_e)['native_pulses']+max(
            code.recovery_resources('E', physical_e, coin, raw)['native_pulses']
            for coin in (0, 1) for raw in circuit.record_values('E', coin))
        encoding = 2*code.encoding_resources()['native_pulses']+code.decoding_resources()['native_pulses']
        total = c+e+encoding+rotation+pointer['native_pulses']
        costs[str(rotations().index(p))] = {'C_repair_recovery_upper': c, 'E_repair_recovery_upper': e,
            'actual_port_rotation_invocation': rotation,
            'encoding_decoding': encoding, 'marker_read': pointer['native_pulses'],
            'native_pulses_upper': total}
    maximum = max(value['native_pulses_upper'] for value in costs.values())
    p = selected_rotation()
    pair = sum(local_code.encoding_resources()['native_pulses'] for local_code in PAIR_CODES.values())
    pair += 2*code._permutation_resources(carrier_permutation(p))['native_pulses']
    for law in ('C', 'E'):
        local_code = PAIR_CODES[law]
        seam = tuple(sorted((p[SEAMS[law][0]], p[SEAMS[law][1]])))
        pair += local_code.repair_resources(law, seam)['native_pulses']+max(
            local_code.recovery_resources(law, seam, coin, raw)['native_pulses']
            for coin in (0, 1) for raw in circuit.record_values(law, coin))
    preparation = preparation_resources()
    pair += preparation['native_pulses']
    return {'per_reached_frame': costs, 'one_history_native_pulses_upper': maximum,
        'required_proper_load_export_service_events_upper': 2*maximum,
        'two_carrier_geometric_preparation': preparation,
        'two_carrier_preparation_persistent_history_native_pulses_upper': pair,
        'two_carrier_preparation_persistent_history_proper_services_upper': 2*pair,
        'two_carrier_decoded_register_qubits': 6, 'two_carrier_native_decoding_invocations': 0,
        'two_carrier_input_decoded_registers_and_blank_state_supplied': True, 'two_carrier_additional_blank_encoding_qubits': 12,
        'geometric_source_preparation_reflections': 6,
        'geometric_source_preparation_interaction_and_phase_services_counted': True,
        'extra_chart_memory_bits_for_record_determined_rotation': 0, 'native_coin_raw_record_and_marker_bits': 8,
        'coherent_chart_register_qubits': 0,
        'source_frame_is_an_admitted_classical_control_not_a_hidden_native_state_read': True,
        'physical_load_export_protection_and_memory_services_executed': False,
        'additional_physical_error_contract': 'add required_service_events*eta for per-service diamond error eta, plus marker/preparation/memory errors'}


def execute():
    exact = exact_covariance()
    arrays, pointer = execute_source()
    arrays.update(pair_source())
    return {'schema': 'oph.native_code_transport.v2',
        'classification': 'source-prepared persistent native code, quintet and geometric pair6 under actual port rotations',
        'compatibility': 'the source-derived axis matching code and its invariant sum-zero sector; all native C/E records have scalar compressed Gram',
        'exact_native_covariance': exact,
        'protocol': {'native_seams_in_original_chart': {k: list(v) for k, v in SEAMS.items()},
            'source_equivariant_coloring_choice': 0,
            'nonidentity_rotation': list(selected_rotation()),
            'chart_selection': 'nonzero committed C current selects the first lexical proper rotation moving the declared marker and reversing the E seam; zero selects identity',
            'clock_parameter': 'ordered source intervention events; no physical-duration interpretation',
            'chart_recovery': 'the same source-derived code recovery on the actual physical seam and raw record',
            'marker': '(color0-color1)/sqrt2, transported by the actual axis permutation; both pointer outcomes retained'},
        'admissions': {'complete_source_dilation': SOURCE_ADMISSION, 'fresh_binary_replay': REPLAY_ADMISSION},
        'execution': metrics(arrays), 'two_carrier_geometric_object': pair_metrics(arrays), 'resources': resource_counts(pointer),
        'scope': {'persistent_quantum_memory_component': True,
            'actual_native_rotation_and_record_instrument_share_one_carrier': True,
            'same_fixed_code_invariant_under_all_proper_rotations': True,
            'continuous_SM_gauge_action_or_particle_identification_derived': False,
            'arbitrary_code_isometry_used_to_identify_M6_generators': False}}, arrays


def source_bindings():
    paths = set(attachment.source_bindings()) | set(object_geometry.source_bindings()) | {
        'model/common_source/native_code_transport.py', 'model/common_source/test_native_code_transport.py',
        'model/common_source/native_gauge_object.py',
        'model/common_source/mass_response/history.py', 'tools/strict_json.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def _compare(left, right, path='receipt'):
    if type(left) is not type(right):
        raise ValueError(path+' original type differs')
    if type(left) is dict:
        if left.keys() != right.keys():
            raise ValueError(path+' complete keys differ')
        for key in left:
            _compare(left[key], right[key], path+'.'+key)
    elif type(left) is list:
        if len(left) != len(right):
            raise ValueError(path+' length differs')
        for i, (a, b) in enumerate(zip(left, right, strict=True)):
            _compare(a, b, path+f'[{i}]')
    elif type(left) is float:
        exact = path.startswith('receipt.admissions.') or path.endswith(('fixed_history_channel_error_allowance',
            'fixed_product_dilation_allowance', 'fixed_product_channel_error_allowance'))
        if not math.isfinite(right) or (left != right if exact else not math.isclose(left, right, rel_tol=2e-9, abs_tol=1e-10)):
            raise ValueError(path+' numeric value differs')
    elif left != right:
        raise ValueError(path+' exact value differs')


def compare_execution(report, arrays, fresh, fresh_arrays):
    if type(arrays) is not dict or type(fresh_arrays) is not dict or set(arrays) != set(fresh_arrays):
        raise ValueError('retained and fresh complete native array menus differ')
    metrics(fresh_arrays)
    own = metrics(arrays)
    pair = pair_metrics(arrays)
    for key, value in own.items():
        if key.endswith('_diagnostic'):
            old = report['execution'].get(key)
            if type(old) is not float or not math.isclose(value, old, rel_tol=2e-10, abs_tol=64*np.finfo(float).eps):
                raise ValueError('diagnostic differs from its own binary source arrays: '+key)
    for key, value in pair.items():
        if key.endswith('_diagnostic'):
            old = report['two_carrier_geometric_object'].get(key)
            if type(old) is not float or not math.isclose(value, old, rel_tol=2e-10, abs_tol=64*np.finfo(float).eps):
                raise ValueError('pair diagnostic differs from its own binary source arrays: '+key)
    for key in arrays:
        if arrays[key].dtype.kind in 'iu':
            if not np.array_equal(arrays[key], fresh_arrays[key]):
                raise ValueError('retained complete labels differ')
        elif np.linalg.norm(arrays[key]-fresh_arrays[key]) > REPLAY_ADMISSION:
            raise ValueError('retained complete native outputs differ: '+key)
    normalized = deepcopy(report)
    for key in own:
        if key.endswith('_diagnostic'):
            normalized['execution'][key] = fresh['execution'][key]
    for key in pair:
        if key.endswith('_diagnostic'):
            normalized['two_carrier_geometric_object'][key] = fresh['two_carrier_geometric_object'][key]
    _compare(fresh, normalized)


def check(path):
    path = Path(path)
    report = strict_json.load(path)
    pins = source_bindings()
    if type(report) is not dict or report.get('source_bindings') != pins:
        raise ValueError('transported code source bindings differ')
    artifact = path.with_suffix('.npz')
    if report.get('artifact') != {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError('transported code binary custody differs')
    with np.load(artifact, allow_pickle=False) as saved:
        if len(saved.files) != len(set(saved.files)):
            raise ValueError('duplicate native source artifact keys')
        arrays = {key: saved[key] for key in saved.files}
    fresh, outputs = execute()
    compare_execution({k: v for k, v in report.items() if k not in ('source_bindings', 'artifact')}, arrays, fresh, outputs)
    if source_bindings() != pins:
        raise ValueError('native code source changed during replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--output', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print('Native code-object chart transport and retained records: PASS')
        return
    artifact = args.output.with_suffix('.npz')
    if args.output.exists() or artifact.exists():
        parser.error('fresh receipt and artifact paths required')
    pins = source_bindings()
    report, arrays = execute()
    if source_bindings() != pins:
        raise ValueError('native code source changed during build')
    report['source_bindings'] = pins
    with artifact.open('xb') as out:
        np.savez_compressed(out, **arrays)
    report['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    with args.output.open('x') as out:
        json.dump(report, out, indent=2, sort_keys=True, allow_nan=False)
        out.write('\n')


if __name__ == '__main__':
    main()
