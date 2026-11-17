"""Unconditional reference cooling with an explicit fresh bath and native clock.

A fresh six-bit bath is prepared by an emitted source word in the ground
state of native_response_biased. A complete native C attempt on a separate
fresh counter selects a SWAP for coin0/raw1. Every other label leaves the
reference alone. The old quantum information remains in the bath, including
its correlations with arbitrary untouched systems. No successful branch is
selected or normalized. Fresh counter/bath preparation and this feedback
coupling are added apparatus choices; passive C/E repair does not cool.

The reference here is the decoded64-row core of two geometric objects.
Its native four-carrier geometric embedding is the upstream interface.
The retained factors specify every counter, reference and bath output of
this new preparation stage, without materializing their tensor product.
"""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

import numpy as np
import sympy as s

from . import native_code_transport as transport
from . import native_response_reference as reference
from . import record_code, native_routing
from ..observer.source_scalar import compiler
from ..observer.source_discrimination import circuit

ROOT = Path(__file__).resolve().parents[2]
SOURCE_ADMISSION = 2e-7
REPLAY_ADMISSION = 4e-9
TRIGGER = (0, 1)
COUNTER_SEAM = transport.SEAMS['C']
BIAS = F(16)


def ground_vector():
    """Owned source-derived ground in the existing8-by8 geometric chart."""
    from . import native_response_biased as biased
    value = np.array(biased.ground64(BIAS), dtype=float, copy=True)
    if value.shape != (64,) or not np.isfinite(value).all() or abs(np.linalg.norm(value)-1) > 2e-12:
        raise ArithmeticError('the source-derived biased ground is not a normalized64-row state')
    return value


def raw_records():
    code = transport.PAIR_CODES['C']
    rows = tuple((coin, raw, code.branch_weight('C', COUNTER_SEAM, coin, raw))
                 for coin in (0, 1) for raw in range(4))
    if (sum(row[2] for row in rows) != 1
            or next(w for coin, raw, w in rows if (coin, raw) == TRIGGER) != F(1, 12)):
        raise ArithmeticError('the complete native clock or cooling-trigger probability changed')
    return rows


def remaining_weight(cycles):
    if type(cycles) is not int or not 0 <= cycles <= 4096:
        raise ValueError('an original cycle count from0 through4096 is required')
    return F(11, 12)**cycles


def clock_data():
    return {'law': 'C', 'seam': COUNTER_SEAM, 'trigger': TRIGGER, 'records': raw_records()}


def trigger_record():
    return TRIGGER


def swap_word():
    """Literal SWAP of low reference bits0..5 and fresh bath bits6..11."""
    return [event for bit in range(6) for control, target in
            ((bit, bit+6), (bit+6, bit), (bit, bit+6))
            for event in compiler.cnot(control, target)]


def _counter():
    code = transport.PAIR_CODES['C']
    incoming = code._encode(np.eye(512, 1, dtype=complex))
    rows = code.repair_branches('C', COUNTER_SEAM, incoming)
    resources = {'encoding': code.encoding_resources(),
        'repair': code.repair_resources('C', COUNTER_SEAM),
        'maximum_recovery_native_pulses': max(row['recovery_resources']['native_pulses']
            for row in rows if not row['abort'])}
    return (np.asarray([row['amplitude'][:, 0] for row in rows]),
            np.asarray([(row['coin'], row['record_code']) for row in rows], np.int64), resources)


def _schema(arrays):
    shapes = {'bath_preparation': (64, 64), 'unchanged': (4096, 64), 'replaced': (4096, 64),
              'counter': (8, 512), 'labels': (8, 2)}
    if type(arrays) is not dict or arrays.keys() != shapes.keys():
        raise ValueError('the complete reference-cooling array menu is required')
    for name, shape in shapes.items():
        value = arrays[name]
        dtype = np.dtype('int64') if name == 'labels' else np.dtype(complex)
        if type(value) is not np.ndarray or value.shape != shape or value.dtype != dtype or not np.isfinite(value).all():
            raise ValueError('original finite reference-cooling array schema differs: '+name)


def output_block(arrays, coin, raw, counter_row, input_columns=None):
    """One coherent counter-output slice, keeping all4096 reference/bath rows.

    Counter rows are coherent output indices, never selected measurement
    outcomes. Collecting all512 slices and all8 raw labels reconstructs the
    complete instrument. Input columns may carry an arbitrary reference.
    """
    _schema(arrays)
    if (type(coin) is not int or coin not in (0, 1) or type(raw) is not int or not 0 <= raw < 4
            or type(counter_row) is not int or not 0 <= counter_row < 512):
        raise ValueError('original in-range native raw labels and counter row required')
    value = np.eye(64, dtype=complex) if input_columns is None else input_columns
    if (type(value) is not np.ndarray or value.ndim not in (1, 2) or value.shape[0] != 64
            or value.dtype.kind not in 'fc' or not np.isfinite(value).all()):
        raise ValueError('finite floating64-row input or reference columns required')
    factor = arrays['replaced'] if (coin, raw) == TRIGGER else arrays['unchanged']
    return arrays['counter'][4*coin+raw, counter_row]*(factor@value)


def ideal_marginal(operator, cycles, *, reference_dimension=1):
    """Ideal unconditional channel on the reference core and untouched ancilla.

    This linear map accepts arbitrary finite operators. It does not repair
    positivity, normalize inputs, or substitute an ideal map in execution.
    """
    weight = float(remaining_weight(cycles))
    if type(reference_dimension) is not int or not 1 <= reference_dimension <= 16:
        raise ValueError('an original untouched-reference dimension from1 through16 required')
    r = reference_dimension
    if (type(operator) is not np.ndarray or operator.shape != (64*r, 64*r)
            or operator.dtype.kind not in 'fc' or not np.isfinite(operator).all()):
        raise ValueError('a finite original operator on the stated reference tensor ancilla is required')
    old_ancilla = np.einsum('iaib->ab', operator.reshape(64, r, 64, r))
    vector = ground_vector()
    return weight*operator+(1-weight)*np.kron(np.outer(vector, vector), old_ancilla)


def exact_certificate():
    records = raw_records()
    # Independent two-qubit Zeno counterexample: monitoring the second
    # system cannot remove an X current acting on the first system.
    x, p = s.Matrix([[0, 1], [1, 0]]), s.diag(1, 0)
    P = s.kronecker_product(s.eye(2), p)
    collective = s.kronecker_product(x, s.eye(2))+s.kronecker_product(s.eye(2), x)
    if P*collective*P != s.kronecker_product(x, p):
        raise ArithmeticError('reference-only monitoring counterexample failed')
    return {'schema': 'oph.native_response_cooling.v1',
        'native_clock': {'law': 'C', 'seam': list(COUNTER_SEAM), 'trigger': list(TRIGGER),
            'trigger_probability': '1/12',
            'complete_raw_records': [{'coin': coin, 'raw': raw, 'weight': str(weight),
                'swap': (coin, raw) == TRIGGER, 'abort': raw not in circuit.record_values('C', coin)}
                for coin, raw, weight in records],
            'fresh_separate_counter_each_cycle': True},
        'unconditional_channel': {'reference_dimension': 64,
            'one_cycle': 'C(rho)=11rho/12+|psi><psi| Tr(rho)/12',
            'n_cycles': 'C^n(rho)=q rho+(1-q)|psi><psi| Tr(rho), q=(11/12)^n',
            'untouched_reference_extension': 'q rho_WR+(1-q)|psi><psi| tensor Tr_W(rho_WR)',
            'half_trace_distance_to_reset_upper': '(11/12)^n for normalized states, including arbitrary untouched references',
            'diamond_distance_to_reset': '2(11/12)^n',
            'no_success_postselection_or_normalization': True,
            'old_quantum_information_is_retained_in_bath': True},
        'shared_clock_correlation': {'no_trigger': 'the original reference remains',
            'one_or_more_triggers': 'the reference equals the prepared ground exactly in the ideal law',
            'first_trigger_moves_original_reference_information_to_its_bath': True,
            'later_ground_to_ground_swaps_discard_no_old_information': True,
            'independent_cooling_and_control_event_samples_required': False},
        'fixed_repair_no_go': {'corrected_native_branch_on_code': 'K_h=sqrt(w_h) I',
            'record_only_unitary_feedback': 'sum_h w_h U_h rho U_h^dagger is unital',
            'maximally_mixed_state_is_fixed': True,
            'pure_global_attractor_without_additional_entropy_sink': False,
            'scope': 'the fixed corrected code and state-independent records; not arbitrary added baths or data-dependent measurements'},
        'monitoring_counterexample': {'projection': 'P=I_data tensor |psi><psi|',
            'compressed_collective_H': 'PHP=(H_data+<psi|H_ref|psi>I_data) tensor |psi><psi|',
            'all_data_currents_survive_reference_only_Zeno_limit': True,
            'reference_survival_cost_alone_filters_the_data_Hamiltonian': False},
        'architecture_inputs': ['fresh bath blanks', 'source-derived biased-ground preparation word',
            'fresh native counter preparation', 'record-controlled SWAP', 'retained bath and counter storage'],
        'scope': {'decoded_two_object_reference_core': True,
            'geometric_four_carrier_embedding_is_upstream': True,
            'ground_source': 'native_response_biased.ground64(Fraction(16))',
            'prepared_ground_is_projected_geometric_ground': True,
            'full1296_color_ground_prepared_exactly': False,
            'passive_C_E_already_cools_or_selects_currents': False,
            'cooling_alone_selects_data_current_algebra': False,
            'new_bath_coupling_is_a_supplied_architecture': True}}


def metrics(arrays):
    _schema(arrays)
    ground = ground_vector().astype(complex)
    identity = np.eye(64, dtype=complex)
    unchanged = np.kron(ground[:, None], identity)
    replaced = np.kron(identity, ground[:, None])
    expected_labels = np.asarray([(coin, raw) for coin, raw, _ in raw_records()], np.int64)
    if not np.array_equal(arrays['labels'], expected_labels):
        raise ValueError('complete cooling counter records differ')
    code_zero = transport.PAIR_CODES['C'].ideal_code_matrix()[:, 0]
    target_counter = np.asarray([math.sqrt(float(w))*code_zero for _, _, w in raw_records()])
    errors = {'bath_preparation_column': float(np.linalg.norm(arrays['bath_preparation'][:, 0]-ground)),
        'bath_preparation_complete_Gram': float(np.linalg.norm(arrays['bath_preparation'].conj().T@arrays['bath_preparation']-identity, 2)),
        'unchanged_isometry': float(np.linalg.norm(arrays['unchanged']-unchanged, 2)),
        'replacement_isometry': float(np.linalg.norm(arrays['replaced']-replaced, 2)),
        'complete_native_counter': float(np.linalg.norm(arrays['counter']-target_counter))}
    if max(errors.values()) > SOURCE_ADMISSION:
        raise ValueError('complete source cooling isometry exceeds fixed component admission')
    if np.linalg.norm(arrays['unchanged']-np.kron(arrays['bath_preparation'][:, :1], identity)) > 2e-12:
        raise ValueError('the unknown input was not tensored with its actual prepared fresh bath')
    weights = np.sum(abs(arrays['counter'])**2, axis=1)
    gram = sum((weight*(arrays['replaced'] if (coin, raw) == TRIGGER else arrays['unchanged']).conj().T@
                (arrays['replaced'] if (coin, raw) == TRIGGER else arrays['unchanged'])
                for (coin, raw, _), weight in zip(raw_records(), weights, strict=True)), np.zeros((64, 64), complex))
    complete_gram = float(np.linalg.norm(gram-identity, 2))
    if complete_gram > 6*SOURCE_ADMISSION:
        raise ValueError('the complete unconditional cooling instrument is not normalized')
    return {'component_operator_errors': errors, 'complete_Gram_defect': complete_gram,
        'source_raw_record_probabilities': weights.tolist(),
        'source_trigger_probability': float(weights[1]),
        'unknown_reference_input_columns': 64, 'retained_reference_bath_rows': 4096,
        'retained_native_counter_rows': 512, 'raw_records': 8,
        'whole_tensor_output_materialized': False,
        'error_domain': 'all64 decoded reference inputs, with arbitrary untouched reference columns'}


def execute():
    prepare = compiler.compile_real_state(ground_vector())
    swap = swap_word()
    counter, labels, native_resources = _counter()
    bath = compiler.execute(prepare, 6)
    unchanged = np.kron(bath[:, :1], np.eye(64, dtype=complex))
    replaced = compiler.execute(swap, 12, unchanged)
    arrays = {'bath_preparation': bath, 'unchanged': unchanged, 'replaced': replaced,
              'counter': counter, 'labels': labels}
    own = metrics(arrays)
    code = transport.PAIR_CODES['C']
    recovery_edges = max(len(native_routing.permutation_edges(code.branch_permutation('C', COUNTER_SEAM, coin, raw)))
        for coin, raw, weight in raw_records() if weight)
    macros = (2*len(native_routing.permutation_edges(code.encoder_permutation()))
              +2*len(native_routing.permutation_edges(record_code._route_permutation(COUNTER_SEAM)))+recovery_edges)
    count = (len(prepare)+len(swap)+native_resources['encoding']['native_pulses']
             +native_resources['repair']['native_pulses']+native_resources['maximum_recovery_native_pulses'])
    export = reference._export_diagnostic([prepare, swap], count)
    projection = 5.1e-11
    if any(native_routing._executed_macro(inverse)[1]['projection_change'] > projection for inverse in (False, True)):
        raise ArithmeticError('raw native macro correction exceeds its fixed ceiling')
    raw_change = math.exp(count*math.log1p(2e-11))*math.expm1(macros*math.log1p(projection))
    success = (1+SOURCE_ADMISSION)**2-1+raw_change
    complete = success+export['first_failure_amplitude_bound']
    if not math.isfinite(complete) or complete > 1e-5:
        raise ArithmeticError('complete cooling history exceeds numerical admission')
    report = {'exact': exact_certificate(), 'execution': own,
        'resources': reference._resource_record({'ground_preparation': compiler.resources(prepare),
            'SWAP': compiler.resources(swap), 'native_counter': native_resources,
            'maximum_native_pulses_per_cycle': count, 'native_macro_invocations': macros,
            'supplied_unknown_reference_qubits': 6,
            'fresh_bath_qubits_per_cycle': 6, 'fresh_counter_work_qubits_per_cycle': 9,
            'fresh_coin_record_qubits_per_cycle': 3, 'committed_native_bits_per_cycle': 3,
            'first_export_failure_position_bits': count.bit_length(),
            'first_export_failure_flag_bits': 1,
            'first_export_failure_archive_bits': count.bit_length()+1,
            'first_export_failure_position_commit': 1,
            'maximum_archive_bits_including_export_failure': 3+count.bit_length()+1,
            'retained_quantum_bath_and_counter_qubits_per_cycle': 15,
            'old_reference_information_destroyed_or_reset_in_place': False}),
        'error_budget': {'component_source_admission': SOURCE_ADMISSION,
            'raw_macro_projection_change_ceiling': projection, 'raw_macro_adjustment': raw_change,
            'proper_code_exports': export, 'complete_flagged_operator_error': complete,
            'one_cycle_channel_error': complete*(2+complete),
            'n_cycle_channel_error': 'at most (1+delta)^n-1 for the floating CP maps, delta=one_cycle_channel_error; n*delta applies to genuine CPTP realizations; add the separate ideal (11/12)^n cooling error',
            'numerical_not_outward_rounded_or_physical': True}}
    return report, arrays


def compare_execution(report, arrays, fresh_report, fresh_arrays):
    if type(report) is not dict or report.keys() != fresh_report.keys():
        raise ValueError('complete cooling report differs')
    own = metrics(arrays)
    reference._owned_compare(own, report['execution'])
    reference._fixed_compare(fresh_report['error_budget'], report['error_budget'])
    reference._fixed_compare(fresh_report['resources'], report['resources'])
    _exact_compare(fresh_report['exact'], report['exact'])
    _schema(fresh_arrays)
    for key, actual in arrays.items():
        if key == 'labels':
            if not np.array_equal(actual, fresh_arrays[key]):
                raise ValueError('cooling raw labels differ')
        elif np.linalg.norm(actual-fresh_arrays[key]) > REPLAY_ADMISSION:
            raise ValueError('complete cooling source arrays differ: '+key)


def _exact_compare(expected, actual):
    if type(expected) is not type(actual):
        raise ValueError('cooling exact claim original type differs')
    if type(expected) is dict:
        if expected.keys() != actual.keys():
            raise ValueError('cooling exact claim keys differ')
        for key in expected:
            _exact_compare(expected[key], actual[key])
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise ValueError('cooling exact claim length differs')
        for a, b in zip(expected, actual, strict=True):
            _exact_compare(a, b)
    elif expected != actual:
        raise ValueError('cooling exact claim differs')


def source_bindings():
    from . import native_response_biased as biased
    paths = set(reference.source_bindings())|set(biased.source_bindings())|{
        'model/common_source/native_response_cooling.py', 'model/common_source/test_native_response_cooling.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}
