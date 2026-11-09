"""Remove the coherent selector register by declaring a direct algebraic dilation.

The success operator is unchanged. Complete failure maps are different, and
the target operator's description and realization costs are not bounded here.
"""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.particle_common_contract import contract
from model.particle_filtered_probe import filter as apparatus
from model.particle_primary_survival import survival
from model.particle_primary_capacity import capacity

CLAIM = 'compact complete readout of the original finite filtered response'
SOURCE_INPUTS = (
    'model/particle_common_contract/contract.py',
    'model/particle_common_contract/PROOF.md',
    'model/particle_filtered_probe/filter.py',
    'model/particle_filtered_probe/PROOF.md',
    'model/particle_primary_survival/survival.py',
    'model/particle_primary_survival/PROOF.md',
    'model/particle_primary_mass_readout/PROOF.md',
    'model/particle_primary_capacity/capacity.py',
    'model/particle_primary_capacity/PROOF.md',
)


def compact_instrument(filtered_probe):
    """Reuse the exact Hermitian-contraction dilation with one neutral qubit."""
    b = apparatus.exact_matrix(filtered_probe)
    unitary = apparatus.hermitian_dilation(b)
    dimension = b.rows
    kraus = (unitary[:dimension, :dimension], unitary[dimension:, :dimension])
    if kraus[0] != b or apparatus.clean(sum((k.H*k for k in kraus), s.zeros(dimension))) != s.eye(dimension):
        raise ArithmeticError('compact filter lost its original success or completeness')
    return {'unitary': unitary, 'kraus': kraus, 'success': b,
            'fine_labels': ('original_probe', 'complement'),
            'selector_register_qubits': 0, 'neutral_probe_qubits': 1}


def compact_ramsey(hamiltonian, filtered_probe, duration='1', steps=16, imaginary=False):
    h = apparatus.exact_matrix(hamiltonian)
    instrument = compact_instrument(filtered_probe)
    if h.shape != instrument['success'].shape:
        raise ValueError('the filtered probe must act on the same source carrier')
    rows = apparatus.ramsey_instrument(h, instrument['kraus'], duration, steps, imaginary)
    if len(rows) != 4:
        raise ArithmeticError('both probe and Ramsey bits must be retained')
    return {'filter': instrument, 'rows': rows,
            'fine_outcomes': 4, 'neutral_apparatus_qubits': 2,
            'all_outcomes_retained': True, 'success_normalized': False}


def complete_error_bound(probe_operator_error='0', echo_unitary_error='0',
                         compiler_halfdiamond='0', preparation_halftrace='0'):
    """Compare two COMPACT instruments; this does not compare selector failures."""
    delta = min(F(2), apparatus.exact(probe_operator_error, 'probe operator error'))
    clock = min(F(2), apparatus.exact(echo_unitary_error, 'echo unitary error'))
    sigma = apparatus.exact(compiler_halfdiamond, 'whole compact-instrument error')
    prep = apparatus.exact(preparation_halftrace, 'preparation error')
    root = apparatus._relative_sqrt_upper(delta*delta+2*delta)
    total = min(F(1), root+2*clock+sigma+prep)
    return {'compact_column_error_formula': 'sqrt(delta^2+2delta)',
            'compact_column_error_rational_upper': root,
            'complete_halftrace_error_upper': total,
            'valid_for_all_untouched_references': True,
            'comparison_instruments': 'two direct compact dilations, with all four Ramsey outcomes',
            'equals_selector_complete_failure_law': False}


def conditional_compiler(system_bits, contrast_coefficient, floor_exponent):
    """Factored inherited four-outcome compiler at sigma=d/16, d=c*2^-N."""
    q = apparatus.integer(system_bits, 'source register bits')
    n = apparatus.integer(floor_exponent, 'density floor exponent', minimum=0)
    c = apparatus.exact(contrast_coefficient, 'contrast coefficient', True)
    sigma = survival.dyadic(c/16, n)
    if not survival.dyadic_leq(sigma, survival.dyadic(F(1, 2), 0)):
        raise ValueError('the declared sign-preserving allowance must be below one')
    total = q+2
    return {'padded_source_bits': q, 'retained_record_bits': 2,
            'padded_unitary_bits': total,
            'pulse_bound_exact_expression': f'20*2^{2*total}+4*2^{total}+10*{q}',
            'pulse_bound_factored_upper': {'coefficient': 34, 'power_of_two': 2*total},
            'port_fraction_bits': 2*total+n+apparatus._ceil_log2(1972*16/c),
            'rank_one_tree_angle_bits': total+n+apparatus._ceil_log2(8*(total+2)*16/c),
            'complete_halfdiamond_allowance': sigma,
            'positive_finite_contrast_after_compilation_lower': survival.dyadic(3*c/8, n),
            'additional_preparation_allowance': 0,
            'coherent_source_code_load_export_service_required': True,
            'gatewise_Gauss_preservation_of_compiled_tape': False,
            'all_compiler_failure_and_code_complement_records_retained': True,
            'compiler_workspace_capacity_bounded': False,
            'source_operator_description_cost_bounded': False,
            'native_spatial_tape_executed': False}


def carrier_budget(member):
    contract.verify_member(member)
    source = member['source']; volume = source['side']**3
    degree = source['per_link_word_degree']
    q = member['preparation_plan']['source']['dimension_bits_upper']
    # Source and its declared full conjugate reference, plus two neutral bits.
    terms = [(2, 120*volume+2), (451, 6*degree*volume)]
    exact = capacity.log_expression([(2, 120*volume+2), (119, 6*volume)]) if degree == 1 else None
    return {'source_plus_reference_dimension': 'D_source^2',
            'neutral_probe_qubits': 1, 'Ramsey_qubits': 1,
            'coherent_selector_register_qubits': 0,
            'retained_fine_records_per_phase': 4,
            'one_shot_joint_carrier_dimension': '4 D_source^2',
            'one_shot_joint_capacity_nats': '2 ln D_source + 2 ln 2',
            'one_shot_joint_capacity_exact': exact,
            'one_shot_joint_capacity_upper': capacity.log_expression(terms),
            'padded_source_qubit_upper': q,
            'padded_source_reference_and_apparatus_qubit_upper': 2*q+2,
            'controller_precision_workspace_and_repeated_archive_included': False,
            'conditional_compiler_extra_service_flags_included': False,
            'total_self_contained_executable_capacity_proved': False}


def common_contract_payload(member):
    contract.verify_member(member)
    witness = survival.selected_witness(member)
    d = witness['continuous_original_contrast_lower']
    budget = carrier_budget(member)
    return contract.serial({
        'program': contract.PROGRAM, 'member_id': member['member_id'], **member['ids'],
        'claim': CLAIM, 'probe_descriptor': witness['probe'],
        'inherited_finite_filter_witness_sha256': contract.digest(witness),
        'target_operator': 'the identical finite Cayley/quadrature B_tilde',
        'quadrature_nodes_in_operator_definition': witness['apparatus']['quadrature_cells'],
        'Cayley_subdivisions_in_operator_definition': witness['apparatus']['Cayley_subdivisions'],
        'success_Kraus_operator': 'B_tilde',
        'complement_Kraus_operator': 'sqrt(I-B_tilde^2)',
        'apparatus_service': 'direct controlled exact algebraic Hermitian-contraction dilation of B_tilde',
        'source_energy_zero_changed': False,
        'reference_idle_H': member['reference']['idle_H'],
        'reference_actively_counter_evolved': False,
        'same_original_F0_ReF_ImF_as_selector_apparatus': True,
        'same_blank_Ramsey_Kraus_operators_on_all_inputs': True,
        'complete_failure_channel_equal_to_selector_apparatus': False,
        'autonomous_repair_compaction_proved': False,
        'finite_apparatus_original_contrast_lower': witness['finite_apparatus_original_contrast_lower'],
        'original_transfer_band': witness['original_positive_band'],
        'original_transfer_and_final_excitation_band_weight_lower': witness[
            'finite_apparatus_original_positive_band_weight_lower'],
        'rest_mass_or_fixed_band_noncollapse_proved': False,
        'carrier_budget': budget,
        'conditional_compiler': conditional_compiler(budget['padded_source_qubit_upper'],
            d['coefficient'], d['negative_exponent']),
        'large_primary_operator_or_root_materialized': False,
        'exact_service_description_precision_or_gate_cost_bounded': False,
        'efficient_implementation_claimed': False,
    })


def common_contract_bound(member):
    result = common_contract_payload(member)
    return {'schema': 'oph.compact_primary_mass_probe_bound.v1', 'program': contract.PROGRAM,
            'member_id': member['member_id'], **member['ids'], 'claim': CLAIM,
            'result': result, 'result_sha256': contract.digest(result)}


def verify_bound(member, value):
    if value != common_contract_bound(member):
        raise ValueError('compact apparatus or original source binding differs')


@lru_cache(maxsize=1)
def quantum_control():
    """Execute separate two-level source instruments with different failures."""
    h = s.diag(0, 4); rho = s.diag(s.Rational(3, 4), s.Rational(1, 4))
    a = s.Matrix([[0, 1], [1, 0]])
    original = apparatus.contraction_instrument(h, a, '1/2', 2, 1)
    b = original['success']
    if b != 3*a/5:
        raise ArithmeticError('the independent finite filtered operator changed')
    compact = compact_instrument(b)
    source_input = s.diag(1, 0)
    old_failure = apparatus.clean(sum((k*source_input*k.H for k in original['kraus'][1:]), s.zeros(2)))
    direct_failure = apparatus.clean(compact['kraus'][1]*source_input*compact['kraus'][1].H)
    if old_failure == direct_failure:
        raise ArithmeticError('the control must distinguish the complete failure laws')
    phases = []
    for imaginary in (False, True):
        old = apparatus.ramsey_instrument(h, original['kraus'], '1/2', 1, imaginary)
        direct = compact_ramsey(h, b, '1/2', 1, imaginary)['rows']
        if any(old[j]['K'] != direct[j]['K'] for j in (0, 1)):
            raise ArithmeticError('original blank Ramsey maps differ')
        probabilities = [s.simplify(s.trace(rho*row['K'].H*row['K'])) for row in direct]
        old_probabilities = [s.simplify(s.trace(rho*row['K'].H*row['K'])) for row in old]
        score = s.simplify(sum(row['original_score']*p for row, p in zip(direct, probabilities, strict=True)))
        if sum(probabilities) != 1 or score != (s.Rational(9, 50) if imaginary else 0):
            raise ArithmeticError('the complete compact response is not the original score')
        phases.append({'imaginary': imaginary, 'original_probabilities': [str(p) for p in old_probabilities],
                       'compact_probabilities': [str(p) for p in probabilities], 'original_score': str(score)})
    return {'H': h, 'rho': rho, 'A': a, 'B': b, 'original': original,
            'compact': compact, 'old_failure_on_zero': old_failure,
            'compact_failure_on_zero': direct_failure,
            'receipt': {'source_dimension': 2, 'reference_dimension': 2,
                'is_primary_source': False, 'old_Ramsey_fine_outcomes': 8,
                'compact_Ramsey_fine_outcomes': 4, 'F0': '9/25', 'ReF': '0', 'ImF': '9/50',
                'all_original_blank_Kraus_maps_equal': True,
                'coarse_failure_quantum_channels_differ': True, 'phases': phases}}


def certificate():
    return contract.serial({'schema': 'oph.compact_primary_mass_probe.v1',
        'selected_members': [common_contract_bound(contract.member(n)) for n in (1, 2)],
        'independent_quantum_control': quantum_control()['receipt'],
        'total_executable_capacity_or_efficient_realization_proved': False})
