"""Complete neutral Fourier pump and its retained source/reference phases."""
from itertools import product
import json

import sympy as s

from ...chiral_higgs.model import higgs_memory
from ...chiral_source.model import PAULI
from ..higgs_composites import model as car
from ..higgs_composites.certificate import validate_pins
from ..neutral_mobility.checks import annihilator, one_body
from ..ground_composite_response.spectrum import projector_coefficients
from ..periodic_family.model import Configuration, cells, cell_index, shift_cell, apply_state
from ..peter_weyl_links import model as source
from ..neutral_mobility.model import composite_column
from .model import ROOT, reference_data, momentum_probability

GROUND_RESPONSE_RECEIPT = 'model/coherent_defect_wave/chiral_joint/ground_composite_response/results.json'


def full_local_operator_checks():
    fields = higgs_memory()[0]
    c = tuple(annihilator(2, j) for j in range(2))
    j = s.ImmutableSparseMatrix(s.kronecker_product(c[1], fields[0])
                               - s.kronecker_product(c[0], fields[1]))
    anti = j.H * j + j * j.H
    assert anti.H == anti and anti.eigenvals() == {s.Integer(0): 10, s.Integer(1): 4, s.Integer(2): 6}
    assert j * j != s.zeros(20)
    assert (j.H * j).eigenvals() == {s.Integer(0): 14, s.Integer(1): 4, s.Integer(2): 2}
    # Both source spins and the actual quantum scalar are retained here.
    four = tuple(annihilator(4, k) for k in range(4))
    js = tuple(s.ImmutableSparseMatrix(s.kronecker_product(four[2 * spin + 1], fields[0])
        - s.kronecker_product(four[2 * spin], fields[1])) for spin in (0, 1))
    parity = s.ImmutableSparseMatrix(s.kronecker_product(
        s.diag(*[(-1) ** mask.bit_count() for mask in range(16)]), s.eye(5)))
    up, down = s.Matrix([1, 0]), s.Matrix([0, 1])
    pump = s.ImmutableSparseMatrix((s.kronecker_product(js[0].H * parity, down)
        - s.kronecker_product(js[1].H * parity, up)) / (2 * s.sqrt(2)))
    effect = pump.H * pump
    assert effect == (js[0] * js[0].H + js[1] * js[1].H) / 8
    eigenvalues = effect.eigenvals()
    assert all(0 <= value <= s.Rational(1, 2) for value in eigenvalues)
    for axis in range(3):
        spin = s.kronecker_product(one_body(four, s.kronecker_product(PAULI[axis] / 2, s.eye(2))), s.eye(5))
        total = s.kronecker_product(spin, s.eye(2)) + s.kronecker_product(s.eye(80), PAULI[axis] / 2)
        assert total * pump == pump * spin
    total_parity = s.kronecker_product(parity, -s.eye(2))
    assert total_parity * pump == pump * parity
    scalar_sum = car.add(*(car.multiply(car.neutral_composite(spin),
        car.adjoint(car.neutral_composite(spin))) for spin in (0, 1)))
    for current in range(12):
        assert not car.commutator(car.gauge_current(current), scalar_sum)
        for spin in (0, 1):
            assert not car.commutator(car.gauge_current(current), car.neutral_composite(spin))
    for axis in range(3):
        assert not car.commutator(car.spin_current(axis), scalar_sum)
    return {'full_active_fixed_spin_CAR_scalar_dimension': 20,
        'actual_J_anticommutator_spectrum': {'0': 10, '1': 4, '2': 6},
        'actual_JdaggerJ_spectrum': {'0': 14, '1': 4, '2': 2},
        'J_squared_is_zero': False, 'uniform_Fourier_anticommutator_bound': '2 I',
        'uniform_Fourier_J_operator_norm_bound': 'sqrt(2)',
        'uniform_Fourier_bound_uses_offsite_graded_anticommutation': True,
        'complete_active_two_spin_source_dimension': 80,
        'spin_singlet_success_map_rows': 160,
        'local_pump_success_effect_spectrum': {str(value): multiplicity for value, multiplicity in eigenvalues.items()},
        'pump_formula': '(J_p,up^dagger r_down^dagger-J_p,down^dagger r_up^dagger)/(2 sqrt(2))',
        'success_map_domain': 'complete source carrier tensor pump-reference vacuum',
        'success_map_codomain': 'complete source carrier tensor pump-reference one-fermion sector',
        'BdaggerB_formula_is_vacuum_input_compression': True,
        'other_reference_input_sectors': 'retained in the block unitary; unused double occupancy receives identity',
        'full_CAR_pump_matrix_convention': 'J_p^dagger times source parity before creating the following pump-reference mode',
        'uniform_pump_effect_bound': '0 <= B^dagger B <= I/2',
        'failure_operator': 'sqrt(I-B^dagger B)',
        'unitary_completion': '[[sqrt(I-BdaggerB),-Bdagger],[B,sqrt(I-BBdagger)]] plus identity on unused reference sectors',
        'three_exact_spin_intertwiners': 3, 'total_even_parity_intertwiner': True,
        'all_twelve_native_Gauss_currents_preserved': True,
        'local_ground_restricted_spin_sum_is_scalar': 'Schur lemma on the unique spin-half ground doublet'}


def actual_composite_translation_checks():
    count = 0
    for side in (2, 4):
        config = Configuration(side=side); base = source.empty_key(config)
        for axis in range(3):
            for site in (0, cell_index(side, tuple(side - 1 if j == axis else 0 for j in range(3)))):
                xyz = cells(side)[site]; moved = cell_index(side, shift_cell(side, xyz, axis))
                seam = -1 if xyz[axis] == side - 1 else 1
                for spin in (0, 1):
                    for occupied in range(4):
                        mask = (1 << 0) | (1 << (30 * ((site + 1) % (side ** 3)) + 29))
                        for j in range(2):
                            if occupied >> j & 1: mask |= 1 << (30 * site + 15 * spin + 10 + j)
                        for scalar in range(5):
                            scalars = list(base[2]); scalars[site] = scalar
                            key = mask, base[1], tuple(scalars)
                            for create in (False, True):
                                left = apply_state(composite_column(key, site, spin, create=create, config=config),
                                    lambda value: source.translation_column(value, axis, config))
                                right = apply_state(source.translation_column(key, axis, config),
                                    lambda value: composite_column(value, moved, spin, create=create, config=config))
                                assert left == {target: seam * value for target, value in right.items()}
                                count += 1
    return {'actual_full_carrier_neutral_operator_translation_controls': count,
        'operator': 'J_s and its adjoint, with the full quantum scalar and CAR signs',
        'Fourier_creation_convention': 'J_theta,s^dagger=N^(-3/2) sum_x exp(-i theta.x) J_x,s^dagger',
        'operator_translation_character': 'T_i J_theta^dagger T_i^dagger=exp(i theta_i) J_theta^dagger',
        'allowed_odd_transfer_phases': 'theta_i=pi*(2m_i+1)/N in every direction',
        'even_or_mixed_parity_grid_substitution_allowed': False}


def success_weight():
    reference_data()
    receipt = json.loads((ROOT / GROUND_RESPONSE_RECEIPT).read_text())
    validate_pins(receipt['input_sha256'])
    actual_source = receipt['source']
    assert actual_source['full_CAR_modes'] == 30 and actual_source['scalar_states'] == 5
    assert actual_source['coefficients'] == {
        'up': '1', 'down': '1', 'charged_lepton': '1', 'scalar_penalty': '1'}
    assert not actual_source['spatial_or_composite_hopping_added']
    data = receipt['actual_ground_composite_spectrum']
    assert data['probe_spin_independence_on_rotation_invariant_reference']
    assert data['exact_charge_conjugate_addition_removal_weights_equal']
    setting = next(row for row in data['recorded_settings'] if row['kind'] == 'add' and row['probe_spin'] == 0)
    other = next(row for row in data['recorded_settings'] if row['kind'] == 'add' and row['probe_spin'] == 1)
    assert setting['setting_probability'] == '1/4' and setting['success_operator'] == 'J_kind,spin/2'
    assert setting['exact_cross_moments_target_rows_source_columns'] == other['exact_cross_moments_target_rows_source_columns']
    x = s.Symbol('x')
    polynomial = s.Poly(s.sympify(setting['exact_source_minimal_polynomial']), x)
    factor = s.Poly((x * (x - 1)) ** 3 - 313 * (x * (x - 1)) ** 2
                   + 31708 * x * (x - 1) - 1041468, x)
    assert polynomial.rem(factor).is_zero
    moments = list(map(s.Rational, setting['exact_cross_moments_target_rows_source_columns'][0]))
    numerator = s.Poly(sum(c.as_expr() * value for c, value in
        zip(projector_coefficients(polynomial, factor), moments, strict=True)), x).rem(factor)
    denominator = 8 * polynomial.diff().rem(factor)
    lo, hi = (4 * s.Rational(value) for value in setting['joint_success_probability_interval'])
    assert 0 < lo <= hi < s.Rational(1, 2)
    return {'exact_success_weight': '(' + str(numerator.as_expr()) + ') / (' + str(denominator.as_expr()) + ')',
        'success_weight_numerator_polynomial': str(numerator.as_expr()),
        'success_weight_denominator_polynomial': str(denominator.as_expr()),
        'ground_root_factor': str(factor.as_expr()), 'ground_root_selection': 'least real root',
        'per_setting_ground_receipt_weight_scale': 'trace(P_ground J Jdagger)/32',
        'coherent_spin_singlet_pump_weight_scale': 'trace(P_ground J Jdagger)/8',
        'explicit_per_setting_to_pump_factor': '4',
        'success_probability_interval': [str(lo), str(hi)],
        'failure_probability_interval': [str(1 - hi), str(1 - lo)],
        'uniform_in_size_and_transfer_phase': True,
        'restricted_atomic_effect': 'P_atomic Bdagger B P_atomic = w P_atomic',
        'offsite_terms_vanish_after_atomic_restriction': 'each J flips the local odd fermion parity',
        'full_effect_assumed_to_preserve_atomic_subspace': False}


def complete_pump_records():
    weight = success_weight()
    lo, hi = map(s.Rational, weight['success_probability_interval'])
    settings = []
    for side in (2, 4, 6):
        transfer = (1, 1, 1)  # Full 2N phase indices, hence AP in every direction.
        rows = []
        total = s.Integer(0)
        for momentum in product(range(side), repeat=3):
            probability = momentum_probability(side, momentum); total += probability
            source_even = tuple(2 * k for k in momentum)
            ref = tuple((-2 * k) % (2 * side) for k in momentum)
            source_odd = tuple((k + ell) % (2 * side) for k, ell in zip(source_even, transfer, strict=True))
            assert all(k % 2 == 1 for k in source_odd)
            assert all((a + b) % (2 * side) == ell for a, b, ell in zip(source_odd, ref, transfer, strict=True))
            rows.append({'background_N_phase_indices': list(momentum),
                'background_probability': str(probability),
                'exact_joint_success_probability': 'w * (' + str(probability) + ')',
                'exact_joint_failure_probability': '(1-w) * (' + str(probability) + ')',
                'source_success_2N_phase_indices': list(source_odd),
                'source_failure_2N_phase_indices': list(source_even),
                'site_reference_2N_phase_indices': list(ref),
                'joint_success_probability_interval': [str(probability * lo), str(probability * hi)],
                'joint_failure_probability_interval': [str(probability * (1 - hi)), str(probability * (1 - lo))]})
        assert total == 1
        settings.append({'side': side, 'pump_transfer_2N_phase_indices': list(transfer),
            'pump_reference_success_2N_phase_indices': [(-ell) % (2 * side) for ell in transfer],
            'pump_reference_failure_phase_indices': [0, 0, 0],
            'source_plus_site_reference_success_phase_indices': list(transfer),
            'including_compensating_pump_reference_total_phase_indices': [0, 0, 0],
            'complete_record_count_including_zero_outcomes': 2 * (2 * side) ** 9,
            'positive_record_count': 2 * side ** 3,
            'complete_positive_records': rows,
            'all_other_outcomes': {'count': 2 * (2 * side) ** 9 - 2 * side ** 3,
                'probability': '0', 'retained': True}})
    return {'ground_response_ancestry': weight,
        'site_reference': 'one static neutral spin-half fermionic mode pair per cell, with one reference fermion',
        'pump_reference': 'retained neutral spin-half fermionic pair, opposite raw lepton charge and inverse transfer character',
        'source_success_parity': 'odd', 'source_failure_parity': 'even',
        'total_source_and_references_parity': 'even',
        'source_success_T_to_N': '-I', 'source_failure_T_to_N': 'I',
        'reference_character_cancels_total_phase_without_changing_source_grid': True,
        'all_probabilities_are_unconditional': True,
        'source_momentum_sharp_after_pump': False,
        'momentum_transfer_relative_to_calibrated_reference_is_sharp': True,
        'source_and_reference_phase_records_are_preserved_by_full_PW_quench': True,
        'energy_or_relativistic_dispersion_inferred': False, 'executed_sizes': settings}
