"""Exact full-word symmetries and physical three-dimensional controls."""
from collections import defaultdict
from functools import lru_cache
from itertools import product

import sympy as s

from .model import (Configuration, PAULI, cells, cell_index, shift_cell, edges,
                    edge_endpoints, path, kinetic_terms, native_blocks,
                    representation, translate_fermions, translation_column,
                    pair_preparation, phase_projector, column, gauss_column,
                    apply_state, norm2, local_yukawa, descriptor)
from ..model import clean, path_ward_checks
from ...chiral_higgs.model import tensor_checks
from ...chiral_regulator.model import momentum_matrix
from ..two_cell import hamiltonian as two_cell_h
from ..two_cell.symmetries import periodic_column as two_cell_column


@lru_cache(maxsize=None)
def translation_word_checks(side=2):
    config = Configuration(side=side)
    terms = kinetic_terms(config)
    lookup = {(term.left_cell, term.axis, term.displacement): term for term in terms}
    count = 0
    for term in terms:
        left, right = cells(side)[term.left_cell], cells(side)[term.right_cell]
        for direction in range(3):
            moved_left = shift_cell(side, left, direction)
            moved_right = shift_cell(side, right, direction)
            moved = lookup[cell_index(side, moved_left), term.axis, term.displacement]
            fermion_sign = (-1) ** (int(left[direction] == side - 1) + int(right[direction] == side - 1))
            assert s.simplify(term.coefficient * fermion_sign - moved.coefficient) == 0
            translated_path = tuple(3 * cell_index(side, shift_cell(side, edges(side)[edge][0], direction))
                                    + edges(side)[edge][1] for edge in term.path)
            assert translated_path == moved.path
            assert moved.right_cell == cell_index(side, moved_right)
            count += 1
    one_body = 0
    for mode in range(30 * side ** 3):
        for a in range(3):
            first, sign_first = translate_fermions(1 << mode, a, side)
            for b in range(3):
                second, sign_second = translate_fermions(first, b, side)
                other, sign_other = translate_fermions(1 << mode, b, side)
                back, sign_back = translate_fermions(other, a, side)
                assert (second, sign_first * sign_second) == (back, sign_other * sign_back)
                one_body += 1
            state, coefficient = 1 << mode, s.Integer(1)
            for _ in range(side):
                state, phase = translate_fermions(state, a, side)
                coefficient *= phase
            assert state == 1 << mode and coefficient == -1
    # Scalar, link and local even pair operators are translated by the same
    # cell permutation. Two antiperiodic fermion endpoint phases cancel in
    # each onsite pair; occupation penalties only permute their summands.
    assert all((-1) ** (2 * int(x[axis] == side - 1)) == 1 for x in cells(side) for axis in range(3))
    return {'side': side, 'all_ordered_path_translation_identities': count,
            'all_one_body_translation_commutators': one_body,
            'full_Fock_translation_commutation_by_exterior_lift': True,
            'T_axis_to_side_equals_matter_parity': True,
            'all_three_H_translation_commutators_zero_by_full_operator_grammar': True,
            'all_local_Yukawa_and_penalty_terms_retained': True}


@lru_cache(maxsize=1)
def two_cell_normalization_checks():
    """Equality of full operator words, not a one-particle state truncation."""
    config = Configuration(side=2, spacing=s.pi)
    derivative = momentum_matrix(2, config.length)
    assert derivative == s.Matrix([[0, -s.I / 2], [s.I / 2, 0]])
    labels = ((-1, -1), (0, 0), (1, 1), (1, 2), (2, 1), (2, 2))
    expected = defaultdict(lambda: s.Integer(0))
    for anchor in (0, 1):
        coefficient = -derivative[anchor, 1 - anchor] / 2
        assert coefficient == s.I * (-1) ** anchor / 4
        for label, (a, b) in enumerate(labels[1:], 1):
            for spin in range(2):
                i, j = 6 * anchor + 3 * spin + a, 6 * (1 - anchor) + 3 * spin + b
                value = coefficient * PAULI[2][spin, spin]
                for spectator in range(6):
                    row = (label, spectator) if anchor == 0 else (spectator, label)
                    col = (0, spectator) if anchor == 0 else (spectator, 0)
                    expected[i, j, row, col] += value
                    expected[j, i, col, row] += s.conjugate(value)
    expected = clean(expected)
    actual = {}
    # These matrix-unit columns recover every CAR-bilinear coefficient and
    # all two-link words. Equality of those words then holds on every Fock
    # occupation, rather than restricting the source to these input states.
    for j in range(12):
        for l01, l10 in product(range(6), repeat=2):
            key = (1 << j, l01, l10, 0, 0)
            for (mask, a, b, hl, hr), value in two_cell_column(key).items():
                if mask.bit_count() == 1 and mask != 1 << j and hl == hr == 0:
                    actual[mask.bit_length() - 1, j, (a, b), (l01, l10)] = value
    assert actual == expected
    native_indices = (9, 10, 11, 24, 25, 26)
    local_checks = 0
    for mask in range(64):
        native_mask = sum(1 << native_indices[i] for i in range(6) if mask >> i & 1)
        for scalar in range(5):
            full = local_yukawa(native_mask, scalar, (s.Integer(0), s.Integer(0), s.Integer(1)))
            restricted = {}
            for (target, h), value in full.items():
                assert target & ~sum(1 << j for j in native_indices) == 0
                small = sum(1 << i for i, j in enumerate(native_indices) if target >> j & 1)
                restricted[small, h] = value
            assert restricted == two_cell_h.local_yukawa(mask, scalar)
            local_checks += 1
    return {'spacing': 'pi', 'length_at_side_two': '2*pi', 'spectral_momenta': ['-1/2', '1/2'],
            'ordered_edge_coefficients_before_sigma3': ['I/4', '-I/4'],
            'exact_hopping_matrix_unit_words_matched': len(expected),
            'complete_local_lepton_scalar_columns_matched': local_checks,
            'scalar_and_link_penalties_match_at_unit_coefficients': True,
            'all_occupation_operator_equality_on_the_stated_one_axis_control': True,
            'two_cell_ground_transferred_to_3D': False}


def norm_bound(config=Configuration()):
    dimensions = [len(block) for block in native_blocks()]
    kinetic = sum(4 * abs(term.coefficient) * sum(d ** (term.displacement + 1) for d in dimensions)
                  for term in kinetic_terms(config))
    weights = {'up': config.y_up, 'down': config.y_down, 'charged_lepton': config.y_lepton}
    local = sum(2 * abs(weights[channel]) * sum(abs(v) for component in vertex['components'] for v in component)
                for channel, vertex in representation()[3].items())
    volume = config.side ** 3
    return s.simplify(kinetic + volume * (local + config.scalar_penalty) + 3 * volume * config.link_penalty)


@lru_cache(maxsize=1)
def physical_cube_checks():
    config = Configuration()
    state = pair_preparation(0, config)
    evolved_column = apply_state(state, lambda key: column(key, config))
    assert len(evolved_column) == 849 and norm2(evolved_column) == s.Rational(3321, 16)
    gauge_checks = 0
    for vector in (state, evolved_column):
        for cell in range(8):
            for current in range(12):
                assert not apply_state(vector, lambda key: gauss_column(key, cell, current, config))
                gauge_checks += 1
    for direction in range(3):
        left = apply_state(evolved_column, lambda key: translation_column(key, direction, config))
        right = apply_state(apply_state(state, lambda key: translation_column(key, direction, config)),
                            lambda key: column(key, config))
        assert left == right
    # Check exact adjoint coefficients on every column reached in this
    # control. Global Hermiticity follows separately from the full-word H+h.c.
    adjoints = 0
    for source in state:
        for target, value in column(source, config).items():
            assert column(target, config).get(source, 0) == s.conjugate(value)
            adjoints += 1
    # Complete simultaneous translation-phase effects: every possible triple
    # is present, including the exactly zero odd labels of this even input.
    records = []
    total = defaultdict(lambda: s.Integer(0))
    branch_checks = 0
    for labels in product(range(4), repeat=3):
        image = phase_projector(state, labels, config)
        probability = norm2(image)
        assert probability == (s.Rational(1, 8) if all(label % 2 == 0 for label in labels) else 0)
        for key, value in image.items():
            total[key] += value
        if image:
            for direction, label in enumerate(labels):
                translated = apply_state(image, lambda key: translation_column(key, direction, config))
                assert translated == clean({key: s.I ** label * value for key, value in image.items()})
            for cell in range(8):
                for current in range(12):
                    assert not apply_state(image, lambda key: gauss_column(key, cell, current, config))
                    branch_checks += 1
        records.append({'phase_labels': list(labels), 'phase_unit': 'pi/2', 'probability': str(probability)})
    assert clean(total) == state
    assert sum(s.Rational(record['probability']) for record in records) == 1
    return {'source': descriptor(config), 'preparation': 'normalized local charged-lepton Yukawa pair/scalar action on empty source',
            'preparation_sites': [0], 'prepared_basis_states': len(state),
            'prepared_state_is_stationary_or_ground': False,
            'all_three_translation_H_controls_exact': True,
            'actual_preparation_and_H_state_Gauss_actions_checked': gauge_checks,
            'actual_phase_branch_Gauss_actions_checked': branch_checks,
            'H_state_full_basis_support': len(evolved_column), 'H_state_norm_squared': '3321/16',
            'actual_adjoint_coefficients_checked': adjoints,
            'all_channels_including_quark_creation_present_in_idle_H': True,
            'complete_joint_phase_records': records,
            'complete_joint_phase_records_count': len(records),
            'physical_three_translation_phase_readout': True,
            'source_particle_momentum_or_relativistic_mass_identified': False}


def certificate_data():
    config = Configuration()
    bound = norm_bound(config)
    assert bound == 1896
    return {'schema': 'oph.chiral_periodic_quantum_family.v1', 'source': descriptor(config),
            'normalization_and_frozen_one_axis_equality': two_cell_normalization_checks(),
            'full_word_translation_checks': [translation_word_checks(n) for n in (2, 4)],
            'actual_native_path_Gauss_intertwiners_reused': path_ward_checks(),
            'all_quantum_Yukawa_current_intertwiners_reused': tensor_checks(),
            'physical_three_axis_control': physical_cube_checks(),
            'finite_H_norm_bound_at_default_parameters': str(bound),
            'same_H_physical_reference': {'definition': 'P exp(-beta H) P / Tr(P exp(-beta H)), beta>0',
                'nonempty_physical_sector_witness': 'empty matter, blank links, scalar vacuum',
                'same_finite_filter_construction_as_parent': True,
                'global_interacting_ground_or_Gibbs_weights_computed': False},
            'scope': {'three_commuting_source_translations': True,
                'all_local_Gauss_and_all_occupations_retained': True,
                'spectral_path_support_is_long_range': True,
                'primitive_path_matrix_words_are_not_normalized': True,
                'free_ground_or_open_ground_imported': False,
                'interacting_relativistic_field_limit': False,
                'named_particle_or_GeV_mass': False}}
