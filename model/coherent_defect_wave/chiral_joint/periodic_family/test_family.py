"""Exact full-word and physical three-dimensional source controls."""
import pytest
import sympy as s

from .model import Configuration, descriptor, kinetic_terms, native_blocks
from .checks import (translation_word_checks, two_cell_normalization_checks,
                     norm_bound, physical_cube_checks)
from .certificate import input_hashes


def test_normalization_reproduces_actual_frozen_one_axis_operator():
    result = two_cell_normalization_checks()
    assert result['spectral_momenta'] == ['-1/2', '1/2']
    assert result['ordered_edge_coefficients_before_sigma3'] == ['I/4', '-I/4']
    assert result['exact_hopping_matrix_unit_words_matched'] == 240
    assert result['complete_local_lepton_scalar_columns_matched'] == 320
    assert not result['two_cell_ground_transferred_to_3D']


def test_three_full_word_translations_preserve_actual_paths_and_AP_signs():
    for side, expected in ((2, 72), (4, 1728)):
        result = translation_word_checks(side)
        assert result['all_ordered_path_translation_identities'] == expected
        assert result['all_three_H_translation_commutators_zero_by_full_operator_grammar']
        assert result['T_axis_to_side_equals_matter_parity']
    assert max(len(term.path) for term in kinetic_terms(Configuration(4))) == 3
    assert [len(block) for block in native_blocks()] == [3, 6, 1, 2, 3]
    assert norm_bound() == 1896
    with pytest.raises(ValueError):
        Configuration(2, spacing=1.0)


def test_actual_physical_cube_has_complete_three_phase_measurement():
    result = physical_cube_checks()
    assert result['source']['full_dimension'] == '2^240 * 60^24 * 5^8'
    assert result['H_state_full_basis_support'] == 849
    assert result['H_state_norm_squared'] == '3321/16'
    assert result['actual_preparation_and_H_state_Gauss_actions_checked'] == 192
    assert result['actual_phase_branch_Gauss_actions_checked'] == 768
    assert result['actual_adjoint_coefficients_checked'] == 852
    records = result['complete_joint_phase_records']
    assert len(records) == 64
    assert sum(s.Rational(row['probability']) for row in records) == 1
    assert sum(row['probability'] == '0' for row in records) == 56
    assert not result['prepared_state_is_stationary_or_ground']
    assert not result['source_particle_momentum_or_relativistic_mass_identified']


def test_source_custody_includes_dynamic_data_and_finite_theorem():
    from model.common_source.chiral_candidate import action
    hashes = input_hashes()
    assert set(action.SOURCE_FILES) <= set(hashes)
    assert 'Lean/QFT/PeriodicMomentum.lean' in hashes
    assert 'Lean/lake-manifest.json' in hashes
    assert not descriptor()['primitive_quantum_path_normalized_or_replaced_by_identity']
