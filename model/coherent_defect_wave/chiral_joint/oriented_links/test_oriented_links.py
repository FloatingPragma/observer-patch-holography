"""Full-carrier identities and deliberately broken orientation controls."""
import pytest
import sympy as s

from ..periodic_family.model import Configuration
from .model import link, empty_key, validate_key, path_transitions
from .checks import (matrix_checks, peter_weyl_check, path_checks,
                     rotation_word_checks, actual_hamiltonian_checks)


def test_complete_native_current_and_orientation_matrices():
    result = matrix_checks()
    assert result['complete_matrix_endpoint_Ward_identities'] == 1416
    assert result['missing_transpose_mutation_rejected']
    assert all(row['transporter_operator_norm'] == '1' for row in result['blocks'])


def test_declared_architecture_is_distinct_from_literal_group_compression():
    result = peter_weyl_check()
    assert result['retained_nonvacuum_Peter_Weyl_transition_rank'] == 4
    assert result['declared_link_transition_rank_on_that_input'] == 0
    assert not result['literal_Peter_Weyl_compression']


def test_paths_retain_quantum_absorption_and_bounded_full_operator_norm():
    result = path_checks()
    assert result['mixed_creation_absorption_word_terms'] == 6
    assert result['hopping_only_response_controls'][1]['exact_blank_one_fermion_response_squared'] == '15/32'
    assert not result['local_star_parity_obstruction_removed']


def test_complete_periodic_words_have_all_proper_cubic_rotations():
    result = rotation_word_checks()
    assert result['proper_spatial_rotations'] == 24
    assert result['missing_antiperiodic_seam_mutation_rejected']
    assert result['rotation_does_not_conjugate_matter_internal_representation']


def test_actual_H_uses_nonblank_flux_and_unequal_native_couplings():
    result = actual_hamiltonian_checks()
    assert result['complete_basis_rotation_controls'] == 33
    assert result['actual_quantum_H_G_controls'] == 3
    assert result['actual_quantum_H_translation_controls'] == 3


def test_carrier_and_block_validation():
    with pytest.raises(ValueError):
        link(0, 9)
    with pytest.raises(ValueError):
        path_transitions(())
    with pytest.raises(ValueError):
        path_transitions((119,))
    mask, links, scalars = empty_key()
    with pytest.raises(ValueError):
        validate_key((mask, (119,) + links[1:], scalars))
    assert Configuration(y_up=s.Rational(2, 3)).y_up == s.Rational(2, 3)


def test_custody_rejects_missing_declared_and_import_inputs(monkeypatch):
    from . import certificate as custody
    with monkeypatch.context() as patch:
        patch.setattr(custody.action, 'SOURCE_FILES', (*custody.action.SOURCE_FILES, 'missing-oriented-seed.json'))
        with pytest.raises(FileNotFoundError, match='missing declared'):
            custody.input_hashes()
    monkeypatch.setattr(custody, 'import_closure', lambda seeds: set(seeds) | {'missing-oriented-import.py'})
    with pytest.raises(FileNotFoundError, match='missing oriented-link import'):
        custody.input_hashes()
