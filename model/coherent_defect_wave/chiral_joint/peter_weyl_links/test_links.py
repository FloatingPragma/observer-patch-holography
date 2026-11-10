"""Exact hardware, complete source propagation and rejected mutations."""
import pytest
import sympy as s

from .clebsch import coordinate, intertwiner, carrier_index
from .model import path_transitions, column
from .checks import (clebsch_checks, matrix_checks, path_checks,
                     rotation_word_checks, actual_hamiltonian_checks)
from .experiment import transfer_checks, preparation_checks, supported_bilinear, inner
from ..periodic_family.model import apply_state
from ..neutral_mobility.model import preparation


def test_exhaustive_current_derived_intertwiners_and_independent_wedges():
    row = clebsch_checks()
    assert row['nonzero_retained_channels'] == 37
    assert row['exact_scalar_commutants'] == 11
    assert len(row['independent_native_exterior_wedge_projectors']) == 5
    assert intertwiner(6, -3, 0)[1]['hypercharge_rejected']
    # Removing this retained transition would preserve Ward covariance but
    # fail the literal multiplication/completeness requirement.
    assert coordinate(6, 0, 0)[carrier_index(3, 1, 1), carrier_index(-3, 0, 0)] == 1


def test_complete_Gauss_orientation_and_norm_matrices():
    row = matrix_checks()
    assert row['complete_endpoint_Ward_identities'] == 1416
    assert [q['exact_Gram_projection_rank'] for q in row['blocks']] == [45, 72, 24, 80, 90]
    assert not row['link_excitation_parity_anticommutation']
    assert row['wrong_orientation_transpose_rejected']


def test_quantum_paths_and_complete_cubic_words():
    assert path_checks()['exhaustive_one_link_columns_both_adjoint_branches'] == 238
    assert rotation_word_checks()['proper_spatial_rotations'] == 24
    with pytest.raises(ValueError):
        path_transitions(())
    with pytest.raises(ValueError):
        coordinate(0, 0, 0)


def test_occupied_full_H_and_physical_reference():
    row = actual_hamiltonian_checks()
    assert row['complete_occupied_basis_rotation_controls'] == 33
    assert row['star_commutator_witness_norm_squared'] == '1/4'
    assert preparation_checks()['total_fermion_parity'] == 'even'


def test_supported_bilinear_is_full_H_action_on_a_small_physical_control():
    state = preparation()
    image = apply_state(state, column)
    actual, rows = supported_bilinear(image, state)
    assert actual == inner(image, apply_state(state, column)) == s.Rational(1651, 8)
    assert len(rows) == 33


def test_bare_transfer_with_complete_outcomes_and_full_H_remainder():
    row = transfer_checks()
    assert row['target_spin_moments'][1]['moments_orders_zero_through_five'] == ['0'] * 5 + ['-151*I/32']
    lower, upper = map(s.Rational, row['complete_arrival_probability_interval'])
    failure_lower, failure_upper = map(s.Rational, row['complete_no_arrival_probability_interval'])
    assert 0 < lower < upper < 1
    assert lower + failure_upper == upper + failure_lower == 1
    assert not row['extra_neutral_hopping_term']


def test_custody_rejects_missing_declared_and_import_inputs(monkeypatch):
    from . import certificate as custody
    with monkeypatch.context() as patch:
        patch.setattr(custody.action, 'SOURCE_FILES', (*custody.action.SOURCE_FILES, 'missing-PW-seed.json'))
        with pytest.raises(FileNotFoundError, match='missing declared'):
            custody.input_hashes()
    monkeypatch.setattr(custody, 'import_closure', lambda seeds: set(seeds) | {'missing-PW-import.py'})
    with pytest.raises(FileNotFoundError, match='missing Peter-Weyl import'):
        custody.input_hashes()
