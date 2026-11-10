"""Physical preparation, full source motion and retained complementary records."""
import pytest
import sympy as s

from .model import (Configuration, empty_key, preparation, arrival_projector,
                    complete_readout, norm2, clean, star_signature)
from .checks import local_preparation_checks, full_domain_checks, rotation_checks, experiment


def test_complete_preparation_preserves_gauge_spin_and_total_parity():
    result = local_preparation_checks()
    assert result['complete_preparation_current_commutators'] == 12
    assert result['complete_preparation_spin_commutators'] == 3
    assert result['total_fermion_parity'] == 'even'
    assert not result['D_is_canonical_fermion']


def test_full_occupation_operator_words():
    result = full_domain_checks()
    assert result['all_occupation_Gauss_identities'] == 24
    assert result['complete_normal_ordered_hopping_words_side_two'] == 768
    assert result['exact_word_Hermiticity']


def test_proper_rotations_include_actual_nonblank_full_H():
    result = rotation_checks()
    assert all(row['all_proper_rotations'] == 24 for row in result['complete_composite_word_rotations'])
    assert result['actual_full_H_nonblank_pi_rotation_controls'] == 3


def test_full_H_odd_arrival_has_a_positive_complete_record_bound():
    result = experiment()
    assert result['full_H_second_action_norm_squared'] == '15948017/128'
    assert result['first_action_terms_with_quark_occupations'] == 348
    assert result['complete_arrival_probability_interval'] == ['1/498436', '9/498436']
    assert result['complete_no_arrival_probability_interval'] == ['498427/498436', '498435/498436']
    assert not result['postselection_or_discarded_dynamical_sector']


def test_detector_keeps_all_flux_states_and_the_no_arrival_branch():
    psi = preparation()
    for label in range(119):
        state = {}
        for key, value in psi.items():
            links = list(key[1]); links[0] = label
            state[key[0], tuple(links), key[2]] = value
        arrived, missed = complete_readout(state, 0)
        assert norm2(arrived) + norm2(missed) == 1
        assert norm2(arrived) == int(label == 0)
        assert arrival_projector(arrived, 0) == arrived
        assert not arrival_projector(missed, 0)
        even_flux = {}
        for key, value in psi.items():
            links = list(key[1]); links[0] = label; links[1] = label
            even_flux[key[0], tuple(links), key[2]] = value
        assert arrival_projector(even_flux, 0) == even_flux
    # Basis inputs require rational half-coefficients, with no float injection.
    basis = {next(iter(psi)): 1}
    projected, complement = complete_readout(basis, 0)
    assert norm2(projected) == norm2(complement) == s.Rational(1, 2)
    assert all(not value.has(s.Float) for value in projected.values())
    # A colored spectator on the measured site is a retained non-arrival input.
    key = next(iter(psi))
    colored = {(key[0] | 3, key[1], key[2]): s.Integer(1)}
    assert star_signature(next(iter(colored)))[0] == -1
    assert not arrival_projector(colored, 0)
    assert complete_readout(colored, 0)[1] == colored


def test_exact_parameter_domain():
    with pytest.raises(ValueError):
        Configuration(neutral_hopping=0.1)
    with pytest.raises(ValueError):
        Configuration(neutral_hopping=s.I)
    assert Configuration(neutral_hopping=0).neutral_hopping == 0
    assert not arrival_projector({empty_key(): s.Integer(1)}, 0)


def test_custody_rejects_missing_declared_and_import_inputs(monkeypatch):
    from . import certificate as custody
    with monkeypatch.context() as patch:
        patch.setattr(custody.action, 'SOURCE_FILES', (*custody.action.SOURCE_FILES, 'missing-mobility-seed.json'))
        with pytest.raises(FileNotFoundError, match='missing declared'):
            custody.input_hashes()
    monkeypatch.setattr(custody, 'import_closure', lambda seeds: set(seeds) | {'missing-mobility-import.py'})
    with pytest.raises(FileNotFoundError, match='missing neutral-mobility import'):
        custody.input_hashes()
