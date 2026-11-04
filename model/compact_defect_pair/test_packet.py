import copy
from fractions import Fraction as F

import pytest

from . import packet


@pytest.fixture(scope='module')
def replay():
    return packet.report()


def test_complete_cold_replay(replay):
    stored = packet.strict_json.load(packet.OUTPUT)
    packet.validate(stored, replay)
    assert len(replay['finite_spectra']) == 4
    assert replay['complete_exact_certificate']['exact_zero_mode_count'] == 0
    assert F(replay['infinite_volume']['core_abs_energy_certified_lower']) > F(1, 400)
    for row in replay['finite_spectra']:
        assert len(row['complete_numerical_spectrum']) == row['one_body_dimension']
        assert len(row['isolated_numerical_energies']) == 4
        assert row['minimum_core_probability'] > .99


def test_all_current_and_negative_controls(replay):
    value = replay['same_current_and_pairing_control']
    assert value['native_current_count'] == 12
    assert value['complex_kernel_full_Nambu_control_passed']
    assert value['zero_core_formula_fails']
    assert value['core_shift_preserves_shell_data_but_changes_spectrum']
    assert replay['canonical_physical_premises_discharged'] == []
    assert replay['physical_comparison_count_added'] == 0


@pytest.mark.parametrize('mutation', ['drop_level', 'fake_zero', 'physical_claim', 'source'])
def test_receipt_mutations_fail(replay, mutation):
    value = copy.deepcopy(replay)
    if mutation == 'drop_level':
        value['finite_spectra'][0]['complete_numerical_spectrum'].pop()
    elif mutation == 'fake_zero':
        value['complete_exact_certificate']['exact_zero_mode_count'] = 2
    elif mutation == 'physical_claim':
        value['scope']['observer_consensus_particle_identity_or_relativistic_dispersion'] = True
    else:
        value['source_bindings'].pop(next(iter(value['source_bindings'])))
    with pytest.raises(ValueError, match='differs'):
        packet.validate(value, replay)


def test_source_manifest_excludes_prose_notes_and_central_churn(replay):
    paths = replay['source_bindings']
    assert not any(p.startswith('docs/') for p in paths)
    assert 'model/verify.py' not in paths
    assert 'model/compact_defect_pair/certification.py' in paths
    assert 'model/common_source/chiral_candidate/action.py' in paths
