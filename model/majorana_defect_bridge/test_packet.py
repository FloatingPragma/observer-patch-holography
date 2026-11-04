"""Cold receipt, full-output and inherited-source rejection controls."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def fresh():
    return packet.report()


def test_stored_receipt_matches_complete_cold_calculation(fresh):
    stored = strict_json.load(packet.OUTPUT)
    packet.validate(stored, fresh)
    assert stored['actual_action_spectrum']['full_K45_nullity'] == 1
    assert stored['native_current_Higgs_dictionary']['native_current_count'] == 12


@pytest.mark.parametrize('mutation', ['massless', 'torsion', 'mass', 'gate', 'source_selection', 'omit_car'])
def test_scientific_omissions_and_promotions_are_rejected(fresh, mutation):
    wrong = deepcopy(fresh)
    if mutation == 'massless': wrong['actual_action_spectrum']['neutrino_rows'][0]['massive_Majorana'] = True
    elif mutation == 'torsion': wrong['additive_topological_constraints']['all_three_positive_with_optional_instanton']['torsion_invariant_factors'] = []
    elif mutation == 'mass': wrong['actual_action_spectrum']['neutrino_rows'][1]['action_mass_exact'] = '1'
    elif mutation == 'gate': wrong['scope']['integrated_gate_accepted'] = True
    elif mutation == 'source_selection': wrong['scope']['source_selects_physical_couplings_or_Majorana_phases'] = True
    else: wrong['full_momentum_symbol']['physical_CAR_annihilation_modes'] = 6
    with pytest.raises(ValueError, match='complete fresh calculation'):
        packet.validate(wrong, fresh)


def test_stale_inherited_source_is_rejected_before_replay(monkeypatch):
    original = strict_json.load

    def corrupted(path):
        result = original(path)
        if str(path).endswith('native_matter_attachment_results.json'):
            result['source_bindings']['model/common_source/native_matter_attachment.py'] = '0'*64
        return result

    monkeypatch.setattr(strict_json, 'load', corrupted)
    with pytest.raises(ValueError, match='stale inherited source'):
        packet.source_bindings()
