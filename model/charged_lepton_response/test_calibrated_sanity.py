from copy import deepcopy
from fractions import Fraction as F
import json

import pytest

from . import calibrated_sanity as cs


def reference():
    return json.loads(cs.REFERENCE.read_text())


def test_chosen_phase_matches_mu_and_holds_out_tau():
    r = reference(); phase = cs.choose_phase(r['muon_electron_ratio']['value'])
    assert phase['tau_used_to_choose_phase'] is False
    result = cs.compare(r, phase)
    assert result['calibrated_mu_ratio_inside_reported_one_sigma_window']
    assert result['held_out_tau_ratio_inside_reported_one_sigma_window']
    assert result['comparison_proves_physical_lepton_identity'] is False
    assert abs(result['tau_MeV_approx']-1776.9690272267634) < 1e-9


def test_mutating_tau_changes_only_comparison():
    r = reference(); phase = cs.choose_phase(r['muon_electron_ratio']['value'])
    original = cs.compare(r, phase)
    changed = deepcopy(r); changed['tau_electron_ratio']['value'] = '4000'
    altered = cs.compare(changed, phase)
    assert not altered['held_out_tau_ratio_inside_reported_one_sigma_window']
    assert altered['source_action_ratios'] == original['source_action_ratios']
    assert cs.choose_phase(changed['muon_electron_ratio']['value']) == phase


def test_common_energy_calibration_leaves_ratios_unchanged():
    r = reference(); phase = cs.choose_phase(r['muon_electron_ratio']['value'])
    before = cs.compare(r, phase)
    r['electron_MeV']['value'] = str(7*F(r['electron_MeV']['value']))
    after = cs.compare(r, phase)
    assert before['source_action_ratios'] == after['source_action_ratios']
    assert [F(x)*7 for x in before['optional_calibrated_mass_MeV'][2]] == list(map(F, after['optional_calibrated_mass_MeV'][2]))


@pytest.mark.parametrize('r', [True, 1, 0, -1, 2.5])
def test_invalid_calibration_inputs(r):
    with pytest.raises(ValueError): cs.choose_phase(r)
