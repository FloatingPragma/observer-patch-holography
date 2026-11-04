"""Actual carrier incompatibilities and campaign promotion countercontrols."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import dynamics_response_campaign as campaign


def packet():
    return strict_json.load(campaign.ROOT/campaign.PATH)


def test_campaign_has_exactly_five_issue_owners_without_closing_a_gate():
    result = campaign.validate(packet())
    assert result['selected_issues'] == [1138, 1142, 1143, 1144, 1145]
    assert result['reviewed_components'] == 4
    assert not result['integrated_gate_accepted']


@pytest.mark.parametrize('damage', ['missing_component', 'wrong_owner', 'canonical_closure',
                                  'physical_closure', 'joined_clock', 'joined_reference',
                                  'altered_state', 'altered_evidence', 'discarded_boundary'])
def test_component_promotion_and_silent_join_fail_closed(damage):
    value = deepcopy(packet())
    if damage == 'missing_component':
        value['component_results'].pop()
    elif damage == 'wrong_owner':
        value['component_results'][0]['issues'] = [1145]
    elif damage == 'canonical_closure':
        value['canonical_premises_discharged'] = ['PR-26']
    elif damage == 'physical_closure':
        value['common_interacting_QFT_and_Einstein_model'] = True
    elif damage == 'joined_clock':
        value['composition_boundary']['common_clock_and_refinement_constructed'] = True
    elif damage == 'joined_reference':
        value['composition_boundary']['common_state_reference_constructed'] = True
    elif damage == 'altered_state':
        value['component_results'][1]['interface']['state'] = value['component_results'][0]['interface']['state']
    elif damage == 'altered_evidence':
        source = next(iter(value['component_results'][0]['evidence']))
        value['component_results'][0]['evidence'][source] = '0'*64
    else:
        value['component_results'][1]['remaining'] = []
    with pytest.raises(ValueError):
        campaign.validate(value)


@pytest.mark.parametrize('field', ['arbitrary_correlated_reference_Gibbs_replacement',
                                 'finite_c_thermal_repair_stationarity',
                                 'infinite_volume_thermal_KMS_or_unique_phase_proved',
                                 'intervention_complete_multitime_stationarity_proved',
                                 'full_QFT_or_Einstein_composition'])
def test_supported_spin_join_keeps_its_reference_and_stationarity_domain(field):
    value = deepcopy(packet())
    value['shared_source_joins'][0][field] = True
    with pytest.raises(ValueError, match='supported spin join'):
        campaign.validate(value)


def test_supported_same_source_join_is_retained_without_postselection():
    value = packet()
    campaign.validate(value)
    join = value['shared_source_joins'][0]
    assert join['shared_unconditional_rotated_spin_pullback']
    assert join['finite_volume_thermal_comparison']
    assert join['ground_and_thermal_selection_are_supplied_inputs']
    assert not join['arbitrary_correlated_reference_Gibbs_replacement']
    assert 'No initial-dressing term' in join['thermal_error']
