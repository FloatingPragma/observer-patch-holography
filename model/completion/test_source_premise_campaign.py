"""Actual carrier incompatibilities and campaign promotion countercontrols."""
from copy import deepcopy
from fractions import Fraction as F

import pytest

from tools import strict_json
from . import source_premise_campaign as campaign


def packet():
    return strict_json.load(campaign.ROOT/campaign.PATH)


def test_campaign_has_exactly_five_issue_owners_without_closing_a_gate():
    result = campaign.validate(packet())
    assert result['selected_issues'] == [1138, 1139, 1142, 1143, 1145]
    assert result['reviewed_components'] == 6
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


def test_native_reference_and_original_one_particle_spectrum_have_different_algebras():
    native = strict_json.load(campaign.ROOT/'model/adaptive_sweep_reference/results.json')
    particle = strict_json.load(campaign.ROOT/'model/one_particle_search/results.json')
    assert native['full_carrier_mixing']['population_dimension'] == 264
    assert native['full_carrier_mixing']['stationary_reference_modular_automorphism'] == 'identity'
    assert particle['scope']['source_number_sector'] == 1
    assert particle['scope']['source_geometry'] == [8, 8, 8]
    assert not particle['scope']['localized_stationary_N1_core_constructed']
    assert not particle['scope']['many_body_defects_excluded']
    assert not particle['scope']['full_interacting_QFT_constructed']


def test_bipartite_estimate_does_not_certify_fixed_exchange_continuum():
    from model.bipartite_spin_limit import bounds
    fixed_exchange = bounds.refinement_exponents(2)
    controlled = bounds.refinement_exponents(F(11, 2))
    assert fixed_exchange['fixed_nonzero_exchange']
    assert not fixed_exchange['global_error_proved_to_vanish']
    assert controlled['global_error_proved_to_vanish']
    assert controlled['exchange_power'] < 0


def test_finite_modular_obstruction_does_not_remove_limiting_inclusion_route():
    value = packet()
    source = next(r for r in value['component_results'] if r['id'] == 'finite_modular_inclusion')
    assert 'No lower bound uniform over continuum scaling' in source['interface']['refinement']
    assert any('limiting representation' in text for text in source['remaining'])
    area = next(r for r in value['component_results'] if r['id'] == 'finite_area_response')
    assert area['interface']['state'] != source['interface']['state']
    assert not value['composition_boundary']['common_state_reference_constructed']
