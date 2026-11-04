"""Conditional domains cannot be silently promoted to a physical composition."""
from copy import deepcopy
import pytest
from tools import strict_json
from . import spectrum_einstein_campaign as campaign


def packet():
    return strict_json.load(campaign.ROOT/campaign.PATH)


def test_exact_five_issue_ownership_and_three_scoped_results():
    result = campaign.validate(packet())
    assert result['selected_issues'] == [1139, 1142, 1143, 1145, 1146]
    assert result['reviewed_components'] == 3
    assert result['canonical_premises_discharged'] == 0


@pytest.mark.parametrize('key', [
    'supplied_texture_has_conditionally_identified_spectrum',
    'massless_bulk_and_gapped_defect_regimes_distinct',
    'all_eight_native_valleys_retained',
    'full_Fock_virtual_core_complement_retained',
    'texture_energy_is_actual_ground_orbit_result',
    'texture_energy_proves_stable_defect_particle',
    'ground_texture_self_consistently_supplies_Dirac_mass',
    'native_repair_protects_defect_spectrum',
    'same_source_regional_Einstein_implication_proved',
    'regional_Einstein_premises_instantiated',
    'finite_thermal_eta_is_universal_Newton_coefficient',
    'single_vacuum_tangent_proves_nonlinear_Einstein',
    'physical_particle_identification',
    'one_common_interacting_QFT_and_responsive_Einstein_source'])
def test_every_conditional_domain_is_load_bearing(key):
    value = deepcopy(packet())
    value['conditional_interfaces'][key] = not value['conditional_interfaces'][key]
    with pytest.raises(ValueError, match='domain promoted or erased'):
        campaign.validate(value)


@pytest.mark.parametrize('damage', ['missing', 'clock', 'reference', 'evidence', 'state', 'gate'])
def test_custody_and_common_source_join_fail_closed(damage):
    value = deepcopy(packet())
    if damage == 'missing':
        value['component_results'].pop()
    elif damage == 'clock':
        value['composition_boundary']['common_clock_and_refinement_constructed'] = True
    elif damage == 'reference':
        value['composition_boundary']['common_state_reference_constructed'] = True
    elif damage == 'evidence':
        name = next(iter(value['component_results'][0]['evidence']))
        value['component_results'][0]['evidence'][name] = '0'*64
    elif damage == 'state':
        value['component_results'][2]['interface']['state'] = 'Supplied scalar labels'
    else:
        value['integrated_gate_accepted'] = True
    with pytest.raises(ValueError):
        campaign.validate(value)
