"""Particle poles, exact index and native geometry retain their source domains."""
from copy import deepcopy
import pytest
from tools import strict_json
from . import defect_geometry_campaign as campaign


def packet():
    return strict_json.load(campaign.ROOT/campaign.PATH)


def test_five_issue_ownership_and_four_source_results():
    result = campaign.validate(packet())
    assert result['selected_issues'] == [1139, 1140, 1142, 1143, 1144]
    assert result['reviewed_components'] == 4
    assert result['canonical_premises_discharged'] == 0


@pytest.mark.parametrize('key', list(campaign.INTERFACES))
def test_each_source_boundary_is_load_bearing(key):
    value = deepcopy(packet())
    value['conditional_interfaces'][key] = not value['conditional_interfaces'][key]
    with pytest.raises(ValueError, match='domain promoted or erased'):
        campaign.validate(value)


@pytest.mark.parametrize('key', [
    'cross_parent_states_or_energy_levels_identified',
    'native_interacting_ground_preparation_executed',
    'native_repair_invariance_inferred_from_Hamiltonian_poles',
    'join_to_bounded_radial_index_or_responsive_Gram'])
def test_typed_finite_pole_join_rejects_extra_identifications(key):
    value = deepcopy(packet())
    value['supported_source_join'][key] = True
    with pytest.raises(ValueError, match='source join broadened or erased'):
        campaign.validate(value)


def test_missing_proof_artifact_is_rejected_even_with_receipt_retained():
    value = deepcopy(packet())
    evidence = value['component_results'][0]['evidence']
    del evidence['model/defect_core_response/PROOF.md']
    with pytest.raises(ValueError, match='incomplete component artifact custody'):
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
        value['component_results'][3]['interface']['state'] = 'Arbitrary correlated reused probe'
    else:
        value['integrated_gate_accepted'] = True
    with pytest.raises(ValueError):
        campaign.validate(value)
