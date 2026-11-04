"""Source-particle component progress cannot promote physical matching."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import particle_derivation_campaign as campaign


def packet():
    return strict_json.load(campaign.ROOT/campaign.PATH)


def test_six_results_cover_five_issues_with_a_bound_source_particle_map():
    result = campaign.validate(packet())
    assert result['reviewed_components'] == 6
    assert result['selected_issues'] == [1139, 1140, 1141, 1142, 1143]
    assert result['source_selection_solved'] is False
    assert result['physical_particle_map_complete'] is False
    assert result['source_map_bound'] == campaign.SOURCE_MAP_PATH


@pytest.mark.parametrize('key', list(campaign.BOUNDARY))
def test_each_finite_mathematical_domain_is_retained(key):
    value = deepcopy(packet())
    value['composition_boundary'][key] = not value['composition_boundary'][key]
    with pytest.raises(ValueError, match='domain changed'):
        campaign.validate(value)


@pytest.mark.parametrize('damage', ['missing', 'duplicate', 'clock', 'hash', 'proof_pin',
                                   'count', 'gate', 'discharge', 'issue', 'boundary_key',
                                   'boundary_type', 'ownership'])
def test_missing_sources_and_physical_scope_promotion_rejected(damage):
    value = deepcopy(packet())
    if damage == 'missing':
        value['component_results'].pop()
    elif damage == 'duplicate':
        value['component_results'][1] = deepcopy(value['component_results'][0])
    elif damage == 'clock':
        value['component_results'][0]['interface']['clock'] = 'a measured physical clock'
    elif damage == 'hash':
        evidence = value['component_results'][0]['evidence']
        evidence[next(iter(evidence))] = '0'*64
    elif damage == 'proof_pin':
        evidence = value['component_results'][0]['evidence']
        proof = next(name for name in evidence if name.endswith('/PROOF.md'))
        evidence.pop(proof)
    elif damage == 'count':
        value['physical_comparison_count_added'] = False
    elif damage == 'gate':
        value['integrated_gate_accepted'] = True
    elif damage == 'discharge':
        value['canonical_premises_discharged'] = ['claimed physical matching']
    elif damage == 'issue':
        value['selected_issues'][-1] = 1145
    elif damage == 'boundary_key':
        value['composition_boundary']['unreviewed_join'] = True
    elif damage == 'boundary_type':
        value['composition_boundary']['simulator_observer_scheduler_repair_memory_specs_admitted'] = 1
    else:
        value['component_results'][0]['issues'] = []
    with pytest.raises(ValueError):
        campaign.validate(value)


@pytest.mark.parametrize('damage', ['missing', 'path', 'hash', 'extra', 'type'])
def test_source_particle_map_is_exactly_bound(damage):
    value = deepcopy(packet())
    if damage == 'missing':
        value.pop('source_map')
    elif damage == 'path':
        value['source_map']['path'] = 'model/completion/README.md'
    elif damage == 'hash':
        value['source_map']['sha256'] = '0'*64
    elif damage == 'extra':
        value['source_map']['physical_matching_accepted'] = True
    else:
        value['source_map']['sha256'] = 0
    with pytest.raises(ValueError):
        campaign.validate(value)


def test_nonmapping_campaign_is_rejected():
    with pytest.raises(ValueError, match='campaign'):
        campaign.validate([])


@pytest.mark.parametrize('damage', ['missing', 'path', 'hash', 'extra'])
def test_downstream_comparison_custody_cannot_be_dropped(damage):
    value = deepcopy(packet())
    if damage == 'missing':
        value.pop('downstream_comparison')
    elif damage == 'path':
        value['downstream_comparison']['path'] = campaign.SOURCE_MAP_PATH
    elif damage == 'hash':
        value['downstream_comparison']['sha256'] = '0'*64
    else:
        value['downstream_comparison']['independent_prediction'] = True
    with pytest.raises(ValueError):
        campaign.validate(value)
