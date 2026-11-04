"""Physical identification cannot be inferred from component progress."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import defect_selection_campaign as campaign


def packet():
    return strict_json.load(campaign.ROOT/campaign.PATH)


def test_four_results_cover_five_issues_without_accepting_physical_selection():
    result = campaign.validate(packet())
    assert result['reviewed_components'] == 4
    assert result['selected_issues'] == [1139, 1141, 1142, 1143, 1145]
    assert result['source_selection_solved'] is False


@pytest.mark.parametrize('key', list(campaign.BOUNDARY))
def test_each_mathematical_domain_is_retained(key):
    value = deepcopy(packet())
    value['composition_boundary'][key] = not value['composition_boundary'][key]
    with pytest.raises(ValueError, match='domain changed'):
        campaign.validate(value)


@pytest.mark.parametrize('damage', ['missing', 'clock', 'hash', 'count', 'gate'])
def test_missing_sources_and_scope_promotion_rejected(damage):
    value = deepcopy(packet())
    if damage == 'missing':
        value['component_results'].pop()
    elif damage == 'clock':
        value['component_results'][0]['interface']['clock'] = 'a measured physical clock'
    elif damage == 'hash':
        evidence = value['component_results'][0]['evidence']
        evidence[next(iter(evidence))] = '0'*64
    elif damage == 'count':
        value['physical_comparison_count_added'] = False
    else:
        value['integrated_gate_accepted'] = True
    with pytest.raises(ValueError):
        campaign.validate(value)
