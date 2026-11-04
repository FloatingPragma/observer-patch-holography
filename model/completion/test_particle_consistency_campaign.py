"""Finite results cannot silently become physical or joint-limit acceptance."""
from copy import deepcopy

import pytest

from tools import strict_json
from .particle_consistency_campaign import BOUNDARY, ISSUES, PATH, ROOT, validate


def packet():
    return strict_json.load(ROOT/PATH)


def test_consistency_result_ownership_and_physical_scope():
    result = validate(packet())
    assert result['selected_issues'] == ISSUES
    assert result['reviewed_components'] == 3
    assert result['canonical_premises_discharged'] == 0
    assert result['physical_particle_map_complete'] is False


@pytest.mark.parametrize('key', BOUNDARY)
def test_incompatible_scientific_scope_is_rejected(key):
    p = deepcopy(packet())
    p['composition_boundary'][key] = not BOUNDARY[key]
    with pytest.raises(ValueError, match='composition claim'):
        validate(p)


@pytest.mark.parametrize('damage', [
    'stale_source', 'missing_source', 'reassigned_issue', 'altered_result',
    'missing_component', 'duplicate_component', 'physical_acceptance', 'canonical_discharge',
])
def test_source_and_acceptance_mutations_are_rejected(damage):
    p = deepcopy(packet())
    if damage == 'stale_source':
        p['source_bindings'][next(iter(p['source_bindings']))] = '0'*64
    elif damage == 'missing_source':
        del p['source_bindings']['model/pixel_mass_closure/results.json']
    elif damage == 'reassigned_issue':
        p['component_results'][0]['issues'] = [1147]
    elif damage == 'altered_result':
        p['component_results'][0]['result'] = 'Nature is derived.'
    elif damage == 'missing_component':
        p['component_results'].pop()
    elif damage == 'duplicate_component':
        p['component_results'].append(p['component_results'][0])
    elif damage == 'physical_acceptance':
        p['integrated_gate_accepted'] = True
    elif damage == 'canonical_discharge':
        p['canonical_premises_discharged'] = ['PR-23']
    with pytest.raises(ValueError):
        validate(p)
