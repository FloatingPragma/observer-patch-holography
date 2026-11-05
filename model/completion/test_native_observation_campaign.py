"""Finite observation joins cannot substitute an input or promote physical scope."""
from copy import deepcopy

import pytest
from tools import strict_json
from .native_observation_campaign import ROOT, PATH, ISSUES, validate


def test_four_source_results_cover_five_issues():
    assert validate(strict_json.load(ROOT/PATH)) == {
        'selected_issues': ISSUES, 'finite_components': 4,
        'canonical_premises_discharged': 0, 'integrated_gate_accepted': False}


@pytest.mark.parametrize('mutation', ['reference', 'clock', 'refinement', 'owner',
                                     'missing', 'source', 'premise', 'acceptance'])
def test_substitution_and_physical_promotion_are_rejected(mutation):
    packet = deepcopy(strict_json.load(ROOT/PATH))
    if mutation in ('reference', 'clock', 'refinement'):
        key = 'state' if mutation == 'reference' else mutation
        packet['component_results'][0]['interface'][key] = 'arbitrary continuum physical source'
    elif mutation == 'owner':
        packet['component_results'][2]['issues'] = [1145]
    elif mutation == 'missing':
        packet['component_results'].pop()
    elif mutation == 'source':
        name = next(iter(packet['source_bindings']))
        packet['source_bindings'][name] = '0'*64
    elif mutation == 'premise':
        packet['canonical_premises_discharged'] = ['PR-23']
    else:
        packet['integrated_gate_accepted'] = True
    with pytest.raises(ValueError):
        validate(packet)
