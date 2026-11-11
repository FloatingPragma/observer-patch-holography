"""Separate finite preparations cannot silently become a joined physical law."""
from copy import deepcopy
import json

import pytest

from .finite_record_dynamics import PATH, ROOT, validate


def packet():
    return json.loads((ROOT / PATH).read_text())


def test_live_component_contract():
    assert validate(packet()) == {
        'finite_components': 3, 'common_source_join_established': False,
        'integrated_gate_accepted': False, 'physical_comparison_count_added': 0}


@pytest.mark.parametrize('mutation', ['join', 'gate', 'comparison', 'premise',
                                     'missing', 'source', 'clock', 'transfer', 'owner'])
def test_unproved_join_and_missing_inputs_rejected(mutation):
    changed = deepcopy(packet())
    if mutation == 'join':
        changed['common_source_join_established'] = True
    elif mutation == 'gate':
        changed['integrated_gate_accepted'] = True
    elif mutation == 'comparison':
        changed['physical_comparison_count_added'] = 1
    elif mutation == 'premise':
        changed['canonical_premises_discharged'] = ['PR-23']
    elif mutation == 'missing':
        changed['components'].pop()
    elif mutation == 'source':
        changed['components'][0]['sha256'] = '0' * 64
    elif mutation == 'clock':
        changed['components'][0]['interface'].pop('clock')
    elif mutation == 'transfer':
        changed['components'][0]['common_source_transfer_proved'] = True
    else:
        changed['components'][0]['supports'] = ['CP-10']
    with pytest.raises(ValueError):
        validate(changed)
