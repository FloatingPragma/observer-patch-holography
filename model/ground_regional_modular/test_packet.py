from copy import deepcopy

import pytest

from . import packet


@pytest.fixture(scope='module')
def value():
    return packet._produce()


def test_bound_source_and_exact_complete_ground_replay(value):
    packet.validate(value)
    assert value['source_bindings']==packet.source_bindings()
    assert 'model/spin_resolved_spectrum/results.json' in value['source_bindings']
    assert 'paper/tex_fragments/GROUND_REGIONAL_MODULAR.tex' in value['source_bindings']


@pytest.mark.parametrize('mutation',['kernel','faithfulness','gravity','partial'])
def test_zero_eigenspaces_and_claim_boundaries_cannot_be_promoted(value,mutation):
    changed=deepcopy(value);result=changed['regional_ground']
    if mutation=='kernel':result['all_site_pairs'][0]['kernel_rank']=0
    elif mutation=='faithfulness':result['faithful_on_full_regional_CAR']=True
    elif mutation=='gravity':result['geometric_cap_or_Einstein_identified']=True
    else:result['all_site_pairs'].pop()
    with pytest.raises(ValueError):packet.validate(changed)
