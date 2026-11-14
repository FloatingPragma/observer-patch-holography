"""Cross-runtime comparison rejects scientific, custody and authority mutations."""
from copy import deepcopy
import pytest
from tools.check_openwave_assessment_replay import compare

@pytest.fixture
def receipt():
    return {'versions':{'python':'3.13','numpy':'2.2'},'sources':{'source.py':'abc'},
            'rows':[{'id':1,'full_criterion_met':False,'measurement':{'error':1e-15,'value':.25}}]}

def test_runtime_and_roundoff_are_explicit(receipt):
    actual=deepcopy(receipt);actual['versions']['python']='3.11';actual['rows'][0]['measurement']['error']=2e-15
    compare(receipt,actual)

@pytest.mark.parametrize('mutation',['source','award','bool','nan','number','missing','runtime'])
def test_mutated_replay_rejected(receipt,mutation):
    actual=deepcopy(receipt)
    if mutation=='source':actual['sources']['source.py']='bad'
    elif mutation=='award':actual['rows'][0]['full_criterion_met']=True
    elif mutation=='bool':actual['rows'][0]['id']=True
    elif mutation=='nan':actual['rows'][0]['measurement']['value']=float('nan')
    elif mutation=='number':actual['rows'][0]['measurement']['value']=.3
    elif mutation=='missing':actual['rows']=[]
    else:actual['versions']={}
    with pytest.raises(ValueError):compare(receipt,actual)
