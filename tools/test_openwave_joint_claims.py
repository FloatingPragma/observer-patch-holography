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

@pytest.mark.parametrize('container', ['row', 'extension', 'observer'])
@pytest.mark.parametrize('parameter', ['mass', 'energy_mass', 'energy_coupling', 'time',
                                     'evolution_time', 'commit_probability'])
def test_nested_declared_parameters_are_exact(container, parameter):
    result={'supplied':{parameter:.23}, 'variants':[{'energy_mass':.23,'energy_coupling':.29,
                                                  'energy_accounting':{'total':.25}}]}
    receipt={'rows':[{'measurement':{'result':deepcopy(result)}}],
             'record_extensions':{'interacting_records':deepcopy(result)},
             'observer_protocol':deepcopy(result)}
    actual=deepcopy(receipt)
    target=(actual['rows'][0]['measurement']['result'] if container=='row' else
            actual['record_extensions']['interacting_records'] if container=='extension' else actual['observer_protocol'])
    target['supplied'][parameter] += 1e-13
    with pytest.raises(ValueError):compare(receipt, actual)

@pytest.mark.parametrize('parameter',['energy_mass','energy_coupling','time','commit_probability'])
def test_variant_parameter_outside_supplied_is_exact(parameter):
    receipt={'rows':[{'measurement':{'result':{'variants':[{parameter:.23}]}}}],
             'record_extensions':{'interacting_records':{'variants':[{parameter:.23}]}}}
    for container in ('row','extension'):
        actual=deepcopy(receipt)
        target=actual['rows'][0]['measurement']['result'] if container=='row' else actual['record_extensions']['interacting_records']
        target['variants'][0][parameter]+=1e-13
        with pytest.raises(ValueError):compare(receipt,actual)

def test_extension_readout_and_energy_roundoff_allowed():
    receipt={'record_extensions':{'interacting_records':{'variants':[
        {'energy_mass':.23,'energy_coupling':.29,
         'readout':{'position_probabilities':[.2,.8]},
         'energy_accounting':{'coherent_energy_change':1e-15}}]}}}
    actual=deepcopy(receipt)
    variant=actual['record_extensions']['interacting_records']['variants'][0]
    variant['readout']['position_probabilities'][0]+=1e-13
    variant['energy_accounting']['coherent_energy_change']=-1e-15
    compare(receipt,actual)

def test_nested_measured_metadata_types_remain_exact():
    receipt={'record_extensions':{'control':{'measurements':{'count':1,'accepted':True}}}}
    for field,value in (('count',True),('accepted',1)):
        actual=deepcopy(receipt);actual['record_extensions']['control']['measurements'][field]=value
        with pytest.raises(ValueError):compare(receipt,actual)
