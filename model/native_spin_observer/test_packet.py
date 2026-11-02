from copy import deepcopy
from fractions import Fraction as F
import numpy as np
import pytest
from . import packet,experiment


@pytest.fixture(scope='module')
def retained():return packet.read_pair(packet.OUTPUT)


def test_complete_native_examples_and_reference(retained):
    value,arrays=retained
    assert len(value['execution'])==6
    assert value['execution']['C_repair_choi']['quantum_reference_dimension']==16
    assert arrays['C_crossing_reference_success'].shape==(2048,4)
    for name,row in value['execution'].items():
        assert row['complete_output_error_diagnostic']<float(F(value['native_error_certificates'][name]['complete_isometry_error']))
    p=value['execution']['observed_Neel_history']['logical_pointer_records']
    observed=sum(r['probability'] for r in p if r['spin_pointer']==1)
    assert abs(observed-np.sin(.5)**2)<1e-10
    for law in 'CE':
        p=value['execution'][law+'_repair_choi']['logical_pointer_records']
        assert abs(sum(r['probability'] for r in p if r['raw']==8)-7/16)<1e-10


@pytest.mark.parametrize('field',['infinite_ground_state_prepared','fixed_binary64_control_error_sent_to_zero',
                                'bare_repair_during_each_constituent_pulse_protected','physical_particle_species_identified'])
def test_scope_promotions_rejected(retained,field):
    value,arrays=retained;changed=deepcopy(value);changed['theorem'][field]=True
    with pytest.raises(ValueError):packet.validate(changed,arrays)


def test_type_and_custody_laundering_rejected(retained):
    value,arrays=retained
    for field,new in [('specified_finite_Neel_state_prepared',1),('infinite_ground_state_prepared',0)]:
        changed=deepcopy(value);changed['theorem'][field]=new
        with pytest.raises(ValueError):packet.validate(changed,arrays)
    changed=deepcopy(value);changed['prior_custody'][packet.PRIORS[0]]['old_expensive_execution_repeated']=True
    with pytest.raises(ValueError):packet.validate(changed,arrays)


def test_outcome_and_reference_tampering_rejected(retained):
    value,arrays=retained
    changed=dict(arrays);changed['C_repair_choi_success']=arrays['C_repair_choi_success'][:,:4]
    with pytest.raises(ValueError):packet.validate(value,changed)
    changed=dict(arrays);changed['C_repair_choi_failure_weights']=arrays['C_repair_choi_failure_weights'].copy()
    changed['C_repair_choi_failure_weights'][0,0]=.1
    with pytest.raises(ValueError):packet.validate(value,changed)
    altered=deepcopy(value);altered['execution']['C_repair_choi']['quantum_reference_dimension']=16.
    with pytest.raises(ValueError):packet.validate(altered,arrays)
