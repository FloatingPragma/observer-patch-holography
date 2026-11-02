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


def test_near_zero_host_diagnostic_is_compared_absolutely():
    # Actual macOS/Linux recomputations from the same retained full array.
    packet.compare({'complete_output_error_diagnostic':2.2471471556841873e-15},
                   {'complete_output_error_diagnostic':2.2461137643882236e-15})
    packet.compare({'normalization_residual':1.7319479184152442e-14},
                   {'normalization_residual':1.6542323066914832e-14})
    with pytest.raises(ValueError):
        packet.compare({'normalization_residual':0.},
                       {'normalization_residual':3*packet.REPLAY_ATOL})
    with pytest.raises(ValueError,match='negative'):
        packet.compare({'normalization_residual':0.}, {'normalization_residual':-1e-30})


def test_failure_replay_uses_one_total_amplitude_budget():
    zero=np.zeros((25,4),dtype=np.float64)
    small=np.full_like(zero,(packet.REPLAY_ATOL/20)**2)
    packet.compare_arrays(zero,small,'example_failure_weights')
    # Each forged slot is individually below the amplitude tolerance; their
    # combined norm exceeds it. A per-slot allclose would incorrectly accept.
    distributed=np.full_like(zero,(packet.REPLAY_ATOL/5)**2)
    with pytest.raises(ValueError,match='total failure amplitude'):
        packet.compare_arrays(zero,distributed,'example_failure_weights')
    with pytest.raises(ValueError):
        packet.compare({'complete_first_failure_probability':0.},
                       {'complete_first_failure_probability':float(distributed.sum())})


def test_replay_tolerance_cannot_replace_source_certificate(retained):
    value,arrays=retained
    altered=deepcopy(value)
    altered['native_error_certificates']['exchange_choi']['complete_isometry_error']='1'
    with pytest.raises(ValueError,match='exact theorem'):
        packet.validate(altered,arrays)
    changed=dict(arrays);changed['exchange_choi_success']=arrays['exchange_choi_success'].copy()
    changed['exchange_choi_success'][0,0]+=1e-4
    with pytest.raises((ValueError,ArithmeticError)):
        packet.validate(value,changed)
