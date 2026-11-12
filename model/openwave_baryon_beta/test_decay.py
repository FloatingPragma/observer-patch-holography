import numpy as np
from .certificate import calculate
from .decay import construction,eigensystem,evolve

def test_complete_native_channels_and_unconditional_records():
    r=calculate();assert r['construction']['complete_charge_zero_dimension']==1400
    assert r['construction']['exact_twice_spin_projection_one_dimension']==398
    assert r['calibrated_threshold']['electron_endpoint_GeV']>0
    assert r['calibrated_threshold']['double_electron_Delta_endpoint_GeV']<0
    row=r['records'][6];assert row['conditional_electron_count_probability'][1]>.95
    assert row['conditional_electron_count_probability'][2]>0
    assert abs(sum(row['original_electron_count_probability'])+row['original_failure_probability']-1)<1e-11
    for row in r['records']:
        for key in ('position_probability','momentum_probability'):assert abs(sum(row[key])-1)<1e-11

def test_full_failure_state_evolves_without_loss():
    d=construction();e=eigensystem(d);failure=evolve(d,64,e,d['failure'])
    assert abs(np.vdot(failure,failure).real+d['success']-1)<1e-11
    assert abs(np.vdot(failure,d['H']@failure).real-np.vdot(d['failure'],d['H']@d['failure']).real)<2e-8
