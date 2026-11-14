import numpy as np
import pytest
from link_phase_probe import execute,ExternalPhaseProbe,JointNativeSource


def test_same_source_probe_current_gauge_and_all_native_ledgers():
    r,states=execute();assert all(r['checks'].values())
    for law,x in r['native'].items():
        assert abs(sum(row['probability'] for row in x['complete_outcomes'])-1)<1e-12
        assert abs(sum(row['unnormalized_energy'] for row in x['complete_outcomes'])-x['energy_after'])<1e-12
        assert abs(sum(row['unnormalized_regional_charge'] for row in x['complete_outcomes'])-x['charge_after'])<1e-12
        assert abs(x['energy_after']-x['energy_before']-x['mean_native_work'])<1e-12
        assert abs(np.trace(states['native_'+law])-1)<1e-12
    assert not r['full_criterion_met']


def test_Kraus_dressing_phase_cancels_in_classical_record_maps_and_inputs_reject():
    model=JointNativeSource();probe=ExternalPhaseProbe(model);rho,_,_=model.reference()
    a=np.array([0,0,0,0,.17,0.])
    sub,ledger=probe.native_step(.4*rho,a)
    assert abs(sum(x['probability'] for x in ledger['complete_outcomes'])-.4)<1e-12
    assert abs(np.trace(sub)-.4)<1e-12
    for left,right in zip(probe.native_operators(np.zeros(6)),probe.native_operators(a)):
        for (_,k),(_,ka) in zip(left,right):assert np.max(abs(k@rho@k.T-ka@rho@ka.conj().T))<1e-12
    with pytest.raises(ValueError):probe.native_step(rho,a,probability=-.1)
    with pytest.raises(ValueError):probe.native_step(rho,a,probability=float('nan'))
    for value in ([0]*5,[0,0,0,0,np.nan,0],np.zeros(6,dtype=complex)):
        with pytest.raises(ValueError):probe.hamiltonian(value)
