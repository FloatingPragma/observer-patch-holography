import numpy as np
from . import decay as d

def test_full_current_form_factor_and_reference_shift():
    v=d.construction();p=v['initial'];r=v['target']
    assert abs(np.vdot(r,v['V']@p)-v['C'][0,1]/np.sqrt(2))<2e-12
    assert abs(np.vdot(p,v['V']@p)+3)<2e-12
    assert abs(np.vdot(r,v['V']@r)-2.5)<2e-12
    assert abs(np.vdot(p,v['H']@p)-np.vdot(r,v['H']@r))<1e-10

def test_complete_colour_spin_parity_reduction():
    v=d.construction();B,H,error=d.spin_parity_reduction(v)
    assert B.shape==(888,150)and error<5e-11
    assert np.linalg.norm(B.conj().T@B-np.eye(150))<1e-12
    assert np.max(abs(H-H.conj().T))<1e-12

def test_complete_original_records_and_readout_work():
    row=d.complete_records(d.construction(),720.)
    assert row['conditional_total_lepton_pair_probability']>.99
    assert row['unconditional_total_lepton_pair_probability']>.24
    assert abs(3*sum(row['complete_position_colour_triplet_probability_each'])+row['preparation_failure_probability']-1)<3e-12
    assert abs(3*sum(row['complete_momentum_colour_triplet_probability_each'])+row['preparation_failure_probability']-1)<3e-12
    assert abs(row['successful_position_readout_mean_work']-row['successful_momentum_readout_mean_work'])>1

def test_native_word_and_source_resonance_controls():
    zero=d.complete_records(d.construction(True),32.)
    assert zero['conditional_emission_by_quark_family'][0]<1e-22
    untuned=d.construction(False,False)
    detuning=np.vdot(untuned['target'],untuned['H']@untuned['target'])-np.vdot(untuned['initial'],untuned['H']@untuned['initial'])
    assert abs(detuning-11*d.G/2)<1e-10
