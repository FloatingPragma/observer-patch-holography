import numpy as np
from . import decay as d

def test_complete_parent_eigenpairs_and_symmetry_across_lapack_drivers():
    native=d.sector('meson')
    ids=[i for i,(q,a)in enumerate(native['labels'])if q[2]==1 and a[2]==0]
    labels=[native['labels'][i]for i in ids]
    shift=1/d.mass_data()[0]-1
    H=native['h'][np.ix_(ids,ids)]+shift*np.eye(16)
    columns=[]
    for xa in range(2):
        v=np.zeros(16)
        for xq,aq in((0,xa),(1,1-xa)):
            for sq,sa in((0,1),(1,0)):
                v[labels.index(((xq,sq,1),(aq,sa,0)))]=(1 if sq==0 else -1)/2
        columns.append(v)
    B=np.array(columns).T
    for driver in('evr','evd','ev','evx'):
        e,U=d._shifted_eigensystem(H,shift,driver=driver)
        assert len(e)==16 and np.linalg.norm(H@U-U*e)<2e-12
        assert np.max(abs(e-np.linalg.eigvalsh(H)))<2e-12
        assert abs(e[0]-(shift+(181-np.sqrt(5473))/48))<2e-12
        assert e[1]-e[0]>1/50
        assert np.linalg.norm(U[:,0]-B@(B.T@U[:,0]))<2e-12

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
