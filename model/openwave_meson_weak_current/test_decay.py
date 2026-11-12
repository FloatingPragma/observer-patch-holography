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

def test_complete_exact_family_blocks_across_lapack_drivers():
    v=d.construction(True);B,H,_=d.spin_parity_reduction(v)
    count,labels=d.sp.csgraph.connected_components(d.sp.csr_matrix(H!=0),directed=False)
    assert count==3 and np.array_equal(np.bincount(labels),[50,50,50])
    assert np.count_nonzero(H*(labels[:,None]!=labels[None,:]))==0
    family=np.array([v['hadron_labels'][int(x)//49][0]for x in v['keep']])
    for driver in('evr','evd','ev','evx'):
        e,U=d._complete_eigensystem(H,driver=driver)
        assert U.shape==(150,150)and len(e)==150
        assert np.linalg.norm(U.conj().T@U-np.eye(150))<2e-12
        assert np.linalg.norm(H@U-U*e)<2e-9
        assert np.max(abs(e-np.linalg.eigvalsh(H)))<2e-9
        initial=B.conj().T@v['initial']
        for time in(0.,32.,720.):
            state=B@U@(np.exp(-1j*time*e)*(U.conj().T@initial))
            assert np.count_nonzero(state[family!=1])==0
        assert d.complete_records(v,32.,(e,U))['conditional_emission_by_quark_family'][0]<1e-22

def test_component_partition_retains_arbitrarily_small_nonzero_edges(monkeypatch):
    H=np.diag([1.,2.,3.]).astype(complex);H[0,1]=H[1,0]=1e-30
    original=d.eigh;orders=[]
    def record(block,**kwargs):
        orders.append(len(block));return original(block,**kwargs)
    monkeypatch.setattr(d,'eigh',record)
    e,U=d._complete_eigensystem(H)
    assert sorted(orders)==[1,2]and U.shape==(3,3)and len(e)==3
    assert np.linalg.norm(H@U-U*e)<1e-14

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
