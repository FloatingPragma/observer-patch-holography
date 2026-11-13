import numpy as np
import pytest
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


def test_family_energy_offset_preserves_native_ground_projector():
    v=d.construction();native=d.sector('meson')
    active=[i for i,(_,anti)in enumerate(native['labels'])if anti[2]==0]
    rows=[i for i,label in enumerate(v['hadron_labels'])if label[0]==1 and label[3]==1]
    local=np.array(rows)%32
    h=native['h'][np.ix_(active,active)][np.ix_(local,local)]
    lookup={int(value):i for i,value in enumerate(v['keep'])}
    prepared=np.array([v['initial'][lookup[49*i]]if 49*i in lookup else 0 for i in rows])
    # NumPy's independent divide-and-conquer solve checks the native projector.
    eigenvalues,eigenvectors=np.linalg.eigh(h);ground=eigenvectors[:,0]
    assert eigenvalues[1]-eigenvalues[0]>1e-5
    assert np.linalg.norm(np.outer(prepared,prepared.conj())-np.outer(ground,ground.conj()))<2e-12
    shift=v['mass'][1]-1;shifted=h+shift*np.eye(len(h))
    rounding=64*np.finfo(float).eps*len(h)*max(1,np.linalg.norm(shifted,ord=2))
    assert abs(v['parent_energy']-(eigenvalues[0]+shift))<rounding
    assert np.linalg.norm(shifted@prepared-v['parent_energy']*prepared)<rounding
    shifted_values,shifted_frame=d.eigh(shifted)
    expected_seed=lookup[49*rows[int(np.argmax(abs(shifted_frame[:,0])**2))]]
    assert v['spectral_seed']==expected_seed
    assert np.array_equal(v['seed'],np.eye(len(v['keep']))[:,expected_seed])


@pytest.mark.parametrize('damage',('source','preparation'))
def test_zero_word_invariant_evolution_rejects_nonzero_complements(damage):
    v=dict(d.construction(True))
    family=np.array([v['hadron_labels'][int(x)//49][0]==1 for x in v['keep']])
    inside=np.flatnonzero(family)[0];outside=np.flatnonzero(~family)[0]
    if damage=='source':
        h=v['H'].tolil();h[inside,outside]=h[outside,inside]=1e-30;v['H']=h.tocsr()
    else:
        v['initial']=v['initial'].copy();v['initial'][outside]=1e-30
    with pytest.raises(AssertionError):
        d.complete_records(v,32.)


def test_shifted_eigensystem_rejects_a_different_native_source():
    native=np.array([[1.,.25],[.25,2.]])
    H=native+7*np.eye(2)
    changed=native.copy();changed[0,1]+=.01
    with pytest.raises(AssertionError):
        d._shifted_eigensystem(H,7,native=changed)
