"""Independent full-source checks of hole coherence, records and program choice."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
from scipy.linalg import expm

from model.defect_field_bridge import full_domain
from model.native_program import interpreter as native
from . import hole_records as h, program_selection as programs


def actual_matrices(sites,edge,law,masks):
    index={mask:i for i,mask in enumerate(masks)}
    out=[]
    for coin in range(2):
        for raw in range(full_domain.RAW_SLOTS):
            a=np.zeros((len(masks),len(masks)))
            for j,mask in enumerate(masks):
                code,target,sign,_=full_domain.branch_action(
                    mask,law,coin,sites=sites,flavors=2,seam=edge)
                if code==raw:a[index[target],j]=sign/np.sqrt(2)
            out.append(a)
    return out


@pytest.mark.parametrize('law',['C','E'])
def test_actual_full_spin_matrix_units_and_all_padding(law):
    sites=3; basis=dict(h.hole_basis(sites));masks=list(basis)
    edge=(0,2);matrices=actual_matrices(sites,edge,law,masks)
    assert len(matrices)==32
    assert np.linalg.norm(sum(a.T@a for a in matrices)-np.eye(len(masks)))<1e-14
    for i,x in enumerate(masks):
        for j,y in enumerate(masks):
            if basis[x]==basis[y]:continue
            unit=np.zeros((len(masks),len(masks)));unit[i,j]=1
            actual=sum(a@unit@a.T for a in matrices)
            assert np.linalg.norm(actual-float(h.cross_factor(basis[x],basis[y],edge,law))*unit)<1e-14


def test_scan_identifies_initial_position_without_resolving_internal_spin():
    sites=4;edges=((0,1),(2,3));basis=dict(h.hole_basis(sites))
    rows=h.scan(sites,edges)
    assert len(rows)==sites*2**len(edges)
    assert all(len(row)==2**(sites-1)for row in rows.values())
    for row in rows.values():
        assert len({basis[mask]for mask in row})==1
        assert set(sign for _,sign in row.values()) <= {-1,1}
    for x,y in product(basis,repeat=2):
        assert h.corrected_coefficient(sites,edges,x,y)==int(basis[x]==basis[y])
    assert h.scan_branch(sites,edges,((0,15),(1,15)))=={}


def test_corrected_complete_detector_preserves_arbitrary_reference_spin_blocks():
    sites=3;edges=((0,1),(1,2));basis=dict(h.hole_basis(sites));masks=list(basis)
    d=len(masks);reference=3;rng=np.random.default_rng(291)
    vector=rng.normal(size=(d,reference))+1j*rng.normal(size=(d,reference))
    vector/=np.linalg.norm(vector)
    original=np.outer(vector.ravel(),vector.ravel().conj())
    total=np.zeros_like(original)
    for row in h.scan(sites,edges).values():
        k=np.zeros((d,d));inverse=np.zeros((d,d))
        for start,(finish,sign)in row.items():
            i,j=masks.index(start),masks.index(finish)
            k[j,i]=sign/2**(len(edges)/2)
            inverse[i,j]=sign
        corrected=np.kron(inverse@k,np.eye(reference))
        total+=corrected@original@corrected.T
    expected=np.zeros_like(original)
    for hole in range(sites):
        p=np.kron(np.diag([int(basis[x]==hole)for x in masks]),np.eye(reference))
        expected+=p@original@p
    assert np.linalg.norm(total-expected)<1e-14
    assert abs(np.trace(total)-1)<1e-14
    # Explicit partial trace avoids relying on the detector's own labels.
    assert np.linalg.norm(sum(total[i*reference:(i+1)*reference,i*reference:(i+1)*reference]for i in range(d))
        -sum(original[i*reference:(i+1)*reference,i*reference:(i+1)*reference]for i in range(d)))<1e-14


@pytest.mark.parametrize('theta',[F(0),F(1,4),F(1)])
def test_complete_poisson_generator_matches_decay_and_diffusion(theta):
    sites=4;edges=((0,1),(1,2),(2,3),(0,3));gamma=F(2,3)
    # The fully aligned internal-spin block is invariant under every actual repair.
    masks=[sum(1<<(2*x)for x in range(sites)if x!=hole)for hole in range(sites)]
    generator=np.zeros((sites*sites,sites*sites));classical=np.zeros((sites,sites))
    for edge in edges:
        for law,weight in (('C',1-theta),('E',theta)):
            matrices=actual_matrices(sites,edge,law,masks)
            generator+=float(gamma*weight)*(sum(np.kron(a,a)for a in matrices)-np.eye(sites*sites))
        x,y=edge;classical[x,y]+=float(gamma/2);classical[y,x]+=float(gamma/2)
        classical[x,x]-=float(gamma/2);classical[y,y]-=float(gamma/2)
    time=.37;channel=expm(time*generator)
    for x,y in product(range(sites),repeat=2):
        unit=np.zeros((sites,sites));unit[x,y]=1
        actual=(channel@unit.ravel()).reshape(sites,sites)
        if x!=y:
            expected=np.exp(-time*float(h.coherence_rate(sites,edges,x,y,gamma,theta)))*unit
        else:
            expected=np.diag(expm(time*classical)[:,x])
        assert np.linalg.norm(actual-expected)<2e-14


def test_complete_program_capacity_does_not_select_coherent_action():
    row=programs.certificate();a,b=row['programs']
    assert a['correctable_record_count']==b['correctable_record_count']==10
    assert a['resources']==b['resources']
    assert (a['source_reference_off_diagonal'],b['source_reference_off_diagonal'])==('1/16','3/16')
    assert row['literal_default_c_weight']=='1/2'
    assert not row['fixed_program_predictions_are_undetermined']


def test_actual_uniform_attempts_match_finite_counting_theorem_specialization():
    sites=4;edges=((0,1),(1,2),(2,3),(0,3))
    masks=[sum(1<<(2*x)for x in range(sites)if x!=hole)for hole in range(sites)]
    average=sum(np.kron(a,a)for edge in edges for law in ('C','E')
                for a in actual_matrices(sites,edge,law,masks))/(2*len(edges))
    for y in (1,2):
        unit=np.zeros((sites,sites));unit[0,y]=1
        row=h.counting_word_multiplier(sites,edges,0,y,7)
        expected=float(row['complete_coherence_multiplier'])*unit
        assert np.linalg.norm((np.linalg.matrix_power(average,7)@unit.ravel()).reshape(sites,sites)-expected)<1e-14
    assert h.counting_word_multiplier(sites,edges,0,1,1)['one_attempt_decrement']==F(5,8)
    assert h.counting_word_multiplier(sites,edges,0,2,1)['one_attempt_decrement']==F(3,4)


@pytest.mark.parametrize('theta',[F(0),F(1,2),F(1)])
def test_complete_record_dilation_Gram_and_reference_stable_bound(theta):
    sites=3;edges=((0,1),(1,2));basis=dict(h.hole_basis(sites));masks=list(basis)
    gram=h.identity_dilation_gram(sites,edges,theta)
    dimension=len(masks);blocks=[];ideal=[]
    for edge in edges:
        for law,weight in(('C',1-theta),('E',theta)):
            for index,a in enumerate(actual_matrices(sites,edge,law,masks)):
                blocks.append(np.sqrt(float(weight)/len(edges))*a)
                ideal.append(np.sqrt(float(weight)/(2*len(edges)))*np.eye(dimension)
                             if index%full_domain.RAW_SLOTS==0 else np.zeros_like(a))
    v=np.concatenate(blocks);v0=np.concatenate(ideal)
    expected=np.diag([float(gram.get((mask,mask),0))for mask in masks])
    assert np.linalg.norm((v-v0).T@(v-v0)-expected)<1e-14
    assert np.linalg.norm(v.T@v-np.eye(dimension))<1e-14
    assert np.linalg.norm(v0.T@v0-np.eye(dimension))<1e-14
    # A nontrivial untouched reference uses the same full record-isometry bound.
    psi=np.arange(1,2*dimension+1,dtype=float).reshape(dimension,2);psi/=np.linalg.norm(psi)
    actual=v@psi;comparison=v0@psi
    overlap=np.vdot(actual,comparison)
    trace_distance=2*np.sqrt(max(0,1-abs(overlap)**2))
    assert trace_distance<=2*np.sqrt(np.max(np.diag(expected)))+1e-14


def test_counted_global_attempt_workaround_and_actual_Hubbard_scope():
    for n in(2,3,8,100):
        row=h.global_attempt_bound(n)
        assert row['global_attempts']==n*n
        assert F(row['one_step_dilation_squared_bound'])<=F(4,n**6)
        assert F(row['complete_source_archive_trace_bound'])==min(F(2),F(4,n))
    source=h.actual_aligned_hubbard_sector()
    assert [r['invariant_aligned_hole_dimension']for r in source['aligned_cases']]==[2,3,4]
    assert abs(source['general_one_hole_leakage']['T_over_c_half_matrix_element'])==1
    assert not source['aligned_sector_identified_as_interacting_ground']


def test_resource_budget_and_unknown_field_command_remain_checked():
    p=programs.selected_program(F(1,4));p['declared_resources']['attempts_upper']+=1
    with pytest.raises(ValueError):native.admit(p)
    p=programs.selected_program(F(1,4));p['commands'][0]={'op':'spin_ground','B':'1/32','repeat':1}
    with pytest.raises(ValueError):native.admit(p)


def test_reject_noncover_and_outside_coherence_domains():
    with pytest.raises(ValueError):h.scan(3,((0,1),))
    with pytest.raises(ValueError):h.coherence_rate(3,((0,1),),0,0)
    with pytest.raises(ValueError):h.coherence_rate(3,((0,1),),0,1,e_weight=2)
    with pytest.raises(ValueError):programs.selected_program(0)


def test_erasing_E_parity_is_detected_as_a_changed_source(monkeypatch):
    actual=h.branch_action
    def changed(*args,**kwargs):
        raw,out,sign,current=actual(*args,**kwargs)
        if args[1]=='E' and 4<=raw<8:raw-=4
        return raw,out,sign,current
    monkeypatch.setattr(h,'branch_action',changed)
    with pytest.raises(ArithmeticError):h.certificate()
