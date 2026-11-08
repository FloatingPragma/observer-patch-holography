"""Independent finite-instrument controls for the local record decoder."""
from collections import defaultdict
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
from math import comb
import copy
import json
import numpy as np
import pytest
from . import readout as r


def trace(ports, incident, edges, coins, mask):
    p=r.plan(ports,incident); state={i:(mask>>i)&1 for e in incident for i in e};events=[]
    for t,(edge,coin) in enumerate(zip(edges,coins,strict=True),1):
        a,b=edge;record=r.record(edge,coin,state[a]-state[b])
        if edge in incident:events.append((t,record))
        state[a]-=record.current;state[b]+=record.current
    return p,tuple(events),sum(state[i]<<k for k,i in enumerate(ports))


@pytest.mark.parametrize('size',[2,3,4,5])
def test_every_number_and_coin_word(size):
    ports=tuple(range(size));incident=tuple((i,i+1) for i in range(size-1))
    p=r.plan(ports,incident)
    for mask in range(1<<size):
        for coins in product((0,1),repeat=size-1):
            _,events,last=trace(ports,incident,p['tree'],coins,mask)
            result=r.decode(p,events,size-1)
            if mask in (0,(1<<size)-1):assert result['status']=='CONSTANT_ERASURE'
            else:assert (result['initial_mask'],result['terminal_mask'])==(mask,last)


def test_boundary_currents_and_first_hits():
    ports=(1,2,3);incident=((0,1),(1,2),(2,3),(3,4))
    words=(((0,1),(3,4),(1,2),(0,1),(2,3)),((2,3),(3,4),(1,2),(2,3)))
    for edges in words:
        for mask in range(32):
            for coins in product((0,1),repeat=len(edges)):
                p,events,last=trace(ports,incident,edges,coins,mask)
                result=r.decode(p,events,len(edges));initial=(mask>>1)&7
                if initial in (0,7):assert result['status']=='CONSTANT_ERASURE'
                else:assert (result['initial_mask'],result['terminal_mask'])==(initial,last)


def test_missing_boundary_record_can_change_the_inferred_initial_mask():
    ports=(1,2);incident=((0,1),(1,2));p,events,_=trace(ports,incident,incident,(0,1),1)
    assert r.decode(p,events,2)['status']=='CONSTANT_ERASURE'
    forged=r.decode(p,events[1:],2)
    assert forged['status']=='DECODED' and forged['initial_mask']==1
    # Source attachment, rather than the decoder alone, certifies completeness.
    actual=r.source_plan(0,tuple(range(12)));bad=r.plan(actual['ports'],tuple(e for e in actual['incident_seams'] if e[1]<12))
    with pytest.raises(ValueError,match='incident-edge'):
        r.raw_query(r.fixed_program(0,1,actual),bad,(0,)*11,{(1,0,1,0):F(1)})


def uniform_blocks(size,steps):
    edges=tuple((i,i+1) for i in range(size-1));p=r.plan(tuple(range(size)),edges)
    blocks={}
    for word in product(edges,repeat=steps):
        for coins in product((0,1),repeat=steps):
            for initial in range(1<<size):
                _,events,last=trace(p['ports'],edges,word,coins,initial)
                key=tuple(record for _,record in events)
                blocks.setdefault(key,[]).append((initial,last))
    # Each full raw label word fixes both selected seams and coins.
    return p,blocks,F(1,(2*len(edges))**steps)


def test_full_quantum_reference_oracle_distance_and_backaction():
    p,blocks,weight=uniform_blocks(3,3);dim=8;ref=3
    rng=np.random.default_rng(56);psi=rng.normal(size=(dim,ref))+1j*rng.normal(size=(dim,ref));psi/=np.linalg.norm(psi)
    rho=np.outer(psi.ravel(),psi.ravel().conj());total=np.zeros_like(rho)
    gamma=F();quantum_distance=0.;const=0.;coherence_survives=False
    for labels,pairs in blocks.items():
        K=np.zeros((dim,dim));
        for initial,last in pairs:K[last,initial]=1
        K*=float(weight)**.5;KR=np.kron(K,np.eye(ref));actual=KR@rho@KR.conj().T;total+=actual
        result=r.decode(p,tuple(enumerate(labels,1)),3)
        oracle=[]
        for initial in range(dim):
            Q=np.zeros((dim,dim));Q[initial,initial]=1
            L=KR@np.kron(Q,np.eye(ref));oracle.append(L@rho@L.conj().T)
        if result['status']=='DECODED':
            assert np.allclose(actual,oracle[result['initial_mask']],atol=1e-14)
            assert sum(np.trace(x).real for j,x in enumerate(oracle) if j!=result['initial_mask'])<1e-14
        else:
            # Distinct classical labels make both positive outputs orthogonal.
            quantum_distance+=(np.trace(actual).real+sum(np.trace(x).real for x in oracle))/2
            if result['status']=='CONSTANT_ERASURE':
                const+=np.trace(actual).real
                coherence_survives |= np.linalg.norm(actual[:ref,-ref:])>1e-8
    eta=np.trace(rho[:ref,:ref]).real+np.trace(rho[-ref:,-ref:]).real
    pcov=r.coverage(2,2,3,exact=True)['exact'];expected=float(pcov)+(1-float(pcov))*eta
    assert abs(np.trace(total).real-1)<1e-13
    assert abs(quantum_distance-expected)<1e-13
    assert abs(const-(1-float(pcov))*eta)<1e-13
    assert coherence_survives


def test_record_only_phase_blindness_is_not_state_identity():
    p,blocks,w=uniform_blocks(3,2)
    plus=np.zeros(8);minus=np.zeros(8);plus[[1,2]]=2**-.5;minus[[1,2]]=[2**-.5,-2**-.5]
    a=np.outer(plus,plus);b=np.outer(minus,minus)
    assert not np.array_equal(a,b)
    for pairs in blocks.values():
        K=np.zeros((8,8))
        for x,y in pairs:K[y,x]=1
        assert np.trace(K@a@K.T)==pytest.approx(np.trace(K@b@K.T))


@pytest.mark.parametrize('M,k,T',[(2,2,3),(3,2,4),(3,3,6),(5,3,8),(1,1,0),(1,1,1)])
def test_coupon_law_against_literal_selector_words(M,k,T):
    result=r.coverage(M,k,T,exact=True)
    failed=sum(not set(range(k))<=set(word) for word in product(range(M),repeat=T))
    assert result['exact']==F(failed,M**T)


@pytest.mark.parametrize('V',[3,8,17,240])
def test_canonical_erasure_all_tested_sectors(V):
    for q in range(V+1):
        for b in (2,min(3,V),min(12,V)):
            result=r.canonical_erasure(V,q,b)
            numerator=(comb(V-b,q) if q<=V-b else 0)+(comb(V-b,q-b) if q>=b else 0)
            assert result['exact']==F(numerator,comb(V,q))
            assert result['exact']<=result['product_upper']


def test_original_continuation_on_local_tree_all_reference_matrix_units():
    p=r.source_plan(0,tuple(range(12)));value=r.fixed_program(0,3,p)
    initial=(1<<0)|(1<<1)|(1<<20);other=(1<<0)|(1<<1)|(1<<21)
    # Same initial local configuration, arbitrary exterior/reference coherence.
    density={(initial,0,other,1):F(2,3),(other,1,initial,0):F(2,3)}
    state=initial;slots=[]
    for edge in p['tree']:
        coin=len(slots)%2;rec=r.record(edge,coin,((state>>edge[0])&1)-((state>>edge[1])&1),'fixed')
        rows=r.source._outcomes(value,value['commands'][len(slots)],{})
        slot=next(i for i,x in enumerate(rows) if r.inherited._same_record(x.record,rec));slots.append(slot)
        state+=(1<<edge[1])*rec.current-(1<<edge[0])*rec.current
    out=r.raw_query(value,p,tuple(slots),density)
    assert out['decoded']['initial_mask']==3
    assert out['density']=={(state,0,state-(1<<20)+(1<<21),1):F(2,3*2**11),
                            (state-(1<<20)+(1<<21),1,state,0):F(2,3*2**11)}
    assert len(out['archive'])==10 and len(out['memory'])==1
    zero=r.raw_query(value,p,(0,)*11,density)
    assert zero['density']=={} and len(zero['archive'])==10


def test_original_one_step_uniform_complete_alphabet():
    p=r.source_plan(0,tuple(range(12)));v=r.uniform_program(0,1,1);rho={(1,0,1,0):F(1)}
    outcomes=r.source._outcomes(v,v['commands'][0],rho)
    assert len(outcomes)==1+24*630
    assert sum(sum(x for (a,i,b,j),x in o.density.items() if a==b and i==j) for o in outcomes)==1
    incident=next(i for i,o in enumerate(outcomes) if o.density and o.record.seam in p['incident_seams'])
    query=r.raw_query(v,p,(incident,),rho)
    assert query['density']==outcomes[incident].density
    assert query['memory']==outcomes[incident].memory


def test_growing_admission_fixed_intercell_and_host_boundaries():
    p=r.source_plan(0,tuple(range(12)));fixed=r.fixed_program(200,80,p)
    assert fixed['declared_resources']['attempts_upper']==11
    stopped=r.source.execute(fixed,{'sentinel':True})
    assert stopped.status=='RESOURCE_STOP_UNEVALUATED_CARRIER'
    w=r.routing.routed_window(4,2);inner=r.routing.routed_window(4,1)
    indices=tuple(i for i,a in enumerate(w['source_ports']) if a in set(inner['source_ports']))
    multicell=r.rooted_plan(w,indices)
    assert len(multicell['ports'])==48
    with pytest.raises(ValueError,match='intercell'):r.fixed_program(4,80,multicell)
    ceiling=r.uniform_ceiling(200,80)
    uniform=r.uniform_program(200,80,ceiling)
    assert uniform['declared_resources']['attempts_upper']==ceiling
    with pytest.raises(ValueError):r.uniform_program(200,80,ceiling+1)
    with pytest.raises(ValueError):r.uniform_program(0,1,2**32+1)
    assert r.coverage(630*4**200,11,2**32)['lower']>F(999,1000)
    with pytest.raises(ValueError,match='boundary'):r.rooted_plan(w,tuple(range(len(w['source_ports']))))


@pytest.mark.parametrize('mutation',[
    lambda p: p.update(tree=tuple(reversed(p['tree']))),
    lambda p: p.update(ports=(False,1,2)),
    lambda p: p.update(extra=True),
    lambda p: p.update(incident_seams=([0,1],(1,2))),
])
def test_plan_rejects_mutations(mutation):
    p=r.plan((0,1,2),((0,1),(1,2)));mutation(p)
    with pytest.raises(ValueError):r.authenticate(p)


@pytest.mark.parametrize('mutation',[
    lambda rec:replace(rec,coin=True),lambda rec:replace(rec,current=1-rec.current),
    lambda rec:replace(rec,seam=(False,1)),lambda rec:replace(rec,selector='fixed'),
    lambda rec:replace(rec,raw=float(rec.raw)),lambda rec:replace(rec,label=float(rec.label)),
])
def test_literal_record_rejects_typed_mutations(mutation):
    p=r.plan((0,1),((0,1),));rec=mutation(r.record((0,1),0,1))
    with pytest.raises(ValueError):r.decode(p,((1,rec),),1)


@pytest.mark.parametrize('args',[(True,2,3),(3,True,3),(3,2,True),(3,4,2),(0,1,1)])
def test_coverage_typed_domains(args):
    with pytest.raises(ValueError):r.coverage(*args)


def test_zero_slots_kept_and_input_protocol_guard():
    p=r.plan((0,1),((0,1),));rec=r.source.Record('repair','uniform',(0,1),'E',0,7,None,None)
    # Find an actual unused address, without guessing its numeric value.
    rec=replace(rec,raw=next(k for k in range(8) if k not in r.source.circuit.record_values('E',0)))
    assert r.decode(p,((1,rec),),1)['status']=='ZERO_RAW_WORD'
    p=r.source_plan(0,tuple(range(12)));v=r.uniform_program(0,1,1)
    with pytest.raises(ValueError):r.raw_query(v,p,(True,),{})
    with pytest.raises(ValueError):r.raw_query(v,p,[1],{})
    with pytest.raises(ValueError):r.raw_query(v,p,(1,),{(1,0,1,0):True})
    hold=r.source.program(federation={'level':0,'total':1},biased=False,commands=[{'op':'hold','repeat':1}])
    with pytest.raises(ValueError):r.raw_query(hold,p,(0,),{})


def test_bound_protocol_and_typed_tree_cannot_be_relabelled():
    p=r.source_plan(0,tuple(range(12)));v=r.fixed_program(0,1,p)
    different=r.source.program(federation={'level':0,'total':1},biased=False,
        commands=[{'op':'repair','selector':{'kind':'fixed','seam':list(p['tree'][0])},'law':'C','repeat':11}],
        capacity=1,overflow='export_oldest')
    with pytest.raises(ValueError,match='protocol'):r.raw_query(different,p,(0,)*11,{})
    altered=copy.deepcopy(p);altered['tree']=list(altered['tree'])
    with pytest.raises(ValueError,match='typed'):r.authenticate(altered)


@pytest.fixture(scope='module')
def packet_value():
    from . import packet
    return packet.produce()


@pytest.mark.parametrize('mutation',[
    lambda p:p['scope'].update(all_local_number_sectors=1),
    lambda p:p['actual_source_kernel_controls'].update(fixed_attempts=11.0),
    lambda p:p['finite_routed_uniform_readout'].update(coverage_upper='0'),
    lambda p:p['source_bindings'].update({'model/native_program/interpreter.py':'0'*64}),
    lambda p:p.update(extra=False),
])
def test_receipt_rejects_typed_or_custody_mutations(tmp_path,packet_value,monkeypatch,mutation):
    from . import packet
    candidate=copy.deepcopy(packet_value);mutation(candidate)
    path=tmp_path/'receipt.json';path.write_text(packet.canonical(candidate))
    monkeypatch.setattr(packet,'produce',lambda:packet_value)
    with pytest.raises(ValueError):packet.verify(path)


def test_receipt_duplicate_keys_rejected(tmp_path,packet_value,monkeypatch):
    from . import packet
    text=packet.canonical(packet_value).replace('"schema":','"schema":"duplicate", "schema":',1)
    path=tmp_path/'receipt.json';path.write_text(text)
    monkeypatch.setattr(packet,'produce',lambda:packet_value)
    with pytest.raises(ValueError):packet.verify(path)


def test_source_binding_closure(packet_value):
    bindings=packet_value['source_bindings']
    required=('model/native_program/interpreter.py','model/observer/source_discrimination/circuit.py',
        'model/observer/source_scalar/compiler.py','model/native_rooted_refinement/routing.py',
        'model/native_sector_readout/readout.py','model/native_fair_local_readout/PROOF.md')
    assert all(p in bindings for p in required)
    assert not any(p.startswith('docs/') for p in bindings)


def test_global_total_resolves_only_feasible_constant_inputs():
    ports=(0,1,2);incident=((0,1),(1,2));p=r.plan(ports,incident)
    for mask in range(64):
        total=mask.bit_count()
        for coins in product((0,1),repeat=2):
            _,events,terminal=trace(ports,incident,incident,coins,mask)
            decoded=r.decode(p,events,2,population=6,total=total)
            both=r.feasibility(6,total,3)['ambiguous_constant'];initial=mask&7
            if both and initial in (0,7):assert decoded['status']=='CONSTANT_ERASURE'
            else:assert (decoded['initial_mask'],decoded['terminal_mask'])==(initial,terminal)
    assert r.decode(p,(),0,population=6,total=0)['initial_mask']==0
    assert r.decode(p,(),0,population=6,total=6)['initial_mask']==7
    assert r.decode(p,(),0,population=6,total=1)['status']=='COVERAGE_ERASURE'


def test_global_total_resolves_constant_with_boundary_backaction():
    p,events,last=trace((1,2),((0,1),(1,2)),((0,1),(1,2)),(0,1),1)
    out=r.decode(p,events,2,population=3,total=1)
    assert out['initial_mask']==0 and out['terminal_mask']==last==1


@pytest.mark.parametrize('V,b',[(6,3),(8,5),(17,12),(240,12)])
def test_tighter_canonical_erasure(V,b):
    for q in range(V+1):
        informed=r.canonical_erasure(V,q,b,use_global_total=True)
        old=r.canonical_erasure(V,q,b)
        feasible=r.feasibility(V,q,b)
        assert informed==(old if feasible['ambiguous_constant'] else {'exact':F(),'product_upper':F()})
        assert informed['exact']<=old['exact']
        assert r.error_bound(F(1,3),informed['exact'],global_constant=feasible['global_constant'])==(F() if q in (0,V) else F(1,3)+F(2,3)*informed['exact'])


def test_global_total_bound_to_literal_source_query():
    p=r.source_plan(0,tuple(range(12)))
    for q,mask,wanted in ((1,1<<20,0),(239,((1<<240)-1)^(1<<20),4095)):
        value=r.fixed_program(0,q,p);slots=[]
        for command in value['commands']:
            rows=r.source._outcomes(value,command,{})
            rec=r.record(tuple(command['selector']['seam']),len(slots)%2,0,'fixed')
            slots.append(next(i for i,row in enumerate(rows) if r.inherited._same_record(row.record,rec)))
        density={(mask,0,mask,0):F(1)}
        uninformed=r.raw_query(value,p,tuple(slots),density)
        informed=r.raw_query(value,p,tuple(slots),density,use_global_total=True)
        assert uninformed['decoded']['status']=='CONSTANT_ERASURE'
        assert informed['decoded']['initial_mask']==wanted
        assert uninformed['density']==informed['density']
        assert uninformed['memory']==informed['memory'] and uninformed['archive']==informed['archive']


@pytest.mark.parametrize('population,total',[(6,None),(None,1),(True,1),(6,True),(6,7),(2,1)])
def test_invalid_global_specification_rejected(population,total):
    p=r.plan((0,1,2),((0,1),(1,2)))
    with pytest.raises(ValueError):r.decode(p,(),0,population=population,total=total)


@pytest.mark.parametrize('total',(0,240))
def test_infeasible_global_sector_raw_word_is_retained(total):
    p=r.source_plan(0,tuple(range(12)));v=r.uniform_program(0,total,1)
    mask=0 if total==0 else (1<<240)-1;density={(mask,0,mask,0):F(1)}
    rows=r.source._outcomes(v,v['commands'][0],density)
    slot=next(i for i,out in enumerate(rows) if out.record.kind=='repair'
              and out.record.seam in p['incident_seams'] and out.record.law=='E'
              and out.record.label==1)
    assert rows[slot].density=={}
    result=r.raw_query(v,p,(slot,),density,use_global_total=True)
    assert result['density']=={} and result['decoded']['status']=='ZERO_RAW_WORD'
    assert result['memory']==rows[slot].memory and result['archive']==rows[slot].archive
    assert len(result['memory'])==1


def test_valid_contrast_word_with_no_global_candidate_is_zero():
    p,events,_=trace((1,2),((0,1),(1,2)),((0,1),(1,2)),(0,1),1)
    assert r.decode(p,events,2,population=3,total=2)['status']=='ZERO_RAW_WORD'
    p,events,_=trace((0,1,2),((0,1),(1,2)),((0,1),(1,2)),(0,1),0)
    assert r.decode(p,events,2,population=3,total=1)['status']=='ZERO_RAW_WORD'


@pytest.mark.parametrize('use_global_total',(False,True))
def test_valid_tree_word_with_no_binary_candidate_retains_every_raw_record(use_global_total):
    p=r.source_plan(0,tuple(range(12)));value=r.fixed_program(0,2,p)
    records=[];slots=[]
    for command in value['commands']:
        rec=r.record(tuple(command['selector']['seam']),1,1,'fixed')
        rows=r.source._outcomes(value,command,{})
        slots.append(next(i for i,row in enumerate(rows) if r.inherited._same_record(row.record,rec)))
        records.append(rec)
    assert slots==[21]*11
    # The (0,1),(1,8) contrasts require x_0-x_8=2; no binary input exists.
    assert records[0].label==records[5].label==1
    assert all(rec.current==0 for rec in records)
    result=r.raw_query(value,p,tuple(slots),{(3,0,3,0):F(1)},use_global_total=use_global_total)
    assert result['density']=={} and result['decoded']['status']=='ZERO_RAW_WORD'
    assert r.inherited._same_records(result['archive']+result['memory'],tuple(records))
    assert len(result['archive'])==10 and len(result['memory'])==1
