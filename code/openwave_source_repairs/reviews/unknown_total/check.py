#!/usr/bin/env python3
"""Retrospective independent source-record and theorem audit; no producer imports."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,deque
import ast,gzip,hashlib,json,sys,time
import numpy as np
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
P=HERE.parents[1]/'work/2026-11-16_codex_unknown_total_extraction'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def advance(q,event):
    a,b,c,j=event['record'];left,right=q[a],q[b]
    if event['law']=='native':
        low=(left+right)//2;high=left+right-low
        pair=(low,high)if c==0 else(high,low)
    else:
        assert event['law']=='swap_hold';pair=(left,right)if c==0 else(right,left)
    assert left-pair[0]==j and pair[1]-right==j
    q[a],q[b]=pair

def reconstruct(census):
    cell=census['cell'];nodes=list(range(12*cell,12*cell+12));delta={v:0 for v in nodes}
    equations=[];tape=census['tape'];cursor=0
    while cursor<len(tape):
        a,b,c,j=tape[cursor];before=delta[a]-delta[b]
        firstj=j;firstc=c;delta[a]-=j;delta[b]+=j;cursor+=1
        while cursor<len(tape):
            x,y,c,j=tape[cursor];assert (x,y)==(a,b)
            delta[a]-=j;delta[b]+=j;cursor+=1
            if c!=firstc:
                equations.append((a,b,2*firstj+j-before));break
        else:raise AssertionError('Census failed to resolve its parity')
    # Solve differences anchored at the declared INITIAL vacuum, not a total.
    known={(0,14)[cell]:0}
    while len(known)<12:
        previous=len(known)
        for a,b,d in equations:
            if a in known:known[b]=known[a]-d
            if b in known:known[a]=known[b]+d
        assert len(known)>previous
    initial=[known[v]for v in nodes];terminal=[known[v]+delta[v]for v in nodes]
    assert census['decoded']['initial']==initial and census['decoded']['terminal']==terminal
    assert census['decoded']['total']==sum(initial)
    return initial

def matrix_pump(q,seed):
    e=np.eye(12)-q.conj().T@q
    superop=np.kron(q,q.conj())+np.outer((np.eye(12)/12).ravel(),e.T.ravel())
    start=np.zeros(144);start[13*seed]=1
    return (np.linalg.matrix_power(superop,512)@start).reshape(12,12)

def main():
    start=time.monotonic();receipt=read(P/'RECEIPT.json');raw=(P/'RECORDS.json.gz').read_bytes()
    assert sha(P/'RECORDS.json.gz')==receipt['records_sha256'];data=json.loads(gzip.decompress(raw))
    for n,h in read(P/'PREREG.json')['sha256'].items():assert sha(P/n)==h,n
    edges={frozenset(x)for x in read(P/'graph.json')['seams']};counts=Counter();outcomes=Counter();censuses={}
    handoffs=[];maxrequests=0;strict=0;physical_reversals=0;zero=0
    for trial in data['trials']:
        initial=[0]*240
        for census in trial['censuses']:
            c=census['cell'];initial[12*c:12*c+12]=reconstruct(census)
        # Audit/source values are compared only AFTER the observer reconstructed them.
        assert initial==trial['initial_audit_only'] and initial[0]==initial[14]==0
        assert max(map(abs,initial))<=64
        k=sum(initial[:12]);assert sum(initial[12:24])==-k and -12<=k<=12
        assert k==trial['harness_K_audit_only']
        censuses[k,trial['law']]=[x['tape']for x in trial['censuses']]
        q=initial[:];states=[q[:]]
        for event in trial['events']:
            a,b,c,j=event['record'];assert frozenset((a,b))in edges
            if event['kind'] in ['retention_move','cleanup_move','orientation_move','frame_move']:
                assert q[a] in(-1,1)and q[b]==0
            if event['kind']=='cleanup_contact':assert q[a]==1 and q[b]==-1
            if event['kind']=='cleanup_move':assert a not in(0,14)and b not in(0,14)
            before=sum(x*x for x in q);advance(q,event);after=sum(x*x for x in q)
            if event['law']=='native':assert after<=before;counts[event['kind']]+=1
            else:assert after==before;counts['swap_hold']+=1
            states.append(q[:])
        assert q==trial['observed_source']==trial['source_audit_only']
        for action in trial['actions']:
            if action['kind']=='strict_repair':
                t=action['event'];a,b=action['request'];before=states[t];after=states[t+1]
                cell=range(0,12)if a<12 else range(12,24)
                assert sum(before[v]**2 for v in cell)==action['potential_before']
                assert sum(after[v]**2 for v in cell)==action['potential_after']
                assert abs(before[a]-before[b])>=2
                if trial['law']=='native':assert action['potential_before']-action['potential_after']>=2;strict+=1
                else:assert action['potential_before']==action['potential_after']
            if 'span'in action:
                lo,hi=action['span'];assert 1<=hi-lo<=64
        outcomes[trial['outcome']]+=1;maxrequests=max(maxrequests,trial['movement_requests'])
        if trial['law']=='native':
            assert all(s['outcome']=='selected'for s in trial['selections'])
            for selection in trial['selections']:assert set(selection['state'])<=({0,1}if sum(selection['state'])>=0 else{-1,0})
            if k:
                expected=[0]*240;expected[0]=1;expected[14]=-1
                assert q==expected and trial['outcome']=='handoff'
                assert trial['cleanup']['pairs']==abs(k)-1
                frame=[e for e in trial['events']if e['kind']=='frame_move'and e['record'][3]]
                assert [e['record'][:2]for e in frame]==[[0,1],[1,0]]
                if k<0:
                    assert any(e['kind']=='orientation_move'and e['record'][3]for e in trial['events']);physical_reversals+=1
                handoffs.append(q)
            else:assert not any(q)and trial['outcome']=='no_unit';zero+=1
        else:assert trial['outcome']=='selection_rejected'
    assert len(handoffs)==24 and physical_reversals==12 and zero==1
    assert all(censuses[k,'native']==censuses[k,'swap_hold']for k in range(-12,13))
    for challenge in data['challenges'].values():
        q=challenge['initial'][:]
        for e in challenge['events']:advance(q,e);counts[e['kind']]+=1
        assert q==challenge['observed']==challenge['audit']
    assert sum(v for k,v in counts.items()if k!='swap_hold')==receipt['literal_native_calls']
    assert counts['swap_hold']==receipt['changed_law_calls']
    # Independently establish connectivity, not just the supplied tree traversal.
    for removed in ({0,14},{0},{1},{14}):
        allowed=set(range(240))-removed;seen={min(allowed)};todo=deque(seen)
        while todo:
            a=todo.popleft()
            for e in edges:
                if a not in e:continue
                b=next(iter(e-{a}))
                if b in allowed and b not in seen:seen.add(b);todo.append(b)
        assert seen==allowed
    # Complete finite endpoint box: the descent and two-coin contrast identity.
    endpoint_cases=0
    for a in range(-64,65):
        for b in range(-64,65):
            for coin in(0,1):
                low=(a+b)//2;high=a+b-low;x,y=(low,high)if coin==0 else(high,low)
                drop=a*a+b*b-x*x-y*y;assert 2*drop==(a-b)**2-((a-b)%2)**2
                j=a-x;opp=x-y;assert 2*j+opp==a-b
                assert min(a,b)<=x<=max(a,b)and min(a,b)<=y<=max(a,b);endpoint_cases+=1
    requests=24*64**2//2*11+22+11*237+3*239+2;assert requests==544020 and maxrequests<=requests
    eta=F(22,2**63)+F(requests,2**64);assert eta==F(receipt['preparation_inconclusive_upper'])<F(3,10**14)
    controller=next(n for n in ast.parse((P/'run.py').read_text()).body if isinstance(n,ast.ClassDef)and n.name=='Controller')
    assert not any(isinstance(n,ast.Attribute)and n.attr=='audit'for n in ast.walk(controller))
    assert not any(isinstance(n,ast.Name)and n.id in('initial','k','K')for n in ast.walk(controller))
    consumer=read(P/'CONSUMER_RESULT.json');maps=np.load(P/'CONSUMER_MAPS.npz')
    native=[np.load(P/'inputs'/n)['certified_Q_proposal']for n in ['first_maps.npz','second_maps.npz']]
    states=[matrix_pump(q,s)for q,s in zip(native,[0,2])]
    difference=max(float(np.max(abs(a-maps[n])))for a,n in zip(states,['first_density','second_density']))
    assert difference<1e-11 and consumer['actual_handoffs_consumed']==24
    common=read(P/'inputs/common_RECEIPT.json')
    assert F(consumer['complete_energy_error_upper'])==F(common['energy_complete_error_upper'])+eta
    assert F(consumer['positive_position_halfturn_TV_upper'])==F(common['spatial_halfturn_TV_upper'])+2*eta
    result={'same_group_only':True,'producer_imports':[],'events_by_kind':dict(counts),'record_only_census_cases':50,
      'strict_native_potential_decreases_checked':strict,'complete_endpoint_cases':endpoint_cases,
      'actual_handoffs':24,'physical_sign_reversals':12,'no_unit_outcomes':1,'controls':25,
      'all_input_fair_request_bound':requests,'complete_inconclusive_upper':str(eta),
      'common_density_independent_power_difference':difference,'all_graph_exclusions_connected':True,
      'bindings':{n:sha(P/n)for n in ['PREREG.json','run.py','observer.py','RECEIPT.json','RECORDS.json.gz','CONSUMER_PREREG.json','join_common.py','CONSUMER_RESULT.json']},
      'seconds':time.monotonic()-start}
    out=HERE/'RESULT.json'
    if out.exists():
        old=read(out);assert {k:v for k,v in old.items()if k!='seconds'}=={k:v for k,v in result.items()if k!='seconds'}
    else:out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k]for k in ['actual_handoffs','physical_sign_reversals','complete_endpoint_cases','common_density_independent_power_difference','seconds']}))
if __name__=='__main__':main()
