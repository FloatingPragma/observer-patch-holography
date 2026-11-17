"""Independent record verification; imports no producer or native module."""
import gzip
import hashlib
import json
from pathlib import Path
from fractions import Fraction

HERE=Path(__file__).resolve().parent


def update(q,record,law):
    a,b,c,j=record; old=q[a]+q[b]
    if law=='native':
        lower=old//2; upper=old-lower
        left,right=(upper,lower) if c else (lower,upper)
    elif law=='swap_hold':
        left,right=(q[b],q[a]) if c else (q[a],q[b])
    else: raise AssertionError(law)
    assert j==q[a]-left and q[b]+j==right
    q[a],q[b]=left,right


def read_probe(tape,initial,cell):
    # Solve independent initial differences while replaying probe incidence.
    base=12*cell; offsets=[0]*240; differences=[];i=0
    while i<len(tape):
        a,b,c,j=tape[i];d0=offsets[a]-offsets[b]
        first=c; firstj=j
        offsets[a]-=j;offsets[b]+=j;i+=1
        while i<len(tape):
            aa,bb,coin,current=tape[i];assert(aa,bb)==(a,b)
            offsets[a]-=current;offsets[b]+=current;i+=1
            if coin!=first:
                differences.append((a,b,2*firstj+current-d0));break
        else:raise AssertionError('incomplete census')
    values={(0,14)[cell]:0}
    while len(values)<12:
        old=len(values)
        for a,b,d in differences:
            if a in values: values[b]=values[a]-d
            if b in values: values[a]=values[b]+d
        assert len(values)>old
    recovered=[values[v] for v in range(base,base+12)]
    assert recovered==initial[base:base+12]
    return recovered,[values[v]+offsets[v] for v in range(base,base+12)]


def main():
    receipt=json.loads((HERE/'RECEIPT.json').read_text())
    raw=(HERE/'RECORDS.json.gz').read_bytes()
    assert hashlib.sha256(raw).hexdigest()==receipt['records_sha256']
    data=json.loads(gzip.decompress(raw)); events=0; native_outcomes={};controls=[]
    handoffs=[]
    for trial in data['trials']:
        q=list(trial['initial_audit_only']); initial=list(q)
        for event in trial['events']:
            update(q,event['record'],event['law']);events+=1
        assert q==trial['source_audit_only']==trial['observed_source']
        for census in trial['censuses']:
            before,after=read_probe(census['tape'],initial,census['cell'])
            assert census['decoded']['initial']==before
            assert census['decoded']['terminal']==after
            assert census['decoded']['total']==sum(before)
        # Reconstruct the complete initial source backward from the terminal source.
        reverse=list(q)
        for event in reversed(trial['events']):
            a,b,c,j=event['record'];reverse[a]+=j;reverse[b]-=j
        assert reverse==initial
        k=trial['harness_K_audit_only']
        assert sum(initial[:12])==k and sum(initial[12:24])==-k
        if trial['law']=='native':
            native_outcomes[k]=trial['outcome']
            assert all(s['outcome']=='selected' for s in trial['selections'])
            for action in trial['actions']:
                if action['kind']=='strict_repair':
                    assert action['potential_after']<=action['potential_before']-2
                if 'span' in action:
                    lo,hi=action['span'];assert 1<=hi-lo<=64
            if k:
                wanted=[0]*240;wanted[0]=1;wanted[14]=-1
                assert q==wanted and trial['pair_code']==14 and trial['positive_frame']==[0,1]
                assert trial['cleanup']['pairs']==abs(k)-1
                assert trial['initial_orientation']==(1 if k>0 else -1)
                if k<0:assert any(e['kind']=='orientation_move' for e in trial['events'])
                handoffs.append(dict(harness_K_audit_only=k,actual_source=q,
                    observed_source=trial['observed_source'],pair_code=14,positive_frame=[0,1],
                    complete_prefix_sha256=hashlib.sha256(json.dumps(trial['events'],sort_keys=True,separators=(',',':')).encode()).hexdigest(),
                    numerical_total_supplied_to_observer=False,source_reset=False))
        else:
            assert trial['outcome']=='selection_rejected'
            controls.append(k)
    assert len(native_outcomes)==25 and native_outcomes[0]=='no_unit'
    assert all(v=='handoff' for k,v in native_outcomes.items() if k)
    assert len(controls)==25 and len(handoffs)==24
    edges={frozenset(e) for e in json.loads((HERE/'graph.json').read_text())['seams']}
    for name,tree in data['graph_spanning_certificates'].items():
        blocked={0,14} if name=='cleanup_238' else {int(name.split('avoid')[1])}
        vertices=set(range(240))-blocked
        reached={min(vertices)}
        for a,b in tree:
            assert a in reached and b not in reached and frozenset((a,b)) in edges
            reached.add(b)
        assert reached==vertices
    for name,challenge in data['challenges'].items():
        q=list(challenge['initial'])
        for event in challenge['events']:update(q,event['record'],event['law'])
        assert q==challenge['observed']==challenge['audit']
    assert len(data['challenges']['forced_native_timeout']['events'])==64
    assert not data['challenges']['invalid_route_requests']['events']
    bound=Fraction(22,2**63)+Fraction(544020,2**64)
    assert Fraction(receipt['preparation_inconclusive_upper'])==bound< Fraction(3,10**14)
    result=dict(independent_events_replayed=events,unknown_totals_checked=25,
        actual_canonical_handoffs=24,reverse_orientation_handoffs=12,
        no_unit_outcomes=1,changed_law_rejections=25,
        complete_preparation_failure_upper=str(bound),producer_imports=False)
    (HERE/'INDEPENDENT_RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    (HERE/'HANDOFFS.json').write_text(json.dumps(dict(candidate=receipt['candidate'],
        producer_records_sha256=receipt['records_sha256'],handoffs=handoffs),indent=2)+'\n')
    print(json.dumps(result,sort_keys=True))


if __name__=='__main__':main()
