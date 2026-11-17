"""Independently reconstruct all controller decisions from public record tapes."""
from pathlib import Path
from collections import Counter,deque
from fractions import Fraction as F
import json

HERE=Path(__file__).resolve().parent
traces=json.loads((HERE/'TRAJECTORIES.json').read_text())
report=json.loads((HERE/'TRAJECTORY_RECEIPT.json').read_text())
graph=json.loads((HERE/'inputs/graph.json').read_text())
edges=[tuple(e) for e in graph['seams'] if max(e)<12]
adj={i:[] for i in range(12)}
for a,b in edges:
    adj[a].append(b)
    adj[b].append(a)


def toward(p,targets):
    d={i:0 for i in targets}
    q=deque(targets)
    while q:
        v=q.popleft()
        for u in adj[v]:
            if u not in d:d[u]=d[v]+1;q.append(u)
    return tuple(sorted((p,min(v for v in adj[p] if d[v]==d[p]-1))))


def check(t):
    initial,sign=t['preparation']
    loads=[sign if i==initial else 0 for i in range(14)]
    h=tail=None
    mode='SEARCH'
    gate=None
    positions=[initial]
    work=[]
    for index,(edge,coin,j,tag,snap) in enumerate(t['events']):
        assert type(coin)is int and coin in (0,1)
        assert type(j)is int
        a,b=edge
        before=positions[-1]
        if tag=='search':
            assert h is None and tuple(edge) in edges
        elif tag=='gate':
            assert h is not None and edge==[12,13] and before in (h,tail)
        elif tag.startswith('hold_'):
            expected=tuple(sorted((h,tail))) if before==tail or gate==0 else (12,13)
            assert tuple(edge)==expected
            assert tag==('hold_idle' if expected==(12,13) else 'hold_active')
        elif tag=='kick':
            assert tuple(edge)==toward(before,[t['kick_target']])
        elif tag=='return':
            assert tuple(edge)==toward(before,[h,tail])
        elif tag=='work':
            assert tuple(edge)==tuple(sorted((h,tail))) and before in (h,tail)
        else:raise AssertionError('Unknown tag')
        if t['held_control']:
            assert j==0
        else:
            total=loads[a]+loads[b]
            low=total//2
            high=total-low
            na,nb=(low,high) if coin==0 else (high,low)
            assert j==loads[a]-na
            loads[a],loads[b]=na,nb
        after,=[i for i,v in enumerate(loads) if v]
        assert loads[after]==sign
        if tag=='search' and j!=0:
            assert h is None
            h,tail=after,before
            mode='HOLD'
        elif tag=='gate':gate=coin
        elif tag=='kick':mode='PERTURBED'
        elif tag=='return':mode='HOLD' if after in (h,tail) else 'RETURN'
        assert snap==[after,h,tail,mode,gate]
        if tag=='work':work.append([index,int(before==tail),int(after==tail),int(after==tail)-int(before==tail)])
        positions.append(after)
    assert work==t['work']
    if not t['held_control'] and 'return_start'in t:
        assert t['return_end']-t['return_start']<=32
        assert positions[t['return_end']] in (h,tail)
    return positions


for t in traces['trials']+traces['held_controls']:check(t)
for p in traces['pairs']:
    pa,pb=check(p['a']),check(p['b'])
    start_a=next(i+1 for i,e in enumerate(p['a']['events']) if e[4][1]is not None)+16
    start_b=next(i+1 for i,e in enumerate(p['b']['events']) if e[4][1]is not None)+16
    for i,(before,after,ea,eb) in enumerate(p['joint']):
        ha,ta=p['a']['events'][ea-1][4][1:3]
        hb,tb=p['b']['events'][eb-1][4][1:3]
        assert after==[int(pa[ea]==ta),int(pb[eb]==tb)]
        previous_a,previous_b=(p['joint'][i-1][2:4] if i else (start_a,start_b))
        assert before==[int(pa[previous_a]==ta),int(pb[previous_b]==tb)]
        if i:assert before==p['joint'][i-1][1]
        delta=(ea-previous_a,eb-previous_b)
        if p['exclusive']:
            c,j=p['selector_records'][i]
            assert type(c)is int and c in (0,1) and j==0
            assert delta==((2,0) if c==0 else (0,2))
        else:assert delta==(2,2)
    assert len(p['joint'])==2048
counts=Counter(w[3] for t in traces['trials'] for w in t['work'])
assert sum(counts.values())==report['work_draws']==768
assert report['work_counts']=={str(k):counts[k] for k in (-1,0,1)}
assert report['empirical_exponential_work_mean']==str(sum(F(counts[k],768)*F(2)**(-k) for k in (-1,0,1)))
assert report['native_event_records']==sum(len(t['events']) for t in traces['trials'])+sum(len(p['a']['events'])+len(p['b']['events'])+len(p['selector_records']) for p in traces['pairs'])
assert len(traces['trials'])==len(traces['held_controls'])==24
assert all(t['failures']==['discovery_timeout'] for t in traces['held_controls'])
assert not report['physical_rest_energy_established']
print(json.dumps({'record_only_controller_reconstructed':True,'work_counts':dict(counts),
                  'native_event_records':report['native_event_records']},sort_keys=True))
