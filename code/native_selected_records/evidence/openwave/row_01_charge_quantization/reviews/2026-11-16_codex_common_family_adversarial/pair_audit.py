"""Independent review of two actual native sources; no producer imports."""
from pathlib import Path
import gzip,hashlib,json,time

B=Path(__file__).resolve().parents[2]/'work/2026-11-16_codex_common_refinement_family'
start=time.monotonic();r=json.loads(gzip.decompress((B/'PAIR_SOURCE.json.gz').read_bytes()))
counts={'native_endpoint_audit_calls':0,'joint_columns':0,'independent_preparations':0}
def event(q,edge,coin):
    assert coin in (0,1);a,b=edge;out=list(q);d=q[a]-q[b]
    current=(d-(2*coin-1)*(abs(d)%2))//2
    out[a]-=current;out[b]+=current;counts['native_endpoint_audit_calls']+=1
    return out,current
def replay(q,tape):
    for a,b,c,j in tape:
        q,actual=event(q,(a,b),c);assert actual==j
    return q
def arc(q,N,m,copy,p,h):
    n=24;target=[copy*N+x for x in (0,n,n+1)]
    if h==10:edge=(copy*N+n+5,copy*N+n+6);coin=0
    else:
        a,b=((0,1),(0,1),(0,2),(0,2),(1,2))[h//2]
        edge=(target[a],target[b]);coin=h%2
    out,j=event(q,edge,coin);pp=next(x for x,v in enumerate(target) if out[v]==1)
    hh=10 if h==10 else 2*(h//2)+(coin^bool(j))
    return out,pp,hh,[*edge,coin,j]
results=[]
for g in r['graphs']:
    m=g['M'];N=31+6*m;assert g['vertices']==2*N
    edges={tuple(sorted(x)) for x in g['edges']};assert (0,N) in edges
    reached={0};queue=[0]
    for vertex in queue:
        for a,b in edges:
            other=b if vertex==a else a if vertex==b else None
            if other is not None and other not in reached:reached.add(other);queue.append(other)
    assert len(reached)==2*N
    preps=g['independent_preparations'];assert [x['seed'] for x in preps]==[72315001,72316001]
    expected_tape=[];initial=[]
    for copy,p in enumerate(preps):
        q=p['initial_audit_only'].copy();assert q[:24]==[2,-1]+[0]*22
        initial+=q;out=replay(q,p['body']);out=replay(out,p['probe'])
        assert p['reading']['after_probe']==out[:24]
        assert p['reading']['initial']==q[:24] and p['reading']['total']==1
        assert out[:24].count(1)==1 and set(out[:24])<={0,1}
        out=replay(out,p['route']);assert p['outcome']=='target_basis_unit' and out[24]==1
        assert sum(out)==3
        assert all(tuple(sorted((a+copy*N,b+copy*N))) in edges for a,b,c,j in p['body']+p['probe']+p['route'])
        expected_tape += [[a+copy*N,b+copy*N,c,j] for a,b,c,j in p['body']+p['probe']+p['route']]
        counts['independent_preparations']+=1
    assert initial==g['initial_audit_only'] and expected_tape==g['global_classical_tape']
    final=replay(initial,expected_tape);assert final==g['final_audit_only'] and sum(final)==6
    assert final[24]==final[N+24]==1
    assert len(g['all_joint_native_columns'])==1089
    inputs=set();outputs=[]
    for row in g['all_joint_native_columns']:
        p,h,q,l=row['input'];inputs.add((p,h,q,l));state=[0]*(2*N)
        target=(0,24,25);state[target[p]]=1;state[N+target[q]]=1
        for copy in (0,1):
            for ring0 in (31,31+3*m):state[copy*N+ring0+2*m]=1
        first,pp,hh,r1=arc(state,N,m,0,p,h);last,qq,ll,r2=arc(first,N,m,1,q,l)
        other,qq2,ll2,s2=arc(state,N,m,1,q,l);reverse,pp2,hh2,s1=arc(other,N,m,0,p,h)
        assert last==reverse and (pp,hh,qq,ll)==(pp2,hh2,qq2,ll2)
        assert r1==s1==row['first_source_record'] and r2==s2==row['second_source_record']
        assert row['output']==[pp,hh,qq,ll]
        assert tuple(sorted(r1[:2])) in edges and tuple(sorted(r2[:2])) in edges
        assert sum(last)==6
        # Tensor factorization is checked directly against isolated native arcs.
        iso0=[0]*N;iso0[target[p]]=1;iso1=[0]*N;iso1[target[q]]=1
        _,a,b,_=arc(iso0,N,m,0,p,h);_,c,d,_=arc(iso1,N,m,0,q,l)
        assert [pp,hh,qq,ll]==[a,b,c,d]
        outputs.append((11*pp+hh)*33+11*qq+ll);counts['joint_columns']+=1
    assert len(inputs)==1089 and sorted(outputs)==list(range(1089))
    assert outputs==g['joint_permutation']
    results.append({'M':m,'vertices':2*N,'actual_selected_targets':2,'spectator_marker_units':4,
                    'complete_joint_columns':1089,'native_R_events_per_pair_word':2})
print(json.dumps({'verdict':'PASS','counts':counts,'graphs':results,
 'input_sha256':hashlib.sha256((B/'PAIR_SOURCE.json.gz').read_bytes()).hexdigest(),
 'same_group_only':True,'producer_imports':False,'seconds':time.monotonic()-start},indent=2))
