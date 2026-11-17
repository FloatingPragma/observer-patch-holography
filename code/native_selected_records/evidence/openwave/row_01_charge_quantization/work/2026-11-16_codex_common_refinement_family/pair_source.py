"""Two independently selected common sources in one native integer graph."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import core


def replay(initial,tape):
    q=tuple(initial)
    for a,b,c,j in tape:
        q,record=core.native(q,(a,b),c)
        assert record==[a,b,c,j]
    return q


def event(q,g,copy,position,choice):
    offset=copy*g['N'];target=[offset+x for x in g['target']]
    if choice==10:
        edge=tuple(offset+x for x in g['blank']);coin=0
    else:
        edges=((0,1),(0,1),(0,2),(0,2),(1,2));a,b=edges[choice//2]
        edge=(target[a],target[b]);coin=choice%2
    out,record=core.native(q,edge,coin);j=record[3]
    new_position=next(i for i,v in enumerate(target) if out[v]==1)
    new_choice=10 if choice==10 else 2*(choice//2)+(coin^int(j!=0))
    return out,new_position,new_choice,record


def run():
    begin=time.perf_counter();start=dict(core.CALLS);graphs=[]
    for m in (4,8):
        g=core.graph(2,m);N=g['N'];source_records=[]
        for seed in (72315001,72316001):
            r=core.charge_route(g,'native',1,seed)
            assert r['outcome']=='target_basis_unit';source_records.append(r)
        initial=source_records[0]['initial_audit_only']+source_records[1]['initial_audit_only']
        tape=[]
        for copy,r in enumerate(source_records):
            offset=copy*N
            tape += [[a+offset,b+offset,c,j] for a,b,c,j in r['body']+r['probe']+r['route']]
        final=replay(initial,tape)
        assert sum(final)==6 and final[g['target'][1]]==final[N+g['target'][1]]==1
        edges=list(g['edges'])+[(a+N,b+N) for a,b in g['edges']]+[(0,N)]
        core.tree(2*N,edges)
        edge_set={tuple(sorted(e)) for e in edges}
        rows=[];permutation=[]
        for p in range(3):
            for h in range(11):
                for q in range(3):
                    for l in range(11):
                        state=[0]*(2*N)
                        state[g['target'][p]]=1;state[N+g['target'][q]]=1
                        for copy in (0,1):
                            for ring in g['rings']:state[copy*N+ring[2*m]]=1
                        first,pp,hh,r1=event(state,g,0,p,h)
                        out,qq,ll,r2=event(first,g,1,q,l)
                        opposite,q2,l2,s2=event(state,g,1,q,l)
                        reversed_out,p2,h2,s1=event(opposite,g,0,p,h)
                        assert out==reversed_out and (pp,hh,qq,ll)==(p2,h2,q2,l2)
                        assert r1==s1 and r2==s2
                        assert all(tuple(sorted(r[:2])) in edge_set for r in (r1,r2))
                        output_index=(11*pp+hh)*33+11*qq+ll;permutation.append(output_index)
                        rows.append({'input':[p,h,q,l],'output':[pp,hh,qq,ll],
                                     'first_source_record':r1,'second_source_record':r2})
        assert sorted(permutation)==list(range(1089))
        R=core.clock(0)['R'];expected=np.argmax(np.kron(R,R),axis=0).tolist()
        assert permutation==expected
        graphs.append({'M':m,'vertices':2*N,'idle_intercopy_bridge':[0,N],
            'edges':edges,'independent_preparations':source_records,'global_classical_tape':tape,
            'initial_audit_only':initial,'final_audit_only':final,'all_joint_native_columns':rows,
            'joint_permutation':permutation,'joint_column_count':len(rows),
            'R_calls_per_pair_word':2,'graph_connected':True})
    return {'candidate':'A-native-common-refinement-family-v1','domain_supplement':'two independently prepared local-unit sources, connected only by idle bridge',
        'graphs':graphs,'passed':True,'literal_native_calls':core.CALLS['native']-start['native'],
        'seconds':time.perf_counter()-begin,'same_group':True,'physical_award':False}


if __name__=='__main__':
    result=run();payload=gzip.compress((json.dumps(result,sort_keys=True,separators=(',',':'))+'\n').encode(),mtime=0)
    Path('PAIR_SOURCE.json.gz').write_bytes(payload)
    summary={k:result[k] for k in ('candidate','domain_supplement','passed','literal_native_calls','seconds','same_group','physical_award')}
    summary['evidence_sha256']=hashlib.sha256(payload).hexdigest();summary['joint_columns']=2178
    Path('PAIR_RECEIPT.json').write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n');print(json.dumps(summary,indent=2))
