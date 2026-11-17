"""Independent global integer check of the two-source preparation and clock."""
import gzip
import hashlib
import json
from pathlib import Path
import time


def endpoint(q,edge,c):
    a,b=edge;out=list(q);total=q[a]+q[b];lo=total//2
    out[a],out[b]=(lo,total-lo) if c==0 else (total-lo,lo)
    return tuple(out),q[a]-out[a]


def run(folder):
    start=time.perf_counter();blob=(folder/'PAIR_SOURCE.json.gz').read_bytes()
    d=json.loads(gzip.decompress(blob));receipt=json.loads((folder/'PAIR_RECEIPT.json').read_text())
    assert hashlib.sha256(blob).hexdigest()==receipt['evidence_sha256'];columns=tapes=0
    for g in d['graphs']:
        n=g['vertices']//2;m=g['M'];target=[0,24,25];state=tuple(g['initial_audit_only'])
        for a,b,c,j in g['global_classical_tape']:
            state,current=endpoint(state,(a,b),c);assert current==j;tapes+=1
        assert state==tuple(g['final_audit_only']) and sum(state)==6
        assert state[24]==state[n+24]==1
        assert [r['seed'] for r in g['independent_preparations']]==[72315001,72316001]
        rows=g['all_joint_native_columns'];assert len(rows)==1089;images=[]
        for row in rows:
            p,h,q,l=row['input'];load=[0]*(2*n);load[target[p]]=load[n+target[q]]=1
            # Four marker units are spectator integer coordinates.
            for shift in (0,n):
                for first in (31,31+3*m):load[shift+first+2*m]=1
            before=tuple(load);outputs=[]
            for rec in (row['first_source_record'],row['second_source_record']):
                a,b,c,j=rec;before,current=endpoint(before,(a,b),c);assert current==j
                outputs.append((a,b,c,j))
            pp,hh,qq,ll=row['output']
            assert before[target[pp]]==before[n+target[qq]]==1
            expected_h=10 if h==10 else 2*(h//2)+((h%2)^int(outputs[0][3]!=0))
            expected_l=10 if l==10 else 2*(l//2)+((l%2)^int(outputs[1][3]!=0))
            assert (hh,ll)==(expected_h,expected_l)
            # Replay reversed global native order and compare the entire load.
            reverse=tuple(load)
            for a,b,c,j in reversed(outputs):reverse,current=endpoint(reverse,(a,b),c);assert current==j
            assert before==reverse
            images.append((11*pp+hh)*33+11*qq+ll);columns+=1
        assert images==g['joint_permutation'] and sorted(images)==list(range(1089))
    return {'passed':True,'same_group':True,'external_G6':False,'joint_columns':columns,
            'global_classical_events':tapes,'evidence_sha256':hashlib.sha256(blob).hexdigest(),
            'seconds':time.perf_counter()-start}


if __name__=='__main__':print(json.dumps(run(Path(__file__).resolve().parent),indent=2,sort_keys=True))
