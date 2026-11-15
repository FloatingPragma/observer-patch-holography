#!/usr/bin/env python3
"""Independent bounded replay of the retained signed compact-source evidence."""
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    start=time.perf_counter()
    prereg=json.loads((HERE/'PREREG.txt').read_text())
    receipt=json.loads((HERE/'RECEIPT.json').read_text())
    assert receipt['prereg_sha256']==digest(HERE/'PREREG.txt')
    for paths in (prereg['files'],json.loads((HERE/'SOURCE_PINS.json').read_text())):
        for name,sha in paths.items():assert digest(ROOT/name)==sha,name
    for name,sha in receipt['artifact_sha256'].items():assert digest(HERE/name)==sha,name
    graph=json.loads((HERE/'evidence/graph.json').read_text());M=len(graph['seams'])
    f={int(k):v for k,v in graph['f'].items()};g={int(k):v for k,v in graph['g'].items()}
    vectors={'positive':{(p,-1):a for p,a in f.items()},'negative':{(-1,r):b for r,b in g.items()},
        'pair':{(p,r):a*b for p,a in f.items() for r,b in g.items()}}
    outputs={arm:defaultdict(int) for arm in vectors};count=0;contacts=0
    with gzip.open(HERE/'evidence/branches.jsonl.gz','rt') as fobj:
        for line in fobj:
            arm,p,r,ei,coin,current,pp,rr=json.loads(line)
            a,b=graph['seams'][ei]
            q={}
            if p>=0:q[p]=1
            if r>=0:q[r]=-1
            d=q.get(a,0)-q.get(b,0)
            j=(d-(2*coin-1)*(abs(d)%2))//2
            assert current==j
            q[a]=q.get(a,0)-j;q[b]=q.get(b,0)+j
            actual={i:v for i,v in q.items() if v}
            expected={}
            if pp>=0:expected[pp]=1
            if rr>=0:expected[rr]=-1
            assert actual==expected
            if arm=='contact':
                assert not actual;contacts+=1
            else:outputs[arm][pp,rr]+=vectors[arm][p,r]
            count+=1
    for arm,vector in vectors.items():
        lam=12 if arm=='pair' else 6
        assert {k:v for k,v in outputs[arm].items() if v}=={k:(2*M-lam)*v for k,v in vector.items()}
        assert outputs[arm].get((-1,-1),0)==0
    arrays=np.load(HERE/'evidence/raw_source_arrays.npz')
    for name in ('f','g'):
        assert np.linalg.norm(arrays['U_'+name]-arrays['ideal_'+name],'fro')<=1e-11
        assert np.linalg.norm(arrays['inverse_'+name]-arrays['ideal_'+name].conj().T,'fro')<=1e-11
    laws=np.load(HERE/'evidence/complete_outcome_laws.npz')
    for law in ('C','E'):
        for mode,wanted in (('coherent',1.),('dephased',1/16),('sign_plant',1/4)):
            array=laws[f'{law}_energy_{mode}']
            assert array.shape==(256,256) and abs(array.sum()-1)<=5e-10
            assert abs(array[0,0]-wanted)<=5e-10
        a=laws[f'{law}_position_coherent'];b=laws[f'{law}_position_dephased']
        assert np.max(abs(a-b))<=5e-11 and abs(a[4:8,16:20].sum()-1)<=5e-10
    elapsed=time.perf_counter()-start
    assert elapsed<60
    summary={'schema':'oph.localized_signed_rest.reader.v1','all_pass':True,
        'native_branch_records_replayed':count,'contact_branches':contacts,
        'energy_clicks':{'coherent':1,'dephased':1/16,'sign_plant':1/4},
        'elapsed_seconds':elapsed,'row_award':False,'full_row_completion_claimed':False}
    (HERE/'READER_RECEIPT.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps(summary,indent=2,sort_keys=True))

if __name__=='__main__':main()
