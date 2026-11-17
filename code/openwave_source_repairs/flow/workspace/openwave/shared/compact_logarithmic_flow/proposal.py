#!/usr/bin/env python3
"""Source-bound full spectral proposal, not an acceptance certificate."""
from pathlib import Path
import hashlib,json,math,platform,sys,time
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PROVIDER=HERE.parent/'localized_rest_join'
sys.path.insert(0,str(PROVIDER/'standalone'/'reverse-engineering-reality'))
from model.observer.native import native

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2,allow_nan=False)+'\n')

def main():
    for name in ('SOURCE_PINS.json','PREREG_PROPOSAL.json'):
        for p,digest in json.loads((HERE/name).read_text())['files'].items():
            assert sha(ROOT/p)==digest,p
    for name in ('proposal.json','source_columns.npz','PROPOSAL_RECEIPT.json'):
        if (HERE/name).exists():raise FileExistsError(name)
    started=time.monotonic()
    graph=json.loads((PROVIDER/'evidence'/'graph.json').read_text())
    seams=[tuple(x) for x in graph['seams']];V=240;K=1260;Q=1<<44
    assert graph['ports']==V and len(seams)==630
    support={4,5,6,7};held=next(i for i,(a,b) in enumerate(seams) if (a in support)!=(b in support))
    dest=np.empty((V,K),dtype=np.uint16);cur=np.empty((V,K),dtype=np.int8)
    changed_dest=np.empty_like(dest);changed_cur=np.empty_like(cur)
    counts=np.zeros((V,V),dtype=np.int64);altered=np.zeros_like(counts);calls=0
    for p in range(V):
        for ei,(a,b) in enumerate(seams):
            for c in (0,1):
                out,(coin,j)=native((int(p==a),int(p==b)),(0,1),c);calls+=1
                assert coin==c and out==(int(p==a)-j,int(p==b)+j)
                y=p if p not in (a,b) else (a if out[0] else b)
                yc=p if ei==held else y;jc=0 if ei==held else j;ec=2*ei+c
                dest[p,ec]=y;cur[p,ec]=j;changed_dest[p,ec]=yc;changed_cur[p,ec]=jc
                counts[y,p]+=1;altered[yc,p]+=1
    f=np.zeros(V,dtype=np.int64);f[[4,5,6,7]]=[1,-1,-1,1]
    for matrix in (counts,altered):
        assert np.array_equal(matrix,matrix.T) and np.all(matrix.sum(axis=0)==K)
    assert int(np.trace(counts))==301140 and np.array_equal(counts@f,1254*f)
    proposals={}
    for name,matrix in [('native',counts),('held_seam',altered)]:
        values,vectors=np.linalg.eigh(matrix.astype(float)/K)
        for j in range(V):
            pivot=int(np.argmax(abs(vectors[:,j])))
            if vectors[pivot,j]<0:vectors[:,j]*=-1
        proposals[name]={'columns':np.rint(vectors.T*Q).astype(np.int64).tolist(),
                         'eigenvalues':np.rint(values*Q).astype(np.int64).tolist()}
    write(HERE/'proposal.json',{'schema':'oph.compact_log_flow.proposal.v1','denominator':Q,
          'source_dimension':V,'laws':proposals,'held_seam_index':held,'held_seam':list(seams[held]),
          'proposer':'numpy.linalg.eigh; maximum-coordinate sign convention; dyadic rounding',
          'numpy_version':np.__version__,'python_version':platform.python_version()})
    np.savez_compressed(HERE/'source_columns.npz',destination=dest,current=cur,
                        changed_destination=changed_dest,changed_current=changed_cur,
                        counts=counts,changed_counts=altered,seams=np.asarray(seams,dtype=np.int16))
    write(HERE/'PROPOSAL_RECEIPT.json',{'native_calls':calls,'retained_primary_columns':V*K,
          'retained_changed_law_columns':V*K,'held_seam_index':held,'held_seam':list(seams[held]),
          'primary_trace_counts':int(np.trace(counts)),'changed_trace_counts':int(np.trace(altered)),
          'source_denominator':K,'eigenvectors_retained_per_law':V,'no_flow_acceptance_test_executed':True,
          'prereg_sha256':sha(HERE/'PREREG_PROPOSAL.json'),
          'artifacts':{n:sha(HERE/n) for n in ('proposal.json','source_columns.npz')},
          'elapsed_seconds':time.monotonic()-started})
    print((HERE/'PROPOSAL_RECEIPT.json').read_text())

if __name__=='__main__':main()
