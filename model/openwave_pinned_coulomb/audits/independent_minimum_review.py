"""Independent local Ward, compressed-sector and finite Green checks."""
import hashlib,json
from pathlib import Path
import numpy as np
import sympy as s
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
from model.coherent_defect_wave.chiral_joint.periodic_family.model import yukawa_words,kinetic_terms,Configuration
ROOT=Path(__file__).resolve().parents[3]
def review():
    receipt=json.loads((ROOT/'model/openwave_pinned_coulomb/results.json').read_text())
    y=representation()[2][11].diagonal();sy=higgs_memory()[1][11].diagonal()
    assert all(x!=0 for x in y)
    ward=0
    for _,r,c,ops,_ in yukawa_words():
        assert sy[r]-sy[c]+sum((1 if add else-1)*y[i%15]for i,add in ops)==0;ward+=1
    counts=[]
    for n in (2,4,6):
        words=kinetic_terms(Configuration(side=n));assert all(w.left_cell!=w.right_cell for w in words);counts.append(len(words))
    # Build the exact Laplacian by undirected edge addition. Repeated
    # oriented n2 links are intentionally retained, as in the source.
    exact=0
    for n in (2,4):
        N=n**3;L=s.zeros(N);points=list(__import__('itertools').product(range(n),repeat=3))
        for x in points:
            i=points.index(x)
            for a in range(3):
                target=list(x);target[a]=(target[a]+1)%n;j=points.index(tuple(target))
                L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
        P=s.ones(N)/N;G=(L+P).inv()-P
        assert L*G==s.eye(N)-P
        for row in receipt['exact_finite_minimum_rows']:
            if row['side']!=n:continue
            value=G[0,0]-G[0,points.index(tuple(row['cell']))]
            assert value==s.Rational(row['unit_kappa_exact_minimum_difference']);exact+=1
    maximum=0.
    for box in receipt['far_field_minimum_energy_and_static_force']:
        n=box['side'];one=2-2*np.cos(2*np.pi*np.arange(n)/n)
        vals=one[:,None,None]+one[None,:,None]+one[None,None,:];vals[0,0,0]=1
        inverse=1/vals;inverse[0,0,0]=0;green=np.fft.ifftn(inverse).real
        for row in box['rows']:
            r=row['radius'];f=(green[r+1,0,0]-green[r-1,0,0])/2
            maximum=max(maximum,abs(f-row['finite_difference_force']));assert f<0
    assert maximum<1e-14
    paths=['model/openwave_pinned_coulomb/'+f for f in ['PROOF.md','source.py','results.json']]
    return {'verdict':'PASS','native_local_hypercharge_Ward_words':ward,'actual_nonself_hopping_word_counts':counts,
        'independent_exact_Green_minimum_differences':exact,'independent_large_box_force_max_error':maximum,
        'proof_review':'Gauge multiplicity factorization and nonempty native character flux saturate the full pinned tensor lower bound. The unknown complete charged minimum cancels. Pinning is a constraint; no invariance, unrestricted screened force or physical alpha follows.',
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths}}
if __name__=='__main__':
    v=review();Path(__file__).with_suffix('.json').write_text(json.dumps(v,indent=2,sort_keys=True)+'\n');print(v['verdict'])
