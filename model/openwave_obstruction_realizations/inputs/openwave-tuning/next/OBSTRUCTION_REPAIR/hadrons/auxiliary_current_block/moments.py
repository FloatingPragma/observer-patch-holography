"""Exact moments for the full-source current-image block, with all complements.

The source is unchanged. Every Cartesian adjoint/vector commutator output and
both-cell current image is retained; no auxiliary output is postselected.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib
import numpy as np
from scipy import sparse
HERE=Path(__file__).resolve().parent;OLD=HERE.parent;PREV=OLD/'auxiliary_spin_refinement'
sys.path.insert(0,str(OLD));sys.path.insert(0,str(PREV))
from baryon_states import mapping,spin_iso_quadratic
from spin_flavour_source import quadratic
from exact_integer_product import exact_product

def integer_generators():
    colour=[]
    for a,b in [(0,1),(0,2),(1,2)]:
        x=np.zeros((3,3),dtype=np.int64);x[a,b]=x[b,a]=1;colour.append(x)
        y=np.zeros((3,3),dtype=np.int64);y[a,b]=-1;y[b,a]=1;colour.append(y)
    colour.extend([np.diag([1,-1,0]).astype(np.int64),np.diag([1,1,-2]).astype(np.int64)])
    spin=[np.array([[0,1],[1,0]],dtype=np.int64),np.array([[0,-1],[1,0]],dtype=np.int64),np.diag([1,-1]).astype(np.int64)]
    return colour,spin

def local_D24(words,cell):
    counts=np.array([((w>>(12*cell))&((1<<12)-1)).bit_count() for w in words],dtype=np.int64)
    U,_=quadratic(words,[[3*(4*cell+a)+c for a in range(4)] for c in range(3)])
    C6=3*U-sparse.diags(counts*counts,format='csr',dtype=np.int64)
    D=sparse.diags(4*counts*(12-counts),format='csr',dtype=np.int64)-C6-3*spin_iso_quadratic(words,cell,'iso')-spin_iso_quadratic(words,cell,'spin')
    D=D.tocsr();D.eliminate_zeros();assert D.dtype==np.int64 and (D-D.T).nnz==0
    return D

def main():
    previous=json.loads((PREV/'BARYON_STATES.json').read_text())['records'];colour,spin=integer_generators();out=[]
    for state in previous:
        name=state['name'];z=np.load(PREV/f'CAR_STATE_{name}.npz');words=z['words'].tolist();R=z['R'];images=[]
        for cell in range(2):
            D=local_D24(words,cell);assert int(np.max(np.asarray(abs(D).sum(axis=1))))*int(abs(R).max())<(1<<62);images.append(D@R)
        assert np.array_equal(images[0]+images[1],z['G24R'])
        L576=exact_product(images[0],images[1]);assert L576.tolist()==L576.T.tolist()
        assert exact_product(images[0]+images[1],images[0]+images[1]).tolist()==state['mixed_current24_square']
        total=np.zeros((20,20),dtype=object);currents=[]
        for a,T in enumerate(colour):
            for i,S in enumerate(spin):
                h=np.kron(np.kron(S,np.eye(2,dtype=np.int64)),T);terms=[]
                for x,y in zip(*np.nonzero(h)):
                    terms.append((int(x),int(y)+12,int(h[x,y])));terms.append((int(x)+12,int(y),-int(h[x,y])))
                dest,B=mapping(words,terms);image=B@R
                assert int(np.max(np.asarray(abs(B).sum(axis=1))))*int(abs(R).max())<(1<<62)
                gram=exact_product(image,image)
                # Actual g_ai has magnitude integer_matrix/4, with one extra
                # sqrt3 in the eighth colour generator. Imaginary phases drop
                # from each norm. Both cell commutators have opposite signs.
                denominator=1536 if a==7 else 512
                total+=gram*(1 if a==7 else 3)
                currents.append(dict(colour_generator=a,spin_generator=i,
                    complete_output_occupation_count=len(dest),output_words_sha256=hashlib.sha256(np.array(dest,dtype=np.int64).tobytes()).hexdigest(),
                    exact_integer_Gram=gram.tolist(),coefficient_in_CdaggerC=f'1/{denominator}'))
                print('ALL HOP-CURRENT OUTPUTS',name,a,i,len(dest),flush=True)
        rec=dict(name=name,N=9,twice_spin=state['twice_spin'],twice_isospin=state['twice_isospin'],
            all24_current_commutator_records=currents,
            CdaggerC1536=total.tolist(),Dleft_Dright576=L576.tolist(),
            exact_Dsum_images_replayed=True,all_auxiliary_one_and_two_outputs_retained=True,
            definitions=dict(C='[H_hop,A], all one-auxiliary outputs',
                CdaggerC='2t² sum_ai B_ai†B_ai, B_ai=cL†g_ai cR-cR†g_ai cL',
                Dleft_Dright='two-auxiliary squared image; no current compression before multiplication'))
        out.append(rec)
    (HERE/'MOMENTS.json').write_text(json.dumps(dict(schema='oph.auxiliary_current_block.moments.v1',source_unchanged=True,records=out),indent=2)+'\n')
if __name__=='__main__':main()
