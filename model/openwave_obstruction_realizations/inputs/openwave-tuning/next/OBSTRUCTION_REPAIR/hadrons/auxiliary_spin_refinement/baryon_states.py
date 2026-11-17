"""Actual frozen-source nucleon and Delta state maps with all current images.

These are states of the complete zero-auxiliary source, used as literal
full625-source trial inputs. They are not eigenstates of the full source.
"""
from pathlib import Path
from itertools import combinations,product
import sys,json
import numpy as np
from scipy import sparse
HERE=Path(__file__).resolve().parent;OLD=HERE.parent;sys.path.insert(0,str(OLD))
from staggered_vacuum import creation,annihilation
from exact_integer_product import exact_product
from spin_flavour_source import quadratic

def words_for(N,s,i):
    choices=list(combinations(range(8),N//3));out=[]
    for sets in product(choices,repeat=3):
        su=0;iu=0
        for row in sets:
            for x in row:su+=int(x%4<2);iu+=int(x%2==0)
        if 2*su-N==s and 2*iu-N==i:out.append(sum(1<<(3*x+c) for c,row in enumerate(sets) for x in row))
    return sorted(out)

def pairs(a,b):return [(3*(4*x+a)+c,3*(4*x+b)+c,1) for x in range(2) for c in range(3)]

def mapping(words,terms,target=None):
    dest={} if target is None else {w:i for i,w in enumerate(target)};rr=[];cc=[];dd=[]
    for j,w in enumerate(words):
        col={}
        for a,b,coef in terms:
            z,sg=annihilation(b,w)
            if z is None:continue
            y,t=creation(a,z)
            if y is not None:col[y]=col.get(y,0)+coef*sg*t
        for y,val in col.items():
            if not val:continue
            if y not in dest:
                assert target is None;dest[y]=len(dest)
            rr.append(dest[y]);cc.append(j);dd.append(val)
    ordered=[0]*len(dest)
    for w,j in dest.items():ordered[j]=w
    B=sparse.csr_matrix((np.array(dd,dtype=np.int64),(rr,cc)),shape=(len(dest),len(words)))
    return ordered,B

def spin_iso_quadratic(words,cell,kind):
    orbitals=range(4*cell,4*cell+4);z=[];plus=[];minus=[]
    for w in words:
        if kind=='spin':z.append(sum((1 if a%4<2 else -1)*((w>>(3*a))&7).bit_count() for a in orbitals))
        else:z.append(sum((1 if a%2==0 else -1)*((w>>(3*a))&7).bit_count() for a in orbitals))
    alphapairs=[(0,2),(1,3)] if kind=='spin' else [(0,1),(2,3)]
    for a,b in alphapairs:
        for c in range(3):plus.append((3*(4*cell+a)+c,3*(4*cell+b)+c,1));minus.append((3*(4*cell+b)+c,3*(4*cell+a)+c,1))
    _,P=mapping(words,plus);_,M=mapping(words,minus)
    return sparse.diags(np.array(z,dtype=np.int64)**2,format='csr',dtype=np.int64)+2*(P.T@P+M.T@M)

def carrier_ops(words):
    n=len(words);N=words[0].bit_count();nl=np.array([(w&((1<<12)-1)).bit_count() for w in words],dtype=np.int64);nr=N-nl
    M=sparse.diags(nl-nr,format='csr',dtype=np.int64);hop=[]
    for a in range(12):hop.extend([(a,a+12,-1),(a+12,a,-1)])
    _,K=mapping(words,hop,words);electric=[]
    for cell,counts in [(0,nl),(1,nr)]:
        U,_=quadratic(words,[[3*(4*cell+a)+c for a in range(4)] for c in range(3)])
        electric.append((3*U-sparse.diags(counts*counts,format='csr',dtype=np.int64)).tocsr())
    G24=sparse.diags(4*(nl*(12-nl)+nr*(12-nr)),format='csr',dtype=np.int64)-electric[0]-electric[1]
    for cell in range(2):G24-=3*spin_iso_quadratic(words,cell,'iso')+spin_iso_quadratic(words,cell,'spin')
    G24=G24.tocsr();G24.eliminate_zeros();assert G24.dtype==np.int64 and (G24-G24.T).nnz==0
    return dict(M=M,K=K,E=electric[0],G24=G24)

def main():
    z=np.load(OLD/'spin_flavour_cache/CAR_frame_3_3_3_0.npz');basewords=z['words'].tolist();baseR=z['R'];old=json.loads((OLD/'spin_flavour_cache/block_3_3_3_0.json').read_text());records=[]
    for name,spin,iso,scale in [('delta',3,3,1),('nucleon',1,1,18)]:
        words=words_for(9,spin,iso);ix={w:i for i,w in enumerate(words)};R=np.zeros((len(words),20),dtype=np.int64)
        if name=='delta':
            for j,w in enumerate(basewords):R[ix[w]]=baseR[j]
        else:
            direct,A=mapping(basewords,pairs(3,0));middle,B=mapping(basewords,pairs(3,2));indirect,C=mapping(middle,pairs(3,1));v=A@baseR;u=C@(B@baseR)
            for j,w in enumerate(direct):R[ix[w]]+=2*v[j]
            for j,w in enumerate(indirect):R[ix[w]]+=u[j]
        assert exact_product(R,R).tolist()==(scale*np.array(old['G'],dtype=object)).tolist()
        _,Sp=mapping(words,pairs(0,2)+pairs(1,3));_,Ip=mapping(words,pairs(0,1)+pairs(2,3))
        assert not np.any(Sp@R) and not np.any(Ip@R)
        colour=[]
        for a,b in [(0,1),(1,2)]:
            _,B=mapping(words,[(3*x+a,3*x+b,1) for x in range(8)]);assert not np.any(B@R);colour.append(B)
        ops=carrier_ops(words);compressed={}
        for tag in ['M','K','E']:
            val=exact_product(R,ops[tag]@R).tolist();assert val==(scale*np.array(old['components'][tag],dtype=object)).tolist();compressed[tag]=val
        assert int(np.max(np.asarray(abs(ops['G24']).sum(axis=1))))*int(abs(R).max())<(1<<62)
        image=ops['G24']@R
        assert not np.any(Sp@image) and not np.any(Ip@image) and all(not np.any(B@image) for B in colour)
        D24=exact_product(R,image).tolist();D24square=exact_product(image,image).tolist()
        np.savez_compressed(HERE/f'CAR_STATE_{name}.npz',words=np.array(words,dtype=np.int64),R=R,G24R=image)
        payload={}
        for tag,X in ops.items():payload.update({tag+'_data':X.data,tag+'_indices':X.indices,tag+'_indptr':X.indptr})
        np.savez_compressed(HERE/f'OPERATORS_{name}.npz',**payload)
        records.append(dict(name=name,N=9,twice_spin=spin,twice_isospin=iso,full_colour_Cartan_spin_isospin_weight_dimension=len(words),
            exact_source_frame_dimension=20,exact_intertwiner_scale=scale,
            actual_spin_isospin_highest_weight_verified=True,actual_colour_raising_kernel_verified=True,
            all_three_source_components_replayed_exactly=True,
            G=(scale*np.array(old['G'],dtype=object)).tolist(),components=compressed,
            mixed_current24=D24,mixed_current24_square=D24square,
            full_current_images_retained_beyond_original_SU4_representation=True))
        print('ACTUAL SOURCE STATE AND ALL CURRENT IMAGES',name,len(words),'rank20',flush=True)
    (HERE/'BARYON_STATES.json').write_text(json.dumps(dict(schema='oph.auxiliary_spin_refinement.actual_baryon_states.v1',
        source='frozen zero-auxiliary two-cell SU3 source',
        nucleon_intertwiner='(2 E30+E31 E32)/sqrt18 from the conjugate Sym3 highest weight',
        actual_full625_eigenstate_claim=False,records=records),indent=2)+'\n')
if __name__=='__main__':main()
