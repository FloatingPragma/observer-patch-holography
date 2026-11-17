"""Exact CAR reduction of a two-cell, spin/flavour SU4 invariant colour source.

Highest-weight multiplicities come from the separately executed exact character
census. Projectors are polynomials in the actual full-CAR Casimirs, not labels.
All operator entries, cyclic bases and compressed congruences are integer.
"""
from pathlib import Path
from itertools import combinations,product
from math import gcd,prod
from functools import reduce
from collections import deque
import json,time
import numpy as np
from scipy import sparse
from staggered_vacuum import creation,annihilation
from exact_integer_product import exact_product
HERE=Path(__file__).resolve().parent
CACHE=HERE/'spin_flavour_cache';CACHE.mkdir(exist_ok=True)
PRIME=1000003

def key(lam):return '_'.join(map(str,lam))

def move(word,i,j):
    v,s=annihilation(j,word)
    if v is None:return None,0
    z,t=creation(i,v)
    return z,s*t

def weight_words(lam):
    N=sum(lam);k=N//3
    choices=list(combinations(range(8),k))
    out=[]
    for sets in product(choices,repeat=3):
        counts=[0]*4
        for sites in sets:
            for a in sites:counts[a%4]+=1
        if counts==list(lam):out.append(sum(1<<(3*a+c) for c,sites in enumerate(sets) for a in sites))
    return sorted(out)

def rectangular(words,pairs):
    dest={};rr=[];cc=[];dd=[]
    for j,w in enumerate(words):
        col={}
        for a,b in pairs:
            z,s=move(w,a,b)
            if z is not None:col[z]=col.get(z,0)+s
        for z,s in col.items():
            if not s:continue
            if z not in dest:dest[z]=len(dest)
            rr.append(dest[z]);cc.append(j);dd.append(s)
    return sparse.csr_matrix((np.array(dd,dtype=np.int64),(rr,cc)),shape=(len(dest),len(words)))

def quadratic(words,groups):
    # groups[label] is a same-size list of CAR modes; all coordinates are kept.
    n=len(words);r=len(groups);counts=np.array([[sum((w>>i)&1 for i in g) for g in groups] for w in words])
    Q=sparse.diags(np.sum(counts*counts,axis=1),format='csr',dtype=np.int64)
    for a,b in product(range(r),repeat=2):
        if a==b:continue
        B=rectangular(words,list(zip(groups[a],groups[b])))
        Q=Q+B.T@B
    return Q.tocsr(),counts

def operators(lam,expected_weight):
    name=key(lam);filename=CACHE/f'operators_{name}.npz'
    wordsfile=CACHE/f'words_{name}.json'
    if filename.exists():
        z=np.load(filename);words=json.loads(wordsfile.read_text());ops={}
        for tag in ['M','K','E','C','F']:
            ops[tag]=sparse.csr_matrix((z[tag+'_data'],z[tag+'_indices'],z[tag+'_indptr']),shape=(len(words),)*2)
        return words,ops
    words=weight_words(lam);assert len(words)==expected_weight,(lam,len(words),expected_weight)
    n=len(words);N=sum(lam);ix={w:i for i,w in enumerate(words)}
    nleft=np.array([(w&((1<<12)-1)).bit_count() for w in words],dtype=np.int64)
    M=sparse.diags(2*nleft-N,format='csr',dtype=np.int64)
    K=sparse.lil_matrix((n,n),dtype=np.int64)
    for j,w in enumerate(words):
        for a in range(12):
            for i,h in [(a,a+12),(a+12,a)]:
                z,s=move(w,i,h)
                if z is not None:K[ix[z],j]-=s
    C,_=quadratic(words,[[3*a+c for a in range(8)] for c in range(3)])
    C=3*C-sparse.eye(n,format='csr',dtype=np.int64)*N*N
    F,_=quadratic(words,[[3*(4*x+a)+c for x in range(2) for c in range(3)] for a in range(4)])
    E,lc=quadratic(words,[[3*a+c for a in range(4)] for c in range(3)])
    E=3*E-sparse.diags(nleft*nleft,format='csr',dtype=np.int64)
    ops=dict(M=M,K=K.tocsr(),E=E.tocsr(),C=C.tocsr(),F=F.tocsr())
    payload={}
    for tag,X in ops.items():
        assert (X-X.T).nnz==0
        payload.update({tag+'_data':X.data,tag+'_indices':X.indices,tag+'_indptr':X.indptr})
    np.savez_compressed(filename,**payload);wordsfile.write_text(json.dumps(words))
    return words,ops

def primitive(v):
    g=reduce(gcd,(abs(int(x)) for x in v if x),0)
    if not g:return np.zeros(len(v),dtype=np.int64)
    vals=[int(x)//g for x in v]
    if max(map(abs,vals))<(1<<61):return np.array(vals,dtype=np.int64)
    return np.array(vals,dtype=object)

def exact_action(A,v):
    if A.shape[0]==0:return np.zeros(0,dtype=np.int64)
    vmax=max((abs(int(x)) for x in v),default=0)
    rowmax=int(np.max(np.array(abs(A).sum(axis=1)).ravel()))
    if vmax*rowmax<(1<<61):return A@np.array(v,dtype=np.int64)
    return np.array([sum(int(A.data[p])*int(v[A.indices[p]]) for p in range(A.indptr[i],A.indptr[i+1])) for i in range(A.shape[0])],dtype=object)

def u_casimir(lam):
    n=len(lam)
    return sum(a*(a+n+1-2*(i+1)) for i,a in enumerate(lam))

def projection(v,ops,lam,all_flavour):
    N=sum(lam);k=N//3
    cs=sorted({3*u_casimir(mu)-N*N for mu in product(range(9),repeat=3)
               if sum(mu)==N and list(mu)==sorted(mu,reverse=True) and mu!=(k,k,k)})
    for ev in cs:
        v=primitive(exact_action(ops['C'],v)-ev*v)
        if not np.any(v):return v
    assert not np.any(exact_action(ops['C'],v))
    target=u_casimir(lam)
    fs=sorted({u_casimir(mu) for mu in all_flavour if tuple(mu)!=tuple(lam)
       and all(sum(mu[:i])>=sum(lam[:i]) for i in range(1,5))})
    assert all(x>target for x in fs)
    for ev in fs:
        v=primitive(exact_action(ops['F'],v)-ev*v)
        if not np.any(v):return v
    assert np.array_equal(exact_action(ops['F'],v),target*v)
    return v

class ModRank:
    def __init__(self):self.rows=[];self.pivots=[]
    def add(self,v):
        q=np.array([int(x)%PRIME for x in v],dtype=np.int64)
        for p,row in zip(self.pivots,self.rows):q=(q-int(q[p])*row)%PRIME
        nz=np.flatnonzero(q)
        if not len(nz):return False
        p=int(nz[0]);q=q*pow(int(q[p]),-1,PRIME)%PRIME
        self.rows.append(q);self.pivots.append(p);return True

def exact_dot(v,w):return sum(int(a)*int(b) for a,b in zip(v,w) if a and b)

def build_block(row,all_rows):
    lam=row['U4_highest_weight'];expected=row['multiplicity'];name=key(lam)
    done=CACHE/f'block_{name}.json'
    if done.exists():return json.loads(done.read_text())
    print('BUILD',lam,'carrier',row['joint_colour_and_flavour_weight_dimension'],'rank',expected,flush=True)
    words,ops=operators(lam,row['joint_colour_and_flavour_weight_dimension']);n=len(words)
    rank=ModRank();basis=[];q=deque();rng=np.random.default_rng(314159+sum((i+1)*a for i,a in enumerate(lam)))
    def add(v):
        v=primitive(v)
        if rank.add(v):basis.append(v);q.append(v);return True
        return False
    for seed in range(expected+1):
        if len(basis)==expected:break
        if seed==0 and n==1:
            v=np.ones(1,dtype=np.int64)
        elif seed==0 and lam==[3,3,3,3]:
            v=np.zeros(n,dtype=np.int64);v[words.index(((1<<12)-1)<<12)]=1
        elif seed==0 and lam==[6,3,3,3]:
            v=np.zeros(n,dtype=np.int64);v[words.index((((1<<12)-1)<<12)|7)]=1
        else:
            v=projection(rng.integers(-1,2,n,dtype=np.int64),ops,lam,[r['U4_highest_weight'] for r in all_rows])
        add(v)
        while q and len(basis)<expected:
            v=q.popleft()
            for tag in ['M','E','K']:
                add(exact_action(ops[tag],v))
                if len(basis)==expected:break
        print('  projected seeds',seed+1,'rank',len(basis),'max coefficient',max((abs(int(x)) for v in basis for x in v),default=0),flush=True)
    assert len(basis)==expected
    # Explicit native highest-weight constraints; no eigenvalue-label substitution.
    for a in range(2):
        B=rectangular(words,[(3*x+a,3*x+a+1) for x in range(8)])
        assert all(not np.any(exact_action(B,v)) for v in basis)
    for a in range(3):
        B=rectangular(words,[(3*(4*x+a)+c,3*(4*x+a+1)+c) for x in range(2) for c in range(3)])
        assert all(not np.any(exact_action(B,v)) for v in basis)
    basis_matrix=np.array(basis,dtype=np.int64).T
    G=exact_product(basis_matrix,basis_matrix).tolist()
    compressed={}
    for tag in ['M','K','E']:
        images=[exact_action(ops[tag],v) for v in basis]
        # Every image is still a colour singlet and flavour highest weight.
        assert all(not np.any(exact_action(ops['C'],v)) for v in images)
        assert all(np.array_equal(exact_action(ops['F'],v),u_casimir(lam)*v) for v in images)
        A=exact_product(basis_matrix,np.array(images,dtype=np.int64).T).tolist()
        assert A==list(map(list,zip(*A)))
        compressed[tag]=A
    result=dict(highest_weight=lam,SU4_dimension=row['SU4_dimension'],multiplicity=expected,
        full_Cartan_weight_dimension=n,rank_certifying_prime=PRIME,rank_pivots=rank.pivots,
        exact_actual_CAR_raising_kernels=True,exact_CAR_component_intertwining=True,
        exact_dimension_from_character_census=True,G=G,components=compressed)
    done.write_text(json.dumps(result,indent=2)+'\n')
    print('DONE',lam,flush=True)
    return result

def main():
    census=json.loads((HERE/'SPIN_FLAVOUR_CENSUS.json').read_text())
    blocks=[]
    for sector in census['sectors']:
        if sector['N']>12:continue
        # All competing irreps retained; N>12 follows explicitly declared CAR
        # particle-hole plus cell-reflection intertwiner on colour singlets.
        for r in sector['SU4_irreps']:blocks.append(build_block(r,sector['SU4_irreps']))
    (HERE/'SPIN_FLAVOUR_BLOCKS.json').write_text(json.dumps(dict(schema='oph.spin_flavour.exact_blocks.v1',
        full_CAR_dimension=1<<24,all_half_filling_and_lower_irreps_retained=True,blocks=blocks),indent=2)+'\n')
if __name__=='__main__':main()
