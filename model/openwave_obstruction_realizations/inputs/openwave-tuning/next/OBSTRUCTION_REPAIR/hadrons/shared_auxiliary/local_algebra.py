"""Exact full12-CAR colour/spin/isospin identity and complete local singlets."""
from pathlib import Path
from itertools import product,permutations
from collections import defaultdict
from fractions import Fraction as F
from math import prod,comb,lcm
import sys,json
import numpy as np
from scipy import sparse
import sympy as sy
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from staggered_vacuum import creation,annihilation
from staggered_chain4 import sector


def current(pairs):
    rr=[];cc=[];dd=[]
    for w in range(4096):
        col={}
        for a,b,c in pairs:
            z,s=annihilation(b,w)
            if z is None:continue
            y,t=creation(a,z)
            if y is not None:col[y]=col.get(y,0)+c*s*t
        for y,v in col.items():
            if v:rr.append(y);cc.append(w);dd.append(v)
    return sparse.csr_matrix((np.array(dd,dtype=np.int64),(rr,cc)),shape=(4096,4096))

def shifts(n):
    rho=list(range(n-1,-1,-1));out=[]
    for p in permutations(range(n)):
        sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        out.append((tuple(rho[i]-rho[p[i]] for i in range(n)),sign))
    return out

def dim(lam):return int(prod(F(lam[i]-lam[j]+j-i,j-i) for i in range(len(lam)) for j in range(i+1,len(lam))))
def ucas(lam):return sum(a*(a+len(lam)+1-2*(i+1)) for i,a in enumerate(lam))

def census():
    weights=defaultdict(int)
    for w in range(4096):
        c=[0]*3;s=[0]*2;i=[0]*2
        for j in range(12):
            if w>>j&1:c[j%3]+=1;s[(j//3)//2]+=1;i[(j//3)%2]+=1
        weights[tuple(c+s+i)]+=1
    rows=[]
    for N in range(13):
        irreps=[]
        for cl in product(range(5),repeat=3):
            if sum(cl)!=N or list(cl)!=sorted(cl,reverse=True):continue
            for sl,il in product(range((N+1)//2,min(N,6)+1),repeat=2):
                sp=(sl,N-sl);ip=(il,N-il);mu=0
                for dc,sc in shifts(3):
                    for ds,ss in shifts(2):
                        for di,si in shifts(2):mu+=sc*ss*si*weights.get(tuple(a+b for a,b in zip(cl,dc))+tuple(a+b for a,b in zip(sp,ds))+tuple(a+b for a,b in zip(ip,di)),0)
                assert mu>=0
                if not mu:continue
                S=F(2*sl-N,2);I=F(2*il-N,2);Cc=(F(ucas(cl))-F(N*N,3))/2
                gg=F(N*(12-N),6)-I*(I+1)/2-Cc/4-S*(S+1)/6
                assert gg>=0
                irreps.append(dict(colour_highest_weight=list(cl),colour_Dynkin=[cl[0]-cl[1],cl[1]-cl[2]],
                    twice_spin=int(2*S),twice_isospin=int(2*I),multiplicity=mu,
                    full_dimension=dim(cl)*int(2*S+1)*int(2*I+1)*mu,
                    colour_Casimir=str(Cc),mixed_current_square=str(gg)))
        assert sum(r['full_dimension'] for r in irreps)==comb(12,N)
        rows.append(dict(N=N,full_CAR_dimension=comb(12,N),irreps=irreps,
                         exact_max_G_squared=str(max(F(r['mixed_current_square']) for r in irreps))))
    return rows


def main():
    N=np.array([w.bit_count() for w in range(4096)],dtype=np.int64)
    sp=current([(3*f+c,3*(f+2)+c,1) for f in range(2) for c in range(3)])
    ip=current([(3*a+c,3*(a+1)+c,1) for a in [0,2] for c in range(3)])
    zs=sp@sp.T-sp.T@sp;zi=ip@ip.T-ip.T@ip
    S4=zs@zs+2*(sp.T@sp+sp@sp.T);I4=zi@zi+2*(ip.T@ip+ip@ip.T)
    C6=-sparse.diags(N*N,format='csr',dtype=np.int64);G24=-S4
    for a,b in product(range(3),repeat=2):
        E=current([(3*alpha+a,3*alpha+b,1) for alpha in range(4)])
        Z=current([(3*alpha+a,3*alpha+b,1 if alpha<2 else -1) for alpha in range(4)])
        P=current([(3*f+a,3*(f+2)+b,1) for f in range(2)])
        C6+=3*(E.T@E);G24+=3*(Z.T@Z)+6*(P.T@P+P@P.T)
    rhs=sparse.diags(4*N*(12-N),format='csr',dtype=np.int64)-3*I4-C6-S4
    diff=(G24-rhs).tocsr();diff.eliminate_zeros();assert diff.nnz==0
    print('FULL4096 CAR IDENTITY PASSED',flush=True)
    records=census();operators=[];total=0
    Fs={(a,b):current([(3*a+c,3*b+c,1) for c in range(3)]) for a,b in product(range(4),repeat=2)}
    for k in range(5):
        sec=sector(4,k);words=sec['words'];R=sec['R'];n=3*k;d=R.shape[1];total+=d
        small=lambda A:A[words,:][:,words]
        assert not np.any(small(C6)@R)
        direct=R.T@(small(G24)@R)
        assert np.array_equal(direct,R.T@((sparse.eye(len(words),dtype=np.int64)*4*n*(12-n)-3*small(I4)-small(S4))@R))
        branches=[]
        for s,i in product(range(n%2,min(n,12-n)+1,2),repeat=2):
            constraints=np.vstack([(small(zs)-s*sparse.eye(len(words),dtype=np.int64))@R,
                                   (small(zi)-i*sparse.eye(len(words),dtype=np.int64))@R,
                                   small(sp)@R,small(ip)@R])
            constraints=constraints[np.any(constraints,axis=1)]
            null=sy.Matrix(constraints.tolist()).nullspace() if len(constraints) else sy.eye(d).columnspace()
            if not null:continue
            columns=[]
            for v in null:
                den=lcm(*(int(x.q) for x in v));columns.append([int(x*den) for x in v])
            Q=R@np.array(columns,dtype=np.int64).T
            eigen=(F(n*(12-n))-3*F(i*(i+2),4)-F(s*(s+2),4))/6
            assert np.array_equal(small(G24)@Q,int(24*eigen)*Q)
            branches.append(dict(twice_spin=s,twice_isospin=i,multiplicity=Q.shape[1],full_dimension=(s+1)*(i+1)*Q.shape[1],mixed_current_square=str(eigen)))
        assert sum(r['full_dimension'] for r in branches)==d
        operators.append(dict(N=n,dimension=d,words=words,R=R.tolist(),G=(R.T@R).tolist(),
            F={f'{a},{b}':(R.T@(small(Fs[a,b])@R)).tolist() for a,b in Fs},
            S4=(R.T@(small(S4)@R)).tolist(),I4=(R.T@(small(I4)@R)).tolist(),G24=direct.tolist(),
            twice_Sz=(R.T@(small(zs)@R)).tolist(),twice_Iz=(R.T@(small(zi)@R)).tolist(),
            Splus=(R.T@(small(sp)@R)).tolist(),Iplus=(R.T@(small(ip)@R)).tolist(),branches=branches))
        print('LOCAL SINGLET',n,d,branches,flush=True)
    assert total==92
    norms=[];g=F(1,8)
    for n in range(0,25,3):
        full=[];physical=[]
        for l in range(max(0,n-12),min(12,n)+1):
            r=n-l
            full.append(max(F(x['mixed_current_square']) for x in records[l]['irreps'])+max(F(x['mixed_current_square']) for x in records[r]['irreps']))
            for a,b in product(records[l]['irreps'],records[r]['irreps']):
                if a['colour_Dynkin']==list(reversed(b['colour_Dynkin'])):physical.append(F(a['mixed_current_square'])+F(b['mixed_current_square']))
        norms.append(dict(N=n,exact_full_quark_max_sum_G_squared=str(max(full)),
            exact_global_colour_singlet_max_sum_G_squared=str(max(physical)),
            Q_one_to_two_aux_norm_squared_bound=str(g*g*max(full)),P_zero_to_one_aux_norm_squared=str(g*g*max(physical))))
    out=dict(schema='oph.shared_auxiliary.local_algebra.v1',full4096_identity_exact=True,
        identity='24G²=4N(12-N)-12I²-6C_colour-4S²',
        complete_joint_character_dimension=4096,full_joint_irreps=records,
        local_colour_singlet_dimension=92,complete_two_cell_local_singlet_dimension=92**2,
        all_local_singlet_spin_isospin_branches_verified=True,norms=norms,
        physical_model_solved=False)
    (HERE/'LOCAL_ALGEBRA.json').write_text(json.dumps(out,indent=2)+'\n')
    (HERE/'LOCAL_OPERATORS.json').write_text(json.dumps(dict(sectors=operators),indent=2)+'\n')
    print('EXACT NORMS',norms,flush=True)
if __name__=='__main__':main()
