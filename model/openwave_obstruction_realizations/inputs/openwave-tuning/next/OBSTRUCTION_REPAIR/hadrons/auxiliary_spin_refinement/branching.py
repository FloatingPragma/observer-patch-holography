"""Exact colour/spin/isospin counts and free Q-sector gates for full625 source."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import product,permutations
from functools import lru_cache
import json
HERE=Path(__file__).resolve().parent;OLD=HERE.parent

@lru_cache(None)
def weights(top):
    if len(top)==1:return Counter({top:1})
    out=Counter()
    for row in product(*(range(top[i+1],top[i]+1) for i in range(len(top)-1))):
        for w,n in weights(row).items():out[w+(sum(top)-sum(row),)]+=n
    return out

def dim(top):
    z=F(1)
    for i in range(len(top)):
        for j in range(i+1,len(top)):z*=F(top[i]-top[j]+j-i,j-i)
    return int(z)

def su2_branch(top):
    char=Counter()
    for w,n in weights(tuple(top)).items():char[w[0]+w[1]-w[2]-w[3],w[0]+w[2]-w[1]-w[3]]+=n
    rows=[];parity=sum(top)%2
    for s,i in product(range(parity,25,2),repeat=2):
        m=char[s,i]-char[s+2,i]-char[s,i+2]+char[s+2,i+2]
        assert m>=0
        if m:rows.append(dict(twice_spin=s,twice_isospin=i,multiplicity=m))
    assert sum((r['twice_spin']+1)*(r['twice_isospin']+1)*r['multiplicity'] for r in rows)==dim(top)
    return rows

@lru_cache(None)
def colour_tensor(a,b,p,q):
    total=sum(a)+sum(b)
    if (total-p-2*q)%3:return 0
    d=(total-p-2*q)//3
    if d<0:return 0
    target=(p+q+d,q+d,d);char=Counter()
    for x,n in weights(a).items():
        for y,m in weights(b).items():char[tuple(u+v for u,v in zip(x,y))]+=n*m
    result=0;rho=(2,1,0)
    for per in permutations(range(3)):
        sign=(-1)**sum(per[i]>per[j] for i in range(3) for j in range(i+1,3))
        result+=sign*char[tuple(target[i]+rho[i]-rho[per[i]] for i in range(3))]
    assert result>=0
    return result

def tri(a,b,c):return abs(a-b)<=c<=a+b and (a+b-c)%2==0

def main():
    assert [colour_tensor((2,1,0),(2,1,0),p,q) for p,q in [(0,0),(1,1),(3,0),(0,3),(2,2)]]==[1,2,1,1,1]
    loc=json.loads((OLD/'shared_auxiliary/LOCAL_ALGEBRA.json').read_text())['full_joint_irreps'];original=json.loads((OLD/'SPIN_FLAVOUR_SPECTRUM.json').read_text());out=[]
    for N in [9,12]:
        sectors=Counter();norm={}
        for nl in range(max(0,N-12),min(N,12)+1):
            for a,b in product(loc[nl]['irreps'],loc[N-nl]['irreps']):
                if a['colour_Dynkin']!=list(reversed(b['colour_Dynkin'])):continue
                ca=tuple(a['colour_highest_weight']);cb=tuple(b['colour_highest_weight']);assert colour_tensor(ca,cb,0,0)==1
                for s in range(abs(a['twice_spin']-b['twice_spin']),a['twice_spin']+b['twice_spin']+1,2):
                    for i in range(abs(a['twice_isospin']-b['twice_isospin']),a['twice_isospin']+b['twice_isospin']+1,2):
                        sectors[s,i]+=a['multiplicity']*b['multiplicity'];value=F(a['mixed_current_square'])+F(b['mixed_current_square']);norm[s,i]=max(norm.get((s,i),F()),value)
        expected=14112 if N==9 else 24696
        assert sum((s+1)*(i+1)*m for (s,i),m in sectors.items())==expected
        print('EXACT GLOBAL COLOUR/SPIN/ISO CENSUS',N,'states',expected,'largest block',max(sectors.values()),flush=True)
        source_branches=[]
        for r in original['results'][0]['all_competing_irrep_records']:
            if r['N']!=N:continue
            source_branches.append(dict(U4_highest_weight=r['U4_highest_weight'],spin_isospin_branches=su2_branch(r['U4_highest_weight'])))
        q_records=[]
        for (s,i),m in sorted(sectors.items()):
            free={}
            for auxnum in [1,2]:
                targets=[(1,1)] if auxnum==1 else [(0,0),(1,1),(3,0),(0,3),(2,2)];auxspin=[2] if auxnum==1 else [0,2,4];winner=None
                for positive in range(max(0,N-12),min(N,12)+1):
                    negative=N-positive;found=False
                    for a,b in product(loc[negative]['irreps'],loc[positive]['irreps']):
                        if not tri(a['twice_isospin'],b['twice_isospin'],i):continue
                        if not any(tri(sq,sa,s) for sq in range(abs(a['twice_spin']-b['twice_spin']),a['twice_spin']+b['twice_spin']+1,2) for sa in auxspin):continue
                        if not any(colour_tensor(tuple(a['colour_highest_weight']),tuple(b['colour_highest_weight']),p,q) for p,q in targets):continue
                        found=True;break
                    if found:winner=positive;break
                free[str(auxnum)]=dict(sector_present=winner is not None,minimum_positive_band_occupation=winner,exact_free_energy_coefficient=None if winner is None else 2*winner-N)
            q_records.append(dict(twice_spin=s,twice_isospin=i,complete_P0_highest_weight_dimension=m,
                exact_P0_max_Gsum_squared=str(norm[s,i]),exact_P0_Q0_coupling_norm_squared=str(norm[s,i]/64),Q0_free_band_bounds=free))
        out.append(dict(N=N,complete_P0_colour_singlet_dimension=expected,largest_spin_isospin_block=max(sectors.values()),
            source_SU4_branching=source_branches,spin_isospin_sectors=q_records))
    (HERE/'BRANCHING.json').write_text(json.dumps(dict(schema='oph.auxiliary_spin_refinement.branching.v1',results=out,
        all_colour_spin_isospin_multiplets_retained=True,free_band_Q_bounds_are_for_full625_source=True),indent=2)+'\n')
    print('ALL BRANCHING AND Q FREE BAND GATES PASSED',flush=True)
if __name__=='__main__':main()
