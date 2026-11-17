"""Next even length of the same staggered SU3 vacuum source.

All colour weights, all physical singlets and all baryon-number sectors are
retained. Rational nullspaces and source intertwiners establish compression.
Floating spectra are only interval proposals and probe diagnostics.
"""
from fractions import Fraction as F
from itertools import combinations,product
from math import lcm,floor
from pathlib import Path
import json
import numpy as np
import scipy.linalg as la
import sympy as sy
import audit
from staggered_vacuum import creation,annihilation
HERE=Path(__file__).resolve().parent


def one_body(word,sites,a,b):
    out={}
    for x in sites:
        w,s=annihilation(3*x+b,word)
        if w is None:continue
        z,t=creation(3*x+a,w)
        if z is not None:out[z]=out.get(z,0)+s*t
    return {w:s for w,s in out.items() if s}


def colour6(words,sites):
    ix={w:i for i,w in enumerate(words)}
    out=np.zeros((len(words),len(words)),dtype=np.int64)
    mask=sum(7<<(3*x) for x in sites)
    for j,w in enumerate(words):
        out[j,j]-=(w&mask).bit_count()**2
        for a,b in product(range(3),repeat=2):
            for mid,s in one_body(w,sites,b,a).items():
                for z,t in one_body(mid,sites,a,b).items():
                    out[ix[z],j]+=3*s*t
    assert np.array_equal(out,out.T)
    return out


def sector(L,k):
    choices=list(combinations(range(L),k))
    words=sorted(sum(1<<(3*x+c) for c,sites in enumerate(sets) for x in sites)
                 for sets in product(choices,repeat=3))
    ix={w:i for i,w in enumerate(words)}
    C=colour6(words,list(range(L)))
    null=sy.Matrix(C.tolist()).nullspace()
    cols=[]
    for v in null:
        den=lcm(*(int(x.q) for x in v))
        cols.append([int(x*den) for x in v])
    R=np.array(cols,dtype=np.int64).T
    assert R.shape[0]==len(words) and not np.any(C@R)
    # nullspace is exact, exhaustive and computed anew rather than assumed.
    G=R.T@R
    assert audit.inertia([list(map(F,-row)) for row in G])==len(cols)
    M=np.diag([sum((-1)**x*((w>>(3*x))&7).bit_count() for x in range(L)) for w in words])
    K=np.zeros_like(C)
    for j,w in enumerate(words):
        for x in range(L-1):
            for c in range(3):
                for i,h in [(3*x+c,3*(x+1)+c),(3*(x+1)+c,3*x+c)]:
                    v,s=annihilation(h,w)
                    if v is None:continue
                    z,t=creation(i,v)
                    if z is not None:K[ix[z],j]-=s*t
    E=sum((colour6(words,list(range(cut))) for cut in range(1,L)),np.zeros_like(C))
    # Every compressed component leaves the complete singlet kernel invariant.
    for component in [M,K,E]:assert not np.any(C@(component@R))
    return dict(L=L,k=k,N=3*k,words=words,R=R,G=G,
                components=[R.T@X@R for X in [M,K,E]],
                weight_dimension=len(words),singlet_dimension=R.shape[1])


def inertia_at(A,G,x):
    return audit.inertia([[F(A[i][j])-x*int(G[i,j]) for j in range(len(G))] for i in range(len(G))])


def certify_lowest(A,G,proposal):
    denominator=1<<30
    center=F(floor(float(proposal)*denominator),denominator)
    radius=F(2,denominator)
    for _ in range(12):
        lo,hi=center-radius,center+radius
        if inertia_at(A,G,lo)==0 and inertia_at(A,G,hi)==1:return lo,hi
        radius*=2
    raise ArithmeticError('ground interval was not isolated')


def main():
    inp=json.loads((HERE/'RAW_COEFFICIENT_INPUTS.json').read_text())
    mass=F(inp['432']['declared_light_midpoint'])
    m=F(round(mass*(1<<36)),1<<36)
    error=12*abs(m-mass)
    sectors=[]
    for k in range(5):
        print('constructing',k,flush=True)
        s=sector(4,k);sectors.append(s)
        print('N',s['N'],'weight',s['weight_dimension'],'singlet',s['singlet_dimension'],flush=True)
    assert [s['singlet_dimension'] for s in sectors]==[1,20,50,20,1]
    results=[]
    for kappa in [F(1,2),F(1),F(2)]:
        numerical=[];certificates=[]
        for s in sectors:
            M,K,E=s['components'];d=s['singlet_dimension'];G=s['G']
            A=[[m*int(M[i,j])+F(1,8)*int(K[i,j])+kappa*F(int(E[i,j]),6) for j in range(d)] for i in range(d)]
            values,vectors=la.eigh(np.array(A,dtype=float),G.astype(float))
            if s['N'] in [0,12]:
                assert A==[[F()]];interval=(F(),F())
            else:
                interval=certify_lowest(A,G,values[0]);interval=(interval[0]-error,interval[1]+error)
            numerical.append((values,vectors))
            certificates.append(dict(N=s['N'],singlet_dimension=d,ground_interval_exact=list(map(str,interval)),
                complete_spectrum_diagnostic=values.tolist()))
            print('kappa',kappa,'N',s['N'],'ground',list(map(float,interval)),flush=True)
        vacuum=certificates[2]
        vlo,vhi=map(F,vacuum['ground_interval_exact'])
        assert all(vhi<F(r['ground_interval_exact'][0]) for j,r in enumerate(certificates) if j!=2)
        gaps=[]
        for j in [1,3]:
            blo,bhi=map(F,certificates[j]['ground_interval_exact'])
            gaps.append((blo-vhi,bhi-vlo))
        # The same local epsilon creator/remover, now from the actual N6 vacuum.
        sv=sectors[2];psi=sv['R']@numerical[2][1][:,0]
        residues=[]
        for sign,j in [(-1,1),(1,3)]:
            target=sectors[j];target_ix={w:i for i,w in enumerate(target['words'])}
            out=np.zeros(len(target['words']))
            for w,amplitude in zip(sv['words'],psi):
                z=w;s=1
                for mode in ([2,1,0] if sign==1 else [0,1,2]):
                    z,t=(creation if sign==1 else annihilation)(mode,z)
                    if z is None:break
                    s*=t
                if z is not None:out[target_ix[z]]+=s*amplitude
            frames=target['R']@numerical[j][1]
            overlaps=frames.T@out
            weights=abs(overlaps)**2
            assert abs(weights.sum()-np.dot(out,out))<1e-10
            residues.append(dict(baryon_change=sign,local_probe_success=float(weights.sum()),
                all_spectral_probabilities_diagnostic=weights.tolist(),ground_line_probability_diagnostic=float(weights[0]),
                failed_probe_probability=float(1-weights.sum())))
        unit=F(inp['electron_anchor_GeV'])/mass
        results.append(dict(kappa=str(kappa),complete_sector_records=certificates,
            ground_vacuum_N=6,baryon_gap_intervals_source=[list(map(str,g)) for g in gaps],
            conditional_GeV_intervals=[list(map(float,[g[0]*unit,g[1]*unit])) for g in gaps],
            complete_local_probe_records=residues))
    out=dict(schema='oph.obstruction_repair.staggered_chain4.v1',L=4,full_CAR_dimension=4096,
        all_baryon_sectors_retained=True,colour_zero_exterior_flux=True,
        actual_complete_singlet_dimensions=[s['singlet_dimension'] for s in sectors],
        actual_mass_exact=str(mass),rational_computational_mass=str(m),uniform_mass_operator_error=str(error),
        declared_t='1/8',declared_kappas=['1/2','1','2'],new_empirical_fit=False,
        results=results,physical_nucleon_identification=False,common_spatial_Dirac_transfer=False,
        native_complete_program_executed=False)
    (HERE/'STAGGERED_CHAIN4_RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in ['kappa','conditional_GeV_intervals']} for r in results],indent=2))


if __name__=='__main__':main()
