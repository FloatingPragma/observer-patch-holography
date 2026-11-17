"""Complete8464-state second-order source with actual spin/isospin currents.

This is the derived effective operator. Its spectrum is not promoted to the
full625-auxiliary source without the separately reported remainder bound.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import product
from math import lcm,floor
import sys,json
import numpy as np
import scipy.linalg as la
import sympy as sy
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from spin_flavour_spectrum import inertia
CACHE=HERE/'effective_blocks';CACHE.mkdir(exist_ok=True)

def local(r):
    d=r['dimension'];G=np.array(r['G'],dtype=np.int64);z=[];w=[]
    for key,out in [('twice_Sz',z),('twice_Iz',w)]:
        Z=np.array(r[key],dtype=np.int64)
        out.extend(int(Z[j,j]//G[j,j]) for j in range(d))
        assert np.array_equal(Z,G*np.array(out)[None,:])
    return dict(N=r['N'],d=d,G=G,z=z,w=w,Sp=np.array(r['Splus'],dtype=np.int64),Ip=np.array(r['Iplus'],dtype=np.int64),
        G24=np.array(r['G24'],dtype=np.int64),F={tuple(map(int,k.split(','))):np.array(v,dtype=np.int64) for k,v in r['F'].items()})

def sub_tensor(A,B,rows,cols):return np.array([[int(A[a,c])*int(B[b,d]) for c,d in cols] for a,b in rows],dtype=np.int64)

def blocks(left,right):
    filename=CACHE/f"N{left['N']}_{right['N']}.json"
    if filename.exists():return json.loads(filename.read_text())
    pairs=list(product(range(left['d']),range(right['d'])));weight={(a,b):(left['z'][a]+right['z'][b],left['w'][a]+right['w'][b]) for a,b in pairs};count=Counter(weight.values());N=left['N']+right['N'];out=[]
    for s,i in product(range(N%2,13,2),repeat=2):
        expected=count[s,i]-count[s+2,i]-count[s,i+2]+count[s+2,i+2]
        assert expected>=0
        if not expected:continue
        cols=[p for p in pairs if weight[p]==(s,i)];sr=[p for p in pairs if weight[p]==(s+2,i)];ir=[p for p in pairs if weight[p]==(s,i+2)]
        SP=sub_tensor(left['Sp'],right['G'],sr,cols)+sub_tensor(left['G'],right['Sp'],sr,cols)
        IP=sub_tensor(left['Ip'],right['G'],ir,cols)+sub_tensor(left['G'],right['Ip'],ir,cols)
        # Empty codomains are retained as exact0-by-dim constraints.
        SP=SP.reshape(len(sr),len(cols));IP=IP.reshape(len(ir),len(cols));constraints=np.vstack([SP,IP]);constraints=constraints[np.any(constraints,axis=1)]
        null=sy.Matrix(constraints.tolist()).nullspace() if len(constraints) else sy.eye(len(cols)).columnspace()
        assert len(null)==expected
        Q=[]
        for v in null:
            den=lcm(*(int(x.q) for x in v));Q.append([int(x*den) for x in v])
        Q=np.array(Q,dtype=np.int64).T
        assert not np.any(SP@Q) and not np.any(IP@Q)
        G=sub_tensor(left['G'],right['G'],cols,cols)
        X=sum((sub_tensor(left['F'][a,b],right['F'][b,a],cols,cols) for a,b in product(range(4),repeat=2)),np.zeros_like(G))
        J=sub_tensor(left['G24'],right['G'],cols,cols)+sub_tensor(left['G'],right['G24'],cols,cols)
        for A in [G,X,J]:assert np.array_equal(A,A.T)
        out.append(dict(N_left=left['N'],N_right=right['N'],twice_spin=s,twice_isospin=i,
            highest_weight_dimension=expected,full_sector_dimension=(s+1)*(i+1)*expected,
            weight_basis_dimension=len(cols),G=(Q.T@G@Q).tolist(),exchange=(Q.T@X@Q).tolist(),mixed_current24=(Q.T@J@Q).tolist()))
    assert sum(r['full_sector_dimension'] for r in out)==left['d']*right['d']
    filename.write_text(json.dumps(out,indent=2)+'\n');print('EXACT EFFECTIVE BLOCKS',left['N'],right['N'],len(out),flush=True);return out

def main():
    data=json.loads((HERE/'LOCAL_OPERATORS.json').read_text())['sectors'];loc=list(map(local,data));allblocks=[]
    for left,right in product(loc,repeat=2):allblocks.extend(blocks(left,right))
    assert sum(r['full_sector_dimension'] for r in allblocks)==8464
    inp=json.loads((HERE.parent/'RAW_COEFFICIENT_INPUTS.json').read_text());m0=F(inp['432']['declared_light_midpoint']);m=F(round(m0*(1<<36)),1<<36);error=24*abs(m-m0);results=[]
    for kappa in [F(1,2),F(1),F(2)]:
        records=[]
        for b in allblocks:
            d=b['highest_weight_dimension'];G=np.array(b['G'],dtype=np.int64);X=b['exchange'];J=b['mixed_current24'];N=b['N_left']+b['N_right']
            A=[[(m*(b['N_left']-b['N_right'])-F(3,256)*N/kappa)*int(G[i,j])+F(1,128)*F(X[i][j])/kappa-F(1,4608)*F(J[i][j])/kappa for j in range(d)] for i in range(d)]
            gf=G.astype(float);sc=1/np.sqrt(np.diag(gf));values=la.eigh(sc[:,None]*np.array(A,float)*sc[None,:],sc[:,None]*gf*sc[None,:],eigvals_only=True);mult=int(sum(abs(x-values[0])<1e-8 for x in values))
            center=F(floor(values[0]*(1<<32)),1<<32);rad=F(2,1<<32)
            for rep in range(20):
                lo,hi=center-rad,center+rad
                if inertia(A,G.tolist(),lo)==0 and inertia(A,G.tolist(),hi)==mult:break
                rad*=2
            else:raise ArithmeticError('effective interval failed')
            records.append(dict(N=N,N_left=b['N_left'],N_right=b['N_right'],twice_spin=b['twice_spin'],twice_isospin=b['twice_isospin'],
                highest_weight_dimension=d,computational_inertia_interval=list(map(str,[lo,hi])),
                true_mass_effective_interval=list(map(str,[lo-error,hi+error])),ground_band_eigenvalue_count=mult,
                all_highest_weight_spectrum_diagnostic=values.tolist()))
        perN=[]
        for N in range(0,25,3):
            rs=[r for r in records if r['N']==N];lo=min(F(r['true_mass_effective_interval'][0]) for r in rs);hi=min(F(r['true_mass_effective_interval'][1]) for r in rs)
            candidates=[{k:r[k] for k in ['N_left','N_right','twice_spin','twice_isospin']} for r in rs if F(r['true_mass_effective_interval'][0])<=hi]
            perN.append(dict(N=N,exact_effective_ground_interval=list(map(str,[lo,hi])),possible_ground_sectors=candidates))
        vlo=min(F(r['exact_effective_ground_interval'][0]) for r in perN);vhi=min(F(r['exact_effective_ground_interval'][1]) for r in perN)
        vacuum=[r for r in perN if F(r['exact_effective_ground_interval'][0])<=vhi]
        results.append(dict(kappa=str(kappa),all_effective_sector_records=records,all_N_ground_records=perN,possible_effective_vacua=vacuum,
            effective_vacuum_interval=list(map(str,[vlo,vhi])),source_remainder_must_be_applied_before_physical_comparison=True))
        print('EFFECTIVE VACUUM',kappa,vacuum,flush=True)
    out=dict(schema='oph.shared_auxiliary.complete_effective_second_order.v1',complete_P_dimension=8464,
        all_spin_isospin_and_occupation_sectors_retained=True,total_highest_weight_blocks=len(allblocks),
        largest_highest_weight_dimension=max(b['highest_weight_dimension'] for b in allblocks),
        exact_source_derived_coefficients=True,t='1/8',g='1/8',full_source_equality_claim=False,
        mass_rounding_operator_bound=str(error),results=results,physical_mass_comparison_performed=False,
        H_eff2_remainder_is_separate_and_may_be_too_broad=True)
    (HERE/'EFFECTIVE_SECOND_ORDER.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
