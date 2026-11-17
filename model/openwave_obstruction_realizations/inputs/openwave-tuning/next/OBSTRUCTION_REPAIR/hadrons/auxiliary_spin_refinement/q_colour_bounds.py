"""Sharper full-Q bounds using actual colour-electric zero sectors.

No effective Hamiltonian is substituted. A second projection bounds the
fixed-auxiliary-number quark Hamiltonian itself by its exact hopping image.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import product,permutations
from functools import lru_cache
import sys,json
from branching import weights,dim,su2_branch,tri
HERE=Path(__file__).resolve().parent;OLD=HERE.parent;sys.path.insert(0,str(OLD/'shared_auxiliary'))
from feshbach_bounds import sqrt_bounds,lower_pair

def ucas(lam):return sum(a*(a+5-2*(j+1)) for j,a in enumerate(lam))

def local_flavour(N,aux):
    if N%3:return None
    k=N//3
    if aux:
        if k not in [1,2,3]:return None
        colour=(k+1,k,k-1)
    else:
        if k not in range(5):return None
        colour=(k,k,k)
    return tuple(sum(a>=j+1 for a in colour) for j in range(4))

@lru_cache(None)
def tensor(a,b):
    total=sum(a)+sum(b);char=Counter()
    for x,n in weights(a).items():
        for y,m in weights(b).items():char[tuple(u+v for u,v in zip(x,y))]+=n*m
    out=[];rho=(3,2,1,0)
    for top in product(range(7),repeat=4):
        if sum(top)!=total or tuple(sorted(top,reverse=True))!=top:continue
        mult=0
        for p in permutations(range(4)):
            sign=(-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))
            mult+=sign*char[tuple(top[i]+rho[i]-rho[p[i]] for i in range(4))]
        assert mult>=0
        if mult:out.append((top,mult,su2_branch(top)))
    assert sum(dim(top)*m for top,m,_ in out)==dim(a)*dim(b)
    return tuple(out)

@lru_cache(None)
def electric_P_bounds(N,auxnum,s,i):
    possibilities=[]
    for al in [0,1]:
        ar=auxnum-al
        if ar not in [0,1]:continue
        for nl in range(max(0,N-12),min(12,N)+1):
            nr=N-nl;a=local_flavour(nl,al);b=local_flavour(nr,ar)
            if a is None or b is None:continue
            auxspins=[0] if auxnum==0 else [2] if auxnum==1 else [0,2,4]
            vals=[]
            for top,mult,branches in tensor(a,b):
                if not any(r['twice_isospin']==i and any(tri(r['twice_spin'],sa,s) for sa in auxspins) for r in branches):continue
                value=F(N)-F(ucas(top)-ucas(a)-ucas(b),3);assert value>=0
                vals.append(value)
            if vals:possibilities.append(dict(N_left=nl,N_right=nr,aux_left=al,aux_right=ar,
                exact_max_P_K_squared=str(max(vals)),quark_flavour_left=list(a),quark_flavour_right=list(b)))
    if not possibilities:return None
    return dict(exact_max_P_K_squared=str(max(F(x['exact_max_P_K_squared']) for x in possibilities)),
                minimum_mass_number=min(x['N_left']-x['N_right'] for x in possibilities),all_allowed_local_occupancies=possibilities)

def main():
    groups=json.loads((HERE/'BRANCHING.json').read_text())['results'];scalars=json.loads((HERE/'SCALAR_SPIN_BOUNDS.json').read_text())['results'];norms={r['N']:r for r in json.loads((OLD/'shared_auxiliary/LOCAL_ALGEBRA.json').read_text())['norms']};inp=json.loads((OLD/'RAW_COEFFICIENT_INPUTS.json').read_text());m=F(inp['432']['declared_light_midpoint']);omega=sqrt_bounds(F(1,64)+m*m);out=[]
    for base in scalars:
        kappa=F(base['kappa']);gamma=F(4,3)*kappa;records=[]
        for group in groups:
            N=group['N']
            for sec in group['spin_isospin_sectors']:
                s,i=sec['twice_spin'],sec['twice_isospin'];diagonals=[];receipts=[]
                for auxnum in [1,2]:
                    free=sec['Q0_free_band_bounds'][str(auxnum)]
                    if not free['sector_present']:continue
                    coef=free['exact_free_energy_coefficient'];flo=coef*(omega[1] if coef<0 else omega[0]);p=electric_P_bounds(N,auxnum,s,i)
                    if p is None:quark=gamma+flo
                    else:quark=lower_pair(m*p['minimum_mass_number'],gamma+flo,F(p['exact_max_P_K_squared'])/64)
                    # The elementary free lower bound remains independently valid.
                    quark=max(quark,flo);diag=auxnum*3*kappa+quark;diagonals.append(diag)
                    receipts.append(dict(auxiliary_number=auxnum,electric_zero_sector=p,
                        free_quark_lower_exact=str(flo),electric_complement_lower_exact=str(gamma+flo),
                        rigorous_quark_electric_lower_exact=str(quark),fixed_auxiliary_diagonal_lower_exact=str(diag)))
                a2=F(sec['exact_P0_Q0_coupling_norm_squared']);qlo=None
                if len(diagonals)==1:qlo=diagonals[0]
                elif len(diagonals)==2:qlo=lower_pair(diagonals[0],diagonals[1],F(norms[N]['Q_one_to_two_aux_norm_squared_bound']))
                row=next(r for r in base['all_P0_present_spin_isospin_sectors'] if (r['N'],r['twice_spin'],r['twice_isospin'])==(N,s,i));plo,phi=map(F,row['P0_ground_interval_exact'])
                lower=plo if qlo is None else lower_pair(plo,qlo,a2)
                records.append(dict(N=N,twice_spin=s,twice_isospin=i,
                    all_fixed_auxiliary_electric_receipts=receipts,Q_lower_exact=None if qlo is None else str(qlo),
                    P0_Q0_norm_squared=str(a2),P0_max_Gsum_squared=sec['exact_P0_max_Gsum_squared'],
                    P0_ground_interval_exact=list(map(str,[plo,phi])),
                    refined_full_source_lower_exact=str(lower),refined_source_lower_diagnostic=float(lower)))
        out.append(dict(kappa=str(kappa),records=records))
        print('COLOUR-RESOLVED Q GATES',kappa,[{k:r[k] for k in ['N','twice_spin','twice_isospin','refined_source_lower_diagnostic']} for r in records if (r['N']==9 and r['twice_spin']==r['twice_isospin'] and r['twice_spin'] in [1,3]) or (r['N']==12 and r['twice_spin']==2 and r['twice_isospin']==0)],flush=True)
    (HERE/'Q_COLOUR_BOUNDS.json').write_text(json.dumps(dict(schema='oph.auxiliary_spin_refinement.q_colour.v1',
        actual_fixed_auxiliary_Hamiltonian_bounded=True,zero_electric_auxiliary_sectors_complete=True,
        local_projected_bilinear_identity='delta_colour/3 F between all total-colour singlets, including quark-adjoint plus auxiliary-adjoint',
        local_quark_U4_Casimir='7N-N²/3-6n_aux on a local total-colour singlet',
        hopping_square='N-(C_U4_total-C_U4_left-C_U4_right)/3',
        complete_U4_tensor_character_dimensions_verified=True,
        no_H_eff_eigenvalue_substitution=True,results=out),indent=2)+'\n')
if __name__=='__main__':main()
