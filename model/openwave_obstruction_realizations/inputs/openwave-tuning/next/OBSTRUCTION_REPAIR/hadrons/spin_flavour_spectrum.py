"""Exact inertia certificates for every competing spin/flavour source block."""
from pathlib import Path
from fractions import Fraction as F
from math import lcm,floor
import json
import numpy as np
import scipy.linalg as la
import audit
HERE=Path(__file__).resolve().parent

def integer_inertia(A):
    n=len(A);M=[list(map(int,r)) for r in A];original=[r[:] for r in M]
    prev=1;neg=0
    for k in range(n):
        pivot=next((i for i in range(k,n) if M[i][i]),None)
        if pivot is None:
            if all(not M[i][j] for i in range(k,n) for j in range(k,n)):return neg
            return audit.inertia([list(map(F,r)) for r in original])
        if pivot!=k:
            M[k],M[pivot]=M[pivot],M[k]
            for r in M:r[k],r[pivot]=r[pivot],r[k]
        d=M[k][k]
        neg+=int(d*prev<0)
        for i in range(k+1,n):
            for j in range(i,n):
                num=d*M[i][j]-M[i][k]*M[k][j]
                q,r=divmod(num,prev);assert r==0
                M[i][j]=M[j][i]=q
        prev=d
    return neg

def inertia(A,G,x):
    n=len(A);D=[[F(A[i][j])-x*G[i][j] for j in range(n)] for i in range(n)]
    den=lcm(*(v.denominator for row in D for v in row))
    return integer_inertia([[int(v*den) for v in row] for row in D])

def tests():
    rng=np.random.default_rng(62731)
    for n in range(1,9):
        for rep in range(8):
            A=rng.integers(-5,6,size=(n,n));A=A+A.T
            assert integer_inertia(A)==audit.inertia([list(map(F,r)) for r in A.tolist()])
    return 64

def numerical(A,G):
    g=np.array(G,dtype=float);a=np.array(A,dtype=float);s=1/np.sqrt(np.diag(g));g=s[:,None]*g*s[None,:];a=s[:,None]*a*s[None,:]
    condition=float(np.linalg.cond(g))
    if condition<1e8:return la.eigh(a,g,eigvals_only=True),dict(method='scaled float generalized eigensolver',Gram_condition_diagnostic=condition)
    import mpmath as mp
    with mp.workdps(90):
        gg=mp.matrix(G);aa=mp.matrix([[mp.mpf(v.numerator)/v.denominator for v in r] for r in A]);L=mp.cholesky(gg);Li=L**-1;H=Li*aa*Li.T
        values=la.eigvalsh(np.array(H.tolist(),dtype=float))
        return values,dict(method='90-digit Gram whitening and float spectral proposal',Gram_condition_diagnostic=condition)

def main():
    checks=tests();inp=json.loads((HERE/'RAW_COEFFICIENT_INPUTS.json').read_text());exactm=F(inp['432']['declared_light_midpoint']);m=F(round(exactm*(1<<36)),1<<36)
    error=24*abs(exactm-m);unit=F(inp['electron_anchor_GeV'])/exactm
    source=json.loads((HERE/'SPIN_FLAVOUR_BLOCKS.json').read_text());results=[]
    for kappa in [F(1,2),F(1),F(2)]:
        records=[]
        for block in source['blocks']:
            lam=block['highest_weight'];N=sum(lam);d=block['multiplicity'];G=block['G'];M=block['components']['M'];K=block['components']['K'];E=block['components']['E']
            A=[[m*M[i][j]+F(1,8)*K[i][j]+kappa*F(E[i][j],6) for j in range(d)] for i in range(d)]
            values,proposal=numerical(A,G)
            mult=int(sum(abs(v-values[0])<1e-7 for v in values))
            if N==0:lo=hi=computational_lo=computational_hi=F()
            else:
                den=1<<29;c=F(floor(float(values[0])*den),den);rad=F(2,den)
                for attempt in range(20):
                    lo,hi=c-rad,c+rad
                    below,above=inertia(A,G,lo),inertia(A,G,hi)
                    if below==0 and above==mult:break
                    rad*=2
                else:raise ArithmeticError(('ground enclosure not isolated',lam,kappa,below,above))
                computational_lo,computational_hi=lo,hi
                lo-=error;hi+=error
            record=dict(N=N,U4_highest_weight=lam,SU4_dimension=block['SU4_dimension'],multiplicity=d,
                lowest_energy_exact_interval=list(map(str,[lo,hi])),
                computational_inertia_interval_exact=list(map(str,[computational_lo,computational_hi])),certified_eigenvalue_count_in_ground_interval=mult,
                complete_spectrum_diagnostic=values.tolist(),floating_proposal=proposal,
                exact_inertia_lower=0 if N else None,exact_inertia_upper=mult if N else None,exact_scalar_zero=N==0)
            records.append(record)
            print('CERTIFIED',kappa,lam,float(lo),float(hi),flush=True)
        vacuum=min(records,key=lambda r:F(r['lowest_energy_exact_interval'][0]));vlo,vhi=map(F,vacuum['lowest_energy_exact_interval'])
        assert all(vhi<F(r['lowest_energy_exact_interval'][0]) for r in records if r is not vacuum)
        baryon=min([r for r in records if r['N']==9],key=lambda r:F(r['lowest_energy_exact_interval'][0]));blo,bhi=map(F,baryon['lowest_energy_exact_interval'])
        assert all(bhi<F(r['lowest_energy_exact_interval'][0]) for r in records if r['N']==9 and r is not baryon)
        gap=[blo-vhi,bhi-vlo]
        ground_by_N=[]
        for N in range(0,13,3):
            rows=[r for r in records if r['N']==N]
            lo=min(F(r['lowest_energy_exact_interval'][0]) for r in rows);hi=min(F(r['lowest_energy_exact_interval'][1]) for r in rows)
            ground_by_N.append(dict(N=N,exact_interval=list(map(str,[lo,hi])),possible_ground_irreps=[r['U4_highest_weight'] for r in rows if F(r['lowest_energy_exact_interval'][0])<=hi]))
        e2lo,e2hi=map(F,next(r['exact_interval'] for r in ground_by_N if r['N']==6))
        binding=[2*blo-e2hi-vhi,2*bhi-e2lo-vlo]
        results.append(dict(kappa=str(kappa),all_competing_irrep_records=records,
            vacuum_N=vacuum['N'],vacuum_U4_highest_weight=vacuum['U4_highest_weight'],vacuum_SU4_dimension=vacuum['SU4_dimension'],
            ground_abs_baryon_one_hole_U4_highest_weight=baryon['U4_highest_weight'],
            ground_abs_baryon_one_SU4_dimension=baryon['SU4_dimension'],
            Bminus1_and_Bplus1_equal_by_CAR_intertwiner=True,
            baryon_gap_exact_interval=list(map(str,gap)),conditional_baryon_gap_GeV_interval=[float(x*unit) for x in gap],
            complete_baryon_number_ground_records=ground_by_N,
            two_baryon_binding_sign_convention='2(E_B1-E_vac)-(E_B2-E_vac); positive means finite-volume binding',
            two_baryon_binding_exact_interval=list(map(str,binding)),
            conditional_two_baryon_binding_GeV_interval=[float(x*unit) for x in binding]))
    out=dict(schema='oph.obstruction_repair.spin_flavour_spectrum.v1',all_competing_SU4_irreps_certified=True,
        all_baryon_number_sectors_certified=True,source_wide_vacuum_certified=True,
        exact_inertia_self_checks=checks,source_mass=str(exactm),computational_mass=str(m),mass_operator_error_bound=str(error),
        t='1/8',common_energy_unit_GeV=float(unit),results=results,
        nucleon_delta_splitting=False,physical_mass_match_claim=False,
        common_spatial_Dirac_transfer=False,native_complete_program_executed=False)
    (HERE/'SPIN_FLAVOUR_SPECTRUM.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in ['kappa','vacuum_N','vacuum_U4_highest_weight','ground_abs_baryon_one_hole_U4_highest_weight','conditional_baryon_gap_GeV_interval']} for r in results],indent=2))
if __name__=='__main__':main()
