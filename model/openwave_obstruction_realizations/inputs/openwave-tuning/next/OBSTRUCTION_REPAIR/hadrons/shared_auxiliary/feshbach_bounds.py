"""Nonperturbative bounds on the actual full625-auxiliary shared-vacuum source.

P0 is both auxiliaries in vacuum, hence the complete frozen interacting source.
Q0 retains one and two auxiliary excitations. Every square root is enclosed by
exact rational arithmetic. The dressed trial is a literal full-source vector.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import json
import numpy as np
import scipy.linalg as la
import sympy as sy
import mpmath as mp
HERE=Path(__file__).resolve().parent;OLD=HERE.parent

def sqrt_bounds(x,bits=100):
    assert x>=0
    d=1<<bits;n=isqrt(x.numerator*d*d//x.denominator);lo=F(n,d)
    return lo,lo if lo*lo==x else F(n+1,d)

def lower_pair(p,q,b2):return (p+q-sqrt_bounds((p-q)**2+4*b2)[1])/2

def upper_pair(p,q,b2):return (p+q-sqrt_bounds((p-q)**2+4*b2)[0])/2

def frac(x):return F(int(x.p),int(x.q))
def mp_matrix(A):return mp.matrix([[mp.mpf(str(v.p))/mp.mpf(str(v.q)) for v in row] for row in A.tolist()])

def vacuum_trials(m):
    b=json.loads((OLD/'spin_flavour_cache/block_3_3_3_3.json').read_text());G=sy.Matrix(b['G']);M,K,E=[sy.Matrix(b['components'][k]) for k in ['M','K','E']]
    Gi=G.inv(method='DM')
    # Between globalSU4 singlets, exact SU4 averaging of the mixed current.
    D=sy.Rational(32,5)*G-sy.Rational(2,45)*M*Gi*M+E/180
    assert D==D.T
    trials=[]
    for kappa in [F(1,2),F(1),F(2)]:
        H=sy.Rational(m.numerator,m.denominator)*M+K/8+sy.Rational(kappa.numerator,kappa.denominator)*E/6
        with mp.workdps(90):
            L=mp.cholesky(mp_matrix(G));Li=L**-1;h=Li*mp_matrix(H)*Li.T
            _,vec=la.eigh(np.array(h.tolist(),dtype=float));proposal=Li.T*mp.matrix(vec[:,0].tolist());scale=max(abs(v) for v in proposal)
            c=sy.Matrix([int(mp.nint(v/scale*(1<<64))) for v in proposal])
        norm=(c.T*G*c)[0];p=frac((c.T*H*c)[0]/norm);d=frac((c.T*D*c)[0]/norm)
        kinetic=frac((c.T*K*c)[0]/(8*norm));cross=frac((c.T*D*Gi*H*c)[0]/norm)
        assert d>0
        Delta=3*kappa;g2=F(1,64);q=Delta+(cross-kinetic)/d;b2=g2*d
        bound=upper_pair(p,q,b2)
        trials.append(dict(kappa=str(kappa),integer_source_coordinates=list(map(int,c)),exact_source_norm=str(frac(norm)),
            exact_bare_Rayleigh=str(p),exact_mixed_current_expectation=str(d),exact_hopping_expectation=str(kinetic),
            exact_G2_H_expectation=str(cross),exact_auxiliary_trial_diagonal=str(q),exact_offdiagonal_squared=str(b2),
            exact_full_source_variational_upper=str(bound),variational_upper_diagnostic=float(bound),
            total_spin=0,total_isospin=0,full625_source_trial=True))
        print('EXACT DRESSED VACUUM TRIAL',kappa,float(bound),'G²',float(d),flush=True)
    return trials

def main():
    local=json.loads((HERE/'LOCAL_ALGEBRA.json').read_text());norms={r['N']:r for r in local['norms']}
    old=json.loads((OLD/'SPIN_FLAVOUR_SPECTRUM.json').read_text());inp=json.loads((OLD/'RAW_COEFFICIENT_INPUTS.json').read_text());m=F(inp['432']['declared_light_midpoint']);unit=F(inp['electron_anchor_GeV'])/m
    omega=sqrt_bounds(F(1,64)+m*m);trials=vacuum_trials(m);out=[]
    for base,trial in zip(old['results'],trials):
        kappa=F(base['kappa']);Delta=3*kappa;records=[]
        for N in range(0,25,3):
            sourceN=min(N,24-N);ground=next(r for r in base['complete_baryon_number_ground_records'] if r['N']==sourceN)
            plo,phi=map(F,ground['exact_interval']);info=norms[N];a2=F(info['P_zero_to_one_aux_norm_squared']);b2=F(info['Q_one_to_two_aux_norm_squared_bound'])
            free_count=min(N,24-N)
            # At half filling the entire negative free band is a unique colour
            # singlet; a quark adjoint needs at least one particle-hole pair.
            one_count=free_count-2 if N==12 else free_count
            d1=Delta-one_count*omega[1];d2=2*Delta-free_count*omega[1]
            qlo=lower_pair(d1,d2,b2)
            lower=lower_pair(plo,qlo,a2)
            upper=F(trial['exact_full_source_variational_upper']) if N==12 else phi
            assert lower<=upper
            records.append(dict(N=N,P0_ground_interval_exact=list(map(str,[plo,phi])),
                Q0_one_aux_lower=str(d1),Q0_two_aux_lower=str(d2),Q0_lower_exact=str(qlo),
                P0_Q0_coupling_norm_squared_bound=str(a2),Q0_one_two_coupling_norm_squared_bound=str(b2),
                full_source_ground_interval_exact=list(map(str,[lower,upper])),
                full_source_ground_interval_diagnostic=[float(lower),float(upper)],
                half_filling_quark_adjoint_free_band_exclusion_used=N==12))
        vacuum=next(r for r in records if r['N']==12);vlo,vhi=map(F,vacuum['full_source_ground_interval_exact'])
        Nstable=all(F(r['full_source_ground_interval_exact'][0])>vhi for r in records if r['N']!=12)
        # Any nonzero total spin or isospin P0 state is orthogonal to every
        # SU4 singlet, so all other irreps furnish a conservative lower bound.
        competing=min(F(r['lowest_energy_exact_interval'][0]) for r in base['all_competing_irrep_records'] if r['N']==12 and r['U4_highest_weight']!=[3,3,3,3])
        otherSI=lower_pair(competing,F(vacuum['Q0_lower_exact']),F(vacuum['P0_Q0_coupling_norm_squared_bound']))
        spin_certified=otherSI>vhi
        gaps=[]
        for N in [9,15]:
            blo,bhi=map(F,next(r['full_source_ground_interval_exact'] for r in records if r['N']==N));gap=[blo-vhi,bhi-vlo]
            gaps.append(dict(N=N,exact_gap_interval=list(map(str,gap)),conditional_GeV_interval=[float(x*unit) for x in gap]))
        # H0 electric/aux strong-coupling expansion: conservative full-source
        # norm bounds. A finite bound may still be far too broad to identify a
        # ground band or physical mass; no promotion follows from its existence.
        v=12*omega[1]+2*sqrt_bounds(F(17,192))[1] #2g sqrt(max_localG²=17/3)
        gamma=F(4,3)*kappa;bPQ2=F(1,64)*(20+F(34,3));energy_abs=max(abs(vlo),abs(vhi));den=gamma-v-energy_abs
        remainder=None if den<=0 else bPQ2*(v+energy_abs)/(gamma*den)
        out.append(dict(kappa=str(kappa),all_N_sector_bounds=records,
            N12_vacuum_sector_certified=Nstable,spin0_isospin0_vacuum_certified=spin_certified,
            other_spin_or_isospin_lower_exact=str(otherSI),vacuum_trial=trial,
            baryon_gap_bounds=gaps,H_eff2_uniform_operator_remainder=dict(H0_gap=str(gamma),perturbation_norm_upper=str(v),
                offdiagonal_norm_squared_upper=str(bPQ2),energy_absolute_upper=str(energy_abs),resolvent_margin=str(den),
                bound_exists=den>0,bound_exact=None if remainder is None else str(remainder),
                source_units_diagnostic=None if remainder is None else float(remainder),
                usable_as_physical_mass_accuracy=False)))
        print('FULL625 BOUNDS',kappa,'N12',Nstable,'spin0I0',spin_certified,'gaps',gaps,'remainder',None if remainder is None else float(remainder),flush=True)
    result=dict(schema='oph.shared_auxiliary.full625_feshbach.v1',
        full_auxiliary_carrier=625,full_CAR_dimension=1<<24,full_total_carrier=(1<<24)*625,
        exact_local_algebra_and_coupling_norms=True,all_N_sectors_retained=True,
        P0='both auxiliary vacua; entire frozen interacting quark source',
        Q0='all one- and two-auxiliary states satisfying total Gauss law',
        full_Q0_resolvent_replaced_by_bounds_not_by_Delta_scalar=True,
        new_empirical_fit=False,t='1/8',g='1/8',Delta='3*kappa',mass=str(m),results=out,
        interacting_spin_resolved_baryon_spectrum_computed=False,
        physical_mass_solution=False,native_complete_program_executed=False)
    (HERE/'FESHBACH_BOUNDS.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
