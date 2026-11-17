"""Nonperturbative full-source spin bounds through a complete low SU4 block.

The 20-dimensional subspace is invariant for H_old, not for the auxiliary
coupling. Its entire mixed-current complement is retained by an exact D²
covariance. Positivity is proved by exact integer inertia, not floating roots.
"""
from pathlib import Path
from fractions import Fraction as F
from math import floor,lcm
import sys,json
import numpy as np
import scipy.linalg as la
import sympy as sy
import mpmath as mp
HERE=Path(__file__).resolve().parent;OLD=HERE.parent
sys.path.insert(0,str(OLD));sys.path.insert(0,str(OLD/'shared_auxiliary'))
from spin_flavour_spectrum import integer_inertia
from feshbach_bounds import frac,mp_matrix,upper_pair

def rational(x):return sy.Rational(x.numerator,x.denominator)
def pd(A):
    rows=[[frac(x) for x in row] for row in A.tolist()];den=lcm(*(x.denominator for row in rows for x in row));n=len(rows)
    count=integer_inertia([[-int(x*den) for x in row] for row in rows]);return count==n,count

def whiten(G,*matrices):
    with mp.workdps(90):
        L=mp.cholesky(mp_matrix(G));Li=L**-1
        return [np.array((Li*mp_matrix(A)*Li.T).tolist(),float) for A in matrices]

def trial(G,H,K,D,Gi,kappa):
    with mp.workdps(90):
        L=mp.cholesky(mp_matrix(G));Li=L**-1;h=Li*mp_matrix(H)*Li.T
        _,v=la.eigh(np.array(h.tolist(),float));proposal=Li.T*mp.matrix(v[:,0].tolist());scale=max(abs(x) for x in proposal)
        c=sy.Matrix([int(mp.nint(x/scale*(1<<64))) for x in proposal])
    norm=(c.T*G*c)[0];p=frac((c.T*H*c)[0]/norm);d=frac((c.T*D*c)[0]/norm)
    kin=frac((c.T*K*c)[0]/(8*norm));cross=frac((c.T*D*Gi*H*c)[0]/norm)
    q=3*kappa+(cross-kin)/d;b2=d/64;upper=upper_pair(p,q,b2)
    return dict(integer_source_coordinates=list(map(int,c)),source_norm_exact=str(frac(norm)),
        bare_Rayleigh_exact=str(p),Gsum_expectation_exact=str(d),hopping_expectation_exact=str(kin),
        Gsum_H_expectation_exact=str(cross),one_auxiliary_diagonal_exact=str(q),coupling_squared_exact=str(b2),
        full_source_upper_exact=str(upper),full_source_upper_diagnostic=float(upper),
        literal_full_source_trial=True,full_source_eigenvector_claim=False)

def main():
    states=json.loads((HERE/'BARYON_STATES.json').read_text())['records'];branch=json.loads((HERE/'BRANCHING.json').read_text())['results'][0]['source_SU4_branching']
    qbounds=json.loads((HERE/'Q_COLOUR_BOUNDS.json').read_text())['results'];old=json.loads((OLD/'SPIN_FLAVOUR_SPECTRUM.json').read_text())['results'];inp=json.loads((OLD/'RAW_COEFFICIENT_INPUTS.json').read_text());m=F(inp['432']['declared_light_midpoint']);outputs=[]
    cached={}
    for state in states:
        G=sy.Matrix(state['G']);Gi=G.inv(method='DM');M,K,E=[sy.Matrix(state['components'][t]) for t in ['M','K','E']]
        D=sy.Matrix(state['mixed_current24'])/24;D2=sy.Matrix(state['mixed_current24_square'])/576
        V=D2-D*Gi*D;assert V==V.T
        # V is exactly the Gram form of (1-P_low) D R by construction.
        cached[state['name']]=(G,Gi,M,K,E,D,D2,V)
    for qb,ob in zip(qbounds,old):
        kappa=F(qb['kappa']);records=[]
        for state in states:
            name=state['name'];s,i=state['twice_spin'],state['twice_isospin'];G,Gi,M,K,E6,D,D2,V=cached[name]
            H=rational(m)*M+K/8+rational(kappa)*E6/6
            qr=next(r for r in qb['records'] if (r['N'],r['twice_spin'],r['twice_isospin'])==(9,s,i))
            # Downward rounding only weakens the Q bound and controls certificate sizes.
            qraw=F(qr['Q_lower_exact']);q=F(floor(qraw*(1<<36)),1<<36)
            dmax=F(qr['P0_max_Gsum_squared']);competitors=[]
            for rep in branch:
                lam=rep['U4_highest_weight']
                if lam==[3,3,3,0]:continue
                if not any((r['twice_spin'],r['twice_isospin'])==(s,i) for r in rep['spin_isospin_branches']):continue
                r=next(r for r in ob['all_competing_irrep_records'] if r['U4_highest_weight']==lam)
                competitors.append(dict(U4_highest_weight=lam,ground_lower_exact=r['lowest_energy_exact_interval'][0]))
            ehigh=min(F(r['ground_lower_exact']) for r in competitors)
            tr=trial(G,H,K,D,Gi,kappa)
            hn,dn,vn=whiten(G,H,D,V);eye=np.eye(len(G.tolist()))
            def numerical_min(x):
                c=1/(64*(float(q)-x));h=float(ehigh)-x-c*float(dmax)
                if h<=0:return -1e100
                return la.eigvalsh(hn-x*eye-c*dn-c*c/h*vn,subset_by_index=[0,0])[0]
            lo=float(F(qr['refined_full_source_lower_exact']))-.01;hi=float(F(tr['full_source_upper_exact']))
            assert numerical_min(lo)>0
            for _ in range(55):
                mid=(lo+hi)/2
                if numerical_min(mid)>0:lo=mid
                else:hi=mid
            # A floating root only proposes a rational endpoint; all accepted
            # endpoints are independently checked using exact positivity.
            candidate=F(floor(lo*(1<<30))-4,1<<30);attempts=[]
            for attempt in range(20):
                c=F(1,64)/(q-candidate);h=ehigh-candidate-c*dmax
                assert q>candidate and h>0
                B=H-rational(candidate)*G-rational(c)*D-rational(c*c/h)*V
                passed,count=pd(B);attempts.append(dict(energy=str(candidate),negative_inertia_of_minus_matrix=count))
                if passed:break
                candidate-=F(1<<attempt,1<<30)
            assert passed
            records.append(dict(name=name,N=9,twice_spin=s,twice_isospin=i,
                low_SU4_highest_weight=[3,3,3,0],low_multiplicity_dimension=20,
                all_other_SU4_competitors=competitors,other_P0_ground_lower_exact=str(ehigh),
                Q_lower_exact_before_rounding=str(qraw),Q_lower_exact=str(q),D_operator_norm_upper_exact=str(dmax),
                current_covariance='D2-D G^-1 D; full CAR current images, not compressed-square replacement',
                full_source_lower_exact=str(candidate),full_source_lower_diagnostic=float(candidate),
                second_Feshbach_margin_exact=str(h),coefficient_c_exact=str(c),
                exact_positive_definite_dimension=20,exact_inertia_receipts=attempts,
                trial=tr))
            print('FULL SOURCE MATRIX BOUND',kappa,name,float(candidate),float(F(tr['full_source_upper_exact'])),flush=True)
        nuc=next(r for r in records if r['name']=='nucleon');nupper=F(nuc['trial']['full_source_upper_exact']);allcomp=[]
        for row in qb['records']:
            if row['N']!=9 or (row['twice_spin'],row['twice_isospin'])==(1,1):continue
            found=next((r for r in records if (r['twice_spin'],r['twice_isospin'])==(row['twice_spin'],row['twice_isospin'])),None)
            lower=F(found['full_source_lower_exact']) if found else F(row['refined_full_source_lower_exact'])
            allcomp.append(dict(twice_spin=row['twice_spin'],twice_isospin=row['twice_isospin'],lower_exact=str(lower),
                above_nucleon_trial=lower>nupper))
        outputs.append(dict(kappa=str(kappa),records=records,all_competing_P0_present_baryon_spin_sectors=allcomp,
            Q_only_sectors_positive_by_prior_bounds=True,
            full_source_baryon_spin_half_isospin_half_certified=all(r['above_nucleon_trial'] for r in allcomp),
            source_wide_vacuum_sector_from_frozen_packet=True))
        print('ACTUAL BARYON SPIN CERTIFIED',kappa,outputs[-1]['full_source_baryon_spin_half_isospin_half_certified'],flush=True)
    (HERE/'MATRIX_FESHBACH.json').write_text(json.dumps(dict(schema='oph.auxiliary_spin_refinement.matrix_feshbach.v1',
        full625_source=True,all_current_complements_retained=True,effective_eigenvector_substitution=False,
        exact_PD_proof='Hlow-EG-cD-(c²/h)(D2-DG^-1D)>0; c=g²/(q-E); h=e_high-E-c||D||',
        physical_mass_solution=False,common_spatial_Dirac_transfer=False,results=outputs),indent=2)+'\n')
if __name__=='__main__':main()
