"""Full-source two-level current block with exact H² leakage and Feshbach bounds."""
from pathlib import Path
from fractions import Fraction as F
from math import floor
import sys,json
import sympy as sy
import numpy as np
import scipy.linalg as la
HERE=Path(__file__).resolve().parent;OLD=HERE.parent;PREV=OLD/'auxiliary_spin_refinement'
sys.path.insert(0,str(OLD));sys.path.insert(0,str(PREV));sys.path.insert(0,str(OLD/'shared_auxiliary'))
from matrix_feshbach import rational,pd,whiten
from feshbach_bounds import lower_pair

def main():
    states=json.loads((PREV/'BARYON_STATES.json').read_text())['records'];mom=json.loads((HERE/'MOMENTS.json').read_text())['records'];previous=json.loads((PREV/'MATRIX_FESHBACH.json').read_text())['results'];inp=json.loads((OLD/'RAW_COEFFICIENT_INPUTS.json').read_text());m=F(inp['432']['declared_light_midpoint']);out=[]
    data={}
    for s,mo in zip(states,mom):
        assert s['name']==mo['name'];G=sy.Matrix(s['G']);Gi=G.inv(method='DM');D=sy.Matrix(s['mixed_current24'])/24;Di=D.inv(method='DM');D2=sy.Matrix(s['mixed_current24_square'])/576
        M,K,E=[sy.Matrix(s['components'][x]) for x in ['M','K','E']];T=sy.Matrix(mo['CdaggerC1536'])/1536;L=sy.Matrix(mo['Dleft_Dright576'])/576
        data[s['name']]=(G,Gi,D,Di,D2,M,K,E,T,L)
    for base in previous:
        kappa=F(base['kappa']);records=[]
        for prior in base['records']:
            name=prior['name'];G,Gi,D,Di,D2,M,K,E,T,L=data[name];Delta=3*kappa;H=rational(m)*M+K/8+rational(kappa)*E/6;J=H+rational(Delta)*G
            Z=(H*Gi*D-D*Gi*H)/2-K/8
            Fone=rational(Delta)*D+(H*Gi*D+D*Gi*H)/2-K/8
            W=J*Gi*D*Gi*J+J*Gi*Z+Z.T*Gi*J+T
            H2one=D2/64+W+L/16
            V=H2one-D*Gi*D/64-Fone*Di*Fone
            assert Fone==Fone.T and V==V.T
            Gram=sy.diag(G,D);B=H.row_join(D/8).col_join((D/8).row_join(Fone));Leak=sy.diag(sy.zeros(20),V)
            # The 40-vector frame consists of the entire old low block and its
            # complete one-auxiliary current image, not selected current outputs.
            q=F(prior['Q_lower_exact']);ehigh=F(prior['other_P0_ground_lower_exact']);dmax=F(prior['D_operator_norm_upper_exact']);alpha=lower_pair(ehigh,q,dmax/64)
            # Strictly downward rational rounding of the complement lower bound.
            alpha=F(floor(alpha*(1<<36)),1<<36)
            Vp=(D2-D*Gi*D)/64;Vq=V-Vp
            LeakP=sy.diag(sy.zeros(20),Vp);LeakQ=sy.diag(sy.zeros(20),Vq)
            bn,vn,vpn,vqn=whiten(Gram,B,Leak,LeakP,LeakQ);vdiag=la.eigvalsh(vn)
            assert min(vdiag[0],la.eigvalsh(vpn)[0],la.eigvalsh(vqn)[0])>-1e-8
            eigen=la.eigvalsh(bn);den=1<<29;upper=F(floor(float(eigen[0])*den)+4,den)
            # Gram D is rational, so use a direct exact negative count.
            from math import lcm
            from spin_flavour_spectrum import integer_inertia
            def count_below(x):
                A=B-rational(x)*Gram;vals=[[F(int(z.p),int(z.q)) for z in row] for row in A.tolist()];scale=lcm(*(z.denominator for row in vals for z in row));return integer_inertia([[int(z*scale) for z in row] for row in vals])
            count=count_below(upper);assert count>=1
            def numerical_min(x):
                if x>=float(alpha):return -1e100
                return la.eigvalsh(bn-x*np.eye(40)-vn/(float(alpha)-x),subset_by_index=[0,0])[0]
            lo=min(-2.,float(F(prior['full_source_lower_exact']))-.1);hi=min(float(upper),float(alpha)-1e-9)
            while numerical_min(lo)<=0:lo*=2
            for _ in range(55):
                mid=(lo+hi)/2
                if numerical_min(mid)>0:lo=mid
                else:hi=mid
            candidate=F(floor(lo*(1<<29))-4,1<<29);attempts=[]
            for attempt in range(20):
                assert candidate<alpha
                S=B-rational(candidate)*Gram-Leak/rational(alpha-candidate)
                passed,neg=pd(S);attempts.append(dict(energy=str(candidate),negative_inertia_of_minus_matrix=neg))
                if passed:break
                candidate-=F(1<<attempt,1<<29)
            assert passed
            # A weighted two-sector complement keeps the current leakage into
            # other quark irreps separate from its leakage into Q. Young's
            # inequality gives exact positive block-diagonal lower margins.
            a2=dmax/64
            def weighted_min(x,r):
                if x>=float(alpha):return -1e100
                edge=float(a2)/(float(q)-x);window=float(ehigh)-x-edge
                if window<=0:return -1e100
                eta=edge+r*window;hp=float(ehigh)-x-eta;hq=float(q)-x-float(a2)/eta
                if min(hp,hq)<=0:return -1e100
                return la.eigvalsh(bn-x*np.eye(40)-vpn/hp-vqn/hq,subset_by_index=[0,0])[0]
            proposals=[]
            for rr in range(1,16):
                r=F(rr,16);wl=min(-2.,float(candidate)-.1);wh=min(float(upper),float(alpha)-1e-9)
                while weighted_min(wl,float(r))<=0:wl*=2
                for _ in range(55):
                    mid=(wl+wh)/2
                    if weighted_min(mid,float(r))>0:wl=mid
                    else:wh=mid
                proposals.append((wl,r))
            proposal,r=max(proposals,key=lambda z:z[0]);weighted=F(floor(proposal*(1<<29))-4,1<<29);weighted_attempts=[]
            for attempt in range(20):
                edge=a2/(q-weighted);window=ehigh-weighted-edge;eta=edge+r*window
                hp=ehigh-weighted-eta;hq=q-weighted-a2/eta;assert min(hp,hq)>0
                test=B-rational(weighted)*Gram-LeakP/rational(hp)-LeakQ/rational(hq)
                good,nneg=pd(test);weighted_attempts.append(dict(energy=str(weighted),negative_inertia_of_minus_matrix=nneg))
                if good:break
                weighted-=F(1<<attempt,1<<29)
            assert good
            record=dict(name=name,N=9,twice_spin=prior['twice_spin'],twice_isospin=prior['twice_isospin'],
                actual_full_source_frame_dimension=40,complement_lower_exact=str(alpha),
                full_source_lower_exact=str(candidate),full_source_lower_diagnostic=float(candidate),
                full_source_variational_upper_exact=str(upper),full_source_variational_upper_diagnostic=float(upper),
                exact_PD_dimension=40,exact_inertia_receipts=attempts,compressed_upper_inertia=count,
                weighted_complement_bound=dict(full_source_lower_exact=str(weighted),full_source_lower_diagnostic=float(weighted),
                    proof_allocation_r_exact=str(r),eta_exact=str(eta),P_complement_margin_exact=str(hp),Q_complement_margin_exact=str(hq),
                    exact_PD_dimension=40,exact_inertia_receipts=weighted_attempts,
                    all_proof_allocation_trials=[dict(r=str(rr),floating_lower_proposal=x) for x,rr in proposals]),
                strongest_full_source_lower_exact=str(max(candidate,weighted,F(prior['full_source_lower_exact']))),
                current_block_ritz_spectrum_diagnostic=eigen.tolist(),
                current_covariance_eigenvalues_diagnostic=vdiag.tolist(),
                all_two_auxiliary_outputs_retained=True,
                old_certified_lower_exact=prior['full_source_lower_exact'],old_trial_upper_exact=prior['trial']['full_source_upper_exact'],
                new_lower_stronger=candidate>F(prior['full_source_lower_exact']),new_upper_stronger=upper<F(prior['trial']['full_source_upper_exact']))
            records.append(record)
            print('FULL SOURCE CURRENT BLOCK',kappa,name,'bounds',float(candidate),float(upper),'old',float(F(prior['full_source_lower_exact'])),float(F(prior['trial']['full_source_upper_exact'])),'alpha',float(alpha),'weighted',float(weighted),flush=True)
        out.append(dict(kappa=str(kappa),records=records))
    (HERE/'CURRENT_BLOCK.json').write_text(json.dumps(dict(schema='oph.auxiliary_current_block.full_source.v1',source_unchanged=True,
        P='40-dimensional frame R plus A R; all20 low source multiplicities and every current output',
        covariance='R40†H²R40-H40 G40^-1 H40; exactly zero in first20 rows and columns',
        complement='all other old-source irreps plus Q orthogonal to A R; bounded by the prior full-Q bound',
        no_effective_eigenvector_substitution=True,no_mass_fit=True,results=out,physical_rows_solved=False),indent=2)+'\n')
if __name__=='__main__':main()
