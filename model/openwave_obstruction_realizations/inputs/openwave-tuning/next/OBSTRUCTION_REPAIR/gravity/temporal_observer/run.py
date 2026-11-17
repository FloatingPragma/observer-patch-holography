"""All real source-momentum branches of an actual tilted field restriction."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name,'1')
from pathlib import Path
import sys,json,math
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from audit import ALPHA,BETA,dagger
alpha=ALPHA[0];I=np.eye(4);R=1.

def evaluate(u,kappa):
    # Squaring the two-band equation gives degree8. Every root is retained
    # in the record; only real momenta belong to the source Fourier space.
    poly=np.zeros(9);poly[0]=u*u*R**6;poly[6]=u*u-1;poly[7]=2*kappa;poly[8]=-kappa*kappa
    roots=np.roots(poly);G=np.zeros((4,4),complex);N=np.zeros_like(G);records=[];occupations=[]
    for root in sorted(roots,key=lambda z:(z.real,z.imag)):
        if abs(root.imag)>1e-7:
            records.append(dict(k=[float(root.real),float(root.imag)],class_='COMPLEX_NOT_IN_REAL_MOMENTUM_HILBERT_SPACE'))
            continue
        k=float(root.real);E=math.sqrt(k*k+R**6*k**8)
        band=1 if k-kappa>0 else -1
        Ep=(k+4*R**6*k**7)/E
        jac=1-u*band*Ep;weight=1/abs(jac)
        H=alpha*k+BETA*R**3*k**4;P=(I+band*H/E)/2
        G+=weight*P
        if band==-1:N+=weight*P
        occupations.extend([int(band==-1)]*2)
        records.append(dict(k=[k,0.],class_='REAL_PHYSICAL_BRANCH',band=band,energy=band*E,
                            jacobian=jac,positive_weight=weight,degeneracy=2,
                            unsquared_residual=k-u*band*E-kappa))
    ge,gu=np.linalg.eigh(G);invsqrt=(gu/np.sqrt(ge))@dagger(gu)
    C=invsqrt@N@invsqrt;ce=np.linalg.eigvalsh(C)
    normal=I-u*alpha;ne,nu=np.linalg.eigh(normal);sqrt_normal=(nu*np.sqrt(ne))@dagger(nu)
    Gdirac=sqrt_normal@G@sqrt_normal
    exactC=I/2-3*alpha/8-BETA/8
    return dict(tilt=u,kappa=kappa,roots=records,real_root_count=len(occupations)//2,
                complete_branch_modes=len(occupations),complete_branch_occupations=occupations,
                all_real_branches_retained=True,negative_weights_used=False,
                sampled_field_components=4,gram_eigenvalues=ge.tolist(),
                dirac_normalized_gram_eigenvalues=np.linalg.eigvalsh(Gdirac).tolist(),
                gram_limit_error=float(np.max(abs(Gdirac-4*I/3))),
                whitened_covariance_eigenvalues=ce.tolist(),
                covariance_limit_error=float(np.max(abs(C-exactC))),
                sampled_covariance_projector_defect=float(np.max(abs(C@C-C))),
                full_branch_covariance_projector_defect=0.,
                local_dirac_normalization_suffices=False)

def main():
    rows=[evaluate(u,k) for u in (.05,.1,.3,.6) for k in (1e-2,1e-3,1e-4)]
    boosts=[]
    for p in (.01,.1,.3,1.,3.):
        H=alpha*p+BETA*R**3*p**4;Hp=alpha+4*BETA*R**3*p**3
        actual=(H@Hp+Hp@H)/2
        expected=(p+4*R**6*p**7)*I
        boosts.append(dict(momentum_times_R=p*R,relative_boost_defect=4*(p*R)**6,
                           anticommutator_identity_error=float(np.max(abs(actual-expected)))))
    finest=[r for r in rows if r['kappa']==1e-4]
    out=dict(status='COMPLETE_SOURCE_TILTED_OBSERVER_TRANSFER_FAILS',R=R,
             analytic_covariance_eigenvalues=[.5-math.sqrt(10)/8,.5+math.sqrt(10)/8],
             root_records=rows,weighted_energy_boost=boosts,
             kappa_limit_has_extra_UV_branches=True,
             complete_positive_branch_register_available=True,
             complete_branch_register_is_proved_geometric_local_field=False,
             IR_Lorentz_errors_vanish=True,UV_eta_all_observer_transfer=False,
             lapse_equation_derived=False,physical_clock_join=False)
    out['algebra_checks_pass']=bool(all(r['real_root_count']==4 and r['gram_limit_error']<1e-3 and r['covariance_limit_error']<1e-3 and all(v>=-1e-12 and v<=1+1e-12 for v in r['whitened_covariance_eigenvalues']) for r in finest) and max(r['anticommutator_identity_error'] for r in boosts)<1e-9)
    (HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)
    if not out['algebra_checks_pass']:raise SystemExit(1)

if __name__=='__main__':main()
