"""Independent exact review of metric, regional-force and capacity identities."""
from pathlib import Path
import hashlib
import json

import numpy as np
import sympy as s

ROOT=Path(__file__).resolve().parents[3]


def run():
    # Independent construction of the ten measurement rows, using SymPy's
    # exact inversion rather than the owning decoder implementation.
    pairs=tuple((i,j) for i in range(4) for j in range(i,4))
    vectors=[s.Matrix([1,0,0,0])]
    for i in range(1,4):
        for sign in (1,-1):
            v=s.zeros(4,1);v[0]=s.Rational(5,3);v[i]=sign*s.Rational(4,3)
            vectors.append(v)
    for i,j in ((1,2),(1,3),(2,3)):
        v=s.zeros(4,1);v[0]=3;v[i]=v[j]=2;vectors.append(v)
    assert all(-v[0]**2+sum(v[i]**2 for i in range(1,4))==-1 for v in vectors)
    design=s.Matrix([[v[i]*v[j]*(1 if i==j else 2) for i,j in pairs] for v in vectors])
    inverse=design.inv()
    norms=[sum(abs(x) for x in inverse.row(i)) for i in range(10)]
    assert max(norms)==s.Rational(17,8)
    sharp=inverse*s.Matrix([-1,1,1,0,0,0,0,0,0,0])
    assert sharp[pairs.index((1,1))]==s.Rational(17,8)
    assert sharp[pairs.index((1,1))]>2  # A reduced decoder constant fails.
    assert all(design[:i,:].col_join(design[i+1:,:]).rank()==9 for i in range(10))

    # Exact metric determinant and response ratios from the displayed density.
    r=s.symbols('r',positive=True)
    density=s.diag(-1/r,1/r,r*r,r*r)
    metric=s.simplify(s.sqrt(-density.det())*density.inv())
    assert metric==s.diag(-r*r,r*r,1/r,1/r)
    lapse=s.sqrt(-metric[0,0]);area=2*s.sqrt(metric[2,2]*metric[3,3])
    r0,r1=s.Rational(9,10),s.Rational(11,10)
    ratio=s.simplify((lapse.subs(r,r1)-lapse.subs(r,r0))/(area.subs(r,r1)-area.subs(r,r0)))
    assert ratio==-s.Rational(99,200)
    transverse=-metric[0,0]/metric[2,2]
    assert s.simplify(transverse.subs(r,r1)-transverse.subs(r,r0))==s.Rational(301,500)

    # Direct original source evolution, with no target response-module import.
    from model.native_equilibration.reciprocal_vacuum import finite_model
    source=finite_model()
    p=source.probabilities
    reset=np.r_[np.full(8,p[:8].sum()/8),np.full(8,p[8:].sum()/8)]
    updated=source.transition@reset
    dp=float(updated[8:].sum()-p[8:].sum())
    measured_lapse=np.repeat([.9,1.1],8)@(updated-reset)
    measured_area=np.repeat([20/9,20/11],8)@(updated-reset)
    assert dp>0 and abs(measured_lapse-float(ratio)*measured_area)<2e-15
    assert abs(measured_lapse-.2*dp)<2e-15

    # Independent constant-Hessian curvature specialization. Varying Phi
    # independently of Psi tests the spatial diagonal residual contribution.
    x,y,z,e=s.symbols('x y z e',real=True)
    coords=(s.symbols('t'),x,y,z)
    phi=e*(x*x+2*y*y+3*z*z)
    psi=e*(4*x*x+5*y*y+6*z*z)
    background=s.diag(-1,1,1,1)
    h=s.diag(-2*phi,-2*psi,-2*psi,-2*psi)
    # Compute Christoffel first; this differs from the target Ricci formula.
    gamma=[[[sum(background[a,b]*(s.diff(h[b,c],coords[d])+s.diff(h[b,d],coords[c])-s.diff(h[c,d],coords[b])) for b in range(4))/2
             for d in range(4)] for c in range(4)] for a in range(4)]
    ricci=s.Matrix(4,4,lambda b,d:sum(s.diff(gamma[a][b][d],coords[a])-s.diff(gamma[a][b][a],coords[d]) for a in range(4)))
    scalar=sum(background[a,a]*ricci[a,a] for a in range(4))
    einstein=ricci-background*scalar/2
    lap=lambda f:sum(s.diff(f,c,c) for c in (x,y,z))
    assert einstein[0,0]==2*lap(psi)
    assert sum(einstein[i,i] for i in range(1,4))==2*lap(phi-psi)
    assert s.simplify((einstein[0,0]+sum(einstein[i,i] for i in range(1,4)))/2)==lap(phi)
    epsilon,R,radius=s.symbols('epsilon R radius',positive=True)
    q=s.symbols('q',nonnegative=True)
    force=s.integrate(2*epsilon*q*q,(q,0,R))/radius**2
    assert force==2*epsilon*R**3/(3*radius**2)
    assert force==4*(epsilon*R**3/(6*radius**2))

    # Capacity laws derived from continuity, with exact independent fractions.
    a,k,N0,c=s.symbols('a k N0 c',positive=True)
    eta=s.symbols('eta',real=True)
    rho=k/(N0*a**eta)
    w=s.simplify(-1-a*s.diff(rho,a)/(3*rho))
    assert w==eta/3-1 and w.subs(eta,-1)==-s.Rational(4,3)
    assert s.simplify((c*k)/(c*N0*a**eta)-rho)==0
    histories=[]
    for exponent in (-2,0,3):
        for scale in (s.Rational(1,3),s.Rational(1),s.Rational(3)):
            n=scale**exponent
            om,orr,od=s.Rational(2,7),s.Rational(1,14),s.Rational(9,14)
            E2=om/scale**3+orr/scale**4+od/n
            recovered=od/(E2-om/scale**3-orr/scale**4)
            assert recovered==n
            histories.append({'exponent':exponent,'a':str(scale),'recovered_capacity_ratio':str(recovered)})
    paths=['model/openwave_gravity_capacity/'+f for f in ('PROOF.md','response.py','capacity.py','Transfer.lean','results.json')]
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    return {'schema':'oph.independent_gravity_capacity_review.v1',
            'reviewed_sha256':hashes,
            'reviewer_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'verdict':'PASS_SCOPED_MATHEMATICAL_CONSEQUENCES',
            'exact_decoder_norm':'17/8','deleted_probe_rank_checks':10,
            'lapse_area_response_ratio':str(ratio),
            'native_complete_thermal_update_lapse_diagnostic':float(measured_lapse),
            'independent_curvature_method':'Christoffel derivatives on independent quadratic Phi/Psi',
            'Poisson_Phi_to_Psi_error_ratio':'4',
            'capacity_histories':histories,
            'countercontrols':{'decoder_constant_two_fails':True,'missing_each_probe_loses_rank':True,
                'dropping_spatial_residuals_misses_factor_four':True,'decreasing_capacity_is_phantom':True,
                'rescaled_capacity_amplitude_cannot_be_identified_from_density':True},
            'scope':{'native_Newton_force_proved':False,'spatial_gravitational_redshift_proved':False,
                'cosmic_capacity_or_numeric_Lambda_selected':False,
                'Lean_kernel_rerun_by_this_audit':False}}


if __name__=='__main__':
    result=run()
    path=ROOT/'model/openwave_particle_readouts/gravity_review.json'
    path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['verdict'])
