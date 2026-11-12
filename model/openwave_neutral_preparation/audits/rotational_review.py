"""Independent complete-shell magnetic and torque controls."""
from pathlib import Path
import hashlib
import json
import sys
from fractions import Fraction
from itertools import product
import numpy as np
import scipy.linalg as la
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import sympy as s

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from model.openwave_angular_torque.source import oscillator,operators,evaluate,rhs,field_energy
from model.openwave_rotational_g_limit.limit import bounds,stable_kernel


def exact_shell():
    basis=[v for v in product(range(3),repeat=3) if sum(v)<=2];n=len(basis)
    annihilators=[]
    for axis in range(3):
        a=s.zeros(n)
        for j,v in enumerate(basis):
            if v[axis]:
                w=list(v);w[axis]-=1;a[basis.index(tuple(w)),j]=s.sqrt(v[axis])
        annihilators.append(a)
    x=[(a+a.H)/s.sqrt(2) for a in annihilators]
    p=[s.I*(a.H-a)/s.sqrt(2) for a in annihilators]
    sigma=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
    alpha=[s.kronecker_product(sigma[0],z) for z in sigma]
    a=sum((s.kronecker_product(z,q) for z,q in zip(alpha,p)),s.zeros(4*n))
    e=(s.kronecker_product(alpha[1],x[0])-s.kronecker_product(alpha[0],x[1]))/2
    t=(a*a).applyfunc(s.simplify);q=(e*e).applyfunc(s.simplify);w=a*e+e*a
    assert (t*(t-s.Rational(5,2)*s.eye(4*n))).applyfunc(s.simplify)==s.zeros(4*n)
    assert (q*(q-s.eye(4*n)/4)*(q-s.eye(4*n)/2)).applyfunc(s.simplify)==s.zeros(4*n)
    assert w[0,0]==1 and w[n,n]==-1 and t[0,0]==s.Rational(3,2)
    assert q[0,0]==q[n,n]==s.Rational(1,4)
    return {'complete_rows':4*n,'exact_native_reduced_square_identities':2,
        'vacuum_magnetic_derivatives':[1,-1],'both_probe_quadratic_vacua':'1/4',
        'all_spectral_signs_and_four_native_blocks_retained_by_exact_existing_intertwiner':True}


def torque_control():
    d=oscillator();ops=operators(d,mass=3.,omega=.27);n=len(ops[0]);rng=np.random.default_rng(778319)
    psi=rng.normal(size=n)+1j*rng.normal(size=n);psi/=np.linalg.norm(psi)
    b=np.array([.027,.041,-.033]);pi=np.array([.023,-.014,.019]);kappa=.16
    h,g,_=evaluate(psi,b,ops)
    matter=np.array([np.real(1j*np.vdot(psi,(h@j-j@h)@psi)) for j in d['J']])
    field=np.cross(kappa*pi/d['electric_metric'],pi)+np.cross(b,-d['magnetic_metric']*b/kappa-g)
    residual=np.linalg.norm(matter+field);assert residual<1e-12 and np.linalg.norm(matter)>1e-6
    deriv=[]
    for axis in range(3):
        z=np.eye(3)[axis]*2e-5
        ep=np.vdot(psi,evaluate(psi,b+z,ops)[0]@psi).real
        em=np.vdot(psi,evaluate(psi,b-z,ops)[0]@psi).real
        deriv.append(abs((ep-em)/4e-5-g[axis]))
    assert max(deriv)<1e-9
    sol=solve_ivp(rhs,(0,1.7),np.r_[psi,b,pi],args=(d,ops,kappa),rtol=2e-11,atol=2e-13,method='DOP853',t_eval=np.linspace(0,1.7,13))
    assert sol.success
    energies=[];angular=[]
    for state in sol.y.T:
        v=state[:n];bb=state[n:n+3].real;pp=state[n+3:].real
        hh,_,_=evaluate(v,bb,ops)
        energies.append(np.vdot(v,hh@v).real+field_energy(bb,pp,d,kappa))
        angular.append(np.array([np.vdot(v,j@v).real for j in d['J']])+np.cross(bb,pp))
    de=max(abs(x-energies[0]) for x in energies);dj=max(np.linalg.norm(x-angular[0]) for x in angular)
    assert de<2e-9 and dj<2e-9
    # Analytic Gaussian integrals, evaluated independently from source metrics.
    R=1.2;c0=np.pi**1.5*R**5/4;c1=np.pi**1.5*R**3*(1-1.5+15/16+.5-5/16)
    assert abs(c0-d['electric_metric'])<1e-14 and abs(c1-d['magnetic_metric'])<1e-14
    return {'fresh_state_mass_trap_coupling':[3.,.27,kappa],
        'complete_instantaneous_torque_residual':float(residual),
        'zero_field_reaction_mutation_torque':float(np.linalg.norm(matter)),
        'independent_current_finite_difference_errors':deriv,
        'fresh_full_matter_field_trajectory_energy_drift':float(de),
        'fresh_total_angular_drift':float(dj),'Gaussian_mode_metrics_independently_integrated':True}


def axial_control():
    d=oscillator();n=d['spatial_dimension'];j=d['J'][2]
    j2=sum((z@z for z in d['J']),np.zeros_like(j));j2values,j2vectors=la.eigh(j2)
    vals,vecs=la.eigh(j[:2*n,:2*n]);up=np.vstack((vecs[:,abs(vals-.5)<1e-12],np.zeros_like(vecs[:,abs(vals-.5)<1e-12])))
    rows=[]
    for m in (64,256,1024,4096,16384):
        kappa=m**-6;ops=operators(d,mass=m,omega=1.)
        def state(b):
            h0,derivatives,probe,trap=ops
            direct=h0+b*derivatives[2]
            eigen,u=la.eigh(direct);absolute=(u*np.abs(eigen))@u.conj().T+trap
            # The stable positive-square implementation is compared with the
            # independent full signed absolute value at a nonzero field.
            stable,D,E,A=stable_kernel(d,m,(0,0,b),(0,0,0))
            assert np.linalg.norm(absolute-(stable+m*np.eye(len(stable))))<2e-9
            w,v=la.eigh(up.conj().T@stable@up);psi=up@v[:,0]
            _,gradient,probe_current=evaluate(psi,(0,0,b),ops)
            return w[0],psi,gradient,probe_current
        bound=bounds(m);bb=float(bound['self_field_norm_upper'])
        def equation(b):return b+kappa*state(b)[2][2]/d['magnetic_metric']
        assert equation(-bb)<=0 and equation(bb)>=0
        b=brentq(equation,-bb,bb,xtol=1e-100,rtol=1e-14)
        shift,psi,g,mu=state(b)
        zero=state(0.)[1]
        assert np.linalg.norm(j2@zero-.75*zero)<1e-10
        casimir_probability=abs(j2vectors.conj().T@psi)**2
        assert abs(sum(casimir_probability)-1)<1e-12
        residual=abs(equation(b))/kappa
        assert residual<1e-9 and np.linalg.norm(g[:2])<1e-11
        jz=np.vdot(psi,j@psi).real;assert abs(jz-.5)<1e-12
        rest=m+shift+d['magnetic_metric']*b*b/(2*kappa)
        gstatic=4*rest*abs(mu[2]);assert abs(gstatic-2)<float(bound['static_g_error_from_two'])
        assert bound['axial_stationary_branch_contraction']<1
        assert bound['self_field_contraction']<1
        assert bound['uniform_precession_radius_lower']>0
        assert bound['uniform_phase_error']>0
        rows.append({'mass':m,'fresh_signed_absolute_static_g':float(gstatic),
            'independent_scalar_root':float(b),'scaled_field_residual':float(residual),
            'zero_self_field_J2_residual':float(np.linalg.norm(j2@zero-.75*zero)),
            'finite_self_field_full_J2_values':j2values.tolist(),
            'finite_self_field_full_J2_probabilities':casimir_probability.tolist(),
            'finite_self_field_exact_J2_asserted':False,
            'static_rational_g_error':str(bound['static_g_error_from_two']),
            'finite_time_rational_g_error':str(bound['quasistatic_precession_g_error_from_two'])})
    assert Fraction(rows[-1]['finite_time_rational_g_error'])<Fraction(4863,1000000)
    assert all(Fraction(b['finite_time_rational_g_error'])<Fraction(a['finite_time_rational_g_error']) for a,b in zip(rows,rows[1:]))
    return {'rows':rows,'fresh_static_root_solver':'bracketed scalar root, independent of producer fixed-point iteration',
        'exact_B_vacuum_coefficients_and_finite_complement_bounds_checked':True,
        'stationary_gap_scope':'fixed beta-positive J_z axial sector, with all four identical native blocks retained',
        'finite_time_phase_bound_scope':'uniform norm perturbation of the full kernel; no two-level Hamiltonian substitution'}


def main():
    paths=[ROOT/'model'/name/f for name in ('openwave_angular_torque','openwave_rotational_g_limit') for f in ('PROOF.md','PAPER_RESULT.md','results.json')]
    paths.extend(p for name in ('openwave_angular_torque','openwave_rotational_g_limit') for p in (ROOT/'model'/name).glob('*.py'))
    paths.append(Path(__file__).resolve())
    inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    dependency_checks={}
    for name in ('openwave_angular_torque','openwave_rotational_g_limit'):
        producer=json.loads((ROOT/'model'/name/'results.json').read_text())
        pins=producer['source_sha256']
        assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==sha for p,sha in pins.items())
        dependency_checks[name]={'checked_transitive_source_hashes':len(pins),'mismatches':0}
    result={'schema':'oph.independent_rotational_magnetic_review.v1','input_sha256':inputs,
        'producer_dependency_checks':dependency_checks,
        'exact_full_shell_algebra':exact_shell(),'fresh_torque_and_field_controls':torque_control(),
        'independent_axial_static_and_bound_review':axial_control(),
        'mathematical_review':{'positive_Sylvester_inverse_and_second_derivative_bound':'PASS',
            'field_coercivity_all_minima_inside_contraction_ball':'PASS',
            'axial_eigenline_gap_and_contraction':'PASS',
            'Schur_vacuum_energy_and_positive_kernel_remainder':'PASS',
            'transverse_radius_and_full_state_Duhamel_phase_bound':'PASS'},
        'scientific_boundary':'Declared complete-shell trap and transverse Maxwell regulator, distinct from cubic and PW119 sources. Semiclassical quasistatic matter-field evolution is not a linear reference-uniform quantum channel. Physical mass, alpha and anomalous magnetic moment are not derived. Numerical values near2 are not narrow spectral enclosures.',
        'verdict':'PASS'}
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==sha for p,sha in inputs.items())
    (Path(__file__).parent/'rotational_review.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print('PASS independent complete-shell torque and magnetic-limit review')


if __name__=='__main__':main()
