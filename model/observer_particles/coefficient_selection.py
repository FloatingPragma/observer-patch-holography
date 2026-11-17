#!/usr/bin/env python3
"""Exact invariant coefficient spaces and conditional gauge-ray closure.

No empirical masses/couplings enter. The action, A3 output-map and physical
scale boundaries are explicit; this is not a pair of full OPH countermodels.
"""
from __future__ import annotations
from evidence_contract import public_receipt
import hashlib, itertools, json, math, sys
from pathlib import Path
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins
sys.dont_write_bytecode = True
import numpy as np
import sympy as s
from scipy.optimize import brentq
from mpmath import iv

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PREVIOUS = HERE
sys.path.insert(0, str(PREVIOUS))
import representation_identification as rep
import scalar_bilinear_route as pair

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def sparse(matrix): return rep.sparse(matrix)
def clean(matrix): return s.SparseMatrix(matrix.applyfunc(s.simplify))
def zero(matrix): return all(s.expand(x)==0 for x in matrix.todok().values())

def commutant(generators, cartans):
    """Solve the full commutator equations after exact weight pruning."""
    n=generators[0].rows
    weights=[tuple(h[i,i] for h in cartans) for i in range(n)]
    entries=[(i,j) for i in range(n) for j in range(n) if weights[i]==weights[j]]
    lookup={ij:k for k,ij in enumerate(entries)}
    rows=[]
    for g in generators:
        equations={}
        for (i,k),v in g.todok().items():
            for j in range(n):
                if (k,j) in lookup:
                    row=equations.setdefault((i,j),{})
                    z=lookup[k,j];row[z]=row.get(z,0)+v
                if (j,i) in lookup:
                    row=equations.setdefault((j,k),{})
                    z=lookup[j,i];row[z]=row.get(z,0)-v
        rows.extend({k:s.expand(v) for k,v in r.items() if v} for r in equations.values())
    rows=[r for r in rows if r]
    matrix=s.SparseMatrix(len(rows),len(entries),{(i,j):v for i,r in enumerate(rows) for j,v in r.items()})
    kernel=matrix.nullspace()
    out=[s.SparseMatrix(n,n,{ij:v[k] for k,ij in enumerate(entries) if v[k]}) for v in kernel]
    assert all(zero(g*x-x*g) for g in generators for x in out)
    return out, {'weight_pruned_unknowns':len(entries),'nonzero_equation_rows':len(rows),
                 'exact_rank':len(entries)-len(out),'complex_nullity':len(out)}

def spectrum_projectors(matrix):
    values=sorted(matrix.eigenvals(), key=lambda x:float(x))
    out=[]
    for v in values:
        p=s.eye(matrix.rows)
        for w in values:
            if w!=v:p=p*(matrix-w*s.eye(matrix.rows))/(v-w)
        out.append((v,clean(p)))
    return out

def entropy_moments(alpha, kind, cutoff):
    if kind==2:
        n=np.arange(cutoff+1,dtype=float)
        dims=n+1;cas=n*(n+2)/4
    else:
        pairs=[(p,q) for p in range(cutoff+1) for q in range(cutoff+1-p)]
        p,q=np.array(pairs,dtype=float).T
        dims=(p+1)*(q+1)*(p+q+2)/2;cas=(p*p+q*q+p*q)/3+p+q
    weights=dims*np.exp(-4*math.pi**2*alpha*cas);prob=weights/weights.sum()
    ell=float(prob@np.log(dims))
    cov=float(prob@(cas*np.log(dims))-(prob@cas)*ell)
    return ell,cov

def entropy_monotonicity_certificate():
    """Interval-check an analytic full-SU3 covariance bound, not a grid test.

    Positive summands decrease between each pair of alpha endpoints. Bounds
    therefore apply throughout each interval, including the omitted tail.
    Above .04, E(C log d)>=log(3)E(C) and E(log d)<log(3) suffice.
    """
    iv.dps=30
    pi=iv.pi;N=20
    def terms(a):
        total=[iv.mpf(0) for _ in range(4)]
        for p in range(N+1):
            for q in range(N+1-p):
                dim=iv.mpf((p+1)*(q+1)*(p+q+2))/2
                cas=iv.mpf(p*p+q*q+p*q)/3+p+q
                weight=dim*iv.exp(-4*pi*pi*a*cas);ld=iv.ln(dim)
                for i,v in enumerate([weight,weight*cas,weight*ld,weight*cas*ld]):total[i]+=v
        return total
    def tails(a):
        n=N+1;z=4*pi*pi*a
        base=iv.mpf((n+1)*(n+2)**3)/8*iv.exp(-z*(iv.mpf(n*n)/4+n))
        ratio=iv.mpf(n+2)/(n+1)*(iv.mpf(n+3)/(n+2))**3*iv.exp(-z*(iv.mpf(2*n+1)/4+1))
        logtail=base*(n+2)/(1-ratio*iv.mpf(n+3)/(n+2))
        c=iv.mpf(n*n)/3+n;cn=iv.mpf((n+1)**2)/3+n+1
        ctail=base*c/(1-ratio*cn/c)
        return ctail,logtail
    boundaries=[iv.mpf(19)/1000+iv.mpf(i)*21/40000 for i in range(41)]
    values=[terms(a) for a in boundaries];bounds=[]
    for i in range(40):
        lo,hi=values[i],values[i+1];tc,tl=tails(boundaries[i])
        lower=hi[3]*hi[0]-(lo[1]+tc)*(lo[2]+tl)
        assert lower.a>0
        bounds.append(float(lower.a))
    at04=terms(iv.mpf(1)/25);_,tl=tails(iv.mpf(1)/25)
    upper=at04[2]+tl
    assert upper.b<iv.ln(3).a
    def mean3(a):
        total=terms(a);_,tail=tails(a)
        remainder=iv.mpf([0,tail.b])
        return (total[2]+remainder)/(total[0]+remainder)
    def mean2(a):
        z=iv.mpf(0);num=iv.mpf(0);factor=4*pi*pi*a
        for n in range(N+1):
            weight=(n+1)*iv.exp(-factor*iv.mpf(n*(n+2))/4)
            z+=weight;num+=weight*iv.ln(n+1)
        n=N+1
        first=(n+1)**2*iv.exp(-factor*iv.mpf(n*(n+2))/4)
        ratio=(iv.mpf(n+2)/(n+1))**2*iv.exp(-factor*iv.mpf(2*n+3)/4)
        tail=first/(1-ratio);remainder=iv.mpf([0,tail.b])
        return (num+remainder)/(z+remainder)
    alpha=(79+48*iv.sqrt(2)+42*iv.sqrt(3)+3*iv.sqrt(6))/31104
    inv=1/alpha;budget=((1+iv.sqrt(5))/2+iv.sqrt(pi)*alpha)/4
    amin=8*alpha/3
    assert amin.a>iv.mpf(19)/1000
    left=mean2(amin)+mean3(amin)-budget
    right=mean2(134*alpha/23)-budget
    assert left.a>0 and right.b<0
    bracket=[iv.mpf('5.002666247'),iv.mpf('5.002666250')]
    residuals=[]
    for t in bracket:
        residuals.append(mean2(1/(3*inv/8-109*t/24))+mean3(1/(3*inv/8-67*t/8))-budget)
    assert residuals[0].a>0 and residuals[1].b<0
    return {'domain':'alpha3>=19/1000',
        'interval_count':40,'finite_shell_cutoff':20,'precision_digits':30,
        'least_interval_lower_bound_for_unnormalized_covariance':min(bounds),
        'at_alpha_point04_logdimension_numerator_upper':str(upper),
        'above_point04_proof':'E(C log d)>=log(3)E(C), E(log d)<=numerator<log(3); hence covariance>0',
        'SU2_proof':'C(n) and log(n+1) strictly increase; positive weighted covariance is strictly positive',
        'full_entropy_endpoint_residual_enclosures':[str(left),str(right)],
        'certified_t_root_bracket':['5.002666247','5.002666250'],
        'root_bracket_residual_enclosures':list(map(str,residuals)),
        'strict_monotonicity_for_full_sums':True}

def gauge_ray():
    alpha=(79+48*math.sqrt(2)+42*math.sqrt(3)+3*math.sqrt(6))/31104
    inv=1/alpha;P=(1+math.sqrt(5))/2+math.sqrt(math.pi)*alpha
    tmax=3*inv/67
    def inverse(t):return np.array([5*inv/8+109*t/24,3*inv/8-109*t/24,3*inv/8-67*t/8])
    def residual(t,cutoff=40):
        x=inverse(t)
        return entropy_moments(1/x[1],2,cutoff)[0]+entropy_moments(1/x[2],3,cutoff)[0]-P/4
    root=brentq(residual,0,tmax*(1-1e-10),xtol=5e-15)
    c=(3*inv-11*root)/16;x=inverse(root);a=1/x
    ell2,cov2=entropy_moments(a[1],2,40);ell3,cov3=entropy_moments(a[2],3,40)
    derivative=-4*math.pi**2*(cov2*(109/24)*a[1]**2+cov3*(67/8)*a[2]**2)
    assert c>0 and 0<root<tmax and derivative<0
    assert abs(x[0]+x[1]-inv)<1e-12 and abs(residual(root))<1e-13
    assert np.max(np.abs(x-(c*np.array([10/3,2,2])+root*np.array([41/6,-19/6,-7]))))<1e-12
    higher=brentq(lambda t:residual(t,60),0,tmax*(1-1e-10),xtol=5e-15)
    # Uniform positive tail majorants at the smallest admitted alpha.
    # SU2: d log d <= d^2. SU3 shell s: count s+1, d<=(s+2)^3/8,
    # log d<=s+2, C>=s^2/4+s. Successive envelope ratios decrease.
    z=4*math.pi**2/(3*inv/8);n=41
    f2=(n+1)**2*math.exp(-z*n*(n+2)/4)
    r2=((n+2)/(n+1))**2*math.exp(-z*(2*n+3)/4)
    f3=(n+1)*(n+2)**4/8*math.exp(-z*(n*n/4+n))
    r3=(n+2)/(n+1)*((n+3)/(n+2))**4*math.exp(-z*((2*n+1)/4+1))
    tail2=f2/(1-r2);tail3=f3/(1-r3)
    # Numerator tails majorize denominator tails; ell <= log(max included d).
    entropy_error=tail2*(1+math.log(41))+tail3*(1+math.log((42**3)/8))
    # Full 3x3 Jacobian of (ray eliminated EM, ray plane, entropy) in xY,x2,x3.
    jac=np.array([[1.,1.,0.],[69.,-333.,218.],
                  [0.,4*math.pi**2*cov2/x[1]**2,4*math.pi**2*cov3/x[2]**2]])
    grid=np.linspace(0,tmax*(1-1e-8),501)
    vals=[residual(float(t)) for t in grid]
    mono=entropy_monotonicity_certificate()
    assert 1/(3*inv/8)>.019
    return {'alpha_Chain_candidate':alpha,'P_candidate':P,'alpha_inverse':inv,
        'equations':['x=c*(10/3,2,2)+t*(41/6,-19/6,-7)','xY+x2=1/alpha_Chain',
                     'ell2(1/x2)+ell3(1/x3)=P/4'],
        'premises':['matter-trace ray used as physical kinetic boundary',
                    'three-family one-elementary-doublet one-loop beta column',
                    'same-scale tree EM matching to Chain scalar without thresholds',
                    'retained OPH representation entropy law'],
        't_log_MU_over_mu_div_2pi':root,'ray_scale_c':c,
        'alphaY':float(a[0]),'alpha2':float(a[1]),'alpha3':float(a[2]),
        'sin2_theta_same_scale_tree':float(x[1]/inv),
        'MU_over_mu':math.exp(2*math.pi*root),
        'entropy_terms':[ell2,ell3],'entropy_residual':residual(root),
        'positive_domain_t':[0,tmax],'endpoint_residuals':[vals[0],vals[-1]],
        'scalar_residual_derivative_at_root':derivative,
        'joint_constraint_jacobian_determinant':float(np.linalg.det(jac)),
        'joint_constraint_jacobian_rank':int(np.linalg.matrix_rank(jac)),
        'grid_sign_change_count':sum(x*y<0 for x,y in zip(vals,vals[1:])),
        'global_uniqueness_claimed':True,'local_isolated_solution':True,
        'global_uniqueness_scope':'positive c and t>=0 with all couplings positive, the stated same-scale EM/ray/entropy equations',
        'monotonicity_certificate':mono,
        'endpoint_scope':'interval enclosures prove positive residual at t=0 and negative limiting residual as t tends to3A/67, including full-sum tail bounds',
        'cutoff40_vs60_root_difference':abs(higher-root),
        'uniform_analytic_entropy_tail_majorant_float_evaluation':entropy_error,
        'tail_scope':'SU2 n<=40; SU3 p+q<=40. Analytic geometric majorants evaluated float64, not interval-certified rounding.',
        'physical_matching_or_thresholds_derived':False,'absolute_scale_selected':False,
        'old_transmutation_law_imposed':False,'empirical_inputs':[]}

def main():
    packet=json.loads((PREVIOUS/'results/representation_identification_operators.json').read_text())
    scalar=json.loads((PREVIOUS/'results/scalar_bilinear_operators.json').read_text())
    _,dg=rep.fock_data();indices=packet['even_indices']
    def matter(h):return s.SparseMatrix(dg(s.SparseMatrix(h)).extract(indices,indices))
    q=pair.load_matrix(packet['even_identification_operators']['q'])
    h1=matter(rep.unit(5,0,0)-rep.unit(5,1,1));h2=matter(rep.unit(5,1,1)-rep.unit(5,2,2))
    t3=matter((rep.unit(5,3,3)-rep.unit(5,4,4))/2)
    raising=[matter(rep.unit(5,i,j)) for i,j in [(0,1),(1,0),(1,2),(2,1),(3,4),(4,3)]]
    mcomm,mrank=commutant(raising,[q,h1,h2,t3]);assert len(mcomm)==5
    norm=[1 if a==b else s.sqrt(2) for a,b in scalar['symmetric_pair_basis']]
    scale=s.diag(*norm);scaleinv=s.diag(*[1/x for x in norm])
    def lift(h):return clean(scaleinv*pair.pair_lift(h,True)[0]*scale)
    pgens=[lift(g) for g in raising];pcartans=[lift(h) for h in [q,h1,h2,t3]]
    pcomm,prank=commutant(pgens,pcartans)
    assert len(pcomm)==26
    # Convert back to the normalized CAR-pair basis, then verify all native currents.
    normalized=[clean(scale*x*scaleinv) for x in pcomm]
    native=[pair.pair_lift(matter(pair.load_matrix(g)),True)[0] for g in packet['spin_lifted_one_body_generators']]
    assert all(zero(g*x-x*g) for g in native for x in normalized)
    cp=pair.load_matrix(scalar['color_C2']);wp=pair.load_matrix(scalar['weak_C2']);qp=pair.load_matrix(scalar['q'])
    rows=[]
    for charge in sorted(set(qp.diagonal())):
        ix=[i for i,v in enumerate(qp.diagonal()) if v==charge]
        cb=cp.extract(ix,ix);wb=wp.extract(ix,ix)
        for cv,pc in spectrum_projectors(cb):
            for wv,pw in spectrum_projectors(wb):
                dim=int(s.trace(pc*pw))
                if dim:
                    color_dim={s.Rational(0):1,s.Rational(4,3):3,s.Rational(3):8,s.Rational(10,3):6}[cv]
                    weak_dim={s.Rational(0):1,s.Rational(3,4):2,s.Rational(2):3}[wv]
                    mult=s.Rational(dim,color_dim*weak_dim);assert mult.is_Integer
                    rows.append({'q6Y':str(charge),'color_C2':str(cv),'weak_C2':str(wv),
                                 'dimension':dim,'irrep_dimension':color_dim*weak_dim,'multiplicity':int(mult)})
    assert sum(r['dimension'] for r in rows)==120
    assert sum(r['multiplicity']**2 for r in rows)==26
    # Actual duplicate positive-charge scalar pair sectors: Le and Qd.
    pd=pair.load_matrix(scalar['Higgs_charge_doublet_projector'])
    pairs=scalar['symmetric_pair_basis'];t3p=pair.pair_lift(t3,True)[0]
    jm=pair.pair_lift(matter(rep.unit(5,4,3)),True)[0]
    columns=[]
    for origins in [(-3,6),(1,2)]:
        ix=[i for i,(a,b) in enumerate(pairs) if tuple(sorted((q[a,a],q[b,b])))==origins and t3p[i,i]==s.Rational(1,2)]
        images=[pd[:,i] for i in ix];v=next(v for v in images if not zero(v))
        v=clean(v/s.sqrt((v.H*v)[0]));low=clean(jm*v)
        assert (low.H*low)[0]==1
        columns.extend([v,low])
    v=s.Matrix.hstack(*columns);assert v.H*v==s.eye(4)
    for g in native:assert zero(g*v-v*(v.H*g*v))
    # Two inequivalent positive kernels in one fixed induced CAR metric.
    k0=s.diag(1,1,2,2);k1=s.diag(1,1,3,3)
    hp0=clean(v*k0*v.H);hp1=clean(v*k1*v.H)
    assert all(zero(g*h-h*g) for g in native for h in [hp0,hp1])
    plus=s.Matrix([1,1])/s.sqrt(2);rho=s.kronecker_product(plus*plus.H,s.eye(2)/2)
    effect=s.kronecker_product(plus*plus.H,s.eye(2))
    p=[]
    for k in [k0,k1]:
        u=s.diag(*[s.exp(-s.I*s.pi*d) for d in k.diagonal()])
        p.append(s.simplify(s.trace(effect*u*rho*u.H)))
    assert p==[0,1]
    assert all(zero((v.H*g*v)*rho-rho*(v.H*g*v)) for g in native)
    # Trace indices and exact family cancellation, on this actual M15.
    y=q/6;color_t=matter((rep.unit(5,0,0)-rep.unit(5,1,1))/2)
    trace_indices=[s.trace(y*y),s.trace(t3*t3),s.trace(color_t*color_t)]
    assert trace_indices==[s.Rational(10,3),s.Rational(2),s.Rational(2)]
    ng,nh,x1,x2,x3=s.symbols('ng nh x1 x2 x3')
    beta=s.Matrix([20*ng/9+nh/6,-s.Rational(22,3)+4*ng/3+nh/6,-11+4*ng/3])
    determinant=s.expand(s.Matrix.hstack(s.Matrix([x1,x2,x3]),s.Matrix(trace_indices),beta).det())
    assert s.diff(determinant,ng)==0
    assert s.expand(determinant.subs(nh,1)-(-s.Rational(23,3)*x1+37*x2-s.Rational(218,9)*x3))==0
    # A full fixed faithful state and common modular clock has no such freedom:
    # [deltaH,Eij]=0 for every Eij forces deltaH to be scalar. Vacuum-relative
    # pair amplitudes additionally distinguish that pair-sector scalar.
    support=sorted(set().union(*(set(x.todok()) for x in normalized)))
    coeff=s.Matrix([[x[i,j] for x in normalized] for i,j in support])
    assert coeff.rank()==26
    output=HERE/'results';output.mkdir(exist_ok=True)
    artifact=output/'coefficient_selection_operators.json'
    artifact.write_text(json.dumps({'normalized_pair_commutant_basis':[sparse(x) for x in normalized],
        'two_positive_doublet_isometry':sparse(v),'kernel_control_A':sparse(k0),'kernel_control_B':sparse(k1)},indent=2)+'\n', encoding='utf-8', newline='\n')
    sources=[PREVIOUS/'representation_identification.py',PREVIOUS/'scalar_bilinear_route.py',
        PREVIOUS/'results/representation_identification_operators.json',PREVIOUS/'results/scalar_bilinear_operators.json',
        ROOT/'claims/axiom_registry.yaml',
        ROOT/'code/P_derivation/paper_math.py']
    sources += [ROOT/'Lean/Screen'/n for n in ['KineticFormDichotomy.lean','KineticFamilyCancellation.lean',
        'RGRepresentationFrontier.lean','GaugeKineticInvariantForms.lean','SMStructureAdequacySurface.lean']]
    result={'status':'EXACT_COEFFICIENT_SPACES_AND_CONDITIONAL_GAUGE_RAY_CLOSURE_PASS',
        'producer_sha256':sha(__file__),'source_sha256':bound_pins({str(p.relative_to(ROOT)):sha(p) for p in sources}),
        'artifacts':{artifact.name:sha(artifact)},'empirical_inputs':[],
        'matter_commutant':mrank,'full_Lorentz_scalar_pair_commutant':prank,
        'pair_isotypic_decomposition':rows,'all12_current_commutators_for_all26_basis_elements_zero':True,
        'invariant_Hermitian_pair_kernel_real_dimension':26,
        'color_singlet_pair_kernel_real_dimension':7,
        'doublet_number_preserving_pair_kernel_real_dimension':5,
        'charge_aligned_effective_doublet_kinetic_Hermitian_dimension':9,
        'effective_field_dimension_caveat':'Herm3 is for independent charge-aligned effective doublet coordinates; not nine independent number-preserving composite pair kernels or a CAR/Fierz-reduced quartic count.',
        'Yukawa_intertwiner_complex_dimensions':{'one_independent_effective_doublet':3,'three_independent_effective_doublets':9,
            'scope':'one matter copy; full color-singlet doublet kernels derived above; conjugate Higgs slot used for two channels. Composite substitution can impose CAR/Fierz relations and is not counted here.'},
        'canonical_trace_metric':{'actual_M15_indices_Y_weak_color':[str(v) for v in trace_indices],
            'physical_action_identification_selected':False,
            'induced_CAR_pair_metric_fixed_in_controls':True,
            'port_vs_matter_SU2_SU3_Killing_ratios':['6','3/2']},
        'naturality':{'under_internal_source_transport':'all26 kernels commute with every source generator, hence its exponentials',
            'under_basis_change':'kernel, metric and detectors transform together; this covariance does not equate changed kernels at fixed detectors',
            'under_all_module_automorphisms_as_extra_symmetry':'only scalar on each of17 isotypic sectors; retains17 real weights',
            'all_A1_A2_refinement_constraints_instantiated':False},
        'gauge_invariant_observer_control':{'pair_domain_dimension':4,'common_metric':'I4',
            'kernel_eigenvalues_A':[1,1,2,2],'kernel_eigenvalues_B':[1,1,3,3],
            'state':'|+><+| in multiplicity tensor I2/2','effect':'|+><+| tensor I2',
            'time_parameter':'declared pair-unitary model coordinate s=pi',
            'complete_binary_probabilities_A':['0','1'],'complete_binary_probabilities_B':['1','0'],
            'same_representation_anomaly_and_gauge_constraints':True,'same_faithful_KMS_reference_and_modular_clock':False,
            'new_native_pulse_execution':False,'all_OPH_axiom_countermodels_claimed':False},
        'entropy_and_clock':{'I4_over4_maxentropy_stationary_for_both':True,
            'stationarity_not_equal_modular_generator':True,
            'fixed_full_faithful_reference_and_full_modular_flow_fix_H_up_to_scalar':True,
            'current_only_commutators_blind_dimension':26,
            'complete_vacuum_relative_pair_response_linear_rank':26},
        'one_loop_family_cancellation':{'determinant':str(determinant),'ng_derivative_exact_zero':True},
        'conditional_matter_trace_Chain_P_entropy_gauge_join':gauge_ray(),
        'source_selection_boundary':'The combined executed gauge/matter consequences constrain representation and conditional gauge parameters. No inspected map from their A3 optimizer/P constraints to one point of the scalar/Yukawa coefficient space is supplied.'}
    (output/'coefficient_selection.json').write_text(json.dumps(public_receipt(result),indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status':result['status'],'pair_commutant':prank,'gauge_join':result['conditional_matter_trace_Chain_P_entropy_gauge_join']},indent=2))

if __name__=='__main__':main()
