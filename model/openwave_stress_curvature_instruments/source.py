"""Ten native observable means acquire a complete finite tensor stencil."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import ceil
import numpy as np
import sympy as s
from model.particle_geometric_response.geometry import native_symbol,native_strain,strain_basis,clean
from model.openwave_staggered_ward import source as staggered
from model.openwave_tensor_record_transport import source as transport
from model.openwave_einstein_error import source as decoder

ETA=s.Rational(1,2)


def rational_row_bound(A):
    # Gaussian-rational l1 majorant avoids converting exact radicals to floats.
    return max(sum(abs(s.re(A[i,j]))+abs(s.im(A[i,j])) for j in range(A.cols)) for i in range(A.rows))


@lru_cache(None)
def operators():
    base=native_symbol(s.eye(3));H=base['total'];ops={(0,0):H}
    for i in range(3):ops[0,i+1]=clean(s.kronecker_product(s.eye(4),base['link_data']['D'][i]))
    for name,q in strain_basis():ops[int(name[0])+1,int(name[1])+1]=native_strain(q)['total']
    assert len(ops)==10 and all(A==A.H and A.shape==(60,60) for A in ops.values())
    v=s.zeros(60,1);v[0]=s.Rational(3,5);v[17]=s.Rational(4,5)
    rho=s.eye(60)/120+v*v.H/2
    assert s.trace(rho)==1
    bounds={k:rational_row_bound(A) for k,A in ops.items()}
    means={k:s.factor(s.trace(rho*A)) for k,A in ops.items()}
    return H,ops,rho,bounds,means


def instrument_records():
    H,ops,rho,bounds,means=operators();Hn=np.array(H).astype(complex);rn=np.array(rho).astype(complex)
    R_H=rational_row_bound(H);rows=[]
    for key,O in ops.items():
        scale=bounds[key]
        if scale==0:scale=s.S.One
        A=O/scale;mean=means[key]/scale;second=s.factor(s.trace(rho*A*A))
        probabilities=[s.factor((1+2*ETA*mean+ETA**2*second)/4),
                       s.factor((1-2*ETA*mean+ETA**2*second)/4),s.factor((1-ETA**2*second)/2)]
        assert all(p.is_nonnegative for p in probabilities) and sum(probabilities)==1
        assert s.factor(scale*(probabilities[0]-probabilities[1])/ETA)==means[key]
        An=np.array(A).astype(complex);ev,U=np.linalg.eigh(An)
        assert max(abs(ev))<=1+2e-14
        K=[(np.eye(60)+float(ETA)*An)/2,(np.eye(60)-float(ETA)*An)/2,
           (U*np.sqrt((1-float(ETA)**2*ev**2)/2))@U.conj().T]
        complete=sum(k.conj().T@k for k in K)
        output=sum(k@rn@k.conj().T for k in K)
        actual=[float(np.trace(k@rn@k.conj().T).real) for k in K]
        residual=float(np.linalg.norm(complete-np.eye(60)))
        assert residual<8e-14 and max(abs(p-float(q)) for p,q in zip(actual,probabilities))<2e-14
        work=float(np.trace(Hn@(output-rn)).real)
        assert abs(work)<=float(R_H*3*ETA**2/2)+1e-12
        rows.append({'component':list(key),'operator_norm_upper_exact':str(scale),
            'source_mean_exact':str(means[key]),'all_three_probabilities_exact':list(map(str,probabilities)),
            'record_labels_exact':[str(scale/ETA),str(-scale/ETA),'0'],
            'complete60_carrier_instrument_residual_diagnostic':residual,
            'source_work_diagnostic':work,'source_work_absolute_upper_exact':str(R_H*3*ETA**2/2)})
    return rows


def empty_tensor():
    return {(a,b):{} for a,b in product(range(4),repeat=2)}


def basis_stencils():
    _,ops,_,bounds,_=operators();keys=list(ops)
    out={k:empty_tensor() for k in keys}
    # Three supplied metric-readout profiles of actual native observables.
    specs=(((0,0),(0,0),(0,0,0,0),F(1,16)),
           ((1,1),(1,1),(0,0,0,0),F(1,32)),
           ((0,1),(0,1),(1,1,0,0),F(1,64)))
    metric=[]
    for op,ab,p,weight in specs:
        h=empty_tensor();coefficient=weight/F(bounds[op] or 1)
        h[ab]={p:coefficient};h[ab[::-1]]=h[ab]
        curvature=staggered.einstein(h)
        assert all(not staggered.divergence(curvature,b) for b in range(4))
        out[op]={k:staggered.add(out[op][k],curvature[k]) for k in curvature}
        metric.append({'source_component':list(op),'metric_component':list(ab),'site':list(p),'coefficient_exact':str(coefficient)})
    # Work-conjugate stress components at their native staggered locations.
    for key in keys:
        p=staggered.parity(key)
        out[key][key]=staggered.add(out[key][key],{p:F(-1,128)})
        out[key][key[::-1]]=out[key][key]
    for tensor in out.values():staggered.validate_tensor(tensor)
    return out,metric


def tensor_from_means(means):
    basis,_=basis_stencils()
    return {ab:staggered.add(*(staggered.scale(F(means[k]),tensor[ab]) for k,tensor in basis.items()))
            for ab in product(range(4),repeat=2)}


def acquisition_budget(epsilon=F(1,1000),failure=F(1,100)):
    if epsilon<=0 or not 0<failure<1:raise ValueError('positive precision and failure in(0,1) required')
    _,ops,_,bounds,_=operators();count=len(ops)
    exponent=0
    while F(2*count,2**exponent)>failure:exponent+=1
    rows=[]
    for key in ops:
        magnitude=F(bounds[key] or 1)/F(ETA)
        n=ceil(2*magnitude*magnitude*exponent/(epsilon*epsilon))
        rows.append({'component':list(key),'independent_copies':n})
    return {'per_native_mean_error_exact':str(epsilon),'simultaneous_failure_upper_exact':str(failure),
        'integer_exponential_majorant':exponent,'per_setting_copies':rows,
        'total_prepared_copies':sum(r['independent_copies'] for r in rows),
        'method':'Hoeffding: 2 exp(-N epsilon^2/(2 M_j^2)) <=2 exp(-L)<2*2^-L; union over10 settings.'}


def tensor_certificate():
    _,ops,_,_,means=operators();basis,metric=basis_stencils();B=tensor_from_means(means)
    sites=sorted(set().union(*(set(transport.centers(t)) for t in basis.values())))
    assert len(sites)>0
    Cmax=F()
    for ab in B:
        support=set().union(*(set(t[ab]) for t in basis.values()))
        for p in support:Cmax=max(Cmax,sum((abs(t[ab].get(p,F())) for t in basis.values()),F()))
    tau,actual,_=transport.exact_transport_check(B)
    examples=[];qmax=s.S.Zero
    # Every coefficient tensor obeys the same exact inverse. The complete
    # expectation stencil is decoded at every scalar site, including halos.
    for x in sites:
        C=transport.colocate(B,x);q=decoder.charges(C)
        recovered=decoder.reconstruct(C)
        assert s.simplify(recovered-C+decoder.metric_coefficient(C)*decoder.ETA)==s.zeros(4)
        qmax=max(qmax,max(abs(v) for v in q))
        examples.append({'site':list(x),'null_charges_exact':list(map(str,q)),
            'metric_coefficient_exact':str(decoder.metric_coefficient(C))})
    epsilon=F(1,1000);component=Cmax*epsilon;null=8*component
    divergence=max((abs(x) for b in range(4) for x in staggered.divergence(B,b).values()),default=F())
    assert qmax>0 and divergence>0
    return {'metric_readout_profiles':metric,'stress_source_coefficient_exact':'1/128',
        'native_operator_basis_size':len(ops),'complete_scalar_support_count':len(sites),
        'all_colocated_null_records':examples,'same_state_mean_tensor_decoder_identity':True,
        'largest_operator_coefficient_l1_exact':str(Cmax),
        'actual_mean_staggered_transport_budget_exact':str(tau),'actual_mean_transport_error_exact':str(actual),
        'actual_mean_maximum_null_residual_exact':str(qmax),'actual_mean_divergence_maximum_exact':str(divergence),
        'acquisition_error':{'native_mean_epsilon_exact':str(epsilon),'every_staggered_or_colocated_component_error_exact':str(component),
            'every_null_charge_error_upper_exact':str(null),
            'tracefree_decoder_max_entry_error_upper_exact':str((2+s.sqrt(3))*s.Rational(null.numerator,null.denominator)),
            'full_colocated_tensor_error_with_measured_metric_coefficient_upper_exact':str((2+s.sqrt(3))*s.Rational(null.numerator,null.denominator)+s.Rational(component.numerator,component.denominator)),
            'transport_budget_estimation_error_upper_exact':str(2*component),
            'divergence_estimation_error_upper_at_spacing1_exact':str(8*component)},
        'conditional_tensor_theorem':'For acquired null residual qhat and variation tauhat: delta=(2+sqrt(3))*(qhat+8*Cmax*epsilon)+tauhat+2*Cmax*epsilon; |B-lambda eta|<=delta. For acquired divergence what at spacing a=1, each connected lambda step <=what+8*Cmax*epsilon+8*delta.',
        'nonzero_residual_is_retained':True,'Einstein_equation_enforced_by_acquisition':False}


def backaction_controls():
    I=s.eye(2);X=s.Matrix([[0,1],[1,0]]);Z=s.diag(1,-1);rho=(I+X)/2
    K=[(I+ETA*Z)/2,(I-ETA*Z)/2,s.sqrt((1-ETA**2)/2)*I]
    assert sum((k.H*k for k in K),s.zeros(2))==I
    after=sum((k*rho*k.H for k in K),s.zeros(2))
    assert s.trace(after*X)==1-ETA**2
    missing=sum((k.H*k for k in K[:2]),s.zeros(2))
    assert missing!=I and missing==s.Rational(5,8)*I
    bell=s.Matrix([1,0,0,1])/s.sqrt(2)
    joint=sum((s.kronecker_product(k,I)*bell*bell.H*s.kronecker_product(k.H,I) for k in K),s.zeros(4))
    assert s.trace(joint)==1
    return {'same_state_X_mean_exact':'1','X_mean_after_prior_Z_instrument_exact':str(s.trace(after*X)),
        'dropped_zero_event_total_effect_exact':'5I/8','reference_joint_trace_exact':str(s.trace(joint)),
        'sequential_reuse_is_not_independent_same_state_acquisition':True}


def run():
    H,ops,rho,bounds,means=operators()
    noncommuting=[]
    keys=list(ops)
    for a,b in zip(keys,keys[1:]):
        comm=clean(ops[a]*ops[b]-ops[b]*ops[a])
        norm=s.simplify(s.trace(comm.H*comm))
        if norm:noncommuting.append({'components':[list(a),list(b)],'commutator_Hilbert_Schmidt_squared_exact':str(norm)})
    assert noncommuting
    return {'schema':'oph.native_stress_curvature_instruments.v1',
        'verdict':'PASS_COMPLETE_NATIVE_FINITE_TENSOR_ACQUISITION_WITH_RESOURCE_AND_ERROR_BOUNDS',
        'source':{'carrier':'All60 spin/internal columns of the declared one-cell native Wilson symbol; no column or measurement outcome discarded.',
            'source_status':'Inherited native-coordinate one-particle symbol with supplied noncommuting links/coframe; no normalizable full spatial-field embedding or native12-port generation is inferred.',
            'H':'native_symbol(I).total, including its declared magnetic scalar',
            'stress':'Six exact coframe derivatives from native_strain; T00=H and three native covariant momenta are derivatives of supplied lapse and shift couplings.',
            'preparation':'rho=I60/120+|psi><psi|/2, psi=(3|0>+4|17>)/5; a full-rank supplied state.',
            'prepared_source_energy_exact':str(means[0,0]),
            'preparation_work_above_ground_interval_exact':['0',str(means[0,0]+rational_row_bound(H))],
            'clock':'Every setting is acquired at the same specified native time0 on an independently prepared copy; arbitrary common time requires evolving each copy with this same H.',
            'metric':'Three explicit compact staggered profiles of native observable means; the readout assignment is declared, not a derived spacetime metric.'},
        'complete_native_instruments':instrument_records(),'noncommuting_source_pairs':noncommuting,
        'measurement_contract':{'kraus':'K+-=(I+-eta A)/2, K0=sqrt((I-eta^2 A^2)/2), A=O/R, eta=1/2',
            'all_outcomes':'+R/eta, -R/eta, 0; expectation exactly Tr(rho O)',
            'external_reference':'The Stinespring isometry tensors with any spectator reference; initialized apparatus and a classical10-setting schedule are supplied.',
            'trace_norm_disturbance_upper_exact':'3/8 for any input and quantum reference',
            'source_work_bound':'|Delta E|<=||H||*(3 eta^2/2); exact rational row bound used for ||H||',
            'backaction':'Independent preparations make all10 means refer to the same premeasurement state; incompatible values are not assigned jointly to a single copy.',
            'numerical_work_scope':'Full60-state Kraus evaluation is diagnostic; exact completeness, probability formulas and analytic work bounds are the certificate.',
            'memory_and_preparation_work':'The source preparation costs its actual energy relative to the chosen reference. Reset and schedule work are not claimed zero.'},
        'same_state_sampling':acquisition_budget(),'tensor_transfer':tensor_certificate(),
        'negative_controls':backaction_controls(),
        'finite_theorems_reused':['Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean','Lean/ObserverPatchHolography/EinsteinResidualBound.lean','Lean/ObserverPatchHolography/TensorRecordTransport.lean'],
        'scope':{'complete_finite_quantum_acquisition_to_tensor_error_interface':True,
            'source_metric_stress_physical_identification_proved':False,'native_12_port_process_compilation_proved':False,
            'nonlinear_Einstein_continuum_or_Clay_proof':False,'new_Nature_validation':False}}
