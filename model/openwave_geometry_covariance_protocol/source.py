"""Exact paired records, finite-sample bounds and common-source unit identities."""
from fractions import Fraction as F
from math import comb, ceil
from pathlib import Path
import json
import numpy as np
from model.openwave_collective_finite_strength.bounds import exp_bounds, root_bounds
from model.openwave_universal_energy_geometry.source import QSIZE as Q, OMEGA as W, profile
from model.openwave_gravity_unit_constraint.source import pi_interval
from model.openwave_collective_gravity_strength.source import inherited_units

ROOT = Path(__file__).resolve().parents[2]


def valid_n(n):
    if type(n) is not int or n < 1:
        raise ValueError('a positive integer geometry count is required')


def pair_records(n):
    """Every count pair, including equal coordinates and zero covariance records."""
    valid_n(n)
    x = [Q*F(2*k-n,n) for k in range(n+1)]
    p = [F(comb(n,k),2**n) for k in range(n+1)]
    v = Q*Q/n
    rows=[]; c_lo=c_hi=second_hi=F(0); mean_lo=mean_hi=F(0)
    variance_mean=variance_second=F(0)
    for k,a in enumerate(x):
        ea=exp_bounds(a)
        mean_lo+=p[k]*ea[0];mean_hi+=p[k]*ea[1]
        for j,b in enumerate(x):
            probability=p[k]*p[j]; difference=a-b
            eb=exp_bounds(b)
            candidates=[difference*(u-w)/2 for u in ea for w in eb]
            lo=max(F(0),min(candidates));hi=max(candidates)
            if a==b:lo=hi=F(0)
            z=difference*difference/2
            c_lo+=probability*lo;c_hi+=probability*hi
            second_hi+=probability*hi*hi
            variance_mean+=probability*z;variance_second+=probability*z*z
            rows.append({'left_count':k,'right_count':j,'probability_exact':str(probability),
                'left_Q_exact':str(a),'right_Q_exact':str(b),
                'covariance_record_interval_exact':[str(lo),str(hi)],
                'variance_record_exact':str(z)})
    assert sum(F(r['probability_exact']) for r in rows)==1
    assert variance_mean==v and variance_second==v*v*(3-F(1,n))
    lower,upper=exp_bounds(-Q)[0]*v,exp_bounds(Q)[1]*v
    assert lower<=c_lo<=c_hi<=upper
    assert second_hi<=3*exp_bounds(2*Q)[1]*v*v
    return {'N':n,'complete_count_pair_records':rows,
        'covariance_interval_exact':[str(c_lo),str(c_hi)],
        'mean_lapse_interval_exact':[str(mean_lo),str(mean_hi)],
        'variance_mean_exact':str(variance_mean),
        'variance_second_moment_exact':str(variance_second),
        'covariance_second_moment_upper_exact':str(second_hi),
        'all_N_second_moment_upper_exact':str(3*exp_bounds(2*Q)[1]*v*v)}


def sampling_plan(relative_error=F(1,10), blocks=33, mean_relative_error=F(1,100)):
    """Median of independent block means, with an exact binomial failure bound."""
    e=F(relative_error);em=F(mean_relative_error)
    if not 0<e<1 or not 0<em<1 or type(blocks)is not int or blocks<1 or blocks%2!=1:
        raise ValueError('relative error in (0,1) and a positive odd block count required')
    variance_ratio=3*exp_bounds(4*Q)[1]
    mean_variance_ratio=Q*Q/(2*(1-Q*Q))
    size=max(ceil(4*variance_ratio/(e*e)),ceil(4*mean_variance_ratio/(em*em)))
    failure=sum((F(comb(blocks,j))*F(1,4)**j*F(3,4)**(blocks-j)
                 for j in range((blocks+1)//2,blocks+1)),F(0))
    assert variance_ratio/(size*e*e)<=F(1,4)
    assert mean_variance_ratio/(size*em*em)<=F(1,4)
    return {'relative_covariance_error_exact':str(e),'independent_blocks':blocks,
        'independent_count_pairs_per_block':size,'total_independent_count_pairs':size*blocks,
        'single_block_failure_upper_exact':'1/4','median_failure_upper_exact':str(failure),
        'relative_mean_lapse_error_exact':str(em),
        'joint_covariance_and_mean_failure_upper_exact':str(min(F(1),2*failure)),
        'joint_coverage_uses_union_bound_not_independence_of_statistics':True,
        'count_N_independent_sample_bound':True,
        'independent_preparation_and_complete_coordinate_resolution_required':True}


def quantization_bias(n,eta):
    """Deterministic coordinate error propagated into the pair covariance mean."""
    valid_n(n);eta=F(eta)
    if not 0<=eta<=F(1,4):raise ValueError('coordinate error must be in [0,1/4]')
    v=Q*Q/n
    e=exp_bounds(Q+eta)[1];shift=exp_bounds(eta)[1]-1
    bound=e*(2*eta*root_bounds(2*v)[1]+2*eta*eta)+shift*exp_bounds(Q)[1]*v
    cov_floor=exp_bounds(-Q)[0]*v
    return {'absolute_coordinate_error_exact':str(eta),
        'covariance_bias_upper_exact':str(bound),
        'relative_covariance_bias_upper_exact':str(bound/cov_floor),
        'mean_lapse_relative_lower_exact':str(exp_bounds(-eta)[0]),
        'mean_lapse_relative_upper_exact':str(exp_bounds(eta)[1]),
        'resolution_is_a_supplied_apparatus_calibration':True}


def normalized_strength_interval(n):
    """Far-reference weak-rest normalization cancels the mean lapse exactly."""
    valid_n(n);a=Q/n;S=profile()['stiffness_exact'];plo,phi=pi_interval()
    # tanh(a) lies in [a-a^3/3,a]. No huge power or full carrier is allocated.
    return S*Q*(a-a**3/3)/(2*phi*W),S*Q*a/(2*plo*W)


def propagate_strength(covariance, source_rest_slope, photon_slope,
                       energy_unit, speed_ratio=(F(1),F(1))):
    """One common clock and length, with all supplied calibration intervals."""
    intervals=[tuple(map(F,x))for x in (covariance,source_rest_slope,photon_slope,energy_unit,speed_ratio)]
    if any(lo<=0 or hi<lo for lo,hi in intervals):
        raise ValueError('ordered strictly positive intervals required')
    (cl,ch),(sl,sh),(pl,ph),(ul,uh),(vl,vh)=intervals
    S=profile()['stiffness_exact'];pi_lo,pi_hi=pi_interval()
    return (vl*S*cl/(2*pi_hi*W*sh*ph*uh*uh),
            vh*S*ch/(2*pi_lo*W*sl*pl*ul*ul))


def complete_instrument_control(n):
    """Full small geometry/matter/photon carrier, with a count-Lueders instrument."""
    valid_n(n)
    if n>5:raise ValueError('the complete matrix control has N<=5')
    d=2**n;counts=np.array([i.bit_count()for i in range(d)])
    q=np.array([float(Q)*((2*k-n)/n)for k in counts])
    hg=np.eye(d)*n*float(W)/2
    for i in range(d):
        for bit in range(n):hg[i,i^(1<<bit)]=-float(W)/2
    # These complete small spectator controls test the operator identity;
    # the proof keeps the full native matter and unbounded photon number space.
    matter=np.diag([0.,1.,3.]);photon=np.diag([0.,.7,1.4,2.1])
    H=np.kron(hg,np.eye(12))+np.kron(np.diag(np.exp(q)),np.kron(matter,np.eye(4)))+np.kron(np.diag(np.exp(q/2)),np.kron(np.eye(3),photon))
    labels=np.repeat(counts,12)
    pinched=np.where(labels[:,None]==labels[None,:],H,0.)
    difference=pinched-H
    expected=np.kron(n*float(W)/2*np.eye(d)-hg,np.eye(12))
    assert np.array_equal(difference,expected)
    psi=(np.arange(d*12)+1)+1j*np.arange(d*12)[::-1];psi/=np.linalg.norm(psi)
    probabilities=[float(np.sum(abs(psi[labels==k])**2))for k in range(n+1)]
    assert abs(sum(probabilities)-1)<1e-14
    work=float(np.vdot(psi,difference@psi).real)
    assert abs(work-float(np.vdot(psi,expected@psi).real))<1e-13
    # The complete cyclic memory dilation is a permutation on every input.
    images={(i,(record+counts[i])%(n+1))for i in range(d)for record in range(n+1)}
    assert len(images)==d*(n+1)
    return {'N':n,'complete_small_control_dimension':d*12,
        'all_count_probabilities':probabilities,'arbitrary_state_measurement_work_source_units_diagnostic':work,
        'count_memory_dilation_is_a_complete_permutation':True,
        'complete_operator_work_identity_error_exact':'0',
        'small_photon_number_control_is_not_a_full_Fock_truncation_claim':True}


def run():
    u,hc,target,_,_=inherited_units()
    parent=json.loads((ROOT/'model/openwave_gravity_noise_normalization/results.json').read_text())
    inherited_n=parent['conditional_calibrated_member']['N']
    plan=sampling_plan();copies=2*plan['total_independent_count_pairs']
    family=[]
    for n in (1,2,4,8,16):
        lo,hi=normalized_strength_interval(n)
        family.append({'N':n,'count_choice':'declared small control size; no observed-G matching rule',
            'far_reference_rest_normalized_weak_strength_GeV_minus2_interval_exact':[str(lo/u**2),str(hi/u**2)],
            'single_common_far_photon_speed_matched_length_fm_exact':str(hc/u),
            'physical_finite_source_Newton_law_established':False})
    for a,b in zip(family,family[1:]):
        assert F(b['far_reference_rest_normalized_weak_strength_GeV_minus2_interval_exact'][1])<F(a['far_reference_rest_normalized_weak_strength_GeV_minus2_interval_exact'][0])
    v=Q*Q/inherited_n;eta=root_bounds(v)[0]/10000
    example=pair_records(4)
    cc=sum(map(F,example['covariance_interval_exact']))/2
    mm=sum(map(F,example['mean_lapse_interval_exact']))/2
    ec=F(plan['relative_covariance_error_exact']);em=F(plan['relative_mean_lapse_error_exact'])
    propagated=propagate_strength((cc/(1+ec),cc/(1-ec)),
        (mm/(1+em),mm/(1-em)),(1,1),(u,u))
    expected=normalized_strength_interval(4)
    assert propagated[0]<=expected[0]/u**2<=expected[1]/u**2<=propagated[1]
    return {'schema':'oph.complete_geometry_covariance_protocol.v1',
        'verdict':'PASS_COMPLETE_COVARIANCE_INSTRUMENT_AND_NORMALIZATION_DOMAIN',
        'source':{'Hamiltonian':'Hg+exp(Q)K_source+exp(F Q)K_probe+exp(f_gamma Q)K_photon',
            'matter':'Full positive native Gaussian Fock energies transported by the exact native intertwiner; all orthogonal sectors retained',
            'photon':'Declared positive two-polarization free-photon Fock energy; complete number space retained',
            'geometry':'Hg=omega/2 sum(I-X), Q=q/N sum Z, all 2^N states',
            'complete_operator_domain':'D(K_photon); finite geometry and native matter factors are bounded',
            'same_time_parameter_and_length_for_all_factors':True,
            'source_couplings_are_declared_architecture':True},
        'complete_instrument':{'effects':'P_k(Q) tensor identity_matter tensor identity_photon, every k=0,...,N',
            'memory_dilation':'sum_k P_k tensor cyclic_shift^k, with N+1 degenerate memory levels',
            'all_count_pairs_retained':True,
            'pair_covariance_record':'D(x,y)=(x-y)(exp(x)-exp(y))/2',
            'pair_variance_record':'Z(x,y)=(x-y)^2/2',
            'unbiased_means':'E D=Cov(Q,exp Q), E Z=q^2/N',
            'work_operator':'pinching(H)-H=omega/2 sum X on geometry, identity on matter/photon',
            'vacuum_work_per_independent_pair_GeV':'N*omega*u',
            'ground_preparation_and_blank_memories_supplied':True,
            'memory_erasure_or_autonomous_detector_work_source_derived':False},
        'all_N_statistical_theorem':{'variance_record_second_moment':'(q^2/N)^2*(3-1/N)',
            'covariance_record_second_moment_upper':'3 exp(2q)*(q^2/N)^2',
            'covariance_lower':'exp(-q)*q^2/N',
            'relative_covariance_record_variance_upper':'3 exp(4q)',
            'relative_pair_mean_lapse_variance_upper':'v/[2(1-v)] <= q^2/[2(1-q^2)]',
            'sample_plan':plan},
        'error_propagation_example':{'N':4,'observation_role':'synthetic model means; no physical noise observation',
            'covariance_estimate_exact':str(cc),'mean_lapse_estimate_exact':str(mm),
            'conditional_weak_G_interval_GeV_minus2_exact':list(map(str,propagated)),
            'joint_failure_upper_exact':plan['joint_covariance_and_mean_failure_upper_exact'],
            'coordinate_error_assumed_exact':'0','independent_G_or_N_inference_from_Nature':False},
        'complete_small_pair_records':[pair_records(n)for n in (1,2,4,8)],
        'complete_small_operator_controls':[complete_instrument_control(n)for n in (1,2,3,4)],
        'unit_and_rest_normalization':{'M(f)':'cosh(f*q/N)^N',
            'weak_source_rest_slope':'epsilon_source_prime(0)=M(1)',
            'far_probe_rest_slope':'M(0)=1',
            'dressed_photon_group_speed_in_source_units':'M(f_gamma)',
            'bare_label_coefficient':'g_N=S Cov(Q,exp Q)/(2 pi omega)',
            'rest_normalized_coefficient':'g_N/M(1)',
            'common_unit_invariant':'G_weak/(c_gamma/c)=g_N/[u^2 M(1) M(f_gamma)]',
            'far_reference_cancellation':'f_gamma=0: G_weak=S q tanh(q/N)/(2 pi omega u^2) when c_gamma=c',
            'near_reference_control':'f_gamma=1 gives one additional divisor M(1)',
            'finite_native_source_masses_require_epsilon(E),not_M(1)*E':True,
            'finite_native_poles_or_clock_dictionary_transferred':False},
        'independent_architecture_counterfamily':family,
        'inherited_G_selected_resource_example':{'N':inherited_n,
            'N_is_selected_from_observed_G':True,'independent_G_prediction':False,
            'source_coordinate_variance_exact':str(v),
            'hypothetical_coordinate_resolution':quantization_bias(inherited_n,eta),
            'independent_geometry_preparations':copies,
            'complete_measurement_work_GeV_exact':str(F(copies)*inherited_n*W*u/2),
            'calibration_target_G_GeV_minus2_exact':str(target)},
        'scope':{'complete_covariance_protocol_on_declared_common_extension':True,
            'source_count_fixed_by_independent_physical_noise_data':False,
            'physical_metric_noise_observable_identified':False,
            'physical_Newton_strength_independently_validated':False,
            'native_twelve_port_process_realization':False,
            'finite_energy_universality_or_Einstein_continuum_proved':False,
            'additional_Nature_validations':0}}
