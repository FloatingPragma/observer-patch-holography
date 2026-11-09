"""Joint clock/address instruments without a conserved momentum assumption."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
import sympy as s

from model.particle_common_contract import contract
from model.particle_filtered_probe import filter as apparatus
from model.particle_primary_momentum import momentum
from model.particle_primary_mass_readout import readout
from model.particle_primary_survival import survival
from model.spin_shell_readout import bounds as clock_bounds

SOURCE_INPUTS=(
    'model/particle_primary_joint_spectroscopy/spectroscopy.py',
    'model/particle_primary_joint_spectroscopy/PROOF.md',
    'model/particle_common_contract/contract.py',
    'model/particle_common_contract/PROOF.md',
    'model/particle_primary_momentum/momentum.py',
    'model/particle_primary_momentum/PROOF.md',
    'model/particle_primary_survival/survival.py',
    'model/particle_primary_survival/PROOF.md',
    'model/particle_primary_mass_readout/readout.py',
    'model/particle_primary_mass_readout/PROOF.md',
    'model/particle_filtered_probe/filter.py',
    'model/particle_filtered_probe/PROOF.md',
    'model/spin_shell_readout/bounds.py',
    'model/spin_shell_readout/PROOF.md',
    'model/spectral_readout/kernels.py',
    'model/particle_primary_variational_readout/variational.py',
    'model/particle_primary_variational_readout/PROOF.md',
    'model/particle_primary_endpoint_bounds/endpoint.py',
    'model/particle_primary_endpoint_bounds/PROOF.md',
)


def energy_instrument(unitary, bins=4):
    """Small exact execution; U is the specified source-clock step."""
    u=apparatus.exact_matrix(unitary);J=apparatus.integer(bins,'clock bins',2)
    if J>16 or u.rows>8:raise ValueError('exact matrix control is limited to16 bins and8 states')
    if apparatus.clean(u.H*u)!=s.eye(u.rows):raise ValueError('a unitary source clock required')
    phase=s.expand_complex(s.exp(2*s.pi*s.I/J))
    powers=[u**r for r in range(J)]
    rows=tuple(apparatus.clean(sum((phase**(a*r)*powers[r] for r in range(J)),s.zeros(u.rows))/J)
               for a in range(J))
    if apparatus.clean(sum((k.H*k for k in rows),s.zeros(u.rows)))!=s.eye(u.rows):
        raise ArithmeticError('complete clock kernel lost probability')
    return rows


def joint_instrument(energy_rows, spatial_rows):
    first=tuple(apparatus.exact_matrix(k) for k in energy_rows)
    if not first or not spatial_rows:raise ValueError('complete nonempty energy and spatial instruments required')
    d=first[0].rows
    for menu in (first,tuple(row['K'] for row in spatial_rows)):
        if any(k.shape!=(d,d) for k in menu) or apparatus.clean(sum((k.H*k for k in menu),s.zeros(d)))!=s.eye(d):
            raise ValueError('each input instrument must be complete on the same source')
    if len(first)**2*len(spatial_rows)>8192:raise ValueError('explicit joint control exceeds8192 records')
    rows=[]
    for a,initial in enumerate(first):
        for spatial in spatial_rows:
            for c,final in enumerate(first):
                rows.append({'initial_energy_bin':a,'branch':spatial['branch'],
                    'wavevector_index':spatial['wavevector_index'],'final_energy_bin':c,
                    'K':apparatus.clean(final*spatial['K']*initial)})
    if apparatus.clean(sum((row['K'].H*row['K'] for row in rows),s.zeros(d)))!=s.eye(d):
        raise ArithmeticError('the complete two-clock/address history lost a record')
    return rows


def clock_plan(member, resolution='1/64', tail_error='1/1024', clock_error='1/1024'):
    """Reuse the complete summed Fejer tail with an alias-free signed energy chart."""
    contract.verify_member(member)
    eta=apparatus.exact(resolution,'energy resolution',True)
    beta=apparatus.exact(tail_error,'one-clock tail allowance',True)
    error=apparatus.exact(clock_error,'one-clock complete implementation allowance',True)
    if max(beta,error)>=1:raise ValueError('normalized allowances below one required')
    source=member['preparation_plan']['source']
    norm=max(abs(F(source['spectral_lower'])),abs(F(source['spectral_upper'])))
    tau=1/(4*norm)
    if tau*eta>1:raise ValueError('signed chart requires tau*resolution<=1')
    margin=tau*eta/7 # 2pi<7 makes this a lower bound on the true turn margin.
    plan=clock_bounds.clock_plan(margin,beta);J=plan['clock_size']
    square=F(J-1,768)/error # tau^3 ||H||^3/12 =1/768.
    steps=isqrt(-(-square.numerator//square.denominator))
    if steps*steps<square:steps+=1
    steps=max(1,steps)
    actual=F(J-1,768*steps*steps)
    if actual>error:raise ArithmeticError('repeated source-clock error exceeds allowance')
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'energy_resolution':eta,'source_norm_upper':norm,'clock_step':tau,
        'signed_energy_center_formula':'2pi*(a if 2a<J else a-J)/(J*tau)',
        'source_energy_zero_changed':False,'source_spectrum_phase_radius_upper':F(1,4),
        'alias_probability_upper':0,'clock_plan':plan,
        'Cayley_substeps_per_base_tick':steps,
        'one_complete_clock_error_upper':actual,
        'two_complete_clock_error_upper':min(F(1),2*actual),
        'maximum_source_time_per_readout':(J-1)*tau,
        'clock_histogram_executed':False})


def common_payload(member,resolution='1/64',tail_error='1/1024',clock_error='1/1024'):
    spatial=momentum.common_payload(member)
    plan=clock_plan(member,resolution,tail_error,clock_error)
    J=plan['clock_plan']['clock_size'];M=spatial['neutral_address_dimension']
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'spatial_plan':spatial,'energy_plan':plan,
        'complete_history_Kraus':'F_c(H) K_(b,k) F_a(H)',
        'complete_fine_outcomes':2*M*J*J,
        'two_energy_register_dimensions':[J,J],
        'spatial_register_dimension':2*M,
        'compact_two_bit_Ramsey_carrier_bound_applies':False,
        'initial_source_marginal_preserved':True,
        'initial_energy_record_preserves_joint_purification':False,
        'initial_energy_record_retained':True,'all_failures_and_reference_retained':True,
        'joint_source_law':'sum_ij p_i |K_(b,k),ji|^2 |F_a(E_i)|^2 |F_c(E_j)|^2',
        'complete_error_to_exact_two_clock_finite_spatial_instrument_upper':plan['two_complete_clock_error_upper'],
        'additional_error_if_comparing_ideal_temporal_filter':spatial['complete_error_to_ideal_filtered_aperture_upper'],
        'energy_label_is_transfer_only_after_subtraction':True,
        'spatial_label_is_conserved_momentum':False,
        'resolved_primary_ridge_or_branch_computed':False,
        'full_joint_primary_histogram_computed':False,
        'continuum_limit_required_for_this_finite_instrument':False})


def observed_transfer_to_excitation(member,transfer_interval,original_weight_lower,
                                   resolution='1/64',tail_error='1/1024',clock_error='1/1024',
                                   initial_excitation_cut='1/8',additional_event_error='0'):
    """Conditional adapter for an actual original joint-record event, not supplied observations."""
    plan=clock_plan(member,resolution,tail_error,clock_error)
    lo,hi=(F(x) for x in transfer_interval)
    if lo>hi:raise ValueError('ordered recorded transfer interval required')
    weight=apparatus.exact(original_weight_lower,'original event lower weight')
    extra=apparatus.exact(additional_event_error,'additional complete/sampling event error')
    cut=apparatus.exact(initial_excitation_cut,'initial excitation cutoff',True)
    if max(weight,extra)>1:raise ValueError('normalized original event data required')
    eta=F(resolution);tail=F(plan['clock_plan']['complete_clock_tail_upper'])
    implementation=F(plan['two_complete_clock_error_upper'])
    cooling=F(member['state']['target_ground_excess'])/cut
    total=2*tail+implementation+extra+cooling
    if hi+2*eta+cut<0:raise ValueError('the event has no admitted nonnegative final excitation')
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'supplied_original_joint_event_lower':weight,
        'final_excitation_interval':[max(F(0),lo-2*eta),hi+2*eta+cut],
        'original_final_excitation_weight_lower':max(F(0),weight-total),
        'two_energy_tail_allowance':2*tail,'complete_clock_allowance':implementation,
        'additional_event_allowance':extra,'initial_excitation_tail_allowance':cooling,
        'supplied_event_probability_executed_here':False,
        'numerical_ground_energy_required':False,'rest_frame_or_particle_pole_proved':False})


def ground_energy_bracket(member, lower_receipt=None, upper_receipt=None):
    """Regenerate actual full-physical-source producers before using their bounds."""
    from model.particle_primary_variational_readout import variational
    contract.verify_member(member)
    upper_receipt=variational.bind_member(member) if upper_receipt is None else upper_receipt
    variational.verify_bound(member,upper_receipt)
    exact_upper=F(upper_receipt['result']['ground_energy_upper'])
    grid=1<<32
    upper=F(-(-(exact_upper*grid).numerator//(exact_upper*grid).denominator),grid)
    if lower_receipt is None:
        lower=F(member['preparation_plan']['source']['spectral_lower'])
        lower_producer='primary full-carrier source norm floor'
        lower_digest=contract.digest(member['preparation_plan']['source'])
    else:
        from model.particle_primary_endpoint_bounds import endpoint
        endpoint.verify_bound(member,lower_receipt)
        lower=F(lower_receipt['result']['physical_ground_lower'])
        lower_producer='actual physical endpoint lower bound at coupling one'
        lower_digest=contract.digest(lower_receipt)
    if lower>upper:raise ArithmeticError('authenticated full physical ground bounds are inconsistent')
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'ground_energy_interval':[lower,upper],
        'upper_outward_rounding_error_upper':F(1,grid),
        'selected_mean_source_energy_interval':[lower,upper+F(member['state']['target_ground_excess'])],
        'lower_bound_producer':lower_producer,'lower_producer_receipt_sha256':lower_digest,
        'upper_bound_producer':'actual physical packet Rayleigh quotient at coupling one',
        'upper_producer_receipt_sha256':contract.digest(upper_receipt),
        'producer_results_recomputed':True,
        'trial_replaces_selected_preparation':False,
        'trial_residual_is_an_excitation_or_gap':False,
        'actual_source_eigenvalue_or_histogram_executed':False})


def observed_final_energy_to_excitation(member,final_energy_interval,original_weight_lower,
                                       lower_receipt=None,upper_receipt=None,
                                       resolution='1/64',tail_error='1/1024',clock_error='1/1024',
                                       additional_event_error='0'):
    """Use a verified physical ground bracket for a complete joint-record event."""
    bracket=ground_energy_bracket(member,lower_receipt,upper_receipt)
    plan=clock_plan(member,resolution,tail_error,clock_error)
    lo,hi=(F(x) for x in final_energy_interval)
    if lo>hi:raise ValueError('ordered final-energy record interval required')
    weight=apparatus.exact(original_weight_lower,'original event lower weight')
    extra=apparatus.exact(additional_event_error,'additional complete/sampling event error')
    if max(weight,extra)>1:raise ValueError('normalized original event data required')
    eta=F(resolution)
    inferred=readout.absolute_to_excitation((lo-eta,hi+eta),bracket['ground_energy_interval'])
    tail=F(plan['clock_plan']['complete_clock_tail_upper'])
    implementation=F(plan['two_complete_clock_error_upper'])
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'verified_ground_bracket':bracket,
        'supplied_final_energy_record_interval':[lo,hi],
        'supplied_original_joint_event_lower':weight,
        'final_excitation_interval':inferred['final_excitation_interval'],
        'original_final_excitation_weight_lower':max(F(0),weight-tail-implementation-extra),
        'final_energy_tail_allowance':tail,
        'complete_two_clock_allowance':implementation,
        'additional_event_allowance':extra,
        'initial_cooling_tail_required_for_this_inference':False,
        'supplied_event_probability_executed_here':False,
        'same_branch_momentum_shell_or_rest_mass_proved':False})


def actual_zero_wavevector_band(member):
    """A proved broad actual-source band with factored sign-resolving clock costs."""
    contract.verify_member(member)
    spatial=momentum.common_payload(member);M=spatial['neutral_address_dimension']
    base=member['preparation_plan']['source'];a=F(member['source']['spacing'])
    width=F(base['spectral_upper'])-F(base['spectral_lower'])
    norm=max(abs(F(base['spectral_lower'])),abs(F(base['spectral_upper'])))
    N=spatial['separate_unfiltered_original_blank_weight_each_label_lower']['negative_exponent']
    K=F(4949,20);coefficient=(K/M-1)/(2*a*a*width*width)
    if coefficient<=0:raise ArithmeticError('the exact aperture commutator does not certify this band')
    endpoint_exp=(N+1)//2;tau=1/(4*norm)
    # eta=delta/4, beta=weight/16; J*turn_margin*beta >=1.
    margin_coefficient=tau/(28*a)
    bits=max(1,N+endpoint_exp+apparatus._ceil_log2(16/(coefficient*margin_coefficient)))
    # Each complete repeated-Cayley clock error <=weight/16.
    needed=bits+N+apparatus._ceil_log2(1/(48*coefficient))
    substep_bits=max(0,(needed+1)//2)
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'probe':'unfiltered k=0 coherent aperture mean; same actual primary state',
        'wavevector_index':0,'wavevector_zero_identified_as_rest_frame':False,
        'norm_squared_commutator_on_physical_blank':K/(M*a*a),
        'original_positive_transfer_band_lower':survival.dyadic(1/a,endpoint_exp),
        'original_positive_transfer_band_upper':width,
        'original_positive_band_weight_lower':survival.dyadic(coefficient,N),
        'energy_resolution':survival.dyadic(1/(4*a),endpoint_exp),
        'reported_transfer_lower':survival.dyadic(1/(2*a),endpoint_exp),
        'reported_transfer_upper':'source_width + lower_reported_endpoint',
        'clock_step':tau,'clock_bins':{'base':2,'exponent':bits},
        'Cayley_substeps_per_tick':{'base':2,'exponent':substep_bits},
        'one_energy_tail_allowance':survival.dyadic(coefficient/16,N),
        'one_complete_clock_allowance':survival.dyadic(coefficient/16,N),
        'proved_reported_original_band_weight_lower':survival.dyadic(3*coefficient/4,N),
        'fine_outcomes':{'coefficient':2*M,'base':2,'exponent':2*bits},
        'spatial_filter_used':'none; exact direct algebraic contraction dilation',
        'additional_preparation_compiler_or_sampling_allowance':0,
        'huge_register_or_pulse_tape_allocated':False,
        'broad_support_is_resolved_ridge_or_particle_branch':False,
        'uniform_refinement_weight_or_mass_proved':False})


@lru_cache(maxsize=1)
def exact_control():
    H=s.diag(0,1);rho=s.diag(s.Rational(3,4),s.Rational(1,4))
    X=s.Matrix([[0,1],[1,0]])
    U=s.diag(1,-s.I) # exact exp(-i*pi*H/2); signed centers(0,1,-2,-1).
    energy=energy_instrument(U,4)
    spatial=momentum.fourier_instrument((X/2,-X/2))
    rows=joint_instrument(energy,spatial['rows'])
    for row in rows:row['p']=s.simplify(s.trace(rho*row['K'].H*row['K']))
    root=apparatus.positive_sqrt(rho);omega=s.Matrix([root[i,j] for j in range(2) for i in range(2)])
    joint_before=omega*omega.H
    joint_after=sum((s.kronecker_product(s.eye(2),f)*joint_before*s.kronecker_product(s.eye(2),f.H)
                     for f in energy),s.zeros(4))
    return {'H':H,'rho':rho,'source_step':U,'clock_step':'pi/2',
        'energy_rows':energy,'spatial_rows':spatial['rows'],'rows':rows,
        'initial_reference_joint_before':joint_before,'initial_reference_joint_after':joint_after,
        'complete_original_probability_sum':sum(row['p'] for row in rows),
        'is_primary_source_execution':False}


def bind_member(member,**kwargs):
    result=common_payload(member,**kwargs)
    return {'schema':'oph.primary_joint_energy_spatial_records.v1','result':result,
            'result_sha256':contract.digest(result)}


def verify_bound(member,receipt,**kwargs):
    if receipt!=bind_member(member,**kwargs):raise ValueError('joint spectroscopy receipt differs from source or domain')


def certificate():
    from model.particle_primary_endpoint_bounds import endpoint
    member=contract.member(2);control=exact_control()
    return contract.serial({'schema':'oph.primary_joint_spectroscopy.v1',
        'actual_source_instrument':bind_member(member),
        'proved_broad_actual_source_band':actual_zero_wavevector_band(member),
        'authenticated_physical_ground_bracket_same_aperture_member':ground_energy_bracket(member,endpoint.bind_member(member)),
        'independent_exact_control':{'all_records':len(control['rows']),
            'nonzero_records':[{k:(str(v) if k=='p' else v) for k,v in row.items() if k!='K'}
                               for row in control['rows'] if row['p']!=0],
            'probability_sum':str(control['complete_original_probability_sum']),
            'initial_energy_step_changes_joint_purification':True,'is_primary':False},
        'actual_primary_resolved_branch_or_mass_GeV':None})
