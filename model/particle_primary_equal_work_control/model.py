"""Same selected source, exact equal work, and a distinct complete apparatus."""
from fractions import Fraction as F
from math import factorial

from model.particle_common_contract import contract, provenance
from model.particle_primary_critical_refinement import model as critical
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS=(
    'model/particle_primary_equal_work_control/__init__.py',
    'model/particle_primary_equal_work_control/__main__.py',
    'model/particle_primary_equal_work_control/model.py',
    'model/particle_primary_equal_work_control/test_model.py',
    'model/particle_primary_equal_work_control/PROOF.md',
    'model/particle_primary_equal_work_control/README.md',
    'model/particle_primary_critical_refinement/PROOF.md',
    'model/particle_primary_weighted_spectral_tail/PROOF.md',
    'model/particle_primary_current_tail/PROOF.md',
    'model/particle_primary_joint_spectroscopy/PROOF.md',
    'model/particle_common_contract/PROOF.md',
)


def _atan_interval(denominator,terms=12):
    x=F(1,denominator)
    partial=sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(terms)),F())
    other=partial+(-1)**terms*x**(2*terms+1)/F(2*terms+1)
    return min(partial,other),max(partial,other)


def pi_interval():
    """Machin's identity and alternating rational remainders."""
    a,b=_atan_interval(5);c,d=_atan_interval(239)
    return 16*a-4*d,16*b-4*c


def _dyadic_bounds(lo,hi,bits=48):
    q=2**bits
    return F((lo*q).numerator//(lo*q).denominator,q),F(-((-hi*q).numerator//(-hi*q).denominator),q)


def coin_interval(level):
    """Outward bounds on the exact algebraic coin, without enumerating V cells."""
    n=ground._i(level,'source level',0);R=2**n
    if n==0:return F(1),F(1)
    plo,phi=pi_interval();lo=F();hi=F()
    # 1-cos x lies between its six-term and five-term alternating sums.
    # x<=pi/R^2<=pi/4<1 for every n>=1.
    for k in range(1,7):
        power_sum=R**(2*k)+2*sum(u**(2*k) for u in range(1,R))
        coefficient=F(R*power_sum**2,2*R**(8*k)*factorial(2*k))
        if k%2:
            lo+=coefficient*plo**(2*k);hi+=coefficient*phi**(2*k)
        else:
            lo-=coefficient*phi**(2*k)
            if k<6:hi-=coefficient*plo**(2*k)
    return _dyadic_bounds(lo,hi)


def coin_recipe(level,widths=(1,1,1),amplitude=1):
    n=ground._i(level,'source level',0)
    profile=critical.profile_geometry(n,widths,amplitude)
    a=profile['spacing'];d2=profile['three_axis_gradient_square_sum']
    upper=F(5,12)*F(amplitude)**2*a*a*d2
    return {'level':n,'spacing':a,'physical_half_widths':list(widths),
        'amplitude':F(amplitude),
        'exact_probability':'p=(1/12) sum_(i,z) [1-cos(pi*s*a*(f(z+e_i)-f(z)))]',
        'algebraic_phase_order_divides':profile['phase_root_of_unity_order_divides'],
        'probability_upper':upper,
        'probability_admitted_by_bound':upper<=1,
        'same_profile_work_exactly_matched_when_probability_at_most_one':True}


def bind_member(member):
    contract.verify_member(member)
    source=member['source'];n=source['level'];a=F(source['spacing']);V=source['side']**3
    inherited=critical.bind_member(member);value=inherited['result'];profile=value['profiles'][0]
    lo,hi=coin_interval(n)
    envelope=F(1) if n==0 else F(45,32)*a
    if not 0<lo<=hi<=envelope<=1:raise ArithmeticError('exact parity coin left its source domain')
    wlo=F(profile['original_signed_work_lower']);whi=F(profile['original_signed_work_upper'])
    uniform=F(profile['uniform_original_signed_work_lower'])
    # Each parity changes only its incident hops; max degree is3 at n=0,6 otherwise.
    degree=3 if n==0 else 6
    m2lo=wlo*wlo/hi;m2hi=3600*degree**2*hi/(a*a)
    uniform_m2lo=F(32,45)*uniform**2/a
    uniform_m2hi=F(182250)/a
    if not uniform_m2lo<=m2lo<=m2hi<=uniform_m2hi:
        raise ArithmeticError('same-source two-sided moment estimates failed')
    clock=value['complete_two_energy_clock'];error=F(value['complete_energy_clock_signed_mean_error_upper'])
    clockloss=2*F(clock['clock_plan']['complete_clock_tail_upper'])+F(clock['two_complete_clock_error_upper'])
    if error>=wlo:raise ArithmeticError('unchanged clock lost positive work')
    result=contract.serial({'schema':'oph.primary_equal_work_control_result.v1',
        'member_id':member['member_id'],**member['ids'],
        'source_critical_bound_sha256':inherited['result_sha256'],
        'source_state_reference_H_and_clock_unchanged':True,
        'apparatus_identical_to_critical_phase_menu':False,
        'apparatus':'K_(x,0)=sqrt((1-p)/V)*I; K_(x,1)=sqrt(p/V)*(-1)^N_x',
        'retained_origin_coin_outcomes':2*V,
        'exact_coin':'p=kbar/4; kbar=2R sum_(u,v=1..R) m_u m_v [1-cos(pi*u*v/R^4)]',
        'transverse_multiplicity':'m_R=1, m_u=2 for1<=u<R; R=2^n',
        'exact_coin_algebraic_root_order':2*(2**n)**4,
        'exact_coin_outward_interval':[lo,hi],
        'uniform_coin_upper':envelope,
        'coin_contains_fitted_mass_or_selected_expectation':False,
        'exact_selected_mean_identity':'W_control=W_critical=-kbar*Tr(rho H_hop)/V',
        'original_work_interval':[wlo,whi],
        'uniform_original_work_lower':uniform,
        'original_second_moment_interval':[m2lo,m2hi],
        'uniform_second_moment_interval':[uniform_m2lo,uniform_m2hi],
        'second_moment_order':'Theta(1/a) on this distinct actual-source apparatus',
        'original_nonelastic_probability_upper':hi,
        'ideal_transfer_measure_formula':'mu=(1-p)*delta_0+p*nu_parity',
        'ideal_transfer_limit':'delta_0',
        'positive_energy_weighted_tail_lower_formula':'integral_(omega>R) omega dmu >= W_critical-R*p',
        'uniform_energy_weighted_integrability':False,
        'complete_CQ_output_distance_from_matching_idle_experiment_upper':hi,
        'coherent_unread_pointer_distance_upper_formula':'min(1,2*sqrt(p)); pointer coherence is not discarded in this comparison',
        'same_critical_complete_clock_sha256':contract.digest(clock),
        'complete_clock_mean_error_upper':error,
        'complete_reported_work_lower':wlo-error,
        'energy_resolution':clock['energy_resolution'],
        'complete_clock_probability_error_upper':clockloss,
        'fixed_positive_reported_band_upper_formula':'p+clock_error when lower_endpoint>2*resolution',
        'all_fine_clock_origin_coin_records':2*V*clock['clock_plan']['clock_size']**2,
        'critical_phase_collapse_proved':False,
        'original_critical_current_covariance_evaluated':False,
        'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_equal_work_control_bound.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('equal-work source or apparatus certificate changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*critical.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_equal_work_control_certificate.v1',
        'members':[bind_member(contract.member(n)) for n in (0,1,2,8)],
        'source_sha256':pins,'physical_mass_table_GeV':[]}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('equal-work inputs changed during replay')
    return result
