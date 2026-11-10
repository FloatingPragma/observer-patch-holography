"""Source-only finite-window packet collision with all errors below signal."""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
RER=HERE.parents[1]
from model.interacting_block_quadrature.model import ceil_fraction,ceil_sqrt_fraction,dec
import mpmath as mp
import numpy as np
from scipy import sparse
from .finite_collision import certificate as finite_collision


PARENT_QUADRATURE='model/interacting_block_quadrature/results.json'
PARENT_BLOCKS='model/interacting_block_propagation/results.json'
PARENTS=[PARENT_QUADRATURE,PARENT_BLOCKS]
NATIVE_PREPARATION_PARENT='model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json'


def source():
    parents={};pins={}
    for receipt in PARENTS:
        path=RER/receipt
        value=json.loads(path.read_text())
        for name,pin in value['input_sha256'].items():
            relative=Path(name)
            if relative.is_absolute() or '..' in relative.parts:
                raise ValueError('parent path must be repository relative: '+name)
            source_path=RER/relative
            if not source_path.is_file() or hashlib.sha256(source_path.read_bytes()).hexdigest()!=pin:
                raise ValueError('immutable parent source mismatch: '+name)
            pins[name]=pin
        pins[receipt]=hashlib.sha256(path.read_bytes()).hexdigest()
        parents[receipt]=value
    if NATIVE_PREPARATION_PARENT not in pins:
        raise ValueError('missing native deterministic-trial ancestry')
    native=json.loads((RER/NATIVE_PREPARATION_PARENT).read_text())['preparation']
    if native['actual_local_preparation']['preparation_success_probability']!='1':
        raise ValueError('initial native local trial is not deterministic')
    p0=F(native['atomic_cube_trial']['exact_interacting_ground_trial_overlap_lower'])
    if F(parents[PARENT_BLOCKS]['preparation_scaling']['initial_excited_odds_upper'])!=(1-p0)/p0:
        raise ValueError('filter tree must start from the unfiltered trial')
    if F(native['finite_energy_agreement_filter']['success_probability_lower'])!=p0*p0:
        raise ValueError('first agreement-filter probability has different ancestry')
    symmetry=parents[PARENT_QUADRATURE]['retained_native_symmetry_reference']
    if symmetry['group']['native_full_cover_order']!=3072:
        raise ValueError('native symmetry cover changed')
    if not symmetry['completion']['complete_record_probabilities_and_norm_errors_preserved']:
        raise ValueError('complete regular-reference lift is required')
    return parents,pins


def inputs():
    parents,pins=source()
    return parents[PARENT_BLOCKS],parents[PARENT_QUADRATURE],pins[PARENT_QUADRATURE]


def stirling(n,k):
    if n==k==0:return 1
    if n==0 or k==0 or k>n:return 0
    return k*stirling(n-1,k)+stirling(n-1,k-1)


def moment_bound(a,V):
    # Exact domination: n+6 <= in distribution 8+2*NB((3V+4)/2,q).
    u=(a+1/a-2)/4;r=F(3*V+4,2)
    total=F(0)
    for j in range(9):
        raw=F(0)
        for k in range(j+1):
            rising=math.prod(r+i for i in range(k))
            raw+=stirling(j,k)*rising*u**k
        total+=math.comb(8,j)*8**(8-j)*2**j*raw
    return total


def menus(side):
    I=[(2,1,0),(2,-1,0)]
    O=[(0,1,2),(0,-1,2)]
    minus=lambda v:tuple(-x for x in v)
    conserved=[]
    for i,j,k,l in itertools.product(I,list(map(minus,I)),O,list(map(minus,O))):
        if all((i[a]+j[a]-k[a]-l[a])%side==0 for a in range(3)):
            conserved.append([i,j,k,l])
    assert len(conserved)==6
    assert set(I).isdisjoint(O)
    return conserved


def preparation(parent,N,V,amplitude):
    prep=parent['preparation_scaling']
    r0=F(prep['initial_excited_odds_upper'])
    c=F(prep['geometric_odds_contraction'])
    count=N*V;depth=0;odds=r0
    target=(amplitude/8)**2
    while 2*count*odds>target:
        odds*=c;depth+=1
    log_success=mp.mpf(0)
    for j in range(depth):
        log_success-=count*2**(depth-j)*mp.log1p(dec(r0*c**j))
    return {'constituent_native_blocks':count,'physical_source_cells':8*count,
            'deterministic_initial_trial':True,'base_trial_success':'1',
            'initial_state':'unfiltered atomic/cube trial',
            'filter_depth_counts_all_rounds':True,
            'binary_filter_depth':depth,'trials_per_constituent':2**depth,
            'retained_purification_vector_error_squared_upper_factored':'2*N*V*r0*c^d',
            'retained_purification_vector_error_squared_upper_decimal':str(dec(2*count*odds)),
            'target_vector_error_squared_exact':str(target),
            'same_frozen_filter_M':prep['same_frozen_filter_M'],
            'initial_odds_upper_exact':str(r0),'odds_contraction_exact':str(c),
            'success_lower_factored':'product_(j=0..d-1)(1+r0*c^j)^(-N*V*2^(d-j))',
            'log_success_lower_decimal':str(log_success),
            'decimal_logs_are_diagnostics_exact_success_bound_is_factored':True,
            'all_preparation_flags_retained':True}


def complete_time_series_preparation(direct):
    from .finite_collision import exact_endpoint
    iv=mp.iv
    prep=direct['complete_native_preparation']
    squared=2*prep['constituent_native_blocks']*F(prep['initial_odds_upper_exact'])*F(prep['odds_contraction_exact'])**prep['binary_filter_depth']
    error=iv.sqrt(iv.mpf(squared.numerator)/squared.denominator)
    def number(value):
        value=F(value)
        return iv.mpf(value.numerator)/value.denominator
    for row in direct['complete_time_series']:
        low,high=map(F,row['ideal_hit_probability_interval_exact'])
        prepared_low=max(F(0),exact_endpoint(iv.sqrt(number(low))-error))**2
        prepared_high=min(F(1),max(F(0),exact_endpoint(iv.sqrt(number(high))+error,True))**2)
        row['native_preparation_hit_probability_interval_exact']=[str(prepared_low),str(prepared_high)]
        row['native_preparation_complement_probability_interval_exact']=[str(1-prepared_high),str(1-prepared_low)]
        row['native_preparation_vector_error_upper_exact']=str(exact_endpoint(error,True))
        free_low,free_high=map(F,row['same_admitted_free_hit_probability_interval_exact'])
        free_prepared_low=max(F(0),exact_endpoint(iv.sqrt(number(free_low))-error))**2
        free_prepared_high=min(F(1),max(F(0),exact_endpoint(iv.sqrt(number(free_high))+error,True))**2)
        row['same_admitted_free_preparation_error_hit_interval_exact']=[str(free_prepared_low),str(free_prepared_high)]
        row['same_admitted_free_preparation_error_complement_interval_exact']=[str(1-free_prepared_high),str(1-free_prepared_low)]
        row['native_minus_matched_free_hit_difference_interval_exact']=[str(prepared_low-free_prepared_high),str(prepared_high-free_prepared_low)]
        if prep['binary_filter_depth']==0:
            row['preparation_success_probability_exact']='1'
            row['preparation_failure_probability_exact']='0'
            row['unconditional_hit_probability_interval_exact']=[str(prepared_low),str(prepared_high)]
            row['unconditional_complement_probability_interval_exact']=[str(1-prepared_high),str(1-prepared_low)]
        for site in row['occupation_marginals']:
            site['component_record_kind']='joint reference-relative native projectors'
            lower,upper=map(F,site['ground_ideal_interval_exact'])
            site['ground_native_preparation_interval_exact']=[
                str(max(F(0),exact_endpoint(number(lower)-error))),
                str(min(F(1),exact_endpoint(number(upper)+error,True)))]
            site['higher_native_sector_probability_upper_exact']=str(squared)
            for color in site['components']:
                lower,upper=map(F,color['ideal_interval_exact'])
                color['native_preparation_interval_exact']=[
                    str(max(F(0),exact_endpoint(number(lower)-error))),
                    str(min(F(1),exact_endpoint(number(upper)+error,True)))]
            twirled_lower=sum(F(color['ideal_interval_exact'][0]) for color in site['components'])/3
            twirled_upper=sum(F(color['ideal_interval_exact'][1]) for color in site['components'])/3
            prepared_twirled_lower=sum(F(color['native_preparation_interval_exact'][0]) for color in site['components'])/3
            prepared_twirled_upper=sum(F(color['native_preparation_interval_exact'][1]) for color in site['components'])/3
            site['reference_traced_source_components']=[{
                'component':color,'ideal_interval_exact':[str(twirled_lower),str(twirled_upper)],
                'display_midpoint':float((twirled_lower+twirled_upper)/2),
                'native_preparation_interval_exact':[str(prepared_twirled_lower),str(prepared_twirled_upper)]}
                for color in range(3)]


def member(a,side,parent):
    a=F(a);V=side**3;z=6;M=3*V
    D=1/a;t=(1-a*a)/(12*a);T=F(1)
    pi_hi=F(22,7)
    A_lo=6-4*pi_hi*pi_hi*5/side**2
    assert 0<A_lo<6
    omega2_hi=D*D-2*D*t*A_lo
    single_vertex_lo=t*A_lo/V*(D*D/omega2_hi+1)
    packet_vertex_lo=F(3,2)*single_vertex_lo
    M8=moment_bound(a,V);B=ceil_sqrt_fraction(M8)
    C=12*z*V*t
    G=16*M*M/a**4+16*M/a**2+4
    R_B=B*(T*C+T*T*C*C*G/2)
    N=max(2,ceil_fraction(8*R_B/(T*packet_vertex_lo)))
    while 8*F(B,(N+7)**4)>T*packet_vertex_lo/(8*N) or F(B,(N+7)**4)>F(1,4):N*=2
    amplitude=T*packet_vertex_lo/N
    born=R_B/N**2;projection=8*F(B,(N+7)**4)
    assert born<=amplitude/8 and projection<=amplitude/8
    prep=preparation(parent,N,V,amplitude)
    probability=(amplitude/2)**2
    energy_halfwidth=F(1,4)
    energy_residual=3*C*B/N
    energy_tail=(energy_residual/energy_halfwidth)**2
    prepared_energy_tail=min(F(1),(energy_residual/energy_halfwidth+amplitude/8)**2)
    # The scalar packet frequencies have the same fixed physical menu.
    actual_A=2*(mp.cos(4*mp.pi/side)+mp.cos(2*mp.pi/side)+1)
    omega=mp.sqrt(dec(D*D)-2*dec(D*t)*actual_A)
    p=[4*mp.pi/dec(a*side),2*mp.pi/dec(a*side),mp.mpf(0)]
    kg=mp.sqrt(1+sum(x*x for x in p)/6)
    dispersion_bound=dec(a*a)*(sum(x*x for x in p)/12+sum(x**4 for x in p)/144)
    assert abs(omega-kg)<=dispersion_bound
    actual_vertex=mp.mpf('1.5')*dec(t)/V*actual_A*(dec(D*D)/(omega*omega)+1)
    assert dec(packet_vertex_lo)<=actual_vertex
    return {'spacing':str(a),'side':side,'macro_sites':V,'physical_box_side':str(a*side),
            'time_Delta_times_T':'1','constituent_blocks_per_macrocell_N':N,
            'source_parameters':{'D_over_Delta':str(D),'t_over_Delta':str(t)},
            'packet_menu':{'incoming_alpha':[[2,1,0],[2,-1,0]],'incoming_beta':'negative incoming_alpha',
                           'outgoing_alpha':[[0,1,2],[0,-1,2]],'outgoing_beta':'negative outgoing_alpha',
                           'coordinate_unit':'2*pi/side','each_packet_amplitude':'1/sqrt(2)',
                           'momentum_conserving_quadruples':menus(side),
                           'free_separate_component_momentum_forbids_channel_exactly':True,
                           'physical_momenta_decimal':[str(x) for x in p]},
            'exact_bounds':{'adjacency_shell_eigenvalue_lower':str(A_lo),
                            'packet_vertex_over_Delta_lower':str(packet_vertex_lo),
                            'eighth_number_moment_upper':str(M8),'moment_norm_B_upper':B,
                            'weighted_difference_constant_C_delta_over_Delta':str(C),
                            'free_weighted_evolution_constant_G':str(G),
                            'Born_remainder_numerator_R_B':str(R_B),
                            'leading_collision_amplitude_lower':str(amplitude),
                            'second_order_Born_amplitude_error_upper':str(born),
                            'capacity_projection_and_detector_amplitude_error_upper':str(projection),
                            'conditional_hit_probability_lower':str(probability),
                            'energy_halfwidth_over_Delta':str(energy_halfwidth),
                            'exact_admitted_packet_energy_residual_over_Delta_upper':str(energy_residual),
                            'exact_admitted_packet_probability_outside_actual_energy_window_upper':str(energy_tail),
                            'native_prepared_packet_probability_outside_actual_energy_window_upper':str(prepared_energy_tail)},
            'decimals':{'omega_over_Delta':str(omega),'KG_omega_over_Delta':str(kg),
                         'KG_frequency_error_bound_over_Delta':str(dispersion_bound),
                         'actual_packet_vertex_over_Delta':str(actual_vertex),
                         'conditional_hit_probability_lower':str(dec(probability)),
                         'exact_admitted_packet_energy_window_tail_upper':str(dec(energy_tail)),
                         'native_prepared_packet_energy_window_tail_upper':str(dec(prepared_energy_tail)),
                         'log_unconditional_hit_probability_lower':str(mp.mpf(prep['log_success_lower_decimal'])+mp.log(dec(probability)))},
            'complete_preparation':prep}


def exact_two_particle_matrix(N,side=5):
    """Direct positional matrix, two colors, independent of Fourier formula."""
    V=side;t=F(1,8)
    H=np.zeros((V*V,V*V))
    free=np.zeros_like(H)
    for x,y in itertools.product(range(V),repeat=2):
        origin=x*V+y
        for component in (0,1):
            pos=(x,y)[component]
            other=(y,x)[component]
            for sign in (-1,1):
                target=(pos+sign)%V
                # Shared capacity matters if the spectator is at either end.
                factor=math.sqrt(1-1/N) if other in (pos,target) else 1.
                pair=[x,y];pair[component]=target
                dest=pair[0]*V+pair[1]
                H[dest,origin]-=float(t)*factor
                free[dest,origin]-=float(t)
    assert np.max(abs(H-H.T))<1e-15
    k=2*np.pi/side
    coordinates=np.asarray(list(itertools.product(range(V),repeat=2)))
    incoming=np.exp(1j*k*(coordinates[:,0]-coordinates[:,1]))/V
    outgoing=np.exp(-1j*k*(coordinates[:,0]-coordinates[:,1]))/V
    element=np.vdot(outgoing,(H-free)@incoming)
    expected=float(t)*(1-math.sqrt(1-1/N))/V*8*math.cos(k)
    assert abs(element-expected)<1e-14
    assert abs(np.vdot(outgoing,free@incoming))<1e-14
    return {'N':N,'side':side,'two_color_position_dimension':V*V,
            'exact_capacity_effect_matrix_element':float(element.real),
            'closed_form_matrix_element':expected,'free_matrix_element_zero':True}


def independent_bog_vertex(side=5):
    """Evaluate the quartic cross term via finite quasiparticle CAR-free matrices."""
    # One component, beta occupation <=2 is enough for the two quadratic factors.
    basis=[n for n in itertools.product(range(3),repeat=side) if sum(n)<=2]
    ids={n:i for i,n in enumerate(basis)}
    lowers=[]
    for j in range(side):
        b=np.zeros((len(basis),len(basis)),complex)
        for col,n in enumerate(basis):
            if n[j]:
                target=list(n);target[j]-=1;b[ids[tuple(target)],col]=math.sqrt(n[j])
        lowers.append(b)
    D=2.;t=.125
    # A cycle has degree two; use its own stable exact dispersion here.
    frequencies=np.sqrt(D*D-4*D*t*np.cos(2*np.pi*np.arange(side)/side))
    u=(np.sqrt(D/frequencies)+np.sqrt(frequencies/D))/2
    v=(np.sqrt(D/frequencies)-np.sqrt(frequencies/D))/2
    bsites=[]
    for x in range(side):
        b=sum((np.exp(2j*np.pi*j*x/side)*(u[j]*lowers[j]+v[j]*lowers[(-j)%side].T.conj()) for j in range(side)),np.zeros_like(lowers[0]))/np.sqrt(side)
        bsites.append(b)
    q=[b+b.T.conj() for b in bsites]
    num=[b.T.conj()@b for b in bsites]
    def state(j):
        n=[0]*side;n[j]=1
        s=np.zeros(len(basis));s[ids[tuple(n)]]=1
        return s
    ia,ib,oa,ob=state(1),state(side-1),state(side-1),state(1)
    value=0j
    for x in range(side):
        y=(x+1)%side
        value+=t/2*(np.vdot(oa,q[x]@q[y]@ia)*np.vdot(ob,(num[x]+num[y])@ib)
                    +np.vdot(ob,q[x]@q[y]@ib)*np.vdot(oa,(num[x]+num[y])@ia))
    A=2*math.cos(2*np.pi/side)
    expected=t*A/side*((D/frequencies[1])**2+1)
    assert abs(value-expected)<1e-13
    return {'side':side,'one_color_quasiparticle_basis_dimension':len(basis),
            'quartic_Bogoliubov_matrix_element':float(value.real),
            'closed_form_matrix_element':expected,'same_free_total_energy':True}


@lru_cache(maxsize=1)
def calculate():
    mp.mp.dps=80
    parent,q,sha=inputs()
    members=[member(a,side,parent) for a,side in ((F(1,4),32),(F(1,8),64))]
    direct=finite_collision()
    direct['retained_symmetry_reference']={
        'native_cover':'Z4^3 semidirect binary octahedral group',
        'native_cover_order':3072,'macro_group':'C5 ring translations','macro_group_order':5,
        'joint_group':'native cover direct product C5',
        'retained_regular_reference_dimension':15360,
        'source_reference_clock':'H_cycle tensor I; H_reference=0',
        'three_dimensional_macro_cubic_action_claim':False,
        'record_marginals_are_joint_reference_relative':True,
        'hit_detector_record_kind':'joint reference-relative target projector',
        'occupation_display_toggle_changes_detector_protocol':False,
        'reference_traced_source_is_group_twirl':True,
        'additional_binary_detector_dimension':2,
        'reference_and_detector_dimension':30720,
        'additional_capacity_nats':'log(15360)+log(2)',
        'native_preparation_references_are_retained_separately':True}
    direct['matched_comparator_full_carrier']={
        'common_carrier':'Fock(15 modes) direct sum complete native H_rest, tensor all retained references',
        'physical_embedding':'canonical capacity embedding on H_low; identity on H_rest',
        'free_comparator_unitary':'U_free direct sum actual U_native restricted to H_rest',
        'native_extension':'actual U_native on embedded physical carrier; identity on unused Fock complement',
        'detector':'same admitted target projector, zero on H_rest, and its full complement',
        'preparation_vector_error_preserved_by_isometry':True,
        'all_native_and_oscillator_complement_outcomes_retained':True,
        'additional_oscillators_are_native_hardware':False}
    direct['logical_capacity']={
        'native_blocks':5,'native_cells':40,
        'one_block_kinematic_dimension':'(5*2^30)^8*119^24',
        'full_source_kinematic_capacity_nats':'5*(240*log(2)+8*log(5)+24*log(119))',
        'numerically_executed_reducing_dimension':1024,
        'numerically_executed_reducing_capacity_nats':'log(1024)',
        'full_native_matrix_numerically_enumerated':False,
        'native_complement_control':'exact reducing-subspace identity and full-source preparation-vector bound',
        'native_trial_spin_reference_capacity_nats':'40*log(4)',
        'retained_regular_reference_capacity_nats':'log(15360)',
        'binary_hit_complement_pointer_capacity_nats':'log(2)',
        'filter_registers_used_at_depth_zero':0,
        'separate_eight_character_register_used':False,
        'total_declared_capacity_nats':'5*(240*log(2)+8*log(5)+24*log(119))+40*log(4)+log(15360)+log(2)',
        'comparison_Fock_space_is_part_of_native_hardware':False}
    direct['complete_native_preparation']=preparation(parent,1,5,F(direct['exact_source_amplitude_lower']))
    complete_time_series_preparation(direct)
    direct['log_unconditional_hit_probability_lower']=str(
        mp.mpf(direct['complete_native_preparation']['log_success_lower_decimal'])
        +mp.log(dec(F(direct['conditional_hit_probability_lower_after_preparation_margin']))))
    result={'schema':'oph.interacting_block_scattering.v1',
            'same_public_quadrature_receipt_sha256':sha,
            'validated_deterministic_trial_ancestor':NATIVE_PREPARATION_PARENT,
            'source_input_sha256':q['input_sha256'],
            'direct_capacity_matrix_checks':[exact_two_particle_matrix(N) for N in (2,3,8,64)],
            'independent_quartic_vertex':independent_bog_vertex(),
            'interval_certified_five_block_collision':direct,
            'finite_members':members,
            'scope':{'same_paired_native_source_and_clock':True,
                     'prepared_two_packet_state_is_actual_interacting_ground':False,
                     'all_preparation_and_detector_records_retained':True,
                     'amplitude_errors_bounded_below_nonzero_vertex':True,
                     'fixed_finite_window_nonzero_interaction_proved':True,
                     'N1_same_admitted_free_channel_forbidden':False,
                     'matched_full_free_preparation_inclusive_response_difference_proved':True,
                     'asymptotic_scattering_matrix_proved':False,
                     'nontrivial_interacting_relativistic_UV_limit_proved':False,
                     'measured_inputs':[]}}
    return result
