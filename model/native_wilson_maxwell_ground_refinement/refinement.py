"""Uniform fixed-cube interacting ground and complete photon-cap refinement."""
from fractions import Fraction as F
from math import factorial,comb
from pathlib import Path
import hashlib,json,sys
import mpmath as mp

HERE=Path(__file__).resolve().parent
RER=HERE.parents[1]
sys.path.insert(0,str(RER))
from model.native_wilson_maxwell_ground import ground
from model.native_wilson_maxwell_refinement import bounds
from model.affine_gauss_mobile.model import rounded

PARENT='model/native_wilson_maxwell_ground/results.json'
PARENT_SHA='0e9edccffa2d8571f0b325feaee7c49408714d7b838e3e1c57ec6a183c15fa17'
CUTOFF_PARENT='model/native_wilson_maxwell_refinement/results.json'
CUTOFF_SHA='49da49cb0745d1f78563f50fd4f6829e45e957d2bec7e4a87064075ce0a202ec'


def admit():
    path=RER/PARENT
    assert hashlib.sha256(path.read_bytes()).hexdigest()==PARENT_SHA
    data=json.loads(path.read_text())
    for name,h in data['input_sha256'].items():
        assert hashlib.sha256((RER/name).read_bytes()).hexdigest()==h,name
    cutoff=RER/CUTOFF_PARENT
    assert hashlib.sha256(cutoff.read_bytes()).hexdigest()==CUTOFF_SHA
    for name,h in json.loads(cutoff.read_text())['input_sha256'].items():
        assert hashlib.sha256((RER/name).read_bytes()).hexdigest()==h,name
    return data,ground.admit()


def ceil(value):return -(-value.numerator//value.denominator)


def exp_interval(value):
    """Outward rational range reduction; no floating transcendental input."""
    if value<0:raise ValueError('nonnegative argument required')
    shifts=0
    while value>1:value/=2;shifts+=1
    lo=bounds.rounded(value,False,90);hi=bounds.rounded(value,True,90)
    tl=th=sl=sh=F(1)
    for n in range(1,257):
        tl=bounds.rounded(tl*lo/n,False,90)
        th=bounds.rounded(th*hi/n,True,90)
        sl+=tl;sh+=th
    sh+=th*hi/F(257)/(1-hi/F(258))
    for _ in range(shifts):
        sl=bounds.rounded(sl*sl,False,90)
        sh=bounds.rounded(sh*sh,True,90)
    return sl,sh


def exp_lower(value):return exp_interval(value)[0]
def exp_upper(value):return exp_interval(value)[1]


def scientific():
    parent,source=admit();old=parent['bounds_exact'];g=F(1,64)
    coeffs=parent['exact_vacuum_contractions']
    cl,cu=map(F,coeffs['Coulomb_vacuum_coefficient_interval'])
    vl=F(1,16)+1/(12*bounds.root(6));vu=F(1,16)+1/(12*bounds.root(6,False))
    sigma=source['independent_current_coefficient']['uniform_polynomial_method']
    sigma_up=bounds.root(F(sigma['pair_one_photon_V1_vacuum_norm_squared_over_m_squared_interval'][1]))
    sigma_low=F(sigma['pair_one_photon_V1_vacuum_norm_over_m_lower'])
    muL=F(old['free_vacuum_shift_over_m_lower']);muU=F(old['free_vacuum_shift_over_m_upper'])
    c=F(old['counterterm_over_m']);mlo,mhi=map(F,source['clock_and_capacity']['source_selected_excitation_interval'])
    A2=96**2*g*g*vu;eps=F(1,10);b=A2/(8*eps)+2*eps
    d0=2*(1-eps)-b;first=(1-eps)*(2*bounds.root(3,False)+1)-b
    Rv=24*bounds.root(3)*vu+bounds.root(8*cu)
    cstar=F(1,10**6)
    weighted=(g*sigma_up/bounds.root(first-cstar,False)+g*g*Rv/bounds.root(d0-cstar,False))**2
    margin=muL-cstar-weighted
    assert margin>F(7,10**6)
    gamma=d0-muU
    theta=(g*sigma_up+g*g*Rv)/gamma
    gap_native=F(parent['complete_native_complement']['actual_neutral_triplet_gap_interval'][0])
    full_gap=min(mlo*gamma,gap_native+mlo*(c-muU))
    assert full_gap>0
    onsite=4*bounds.root(3,False)+8*bounds.root(7,False)+4*bounds.root(13,False)-40
    assert onsite>0
    B=F(19397);native_norm=480+F(360,2**34)
    raw_ground_bound=32*native_norm+mhi*(3+bounds.root(6)+muU)
    assert raw_ground_bound<B
    M=ceil(4*B/(F(1,1024)*full_gap));epsilon=4*B/(M*full_gap)
    pcode=F(source['actual_preparation']['global_code_probability_lower'])
    p=bounds.rounded(pcode/(1+theta*theta),False,60)
    odds=(1-p)/p;qtree=F(1);steps=[]
    for level in range(2):
        pj=1/(1+odds);qtree*=pj**(2**(2-level))
        odds=odds*odds+2*epsilon*epsilon*odds
        steps.append({'level':level+1,'source_ground_probability_lower':str(bounds.rounded(pj,False,60)),
            'next_ground_odds_upper':str(bounds.rounded(odds,True,60))})
    fidelity=1/(1+odds);trace=bounds.root(1-fidelity)
    pur=bounds.root(2*(1-bounds.root(fidelity,False)))
    odd=4*bounds.root(15)*vu*bounds.root(vu)*g**3
    target_V=bounds.root(2*A2)+8*g*g
    amplitude=((g*sigma_low-odd)/bounds.root(1+theta*theta)-target_V*theta)/(11+muU)
    assert amplitude>pur>0
    r=bounds.displacement_radius(g)
    # Heat-time comparison is a proof device, not a physical clock.
    heat_s=F(12);real_time=F(1,64);radius=F(1,64)*(bounds.root(F(1,2))+bounds.root(2))
    heat_K=1024;J=1700;x=48*(heat_s+real_time);R=J*r+radius
    assert 3*R*R<=F(heat_K+1,4)
    prefactor=1+bounds.root(heat_K+1)*R
    heat_word= exp_upper(x+R*R)*prefactor/F(2**(heat_K+1))
    heat_tail=bounds.dyson_tail(x,J)
    heat_delta=heat_word+heat_tail
    heat_ground_error=2*theta/exp_lower(gamma*heat_s)
    heat_record=heat_ground_error+2*bounds.root(1+theta*theta)*exp_upper(muU*heat_s)*heat_delta
    assert heat_record<F(1,10**8)
    # Same finite Cayley pointer and complete tree record alphabet at all caps.
    K=2**96;L=K//2
    threshold=ceil(4*B/mlo)
    assert K>=threshold and r<F(1,200)
    neumann=48*mhi/B
    assert neumann<1
    Cayley_constant=2/(1-neumann)*(48*mhi/B/F(threshold+1)+192)
    assert Cayley_constant<400
    dC=F(400,K+1)
    # Move the three common source-free filter factors before the filters;
    # every filter commutes with every individual source Hamiltonian.
    tau=3*mhi
    dU=48*tau/F((K+1)**2)+2304*r*tau*tau*bounds.root(F(2,K+1))
    full_tree_error=4*dU+6*(M-1)*dC
    assert full_tree_error<F(1,10**8)
    # Mean photon number after the complete clock program is bounded by
    # conserved source energy, including every accepted/failure branch.
    energy_mean=muU+(1-pcode)*(64*native_norm/mlo+8*c)
    number_mean=(energy_mean+b)/(2*(1-eps))
    # An additional two-field-test setting has a larger displacement radius.
    # For r_test <= 1/32, the low-subspace exponential estimate below still
    # gives delta_W <= (K+1)^-2 once K>=threshold (proof states the guard).
    assert radius<F(1,16)
    test_radius=F(1,32)
    assert F(6)*test_radius*bounds.root(K+1)/(K+1-L)<F(1,16)
    Weyl_error=F(1,(K+1)**2)+2*bounds.root((bounds.root(number_mean)+radius)**2/F(L+1))
    field_history_error=full_tree_error+2*Weyl_error
    assert field_history_error<F(1,10**8)
    raw_cap=32*(600+F(360,2**34))+(185+(K+1)*bounds.root(6)+224*g*g+8*c)*mhi+128
    photon_dim=comb(K+5,5)
    mp.iv.dps=90;log=mp.iv.ln
    base=32*(256*log(2)+8*log(5)+24*log(119)+log(2**31+1))+7*log(9)+log(65)+log(photon_dim)
    capacity=rounded(4*base+3*log(M),places=60)
    vals={'relative_Peierls_A_squared_upper':A2,'Young_epsilon':eps,'Young_b_upper':b,
        'vacuum_complement_form_lower':d0,'first_current_vector_form_lower':first,
        'uniform_code_energy_shift_lower_over_m':cstar,'weighted_Schur_term_upper':weighted,
        'strict_Schur_margin_lower':margin,'uniform_internal_code_gap_lower_over_m':gamma,
        'uniform_full_physical_gap_lower':full_gap,'uniform_ground_angle_upper':theta,
        'zero_phase_removed_matter_baseline_lower_over_m':onsite,
        'uniform_raw_ground_absolute_bound':raw_ground_bound,
        'uniform_initial_physical_trial_ground_probability_lower':p,
        'uniform_complete_filter_tree_success_lower':qtree,
        'uniform_conditional_ground_fidelity_lower':fidelity,
        'uniform_conditional_ground_trace_distance_upper':trace,
        'uniform_conditional_ground_purification_error_upper':pur,
        'uniform_actual_ground_T111_amplitude_lower':amplitude,
        'uniform_actual_ground_T111_probability_lower':amplitude*amplitude,
        'uniform_physical_preparation_and_T111_joint_probability_lower':qtree*(amplitude-pur)**2,
        'heat_combined_word_and_Dyson_error_upper':heat_delta,
        'heat_ground_projection_error_upper':heat_ground_error,
        'ideal_ground_two_field_record_error_at_K1024_upper':heat_record,
        'fixed_pointer_complete_prep_and_count_error_upper':full_tree_error,
        'fixed_pointer_complete_prep_and_two_field_tests_error_upper':field_history_error,
        'all_branch_photon_number_mean_after_filters_upper':number_mean,
        'all_cap_raw_H_norm_at_explicit_K_upper':raw_cap}
    exact={k:str(bounds.rounded(v,'upper' in k or k=='uniform_raw_ground_absolute_bound'))
        for k,v in vals.items()}
    return {'schema':'oph.fixed_cube_interacting_ground_refinement.v1','parent_ground_sha256':PARENT_SHA,
        'parent_photon_refinement_sha256':CUTOFF_SHA,
        'same_source_family':{'photon_caps':'every integer K>=3','changed_parameter':'total photon cap only',
            'source_local_completion':'the same +m c(g) sum_x(I-P_x) as the parent, with no scalar subtraction',
            'native_mass_reference_Gauss_density_current_and_clock_unchanged':True,
            'unbounded_five_oscillators':'mathematical comparison, no finite member contains this carrier',
            'spatial_volume_or_spacing_refinement':False},
        'uniform_ground_bounds_exact':exact,'diagnostic_binary64':{k:float(F(v)) for k,v in exact.items()},
        'fixed_complete_energy_clock':{'B':str(B),'pointer_M':M,'filter_count':3,'source_copies':4,
            'ground_to_all_excited_chord_lower':'Delta/(2B), with only |E_ground|<=B and Delta<=B',
            'all_M_basis_records_and_all_noop_failure_branches_retained':True,
            'same_pointer_alphabet_for_all_K_and_mathematical_comparison':True,
            'additional_apparatus_setting_relative_to_parent_K3_clock':True,
            'parent_clock_or_receipt_modified':False,'stages':steps,
            'source_elapsed_time_before_count':'4, including one preparation and three filter slots',
            'source_elapsed_time_after_count':'5, including final common source factor',
            'imaginary_heat_time_is_physical_clock':False},
        'constructive_ideal_ground_history':{'K':heat_K,'photon_dimension':comb(heat_K+5,5),
            'heat_parameter':str(heat_s),'autonomous_dimensionless_wait':str(real_time),
            'two_complete_E_B_test_strengths':['1/64','1/64'],'Dyson_order':J,
            'actual_full_source_matrix_enumerated':False},
        'constructive_complete_preparation_history':{'K':K,'K_dyadic_exponent':96,
            'photon_dimension':photon_dim,'Cayley_norm_error_upper':str(dC),
            'complete_count_alphabet_at_K':290*(K+1)+1,
            'full_source_reference_trial_filter_capacity_nats_interval':capacity,
            'capacity_formula':'4*(32*(256ln2+8ln5+24ln119+ln(2^31+1))+7ln9+ln65+ln binom(K+5,5))+3ln M',
            'additional_complete_count_pointer_capacity':'ln[290(K+1)+1]',
            'full_controller_H_norm_upper':'4 times the displayed all-cap raw source bound',
            'each_duration_one_control_norm_upper':'4 B_raw(K)+pi',
            'gross_preparation_tree_work_upper':'8*(4 B_raw(K)+pi)',
            'matrices_enumerated_or_resource_feasibility_claimed':False},
        'boundaries':{'fixed_cube_interacting_ground_and_complete_record_limit':True,
            'uniform_spatial_or_relativistic_cone':False,'interacting_UV_limit':False,
            'ground_radiation_or_SM_mass_identification':False,
            'all_cap_exact_link_CCR_or_commutativity':False}}


if __name__=='__main__':print(json.dumps(scientific(),indent=2,sort_keys=True))
