"""Exact geometry and rational errors for selected critical-phase work.

The work coefficient is a selected-state expectation, bounded by an actual
trial-energy inequality. It is not evaluated in that trial state. Energy
moments can escape to the cutoff despite this nonzero work limit.
"""
from fractions import Fraction as F

from model.particle_common_contract import contract, provenance
from model.particle_primary_energy_tightness import model as energy
from model.particle_primary_current_tail import model as form
from model.particle_primary_spectral_scale import scale
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS=(
 'model/particle_primary_critical_refinement/__init__.py',
 'model/particle_primary_critical_refinement/__main__.py',
 'model/particle_primary_critical_refinement/model.py',
 'model/particle_primary_critical_refinement/test_model.py',
 'model/particle_primary_critical_refinement/PROOF.md',
 'model/particle_primary_critical_refinement/README.md',
 'model/particle_primary_energy_tightness/PROOF.md',
 'model/particle_primary_local_fluctuation/PROOF.md',
 'model/particle_primary_current_tail/PROOF.md',
 'model/particle_primary_spectral_scale/PROOF.md',
 'model/particle_primary_channel_labels/PROOF.md',
 'model/particle_primary_joint_spectroscopy/PROOF.md',
 'model/particle_common_contract/PROOF.md',
)


def _widths(widths):
    if len(widths)!=3:raise ValueError('three positive physical tent widths are required')
    return tuple(ground._i(x,'physical tent half-width',1) for x in widths)


def tent_moments(radius):
    m=ground._i(radius,'integer tent radius',1)
    return {'square_sum':F(2*m*m+1,3*m),
        'fourth_sum':F(6*m**4+10*m*m-1,15*m**3),
        'gradient_square_sum':F(2,m),'gradient_fourth_sum':F(2,m**3)}


def profile_geometry(level,widths=(1,1,1),amplitude=1):
    """A full finite profile, without enumerating the growing carrier."""
    n=ground._i(level,'source level',0);R=2**n;a=F(1,R)
    b=_widths(widths);s=F(amplitude)
    if not s:raise ValueError('the comparison profile must have nonzero amplitude')
    if max(b)>R:raise ValueError('the profile is wider than half the finite address list')
    moments=[tent_moments(x*R) for x in b];D2=[];D4=[]
    for i in range(3):
        j,k=[x for x in range(3) if x!=i]
        D2.append(moments[i]['gradient_square_sum']*moments[j]['square_sum']*moments[k]['square_sum'])
        D4.append(moments[i]['gradient_fourth_sum']*moments[j]['fourth_sum']*moments[k]['fourth_sum'])
    d2=sum(D2,F());d4=sum(D4,F());shape=sum((F(b[(i+1)%3]*b[(i+2)%3],b[i]) for i in range(3)),F())
    # k_bar/a = pi²*s²*a*D2/6 minus a nonnegative cosine remainder.
    quadratic_pi2=s*s*a*d2/6
    limit_pi2=F(4,27)*s*s*shape
    discretization_pi2=quadratic_pi2-limit_pi2
    cosine_error=F(25,18)*s**4*a**3*d4  # pi^4 < 100
    if discretization_pi2<0:raise ArithmeticError('tent discretization sign failed')
    relative_cosine_error=F(5,6)*s*s*a*a*d4/d2  # pi² < 10
    small_gradient=abs(s)*a*a/min(b)<=1
    work_lower=F(2,3)*s*s*a*d2 if small_gradient else F()
    work_upper=F(5,3)*s*s*a*d2
    return {'level':n,'spacing':a,'physical_half_widths':list(b),'amplitude':s,
        'integer_tent_radii':[x*R for x in b],
        'nonzero_profile_sites':(2*b[0]*R-1)*(2*b[1]*R-1)*(2*b[2]*R-1),
        'phase_root_of_unity_order_divides':2*s.denominator*b[0]*b[1]*b[2]*R**4,
        'axis_gradient_square_sums':D2,'axis_gradient_fourth_sums':D4,
        'three_axis_gradient_square_sum':d2,'three_axis_gradient_fourth_sum':d4,
        'quadratic_work_over_selected_density_pi2_coefficient':quadratic_pi2,
        'limiting_work_over_selected_density_pi2_coefficient':limit_pi2,
        'discretization_pi2_coefficient':discretization_pi2,
        'cosine_nonlinearity_error_upper':cosine_error,
        'absolute_error_from_pi2_limit_upper':10*discretization_pi2+cosine_error,
        'relative_cosine_error_upper':relative_cosine_error,
        'small_gradient_lower_bound_applies':small_gradient,
        'work_over_selected_density_lower':work_lower,
        'work_over_selected_density_upper':work_upper,
        'uniform_work_over_selected_density_lower':F(16,27)*s*s*shape if small_gradient else F(),
        'pointwise_double_commutator_upper':300*s*s*a*d2,
        'pointwise_sign_averaged_work_upper':150*s*s*a*d2}


def profile_ratio(level,widths=(1,1,1),amplitude=1):
    """The unknown selected hopping coefficient cancels exactly."""
    p=profile_geometry(level,widths,amplitude);r=profile_geometry(level)
    eta=p['relative_cosine_error_upper'];eta0=r['relative_cosine_error_upper']
    q=p['quadratic_work_over_selected_density_pi2_coefficient']/r['quadratic_work_over_selected_density_pi2_coefficient']
    if eta0>=1:raise ValueError('reference profile lacks this relative-error enclosure')
    return {'reference_widths':[1,1,1],'reference_amplitude':1,
        'quadratic_discrete_work_ratio':q,
        'actual_original_work_ratio_lower':max(F(),q*(1-eta)),
        'actual_original_work_ratio_upper':q/(1-eta0),
        'limiting_work_ratio':p['limiting_work_over_selected_density_pi2_coefficient']/r['limiting_work_over_selected_density_pi2_coefficient'],
        'same_unknown_selected_density_cancels_exactly':True,
        'ratio_is_energy_weighted_response_not_particle_dispersion':True}


def spatial_frequency_ratio(level,frequency_multiple):
    """Real cosine/sine instruments use q=pi*r and a common normalization."""
    n=ground._i(level,'source level',0);a=F(1,2**n);r=ground._q(frequency_multiple,'positive frequency multiple',True)
    eta0=F(5,6)*a*a;eta=eta0*r*r
    return {'level':n,'frequency_over_pi':r,
        'exact_ratio_formula':'(1-cos(pi*a*r))/(1-cos(pi*a))',
        'ratio_lower':max(F(),r*r*(1-eta)),'ratio_upper':r*r/(1-eta0),
        'limiting_ratio':r*r,'original_common_instrument_normalization':'30V',
        'same_selected_density_cancels':True,
        'fixed_positive_energy_tail_survives':False,
        'phase_setting_is_conserved_momentum':False}


def selected_bound(member):
    contract.verify_member(member);source=member['source'];n=source['level']
    a=F(source['spacing']);R=2**n;L=source['side'];V=L**3
    if a*R!=1 or L!=2*R*R:raise ArithmeticError('primary refinement changed')
    inherited=energy.selected_density_sum_rule(member)
    h=F(inherited['selected_negative_dimensionless_hopping_density_lower'])
    inherited_uniform=F(inherited['uniform_negative_dimensionless_hopping_density_lower'])
    hstar=inherited_uniform+F(1,128)
    finite_correction=F(3,7)*a*a-F(3,28)*a**4-a**8/128
    if h!=hstar+finite_correction or finite_correction<F(281,896)*a*a:
        raise ArithmeticError('actual finite-size hopping-density correction failed')
    hupper=F(30*source['interior_hopping_bonds'],V)
    preparation_source=member['preparation_plan']['source']
    norm=max(abs(F(preparation_source['spectral_lower'])),abs(F(preparation_source['spectral_upper'])))
    mean_allowance=a*a*(F(16,9)*hstar)/128
    clock=spectroscopy.clock_plan(member,mean_allowance/16,mean_allowance/(1024*norm),mean_allowance/(1024*norm))
    clock_error=mean_allowance/8+56*norm*(2*F(clock['clock_plan']['complete_clock_tail_upper'])+F(clock['two_complete_clock_error_upper']))
    if clock_error>F(11,32)*mean_allowance:raise ArithmeticError('complete-clock work error exceeded its budget')
    profiles=[((1,1,1),F(1)),((1,1,1),F(1,2)),((1,1,1),F(2))]
    if n>=1:profiles.extend([((2,1,1),F(1)),((2,2,2),F(1))])
    rows=[]
    for b,s in profiles:
        row=profile_geometry(n,b,s)
        row['original_signed_work_lower']=h*row['work_over_selected_density_lower']
        row['original_signed_work_upper']=hupper*row['work_over_selected_density_upper']
        row['uniform_original_signed_work_lower']=hstar*row['uniform_work_over_selected_density_lower']
        row['liminf_original_signed_work_lower_pi2_coefficient']=hstar*row['limiting_work_over_selected_density_pi2_coefficient']
        row['limsup_original_signed_work_upper_pi2_coefficient']=90*row['limiting_work_over_selected_density_pi2_coefficient']
        row['original_absolute_transfer_moment_upper']=row['original_signed_work_upper']+2*F(member['state']['target_ground_excess'])
        row['finite_reported_signed_work_interval']=[row['original_signed_work_lower']-clock_error,row['original_signed_work_upper']+clock_error]
        if n>=1:row['relative_to_unit_tent']=profile_ratio(n,b,s)
        rows.append(row)
    return contract.serial({'schema':'oph.primary_critical_refinement_result.v1',
        'member_id':member['member_id'],**member['ids'],
        'selected_hopping_producer_sha256':contract.digest(contract.serial(inherited)),
        'selected_density_symbol':'h_n=-a*Tr(rho_n H_hop)/V',
        'selected_density_lower':h,'uniform_selected_density_lower':hstar,
        'inherited_uniform_selected_density_lower':inherited_uniform,
        'exact_finite_size_density_floor_correction':finite_correction,
        'finite_size_density_floor_correction_lower':F(281,896)*a*a,
        'selected_density_upper':hupper,'uniform_selected_density_upper':90,
        'profiles':rows,'spatial_frequency_ratios':[spatial_frequency_ratio(n,r) for r in (F(1,2),F(2),F(3))],
        'complete_retained_origin_and_sign_outcomes':2*V,
        'complete_two_energy_clock':clock,
        'complete_energy_clock_signed_mean_error_upper':clock_error,
        'complete_fine_record_count':2*V*clock['clock_plan']['clock_size']**2,
        'complete_clock_mean_error_is_O_a_squared':True,
        'complete_clock_histogram_numerically_executed':False,
        'same_source_state_idle_reference_and_clock':True,
        'intervention_preserves_source_state':False,
        'initial_energy_record_preserves_joint_purification':False,
        'original_unitary_spectral_mass':1,
        'work_ratio_is_independent_of_unknown_selected_density':True,
        'subsequence_positive_Dirichlet_work_form':True,
        'work_form_transferred_to_limiting_probability_measure':False,
        'inelastic_fixed_energy_probability_proved':False,
        'actual_uniform_dimensionless_gap_proved':False,
        'pure_relative_bond_remainder_squared_norm_lower':15,
        'pure_relative_full_star_remainder_squared_norm_lower':15*(3 if n==0 else 6),
        'relative_form_boundary_applies_to_fixed_product_anchor':True,
        'all_gap_stability_theorems_excluded':False,
        'physical_mass_table_GeV':[]})


def bind_member(member):
    result=selected_bound(member)
    return {'schema':'oph.primary_critical_refinement_bound.v1','result':result,
        'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('critical-profile source or bound changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*energy.SOURCE_INPUTS,*form.SOURCE_INPUTS,*scale.SOURCE_INPUTS,*spectroscopy.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_critical_refinement_certificate.v1',
        'members':[bind_member(contract.member(n)) for n in (0,1,2,8)],
        'source_sha256':pins,'physical_mass_table_GeV':[]}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('critical-profile inputs changed during replay')
    return result
