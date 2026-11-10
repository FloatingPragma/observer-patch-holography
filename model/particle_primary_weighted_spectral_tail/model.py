"""Close an actual selected-state second moment through local third moments."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

from model.particle_common_contract import contract, provenance
from model.particle_primary_work_spectrum import model as work
from model.particle_primary_locality import locality
from model.particle_primary_stress import stress

SOURCE_INPUTS=tuple('model/particle_primary_weighted_spectral_tail/'+name for name in
    ('__init__.py','__main__.py','model.py','test_model.py','PROOF.md','README.md'))+(
    'model/particle_primary_work_spectrum/PROOF.md',
    'model/particle_primary_critical_refinement/PROOF.md',
    'model/particle_primary_locality/PROOF.md',
    'model/particle_primary_stress/PROOF.md',
    'model/particle_gauge_characteristics/PROOF.md',
    'model/particle_primary_current_tail/PROOF.md')


def sqrt_upper(value,bits=32):
    """A rational outward square root, with error at most2^-bits."""
    value=F(value)
    if value<0 or type(bits) is not int or bits<0:raise ValueError('nonnegative radicand and bit count required')
    scale=1<<bits;scaled=value*scale*scale
    integer=(scaled.numerator+scaled.denominator-1)//scaled.denominator
    root=isqrt(integer)
    if root*root<integer:root+=1
    result=F(root,scale)
    if result*result<value:raise ArithmeticError('square root rounded inward')
    return result


def quadratic_moment_cap(positive_first_upper,third_constant,third_M2_coefficient,
                         negative_second_upper,negative_third_upper):
    """(M2-N2)_+² <= P1*(A+B*M2+N3), with no success normalization."""
    P,A,B,N2,N3=map(F,(positive_first_upper,third_constant,third_M2_coefficient,
                       negative_second_upper,negative_third_upper))
    if min(P,A,B,N2,N3)<0 or P==0:raise ValueError('positive first moment and nonnegative caps required')
    center=N2+B*P/2
    radicand=(B*P/2)**2+P*(A+N3+B*N2)
    result=(center+sqrt_upper(radicand))*(1<<32)
    return F((result.numerator+result.denominator-1)//result.denominator,1<<32)


@lru_cache(maxsize=1)
def native_derivative_control():
    data=contract.characteristics.native_generators()
    diagonal=tuple(F(str(x)) for x in data['casimir'].diagonal())
    gram=tuple(F(str(data['gram'][i,i])) for i in range(12))
    if len(data['generators'])!=12 or sum(diagonal)!=12 or min(gram)<=0:
        raise ArithmeticError('actual electric trace-metric generator normalization changed')
    return {'native_Casimir_diagonal':diagonal,'native_generator_Gram_diagonal':gram,
        'sum_normalized_generator_squared_norm_upper':12,
        'sum_normalized_generator_squared_norm_bound_source':'sum_A ||T_A||² <=sum_A Tr(T_A²)=Tr(Casimir)=12',
        'word_projection_commutes_with_each_electric_generator':True}


def local_current_caps(member):
    contract.verify_member(member);source=member['source'];n=source['level']
    a=F(source['spacing']);ell=source['per_link_word_degree'];V=source['side']**3
    bonds=source['interior_hopping_bonds'];D2=work.critical.profile_geometry(n)['three_axis_gradient_square_sum']
    # Full product-cutoff owner supports: matter1, electric1, hopping2,
    # magnetic3. A site meets at most6 hopping and9 magnetic terms.
    kappa=90+3*F(source['electric_coefficient'])*F(29,10)*ell**2+6*30+9*F(source['magnetic_coefficient'])
    if kappa>=1200*ell**2:raise ArithmeticError('source local norm envelope changed')
    c_norm=F(180) if n==0 else 150*a*D2
    if c_norm>F(2025,4):raise ArithmeticError('even phase change lacks its uniform cap')
    # Individual cosine bond norm <=150a³; each cell meets <=6 such bonds.
    kappa_C=900*a**3
    squares=3000*F(bonds,V)*D2
    double_H=240*kappa/a*squares
    double_C=240*kappa_C*squares
    norm_A=double_H/2+double_C;B=3*c_norm
    derivative=native_derivative_control()
    # The original state's mean electric energy is <=91V/a. Combining
    # its Casimir expectation with the exact multiplication derivatives
    # removes the link-degree factor from this selected expectation.
    kappa_no_electric=270+9*F(source['magnetic_coefficient'])
    electric_density=(30*F(bonds,V)+F(source['magnetic_coefficient'])*F(source['interior_plaquettes'],V))/a+F(member['state']['target_ground_excess'])/V
    if electric_density>=91/a:raise ArithmeticError('actual selected electric density exceeds its uniform cap')
    selected_electric=F(source['electric_coefficient'])/(2*a)*(24*squares+F(6912000)/a+2880000*a**4*D2)
    selected_A=120*kappa_no_electric/a*squares+selected_electric/2+double_C
    A=min(norm_A,selected_A)
    if norm_A>19440000000*ell**2/a**2 or selected_A>14100000000/a**2:
        raise ArithmeticError('family third-moment envelope changed')
    return {'spacing':a,'link_degree':ell,'owner_H_norm_sum_times_a_upper':kappa,
        'cosine_operator_norm_upper':c_norm,
        'cosine_owner_norm_sum_upper':kappa_C,
        'origin_averaged_cosine_plus_sine_bond_norm_squares_upper':squares,
        'origin_averaged_H_double_commutator_sum_upper':double_H,
        'origin_averaged_sine_sine_cosine_double_commutator_upper':double_C,
        'native_electric_derivative_control':derivative,
        'selected_electric_energy_upper_per_cell_from_actual_box':electric_density,
        'selected_electric_energy_upper_per_cell':91/a,
        'selected_total_Casimir_expectation_upper_per_cell':F(182,105),
        'selected_electric_double_commutator_expectation_sum_upper':selected_electric,
        'cutoff_norm_third_constant_upper':norm_A,
        'selected_electric_third_constant_upper':selected_A,
        'signed_third_constant_upper':A,'signed_third_M2_coefficient_upper':B,
        'third_moment_operator_identity':'m3=half average Tr rho([C,[H,C]]+[S,[H,S]])+average Tr rho(C^3+3SCS+[S,[S,C]])',
        'product_link_cutoff_and_all_physical_sectors_retained':True,
        'parent_or_Haar_state_substituted':False}


def bind_member(member):
    contract.verify_member(member);inherited=work.bind_member(member);old=inherited['result']
    caps=local_current_caps(member);exact_low=F(old['actual_modular_slope']['modular_slope_lower'])
    low=F(exact_low.numerator//exact_low.denominator)
    if low<=0 or low>exact_low:raise ArithmeticError('positive compact modular slope required')
    n1=1/(3*low);n2=2/low**2;n3=6/low**3
    positive=F(old['original_signed_work_interval'][1])+n1
    if positive>=444 or max(n2,n3)>1:raise ArithmeticError('source positive/negative moment envelope changed')
    weighted=quadratic_moment_cap(positive,caps['signed_third_constant_upper'],
                                  caps['signed_third_M2_coefficient_upper'],n2,n3)
    old_cap=F(old['second_moment']['selected_second_moment_upper']);selected=min(old_cap,weighted)
    a=caps['spacing'];ell=caps['link_degree'];uniform=4000000/a
    if selected>uniform:raise ArithmeticError('selected O(1/a) bound failed')
    mean=F(old['original_signed_work_interval'][0]);ideal=work.local.positive_window(mean,selected)
    # The unchanged critical clock was already authenticated for the wider
    # inherited window. Its error is applied to the original event directly.
    clock=old['finite_reported_positive_window'];eta=F(clock['energy_resolution'])
    loss=F(clock['two_complete_clock_probability_loss_upper'])
    left,right=ideal['original_transfer_interval'];probability=ideal['original_weight_lower']
    if left<=2*eta or loss>=probability/2:raise ArithmeticError('same complete clock lost the positive event')
    result=contract.serial({'schema':'oph.primary_weighted_spectral_tail_result.v1',
        'member_id':member['member_id'],**member['ids'],
        'source_work_bound_sha256':inherited['result_sha256'],
        'actual_menu':'all critical unit-tent origins and both signs with weight1/(2V)',
        'positive_integer_modular_slope_lower_used':low,
        'local_current_caps':caps,'negative_moment_uppers':{'first':n1,'second':n2,'third':n3},
        'positive_first_moment_upper':positive,
        'selected_weighted_quadratic_second_moment_upper':weighted,
        'inherited_operator_second_moment_upper':old_cap,
        'selected_second_moment_upper':selected,
        'second_moment_cap_improvement_factor':old_cap/selected,
        'family_second_moment_upper':uniform,
        'family_second_moment_formula':'M2 <=4000000/a; a=2^-n, independent of ell=n+1',
        'weighted_tail_bound_formula':'integral_{|omega|>R}|omega|dmu <=selected_second_moment_upper/R',
        'ideal_original_positive_window':ideal,
        'finite_reported_positive_window':{'reported_transfer_interval':[left-2*eta,right+2*eta],
            'original_probability_lower':probability-loss,
            'two_complete_clock_probability_loss_upper':loss,
            'same_critical_two_clock_plan_sha256':clock['same_critical_two_clock_plan_sha256'],
            'energy_resolution':eta,'all_fine_records_retained':clock['all_fine_records_retained'],
            'conditioned_or_postselected_probability':False},
        'source_state_reference_clock_and_action_unchanged':True,
        'finite_window_floor_order':'a',
        'window_upper_endpoint_order':'1/a',
        'uniform_selected_second_moment_proved':False,
        'uniform_selected_energy_weighted_tail_proved':False,
        'fixed_positive_energy_probability_limit_proved':False,
        'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_weighted_spectral_tail_bound.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('weighted spectral moment or selected source changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*work.SOURCE_INPUTS,*locality.SOURCE_INPUTS,*stress.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_weighted_spectral_tail_certificate.v1',
        'members':[bind_member(contract.member(n)) for n in (0,1,2,8)],
        'source_sha256':pins,'physical_mass_table_GeV':[]}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('weighted spectral inputs changed during replay')
    return result
