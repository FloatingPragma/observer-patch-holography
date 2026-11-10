"""Full-carrier resolvent and unitary Taylor bounds with complete outcomes."""
import sympy as s
from ..periodic_family.checks import norm_bound
from .model import CONFIG,BLOCKS,leading_histories


def feshbach_bound(setting):
    d=len(BLOCKS[setting]);b=norm_bound(CONFIG)
    moment=s.Rational(leading_histories(setting)['full_H_moments_0_through_5'][4])
    radius=1024*d*b**5;assert radius>b
    leading=-moment/radius**3
    error=b**5/(radius**4*(1-b/radius))
    assert 0<error<leading
    return {'setting':setting,'P':'fixed gauge/spin-singlet matter-reference preparation, all scalar/link states',
            'full_unconstrained_P_rank':'119^24 * 5^8',
            'full_unconstrained_Q_rank':'(2^270-1) * 119^24 * 5^8',
            'Gauss_restriction':'P and Q commute with every current; apply the same identity on each complete Gauss sector',
            'P_plus_Q_equals_full_carrier_identity':True,
            'P_H_P':'scalar and link number energies; no off-diagonal loop term',
            'elimination':'H_eff(z)=PHP-PHQ(QHQ-z)^(-1)QHP',
            'z':str(-radius),'certified_Q_resolvent_distance_lower_bound':str(radius-b),
            'Q_irreducible_two_and_three_factor_loop_moments':['0','0'],
            'Q_irreducible_four_factor_loop_moment':str(moment),
            'leading_loop_coefficient_at_z':str(leading),
            'all_higher_full_H_orders_error_upper_bound':str(error),
            'exact_loop_kernel_real_part_enclosure':[str(leading-error),str(leading+error)],
            'exact_loop_kernel_imaginary_part_absolute_upper_bound':str(error),
            'exact_loop_kernel_nonzero_by_strictly_positive_real_part':True,
            'energy_dependence':'leading moment/z^3 plus the complete resolvent remainder',
            'autonomous_energy_independent_pure_gauge_H_derived':False}


def complete_records():
    b=norm_bound(CONFIG);time=1/(64*max(map(len,BLOCKS.values()))*b**3)
    remainder=b**6*time**6/s.factorial(6);rows=[]
    for setting,block in BLOCKS.items():
        moments=leading_histories(setting)['full_H_moments_0_through_5']
        m4,m5=map(s.Rational,moments[4:6]);assert m4<0 and m5==10*m4
        leading=-m4*time**4/s.factorial(4)
        fifth=abs(m5)*time**5/s.factorial(5)
        lower=leading-remainder;upper=leading+fifth+remainder
        assert 0<lower<upper<1
        success=(lower**2/2,upper**2/2)
        failure=(s.Rational(1,2)-success[1],s.Rational(1,2)-success[0])
        assert failure[0]>0
        rows.append({'setting':setting,'choice_probability':'1/2',
            'time_in_same_source_clock':str(time),'leading_amplitude_absolute_value':str(leading),
            'fifth_order_imaginary_amplitude_absolute_value':str(fifth),
            'full_unitary_sixth_order_remainder_upper_bound':str(remainder),
            'success_effect':'rank-one projector onto the declared normalized source-reference singlet times loop trace state',
            'failure_effect':'identity minus success on the complete source and reference carrier',
            'original_setting_success_probability_interval':[str(v) for v in success],
            'original_setting_failure_probability_interval':[str(v) for v in failure],
            'postmeasurement_states':'retain E exp(-itH) initial/sqrt(2) for both complementary outcomes',
            'success_conditioning_applied':False})
    lower=sum(s.Rational(row[key][0]) for row in rows for key in
              ('original_setting_success_probability_interval','original_setting_failure_probability_interval'))
    upper=sum(s.Rational(row[key][1]) for row in rows for key in
              ('original_setting_success_probability_interval','original_setting_failure_probability_interval'))
    assert lower<=1<=upper
    return {'same_complete_PW_Hamiltonian_used':True,'same_clock_for_both_settings':True,
            'Taylor_method':'unitary integral remainder after degree five; norm <= ||H||^6 t^6/720',
            'real_leading_and_imaginary_fifth_term_separated_exactly':True,
            'settings':rows,'complete_four_outcome_probability_enclosure':[str(lower),str(upper)],
            'all_complementary_states_and_reference_occupations_retained':True,
            'resolvent_kernel_used_as_time_generator':False}
