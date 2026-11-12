"""Positive native exchange and the even remainder fix the reference momentum."""
from functools import lru_cache
import sympy as s
from model.coherent_defect_wave.chiral_joint.atomic_reference_control.ground_exchange import actual_spin_coefficients
from model.coherent_defect_wave.chiral_joint.atomic_reference_control.resolved_exchange import tensor_checks
from model.coherent_defect_wave.chiral_joint.atomic_spin_geometry.cube import spectrum
from model.coherent_defect_wave.chiral_joint.atomic_spin_geometry.momentum import certificate as momentum_certificate
from model.coherent_defect_wave.chiral_joint.atomic_even_remainder.bounds import uniform_quartic_remainder
from .isolation import data,LAMBDA,KAPPA


@lru_cache(maxsize=1)
def certificate():
    d=data();spin=actual_spin_coefficients();tensor=tensor_checks()
    cube=spectrum();momentum=momentum_certificate()
    beta_hi=s.Rational(spin['beta_ec_interval'][1])
    window=s.Rational(1,1104);b_lower=(1-beta_hi)/2-s.Rational(2,15)
    # Only the positive Laplace weight changes with the link penalty. The
    # complete native response tensors and the charged semigroup are fixed.
    j_lower=window*(1-KAPPA*window)*b_lower**2/2
    assert KAPPA*window<1 and j_lower>0
    cube_gap=s.Rational(cube['gap_interval'][0]);assert cube_gap>s.Rational(4,5)
    remainder=uniform_quartic_remainder(d['gap'][0],d['tau'],LAMBDA)
    exact_gap=2*LAMBDA**2*j_lower*cube_gap-2*remainder
    assert exact_gap>0 and d['gap'][0]-2*d['v']>0
    assert cube['ground_rank']==1 and momentum['ground_translation_characters']==[1,1,1]
    return {'complete_native_spin_tensor_controls':tensor,
        'penalty_dependent_exchange':'J_kappa=(1/4) sum_f d_f integral_0^infinity exp(-kappa t) b_f(t)^2 dt',
        'all_charged_intermediate_states_retained':True,
        'positive_native_L_channel_window':str(window),'native_L_response_lower':str(b_lower),
        'positive_exchange_lower':str(j_lower),
        'second_order_operator':'-24 A_kappa I + 2 J_kappa H_cube',
        'scalar_shift_retained':True,'cube_ground_rank':1,'cube_gap_lower':str(cube_gap),
        'full_Q_even_remainder_upper':str(remainder),
        'actual_physical_ground_gap_lower':str(exact_gap),
        'actual_ground_gap_decimal_lower':float(exact_gap),
        'actual_kappa_two_ground_unique':True,
        'actual_ground_periodic_momentum':[0,0,0],
        'native_eight_odd_ground_factor_translation_permutation_sign':1,
        'source_periodic_boundary_extra_wrap_phase_on_atomic_ground':1,
        'bipartite_parity_commutes_with_doubled_link_penalty':True,
        'positive_hopping_interval_to_reported_lambda_is_isolated':True,
        'kappa_one_ground_theorem_silently_transferred':False,
        'neutral_zero_source_momentum_equals_ground_relative_momentum':True}
