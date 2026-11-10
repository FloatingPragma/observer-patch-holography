"""Even neutral scalar controls, full-carrier failures and uniform evaluation."""
import sympy as s
from .model import source_identities


def square_root_series(order):
    if type(order) is not int or order<0:raise ValueError('nonnegative exact polynomial order required')
    ratio=s.Rational(165,256);coefficients=[s.Integer(1)];a=s.Rational(1,2)
    for k in range(1,order+1):
        coefficients.append(-a)
        a*=s.Rational(2*k-1,2*k+2)
    tail=a*ratio**(order+1)/(1-ratio)
    return coefficients,tail


def certificate():
    identities=source_identities()
    assert identities['strict_rational_Y_squared_upper_bound']=='165'
    assert s.Rational(165,256)<1
    estimates=[]
    for order in (0,16,64,128):
        coefficients,tail=square_root_series(order)
        estimates.append({'order':order,'failure_operator_error_upper_bound':str(tail),
            'finite_polynomial_unitarity_defect_upper_bound':str(2*tail),
            'finite_polynomial_claimed_exact_unitary':False})
    return {'source_preparation':'selected normalized global physical ground projector P_ground/2',
        'ground_selection_is_declared_preparation_law':True,
        'source_total_fermion_parity':'odd','both_probe_operators_are_even':True,
        'both_probe_operators_have_physical_spin':'0','both_probe_operators_preserve_raw_L_and_B3':True,
        'both_probe_operators_commute_with_all_twelve_Gauss_currents':True,
        'settings_probability':'1/2 each, recorded',
        'scalar_number_success':'N/2','scalar_number_full_failure':'I-N+(sqrt(3)/2)N',
        'total_Yukawa_success':'Y/16','total_Yukawa_full_failure':'positive sqrt(I-Y^2/256)',
        'full_CAR_Kraus_identity':'K^dagger K+F^dagger F=I',
        'even_pointer_unitary':'U=[[F,-K],[K,F]], K=K^dagger, F=sqrt(I-K^2)',
        'both_unitary_identities_follow_from_commuting_positive_square_root':True,
        'setting_and_outcome_registers_retained':True,'failure_has_no_energy_postselection':True,
        'Y_failure_exact_ground_action':'sqrt(1-xi/256) P_ground',
        'Y_failure_is_defined_on_full_carrier':True,
        'uniform_square_root_algorithm':'F_n=I-sum_(k=1)^n a_k (Y^2/256)^k; a1=1/2, a_(k+1)=a_k(2k-1)/(2k+2)',
        'uniform_square_root_tail':'a_(n+1)(165/256)^(n+1)/(1-165/256)',
        'uniform_full_carrier_approximation_bounds':estimates,
        'exact_instrument_uses_unique_positive_finite_matrix_square_root':True,
        'controller_is_declared_finite_apparatus':True,'controller_is_idle_Hamiltonian':False,
        'arbitrary_input_ground_filter_probability_supplied':False}
