"""On-shell second-order exchange, with every charged denominator retained."""
from functools import lru_cache
import sympy as s
from ..higgs_composites import model as car
from ..neutral_local_spectrum.symmetries import conjugate_polynomial
from ..periodic_family.model import Configuration,kinetic_terms,cells
from ...chiral_source.model import PAULI
from .ground_exchange import second_moment_certificate
from .lattice import gap_data,hopping_norm,local_norm
from .perturbation import parameters,iteration_bounds


@lru_cache(maxsize=1)
def tensor_checks():
    a,b=s.symbols('a b',real=True);theta=s.I*PAULI[1]
    def response(spin,other,sign):
        # sign + is Fminus_(spin,other); sign - is Fplus_(other,spin).
        return a*int(spin==other)*s.eye(2)+sign*b*sum((sigma[other,spin]*sigma/2 for sigma in PAULI),s.zeros(2))
    for spin in range(2):
        for other in range(2):
            minus=response(spin,other,1)
            assert minus.H==response(other,spin,1)
            assert (theta*minus.conjugate()*theta.H-response(other,spin,-1)).applyfunc(s.expand)==s.zeros(2)
    for mode in range(30):
        assert conjugate_polynomial(car.monomial(((mode,False),)))==car.scale(car.monomial(((mode,True),)),s.I)
    for axis in range(3):
        spin=car.spin_current(axis)
        assert conjugate_polynomial(spin)==car.scale(spin,-1)
    rows=[]
    for axis,sigma in enumerate(PAULI):
        result=s.zeros(4)
        for (spin,t),left in sigma.todok().items():
            for (u,w),right in sigma.todok().items():
                result+=s.conjugate(left)*right*s.kronecker_product(response(u,spin,-1),response(t,w,1))
                result+=left*s.conjugate(right)*s.kronecker_product(response(spin,u,1),response(w,t,-1))
        exchange=sum(((1 if i==axis else -1)*s.kronecker_product(p/2,p/2) for i,p in enumerate(PAULI)),s.zeros(4))
        expected=4*(a*a*s.eye(4)-b*b*exchange)
        assert (result-expected).applyfunc(s.expand)==s.zeros(4)
        x=s.Symbol('x')
        assert s.factor(exchange.charpoly(x).as_expr())==(4*x-1)**3*(4*x+3)/256
        rows.append({'axis':axis,'complete_two_orientation_spin_contraction':str(expected),
                     'scalar_coefficient':'4 a(t)^2','spin_coefficient':'-4 b(t)^2',
                     'exchange_eigenvalues':['-3/4','1/4','1/4','1/4']})
    terms=kinetic_terms(Configuration());assert len(terms)==24
    assert {term.path for term in terms}=={(e,) for e in range(24)}
    assert all(abs(term.coefficient)**2==s.Rational(1,16) for term in terms)
    return {'actual_CAR_charge_conjugation_mode_identities':30,
            'actual_ground_spin_charge_conjugation_identities':3,
            'ground_antiunitary_matrix':'i sigma_y times complex conjugation',
            'pointwise_response_covariance':'Fminus_st=P c_s^dagger E(t)c_t P=a_f delta_st I+b_f sum_i sigma_i,ts S_i; Fplus_st=P c_s E(t)c_t^dagger P=a_f delta_st I-b_f sum_i sigma_i,st S_i',
            'a_and_b_real':'diagonal response blocks are Hermitian; their two spin expectations fix a and b',
            'Hermitian_and_antiunitary_response_tensor_checks':True,'exact_spin_contractions':rows,
            'full_PW_blank_Gram_factor':'delta_ac delta_bd/d_f, including both orthogonal flux orientations',
            'internal_index_sum_with_Gram_factor':'d_f^2/d_f=d_f',
            'all_24_actual_paths_are_unique_single_edges':True,
            'cross_edge_denominator_terms':'0, distinct excited-link sectors are orthogonal and H0 preserves them'}


@lru_cache(maxsize=1)
def certificate():
    tensors=tensor_checks();moment=second_moment_certificate();atomic=gap_data()
    spin=moment['actual_local_ground_spin_coefficients']
    beta_lo,beta_hi=map(s.Rational,spin['beta_ec_interval'])
    chi_lo,chi_hi=map(s.Rational,spin['chi_interval'])
    centered_local=local_norm()+12;assert centered_local==69
    window=1/(16*centered_local);assert window==s.Rational(1,1104)
    response_error=s.Rational(2,15)
    bL_lower=(1-beta_hi)/2-response_error;assert bL_lower>0
    exchange_lower=s.cancel(window*(1-window)*bL_lower**2/2)
    moment_norm=(15+3*chi_hi)/16
    exchange_upper=s.cancel(moment_norm/atomic['fluxgap'])
    assert 0<exchange_lower<exchange_upper
    tau=hopping_norm();gap=atomic['gap']
    splitting_lambda=min(gap/(5*tau),s.cancel(gap**2*exchange_lower/(24*tau**3)))
    assert splitting_lambda>0
    p=parameters(coupling=splitting_lambda)
    v=p['v'];gamma=p['separated']
    exact_remainder=s.cancel(v**3/(gap*gamma)+12*v**4/gamma**3)
    simple_remainder=s.cancel(6*v**3/gap**2)
    assert exact_remainder<=simple_remainder
    split=s.cancel(splitting_lambda**2*exchange_lower-2*simple_remainder)
    assert split>=splitting_lambda**2*exchange_lower/2>0
    # The rational error coefficient increases for 0<=v/g<=1/5.
    u=s.Symbol('u',nonnegative=True)
    factor=1/(1-u)+12*u/(1-u)**3
    assert factor.subs(u,s.Rational(1,5))==s.Rational(95,16)<6
    assert s.factor(s.diff(factor,u))==(u*u+22*u+13)/(u-1)**4
    # On the two-cell-per-axis cube, flipping one spin changes six incident
    # bond expectations. The exact total is one exchange coefficient.
    cells2=cells(2);difference=s.Integer(0)
    for term in kinetic_terms(Configuration()):
        if 0 not in (term.left_cell,term.right_cell):continue
        sign=1 if term.axis==2 else -1
        difference-=s.Rational(sign,2)
    assert difference==1
    rotations=[s.diag(*(1 if i==axis else -1 for i in range(3))) for axis in range(3)]
    onsite=[]
    for xyz in cells2:
        onsite.append(rotations[0]**xyz[0]*rotations[1]**xyz[1]*rotations[2]**xyz[2])
    for term in kinetic_terms(Configuration()):
        assert onsite[term.left_cell].T*rotations[term.axis]*onsite[term.right_cell]==s.eye(3)
    return {'tensor_and_source_checks':tensors,
            'exact_second_order_operator':'K2=-P0 T Q [Q(H0-volume*E0)Q]^(-1) Q T P0',
            'Q_domain':'complete physical lattice Gauss kernel orthogonal to P0',
            'full_charged_denominator_positive_lower_bound':str(atomic['fluxgap']),
            'denominator_is_evaluated_at_actual_atomic_ground_energy':True,
            'complete_edge_decomposition':'K2=sum_e [-A_edge I+J_edge sum_i R_axis,ii S_x,i S_y,i]',
            'scalar_part_retained':'A_edge=4|c_edge|^2 sum_f d_f integral_0^infinity exp(-t) a_f(t)^2 dt',
            'exchange_coefficient':'J_edge=4|c_edge|^2 sum_f d_f integral_0^infinity exp(-t) b_f(t)^2 dt',
            'all_native_families_in_both_integrals':True,
            'actual_complete_onsite_semigroup':'exp[-t(H_local-E0)], without charged-sector truncation',
            'uniform_edge_coefficient_absolute_squared':'1/16',
            'certified_J_edge_interval':[str(exchange_lower),str(exchange_upper)],
            'scalar_A_edge_interval':['0',str(exchange_upper)],
            'strict_lower_bound_channel':'actual L doublet, b_L(0)=(1-beta_ec)/2',
            'local_centered_H_norm_upper_bound':str(centered_local),
            'positive_integration_window':[str(0),str(window)],
            'uniform_bL_deviation_bound_on_window':str(response_error),
            'uniform_bL_lower_bound_on_window':str(bL_lower),
            'nonnegative_other_channels_and_later_times_retained':True,
            'all_24_oriented_edges_have_same_positive_J':True,
            'site_SO3_transformation':'O_x=R_xaxis^x R_yaxis^y R_zaxis^z',
            'exact_transformed_second_order_spin_operator':'J sum_24_oriented_edges S_x dot S_y = 2J sum_12_cube_edges S_x dot S_y',
            'constant_scalar_shift':'-24 A_edge',
            'canonical_graph_effective_H_remainder':'||K_lambda-volume*E0-lambda^2 K2|| <= 6 (|lambda| tau)^3/gap^2',
            'full_state_remainder_bound_uses_all_Q':True,
            'admitted_positive_lambda_for_nonconstant_exact_band':str(splitting_lambda),
            'exact_band_energy_expectation_contrast_lower_bound':str(split),
            'contrast_witness':'canonical dressed all-up atomic spins versus the same state with origin spin flipped',
            'canonical_graph_iteration_control_at_splitting_lambda':iteration_bounds(2,coupling=splitting_lambda),
            'unit_hopping_claim_transferred':False,
            'exact_low_band_ground_spin_spectrum_identified':False}
