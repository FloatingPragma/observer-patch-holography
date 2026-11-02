"""Rational majorants for all nine local state and modular responses."""
from fractions import Fraction as F
from math import factorial

GAMMA_PI3 = F(4249388,150088295)
INVERSE_NORM = F(68,5)
H = tuple(16*factorial(m)*64**m for m in range(4))
HESSIAN = 243*(3*H[1]*H[2]+H[0]*H[3])
# pi<4: gamma_L > GAMMA_PI3/(64 L^4).
BALL_ENTRY_ERROR = F(3,14)*HESSIAN*64/GAMMA_PI3
ENTROPY_ENTRY_ERROR = 5/GAMMA_PI3
MAX_RADIUS_RATIO = F(1,2**28)


def rational(value,name):
    if type(value) not in (int,F):
        raise ValueError('exact rational '+name+' required')
    return F(value)


def neumann_budget(radius_ratio,*,finite_entropy_difference=False):
    q = rational(radius_ratio,'radius ratio')
    if not 0<q<=F(1,8) or type(finite_entropy_difference) is not bool:
        raise ValueError('0<R/L<=1/8 and a Boolean estimator selector required')
    entry = BALL_ENTRY_ERROR+(ENTROPY_ENTRY_ERROR if finite_entropy_difference else 0)
    return 9*INVERSE_NORM*entry*q*q


def entropy_remainder(l1_parameter,radius=F(1,2)):
    t=rational(l1_parameter,'parameter norm');r=rational(radius,'parameter bound')
    if not 0<=t<=r<1:raise ValueError('0<=parameter norm<=bound<1 required')
    return t*t/(2*(1-r*r))


def certificate():
    if not (H==(16,1024,131072,25165824)):
        raise ArithmeticError("exact response certificate failed")
    if not (HESSIAN==195689447424):
        raise ArithmeticError("exact response certificate failed")
    if not (neumann_budget(MAX_RADIUS_RATIO,finite_entropy_difference=True)<F(1,4)):
        raise ArithmeticError("exact response certificate failed")
    return {
        'spacetime_derivative_majorants':list(H),
        'derivative_bound':'||partial^alpha F_normalized|| <= H_m L^(-3/2-m), |alpha|=m<=3',
        'domain':'|x0|+|x|<=L/4; shell L<=|y|<=2L',
        'complex_Cauchy_radius_over_L':'1/64',
        'boosted_spatial_Hessian_majorant':HESSIAN,
        'boosted_Hessian_scaling':'HESSIAN*||theta||_1/L^6',
        'normalized_ball_entry_error_coefficient':str(BALL_ENTRY_ERROR),
        'normalized_entropy_difference_extra_coefficient':str(ENTROPY_ENTRY_ERROR),
        'certified_radius_ratio':str(MAX_RADIUS_RATIO),
        'modular_response_Neumann_budget':str(neumann_budget(MAX_RADIUS_RATIO)),
        'finite_entropy_response_Neumann_budget':str(neumann_budget(MAX_RADIUS_RATIO,finite_entropy_difference=True)),
        'both_budgets_strictly_below':'1/4',
        'normalized_response_inverse_infinity_norm_upper':'272/15',
        'finite_entropy_difference_step':'epsilon=(R/L)^6',
        'relative_entropy_upper':'||theta||_1^2/[2(1-r0^2)] for ||theta||_1<=r0<1',
        'global_pure_vacuum_affine_two_sided_perturbation_asserted':False,
        'joint_lattice_collar_limit_proved':False,
    }
