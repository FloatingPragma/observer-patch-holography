"""Exact domain, entropy and symbolic refinement bounds; no mesh underflow."""
from fractions import Fraction as F
from model.shrinking_modular_entropy.bounds import modulus as linear_modulus
from model.shell_entropy_response.bounds import modulus as hessian_modulus

C_GAMMA = F(4249388, 150088295)
C_B = F(3,14)*195689447424*64/C_GAMMA


def rational(value, name):
    if type(value) not in (int,F):
        raise ValueError(name+' requires an exact rational, not bool/float')
    return F(value)


def domain(inner=F(1,2), outer=4):
    d=rational(inner,'inner');D=rational(outer,'outer')
    if not 0<d<D:raise ValueError('0<inner<outer required')
    return {'inner':d,'outer':D,'maximum_rho':d/4,'exponential_norm_coefficient':8192*D**3/d**6}


def estimates(rho, amplitude, radius=F(1,2), inner=F(1,2), outer=4):
    r=rational(rho,'rho');t=rational(amplitude,'amplitude');b=rational(radius,'radius')
    g=domain(inner,outer)
    if not 0<r<=g['maximum_rho'] or not abs(t)<=b<1:
        raise ValueError('0<rho<=inner/4 and |amplitude|<=radius<1 required')
    K=g['exponential_norm_coefficient'];s=1-b
    return {'norm_squared_upper':K*r**3,
            'fractional_moment_upper':2*K*r**3,
            'hessian_zero_upper':26*K*K*r**6,
            'relative_entropy_upper':13*K*K*t*t*r**6/(s*s),
            'quadratic_remainder_upper':26*K*K*abs(t)**3*r**6/(3*s**3)}


def response_coefficient(radius=F(1,2)):
    b=rational(radius,'radius')
    if not 0<=b<1:raise ValueError('0<=radius<1 required')
    K=domain()['exponential_norm_coefficient']
    inverse_n0=F(165,28)/C_GAMMA
    return C_B+inverse_n0*(1+b/2+13*K*K*b/(1-b)**2)


def neumann_budget(rho, radius=F(1,2)):
    r=rational(rho,'rho')
    if not 0<r<=F(1,8):raise ValueError('0<rho<=1/8 required')
    return 9*F(68,5)*response_coefficient(radius)*r*r


def joint_mesh(rho, radius=F(1,2), leading_coefficient=False):
    r=rational(rho,'rho');b=rational(radius,'radius')
    if not 0<r<=F(1,8) or not 0<=b<1 or type(leading_coefficient) is not bool:
        raise ValueError('valid radius/rho and strict Boolean mode required')
    epsilon=r**(7 if leading_coefficient else 6)
    return {'epsilon':str(epsilon),'B':6,
            'linear':linear_modulus(epsilon,6),
            'hessian':hessian_modulus(epsilon,b,6),
            'mesh_is_minimum_of_both_positive_schedules':True,
            'leading_coefficient_error_after_rho6_division':str(r if leading_coefficient else F(1)),
            'mesh_float_materialized':False}


def certificate():
    r=F(1,2**36);b=F(1,2);K=domain()['exponential_norm_coefficient']
    budget=neumann_budget(r,b)
    if not budget<F(1,4):raise ArithmeticError('finite response inversion budget failed')
    return {'imaginary_modular_parameter':'1/3','fractional_inverse_power':'2/3',
            'all81_exponential_norm_coefficient':str(K),
            'hessian_upper':'26 K^2 rho^6',
            'relative_entropy_upper':'13 K^2 theta^2 rho^6/(1-b)^2',
            'constant_spinor_hessian_pi4_coefficient':'64/35',
            'original_shell_hessian_inverse_pi2_coefficient':str(F(64,35)*C_GAMMA**2),
            'original_shell_entropy_theta2_inverse_pi2_coefficient':str(F(32,35)*C_GAMMA**2),
            'response_error_coefficient_at_b_half':str(response_coefficient(b)),
            'neumann_example_rho':str(r),'neumann_example_budget':str(budget),
            'inverse_norm_upper':'272/15',
            'sample':{k:str(v) for k,v in estimates(F(1,64),F(1,4)).items()},
            'symbolic_mesh':joint_mesh(F(1,8),b,True),
            'fixed_L_required_for_o_R4':True,
            'varying_L_sufficient_condition':'R^2/L^6 -> 0'}
