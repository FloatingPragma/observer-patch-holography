"""Exact entropy bounds and symbolic positive lattice refinement schedules."""
from fractions import Fraction as F

C_GAMMA = F(4249388, 150088295)
K_ZERO = 2048


def rational(value, name):
    if type(value) not in (int, F):
        raise ValueError(name+' must be an exact rational, not bool/float')
    return F(value)


def integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(name+' must be an integer >= '+str(minimum))
    return value


def entropy_bounds(mass, amplitude, radius=F(1, 2)):
    m=rational(mass,'mass');t=rational(amplitude,'amplitude');b=rational(radius,'radius')
    if not 0<=m<=F(1,4) or not abs(t)<=b<1:
        raise ValueError('0<=mass<=1/4 and |amplitude|<=radius<1 required')
    s=1-b
    return {'hessian_zero_lower':4*m*m, 'hessian_zero_upper':m,
            'hessian_amplitude_lower':4*m*m,
            'relative_entropy_lower':2*m*m*t*t,
            'relative_entropy_upper':m*t*t/(2*s*s),
            'quadratic_remainder_upper':m*abs(t)**3/(3*s**3)}


def regularity(depth):
    n=integer(depth,'depth')
    return {'K_prefactor':2048,'K_8192_exponent':n,
            'K_12_exponent':n*(n+1)//2,
            'alpha_denominator_factor':6,'alpha_denominator_base':12,
            'alpha_denominator_exponent':n+1}


def modulus(epsilon, radius=F(1, 2), sobolev_bound=6):
    e=rational(epsilon,'epsilon');b=rational(radius,'radius');B=rational(sobolev_bound,'B')
    if not 0<e<=1 or not 0<=b<1 or B<1:
        raise ValueError('0<epsilon<=1, 0<=radius<1 and B>=1 required')
    s=1-b;tau=e*s*s/16;upper=2/(e*s*s)
    required=32*upper/(e*s**3*tau*tau)
    n=(required.numerator+required.denominator-1)//required.denominator
    endpoint=(2*tau+1/(4*upper))/(s*s)
    tail=4*upper/(s**3*n*tau*tau)
    spatial=e/8
    total=endpoint+tail+spatial
    if not total<=e/2:
        raise ArithmeticError('Hessian modulus exceeds its exact budget')
    return {'epsilon':str(e),'amplitude_radius':str(b),'B':str(B),
            'integral_lower_cut':str(tau),'integral_upper_cut':str(upper),
            'series_terms':n,'regularity':regularity(n),
            'mesh_normalization_bound':str(1/(6*B)),
            'mesh_power_numerator_over_K':str(e*s**3*tau*tau/(2048*upper*B*B)),
            'mesh_power_factor':6,'mesh_power_base':12,'mesh_power_exponent':n+1,
            'mesh_formula':'min(1/(6B),(numerator_over_K/K_N)^(6*12^(N+1)))',
            'mesh_strictly_positive':True,'mesh_float_materialized':False,
            'endpoint_error_upper':str(endpoint),'series_error_upper':str(tail),
            'spatial_error_upper':str(spatial),'hessian_error_upper':str(total),
            'relative_entropy_error_coefficient':str(e/2)}


def shell_bounds(radius_ratio):
    rho=rational(radius_ratio,'radius ratio')
    if not 0<rho<=F(1,2**18):
        raise ValueError('0<R/L<=2^-18 required for this unboosted lower bound')
    m=C_GAMMA*rho**3/48
    if not 2048*rho<=F(1,128) or not C_GAMMA>F(1,64):
        raise ArithmeticError('shell central-value lower bound failed')
    return {'radius_ratio':str(rho),'mass_lower':str(m),
            'continuum_hessian_lower':str(4*m*m),
            'continuum_entropy_theta_squared_lower':str(2*m*m),
            'finite_hessian_lower_with_extra_mesh':str(m*m),
            'finite_entropy_theta_squared_lower_with_extra_mesh':str(m*m/2),
            'positive_mesh_numerator_over_B_squared':str(m/(48*K_ZERO)),
            'positive_mesh_power':72,
            'applies_to':'nine original unboosted directional shells'}


def certificate():
    return {'natural_mass_upper':'1/4','hessian_zero_lower':'4m^2',
            'hessian_zero_upper':'m','relative_entropy_lower':'2m^2 theta^2',
            'remainder_bound':'chi(0)|theta|^3/[3(1-b)^3]',
            'resolvent_comparison_constant':64,
            'resolvent_bound':'64B^2 K_N eta^alpha_N/tau^2+1/(N tau^2)',
            'modulus':modulus(F(1,2)),
            'shell':shell_bounds(F(1,2**18)),
            'sample_entropy_bounds':{k:str(v) for k,v in entropy_bounds(F(1,8),F(1,4)).items()}}
