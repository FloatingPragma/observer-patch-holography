"""Rational uniform-smearing, all-valley, cone and response calculators."""
from fractions import Fraction as F
from math import factorial


def _q(value,name,positive=False):
    if type(value) not in (int,F) or value<0 or (positive and value==0):
        raise ValueError('exact '+('positive ' if positive else 'nonnegative ')+name+' required')
    return F(value)


def time_lipschitz(gradient_norm_sum,norm=F(1),speed=F(1),interaction_upper=F(1)):
    return (_q(speed,'speed',True)*_q(gradient_norm_sum,'gradient norm sum')+
            _q(interaction_upper,'interaction bound')*_q(norm,'smearing norm'))


def valley_overlap(spacing,first_norm,second_norm,first_derivative,second_derivative):
    a=_q(spacing,'spacing',True)
    x,y=_q(first_norm,'first norm'),_q(second_norm,'second norm')
    dx,dy=_q(first_derivative,'first derivative'),_q(second_derivative,'second derivative')
    return a*(dx*y+x*dy)/2


def central_difference_error(spacing,third_derivative_norm):
    return _q(spacing,'spacing',True)**2*_q(third_derivative_norm,'third derivative norm')/6


def free_dynamics_error(spacing,time,third_derivative_norm_sum,norm=F(1),
                        speed=F(1),interaction=F(1)):
    return _q(time,'absolute time')*(_q(interaction,'interaction')*_q(norm,'norm')+
        _q(speed,'speed',True)*central_difference_error(spacing,third_derivative_norm_sum))


def window_response(eta,frequency_radius,momentum_radius,speed=F(1),interaction_upper=F(1)):
    eta=_q(eta,'imaginary part',True);z=_q(frequency_radius,'frequency radius',True)
    p=_q(momentum_radius,'momentum radius');v=_q(speed,'speed',True)
    u=_q(interaction_upper,'interaction bound')
    if z<eta:raise ValueError('frequency domain is empty')
    return {'inverse_norm_upper':z+v*p+u/2+u*u/(4*eta),
            'negative_imaginary_lower':eta/(z+v*p+u)**2,
            'Green_norm_upper':1/eta,
            'full_source_second_moment_upper':(v*p+u)**2}


def cone_bound(spacing,inner_radius,outer_radius,time,speed=F(1),
               first_norm=F(1),second_norm=F(1)):
    a=_q(spacing,'spacing',True);r=_q(inner_radius,'inner radius')
    R=_q(outer_radius,'outer radius',True);t=_q(time,'time')
    v=_q(speed,'speed',True);x=_q(first_norm,'first norm');y=_q(second_norm,'second norm')
    if R<=r or a>min(F(1),(R-r)/8):raise ValueError('spacing does not resolve separated cubes')
    delta=(R-r)/2-48*v*t
    if delta<=0:raise ValueError('supports are not outside the proved propagation cone')
    r_lattice=-int((-r)//a)+1
    R_lattice=int(R//a)-1
    if R_lattice<=r_lattice:raise ArithmeticError('lattice buffer vanished')
    coefficient=192*v*t*(2*r+5)**3*(2*R+1)**2*x*y
    return {'inner_lattice_radius':r_lattice,'outer_lattice_radius':R_lattice,
            'exponential_margin':delta,'prefactor':coefficient,
            'linear_in_spacing_upper':min(2*x*y,coefficient*factorial(7)*a/delta**7),
            'exponential_power_of_spacing':-6,'Taylor_order':7}


def certificate():
    response=window_response(F(1),F(2),F(3),F(2),F(1))
    rows=[]
    for a in (F(1,16),F(1,32),F(1,64)):
        row=cone_bound(a,F(1),F(10),F(1,100),F(1))
        rows.append({'spacing':str(a),**{k:str(v) for k,v in row.items()}})
    return {'all_eight_valleys_retained':True,
            'kinetic_limits':['i v sum_j (-1)^b_j sigma_j partial_j, b in {0,1}^3'],
            'time_Lipschitz':'v sum_j ||partial_j f||+Umax||f||',
            'propagation_speed_sup_norm':'96v',
            'cone_examples':rows,
            'window_parameters':{'eta':'1','Z':'2','P':'3','v':'2','Umax':'1'},
            'window_bounds':{k:str(v) for k,v in response.items()},
            'spacing_uniform_for_bounded_U':True,
            'full_CAR_Duhamel_bound':'|t| (U(a)||f||+v a^2 sum_j||partial_j^3 f||/6)',
            'vanishing_U_free_eight_valley_continuum_proved':True,
            'real_frequency_boundary_control':False,
            'inverse_weak_limit_interchanged':False,
            'Lorentz_metric_or_Einstein_premise_completed':False}
