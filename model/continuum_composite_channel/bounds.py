"""Exact degree, Dyson, centrality and ground-fluctuation bounds."""
from fractions import Fraction as F
from math import ceil, factorial, isqrt


def rational(value, *, positive=False):
    if type(value) not in (int,F) or value<0 or (positive and value==0):
        raise ValueError('nonnegative exact rational required')
    return F(value)


def integer(value, minimum=1):
    if type(value) is not int or value<minimum:raise ValueError('integer outside proved domain')
    return value


def sqrt_upper(value,bits=80):
    x=rational(value);integer(bits)
    if bits>256:raise ValueError('precision exceeds calculator domain')
    scale=1<<bits;n=isqrt(x.numerator*scale*scale//x.denominator)
    if F(n*n,scale*scale)<x:n+=1
    return F(n,scale)


def dyson_remainder(degree,order,coupling,time):
    m=integer(degree);N=integer(order);U=rational(coupling);t=rational(time)
    rising=F(1)
    for j in range(N):rising*=F(m,2)+j
    return (2*U*t)**N*rising/factorial(N)


def approximation(degree,coupling,time,precision):
    m=integer(degree);U=rational(coupling);T=rational(time);ell=integer(precision)
    if not U or not T:return {'steps':0,'degree':m,'norm_upper':F(1),'error_upper':F(0)}
    K=max(1,ceil(8*U*T));norm=(1+F(1,2**ell))**K
    return {'steps':K,'degree':3**K*(m+ell-1)-(ell-1),
            'norm_upper':norm,'error_upper':norm-1}


def commutator_upper(degree,coupling,time,precision,coefficient_sup,operator_norm):
    a=approximation(degree,coupling,time,precision)
    small=rational(coefficient_sup);norm=rational(operator_norm)
    return 2*a['degree']*small*a['norm_upper']+2*norm*a['error_upper']


def mixed_defect_ground_norm_squared(defect,l2_product,I4):
    D=rational(defect);L=rational(l2_product);four=rational(I4)
    if D>2:raise ValueError('occupation defect exceeds its universal range')
    return 4*(four+(8*sqrt_upper(D)+D/2)*L*L)


def smearing_bounds(spacing,f_l2,g_l2,f_sup,g_sup):
    a=rational(spacing,positive=True)
    nf,ng,sf,sg=map(rational,(f_l2,g_l2,f_sup,g_sup))
    return {'l2_product':nf*ng,'I4_upper':a**3*min(sf*sf*ng*ng,sg*sg*nf*nf),
            'quadratic_coefficient_sup':a**3*sf*sg}


def precision_for_spacing(spacing,steps):
    a=rational(spacing,positive=True);K=integer(steps)
    ell=1
    while F(1,2**ell)>a**3/(2*K):ell+=1
    return ell


def certificate():
    examples=[]
    for power in (4,8,12,16):
        a=F(1,2**power);K=8;ell=precision_for_spacing(a,K)
        smearing=smearing_bounds(a,1,1,1,1)
        app=approximation(1,1,1,ell)
        examples.append({'spacing':str(a),'precision':ell,
                         **{k:str(v) for k,v in app.items()},
                         'mixed_defect_commutator_upper':str(commutator_upper(1,1,1,ell,2*smearing['quadratic_coefficient_sup'],3))})
    return {'Dyson_remainder':'(2 Umax |t|)^N (m/2)_N/N!',
            'Bernstein_degree_bound':'||[dGamma(q),A]||<=m||q||||A|| for q=q*, degree(A)<=m',
            'contact_derivation_degree_increment':2,
            'contact_derivation_norm_bound':'||[sum n_up n_down,A]||<=m||A||',
            'quadratic_fluctuation_upper':'4[I4+(8 sqrt(D)+D/2)||u||^2||v||^2]',
            'fixed_time_commutator_rate':'O(a^3 log(1/a)) for smooth bounded physical smears',
            'examples':examples,
            'positive_quadratic_energy_escape_control':{str(L):{k:str(v) for k,v in ultraviolet_channel(1,F(1,2),L,1,2).items()} for L in (8,32,128)},
            'finite_degree_approximation_uniform_in_spacing_and_volume':True,
            'composite_time_continuity_proved':False,
            'energy_tightness_of_composite_proved':False}


def ultraviolet_channel(system_energy,coupling,bath_energy,time,frequency_radius):
    """Exact bounds for a positive quadratic two-mode escape control."""
    E=rational(system_energy,positive=True);k=rational(coupling,positive=True)
    L=rational(bath_energy,positive=True);T=rational(time);Z=rational(frequency_radius)
    if L<=E or E*L<=k*k or L<=Z:
        raise ValueError('control requires Lambda>E0, E0*Lambda>kappa^2, and Lambda>Z')
    return {'memory_Gram':k*k,
            'selfenergy_upper_on_disk':k*k/(L-Z),
            'finite_time_field_error_upper':(2*k+k*k*T)/(L-E),
            'positive_one_particle_determinant':E*L-k*k}
