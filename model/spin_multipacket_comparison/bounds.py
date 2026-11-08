"""Exact bounds for local source creators and positive-transfer multipackets."""
from fractions import Fraction as F
from math import comb, factorial

from model.spin_charge_dispersion.order import CERTIFIED_M0
from model.spin_energy_polarization.bounds import charged_truncation


def rational(value, name, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and value == 0):
        raise ValueError(name+' must be an exact nonnegative rational')
    return F(value)


def cutoff(amplitude):
    c=rational(amplitude,'analytic Schwartz seminorm',True)
    def admissible(n):
        row=charged_truncation(F(n),6,c)
        return row['tail_L1_upper']<=F(1,4) and row['tail_L2_squared_upper']<=F(1,16)
    upper=1
    while not admissible(upper):upper*=2
    lower=0
    while lower+1<upper:
        mid=(lower+upper)//2
        if admissible(mid):upper=mid
        else:lower=mid
    row=charged_truncation(F(upper),6,c)
    return {'analytic_Schwartz_power':6,'analytic_Schwartz_seminorm':c,
            'cutoff_radius':upper,'tail_bounds':row,'overlap_lower':F(15,16),
            'configuration_variance_upper':F(5547,576),
            'variance_strictly_below':12,
            'sampled_profile_claimed_to_certify_seminorm':False}


def local_band(radius=F(1),m0=CERTIFIED_M0):
    r=rational(radius,'momentum radius',True);m=rational(m0,'order lower bound',True)
    if m>F(1,2):raise ValueError('spin order is at most one half')
    d=3*m/1000;a=3*d*d/2;b=2/(d*d)
    if a/12+1/(4*b)!=d*d/4:raise ArithmeticError('conjugate-field tail budget differs')
    return {'order_lower_bound':m,'commutator_expectation_lower':d,
            'momentum_annulus':[r/2,r],'energy_window':[r*a,r*b],
            'dimensionless_energy_window':[a,b],
            'each_charge_original_joint_band_probability_lower':d*d/2,
            'each_charge_full_first_energy_moment_upper':3*r,
            'complete_preparation_records':3,'nonlinear_creator_all_momenta_annular':False}


def smoothstep(degree):
    if type(degree)is not int or degree<2:raise ValueError('integer smoothness at least two required')
    factor=F(factorial(2*degree+1),factorial(degree)**2)
    coeff=[F(0)]*(2*degree+2)
    for j in range(degree+1):coeff[degree+j+1]=factor*(-1)**j*comb(degree,j)/F(degree+j+1)
    return tuple(coeff)


def derivative(coeff,order):
    result=tuple(coeff)
    for _ in range(order):result=tuple(i*result[i]for i in range(1,len(result)))
    return result


def derivative_l1(coeff,order):
    return sum((abs(c)/F(i+1)for i,c in enumerate(derivative(coeff,order))),F())


def filter_bounds(a,b,smoothness=8):
    a=rational(a,'inner lower energy',True);b=rational(b,'inner upper energy',True)
    if a>=b:raise ValueError('ordered positive energy band required')
    coeff=smoothstep(smoothness)
    # The compact piecewise polynomial equals one on[a,b], and zero
    # outside[a/2,2b]. Endpoint derivatives through smoothness vanish.
    hj=[2*b-a/2]
    for j in range(1,smoothness+2):
        hj.append(((a/2)**(1-j)+b**(1-j))*derivative_l1(coeff,j))
    eh=[2*b*hj[0]]+[2*b*hj[j]+j*hj[j-1]for j in range(1,len(hj))]
    mass=max(F(1),(hj[0]+hj[2])/3)
    return {'smoothness':smoothness,'transition_coefficients':coeff,
            'inner_band':[a,b],'spectral_transfer_interval':[a/2,2*b],
            'frequency_derivative_L1_bounds':hj,
            'energy_times_frequency_derivative_L1_bounds':eh,
            'time_filter_L1_upper':mass,
            'time_filter_derivative_L1_upper_per_radius':(eh[0]+eh[2])/3,
            'all_filter_bounds_analytic':True,'Fourier_samples_used_as_bounds':False}


def tail(filter_data,radius,time_cutoff,order=None,derivative_time=False):
    r=rational(radius,'momentum radius',True);t=rational(time_cutoff,'time cutoff',True)
    k=filter_data['smoothness']+1 if order is None else order
    if type(k)is not int or not 2<=k<=filter_data['smoothness']+1:
        raise ValueError('tail derivative order outside certified filter smoothness')
    key='energy_times_frequency_derivative_L1_bounds'if derivative_time else'frequency_derivative_L1_bounds'
    value=filter_data[key][k]/(3*(k-1)*(r*t)**(k-1))
    return r*value if derivative_time else value


def commutator(filter_data,radius,distance,order=None):
    r=rational(radius,'momentum radius',True);d=rational(distance,'support distance',True)
    t=d/(4*72);m=filter_data['time_filter_L1_upper']
    e=tail(filter_data,r,t,order);ep=tail(filter_data,r,t,order,True)
    derivative_norm=r*filter_data['time_filter_derivative_L1_upper_per_radius']
    return {'support_distance':d,'time_cutoff':t,
            'commutator_norm_upper':min(F(1,2),e/m),
            'one_time_derivative_commutator_upper':ep/(2*m)+derivative_norm*e/(2*m*m)}


def configuration(n,radius,profile_radius,filter_data,probability,epsilon):
    if type(n)is not int or n<1:raise ValueError('positive integer packet count required')
    r=rational(radius,'momentum radius',True);r0=rational(profile_radius,'profile cutoff radius',True)
    p=rational(probability,'original single-preparation floor',True)
    eps=rational(epsilon,'commutator error',True)
    if p>F(1,4):raise ValueError('creator probability cannot exceed one quarter')
    allowed=p**n*4**(n-1)/(4*n*(n-1))if n>1 else F(1,2)
    if eps>allowed:raise ValueError('commutator allowance cannot retain the multipacket probability')
    k=filter_data['smoothness']+1;m=filter_data['time_filter_L1_upper']
    # Choose an integer source support separation from the exact tail bound.
    rhs=filter_data['frequency_derivative_L1_bounds'][k]*(4*72/r)**(k-1)/(3*(k-1)*m*eps)
    distance=integer_root_ceiling(rhs,k-1)
    centre_distance=F(distance)+4*r0/r
    ratio=4*comb(n,2)/(4*p)**n
    domain_multiplier=integer_root_ceiling(ratio,3)
    reordering=n*(n-1)*F(4)**(1-n)*eps
    near_error=F(comb(n,2),domain_multiplier**3)*F(4)**(-n)
    lo,hi=filter_data['spectral_transfer_interval']
    return {'packet_count':n,'original_preparation_record_count':3**n,
            'all_equal_source_charges':[n,-n],
            'single_packet_original_probability_lower':p,
            'source_support_distance_required':distance,
            'centre_max_norm_separation_required':centre_distance,
            'search_box_half_width':centre_distance*domain_multiplier,
            'near_tuple_integral_error_upper':near_error,
            'commutator_reordering_error_upper':reordering,
            'rational_good_tuple_original_weight_lower':p**n/4,
            'complete_created_energy_support':[n*r*lo,n*r*hi],
            'source_charge_zero_included':False,
            'effective_good_location_or_ground_evaluation_algorithm':False,
            'particle_count_or_Fock_factorization_inferred':False}


def integer_root_ceiling(value,degree):
    value=rational(value,'root argument')
    if type(degree)is not int or degree<1:raise ValueError('positive integer root degree required')
    upper=1
    while upper**degree<value:upper*=2
    lower=0
    while lower+1<upper:
        mid=(lower+upper)//2
        if mid**degree>=value:upper=mid
        else:lower=mid
    return upper


def cook_bound(norms,vacuum_defects,pair_commutators):
    norms=tuple(rational(x,'operator norm')for x in norms)
    errors=tuple(rational(x,'vacuum derivative norm')for x in vacuum_defects)
    if len(norms)!=len(errors):raise ValueError('one vacuum derivative defect per creator required')
    n=len(norms)
    if set(pair_commutators)!={(i,j)for i in range(n)for j in range(i+1,n)}:
        raise ValueError('complete ordered derivative-commutator family required')
    def product(excluded):
        result=F(1)
        for k,v in enumerate(norms):
            if k not in excluded:result*=v
        return result
    one=sum((errors[i]*product({i})for i in range(n)),F())
    cross=sum((rational(value,'derivative commutator')*product({i,j})for(i,j),value in pair_commutators.items()),F())
    return {'one_packet_defect_term':one,'inter_packet_commutator_term':cross,'full_product_derivative_norm_upper':one+cross}
