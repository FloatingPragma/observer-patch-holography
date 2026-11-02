"""Exact constants for normalized infrared packets and a source diagonal."""
from fractions import Fraction as F


def integer(n):
    if type(n) is not int or n<1:raise ValueError('positive integer scale required')
    return n


def positive(x,name):
    if type(x) not in (int,F) or x<=0:raise ValueError(name+' must be positive exact rational')
    return F(x)


def diagonal(n,rate_over_exchange=F(1),magnetization=F(1,1000)):
    n=integer(n);rate=positive(rate_over_exchange,'rate/exchange');m=positive(magnetization,'magnetization')
    if m>F(1,2):raise ValueError('spin half magnetization required')
    R=n*n;L=128*n**3;E=6*L*(2*L+1)**2
    # Two filter times plus one scaled characteristic time |s|<=n.
    # Total time <=3 n²/J, a uniform compact-s rescaled characteristic test.
    jt=3*n*n
    # Three free segments may be separated by the two local insertions.
    delta=F(24,E**3)+F(8*rate*jt,E**2)+F(6*jt,E)
    c=m*m/96
    # ||A|| <= K n^(3/2). Bounds below are coefficients of K²/c.
    normalized_source=n**4*delta
    exponent=-L+2*R+72*jt
    prefactor=6*jt*(2*L+1)**2
    # n=1 can have a positive exponent; this is a bound, never a convergence claim at n=1.
    return {'scale':n,'momentum_radius':str(F(1,n)), 'spatial_cutoff':R,
            'boundary_support_radius':2*R,'box_radius':L,'edges':E,'total_J_time_bound':jt,
            'source_trace_error_bound':str(delta),
            'source_normalized_error_coefficient_K2_over_c':str(normalized_source),
            'spatial_time_tail_normalized_coefficient':str(F(n**4,(1+n)**8)),
            'boundary_prefactor':prefactor,'boundary_power_of_two':exponent,
            'spectral_weight_lower':str(c/n),
            'spectral_weight_upper':str(F(1,4*n)),
            'rescaled_energy_support':['1/16000','64000000'] if m==F(1,1000) else [str(m/16),str(64/m**2)],
            'source_c_over_J':E**4,'source_U_over_J':E**8,
            'finite_ground_selection_cost_known':False}


def normalized_ratio_error(numerator_error,denominator_error,denominator_lower):
    a=positive(numerator_error,'numerator error');b=positive(denominator_error,'denominator error');w=positive(denominator_lower,'weight lower')
    if b>=w:raise ValueError('normalization denominator not certified positive')
    # |N|<=w for the true characteristic function. Upper uses w>=lower.
    return (a+b)/(w-b)


def certificate():
    # Cube integration: (3/2)^8 * integral_R3 (1+||x||inf)^-8 dx.
    integral=F(8,35);constant=F(3,2)**8*integral
    if not constant<6:raise ArithmeticError('lattice L1 constant')
    # sqrt((Jp²/2)/(10J))=p/sqrt(20)<p/4.
    if not F(1,20)<F(1,16):raise ArithmeticError('upper spectral weight bound')
    # ln2>=2/3 implies exp(48Jt)<=2^(72Jt).
    # With total Jt=3n², L=128n³ the boundary exponent is <=-19n³ for n>=2.
    if any(-128*n**3+218*n*n > -19*n**3 for n in range(2,20)):
        raise ArithmeticError('boundary decay exponent')
    data=[diagonal(n) for n in (2,4,8,16,32)]
    if not all(F(a['source_normalized_error_coefficient_K2_over_c'])>F(b['source_normalized_error_coefficient_K2_over_c']) for a,b in zip(data,data[1:])):
        raise ArithmeticError('normalized source error does not decay')
    return {'lattice_integral_exact':str(integral),'lattice_sum_constant_exact':str(constant),
            'lattice_sum_constant_upper':'6','spatial_filter_norm_power':'p^(-3/2)',
            'norm_square_divided_by_weight_power':'p^(-4)',
            'time_filter_L1_independent_of_p':True,
            'filtered_weight_two_sided':'m^2*p/96 <= w_p <= p/4',
            'boundary_power_upper_n_at_least_2':'-19 n^3',
            'source_normalized_error_order':'O(n^(-3))',
            'truncation_normalized_error_order':'O(n^(-4))',
            'examples':data,
            'normalization_requires_error_below_weight':True}
