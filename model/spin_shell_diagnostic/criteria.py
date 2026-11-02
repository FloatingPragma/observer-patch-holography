"""Exact positive-measure tests; their hypotheses are not asserted for OPH."""
from fractions import Fraction as F


def positive(value,name):
    if isinstance(value,bool) or not isinstance(value,(int,F)) or value<=0:
        raise ValueError(name+' must be a positive rational')
    return F(value)


def ratio_moments(values,weights):
    if not values or len(values)!=len(weights):raise ValueError('equal nonempty lists required')
    x=[positive(v,'ratio') for v in values]
    w=[positive(v,'weight') for v in weights]
    total=sum(w);w=[v/total for v in w]
    mean=sum(a*b for a,b in zip(x,w));inverse=sum(b/a for a,b in zip(x,w))
    variance=sum(b*(a-mean)**2 for a,b in zip(x,w));defect=mean*inverse-1
    paired=sum(w[i]*w[j]*(x[i]-x[j])**2/(x[i]*x[j])
               for i in range(len(x)) for j in range(len(x)))/2
    if defect!=paired or not min(x)**2*defect<=variance<=max(x)**2*defect:
        raise ArithmeticError('positive ratio-moment identity failed')
    return {'mean':mean,'inverse_mean':inverse,'variance':variance,'defect':defect}


def broad_band_constants(lower,upper,speed):
    a=positive(lower,'lower');b=positive(upper,'upper');c=positive(speed,'speed')
    if a>=b:raise ValueError('strictly ordered energy band required')
    return {'ratio_lower':a,'ratio_upper':F(22,7)*b,
            'variance_over_defect_upper':(F(22,7)*b)**2,
            'lattice_to_euclidean_error_coefficient':c/12,
            'first_residual_over_second_residual_upper':1/a,
            'second_residual_over_first_residual_upper':b+c}


def uniform_continuum_control():
    """A moment-level control, not a state of the Heisenberg algebra."""
    # Uniform unit probability on[1,2]. Its inverse moment is log2 in[56/81,7/10].
    lo,hi=F(56,81),F(7,10)
    assert lo+F(1,540)<hi
    mean=F(3,2);second=F(7,3);fourth=F(31,5)
    # mu_p has weight p/10; nu_p weight (5m^2/2)/p, any0<m<=1/2.
    # 1-cos(p)>=11p^2/24 for0<p<=1, and6-(1-cos(p))<=6.
    assert F(1,10)*mean<F(11,24)
    assert F(1,10)*hi<F(1,12)
    assert F(5,8)*hi<1
    assert F(5,8)*mean<F(11,2)
    assert F(1,10)<F(1,4)
    return {'uniform_scaled_energy_support':['1','2'],
            'mu_weight_over_p':'1/10','nu_weight_times_p':'5*m^2/2',
            'ward_magnetization_range':'0<m<=1/2',
            'cross_measure_density':'i*m/2',
            'positive_matrix_measure_rank_one':True,
            'known_first_inverse_moment_and_ward_bounds_satisfied':True,
            'scaled_energy_variance':str(second-mean**2),
            'minimum_squared_energy_shell_residual':str(fourth-second**2),
            'ratio_defect_lower':str(mean*lo-1),
            'spectral_atoms':False,'is_an_actual_Heisenberg_ground_state':False}
