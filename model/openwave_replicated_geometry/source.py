"""Exact binomial spectral sums, without discarding any register sector."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from model.finite_area_response.response import Interval,exponential
from model.openwave_universal_energy_geometry.source import QSIZE as Q, OMEGA, sinhcosh,root_interval,native

ONE=Interval.point(1)

def power(x,n):
    out=ONE
    while n:
        if n%2:out=out*x
        x=x*x;n//=2
    return out

def sum_weights(n,x):
    term=ONE;total=Interval.point(0);rows=[]
    for r in range(1,n+1):
        term=term*x*F(n-r+1,r)
        value=term/r;total=total+value
        rows.append({'excited_registers':r,'multiplicity':comb(n,r),'response_weight_interval':value.encoded()})
    return total,rows

@lru_cache(None)
def frequency(n):
    if n<1:raise ValueError('positive register count required')
    a=Q/n;sa,ca=sinhcosh(a);sq,_=sinhcosh(Q)
    # This analytic source-to-source normalization matches the original
    # finite-q far-field susceptibility, without using observed Newton G.
    return OMEGA*n*power(ca,n)*(sa/ca)/sq

def susceptibility(n,f):
    f=F(f)
    if not 0<=f<=1:raise ValueError('native profile must lie in[0,1]')
    a,b=Q/n,Q*f/n;sa,ca=sinhcosh(a);sb,cb=sinhcosh(b)
    x=sa/ca*sb/cb;total,rows=sum_weights(n,x);sq,_=sinhcosh(Q)
    # Algebraically cancel the common source normalization before enclosing;
    # this is the actual zero-source log-lapse derivative.
    slope=-2*sq*total/(OMEGA*(sa/ca))
    return slope,rows

def vacuum_records(n,f):
    f=F(f);b=Q*f/n;_,cb=sinhcosh(b);_,c2b=sinhcosh(2*b)
    mean=power(cb,n);relative=power(c2b,n)/(mean*mean)-1
    # The exact coordinate readout has n+1 retained multiplicity outcomes.
    coordinates=[{'up_registers':k,'Q_value_exact':str(Q*F(2*k-n,n)),
                  'probability_exact':str(F(comb(n,k),2**n))}for k in range(n+1)]
    assert sum(F(r['probability_exact'])for r in coordinates)==1
    return {'lapse_mean_interval':mean.encoded(),'relative_lapse_variance_interval':relative.encoded(),
            'relative_variance_upper_interval':(exponential(Q*Q*f*f/n)-1).encoded(),
            'coordinate_outcomes':coordinates,'nonselective_coordinate_measurement_work_interval':(frequency(n)/2).encoded()}

def classical_value(z,omega,E,P,f):
    z=F(z)
    return omega*(ONE-root_interval(Interval.point(1-z*z)))/2+E*exponential(Q*z)+P*exponential(Q*f*z)

def classical_derivative(z,omega,E,P,f):
    z=F(z)
    return omega*z/(2*root_interval(Interval.point(1-z*z)))+E*Q*exponential(Q*z)+P*(Q*f)*exponential(Q*f*z)

def classical_minimum(omega,E,P,f):
    lo,hi=F(-1),F(0)
    for _ in range(60):
        mid=(lo+hi)/2;derivative=classical_derivative(mid,omega,E,P,f)
        if derivative.hi<0:lo=mid
        elif derivative.lo>0:hi=mid
        else:break
    mid=(lo+hi)/2;value=classical_value(mid,omega,E,P,f)
    derivative=classical_derivative(mid,omega,E,P,f)
    error=max(abs(derivative.lo),abs(derivative.hi))*(hi-lo)/2
    return Interval(lo,hi),Interval(value.lo-error,value.hi)

def classical_limit_record(n,E,P,f):
    omega=frequency(n);sq,_=sinhcosh(Q);limit_omega=Interval.point(OMEGA*Q)/sq
    location,energy=classical_minimum(omega,E,P,f)
    _,limit_energy=classical_minimum(limit_omega,E,P,f)
    error=(E*exponential(Q)+P*exponential(Q*f))*(exponential(Q*Q/(2*n))-1)
    # Uniform strong convexity of E exp(q z)+P exp(qF z).
    curvature=E*Q*Q*exponential(-Q)+P*(Q*f)**2*exponential(-Q*f)
    variance=2*error/curvature
    displacement=root_interval(4*error/omega)
    frequency_error=max(abs(omega.lo-limit_omega.hi),abs(omega.hi-limit_omega.lo))/2
    return {'registers':n,'source_energy_interval':E.encoded(),'probe_energy_interval':P.encoded(),
        'profile_exact':str(f),'classical_minimum_coordinate_interval':location.encoded(),
        'finite_N_classical_minimum_energy_interval':energy.encoded(),
        'complete_quantum_ground_energy_enclosure':[str(energy.lo),str(energy.hi+error.hi)],
        'coherent_product_trial_energy_excess_upper':str(error.hi),
        'ground_mean_coordinate_displacement_upper':str(displacement.hi),
        'ground_coordinate_variance_upper':str(min(F(1),variance.hi)),
        'limiting_classical_minimum_energy_interval':limit_energy.encoded(),
        'energy_distance_from_limiting_classical_minimum_upper':str(frequency_error+error.hi),
        'all_2_to_N_register_states_and_native_energy_sectors_retained':True}

def run():
    records=[];sq,_=sinhcosh(Q)
    target=-2*Q*sq/OMEGA
    for n in (1,2,4,8,16,32,64,128,256):
        omega=frequency(n);responses=[]
        for f in (F(0),F(1,16),F(1,4),F(1,2),F(1)):
            slope,weights=susceptibility(n,f)
            if f:
                ratio=slope/(target*f)
                lower=1-Q*Q*f*f/(3*n*n);upper=exponential(Q*Q*f/n)
                assert ratio.hi>=lower and ratio.lo<=upper.hi
                # Intersect two independently justified rational enclosures.
                ratio=Interval(max(ratio.lo,lower),min(ratio.hi,upper.hi))
                relative_error=Interval(F(0),max(1-lower,upper.hi-1))
            else:ratio=ONE;relative_error=Interval.point(0)
            responses.append({'source_profile_exact':str(f),'native_log_lapse_source_derivative_interval':slope.encoded(),
                'ratio_to_preserved_linear_profile_interval':ratio.encoded(),
                'analytic_relative_error_bound_interval':relative_error.encoded(),
                'all_excitation_weights':weights})
        spectrum=[{'excited_registers':r,'energy_interval':(omega*F(r,n)).encoded(),'multiplicity':comb(n,r)}for r in range(n+1)]
        assert sum(x['multiplicity']for x in spectrum)==2**n
        records.append({'registers':n,'complete_geometry_dimension':2**n,'full_native_source_probe_dimension':64*2**n,
            'geometry_frequency_numerator_interval':omega.encoded(),'geometry_gap_interval':(omega/n).encoded(),
            'complete_free_geometry_spectrum':spectrum,'responses':responses,'vacuum_lapse':vacuum_records(n,F(1))})
    energies=native()[2];classical=[]
    for n in (1,4,16,64,256):
        for i in (1,3,7):
            for probe in (Interval.point(0),energies[1]):
                row=classical_limit_record(n,energies[i],probe,F(1,2));row['native_source_level']=i
                row['native_probe_level']=0 if probe.hi==0 else 1;classical.append(row)
    return {'schema':'oph.replicated_native_geometry.v1','verdict':'PASS_COMPLETE_REPLICATED_GEOMETRY_RESPONSE_AND_FLUCTUATION_BOUNDS',
      'architecture':{'native_energy_factors':'all eight parent source levels and all eight probe levels retained',
        'Q_N':'(q/N) sum_i Z_i','H_geometry':'omega_N/(2N) sum_i(I-X_i)',
        'complete_coupling':'H_geometry+exp(Q_N)K_source+exp(F Q_N)K_probe',
        'q_exact':str(Q),'inherited_omega_exact':str(OMEGA),
        'omega_N':'omega*N*cosh(q/N)^N*tanh(q/N)/sinh(q)',
        'normalization_ancestry':'analytic preservation of the original finite-q far-field response; no observed G or new empirical target',
        'positive_complete_Hamiltonian':True,'vacuum':'native source/probe ground times all-register X-positive product'},
      'preserved_far_field_log_lapse_slope_interval':target.encoded(),
      'complete_register_records':records,
      'complete_ground_to_classical_limit':{'classical_energy':'omega_N/2*(1-sqrt(1-z^2))+E*exp(q*z)+P*exp(q*F*z)',
        'exact_lower_bound':'quantum ground energy >= minimum of the classical energy; mean Bloch norm<=1 and spectral Jensen',
        'upper_error':'(E*exp(q)+P*exp(q*F))*(exp(q^2/(2N))-1), from a complete spin-product trial and Hoeffding',
        'strict_convexity':'classical energy second derivative >= omega_N/2 on(-1,1)',
        'limiting_frequency':'omega*q/sinh(q)',
        'records':classical,'arbitrary_finite_source_and_probe_energies_admitted':True,
        'uniform_in_unbounded_energy':False,'dynamical_or_spatial_classical_limit':False},
      'uniform_theorems':{'ratio_lower':'1-q^2*F^2/(3*N^2)',
        'ratio_upper':'exp(q^2*F/N)',
        'vacuum_relative_lapse_variance_upper':'exp(q^2*F^2/N)-1',
        'asymptotic_frequency_numerator':'omega*q/sinh(q)',
        'asymptotic_geometry_gap':'omega*q/(N*sinh(q)); tends to zero in the retained clock',
        'unscaled_geometry_control':'without the1/N geometry-energy scaling, susceptibility vanishes as1/N',
        'complete_coordinate_readout_work':'omega_N/2; readout is not energy-free'},
      'scope':{'uniform_weak_source_and_weak_probe_response':True,'vacuum_lapse_noise_vanishes':True,
        'finite_strength_universal_clock_response_proved':False,'dynamical_classical_metric_limit_proved':False,
        'common_geometry_preparation_cost_uniform_in_N':False,'physical_Newton_strength_matched':False,
        'physical_Einstein_spacetime':False,'new_empirical_calibration':False}}
