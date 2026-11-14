"""Integer-outward bounds for prepared Gaussian/rotor comparisons."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import json
import math

import quantum_state as qs

iv=qs.nf.iv
HERE=Path(__file__).resolve().parent


def upper(x):return F(x.hi,iv.SCALE)
def lower(x):return F(x.lo,iv.SCALE)


def reciprocal(x):
    if x.lo<=0:raise ValueError('positive interval required')
    return iv.Interval(iv.SCALE**2//x.hi,iv.ceiling(iv.SCALE**2,x.lo))


def square_root(x):
    if x.lo<0:raise ValueError('nonnegative interval required')
    return iv.Interval(iv.sqrt_fraction(lower(x)).lo,iv.sqrt_fraction(upper(x)).hi)


def exponential_point(value):
    value=F(value);halvings=0
    while abs(value)>F(1,4):value/=2;halvings+=1
    term=iv.Interval.point(1);total=term
    for k in range(1,65):term=term*value/k;total+=term
    remainder=F(1,4)**65/factorial(65)/(1-F(1,4*66))
    error=iv.ceiling(remainder.numerator*iv.SCALE,remainder.denominator)
    total=iv.Interval(total.lo-error,total.hi+error)
    for _ in range(halvings):total=total*total
    return total


def exponential(x):
    return iv.Interval(exponential_point(lower(x)).lo,exponential_point(upper(x)).hi)


def parameters(eps2,N):
    eps2=F(eps2);delta=F(1,N);A=F(9,2)
    sqA=iv.sqrt_fraction(A);scale=square_root(sqA)
    kappa=iv.sqrt_fraction(1-A*delta*delta/4)
    invk=reciprocal(kappa)
    diag=square_root((invk+1)/2)
    off=square_root((invk-1)/2)
    return {'eps':iv.sqrt_fraction(eps2),'vq':sqA*invk/2,
            'vp':invk*reciprocal(sqA)/2,'cov':delta*sqA*invk/4,
            'mq':scale*diag,'mp':off*reciprocal(scale),'kappa':kappa,
            's':delta*sqA,'eps2':eps2,'delta':delta}


def preparation_comparison(q,eps2,N,cap):
    """Finite integer target versus normalized infinite Gaussian samples."""
    target,values,phases,_=qs.gaussian_cycle(q,eps2,N,cap=cap)
    p=parameters(eps2,N);error2=F();normal_lower=F()
    for z,value,phase in zip(qs.centered(q),values,phases):
        z=int(z)
        if abs(z)>cap:continue
        offset=z*p['eps']-p['mp']
        amplitude=exponential(-(offset*offset)*reciprocal(p['vp'])/4)
        phase_exact=-(p['cov']*offset*offset*reciprocal(p['vp'])/2)-p['mq']*offset
        amplitude_error=(iv.Interval.point(F(value,1<<44))-amplitude).absolute_upper()
        phase_error=(iv.Interval.point(F(float(phase)))-phase_exact).absolute_upper()
        error=amplitude_error+upper(amplitude)*phase_error
        error2+=error*error
        normal_lower=max(normal_lower,lower(amplitude))
    # A monotone Gaussian sum is bounded by its two continuous tails.
    # No floating erfc or underflow is used in this positive upper bound.
    a=F(eps2)/(2*upper(p['vp']))
    shift=upper(p['mp']*reciprocal(p['eps']))
    distance=F(cap)-abs(shift)
    if distance<=0:raise ValueError('support must contain the Gaussian center')
    tail=upper(exponential_point(-a*distance*distance))/(a*distance)
    bound=2*upper(iv.sqrt_fraction(error2+tail))/normal_lower
    return {'norm_error_upper_exact':str(bound),'norm_error_upper':math.nextafter(float(bound),math.inf),
            'unscaled_sample_tail_squared_norm_upper_exact':str(tail),
            'interior_amplitude_phase_error_squared_upper_exact':str(error2),
            'source_target':'exact normalized integer amplitudes and literal phase array; no rounded Gaussian assumed exact'}


def analytic_comparison(eps2,N,macros=1):
    p=parameters(eps2,N);eps2=F(eps2);eta=F(1,16)
    L2=4*iv.pi()*iv.pi()/eps2
    c=lower(p['vp']*L2/2)
    ec=upper(exponential_point(-c));e3c=upper(exponential_point(-3*c))
    tau=2*ec/(1-e3c)
    Nmin=lower(iv.sqrt_fraction(1-tau))
    geom=1+2*upper(exponential(-eta*L2/4))/(1-upper(exponential(-2*eta*L2)))
    v=upper(p['vq']);D=1-2*eta*v
    assert D>0
    mq,mp=p['mq'],p['mp'];moment_sum=F()
    for _ in range(N*macros):
        mean=mq.absolute_upper()/D;variance=v/D
        moment=mean**8+28*mean**6*variance+210*mean**4*variance**2+420*mean**2*variance**3+105*variance**4
        moment*=upper(exponential_point(eta*mq.absolute_upper()**2/D))/lower(iv.sqrt_fraction(D))
        moment_sum+=upper(iv.sqrt_fraction(geom*moment))
        mq,mp=(1-F(9,2)*p['delta']**2)*mq+F(9,2)*p['delta']*mp,mp-p['delta']*mq
    nonlinear=p['delta']*eps2*moment_sum/(24*Nmin)
    norm_drift=2*tau/((1+Nmin)*Nmin)
    s=upper(p['s']);s2=F(9,2)*p['delta']**2
    r=s/(4*(1-s2/4))
    phase=N*macros*s*s2/(24*lower(p['kappa'])**3)
    splitting=r*upper(iv.sqrt_fraction(F(17,4)))+phase*upper(iv.sqrt_fraction(F(3,2)))
    return {'cosine_to_normalized_periodized_Floquet_Gaussian_exact':str(nonlinear+norm_drift),
            'cosine_to_normalized_periodized_Floquet_Gaussian_upper':math.nextafter(float(nonlinear+norm_drift),math.inf),
            'substep_and_squeeze_to_continuous_oscillator_exact':str(splitting),
            'substep_and_squeeze_to_continuous_oscillator_upper':math.nextafter(float(splitting),math.inf),
            'periodization_norm_squared_error_exact':str(tau),'weighted_periodization_eta':str(eta)}


def run():
    paths=[Path(__file__),Path(qs.__file__),Path(iv.__file__),HERE/'QUANTUM_CYCLE.json']
    pins={str(p.relative_to(qs.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    rows=[]
    for q,eps2,N in ((16,F(1,3),8),(128,F(1,12),8),(512,F(1,48),8),
                     (2048,F(1,192),8),(2048,F(1,192),64)):
        cap=min(q//2-1,math.ceil(6/math.sqrt(float(eps2))))
        prep=preparation_comparison(q,eps2,N,cap)
        continuum=analytic_comparison(eps2,N)
        par=qs.parameters(eps2,N)
        register=qs.cutoff.transfer_bound(q,N,par['beta'],cap,0)
        total=F(prep['norm_error_upper_exact'])+register+F(continuum['cosine_to_normalized_periodized_Floquet_Gaussian_exact'])
        rows.append({'modulus':q,'epsilon_squared':str(eps2),'substeps':N,'initial_cap':cap,
                     'preparation':prep,'analytical_state_transfer':continuum,
                     'register_transfer_exact':str(register),
                     'complete_finite_state_to_periodized_Gaussian_norm_upper':math.nextafter(float(min(F(2),total)),math.inf),
                     'complete_finite_state_to_periodized_Gaussian_norm_exact':str(min(F(2),total))})
    assert rows[-1]['complete_finite_state_to_periodized_Gaussian_norm_upper']<.04
    assert all(F(row['preparation']['norm_error_upper_exact'])<F(1,10**5) for row in rows)
    # Exponential enclosures must contain exact x=0 and obey e^-x e^x=1.
    for x in (F(),F(1,3),F(6),F(-20)):
        product=exponential_point(x)*exponential_point(-x)
        assert lower(product)<=1<=upper(product)
    assert pins=={str(p.relative_to(qs.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result={'schema':'oph.local12.prepared_quantum_field_bounds.v1','checks':{
        'integer_outward_exponential_and_algebraic_parameters':True,
        'actual_integer_amplitude_and_phase_target_compared':True,
        'infinite_electric_Gaussian_tails_retained':True,
        'periodization_not_assumed_isometric':True,
        'cosine_global_phase_retained':True,
        'separate_Floquet_and_continuous_time_errors':True,
        'same_source_finite_register_bound_composed':True,'source_bytes_unchanged':True},
        'cases':rows,
        'norm_and_measurement_scope':'Norm transfer to normalized periodized Gaussian in the common integer electric basis; arbitrary common terminal instruments contract trace distance. Continuous angle-bin comparison additionally uses the quadrature/folding estimate in PROOF.md.',
        'nonclaims':'The prepared mode is a single optical face loop. Weak-coupling epsilon^2 tends to zero. No fixed-coupling homogeneous photon state, interacting Maxwell field, or native compiler bound at large modulus is supplied.',
        'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'GAUSSIAN_BOUNDS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'cases':[{'q':r['modulus'],'eps2':r['epsilon_squared'],'N':r['substeps'],
        'total_state':r['complete_finite_state_to_periodized_Gaussian_norm_upper'],
        'split':r['analytical_state_transfer']['substep_and_squeeze_to_continuous_oscillator_upper']} for r in rows]}),flush=True)


if __name__=='__main__':run()
