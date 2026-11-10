"""Exact rational bounds under the proved even-source graph hypotheses."""
import sympy as s


def _parameters(gap,tau,coupling):
    gap,tau,coupling=map(s.sympify,(gap,tau,coupling))
    if any(value.is_Rational is not True for value in (gap,tau,coupling)):
        raise ValueError('exact rational gap, norm and hopping coefficient required')
    if gap<=0 or tau<=0:raise ValueError('positive gap and hopping norm required')
    v=abs(coupling)*tau
    if v>gap/5:raise ValueError('hopping exceeds the proved graph interval v<=gap/5')
    return gap,tau,coupling,v


def quartic_remainder(gap,tau,coupling):
    """Return exact full-Q sharp and uniform bounds for the even effective H."""
    gap,tau,coupling,v=_parameters(gap,tau,coupling)
    gamma=gap-v
    resolvent=s.cancel(v**4/(gap*(gap**2-v**2)))
    graph=s.cancel(4*v**4/gamma**3)
    normalization=s.cancel(8*v**4/gamma**3)
    sharp=s.cancel(resolvent+graph+normalization)
    uniform=s.cancel(25*v**4/gap**3)
    assert 0<=sharp<=uniform
    return {'gap':gap,'tau':tau,'coupling':coupling,'v':v,'separated':gamma,
        'even_resolvent_error':resolvent,'graph_nonlinearity_error':graph,
        'canonical_normalization_error':normalization,'sharp':sharp,'uniform':uniform}


def uniform_quartic_remainder(gap,tau,coupling):
    """Stable consumer API: 25*(abs(coupling)*tau)^4/gap^3, v<=gap/5."""
    return quartic_remainder(gap,tau,coupling)['uniform']


def choose_quartic_coupling(gap,tau,second_order_gap):
    """Largest power-of-two coefficient with R4<=lambda^2*second_order_gap/4."""
    gap,tau,_,_=_parameters(gap,tau,s.Integer(0))
    delta=s.sympify(second_order_gap)
    if delta.is_Rational is not True or delta<=0:
        raise ValueError('positive exact rational second-order gap required')
    candidate=s.Integer(1)
    threshold=s.cancel(gap**3*delta/(100*tau**4))
    while candidate*tau>gap/5 or candidate**2>threshold:candidate/=2
    error=uniform_quartic_remainder(gap,tau,candidate)
    assert 0<candidate and error<=candidate**2*delta/4
    assert candidate**2*delta-2*error>=candidate**2*delta/2>0
    return candidate


def coefficient_certificate():
    u=s.Symbol('u',nonnegative=True)
    coefficient=1/(1-u**2)+12/(1-u)**3
    assert s.diff(coefficient,u)==2*u/(1-u**2)**2+36/(1-u)**4
    endpoint=s.factor(coefficient.subs(u,s.Rational(1,5)))
    assert endpoint==s.Rational(1175,48)<25
    return {'dimensionless_ratio_domain':['0','1/5'],
        'sharp_quartic_coefficient':'1/(1-u^2)+12/(1-u)^3',
        'positive_derivative':'2u/(1-u^2)^2+36/(1-u)^4',
        'endpoint_coefficient':str(endpoint),'uniform_coefficient':'25',
        'all_complement_states_retained':True,'formal_series_assumed':False}
