"""Rational error bounds for complete position records of one native link."""
from dataclasses import dataclass
from fractions import Fraction as F
from math import factorial


@dataclass(frozen=True)
class Interval:
    low: F
    high: F

    def __post_init__(self):
        if type(self.low) is not F or type(self.high) is not F or self.low > self.high:
            raise ValueError('ordered exact Fraction endpoints required')

    def expanded(self, error):
        if type(error) is not F or error < 0:
            raise ValueError('nonnegative exact error required')
        return Interval(self.low-error, self.high+error)

    def contains(self, value):
        return self.low <= value <= self.high

    def encoded(self):
        return [str(self.low), str(self.high)]


def probability(d, hopping, electric, time, terms=12):
    """Certified p(t) via the alternating cosine expansion, no float input."""
    if (type(d) is not int or d < 1 or type(terms) is not int or terms < 1
            or any(type(x) is not F for x in (hopping,electric,time))
            or hopping <= 0 or electric < 0 or time <= 0):
        raise ValueError('positive dimension, hopping, time and nonnegative electric energy required')
    a = d*hopping*hopping; omega2 = electric*electric+4*a
    z = omega2*time*time
    if z > 1:
        raise ValueError('the certified short-time evaluator requires Omega^2 t^2<=1')
    contrast = 4*a/omega2
    value = sum(((-1)**(k+1)*contrast*z**k/(2*factorial(2*k)) for k in range(1,terms+1)), F())
    error = contrast*z**(terms+1)/(2*factorial(2*terms+2))
    return Interval(value-error,value+error)


def infer_coefficients(d, time, p1, p2, hopping_upper, electric_upper):
    """Infer J² and epsilon² from all-outcome p(t),p(2t) error intervals.

    The a priori bounds are source specifications. Probability intervals may
    include analytical, sampling and implementation errors. No accepted-only
    frequency or laboratory clock is substituted for the complete records.
    """
    if (type(d) is not int or d < 1 or any(type(x) is not F for x in
            (time,hopping_upper,electric_upper)) or time <= 0 or hopping_upper <= 0
            or electric_upper < 0 or not isinstance(p1,Interval) or not isinstance(p2,Interval)
            or p1.low < 0 or p2.low < 0 or p1.high > 1 or p2.high > 1):
        raise ValueError('admitted source bounds and complete probability intervals required')
    amax = d*hopping_upper**2; omega2max = electric_upper**2+4*amax
    error_a = amax*omega2max**2*time**4/54
    error_b = 17*amax*omega2max**2*time**2/1080
    a = Interval((16*p1.low-p2.high)/(12*time**2),
                 (16*p1.high-p2.low)/(12*time**2)).expanded(error_a)
    b = Interval((4*p1.low-p2.high)/(12*time**4),
                 (4*p1.high-p2.low)/(12*time**4)).expanded(error_b)
    if a.low <= 0 or b.low < 0:
        raise ValueError('record precision or time resolution does not certify positive coefficients')
    omega = Interval(12*b.low/a.high,12*b.high/a.low)
    raw_e = Interval(omega.low-4*a.high,omega.high-4*a.low)
    if raw_e.high < 0:
        raise ValueError('records contradict the specified source Hamiltonian')
    energy = Interval(max(F(),raw_e.low),raw_e.high)
    j2 = Interval(a.low/d,a.high/d)
    ratio = Interval(energy.low/j2.high,energy.high/j2.low)
    return {'hopping_squared': j2, 'electric_squared': energy,
            'electric_over_hopping_squared': ratio, 'omega_squared': omega,
            'a_taylor_error': error_a, 'b_taylor_error': error_b}


def witness():
    rows=[]
    for d in (1,2,3,6):
        j, energy, time = F(1), F(1), F(1,1024)
        a, b = probability(d,j,energy,time),probability(d,j,energy,2*time)
        inferred = infer_coefficients(d,time,a,b,j,energy)
        if not all(inferred[name].contains(F(1)) for name in
                   ('hopping_squared','electric_squared','electric_over_hopping_squared')):
            raise ArithmeticError('exact source coefficient is outside the detector enclosure')
        rows.append({'block_dimension': d, 'time':str(time),
                     'complete_probabilities':[a.encoded(),b.encoded()],
                     'intervals':{k:v.encoded() if isinstance(v,Interval) else str(v)
                                  for k,v in inferred.items()},
                     'true_J_squared':'1','true_epsilon_squared':'1'})
    return rows
