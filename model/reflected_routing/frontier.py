"""Exact finite tail records for the equal-slot detector in an axial preparation."""
import math

from .kernel import require


def rate(speed):
    """Asymptotic minus log probability per free update outside the scalar cone."""
    require(0<=speed<=1,'speed outside local support')
    if speed<=0.5:
        return 0.
    if speed==1:
        return math.log(16)
    s = math.sqrt((4*speed**2-1)/(1-speed**2))
    return 4*(speed*math.asinh(s)-math.asinh(s/2))


def exact_tails(waits,speeds):
    # Coefficients of 2*4^t T_t((6+z+1/z)/8), in increasing powers.
    old = [2]
    p = [1,6,1]
    result = {}
    for t in range(1,max(waits)+1):
        if t in waits:
            denominator = 4*16**t
            rows = {}
            for a,b in speeds:
                first = (a*t+b-1)//b
                numerator = sum(x*x for x in p[t+first:])
                require(0<numerator<=denominator,'invalid equal-slot tail')
                rows[f'{a}/{b}'] = dict(numerator=str(numerator),
                        measured_rate=(math.log(denominator)-math.log(numerator))/t,
                        limiting_rate=rate(a/b))
            result[str(t)] = rows
        following = [0]*(len(p)+2)
        for j,x in enumerate(p):
            following[j] += x
            following[j+1] += 6*x
            following[j+2] += x
        for j,x in enumerate(old):
            following[j+2] -= 16*x
        old,p = p,following
    return result
