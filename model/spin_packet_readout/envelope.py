"""One explicit band-limited packet family, including computable tail bounds."""
from fractions import Fraction as F
from math import comb, factorial
from .coefficients import div, mul, pi_interval, sqrt_interval, trig_interval
from .filters import exact


def certificate():
    delta = F(1, 64)
    # Integral over the real line of (sin y/y)^16 equals pi times this rational.
    integral_over_pi = sum(F((-1)**k*comb(16, k)*(8-k)**15, factorial(15))
                           for k in range(8))
    if F(23, 24)**16 <= F(1, 2):
        raise ArithmeticError('explicit L2 normalization lower bound failed')
    # C=delta^(3/2)/(pi*integral_over_pi)^(3/2)<1/128.
    c = F(1, 128)
    l1 = 16/(7*delta)
    derivative = F(208, 21)
    lam = c*((F(7, 4)*l1+derivative)*(l1+derivative)**2)/2
    weighted_l1 = 128/(3*delta**5)
    weighted_derivative = 256/delta**4
    tail = c*((F(7, 4)*weighted_l1+weighted_derivative)
              *(weighted_l1+weighted_derivative)**2)/2
    return {'profile': 'C exp(3 i X1/4) product_j sinc(Xj/64)^8',
            'sinc_convention': 'sinc y=sin(y)/y, sinc(0)=1',
            'sinc16_integral_over_pi': str(integral_over_pi),
            'normalization': 'C=(1/64)^(3/2)/(pi*sinc16_integral_over_pi)^(3/2)',
            'L2_norm_exactly_one': True,
            'Fourier_support_radius_lower': '5/8',
            'Fourier_support_radius_squared_upper': '51/64',
            'sampling_LCU_bound_Lambda': str(lam),
            'weighted_sampling_LCU_bound_Lambda4': str(tail),
            'tail_operator_bound': 'Lambda4*p^(-3/2)/(R-1)^4 for p<=1,R>=2',
            'coefficient_data_are_computable': True,
            'arbitrary_envelope_oracle_used': False}


def coefficient_rectangle(point, p, bits=64):
    """Enclose f_p(x)/(2 Lambda p^-3/2) from the explicit envelope."""
    if len(point) != 3 or any(type(x) is not int for x in point):
        raise ValueError('an integer cubic-lattice site is required')
    p = exact(p, positive=True)
    if p > 1 or type(bits) is not int or bits < 8:
        raise ValueError('require p<=1 and at least eight precision bits')
    data = certificate()
    q = F(data['sinc16_integral_over_pi'])
    pi = pi_interval(bits+16)
    pq = q*pi[0], q*pi[1]
    root = sqrt_interval(pq[0], bits+16)[0], sqrt_interval(pq[1], bits+16)[1]
    c = div((F(1, 512), F(1, 512)), mul(pq, root))
    amplitude = mul(c, (p**3/(2*F(data['sampling_LCU_bound_Lambda'])),)*2)
    for x in point:
        argument = p*x/64
        sinc = (F(1), F(1)) if argument == 0 else div(trig_interval(argument, bits)[1], (argument, argument))
        for _ in range(3):
            sinc = mul(sinc, sinc)
        amplitude = mul(amplitude, sinc)
    phase = trig_interval(3*p*point[0]/4, bits)
    return mul(amplitude, phase[0]), mul(amplitude, phase[1])
