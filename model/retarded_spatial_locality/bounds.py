"""Exact constants for upper-half-plane inverse locality and Taylor bounds."""
from fractions import Fraction as F
from math import factorial


def _rational(value,name,*,positive=False):
    if type(value) not in (int,F) or value<0 or (positive and value==0):
        raise ValueError('exact '+('positive ' if positive else 'nonnegative ')+name+' required')
    return F(value)


def _order(value,maximum=12):
    if type(value) is not int or not 0<=value<=maximum:
        raise ValueError('integer order0..'+str(maximum)+' required')
    return value


def geometric_power_sum(q,order):
    """Exact sum_(d>=1) d^order q^d, using Eulerian polynomials."""
    q=_rational(q,'geometric ratio');order=_order(order)
    if q>=1:raise ValueError('geometric ratio must be below one')
    coefficients=[1]
    for n in range(order):
        updated=[0]*(len(coefficients)+1)
        # A_(n+1)=(1+n q) A_n+q(1-q) A_n'.
        for k,coefficient in enumerate(coefficients):
            updated[k]+=(k+1)*coefficient
            updated[k+1]+=(n-k)*coefficient
        while len(updated)>1 and updated[-1]==0:updated.pop()
        coefficients=updated
    numerator=sum(F(value)*q**k for k,value in enumerate(coefficients))
    return q*numerator/(1-q)**(order+1)


def shell_moment_upper(rate,order):
    """Bound sum_Z3 ||x||_inf^order exp(-rate||x||_inf).

    For order zero the origin contributes one. For positive order its
    contribution is zero. exp(-rate)<=1/(1+rate) is a rigorous bound.
    """
    rate=_rational(rate,'spatial decay rate',positive=True);order=_order(order,10)
    q=1/(1+rate)
    return (1 if order==0 else 0)+24*geometric_power_sum(q,order+2)+2*geometric_power_sum(q,order)


def weighted_perturbation_sum_upper(rate):
    """Bound sum_(d>=1)(24d²+2)d(1+d)² exp(-rate d)."""
    rate=_rational(rate,'spatial decay rate',positive=True)
    q=1/(1+rate)
    return sum(F(coefficient)*geometric_power_sum(q,order)
               for order,coefficient in ((5,24),(4,48),(3,26),(2,4),(1,2)))


def constants(eta,modulus,speed=F(1,8),interaction=F(1)):
    eta=_rational(eta,'imaginary-frequency lower bound',positive=True)
    z=_rational(modulus,'frequency-modulus upper bound',positive=True)
    c=_rational(speed,'source speed',positive=True)
    u=_rational(interaction,'contact coefficient')
    if z<eta:raise ValueError('frequency domain is empty: modulus is below eta')
    mu=min(F(1,8),eta/(384*c))
    a=192*c/eta**2+8/eta
    b=(z+2*c+u)**2/eta
    weighted_sum=weighted_perturbation_sum_upper(mu)
    alpha=min(mu,1/(2*b*a*weighted_sum))
    if alpha*a*weighted_sum>1/(2*b):
        raise ArithmeticError('weighted perturbation exceeds its Neumann budget')
    return {'eta':eta,'modulus':z,'speed':c,'interaction':u,
            'G_prefactor':a,'G_half_decay':mu,'inverse_norm':b,
            'weighted_sum':weighted_sum,'inverse_decay':alpha,
            'inverse_prefactor':2*b,
            'memory_prefactor':z+u+2*b+2*c,
            'memory_operator_norm':u*u/(4*eta)}


def spatial_taylor_remainder(eta,modulus,increment_l1,order,
                              speed=F(1,8),interaction=F(1)):
    """Absolute operator error for real spatial momentum increment q.

    The polynomial is in total degree <=order; all coefficients are exact
    derivatives of the infinite-ground inverse response. No coefficient
    values or real-axis continuation are supplied by this bound.
    """
    q=_rational(increment_l1,'real momentum increment l1 norm')
    order=_order(order,8);value=constants(eta,modulus,speed,interaction)
    return (value['inverse_prefactor']*q**(order+1)/factorial(order+1)
            *shell_moment_upper(value['inverse_decay'],order+1))


def frequency_taylor_remainder(eta,modulus,increment_modulus,order,
                                speed=F(1,8),interaction=F(1)):
    """Cauchy remainder at real k, within the same upper-half-plane domain."""
    eta=_rational(eta,'imaginary-frequency lower bound',positive=True)
    modulus=_rational(modulus,'frequency-modulus upper bound',positive=True)
    if modulus<eta:raise ValueError('frequency domain is empty: modulus is below eta')
    w=_rational(increment_modulus,'complex frequency increment modulus')
    order=_order(order,8)
    radius=eta/2
    if w>=radius:raise ValueError('frequency increment must be inside eta/2')
    larger=constants(eta/2,F(modulus)+eta/2,speed,interaction)
    ratio=w/radius
    return larger['inverse_norm']*ratio**(order+1)/(1-ratio)


def certificate():
    examples=[]
    for eta,z in ((F(1),F(2)),(F(1,2),F(2)),(F(1,4),F(2))):
        values=constants(eta,z)
        examples.append({k:str(v) for k,v in values.items()})
    return {'domain':'Im z>=eta>0 and |z|<=Z at fixed spacing, c>0,U>=0',
            'distance':'d=||x-y||_infinity in lattice sites',
            'direct_kernel':'||G_xy(z)||<=A(1+d)^2 exp(-2 mu d)',
            'inverse_kernel':'||(G(z)^-1)_xy||<=2B exp(-alpha d)',
            'memory_kernel':'||Sigma_xy(z)||<=(Z+U+2B+2c) exp(-alpha d)',
            'moment_bound':'||L V||<=sqrt3 c+U<=2c+U',
            'inverse_bound':'-Im G(z)>=eta/(Z+2c+U)^2 I and ||G(z)^-1||<=B',
            'surjectivity_proof':'The same strictly signed imaginary part bounds G and G^dagger below; the range is closed and has zero orthogonal complement.',
            'bounded_weight_proof':'Conjugate with exp(alpha f) for bounded one-Lipschitz f; Schur norm of the difference<=alpha A S(mu)<=1/(2B). Choose f(x)=min(||x-y||_infinity,R), then let R exceed the inspected distance.',
            'Fourier_symbol_domain':'Translation-invariant kernels give analytic 2x2 inverse-response symbols for ||Im k||_1<alpha.',
            'real_spatial_remainder':'2B ||q||_1^(p+1) shell_moment_upper(alpha,p+1)/(p+1)!',
            'complex_frequency_remainder':'B_star (2|w|/eta)^(p+1)/(1-2|w|/eta), |w|<eta/2, B_star at(eta/2,Z+eta/2)',
            'rational_examples':examples,
            'constants_are_conservative':True,
            'eta_zero_uniformity_claimed':False,
            'limiting_state_or_symbol_values_computed':False,
            'common_Lorentz_metric_or_Einstein_premise_completed':False,
            'nontrivial_continuum_scattering_claimed':False}
