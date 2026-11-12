"""Outward interval Ritz and current bounds for exact algebraic source matrices."""
from functools import lru_cache
from fractions import Fraction as F
from math import factorial
import numpy as np
import sympy as sp
from mpmath import iv
from .source import matter_exact,ladders,source,MASS,cross


def rational_bounds(cap):
    sqrt_upper={1:F(1),2:F(3,2),3:F(7,4),4:F(2)}[cap]
    a=F(8,5)+F(3,5)*sqrt_upper;x=a*a/MASS**2
    epsilon=a*a/(2*MASS)
    remainder=F(33,2048)*a**14/MASS**13
    derivative=F(2)*a*F(4,5)/MASS*F(33,2048)*x**6*(7-6*x)/(1-x)**2
    radius=F(5);series=sum((F(25**k,factorial(k))for k in range(101)),F(0))
    outside=(F(16,7)*(radius**5+F(5,2)*radius**3+5*radius)+F(40,7)/radius)/series
    return {'kinetic_norm_upper':a,'positive_kernel_shift_upper':epsilon,'excited_band_floor':F(MASS+1),
      'sixth_order_square_root_remainder':remainder,'sixth_order_current_remainder':derivative,
      'original_seed_to_ground_success_lower':1-(epsilon/(1-epsilon))**2,
      'all_state_outside_radius_five_upper':outside}


@lru_cache(maxsize=4)
def exact_sparse(cap):
    m=matter_exact();basis,ann=ladders(cap,True);d=len(basis);eye=sp.eye(d)
    b=sp.kronecker_product(m['B'],eye)
    for derivative,a in zip(m['derivative'],ann):b+=sp.kronecker_product(derivative,(a+a.H)/8)
    e=sp.kronecker_product(m['external'][2],eye)
    jf=-sp.I*(ann[0].H*ann[1]-ann[1].H*ann[0])
    jz=sp.kronecker_product(m['J'][2],eye)+sp.kronecker_product(sp.eye(20),jf)
    number=[sum(m['basis'][(i//d)%10])+sum(basis[i%d])for i in range(20*d)]
    return b.todok(),e.todok(),jz.todok(),number


@lru_cache(maxsize=2048)
def scalar(expr):
    expr=sp.sympify(expr)
    if expr.is_Rational:return iv.mpf(int(expr.p))/int(expr.q)
    if expr.is_Add:return sum((scalar(a)for a in expr.args),iv.mpf(0))
    if expr.is_Mul:
        result=iv.mpf(1)
        for a in expr.args:result*=scalar(a)
        return result
    if expr.is_Pow:
        base,power=expr.args
        return scalar(base)**(iv.mpf(int(power.p))/int(power.q))
    raise TypeError(f'nonalgebraic coefficient {expr}')


def interval_entries(entries):
    return tuple((i,j,iv.mpc(scalar(sp.re(v)),scalar(sp.im(v))))for(i,j),v in entries.items())


def multiply(entries,v):
    out=[iv.mpc(0)for _ in v]
    for i,j,x in entries:out[i]+=x*v[j]
    return out


def inner(a,b):return sum((iv.mpc(x.real,-x.imag)*y for x,y in zip(a,b)),iv.mpc(0))

def norm_upper(v):
    square=sum((x.real**2+x.imag**2 for x in v),iv.mpf(0))
    return iv.sqrt(square)


def endpoints(v):return [float(np.nextafter(float(v.a),-np.inf)),float(np.nextafter(float(v.b),np.inf))]

def number(x):return iv.mpf(x.numerator)/x.denominator if isinstance(x,F)else iv.mpf(x)


def normalized_dyadic(v):
    integers=[[int(round(float(x.real)*2**40)),int(round(float(x.imag)*2**40))]for x in v]
    denominator=sum(a*a+b*b for a,b in integers)
    scale=iv.sqrt(iv.mpf(denominator))
    return [iv.mpc(a,b)/scale for a,b in integers],integers,denominator


def certify(cap):
    iv.dps=50;scalar.cache_clear()
    d=source(cap);b,e,j,diag=exact_sparse(cap);B=interval_entries(b);E=interval_entries(e);J=interval_entries(j)
    phi,integers,denominator=normalized_dyadic(d['psi'])
    powers=[phi]
    for _ in range(6):powers.append(multiply(B,multiply(B,powers[-1])))
    coefficients=[F(sp.binomial(sp.Rational(1,2),k))/MASS**(2*k-1)for k in range(1,7)]
    h=[(MASS+n)*v for n,v in zip(diag,phi)]
    for coefficient,v in zip(coefficients,powers[1:]):
        h=[u+number(coefficient)*x for u,x in zip(h,v)]
    rayleigh=inner(phi,h).real
    residual=norm_upper([x-rayleigh*v for x,v in zip(h,phi)])
    bounds=rational_bounds(cap);delta=number(bounds['sixth_order_square_root_remainder'])
    # A positive Taylor remainder bounds actual Rayleigh energy and variance.
    rup=rayleigh.b+delta;rlow=rayleigh.a;rnorm=residual.b+delta
    separation=MASS+1-rup
    if float(separation.a)<=0:raise ValueError('trial vector does not lie below the excited-band floor')
    low=rlow-rnorm*rnorm/separation;high=rup
    # Exact Jz=+1/2 projection plus the isolated energy band select one line.
    jphi=multiply(J,phi);jerror=norm_upper([x-v/2 for x,v in zip(jphi,phi)])
    bad_weight=(rnorm/separation)**2+jerror.b**2
    trace_error=2*iv.sqrt(bad_weight)
    # Polynomial Frechet derivative expectation uses all k ordered products.
    current=iv.mpc(0);Tp=[]
    for v in powers[:6]:
        be=multiply(B,multiply(E,v));eb=multiply(E,multiply(B,v));Tp.append([x+y for x,y in zip(be,eb)])
    for k,coefficient in enumerate(coefficients,1):
        for l in range(k):current+=number(coefficient)*inner(powers[l],Tp[k-1-l])
    current_error=number(bounds['sixth_order_current_remainder'])+trace_error*number(bounds['kinetic_norm_upper']*F(4,5)/MASS)
    current_interval=iv.mpf([(current.real-current_error).a,(current.real+current_error).b])
    energy_interval=iv.mpf([low.a,high.b])
    g=4*energy_interval*abs(current_interval)
    return {'cap':cap,'precision_decimal_digits':50,'trial_vector_definition':'integer real/imaginary coefficients divided by square root of their exact sum of squares',
        'trial_integer_coefficients':integers,'trial_squared_integer_norm':str(denominator),
        'rational_bounds':{k:str(v)for k,v in bounds.items()},
        'lower_polynomial_Rayleigh_interval':endpoints(rayleigh),'lower_polynomial_residual_upper':float(np.nextafter(float(residual.b),np.inf)),
        'actual_ground_energy_interval':endpoints(energy_interval),'exact_dressed_line_trace_distance_upper':float(np.nextafter(float(trace_error.b),np.inf)),
        'native_probe_derivative_interval':endpoints(current_interval),'rest_energy_normalized_static_g_interval':endpoints(g),
        'scope':'Algebraic matrix entries and normalized rational trial vectors; outward interval arithmetic, Taylor remainders and the proved band gap. These are finite-source enclosures, without a physical electron error budget.'}
