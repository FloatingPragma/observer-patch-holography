"""Arbitrary-precision positive amplitude-tree angles with rational witnesses."""
from fractions import Fraction as F
from math import isqrt
from model.native_spin_observer import precision


def sqrt_interval(value,bits):
    if type(value) is not F or value<0 or type(bits) is not int or bits<1:
        raise ValueError('nonnegative rational and positive bit precision required')
    scale=1<<bits;n=isqrt(value.numerator*scale*scale//value.denominator)
    lo=F(n,scale);return lo,lo if lo*lo==value else F(n+1,scale)


def divide(a,b):
    if b[0]<=0:raise ValueError('positive denominator interval required')
    return precision.mul(a,(1/b[1],1/b[0]))


def arctan_interval(x,bits):
    if type(x) is not F or not 0<=x<=F(1,2):raise ValueError('reduced rational arctangent in [0,1/2] required')
    total=F();power=x;n=0
    while True:
        total+=(-1)**n*power/F(2*n+1);n+=1;power*=x*x
        following=(-1)**n*power/F(2*n+1)
        if abs(following)<=F(1,1<<bits):return min(total,total+following),max(total,total+following)


def node(a,b,bits):
    """Return a dyadic approximation to -atan2(sqrt(b),sqrt(a))."""
    if type(a) is not int or type(b) is not int or min(a,b)<0 or a+b==0:
        raise ValueError('nonzero nonnegative integer subtree squared norms required')
    if type(bits) is not int or bits<1:raise ValueError('positive angle precision required')
    working=bits+12
    while True:
        # Scale removes dependence on the potentially huge PF integer amplitudes.
        sa=sqrt_interval(F(a,a+b),working);sb=sqrt_interval(F(b,a+b),working)
        x=divide(sb,(1+sa[0],1+sa[1]))
        root_lo=sqrt_interval(1+x[0]*x[0],working)[0]
        root_hi=sqrt_interval(1+x[1]*x[1],working)[1]
        t=divide(x,(1+root_lo,1+root_hi))
        if t[1]>F(1,2):working*=2;continue
        lo=-4*arctan_interval(t[1],working)[1]
        hi=-4*arctan_interval(t[0],working)[0]
        # Keep retained enclosures short without losing the rational proof.
        grid=1<<working
        lo=F((lo*grid).numerator//(lo*grid).denominator,grid)
        hi=-F((-hi*grid).numerator//(-hi*grid).denominator,grid)
        scale=1<<bits;mid=(lo+hi)/2;integer=(mid*scale+F(1,2)).numerator//(mid*scale+F(1,2)).denominator
        value=F(integer,scale);error=max(abs(value-lo),abs(value-hi))
        if error<=F(1,1<<bits):return {'angle':value,'lower':lo,'upper':hi,'error':error,'working_bits':working}
        working*=2


def levels(values,bits):
    from model.centered_native_preparation.tree import walsh
    if type(values) is not tuple or not values or len(values)&(len(values)-1) or any(type(v) is not int or v<=0 for v in values):
        raise ValueError('complete strictly positive integer amplitudes required')
    n=len(values).bit_length()-1;prefix=[0]
    for v in values:prefix.append(prefix[-1]+v*v)
    out=[]
    for bit in reversed(range(n)):
        nodes=[]
        for label in range(1<<(n-bit-1)):
            start=label*(1<<(bit+1));half=1<<bit
            nodes.append(node(prefix[start+half]-prefix[start],prefix[start+2*half]-prefix[start+half],bits))
        coefficients=tuple(v/len(nodes) for v in walsh([x['angle'] for x in nodes]))
        out.append({'target':bit,'nodes':nodes,'coefficients':coefficients})
    return out
