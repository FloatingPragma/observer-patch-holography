"""Rational upper Darboux bound for the 3D KLS infrared integral."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import factorial,isqrt


def ceil_fraction(x):return -((-x.numerator)//x.denominator)


def atan_interval(inverse,terms):
    q=F(1,inverse)
    a=sum(((-1)**i*q**(2*i+1)/F(2*i+1) for i in range(terms)),F(0))
    b=a+(-1)**terms*q**(2*terms+1)/F(2*terms+1)
    return min(a,b),max(a,b)


def pi_interval():
    # Machin's identity; alternating arctangent remainders have fixed signs.
    a,b=atan_interval(5,32);c,d=atan_interval(239,10)
    return 16*a-4*d,16*b-4*c


def cosine_upper(x,scale):
    # Taylor through degree40 with Lagrange remainder <=x^42/42!.
    value=sum(((-1)**i*x**(2*i)/F(factorial(2*i)) for i in range(21)),F(0))
    value+=x**42/F(factorial(42))
    return min(scale,ceil_fraction(value*scale))


def certificate(divisions=64):
    if type(divisions) is not int or not 8<=divisions<=256:
        raise ValueError('divisions must be an integer in [8,256]')
    n=divisions; scale=1<<40
    plo,phi=pi_interval()
    assert F(3)<plo<phi<F(22,7)
    values=[cosine_upper(plo*i/n,scale) for i in range(n)]
    total=0; cells=0
    for i,j,k in combinations_with_replacement(range(n),3):
        if i==j==k==0:continue
        multiplicity=1 if i==k else (3 if i==j or j==k else 6)
        cells+=multiplicity
        c=values[i]+values[j]+values[k]
        if c<=0:continue
        if c>=3*scale:raise ArithmeticError('unisolated infrared singularity')
        # F(C)=(C/3)*sqrt((3+C)/(3-C)); its derivative is positive for C>0.
        # ceil(scale*F(C)) is found by one exact integer square root.
        numerator=c*c*(3*scale+c)
        denominator=9*(3*scale-c)
        root=isqrt(numerator//denominator)
        if root*root*denominator<numerator:root+=1
        total+=multiplicity*root
    assert cells==n**3-1
    regular=F(total,scale*n**3)
    # On origin cube, F<=pi*sqrt3/|x|. Enclosing it in its radius-sqrt3*h
    # ball gives normalized integral <=6*pi*sqrt3/n² <33/n².
    assert 6*F(22,7)*F(7,4)<=33 and F(7,4)**2>3
    singular=F(33,n*n)
    upper=regular+singular
    return {'divisions':n,'regular_upper':str(regular),'origin_cell_upper':str(singular),
            'integral_upper':str(upper),'integral_below_two_fifths':upper<F(2,5),
            'cosine_grid_scale':scale,'nonorigin_cells':cells,
            'pi_lower':str(plo),'pi_upper':str(phi)}
