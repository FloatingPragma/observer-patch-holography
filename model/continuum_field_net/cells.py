"""Exact one-dimensional cell-average witnesses for the tensor-product map."""
from fractions import Fraction as F


def exact(value):
    if type(value) not in (int,F):raise ValueError('exact rational required')
    return F(value)


def polynomial_integral(coefficients,left,right):
    left,right=exact(left),exact(right)
    return sum((exact(c)*(right**(j+1)-left**(j+1))/(j+1)
                for j,c in enumerate(coefficients)),F(0))


def product(left,right):
    result=[F(0)]*(len(left)+len(right)-1)
    for i,a in enumerate(left):
        for j,b in enumerate(right):result[i+j]+=exact(a)*exact(b)
    return result


def derivative(coefficients):
    return [j*exact(c) for j,c in enumerate(coefficients)][1:] or [F(0)]


def cell_integrals(coefficients,left,right,spacing):
    """Integrals of a polynomial restricted to [left,right], not normalized samples."""
    left,right,spacing=exact(left),exact(right),exact(spacing)
    if left>=right or spacing<=0 or not coefficients:
        raise ValueError('positive spacing and nonempty bounded interval required')
    first=int(left//spacing)-1;last=int(right//spacing)+1
    result={}
    for n in range(first,last+1):
        lo=max(left,(F(n)-F(1,2))*spacing)
        hi=min(right,(F(n)+F(1,2))*spacing)
        if lo<hi:
            value=polynomial_integral(coefficients,lo,hi)
            if value:result[n]=value
    return result


def pairing(left,right,spacing,alternating=False):
    spacing=exact(spacing)
    if spacing<=0:raise ValueError('positive spacing required')
    return sum(((-1 if alternating and n%2 else 1)*x*right.get(n,F(0))/spacing
                for n,x in left.items()),F(0))


def certificate():
    # Continuous compact H1 function f(x)=1-x^2 on [-1,1], zero elsewhere.
    f=(1,0,-1);norm=polynomial_integral(product(f,f),-1,1)
    grad=derivative(f);gradient=polynomial_integral(product(grad,grad),-1,1)
    rows=[]
    for a in (F(1),F(1,2),F(1,4),F(1,8),F(1,16)):
        cells=cell_integrals(f,-1,1,a)
        projected=pairing(cells,cells,a)
        cross=pairing(cells,cells,a,True)
        if not 0<=norm-projected<=a*a*gradient:
            raise ArithmeticError('cell contraction or projection error failed')
        if cross*cross>a*a*norm*gradient:
            raise ArithmeticError('alternating valley overlap estimate failed')
        rows.append({'spacing':str(a),'norm_squared':str(projected),
                     'projection_error_squared':str(norm-projected),
                     'cross_valley_pairing':str(cross)})
    return {'dimension_of_finite_witness':1,'continuum_norm_squared':str(norm),
            'continuum_derivative_norm_squared':str(gradient),'refinements':rows,
            'three_dimensional_argument':'commuting coordinate conditional expectations and telescoping',
            'samples_are_normalized_cell_integrals':True,
            'finite_witness_is_continuum_state':False}
