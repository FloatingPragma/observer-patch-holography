"""Source-exact matrix spectral algebra, residues and finite cyclic bounds.

The source is unchanged. Python integers admit the polynomial identities;
the earlier scalar spectrum's exact intervals admit the spectral bounds.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

import numpy as np

from model.interacting_vacuum import ground
from model.retarded_matter_response import response, spectrum


DEGREE = 78
CYCLIC_GAP = F(299, 10000)
CYCLIC_UPPER = F(37, 20)


def _add(a, b):
    return a[0]+b[0], a[1]+b[1]


def _mul(a, b):
    return a[0]*b[0]+3*a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def _scale(a, value):
    return value*a[0], value*a[1]


def _divide(a, b):
    denominator=b[0]**2-3*b[1]**2
    if not denominator:raise ValueError('zero quadratic-field denominator')
    return (F(a[0]*b[0]-3*a[1]*b[1], denominator),
            F(a[1]*b[0]-a[0]*b[1], denominator))


def _evaluate(coefficients, point):
    value=(0,0)
    for coefficient in coefficients:value=_add(_mul(value,point),coefficient)
    return value


def _convolution(a,b):
    result=[(0,0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):result[i+j]=_add(result[i+j],_mul(x,y))
    return tuple(result)


def _remainder(polynomial, monic):
    if monic[0]!=(1,0):raise ValueError('monic quadratic-field divisor required')
    result=list(polynomial)
    for i in range(len(result)-len(monic)+1):
        coefficient=result[i]
        for j,x in enumerate(monic):
            product=_mul(coefficient,x)
            result[i+j]=(result[i+j][0]-product[0],result[i+j][1]-product[1])
    return tuple(result[-len(monic)+1:])


def _apply(rows,value):
    result=np.zeros(value.shape,dtype=object)
    for i,row in enumerate(rows):
        for j,coefficient in row:result[i]+=coefficient*value[j]
    return result


def annihilates_columns(polynomial):
    data=response.dilation('ground','down')
    _,_,_,rows=spectrum._source()
    a,b=(x.astype(object) for x in data.vectors)
    x,y=np.zeros(a.shape,dtype=object),np.zeros(b.shape,dtype=object)
    for c,d in polynomial:
        x,y=_apply(rows,x),_apply(rows,y)
        x+=c*a+3*d*b
        y+=c*b+d*a
    return not np.any(x) and not np.any(y)


@lru_cache(maxsize=1)
def exact_data():
    scalar=spectrum.exact_data()
    polynomial=scalar['polynomial']
    if len(polynomial)!=DEGREE+1 or not annihilates_columns(polynomial):
        raise ArithmeticError('scalar source polynomial does not annihilate all eight columns')
    data=response.dilation('ground','down')
    _,_,_,rows=spectrum._source()
    a,b=(x.astype(object) for x in data.vectors)
    adjacency=ground.band_data()[0].astype(object)
    identity=np.eye(8,dtype=object)
    if not np.array_equal(adjacency@adjacency,3*identity) or any(adjacency.diagonal()):
        raise ArithmeticError('the original signed cube algebra changed')
    i,j=next((i,j) for i in range(8) for j in range(8) if adjacency[i,j])
    x,y=a.copy(),b.copy();moments=[]
    for order in range(DEGREE):
        ma,mb=a.T@x+3*b.T@y,a.T@y+b.T@x
        diagonal=(int(ma[0,0]),int(mb[0,0]))
        edge=(int(ma[i,j]*adjacency[i,j]),int(mb[i,j]*adjacency[i,j]))
        if (not np.array_equal(ma,diagonal[0]*identity+edge[0]*adjacency) or
                not np.array_equal(mb,diagonal[1]*identity+edge[1]*adjacency)):
            raise ArithmeticError('full source matrix moment left span(I,A)')
        if diagonal!=scalar['moments'][order]:
            raise ArithmeticError('full matrix and local scalar moments disagree')
        moments.append((diagonal,edge))
        x,y=_apply(rows,x),_apply(rows,y)
    edge_numerator=[]
    for j in range(DEGREE):
        value=(0,0)
        for k in range(j+1):value=_add(value,_mul(polynomial[k],moments[j-k][1]))
        edge_numerator.append(value)
    numerator=scalar['numerator']
    square_i=_convolution(numerator,numerator)
    square_a=_convolution(edge_numerator,edge_numerator)
    difference=[_add(i,_scale(a,-3)) for i,a in zip(square_i,square_a,strict=True)]
    if any(x or y for x,y in _remainder(difference,polynomial)):
        raise ArithmeticError('rank-four residue polynomial identity failed')
    return {'polynomial':polynomial,'identity_numerator':numerator,
            'adjacency_numerator':tuple(edge_numerator),'moments':tuple(moments),
            'derivative':scalar['derivative']}


def _field_interval(value,bits=256):
    scale=1<<bits;root=isqrt(3*scale*scale)
    a,b=value
    endpoints=(a+b*F(root,scale),a+b*F(root+1,scale))
    return min(endpoints),max(endpoints)


def _encoded(value):return [str(x) for x in value]


def static_certificate(data=None):
    data=exact_data() if data is None else data
    point=(0,-4)
    denominator=_evaluate(data['polynomial'],point)
    c0=_scale(_divide(_evaluate(data['identity_numerator'],point),denominator),-F(16,1296))
    c1=_scale(_divide(_evaluate(data['adjacency_numerator'],point),denominator),-F(16,1296))
    determinant=_add(_mul(c0,c0),_scale(_mul(c1,c1),-3))
    y0,y1=_divide(c0,determinant),_scale(_divide(c1,determinant),-1)
    if (_add(_mul(c0,y0),_scale(_mul(c1,y1),3))!=(1,0) or
            _add(_mul(c0,y1),_mul(c1,y0))!=(0,0)):
        raise ArithmeticError('exact static inverse failed in the cube algebra')
    fixed=(F('12.9058381883'),F('-4.9026412320'))
    errors=[_field_interval(_add(v,(-w,0))) for v,w in zip((c0,c1),fixed,strict=True)]
    norm_error=max(map(abs,errors[0]))+2*max(map(abs,errors[1]))
    if norm_error>=F(1,10**7):raise ArithmeticError('exact static value disagrees with prior enclosure')
    return {'coefficient_basis':'[rational, coefficient of sqrt(3)]',
            'positive_static_response':'C=-G(0)=c0 I+c1 A',
            'c0':_encoded(c0),'c1':_encoded(c1),
            'inverse_c0':_encoded(y0),'inverse_c1':_encoded(y1),
            'memory_c0':_encoded(_add(y0,(-F(1,2),0))),
            'memory_c1':_encoded(_add(y1,(-F(1,16),0))),
            'exact_static_inverse_identity':True,
            'prior_static_enclosure_verified':'1/10000000',
            'floating_linear_solver_used':False}


def memory_remainder(radius,order):
    if type(radius) not in (int,F) or not 0<=radius<CYCLIC_GAP:
        raise ValueError('exact nonnegative radius below299/10000 required')
    if type(order) is not int or not 0<=order<=32:
        raise ValueError('integer order0..32 required')
    radius=F(radius)
    return radius**(order+1)/(4*CYCLIC_GAP**(order+1)*(CYCLIC_GAP-radius))


@lru_cache(maxsize=1)
def certificate():
    data=exact_data()
    local=spectrum.certificate()
    coefficients=tuple(spectrum._coefficient(x) for x in data['polynomial'])
    edge=tuple(spectrum._coefficient(x) for x in data['adjacency_numerator'])
    derivative=tuple(spectrum._coefficient(x) for x in data['derivative'])
    brackets,bits=spectrum._brackets(coefficients)
    bands=[]
    for bracket,row in zip(brackets,local['rows'],strict=True):
        root=spectrum._refine(coefficients,bracket,bits)
        beta=spectrum._ratio(spectrum._evaluate(edge,root),spectrum._evaluate(derivative,root))
        sign=1 if beta[0]>0 else -1 if beta[1]<0 else 0
        if not sign:raise ArithmeticError('rank-four residue band sign was not certified')
        bands.append({'line':row['line'],'band':sign,
                      'gap_above_ground_interval':row['gap_above_ground_interval'],
                      'local_weight_interval':row['normalized_weight_interval'],
                      'residue_formula':'2 alpha P_plus' if sign==1 else '2 alpha P_minus',
                      'rank':4})
    lower=F(local['rows'][0]['gap_above_ground_interval'][0])
    upper=F(local['rows'][-1]['gap_above_ground_interval'][1])
    if not CYCLIC_GAP<lower<=upper<CYCLIC_UPPER:
        raise ArithmeticError('chosen rational cyclic bounds do not contain the certified spectrum')
    plus=sum(row['band']==1 for row in bands)
    return {'source':'original grounded cube, c=1/8,U=1; full eight-site down-addition response',
            'all_eight_columns_annihilated_exactly':True,
            'moment_matrices_in_span_I_A':DEGREE,
            'adjacency_identity':'A^2=3I, tr(A)=0; P_plus/minus=(I +/- A/sqrt3)/2 have rank4',
            'polynomial_descending_Qsqrt3':[list(x) for x in data['polynomial']],
            'identity_numerator_descending_Qsqrt3':[list(x) for x in data['identity_numerator']],
            'adjacency_numerator_descending_Qsqrt3':[list(x) for x in data['adjacency_numerator']],
            'full_resolvent':'G(z)=16[p_I(x)I+p_A(x)A]/[1296 q(x)], x=16z-4sqrt3',
            'rank_polynomial_identity':'q divides p_I^2-3p_A^2 exactly over Q(sqrt3)',
            'matrix_residue_rows':bands,'band_line_counts':{'plus':plus,'minus':DEGREE-plus},
            'residue_ranks_sum':4*DEGREE,'observed_cyclic_dimension':4*DEGREE,
            'full_source_sector_dimension':560,'orthogonal_unobserved_dimension':560-4*DEGREE,
            'all_sites_have_same_local_spectral_measure':True,
            'band_response':'G=g_plus P_plus+g_minus P_minus; g_sign=sum_in_band 2alpha/(z-E)',
            'band_weight_sum':'1 for each band, from M0=I and rank-polynomial identity',
            'exact_inverse':'G^-1=[F I-J A]/[F^2-3J^2] when G=FI+JA',
            'exact_finite_inverse_and_memory_spatial_range':'onsite plus existing signed nearest-neighbor adjacency',
            'strict_observed_cyclic_generator_interval':[str(CYCLIC_GAP),str(CYCLIC_UPPER)],
            'same_interval_for_cyclic_compressed_generator':True,
            'larger_memory_disk':'|z|<299/10000 on the observed cyclic space; the orthogonal dark space is omitted only after exact reducing-space proof',
            'memory_remainder_formula':'r^(p+1)/[4g^(p+1)(g-r)],g=299/10000',
            'examples':[{'radius':str(r),'order':p,'memory_remainder':str(memory_remainder(r,p))}
                        for r,p in ((F(1,100),4),(F(1,100),8),(F(1,50),8))],
            'static_exact':static_certificate(data),
            'all_560_eigenvalues_or_degeneracies_resolved':False,
            'gap_bound_applied_to_unobserved_sector':False,
            'infinite_volume_locality_or_gap_inferred_from_cube':False,
            'common_Lorentz_metric_or_physical_masses_identified':False,
            'native_measurement_executed':False}
