"""Source-generated dressed CAR response and its complete residual.

Each canonical family acts on the same full occupation source. Their mixed
CAR are not zero operators. Orthogonality belongs to the invariant-ground
retarded response space, not to additional physical modes or species.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import lcm

import numpy as np
import sympy as sy

from model.full_fock_field.field import FullFock
from model.strong_coupling_refinement import algebra


def _mul(a,b):
    a=np.asarray(a,dtype=np.int64);b=np.asarray(b,dtype=np.int64)
    bound=a.shape[1]*max(abs(int(x))for x in a.flat)*max(abs(int(x))for x in b.flat)
    if bound>=2**62:raise OverflowError('exact integer product exceeds its arithmetic domain')
    return a@b


def _comm(a,b):return _mul(a,b)-_mul(b,a)
def _anti(a,b):return _mul(a,b)+_mul(b,a)


def source(sides,c,u):
    """Literal full-number-block signed source, without ground projection."""
    sites,edges=algebra.geometry(sides)
    if sites>3:raise ValueError('complete exact oracle admits two or three sites')
    c,u=F(c),F(u)
    if c<=0 or u<=0:raise ValueError('positive declared hopping and interaction required')
    scale=lcm((c/2).denominator,(u/2).denominator)
    size=1<<(2*sites);h=np.zeros((size,size),dtype=np.int64)
    for number in range(2*sites+1):
        masks,t2,d=algebra.operators(sites,edges,number)
        block=int(scale*c/2)*t2.toarray()+np.diag(int(scale*u)*d-int(scale*u/2)*number)
        h[np.ix_(masks,masks)]=block
    space=FullFock(2*sites)
    a=tuple(space.annihilator(i).toarray().real.astype(np.int64)for i in range(2*sites))
    identity=np.eye(size,dtype=np.int64)
    number=tuple(_mul(x.T,x)for x in a)
    b=tuple(_mul(2*number[i^1]-identity,x)for i,x in enumerate(a))
    phase=np.diag([(-1)**sum((mask>>(2*x))&3==3 for x in range(sites))for mask in range(size)])
    one=np.zeros((2*sites,2*sites),dtype=np.int64)
    for x,y,sign in edges:
        for spin in(0,1):one[2*x+spin,2*y+spin]=one[2*y+spin,2*x+spin]=int(scale*c/2)*sign
    for i,x in enumerate(b):
        if np.any(x+_mul(_mul(phase,a[i]),phase)):
            raise ArithmeticError('actual dressed CAR differs from the onsite source phase')
        for j,y in enumerate(b):
            if np.any(_anti(x,y))or np.any(_anti(x,y.T)-(identity if i==j else 0)):
                raise ArithmeticError('dressed operators fail exact CAR')
        expected=sum((one[i,j]*a[j]for j in range(len(a))),np.zeros_like(h))+int(scale*u/2)*b[i]
        if np.any(_comm(a[i],h)-expected):
            raise ArithmeticError('actual first source response equation failed')
    if not np.any(_anti(a[0],b[0].T)):
        raise ArithmeticError('the two families were incorrectly made independent')
    return {'sites':sites,'points':tuple(product(*(range(v)for v in sides))),
            'edges':edges,'scale':scale,'H':h,'a':a,'b':b,
            'number':number,'phase':phase,'one':one,'c':c,'u':u}


def psd_rank(matrix):
    """Exact rational elimination with positive pivots and zero-row checks."""
    matrix=np.asarray(matrix)
    if (matrix.ndim!=2 or matrix.shape[0]!=matrix.shape[1]
            or matrix.dtype.kind not in 'iu' or not np.array_equal(matrix,matrix.T)):
        raise ValueError('an exact symmetric integer matrix required')
    a=[[F(int(x))for x in row]for row in matrix];size=len(a);rank=0
    for k in range(size):
        pivot=a[k][k]
        if pivot<0:raise ArithmeticError('claimed ground has a negative exact pivot')
        if not pivot:
            if any(a[k][j]or a[j][k]for j in range(k,size)):
                raise ArithmeticError('zero pivot has a nonzero row')
            continue
        rank+=1
        for i in range(k+1,size):
            if not a[i][k]:continue
            factor=a[i][k]/pivot
            for j in range(i,size):
                a[i][j]-=factor*a[k][j];a[j][i]=a[i][j]
    return rank


def ground(data,energy):
    energy=F(energy);scale=data['scale'];h=data['H']
    if (energy*scale).denominator!=1:raise ValueError('rational ground must fit source integer scale')
    shifted=h-int(energy*scale)*np.eye(len(h),dtype=np.int64)
    rank=psd_rank(shifted)
    kernel=sy.Matrix(shifted).nullspace()
    if not kernel or len(kernel)!=len(h)-rank:
        raise ArithmeticError('complete exact ground kernel and PSD rank disagree')
    v=sy.Matrix.hstack(*kernel);projector=v*(v.T*v).inv()*v.T
    rho=projector/len(kernel);den=int(sy.ilcm(*[x.q for x in rho]))
    num=np.array(rho*den,dtype=np.int64)
    if np.trace(num)!=den or np.any(_mul(h,num)-int(energy*scale)*num):
        raise ArithmeticError('normalized full source ground is inconsistent')
    return num,den,rank


def _trace(rho,operator):
    numerator,denominator=rho
    return F(sum(int(x)*int(y)for x,y in zip(numerator.T.flat,operator.flat,strict=True)),denominator)


def _square(matrix):
    return tuple(tuple(sum((matrix[i][k]*matrix[k][j]for k in range(len(matrix))),F())
                       for j in range(len(matrix)))for i in range(len(matrix)))


def bond_observable(data,x,y):
    """Integer numerator 4(Sx.Sy+eta_x.eta_y), with actual bipartite signs."""
    a,n=data['a'],data['number'];identity=np.eye(len(data['H']),dtype=np.int64)
    points=data['points']
    spin_plus={z:_mul(a[2*z].T,a[2*z+1])for z in(x,y)}
    spin_z={z:n[2*z]-n[2*z+1]for z in(x,y)}
    eta_plus={z:(-1)**sum(points[z])*_mul(a[2*z].T,a[2*z+1].T)for z in(x,y)}
    eta_z={z:n[2*z]+n[2*z+1]-identity for z in(x,y)}
    result=np.zeros_like(identity)
    for plus,z in((spin_plus,spin_z),(eta_plus,eta_z)):
        result+=_mul(z[x],z[y])+2*(_mul(plus[x],plus[y].T)+_mul(plus[x].T,plus[y]))
    return result


@lru_cache(maxsize=4)
def finite_case(sites):
    if sites==2:c,u,e0=F(3,2),F(4),F(-9,2)
    elif sites==3:c,u,e0=F(2),F(3,7),F(-45,14)
    else:raise ValueError('actual rational dimer or three-site control required')
    data=source((sites,1,1),c,u);num,den,rank=ground(data,e0);rho=(num,den)
    h,scale,a,b=data['H'],data['scale'],data['a'],data['b'];modes=len(a)
    alpha=u/2;identity=np.eye(len(h),dtype=np.int64)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if _trace(rho,_anti(x,y.T)):
                raise ArithmeticError('the actual invariant ground embeddings are not orthogonal')
    comm=tuple(_comm(h,x.T)for x in b)
    C=tuple(tuple(_trace(rho,_anti(x,y))/scale for y in comm)for x in b)
    G=tuple(tuple(_trace(rho,_anti(x.T,y))/(scale*scale)for y in comm)for x in comm)
    C2=_square(C)
    residual=tuple(tuple(G[i][j]-C2[i][j]-(alpha*alpha if i==j else 0)
                         for j in range(modes))for i in range(modes))
    one=tuple(tuple(F(int(x),scale)for x in row)for row in data['one'])
    corr={}
    for x,y,_ in data['edges']:
        corr[x,y]=_trace(rho,bond_observable(data,x,y))
    for i in range(modes):
        for j in range(modes):
            expected=one[i][j]*corr.get(tuple(sorted((i//2,j//2))),F())
            if C[i][j]!=expected:
                raise ArithmeticError('source composite block differs from its spin/eta bond readout')
    # Exact full-source spin/eta symmetries, not presumed symmetry of a fit.
    spin_plus=sum((_mul(a[2*x].T,a[2*x+1])for x in range(sites)),np.zeros_like(h))
    eta_plus=sum(((-1)**x*_mul(a[2*x].T,a[2*x+1].T)for x in range(sites)),np.zeros_like(h))
    for symmetry in(spin_plus,eta_plus):
        if np.any(_comm(h,symmetry))or np.any(_comm(num,symmetry)):
            raise ArithmeticError('full source or its ground projector lost the required symmetry')
    maximum=max(sum(abs(x)for x in row)for row in residual)
    trace=sum(residual[i][i]for i in range(modes))
    if sites==2:
        if trace or any(C[i][j]!=-3*one[i][j]for i in range(modes)for j in range(modes)):
            raise ArithmeticError('actual dimer composite response failed exact closure')
    elif trace!=F(542232,13225) or maximum!=F(173508,13225):
        raise ArithmeticError('actual larger-source residual changed')
    return {'sites':sites,'whole_Fock_dimension':len(h),'c':c,'U':u,'ground_energy':e0,
            'positive_shift_exact_rank':rank,'ground_degeneracy':len(h)-rank,
            'ground_density_integer_denominator':den,'dressed_CAR_exact':True,
            'mixed_CAR_operators_are_independent':False,'ground_response_embeddings_orthogonal':True,
            'source_spin_eta_symmetries_exact':True,'C':C,'residual_Gram':residual,
            'residual_squared_norm_upper':maximum,'residual_squared_trace':trace,
            'four_spin_eta_bond_expectations':{f'{x},{y}':value for(x,y),value in corr.items()},
            'physical_mass':None,'source_particle_species_added':0}


def approximation_budget(residual_square,interaction,time=F(0),imaginary_part=F(1)):
    residual_square,interaction,time,imaginary_part=map(F,(residual_square,interaction,time,imaginary_part))
    if residual_square<0 or interaction<0 or imaginary_part<=0:
        raise ValueError('nonnegative residual/interaction and positive resolvent height required')
    alpha=interaction/2
    # Squared form avoids rounding a square root in the exact receipt.
    squared=residual_square*min(abs(time),alpha*time*time/2)**2
    return {'original_response_vector_error_squared_upper':squared,
            'compressed_resolvent_norm_error_upper':residual_square/imaginary_part**3,
            'lorentzian_spectral_density_error_times_pi_upper':residual_square/imaginary_part**3,
            'source_success_probability_division':False}


def native_twist():
    """Literal onsite double-occupancy phase, including every export outcome."""
    from model.interacting_field import source as native
    from model.field_thermal.thermal import execute_with_exports
    from model.field_symmetry.readout import pulse_budget
    _,info=native.contact(np.pi,1.)
    word=info['word'];actual,events=execute_with_exports(word,2,np.eye(4,dtype=complex))
    if len(events)!=len(word):raise ArithmeticError('onsite phase omitted source failure events')
    complements=[]
    for event,outputs in zip(word,events,strict=True):
        dimension=1<<len(event['targets'])
        if len(outputs)!=6-dimension or any(x.shape!=(4//dimension,4)for x in outputs):
            raise ArithmeticError('onsite phase omitted a source complement')
        complements.extend(outputs)
    target=np.diag([1,1,1,-1]);joined=np.concatenate([actual-target,*complements])
    error=float(np.linalg.norm(joined,2));budget,_=pulse_budget(word)
    if error>float(2*budget):raise ArithmeticError('native onsite twist exceeds complete dilation budget')
    counts={key:value for key,value in info['resources'].items()if type(value)in(int,bool)}
    return {'native_pulses':len(word),'full_input_columns':4,'failure_slots':len(complements),
            'complete_dilation_operator_error_interval':['0',str(2*budget)],
            'numerical_admission_passed':True,
            'resources':counts,'floating_diagnostic_is_outward_certificate':False,
            'odd_transfer_or_reservoir_service_derived':False}


def certificate():
    rows=[finite_case(2),finite_case(3)]
    return {'canonical_operator_families':2,'independent_CAR_flavours_added':0,
            'physical_species_added':0,'all_volume_first_composite_equation':True,
            'all_volume_residual_discarded':False,'finite_full_ground_controls':rows,
            'native_onsite_twist':native_twist(),
            'three_site_response_error':approximation_budget(rows[1]['residual_squared_norm_upper'],
                                                              rows[1]['U'],F(1,8),F(1)),
            'dimer_or_aligned_reference_transferred_to_full_continuum':False}
