"""Magnetic Bloch moments of the actual centered signed-CAR source.

The local polynomials use exact Gaussian rationals in the Majorana algebra.
Their ground expectations remain source observables, not fitted parameters.
"""
from fractions import Fraction as F
from itertools import permutations, product

import numpy as np

from model.field_symmetry.symmetry import link_sign
from model.source_family_compatibility.hole_dispersion import gamma_matrices

_ZERO=(F(0),F(0)); _ONE=(F(1),F(0)); _I=(F(0),F(1))


def _q(value):
    return value if type(value)is tuple else (F(value),F(0))


def _qadd(a,b):return (a[0]+b[0],a[1]+b[1])
def _qmul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def _qconj(a):return(a[0],-a[1])


def add(*operators):
    out={}
    for op in operators:
        for mask,value in op.items():
            out[mask]=_qadd(out.get(mask,_ZERO),value)
            if out[mask]==_ZERO:del out[mask]
    return out


def scale(op,value):
    value=_q(value)
    return {m:z for m,c in op.items()if(z:=_qmul(c,value))!=_ZERO}


def _product_sign(a,b):
    parity=0
    while b:
        bit=b&-b;index=bit.bit_length()-1
        parity^=(a>>(index+1)).bit_count()&1;b^=bit
    return -1 if parity else 1


def multiply(a,b):
    out={}
    for x,c in a.items():
        for y,d in b.items():
            key=x^y;value=_qmul(_qmul(c,d),_q(_product_sign(x,y)))
            out[key]=_qadd(out.get(key,_ZERO),value)
            if out[key]==_ZERO:del out[key]
    return out


def dagger(op):
    return {m:_qmul(_qconj(c),_q((-1)**(m.bit_count()*(m.bit_count()-1)//2)))
            for m,c in op.items()}


def commutator(a,b):return add(multiply(a,b),scale(multiply(b,a),-1))
def anticommutator(a,b):return add(multiply(a,b),multiply(b,a))


def annihilator(mode):
    if type(mode)is not int or mode<0:raise ValueError('nonnegative integer mode required')
    return {1<<(2*mode):(F(1,2),F(0)),1<<(2*mode+1):(F(0),F(1,2))}


def number(mode):
    a=annihilator(mode);return multiply(dagger(a),a)


def dressed(mode):
    return multiply(add(scale(number(mode^1),2),{0:(F(-1),F(0))}),annihilator(mode))


def bond_derivative(x,y,spin=0):
    """D_xy,s with [b_xs,T]=sum_y h_xy D_xy,s, on the literal CAR."""
    if x==y or spin not in(0,1):raise ValueError('distinct sites and spin zero or one required')
    a=annihilator(2*x+spin);other=annihilator(2*x+(spin^1))
    neighbor=annihilator(2*y+spin);other_neighbor=annihilator(2*y+(spin^1))
    z=add(scale(number(2*x+(spin^1)),2),{0:(F(-1),F(0))})
    current=add(multiply(dagger(other),other_neighbor),scale(multiply(dagger(other_neighbor),other),-1))
    return add(multiply(z,neighbor),scale(multiply(current,a),2))


def hopping(edges):
    """Exact declared real hopping coefficients, including both flavours."""
    out={}
    for x,y,t in edges:
        for spin in(0,1):
            term=multiply(dagger(annihilator(2*x+spin)),annihilator(2*y+spin))
            out=add(out,scale(add(term,dagger(term)),F(t)))
    return out


def hermitian_terms(op):
    """Expansion in Hermitian norm-one Majorana words, a complete Pauli readout."""
    if dagger(op)!=op:raise ValueError('Hermitian local source observable required')
    rows=[]
    for mask,value in sorted(op.items()):
        indices=tuple(j for j in range(mask.bit_length())if mask>>j&1)
        phase=(_ONE,_I,(-F(1),F(0)),(F(0),-F(1)))[(len(indices)*(len(indices)-1)//2)%4]
        coefficient=_qmul(value,_qconj(phase))
        if coefficient[1]:raise ArithmeticError('Hermitian word has nonreal coefficient')
        rows.append({'majorana_indices':indices,'coefficient':coefficient[0]})
    return tuple(rows)


def score_plan(op):
    rows=hermitian_terms(op);bound=sum((abs(r['coefficient'])for r in rows),F())
    return {'terms':rows,'score_absolute_bound':bound,
            'setting_probabilities':tuple(abs(r['coefficient'])/bound for r in rows)if bound else (),
            'all_plus_minus_outcomes_retained':True,'zero_observable_identity_record':not bool(bound)}


def local_observables():
    """Dimensionless g0/c² and g2/c² as exact seven-/three-site observables."""
    star=add(*(bond_derivative(0,y)for y in range(1,7)))
    g0=scale(anticommutator(star,dagger(star)),F(1,4))
    cross=anticommutator(bond_derivative(0,1),dagger(bond_derivative(2,1)))
    g2=scale(add(cross,dagger(cross)),F(1,8))
    if dagger(g0)!=g0 or dagger(g2)!=g2:raise ArithmeticError('source coefficients are not observables')
    return {'g0_over_c_squared':g0,'g2_over_c_squared':g2}


def _bilinears(up,down,epsilon,norm_square=1):
    norm_square=F(norm_square)
    nu=multiply(dagger(up),up);nd=multiply(dagger(down),down)
    return ((scale(multiply(dagger(up),down),1/norm_square),
             scale(add(nu,scale(nd,-1)),1/(2*norm_square))),
            (scale(multiply(dagger(up),dagger(down)),F(epsilon)/norm_square),
             scale(add(scale(add(nu,nd),1/norm_square),{0:(F(-1),F(0))}),F(1,2))))


def _invariant_dot(left,right):
    out={}
    for (p,z),(q,w) in zip(left,right,strict=True):
        out=add(out,multiply(z,w),
                scale(add(multiply(p,dagger(q)),multiply(dagger(p),q)),F(1,2)))
    return out


def spin_eta_dot(x,y,epsilon_x=1,epsilon_y=1):
    return _invariant_dot(_bilinears(annihilator(2*x),annihilator(2*x+1),epsilon_x),
                          _bilinears(annihilator(2*y),annihilator(2*y+1),epsilon_y))


def invariant_observables():
    """Actual q1,q2,qstar polynomials; the neighbour orbital is not an extra mode."""
    q1=spin_eta_dot(0,1,1,-1);q2=spin_eta_dot(0,2,1,1)
    up=add(*(annihilator(2*y)for y in range(1,7)))
    down=add(*(annihilator(2*y+1)for y in range(1,7)))
    star=_invariant_dot(_bilinears(annihilator(0),annihilator(1),1),
                        _bilinears(up,down,-1,6))
    for value in(q1,q2,star):
        if dagger(value)!=value:raise ArithmeticError('invariant correlator is not Hermitian')
    # Spin averaging is an operator identity, before taking a ground expectation.
    force=[];cross=[]
    for spin in(0,1):
        y=add(*(bond_derivative(0,z,spin)for z in range(1,7)))
        force.append(scale(anticommutator(y,dagger(y)),F(1,4)))
        d=anticommutator(bond_derivative(0,1,spin),dagger(bond_derivative(2,1,spin)))
        cross.append(scale(add(d,dagger(d)),F(1,8)))
    if scale(add(*force),F(1,2))!=add({0:(F(9,2),F(0))},scale(star,-12)):
        raise ArithmeticError('source star coefficient does not equal its invariant orbital readout')
    if scale(add(*cross),F(1,2))!=q2:
        raise ArithmeticError('source axial coefficient does not equal the actual next-neighbour correlator')
    return {'q1':q1,'q2':q2,'qstar':star}


def signed_link(x,axis):
    """Periodic parity extension of the actual nonnegative-coordinate source API."""
    return link_sign(tuple(v%2 for v in x),axis)


def magnetic_phase(x,axis):return (-1)**sum(x[:axis])


def permutation_sign(perm):return (-1)**sum(perm[i]>perm[j]for i in range(3)for j in range(i+1,3))


def cubic_rotations():
    return tuple((perm,signs)for perm in permutations(range(3))for signs in product((-1,1),repeat=3)
                 if permutation_sign(perm)*signs[0]*signs[1]*signs[2]==1)


def rotate(x,perm,signs):
    out=[0,0,0]
    for j in range(3):out[perm[j]]=signs[j]*x[j]
    return tuple(out)


def rotation_phase(x,perm):
    return (-1)**sum(x[i]*x[j]for i in range(3)for j in range(i+1,3)if perm[i]>perm[j])


def source_symmetry_audit():
    """Check the literal link/gauge identities, including diagonal annihilation."""
    translations=rotations=normalizations=0
    for x in product(range(-2,3),repeat=3):
        for axis in range(3):
            y=list(x);y[axis]+=1;y=tuple(y)
            for j in range(3):
                xx=list(x);xx[j]+=1;yy=list(y);yy[j]+=1
                if signed_link(tuple(xx),axis)!=magnetic_phase(tuple(xx),j)*magnetic_phase(tuple(yy),j)*signed_link(x,axis):
                    raise ArithmeticError('source magnetic translation does not preserve its links')
                translations+=1
            for perm,signs in cubic_rotations():
                rx=rotate(x,perm,signs);ry=rotate(y,perm,signs);a=perm[axis]
                start=rx if signs[axis]>0 else ry
                if signed_link(start,a)!=rotation_phase(x,perm)*rotation_phase(y,perm)*signed_link(x,axis):
                    raise ArithmeticError('proper cubic lift does not preserve the actual source')
                rotations+=1
                previous=list(x);previous[axis]-=1
                left=rotation_phase(x,perm)*magnetic_phase(x,axis)
                right=magnetic_phase(rx,perm[axis])*rotation_phase(tuple(previous),perm)
                if left!=right:
                    raise ArithmeticError('cubic lift fails to normalize source magnetic translations')
                normalizations+=1
    diagonals=[]
    for delta in product((-1,0,1),repeat=3):
        if sum(abs(x)for x in delta)!=2:continue
        fixes=[(p,s)for p,s in cubic_rotations()if rotate(delta,p,s)==delta and rotation_phase(delta,p)==-1]
        if not fixes:raise ArithmeticError('mixed-axis coefficient lacks its cancelling source symmetry')
        diagonals.append({'displacement':delta,'permutation':fixes[0][0],'axis_signs':fixes[0][1]})
    return {'proper_rotations':len(cubic_rotations()),'exact_link_translation_checks':translations,
            'exact_link_rotation_checks':rotations,'exact_rotation_translation_normalizations':normalizations,
            'mixed_axis_cancellations':diagonals}


def endpoint_bounds(kappa,g0,g2,coupling=F(1)):
    """Exact endpoints of the actual range-two symbols after invariant averaging."""
    kappa,g0,g2,c=map(F,(kappa,g0,g2,coupling))
    if c<=0 or abs(kappa)>3:raise ValueError('positive source hopping and bond coefficient in [-3,3] required')
    a=g0-6*g2;b=g0+6*g2-3*c*c*kappa*kappa
    return {'G_minimum':g0-6*abs(g2),'G_maximum':g0+6*abs(g2),
            'residual_endpoint_S0':a,'residual_endpoint_S3':b,
            'residual_minimum':min(a,b),'residual_maximum':max(a,b),
            'full_momentum_positivity':a>=0 and b>=0}


def symbols(k,kappa,g0,g2,coupling=1.):
    k=np.asarray(k,dtype=float);values=np.asarray([kappa,g0,g2,coupling],dtype=float)
    if k.shape!=(3,)or not np.isfinite(k).all()or not np.isfinite(values).all()or coupling<=0:
        raise ValueError('finite three-momentum and source coefficients required')
    h=float(coupling)*sum(np.cos(k[j])*g for j,g in enumerate(gamma_matrices()))
    g=float(g0)+2*float(g2)*sum(np.cos(2*k));C=float(kappa)*h
    return {'h':h,'C':C,'Y_Gram':g*np.eye(8),'residual_Gram':g*np.eye(8)-C@C}


def uniform_error(kappa_error,g0_error,g2_error,coupling=F(1)):
    dk,d0,d2,c=map(F,(kappa_error,g0_error,g2_error,coupling))
    if min(dk,d0,d2)<0 or c<=0:raise ValueError('nonnegative errors and positive hopping required')
    return {'G_operator_error':d0+6*d2,'residual_operator_error':d0+6*d2+18*c*c*dk,
            'both_bond_estimates_clipped_to_minus3_plus3':True}


def exact_algebra_audit():
    identity={0:_ONE};checks=0
    for x,y in((0,1),(1,2)):
        t=hopping(((x,y,F(1)),))
        for spin in(0,1):
            if commutator(dressed(2*x+spin),t)!=bond_derivative(x,y,spin):
                raise ArithmeticError('local derivative does not equal the actual hopping commutator')
            checks+=1
    # A disjoint odd pair must anticommute, even with Jordan-Wigner signs.
    if anticommutator(bond_derivative(0,1),dagger(bond_derivative(2,3))):
        raise ArithmeticError('odd source locality failed')
    for mode in range(8):
        a=annihilator(mode)
        if anticommutator(a,dagger(a))!=identity:raise ArithmeticError('Majorana codec failed CAR')
    obs=local_observables()
    rows={name:score_plan(op)for name,op in obs.items()}
    return {'local_commutators_checked':checks,'disjoint_odd_anticommutator_zero':True,
            'observable_readouts':rows,
            'tracial_control_g0_over_c_squared':obs['g0_over_c_squared'].get(0,_ZERO)[0],
            'tracial_control_g2_over_c_squared':obs['g2_over_c_squared'].get(0,_ZERO)[0],
            'tracial_control_is_a_ground_state':False}


def certificate():
    return {'source_symmetry':source_symmetry_audit(),'local_algebra':exact_algebra_audit(),
            'invariant_correlator_readouts':{k:score_plan(v)for k,v in invariant_observables().items()},
            'coefficient_identities':{'kappa':'4 q1','g0_over_c_squared':'9/2-12 qstar','g2_over_c_squared':'q2'},
            'uniform_Y_Gram_over_c_squared_upper':F(18),
            'C_is_kappa_times_actual_h_on_invariant_ground_average':True,
            'Y_Gram_symbol':'[g0+2 g2 sum_j cos(2 k_j)] I8 tensor Ispin',
            'residual_symbol':'[g0-6g2+(4g2-kappa^2 c^2) sum_j cos(k_j)^2] I',
            'coefficient_ground_values_assigned_from_finite_cube':False,
            'momentum_is_genuine_two_step_magnetic_translation':True,
            'complete_finite_phase_profiles_are_translation_eigenstates':False,
            'continuum_poles_or_physical_species_derived':False}
