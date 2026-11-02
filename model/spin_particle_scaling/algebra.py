"""Exact affine spatial action in the source and hidden spin frames."""
from itertools import permutations, product
import numpy as np
import sympy as sp


def coordinate(x):
    if type(x) is not tuple or len(x) != 3 or any(type(v) is not int for v in x):
        raise ValueError('three exact integer coordinates required')
    return x


def frame(x):
    x=coordinate(x)
    return np.diag([(-1)**((sum(x)-x[a])%2) for a in range(3)]).astype(np.int64)


def order_translation(x):
    return np.diag([(-1)**(v%2) for v in coordinate(x)]).astype(np.int64)


def rotations():
    out=[]
    for p in permutations(range(3)):
        for s in product((-1,1),repeat=3):
            r=np.zeros((3,3),dtype=np.int64)
            for a in range(3):r[a,p[a]]=s[a]
            det=(int(r[0,0])*(int(r[1,1])*int(r[2,2])-int(r[1,2])*int(r[2,1]))
                 -int(r[0,1])*(int(r[1,0])*int(r[2,2])-int(r[1,2])*int(r[2,0]))
                 +int(r[0,2])*(int(r[1,0])*int(r[2,1])-int(r[1,1])*int(r[2,0])))
            if det==1:out.append(r)
    return tuple(out)


def local_car_frame():
    """Both empty/double states are retained in each local spin bilinear."""
    I=sp.eye(2);sigma=(sp.Matrix([[0,1],[1,0]]),sp.Matrix([[0,-sp.I],[sp.I,0]]),sp.diag(1,-1))
    checked=0
    for x in product(range(4),repeat=3):
        R=I
        for a in range(3):R=R*(sp.I*sigma[a])**x[a]
        F=frame(x)
        U=sp.diag(1,1,1,1);U[1:3,1:3]=R
        for a in range(3):
            S=sp.zeros(4);S[1:3,1:3]=sigma[a]/2
            if sp.simplify(U.H*S*U-int(F[a,a])*S)!=sp.zeros(4):
                raise ArithmeticError('full onsite frame mismatch')
            checked+=1
    return checked


def certificate():
    e=np.array([1,0,0],dtype=np.int64);checked=0
    for r in rotations():
        for t in product(range(2),repeat=3):
            for y in product((-1,0,1),repeat=3):
                x=tuple(int(v) for v in r@np.array(y)+t)
                if not np.array_equal(frame(x)@r@frame(y),frame(t)@r):
                    raise ArithmeticError('affine hidden action failed')
                if not np.array_equal((-1)**(sum(t)%2)*frame(t)@r,order_translation(t)@r):
                    raise ArithmeticError('staggered order action failed')
                checked+=1
    stabilizer=[{'rotation':r.tolist(),'translation_parity':list(t)}
                for r in rotations() for t in product(range(2),repeat=3)
                if np.array_equal(order_translation(t)@r@e,e)]
    if len(stabilizer)!=32:raise ArithmeticError('wrong ordered spatial stabilizer')
    C4=np.array([[1,0,0],[0,0,-1],[0,1,0]],dtype=np.int64)
    c=sp.Matrix(C4);vp=sp.Matrix([0,1,sp.I]);vm=sp.conjugate(vp)
    if c*vp!=-sp.I*vp or c*vm!=sp.I*vm:raise ArithmeticError('wrong C4 circular characters')
    # Source s^a Fourier carrier = hidden carrier + M_a, modulo 2pi.
    carriers=[tuple(0 if a==j else 1 for j in range(3)) for a in range(3)]
    ordered=[tuple((1+v)%2 for v in k) for k in carriers]
    if ordered!=[(1,0,0),(0,1,0),(0,0,1)]:raise ArithmeticError('wrong order wavevectors')
    return {'affine_cases_exact':checked,'full_onsite_CAR_frame_cases':local_car_frame(),
            'proper_cubic_rotations':24,'ordered_spatial_stabilizer_mod2_count':len(stabilizer),
            'ordered_spatial_stabilizer_mod2':stabilizer,
            'order_translation_signs':'diag((-1)^t1,(-1)^t2,(-1)^t3)',
            'hidden_translation_signs':'diag((-1)^(t2+t3),(-1)^(t1+t3),(-1)^(t1+t2))',
            'source_order_carriers_in_pi_units':[list(k) for k in ordered],
            'source_hidden_zero_carriers_in_pi_units':[list(k) for k in carriers],
            'residual_translation_periods':[2,1,1],
            'transverse_residual_zone_centers_in_pi_units':[[0,0,1],[0,1,0]],
            'circular_C4_characters':['-i','+i'],
            'circular_channels_are_single_translation_eigenstates':False,
            'cubic_character_is_relativistic_helicity':False}
