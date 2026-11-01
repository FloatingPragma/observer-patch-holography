"""Actual-source CAR ladders and coherent ground-family intertwiners.

All coefficients are integer pairs a+b*sqrt(3), with common amplitude
normalization 36. No Clebsch-Gordan table or numerical eigenvector is used.
"""
from functools import lru_cache
from math import comb, lcm
from fractions import Fraction as F

import numpy as np
from scipy import sparse

from ..cube_eta_protection import cube
from ..interacting_vacuum import ground
from ..retarded_matter_response import response


NORMS=(1,4,24,144,576)


def combine(*terms):
    out={}
    for coefficient,state in terms:
        for mask,(a,b) in state.items():
            x,y=out.get(mask,(0,0));out[mask]=(x+coefficient*a,y+coefficient*b)
    return {mask:value for mask,value in out.items() if value!=(0,0)}


def car(state,mode,creation):
    operation=ground._create if creation else ground._annihilate
    out={}
    for mask,(a,b) in state.items():
        value=operation(mask,mode)
        if value is not None:out[value[0]]=(value[1]*a,value[1]*b)
    return out


def ladder(state,lower=True):
    return combine(*[(1,car(car(state,2*x+(0 if lower else 1),False),
                           2*x+(1 if lower else 0),True)) for x in range(8)])


def inner(left,right):
    a=b=0
    for mask,(x,y) in left.items():
        u,v=right.get(mask,(0,0));a+=x*u+3*y*v;b+=x*v+y*u
    return a,b


@lru_cache(maxsize=1)
def ground_family():
    states=[ground.reference()]
    for _ in range(4):states.append(ladder(states[-1]))
    if ladder(states[0],False) or ladder(states[-1]):
        raise ArithmeticError('full source ground highest/lowest weight failed')
    for k,state in enumerate(states):
        if inner(state,state)!=(1296*NORMS[k],0):
            raise ArithmeticError('source ground ladder norm changed')
        if k and ladder(state,False)!=combine((k*(5-k),states[k-1])):
            raise ArithmeticError('full source spin-two raising coefficient changed')
    masks=cube.weights(4);lookup={int(mask):j for j,mask in enumerate(masks)}
    h=cube.restrict(cube.operators().hubbard_H.real,4)
    for state in states:
        a=np.zeros(len(masks),dtype=np.int64);b=np.zeros_like(a)
        for mask,(x,y) in state.items():a[lookup[mask]]=x;b[lookup[mask]]=y
        if not np.array_equal(h@a,-12*b) or not np.array_equal(h@b,-4*a):
            raise ArithmeticError('spin-lowered reference is not an exact original ground')
    return tuple(states)


@lru_cache(maxsize=None)
def raising(up,down):
    if down<1:raise ValueError('raising map requires a down fermion')
    source=response.sector_masks(up,down)
    target=response.sector_masks(up+1,down-1)
    lookup={int(mask):j for j,mask in enumerate(target)}
    rows=[];columns=[];values=[]
    for j,mask in enumerate(source):
        for x in range(8):
            first=ground._annihilate(int(mask),2*x+1)
            if first is None:continue
            second=ground._create(first[0],2*x)
            if second is None:continue
            rows.append(lookup[second[0]]);columns.append(j);values.append(first[1]*second[1])
    return sparse.csr_matrix((values,(rows,columns)),shape=(len(target),len(source)),dtype=np.int64)


def representation_certificate(number):
    if number not in (3,5):raise ValueError('retained charged particle number is three or five')
    hs=[response.source_block(up,number-up) for up in range(number+1)]
    dimensions=[h.shape[0] for h in hs]
    for up,h in enumerate(hs):
        minus_plus=sparse.csr_matrix(h.shape,dtype=np.int64)
        plus_minus=sparse.csr_matrix(h.shape,dtype=np.int64)
        if up<number:
            t=raising(up,number-up)
            if (hs[up+1]@t-t@h).nnz:raise ArithmeticError('source H fails spin raising intertwiner')
            minus_plus=t.T@t
        if up:
            t=raising(up-1,number-up+1);plus_minus=t@t.T
        if (plus_minus-minus_plus-(2*up-number)*sparse.eye(h.shape[0],dtype=np.int64)).nnz:
            raise ArithmeticError('complete number-sector spin commutator failed')
    spins=range(number,0,-2)
    multiplicities={str(twice)+'/2':comb(8,(number+twice)//2)*comb(8,(number-twice)//2)
                    -(comb(8,(number+twice)//2+1)*comb(8,(number-twice)//2-1)
                      if (number-twice)//2 else 0) for twice in spins}
    if sum((int(key.split('/')[0])+1)*value for key,value in multiplicities.items())!=comb(16,number):
        raise ArithmeticError('spin multiplicities do not exhaust the full source')
    return {'N':number,'full_dimension':comb(16,number),'weight_space_dimensions':dimensions,
            'spin_multiplicities':multiplicities,'all_source_spin_intertwiners_exact':True,
            'spin_commutator_exact_on_entire_number_sector':True}


def projection_certificate():
    t=raising(4,1)
    down=response.dilation('ground','down').vectors
    up=tuple(v[:56] for v in response.dilation('ground','up').vectors)
    if (t@t.T-5*sparse.eye(56,dtype=np.int64)).nnz:
        raise ArithmeticError('spin-five-halves ladder normalization changed')
    if any(not np.array_equal(t@d,u) for d,u in zip(down,up,strict=True)):
        raise ArithmeticError('source S+ c_down* G differs from c_up* G')
    low5=tuple(5*d-t.T@u for d,u in zip(down,up,strict=True))
    if any(np.any(t@x) for x in low5):raise ArithmeticError('claimed spin-three-halves seed is not highest weight')
    a=ground.band_data()[0];gram=response.gram(up,up)
    if (not np.array_equal(gram[0],648*np.eye(8,dtype=np.int64)) or
            not np.array_equal(gram[1],216*a)):
        raise ArithmeticError('full spatial up-addition Gram is not P_plus')
    h=response.source_block(5,0)
    if (not np.array_equal(h@up[0],-9*up[1]) or
            not np.array_equal(h@up[1],-3*up[0])):
        raise ArithmeticError('source spin-five-halves seed is not at gap lambda')
    removed=tuple(v[56:] for v in response.dilation('ground','up').vectors)
    h3=response.source_block(3,0);removed_gram=response.gram(removed,removed)
    if (not np.array_equal(h3@removed[0],-9*removed[1]) or
            not np.array_equal(h3@removed[1],-3*removed[0]) or
            not np.array_equal(removed_gram[0],648*np.eye(8,dtype=np.int64)) or
            not np.array_equal(removed_gram[1],-216*a)):
        raise ArithmeticError('source removal seed energy or full spatial Gram changed')
    return {'raising_shape':[56,560],'T_T_adjoint':'5 I_56','T_D':'U',
            'P_high_on_M_three_halves':'T_adjoint T/5',
            'D_high':'T_adjoint U/5','D_low':'D-T_adjoint U/5',
            'T_D_low':'0','U_Gram':'P_plus=(I+A/sqrt3)/2',
            'D_high_Gram':'P_plus/5','U_grounded_energy':'sqrt3/16',
            'removal_seed_Gram':'P_minus','removal_energy_cost':'sqrt3/16',
            'all_spatial_coherences_retained':True}


def family_certificate():
    states=ground_family();checked_add=checked_remove=0
    for x in range(8):
        up=car(states[0],2*x,True);down=car(states[0],2*x+1,True)
        low5=combine((5,down),(-1,ladder(up)))
        if ladder(low5,False):raise ArithmeticError('low-spin source seed fails S+ annihilation')
        upper=[up];lower=[low5];removed=[car(states[0],2*x,False)]
        for _ in range(5):upper.append(ladder(upper[-1]))
        for _ in range(4):lower.append(ladder(lower[-1]))
        for _ in range(4):removed.append(ladder(removed[-1]))
        if ladder(upper[-1]) or lower[4] or removed[4]:
            raise ArithmeticError('source spin ladder does not terminate')
        for k,state in enumerate(states):
            expected_down=combine((1,upper[k+1]),(1,lower[k]))
            expected_up=combine((5-k,upper[k]),(-k,lower[k-1] if k else {}))
            if (combine((5,car(state,2*x+1,True)))!=expected_down or
                    combine((5,car(state,2*x,True)))!=expected_up):
                raise ArithmeticError('complete coherent ground creator identity failed')
            expected_removed_down=combine((k,removed[k-1] if k else {}))
            if (car(state,2*x,False)!=removed[k] or
                    car(state,2*x+1,False)!=expected_removed_down):
                raise ArithmeticError('complete coherent ground annihilator identity failed')
            checked_add+=2;checked_remove+=2
    return {'ground_dimension':5,'source_ground_norms_squared_relative_to_G':list(NORMS),
            'creator_columns_checked':checked_add,'annihilator_columns_checked':checked_remove,
            'creator_identities':[
                '5 c_down* G_k = S_minus^(k+1) U + S_minus^k (5D-S_minus U)',
                '5 c_up* G_k = (5-k) S_minus^k U - k S_minus^(k-1) (5D-S_minus U)'],
            'annihilator_identities':['c_up G_k=S_minus^k B',
                                      'c_down G_k=k S_minus^(k-1) B'],
            'arbitrary_coherent_ground_and_linear_CAR_probe_covered':True,
            'normalized_squared_ladder_factors':[
                {'k':k,'up_high':str(F(5-k,5)),'up_low':str(F(k,4)),
                 'down_high':str(F(k+1,5)),'down_low':str(F(4-k,4))} for k in range(5)]}


def _columns(state,number,creation):
    masks=cube.weights(number);lookup={int(mask):j for j,mask in enumerate(masks)}
    a=np.zeros((len(masks),16),dtype=np.int64);b=np.zeros_like(a)
    for mode in range(16):
        for mask,(u,v) in car(state,mode,creation).items():a[lookup[mask],mode]=u;b[lookup[mask],mode]=v
    return a,b


def twirled_moments():
    common=lcm(*NORMS);denominator=5*1296*common
    sums={kind:[np.zeros((16,16),dtype=np.int64) for _ in range(4)] for kind in ('addition','removal')}
    for kind,number,creation,sign in (('addition',5,True,1),('removal',3,False,-1)):
        h=cube.restrict(cube.operators().hubbard_H.real,number)
        for norm,state in zip(NORMS,ground_family(),strict=True):
            a,b=_columns(state,number,creation)
            ka,kb=sign*(h@a+12*b),sign*(h@b+4*a)
            m0=response.gram((a,b),(a,b));m1=response.gram((a,b),(ka,kb))
            for target,value in zip(sums[kind],(*m0,*m1),strict=True):target+=(common//norm)*value
    a=np.kron(ground.band_data()[0],np.eye(2,dtype=np.int64));eye=np.eye(16,dtype=np.int64);d=denominator
    expected={'addition':(3*d//4*eye,d//12*a,4*d*eye+3*d//4*a,d//4*eye),
              'removal':(d//4*eye,-d//12*a,d//4*a,-d//4*eye)}
    for kind in sums:
        if any(not np.array_equal(v,w) for v,w in zip(sums[kind],expected[kind],strict=True)):
            raise ArithmeticError('direct full-ground source twirl disagrees with exact coherent moments')
    return {'all_16_by_16_spatial_and_spin_entries_checked':True,
            'normalized_ground_projector':'P_N4_ground/5',
            'addition_zeroth':'(3/4)I + (sqrt3/12)(A tensor I_spin)',
            'removal_zeroth':'(1/4)I - (sqrt3/12)(A tensor I_spin)',
            'retarded_zeroth':'I_16','retarded_first':'(A tensor I_spin)/16+I_16/4',
            'local_occupation_per_spin':'1/4','local_addition_probability':'3/4'}


def certificate():
    return {'projection':projection_certificate(),'ground_family':family_certificate(),
            'number_sectors':[representation_certificate(n) for n in (5,3)],
            'direct_twirled_moments':twirled_moments(),
            'hidden_total_spin_is_physical_Lorentz_spin':False}
