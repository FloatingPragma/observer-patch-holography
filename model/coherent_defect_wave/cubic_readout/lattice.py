"""Sparse one-occupation execution in the exact positive-Gauss encoding.

Every cell retains its own five physical link-orbit values. Only values that
differ from the fixed background are stored. This is an isometric encoding
of the declared spatial CAR circuit, not a replacement background potential.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

from ..model import JointObserver, SparseState
from ..species.fock import mode, coordinates

STAR=6553634
ZERO=(F(0),F(0))


@lru_cache(maxsize=1)
def observer():
    return JointObserver()


@lru_cache(maxsize=1)
def representatives():
    return tuple(g for g in observer().link_words if g<(g^STAR))


def chi(index):
    p,j,_,_=coordinates(index)
    return -1 if p==5 and j>=3 else 1


@lru_cache(maxsize=8192)
def local_column(word,index,inverse=False):
    """J* U J in the five-link-orbit physical basis, checked against U."""
    if word not in representatives():
        raise ValueError('canonical positive-Gauss link representative required')
    w=coordinates(index)
    source=SparseState({(word,*w):F(1),(word^STAR,*w):F(chi(index))},2)
    out=observer().inverse_step(source) if inverse else observer().step(source)
    values=out.as_dict(); reduced=[]
    for (g,p,j,r,a),value in out.amplitudes:
        target=mode(p,j,r,a)
        if g<(g^STAR):
            if values.get((g^STAR,p,j,r,a),F(0))!=chi(target)*value:
                raise ArithmeticError('literal U leaves the encoded positive-Gauss sector')
            reduced.append(((g,target),F(value)))
    if sum(v*v for _,v in reduced)!=1:
        raise ArithmeticError('encoded source column is not normalized')
    return tuple(reduced)


def plus(a,b):
    return a[0]+b[0],a[1]+b[1]


def scale(a,c):
    return a[0]*c,a[1]*c


def inner_component(a,b):
    return a[0]*b[0]+a[1]*b[1],a[0]*b[1]-a[1]*b[0]


@dataclass(frozen=True)
class State:
    amplitudes:tuple
    norm_denominator:int=1

    def __post_init__(self):
        if type(self.norm_denominator) is not int or self.norm_denominator<=0:
            raise ValueError('positive integer normalization required')
        entries=self.amplitudes.items() if isinstance(self.amplitudes,dict) else self.amplitudes
        out={}
        for (position,index,changes),amplitude in entries:
            position=tuple(position); changes=tuple(changes)
            coordinates(index)
            if len(position)!=3 or any(type(v) is not int or v<0 for v in position):
                raise ValueError('nonnegative three-dimensional position required')
            if tuple(sorted(changes))!=changes or len({x for x,_ in changes})!=len(changes):
                raise ValueError('sorted distinct changed-cell keys required')
            if any(len(x)!=3 or any(type(v) is not int or v<0 for v in x)
                   or g not in representatives() or g==1 for x,g in changes):
                raise ValueError('store only admitted nonbackground link orbits')
            if len(amplitude)!=2 or any(type(v) not in (int,F) for v in amplitude):
                raise ValueError('exact Gaussian rational amplitudes required')
            value=tuple(map(F,amplitude)); key=position,index,changes
            out[key]=plus(out.get(key,ZERO),value)
        object.__setattr__(self,'amplitudes',tuple(sorted((k,v) for k,v in out.items() if v!=ZERO)))

    def as_dict(self):
        return dict(self.amplitudes)

    def norm_squared(self):
        return sum((r*r+i*i for _,(r,i) in self.amplitudes),F(0))/self.norm_denominator


def validate(state,length):
    if type(length) is not int or length<4 or length%2:
        raise ValueError('even periodic side length at least four required')
    if not isinstance(state,State):
        raise ValueError('encoded one-occupation state required')
    if any(any(v>=length for v in x) for (p,_,changes),_ in state.amplitudes
           for x in (p,*(x for x,_ in changes))):
        raise ValueError('position outside the declared torus')


def onsite(state,inverse=False):
    output={}
    for (position,index,changes),value in state.amplitudes:
        field=dict(changes); word=field.get(position,1)
        for (g,target),weight in local_column(word,index,inverse):
            changed=dict(field)
            if g==1:
                changed.pop(position,None)
            else:
                changed[position]=g
            key=position,target,tuple(sorted(changed.items()))
            output[key]=plus(output.get(key,ZERO),scale(value,weight))
    return State(output,state.norm_denominator)


def hopping_layer(state,length,axis,parity,inverse=False):
    validate(state,length)
    if axis not in (0,1,2) or parity not in (0,1):
        raise ValueError('one of six cubic matching layers required')
    output={}
    for (position,index,changes),(real,imag) in state.amplitudes:
        target=list(position)
        target[axis]=(target[axis]+(1 if target[axis]%2==parity else -1))%length
        if inverse:
            stay=((real-imag)/2,(imag+real)/2)
            hop=((real+imag)/2,(imag-real)/2)
        else:
            stay=((real+imag)/2,(imag-real)/2)
            hop=((real-imag)/2,(imag+real)/2)
        for key,value in (((position,index,changes),stay),((tuple(target),index,changes),hop)):
            output[key]=plus(output.get(key,ZERO),value)
    return State(output,state.norm_denominator)


def tick(state,length,inverse=False):
    validate(state,length)
    layers=tuple(product(range(3),range(2)))
    if inverse:
        for axis,parity in reversed(layers):
            state=hopping_layer(state,length,axis,parity,True)
        return onsite(state,True)
    state=onsite(state)
    for axis,parity in layers:
        state=hopping_layer(state,length,axis,parity)
    return state


def translate(state,length,shift):
    validate(state,length)
    shift=tuple(shift)
    if len(shift)!=3 or any(type(x) is not int for x in shift):
        raise ValueError('three integer cell displacements required')
    def moved(x):return tuple((a+b)%length for a,b in zip(x,shift,strict=True))
    return State({(moved(p),i,tuple(sorted((moved(x),g) for x,g in field))):v
                  for (p,i,field),v in state.amplitudes},state.norm_denominator)


def packet(index=None,origin=(0,0,0)):
    index=mode(7,3,0,0) if index is None else index
    return State({(tuple(a+b for a,b in zip(origin,epsilon,strict=True)),index,()):(1,0)
                  for epsilon in product(range(2),repeat=3)},8)


def overlap(left,right):
    if left.norm_denominator!=right.norm_denominator:
        raise ValueError('one unchanged normalization required')
    lookup=right.as_dict(); result=ZERO
    for key,value in left.amplitudes:
        result=plus(result,inner_component(value,lookup.get(key,ZERO)))
    return scale(result,F(1,left.norm_denominator))


def projector(state,momentum):
    """Complete Z2^3 cell-translation projector on the side-four torus."""
    if len(momentum)!=3 or any(k not in (0,1) for k in momentum):
        raise ValueError('three momentum labels in {0,pi} required')
    out={}
    for d in product(range(2),repeat=3):
        sign=(-1)**sum(k*x for k,x in zip(momentum,d,strict=True))
        shifted=translate(state,4,tuple(2*x for x in d))
        for key,value in shifted.amplitudes:
            out[key]=plus(out.get(key,ZERO),scale(value,F(sign,8)))
    return State(out,state.norm_denominator)
