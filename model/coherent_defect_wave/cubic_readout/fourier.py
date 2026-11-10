"""Exact translation-orbit compression of the side-four momentum projectors."""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product

from .lattice import State,ZERO,plus,scale,validate


@dataclass(frozen=True)
class FourierState:
    """Coefficients in an orthonormal momentum-orbit basis, not cell kets."""
    amplitudes:tuple
    norm_denominator:int
    momentum:tuple

    def __post_init__(self):
        momentum=tuple(self.momentum)
        if len(momentum)!=3 or any(type(k) is not int or k not in (0,1) for k in momentum):
            raise ValueError('three exact binary coarse-momentum labels required')
        if type(self.norm_denominator) is not int or self.norm_denominator<=0 or self.norm_denominator%8:
            raise ValueError('positive normalization divisible by the orbit size required')
        normalized=State(self.amplitudes,self.norm_denominator)
        if any(any(x>=2 for x in position) or any(any(x>=4 for x in cell) for cell,_ in field)
               for (position,_,field),_ in normalized.amplitudes):
            raise ValueError('canonical occupied position and side-four field required')
        object.__setattr__(self,'amplitudes',normalized.amplitudes)
        object.__setattr__(self,'momentum',momentum)

    def norm_squared(self):
        return sum((a*a+b*b for _,(a,b) in self.amplitudes),F(0))/self.norm_denominator

    def as_dict(self):
        return dict(self.amplitudes)


def compress(state,momentum):
    validate(state,4)
    momentum=tuple(momentum)
    if len(momentum)!=3 or any(type(k) is not int or k not in (0,1) for k in momentum):
        raise ValueError('three exact binary coarse-momentum labels required')
    out={}
    for (position,index,field),value in state.amplitudes:
        displacement=tuple(x//2 for x in position)
        canonical=tuple(x%2 for x in position)
        changed=tuple(sorted((tuple((x-2*d)%4 for x,d in zip(cell,displacement,strict=True)),g)
                             for cell,g in field))
        key=canonical,index,changed
        sign=(-1)**sum(k*d for k,d in zip(momentum,displacement,strict=True))
        out[key]=plus(out.get(key,ZERO),scale(value,sign))
    return FourierState(tuple(sorted((k,v) for k,v in out.items() if v!=ZERO)),
                        8*state.norm_denominator,momentum)


def expand(state):
    """Embed a projected orbit vector back into the full one-wave cell basis."""
    if not isinstance(state,FourierState) or state.norm_denominator%8:
        raise ValueError('Fourier projection with its original normalization required')
    out={}
    for (position,index,field),value in state.amplitudes:
        for displacement in product(range(2),repeat=3):
            def moved(x):return tuple((a+2*b)%4 for a,b in zip(x,displacement,strict=True))
            sign=(-1)**sum(k*d for k,d in zip(state.momentum,displacement,strict=True))
            key=moved(position),index,tuple(sorted((moved(cell),g) for cell,g in field))
            out[key]=scale(value,F(sign,8))
    return State(out,state.norm_denominator//8)


def phase_branch(prepared,advanced,quadrature,sign):
    if (prepared.momentum!=advanced.momentum or
            prepared.norm_denominator!=advanced.norm_denominator):
        raise ValueError('one momentum and preparation normalization required')
    if quadrature not in (0,1) or sign not in (-1,1):
        raise ValueError('two quadratures and two signs required')
    out={key:scale(value,F(1,2)) for key,value in prepared.amplitudes}
    for key,(real,imag) in advanced.amplitudes:
        value=(real,imag) if quadrature==0 else (-imag,real)
        out[key]=plus(out.get(key,ZERO),scale(value,F(sign,2)))
    return FourierState(tuple(sorted((k,v) for k,v in out.items() if v!=ZERO)),
                        2*prepared.norm_denominator,prepared.momentum)
