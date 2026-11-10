"""Declared common-link CAR spin extension and native Dirac propagation.

The additional four-component memory is architectural. Primitive carrier
port, address and internal labels are retained, not relabelled as this spin.
"""
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
import sympy as s

from ..native_routing_join.model import (State,label_word,label_memory,
    word_label,native_sign,seed,mode)
from ..mobile_flux.model import observer,STAR,representative
from ..species.fock import MODES,wedge_insert,coordinates
from ..spatial.dirac_route import route_step,clifford
from .model import flux_label_column

EXT_MODES=4*MODES
PROGRAM_ID='coherent-spinful-propagating-native-flux-v1'


@dataclass(frozen=True)
class SpinFock:
    amplitudes:tuple
    norm_denominator:int=1

    def __post_init__(self):
        if type(self.norm_denominator) is not int or self.norm_denominator<=0:
            raise ValueError('positive exact normalization required')
        values=self.amplitudes.items() if isinstance(self.amplitudes,dict) else self.amplitudes
        out={}
        for (word,occupied),coefficient in values:
            occupied=tuple(occupied)
            if word not in observer().base.link_words or tuple(sorted(set(occupied)))!=occupied:
                raise ValueError('literal link word and ordered CAR occupation required')
            if any(type(i) is not int or i not in range(EXT_MODES) for i in occupied):
                raise ValueError('native mode with four-component spin memory required')
            coefficient=s.sympify(coefficient)
            if coefficient.has(s.Float) or coefficient.is_finite is not True or coefficient.is_number is not True:
                raise ValueError('exact finite coefficient required')
            out[word,occupied]=s.expand(out.get((word,occupied),0)+coefficient)
        object.__setattr__(self,'amplitudes',tuple(sorted((k,v) for k,v in out.items() if v!=0)))

    def norm_squared(self):
        return s.simplify(sum((s.conjugate(v)*v for _,v in self.amplitudes),s.Integer(0))/self.norm_denominator)


def spin_sign(occupation):
    return native_sign(tuple(i//4 for i in occupation))


def gamma(state):
    return SpinFock({(word^STAR,n):spin_sign(n)*v for (word,n),v in state.amplitudes},state.norm_denominator)


def wave(state,inverse=False):
    out={}
    for (word,occupied),value in state.amplitudes:
        terms={():value}
        for index in occupied:
            primitive,spin=divmod(index,4);following={}
            for wedge,amplitude in terms.items():
                for native,weight in observer().column(word,primitive,inverse):
                    inserted,sign=wedge_insert(wedge,4*native+spin)
                    if sign:following[inserted]=s.expand(following.get(inserted,0)+sign*amplitude*weight)
            terms={k:v for k,v in following.items() if v!=0}
        for wedge,value in terms.items():out[word,wedge]=s.expand(out.get((word,wedge),0)+value)
    return SpinFock(out,state.norm_denominator)


def controls(state,inverse=False):
    out={};law=observer().base.law
    edges=reversed(law.interior) if inverse else law.interior
    edges=tuple(edges)
    for (word,occupied),value in state.amplitudes:
        target=word
        for edge in edges:
            endpoints=law.edges[edge]
            if sum(coordinates(index//4)[0] in endpoints for index in occupied)%2:
                target=law.transport(target,edge)[0]
        out[target,occupied]=s.expand(out.get((target,occupied),0)+value)
    return SpinFock(out,state.norm_denominator)


def native_step(state,inverse=False):
    return wave(controls(state,True),True) if inverse else controls(wave(state))


def encode(label,occupation):
    word=label_word(label)
    return SpinFock({(word,occupation):1,(word^STAR,occupation):spin_sign(occupation)},2)


@lru_cache(maxsize=4096)
def native_column(label,occupation,inverse=False):
    source=encode(label,occupation);image=native_step(source,inverse)
    out={}
    for (word,wedge),value in image.amplitudes:
        if word==representative(word):out[word_label(word,label_memory(label)),wedge]=value
    reconstructed={}
    for (target,wedge),value in out.items():
        assert bool(target)==bool(label) and len(wedge)==len(occupation)
        for key,weight in encode(target,wedge).amplitudes:
            reconstructed[key]=s.expand(reconstructed.get(key,0)+value*weight)
    assert SpinFock(reconstructed,2)==image and gamma(image)==image
    return tuple(sorted(out.items()))


def norm(state):return s.simplify(sum((s.conjugate(v)*v for v in state.values()),s.Integer(0)))


def clean(state):return {k:s.expand(v) for k,v in state.items() if s.expand(v)!=0}


def add(states,coefficients):
    out={}
    for state,coefficient in zip(states,coefficients,strict=True):
        for key,value in state.items():out[key]=s.expand(out.get(key,0)+coefficient*value)
    return clean(out)


@lru_cache(maxsize=256)
def flux_column(labels,length,spacing,inverse=False):
    return flux_label_column(labels,length,spacing,inverse)


def flux_route(state,length,spacing,inverse=False):
    out={}
    for (labels,natives),value in state.items():
        for target,weight in flux_column(labels,length,spacing,inverse):
            out[target,natives]=s.expand(out.get((target,natives),0)+value*weight)
    return clean(out)


@lru_cache(maxsize=8192)
def kinetic_column(index,length,inverse=False):
    cell_lane,local=divmod(index,EXT_MODES);primitive,spin=divmod(local,4)
    cells=tuple(product(range(length),repeat=3));lookup={x:i for i,x in enumerate(cells)}
    if cell_lane not in range(2*len(cells)):raise ValueError('native index on cube required')
    position=cells[cell_lane//2];owner=cell_lane%2
    image=route_step({(position,spin):1},length,0,inverse)
    return tuple(sorted(((2*lookup[x]+owner)*EXT_MODES+4*primitive+target,value)
                        for (x,target),value in image.items()))


def native_route(state,length,inverse=False):
    """Full exterior native Dirac step; no occupation sector is truncated."""
    out={}
    for (labels,occupied),value in state.items():
        if len(labels)!=2*length**3 or tuple(sorted(set(occupied)))!=occupied:
            raise ValueError('complete cube labels and ordered native CAR occupation required')
        terms={():value}
        for index in occupied:
            following={}
            for wedge,amplitude in terms.items():
                for target,weight in kinetic_column(index,length,inverse):
                    inserted,sign=wedge_insert(wedge,target)
                    if sign:following[inserted]=s.expand(following.get(inserted,0)+sign*amplitude*weight)
            terms=clean(following)
        for wedge,value in terms.items():out[labels,wedge]=s.expand(out.get((labels,wedge),0)+value)
    return clean(out)


def local_native(state,lane,inverse=False,flat=False):
    out={}
    for (labels,occupied),value in state.items():
        local=tuple(i%EXT_MODES for i in occupied if i//EXT_MODES==lane)
        spectators=tuple(i for i in occupied if i//EXT_MODES!=lane)
        for (target,wedge),weight in native_column(0 if flat else labels[lane],local,inverse):
            if flat:target=labels[lane]
            moved=list(labels);moved[lane]=target
            combined=tuple(sorted((*spectators,*(lane*EXT_MODES+i for i in wedge))))
            key=tuple(moved),combined;out[key]=s.expand(out.get(key,0)+value*weight)
    return clean(out)


def generator(state,flat=False):
    lanes={i//EXT_MODES for (_,n),_ in state.items() for i in n}
    out={}
    for lane in sorted(lanes):
        out=add((out,state,local_native(state,lane,flat=flat),local_native(state,lane,True,flat)),
                (1,1,-s.Rational(1,2),-s.Rational(1,2)))
    return out


def contact(state):
    """Declared flat-background subtraction, retaining the literal link response."""
    return add((generator(state),generator(state,True)),(1,-1))


def focused_preparation(length=4):
    """Inverse native propagation prepares an incoming wave that focuses at receiver."""
    cells=tuple(product(range(length),repeat=3));lookup={x:i for i,x in enumerate(cells)}
    labels=[0]*(2*length**3);labels[2*lookup[1,1,1]]=1
    receiver=2*lookup[2,2,2]
    target={ (tuple(labels),(receiver*EXT_MODES+4*mode(7,3,0,0),)):s.Integer(1)}
    return native_route(target,length,True),target,receiver


def changed_face(key):return any(label>2 for label in key[0])


def complete_records(state,receiver):
    return { (arrived,changed):norm({key:value for key,value in state.items()
                if any(i//EXT_MODES==receiver for i in key[1])==arrived
                and changed_face(key)==changed})
            for arrived in (False,True) for changed in (False,True)}
