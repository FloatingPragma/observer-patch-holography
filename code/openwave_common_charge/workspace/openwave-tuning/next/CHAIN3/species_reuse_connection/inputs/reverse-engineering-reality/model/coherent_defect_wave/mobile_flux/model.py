"""Vacuum/paired-flux link sectors and gauge-compatible spatial transport.

The intercell gate changes individual carrier boundaries. Its conserved
charge counts complete inside/exterior flux pairs in the declared domain.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from ..model import JointObserver,SparseState
from ..species.fock import FockObserver,FockState,MODES,coordinates,mode

PROGRAM_ID='coherent-mobile-paired-flux-v1'
STAR=6553634
A=(1-s.I)/2
B=(1+s.I)/2


class MobileObserver(JointObserver):
    def __init__(self):
        super().__init__()
        self.original_link_words=self.link_words
        self.link_words=tuple(sorted({*self.link_words,0,STAR}))
        self.dimension=len(self.link_words)*7200
        self._transports={(g,e):self.law.transport(g,e)[0]
                          for g in self.link_words for e in self.law.interior}
        if len(self.link_words)!=12 or any(g not in self.link_words for g in self._transports.values()):
            raise ArithmeticError('the enlarged source domain is not invariant')


class MobileFockObserver(FockObserver):
    def __init__(self):
        super().__init__(control='parity')
        self.base=MobileObserver()

    def face_measurement(self,state):
        self._validate(state)
        flat=FockState({key:v for key,v in state.amplitudes
                       if self.base.law.region_charge(key[0])==0},state.norm_denominator)
        return (('flat',flat),*super().face_measurement(state))


@lru_cache(maxsize=1)
def observer():
    return MobileFockObserver()


def charge(word):
    if type(word) is not int or word not in observer().base.link_words:
        raise ValueError('mobile-flux source word required')
    return observer().base.law.region_charge(word)


def representative(word):
    return min(word,word^STAR)


def chi(occupied,cell):
    return (-1)**sum(i//MODES==cell and coordinates(i%MODES)[0]==5
                     and coordinates(i%MODES)[1]>=3 for i in occupied)


def frame(word,occupied,cell):
    return chi(occupied,cell) if word>(word^STAR) else 1


@dataclass(frozen=True)
class PairState:
    amplitudes:tuple
    norm_denominator:int=1

    def __post_init__(self):
        if type(self.norm_denominator) is not int or self.norm_denominator<=0:
            raise ValueError('positive integer normalization required')
        entries=self.amplitudes.items() if isinstance(self.amplitudes,dict) else self.amplitudes
        out={}
        for (links,occupied),value in entries:
            links,occupied=tuple(links),tuple(occupied)
            if len(links)!=2 or any(type(g) is not int or g not in observer().base.link_words for g in links):
                raise ValueError('two words of the enlarged link domain required')
            if tuple(sorted(set(occupied)))!=occupied or any(type(i) is not int or i not in range(2*MODES) for i in occupied):
                raise ValueError('ordered full two-cell CAR occupation required')
            value=s.sympify(value)
            if value.is_number is not True or value.is_finite is not True or value.has(s.Float):
                raise ValueError('exact finite amplitudes required')
            key=links,occupied;out[key]=s.expand(out.get(key,0)+value)
        object.__setattr__(self,'amplitudes',tuple(sorted((k,v) for k,v in out.items() if v!=0)))

    def as_dict(self):return dict(self.amplitudes)

    def norm_squared(self):
        return s.expand(sum((s.conjugate(v)*v for _,v in self.amplitudes),s.Integer(0))/self.norm_denominator)


def gamma(state,cell):
    if cell not in (0,1):raise ValueError('two-cell index required')
    out={}
    for (links,occupied),value in state.amplitudes:
        moved=list(links);moved[cell]^=STAR
        out[tuple(moved),occupied]=chi(occupied,cell)*value
    return PairState(out,state.norm_denominator)


def positive(state):
    """Normalized plus orbit from an input with definite gauge representatives."""
    for cell in (0,1):
        out=state.as_dict()
        for key,value in gamma(state,cell).amplitudes:
            out[key]=s.expand(out.get(key,0)+value)
        state=PairState(out,2*state.norm_denominator)
    if state.norm_squared()!=1:
        raise ValueError('this preparation requires a normalized definite-word input')
    return state


def constraint_branches(state):
    result=[]
    for left in (1,-1):
        for right in (1,-1):
            out={}
            for coefficient,part in ((1,state),(left,gamma(state,0)),
                                     (right,gamma(state,1)),(left*right,gamma(gamma(state,0),1))):
                for key,value in part.amplitudes:
                    out[key]=s.expand(out.get(key,0)+s.Rational(coefficient,4)*value)
            result.append(((left,right),PairState(out,state.norm_denominator)))
    return tuple(result)


def swap_links(state,keep_gauge_bits=True):
    """Dress either orbit exchange or complete link-word exchange by F."""
    out={}
    for (links,occupied),value in state.amplitudes:
        before=frame(links[0],occupied,0)*frame(links[1],occupied,1)
        if keep_gauge_bits:
            moved=tuple(representative(links[1-cell])^(STAR if links[cell]>(links[cell]^STAR) else 0)
                        for cell in (0,1))
        else:
            moved=tuple(reversed(links))
        after=frame(moved[0],occupied,0)*frame(moved[1],occupied,1)
        key=moved,occupied
        out[key]=s.expand(out.get(key,0)+before*after*value)
    return PairState(out,state.norm_denominator)


def exchange(state,inverse=False,keep_gauge_bits=True):
    aa,bb=(s.conjugate(A),s.conjugate(B)) if inverse else (A,B)
    out={key:aa*value for key,value in state.amplitudes}
    for key,value in swap_links(state,keep_gauge_bits).amplitudes:
        out[key]=s.expand(out.get(key,0)+bb*value)
    return PairState(out,state.norm_denominator)


def onsite(state,cell,inverse=False):
    if cell not in (0,1):raise ValueError('two-cell index required')
    out={}
    for (links,occupied),value in state.amplitudes:
        active=tuple(i%MODES for i in occupied if i//MODES==cell)
        spectators=tuple(i for i in occupied if i//MODES!=cell)
        local=FockState({(links[cell],active):value})
        advanced=observer().inverse_step(local) if inverse else observer().step(local)
        for (word,wedge),coefficient in advanced.amplitudes:
            words=list(links);words[cell]=word
            joined=tuple(sorted((*spectators,*(cell*MODES+i for i in wedge))))
            key=tuple(words),joined
            out[key]=s.expand(out.get(key,0)+coefficient)
    return PairState(out,state.norm_denominator)


def tick(state,inverse=False):
    if inverse:return onsite(onsite(exchange(state,True),1,True),0,True)
    return exchange(onsite(onsite(state,0),1))


def charge_measurement(state):
    return tuple(((left,right),PairState({key:v for key,v in state.amplitudes
                    if (charge(key[0][0]),charge(key[0][1]))==(left,right)},state.norm_denominator))
                 for left in (0,1) for right in (0,1))


def insert_pair(state,cell):
    """Complete reversible preparation: exchange the flat and initial-pair orbits."""
    if cell not in (0,1):raise ValueError('two-cell index required')
    out={}
    for (links,occupied),value in state.amplitudes:
        words=list(links)
        if representative(words[cell]) in (0,1):words[cell]^=1
        out[tuple(words),occupied]=value
    return PairState(out,state.norm_denominator)
