"""Local Gauss-preserving exchange and a finite three-dimensional family.

The exact two-cell interpreter uses the owning joint U inside each cell.
The family specifies finite-depth nearest-neighbor gates on every cell.
"""
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import product

import sympy as s

from ..model import JointObserver, SparseState
from ..gauss import gamma_amplitudes

A = (1-s.I)/2
B = (1+s.I)/2
STAR = 6553634


@lru_cache(maxsize=4096)
def expanded(value):
    return s.expand(value)


@lru_cache(maxsize=4096)
def weight(value):
    return s.expand(s.conjugate(value)*value)


@dataclass(frozen=True)
class PairState:
    amplitudes: tuple
    norm_denominator: int = 1

    def __post_init__(self):
        if type(self.norm_denominator) is not int or self.norm_denominator <= 0:
            raise ValueError('positive integer normalization required')
        entries = self.amplitudes.items() if isinstance(self.amplitudes, dict) else self.amplitudes
        out = {}
        observer = JointObserver()
        for (left, right), value in entries:
            observer._key(*left)
            observer._key(*right)
            exact = s.sympify(value)
            if (isinstance(value, (float, complex)) or exact.has(s.Float)
                    or exact.is_number is not True or exact.is_finite is not True):
                raise ValueError('exact amplitudes required')
            out[(left,right)] = expanded(out.get((left,right),0)+value)
        object.__setattr__(self, 'amplitudes', tuple(sorted((key,value) for key,value in out.items() if value != 0)))

    def as_dict(self):
        return dict(self.amplitudes)

    def norm_squared(self):
        return expanded(sum((weight(value) for _,value in self.amplitudes),s.Integer(0))
                        /self.norm_denominator)


def tensor(left, right):
    return PairState({(a,b): x*y for a,x in left.amplitudes for b,y in right.amplitudes},
                     left.norm_denominator*right.norm_denominator)


def frame_sign(key):
    g,p,j,_,_ = key
    return -1 if g > (g ^ STAR) and p == 5 and j >= 3 else 1


def dressed_swap(state):
    """Exchange wave registers, retaining local links and transporting frames."""
    out = {}
    for (left,right),value in state.amplitudes:
        a,b = (left[0],*right[1:]),(right[0],*left[1:])
        sign = frame_sign(left)*frame_sign(right)*frame_sign(a)*frame_sign(b)
        out[(a,b)] = sign*value
    return PairState(out,state.norm_denominator)


def exchange(state, inverse=False):
    """Square root of the dressed SWAP; both Gauss constraints are preserved."""
    if type(inverse) is not bool:
        raise ValueError('inverse flag must be boolean')
    swapped = dressed_swap(state).as_dict()
    source = state.as_dict()
    a,b = (s.conjugate(A),s.conjugate(B)) if inverse else (A,B)
    return PairState({key: expanded(a*source.get(key,0)+b*swapped.get(key,0))
                      for key in source.keys() | swapped.keys()},state.norm_denominator)


def gamma(state, cell):
    if type(cell) is not int or cell not in (0,1):
        raise ValueError('two-cell index required')
    out = {}
    for keys,value in state.amplitudes:
        key = keys[cell]
        g,p,j,r,a = key
        pair = list(keys)
        pair[cell] = (g ^ STAR,p,j,r,a)
        out[tuple(pair)] = (-value if p == 5 and j >= 3 else value)
    return PairState(out,state.norm_denominator)


def positive_local(state):
    """Normalized positive Gauss preparation for an input with one link word."""
    if len({key[0] for key,_ in state.amplitudes}) != 1:
        raise ValueError('the source preparation requires one definite link word')
    out = state.as_dict()
    for key,value in gamma_amplitudes(state.as_dict()).items():
        out[key] = out.get(key,0)+value
    return SparseState(out,2*state.norm_denominator)


def local_step(state, cell, inverse=False):
    """Execute the literal owning U on one cell, retaining the other as reference."""
    if type(cell) is not int or cell not in (0,1):
        raise ValueError('two-cell index required')
    if type(inverse) is not bool:
        raise ValueError('inverse flag must be boolean')
    observer = JointObserver()
    groups = defaultdict(dict)
    for keys,value in state.amplitudes:
        groups[keys[1-cell]][keys[cell]] = value
    out = {}
    for reference,amplitudes in groups.items():
        local = SparseState(amplitudes,state.norm_denominator)
        result = observer.inverse_step(local) if inverse else observer.step(local)
        for key,value in result.amplitudes:
            pair = (key,reference) if cell == 0 else (reference,key)
            out[pair] = value
    return PairState(out,state.norm_denominator)


def face_measurement(state, cell):
    """All five regional face outcomes; no raw link representative is measured."""
    if type(cell) is not int or cell not in (0,1):
        raise ValueError('two-cell index required')
    law = JointObserver().law
    return tuple(({'cell':cell,'face_index':face,'face':list(law.faces[face])},
                  PairState({keys:value for keys,value in state.amplitudes
                             if law.syndrome(keys[cell][0]) >> face & 1},state.norm_denominator))
                 for face in law.region)


def charge_measurement(state, cell):
    if type(cell) is not int or cell not in (0,1):
        raise ValueError('two-cell index required')
    return tuple(({'cell':cell,'charge':charge},
                  PairState({keys:value for keys,value in state.amplitudes
                             if (2 if keys[cell][2] < 3 else -3) == charge},state.norm_denominator))
                 for charge in (2,-3))


class CubicFamily:
    """A declared finite cubic placement and six fixed matching layers."""
    def __init__(self, length, periodic=False):
        if type(length) is not int or length < 2:
            raise ValueError('integer side length at least two required')
        if type(periodic) is not bool:
            raise ValueError('periodic flag must be boolean')
        if periodic and (length < 4 or length % 2):
            raise ValueError('periodic matching family needs an even side at least four')
        self.length,self.periodic = length,periodic
        self.sites = tuple(product(range(length),repeat=3))
        self.layers = tuple(self._matching(axis,parity) for axis in range(3) for parity in (0,1))

    def _matching(self,axis,parity):
        edges = []
        for x in self.sites:
            if x[axis] % 2 != parity:
                continue
            y = list(x)
            y[axis] += 1
            if y[axis] == self.length:
                if not self.periodic:
                    continue
                y[axis] = 0
            edges.append((x,tuple(y)))
        vertices = [x for edge in edges for x in edge]
        assert len(vertices) == len(set(vertices))
        return tuple(edges)

    def cone(self,roots,ticks):
        roots = frozenset(roots)
        if type(ticks) is not int or ticks < 0 or not roots <= set(self.sites):
            raise ValueError('finite supported roots and nonnegative tick count required')
        support = set(roots)
        for _ in range(ticks):
            for layer in reversed(self.layers):
                touched = [edge for edge in layer if set(edge) & support]
                support.update(x for edge in touched for x in edge)
        return frozenset(support)

    def one_charge_step(self,amplitudes,inverse=False):
        """Exact invariant triplet-color sector above the product background.

        U acts identically on each uniform native cell state in this sector.
        The spatial exchange is evaluated on the literal complete matchings.
        """
        if not set(amplitudes) <= set(self.sites):
            raise ValueError('excitation lies outside the finite box')
        if type(inverse) is not bool:
            raise ValueError('inverse flag must be boolean')
        for value in amplitudes.values():
            exact = s.sympify(value)
            if (isinstance(value,(float,complex)) or exact.has(s.Float)
                    or exact.is_number is not True or exact.is_finite is not True):
                raise ValueError('exact finite amplitudes required')
        state = dict(amplitudes)
        a,b = (s.conjugate(A),s.conjugate(B)) if inverse else (A,B)
        for layer in (reversed(self.layers) if inverse else self.layers):
            out = dict(state)
            for x,y in layer:
                left,right = state.get(x,0),state.get(y,0)
                out[x],out[y] = expanded(a*left+b*right),expanded(b*left+a*right)
            state = {x:value for x,value in out.items() if value != 0}
        return state
