"""Exterior lift of the same link-dependent wave with a shared field register.

The full finite Fock space supplies canonical anticommutation. Executions use
sparse fixed-number vectors without truncating the underlying operator algebra.
The parity-controlled link extension is declared; it agrees with the parent on
one occupation. No spatial spin or physical species assignment is inferred.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from ..model import JointObserver

MODES = 7200


def mode(port, internal, parity, address):
    if (any(type(x) is not int for x in (port, internal, parity, address))
            or port not in range(12) or internal not in range(5)
            or parity not in (0, 1) or address not in range(60)):
        raise ValueError('mode outside the owned wave register')
    return ((port*5+internal)*2+parity)*60+address


def coordinates(index):
    if type(index) is not int or index not in range(MODES):
        raise ValueError('mode index outside the owned wave register')
    block, address = divmod(index, 60)
    block, parity = divmod(block, 2)
    port, internal = divmod(block, 5)
    return port, internal, parity, address


def wedge_insert(occupied, added):
    """Append one vector to an ordered wedge, retaining the permutation sign."""
    if added in occupied:
        return None, 0
    sign = -1 if sum(i > added for i in occupied) % 2 else 1
    return tuple(sorted((*occupied, added))), sign


@dataclass(frozen=True)
class FockState:
    amplitudes: tuple
    norm_denominator: int = 1

    def __post_init__(self):
        if type(self.norm_denominator) is not int or self.norm_denominator <= 0:
            raise ValueError('positive integer normalization required')
        items = self.amplitudes.items() if isinstance(self.amplitudes, dict) else self.amplitudes
        result = {}
        for (g, occupied), value in items:
            occupied = tuple(occupied)
            if (type(g) is not int or tuple(sorted(set(occupied))) != occupied
                    or any(type(i) is not int or i not in range(MODES) for i in occupied)):
                raise ValueError('an ordered occupation set is required')
            if type(value) is int:
                value = F(value)
            if not isinstance(value, F) and not (isinstance(value, s.Expr) and value.is_number
                    and value.is_finite and not value.has(s.Float)):
                raise ValueError('exact finite amplitudes required')
            key = g, occupied
            result[key] = result.get(key, F(0))+value
        object.__setattr__(self, 'amplitudes', tuple(sorted((k,v) for k,v in result.items() if v != 0)))

    def as_dict(self):
        return dict(self.amplitudes)

    def norm_squared(self):
        return s.simplify(sum((v.conjugate()*v for _,v in self.amplitudes), F(0))/self.norm_denominator)


class FockObserver:
    """One shared link register, with an exterior wave and conditional repair.

    The default rule applies T_e when endpoint occupation is odd. The explicit
    presence control applies it for any nonzero endpoint occupation. Both
    extend the parent one-wave law; their different pair outputs witness the
    additional architecture choice required by an occupation extension.
    """
    def __init__(self, control='parity'):
        if control not in ('parity', 'presence'):
            raise ValueError('unknown declared occupation control')
        self.base = JointObserver()
        self.control = control

    def _validate(self, state):
        if not isinstance(state, FockState):
            raise ValueError('Fock state required')
        if any(g not in self.base.link_words for (g, _), _ in state.amplitudes):
            raise ValueError('link word outside the shared field sector')

    def basis_state(self, g, occupied=()):
        result = FockState({(g, tuple(occupied)): F(1)})
        self._validate(result)
        return result

    @lru_cache(maxsize=8192)
    def column(self, g, index, inverse=False):
        p,j,r,a = coordinates(index)
        state = self.base.basis_state(g,p,j,r,a)
        result = self.base.inverse_wave(state) if inverse else self.base.wave_step(state)
        return tuple((mode(p,j,r,a),value) for (word,p,j,r,a),value in result.amplitudes if word == g)

    def wave(self, state, inverse=False):
        self._validate(state)
        output = {}
        for (g, occupied), coefficient in state.amplitudes:
            terms = {(): coefficient}
            for index in occupied:
                following = {}
                for partial, amplitude in terms.items():
                    for image, weight in self.column(g,index,inverse):
                        combined, sign = wedge_insert(partial,image)
                        if sign:
                            following[combined] = following.get(combined,F(0))+sign*amplitude*weight
                terms = {key:value for key,value in following.items() if value != 0}
            for wedge, value in terms.items():
                key = g,wedge
                output[key] = output.get(key,F(0))+value
        return FockState(output,state.norm_denominator)

    def controls(self, state, inverse=False):
        self._validate(state)
        output = {}
        edges = tuple(reversed(self.base.law.interior)) if inverse else self.base.law.interior
        for (g, occupied), value in state.amplitudes:
            word = g
            for edge in edges:
                endpoints = self.base.law.edges[edge]
                number = sum(coordinates(index)[0] in endpoints for index in occupied)
                enabled = number % 2 == 1 if self.control == 'parity' else number > 0
                if enabled:
                    word = self.base.law.transport(word,edge)[0]
            key = word,occupied
            output[key] = output.get(key,F(0))+value
        return FockState(output,state.norm_denominator)

    def step(self,state):
        return self.controls(self.wave(state))

    def inverse_step(self,state):
        return self.wave(self.controls(state,inverse=True),inverse=True)

    def gamma(self,state):
        self._validate(state)
        star = self.base.law.vertex_stars[5]
        return FockState({(g ^ star,occupied):
                          (-value if sum(coordinates(i)[0] == 5 and coordinates(i)[1] >= 3
                                         for i in occupied) % 2 else value)
                          for (g,occupied),value in state.amplitudes},state.norm_denominator)

    def face_measurement(self,state):
        self._validate(state)
        result = []
        for face in self.base.law.region:
            values = {key:value for key,value in state.amplitudes
                      if self.base.law.syndrome(key[0]) >> face & 1}
            result.append((face,FockState(values,state.norm_denominator)))
        return tuple(result)
