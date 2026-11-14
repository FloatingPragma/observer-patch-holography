"""Executable finite observer with coherent seam records and optional probes.

The program uses exact rational amplitudes, with a common square-root
normalization carried separately. It executes the one-load finite construction;
the general local double-load completion is in common.local_repair_unitary.
"""
from dataclasses import dataclass
from fractions import Fraction as F
import sympy as s
from model.common_source import slow_geometry

PROGRAM_ID = 'coherent-seam-repair-v1'


@dataclass(frozen=True)
class State:
    amplitudes: tuple
    norm_denominator: int = 1

    def __post_init__(self):
        values = []
        for value in self.amplitudes:
            if type(value) is int:
                value = F(value)
            if not isinstance(value, F) and not (isinstance(value, s.Expr) and
                    value.is_number and value.is_finite and not value.has(s.Float)):
                raise ValueError('amplitudes must be exact finite numbers')
            values.append(value)
        object.__setattr__(self, 'amplitudes', tuple(values))

    def norm_squared(self):
        return s.simplify(sum(x.conjugate() * x for x in self.amplitudes) / self.norm_denominator)


class Observer:
    """One explicit program; fibre=5 admits a declared static central link.

    Idle events are coherent step applications. Quarter-turns and register
    measurements are optional apparatus interventions, with separate records.
    """

    def __init__(self, fibre=1, central_edge=None):
        if type(fibre) is not int or fibre not in (1, 5):
            raise ValueError('fibre must be 1 or 5')
        self.seams = tuple(tuple(edge) for edge in slow_geometry.source().seams())
        self.addresses = tuple(pair for edge in self.seams for pair in (edge, edge[::-1]))
        if central_edge is not None:
            central_edge = tuple(sorted(central_edge))
            if fibre != 5 or central_edge not in self.seams:
                raise ValueError('central link requires the five-level fibre and a source seam')
        self.fibre, self.central_edge = fibre, central_edge
        self.address_count = len(self.addresses)
        self.dimension = 12 * fibre * 2 * self.address_count

    def index(self, port, parity, address, internal=0):
        values = (port, parity, address, internal)
        if any(type(x) is not int for x in values):
            raise ValueError('register addresses must be integers')
        if not (0 <= port < 12 and parity in (0, 1) and 0 <= address < self.address_count and 0 <= internal < self.fibre):
            raise ValueError('register address outside the declared space')
        return ((port * self.fibre + internal) * 2 + parity) * self.address_count + address

    def _validate(self, state):
        if not isinstance(state, State) or len(state.amplitudes) != self.dimension:
            raise ValueError('state does not belong to this observer')
        if type(state.norm_denominator) is not int or state.norm_denominator <= 0:
            raise ValueError('invalid amplitude normalization')

    def embedded_port(self, port, internal=0):
        """A blank parity record and uniform coherent seam address."""
        result = [F(0)] * self.dimension
        for address in range(self.address_count):
            result[self.index(port, 0, address, internal)] = F(1)
        return State(tuple(result), self.address_count)

    def basis_state(self, port, parity, address, internal=0):
        """An owned computational-register preparation.

        The uniform ensemble of these dimension-many states has reduced
        density I/d. Reading its retained preparation label conditions it.
        """
        result = [F(0)] * self.dimension
        result[self.index(port, parity, address, internal)] = F(1)
        return State(tuple(result))

    def reflect(self, state):
        """The declared address reflection R=2AA*−I."""
        self._validate(state)
        x, y, count = state.amplitudes, list(state.amplitudes), self.address_count
        for block in range(24 * self.fibre):
            start = block * count
            mean = sum(x[start:start+count], F(0)) / count
            for k in range(count):
                y[start+k] = 2 * mean - x[start+k]
        return State(tuple(y), state.norm_denominator)

    def repair(self, state):
        """Controlled reversible seam correction Q, including its owned memory."""
        self._validate(state)
        x, y = state.amplitudes, list(state.amplitudes)
        for address, (a, b) in enumerate(self.addresses):
            for internal in range(self.fibre):
                link = -1 if self.central_edge == tuple(sorted((a, b))) and internal >= 3 else 1
                sites = [self.index(a, 0, address, internal), self.index(b, 0, address, internal),
                         self.index(a, 1, address, internal), self.index(b, 1, address, internal)]
                signs = [1, -link, -1, -link]
                coefficient = sum(sign * x[i] for sign, i in zip(signs, sites)) / 2
                for sign, i in zip(signs, sites):
                    y[i] -= sign * coefficient
        return State(tuple(y), state.norm_denominator)

    def step(self, state):
        return self.repair(self.reflect(state))

    def inverse_step(self, state):
        return self.reflect(self.repair(state))

    def quarter_turn(self, state, address):
        """Optional basis-selective apparatus phase; absent from idle evolution."""
        self._validate(state)
        if type(address) is not int or not 0 <= address < self.dimension:
            raise ValueError('probe address outside the declared space')
        y = list(state.amplitudes)
        y[address] *= s.I
        return State(tuple(y), state.norm_denominator), {'kind': 'quarter_turn', 'address': address}

    def register_measurement(self, state, address):
        """Both unnormalized branches of {P_address,I−P_address}, retaining labels."""
        self._validate(state)
        if type(address) is not int or not 0 <= address < self.dimension:
            raise ValueError('detector address outside the declared space')
        yes, no = [F(0)] * self.dimension, list(state.amplitudes)
        yes[address], no[address] = state.amplitudes[address], F(0)
        return tuple(({'kind': 'register_readout', 'address': address, 'outcome': label},
                      State(tuple(values), state.norm_denominator))
                     for label, values in ((True, yes), (False, no)))

    def preparation_filter(self, state):
        """Complete heralded filter, K+=(2I+W)/3 and K−=sqrt(2)(I−W)/3."""
        self._validate(state)
        evolved = self.step(state).amplitudes
        success = tuple((2*x+y)/3 for x, y in zip(state.amplitudes, evolved))
        failure = tuple(s.sqrt(2)*(x-y)/3 for x, y in zip(state.amplitudes, evolved))
        return (({'kind': 'preparation', 'outcome': 'success'}, State(success, state.norm_denominator)),
                ({'kind': 'preparation', 'outcome': 'failure'}, State(failure, state.norm_denominator)))

    def source_populations(self, state):
        self._validate(state)
        width = self.fibre * 2 * self.address_count
        return tuple(s.simplify(sum(z.conjugate()*z for z in state.amplitudes[a*width:(a+1)*width]) / state.norm_denominator)
                     for a in range(12))

    def embedded_return_amplitude(self, state, port, internal=0):
        """Overlap with the normalized blank-record preparation A|port,internal>."""
        self._validate(state)
        numerator = sum(state.amplitudes[self.index(port, 0, k, internal)] for k in range(self.address_count))
        return s.simplify(numerator / s.sqrt(self.address_count * state.norm_denominator))

    def embedded_return_measurement(self, state, port, internal=0):
        """Complete return detector for the specified blank-record preparation.

        The projector is onto A|port,internal>. The complement is retained;
        this operation belongs to the optional apparatus, not idle repair.
        """
        self._validate(state)
        sites = [self.index(port, 0, k, internal) for k in range(self.address_count)]
        mean = sum((state.amplitudes[i] for i in sites), F(0)) / self.address_count
        yes, no = [F(0)] * self.dimension, list(state.amplitudes)
        for i in sites:
            yes[i] = mean
            no[i] -= mean
        return tuple(({'kind': 'source_return', 'port': port, 'internal': internal, 'outcome': outcome},
                      State(tuple(values), state.norm_denominator))
                     for outcome, values in ((True, yes), (False, no)))
