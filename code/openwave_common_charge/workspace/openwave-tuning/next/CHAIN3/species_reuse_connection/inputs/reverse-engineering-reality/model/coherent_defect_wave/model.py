"""One exact finite clock joining coherent seam waves and movable link defects.

The idle clock has no measured moved/held flag: such a flag would dephase U.
Optional detectors return every unnormalized branch, including zero branches.
"""
from collections import defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction as F
from collections.abc import Mapping
import sympy as s

from model.coherent_seam_repair.protection import BoundaryFlux

PROGRAM_ID = 'coherent-defect-wave-v1'


@dataclass(frozen=True)
class SparseState:
    """Exact coefficient/sqrt(norm_denominator), keyed by (g,p,j,r,a)."""
    amplitudes: tuple
    norm_denominator: int = 1

    def __post_init__(self):
        if type(self.norm_denominator) is not int or self.norm_denominator <= 0:
            raise ValueError('normalization denominator must be a positive integer')
        items = self.amplitudes.items() if isinstance(self.amplitudes, Mapping) else self.amplitudes
        out = {}
        for key, value in items:
            key = tuple(key)
            if len(key) != 5 or any(type(x) is not int for x in key):
                raise ValueError('basis key must contain five integer registers')
            if type(value) is int:
                value = F(value)
            if not isinstance(value, F) and not (isinstance(value, s.Expr) and value.is_number
                    and value.is_finite and not value.has(s.Float)):
                raise ValueError('amplitudes must be exact finite numbers')
            out[key] = out.get(key, F(0)) + value
        object.__setattr__(self, 'amplitudes', tuple(sorted((k, v) for k, v in out.items() if v != 0)))

    def as_dict(self):
        return dict(self.amplitudes)

    def norm_squared(self):
        return s.simplify(sum((v.conjugate()*v for _, v in self.amplitudes), F(0)) / self.norm_denominator)


def combine(left, right, a=1, b=1):
    """Linear combination with one unchanged common normalization."""
    if left.norm_denominator != right.norm_denominator:
        raise ValueError('linear combination requires a common normalization')
    result = {k: a*v for k, v in left.amplitudes}
    for key, value in right.amplitudes:
        result[key] = result.get(key, F(0)) + b*value
    return SparseState(result, left.norm_denominator)


class JointObserver:
    """The fixed 72000-dimensional one-defect/one-wave finite program.

    Tick order: B (all-link-dependent coherent wave), then the five interior
    occupation-controlled transports in sorted carrier-seam order.
    """
    def __init__(self):
        self.law = BoundaryFlux()
        self.addresses = tuple(pair for edge in self.law.edges for pair in (edge, edge[::-1]))
        self.initial_word = 1 << self.law.edges.index((0, 1))
        reached, queue = {self.initial_word}, deque([self.initial_word])
        while queue:
            g = queue.popleft()
            for edge in self.law.interior:
                h = self.law.transport(g, edge)[0]
                if h not in reached:
                    reached.add(h)
                    queue.append(h)
        self.link_words = tuple(sorted(reached))
        self.dimension = len(self.link_words)*12*5*2*len(self.addresses)
        self._transports = {(g, e): self.law.transport(g, e)[0]
                            for g in self.link_words for e in self.law.interior}
        if self.dimension != 72000:
            raise ValueError('the pinned one-defect domain changed')

    def _key(self, g, port, internal, parity, address):
        key = (g, port, internal, parity, address)
        if (any(type(x) is not int for x in key) or g not in self.link_words
                or port not in range(12) or internal not in range(5)
                or parity not in (0, 1) or address not in range(60)):
            raise ValueError('basis register outside the declared joint domain')
        return key

    def _validate(self, state):
        if not isinstance(state, SparseState):
            raise ValueError('expected a joint sparse state')
        for key, _ in state.amplitudes:
            self._key(*key)

    def basis_state(self, g, port, internal=3, parity=0, address=0):
        return SparseState({self._key(g, port, internal, parity, address): F(1)})

    def embedded_port(self, g, port, internal=3):
        return self.embedded_ports(g, (port,), internal)

    def embedded_ports(self, g, ports, internal=3):
        ports = tuple(ports)
        if not ports or len(set(ports)) != len(ports):
            raise ValueError('preparation requires distinct nonempty ports')
        return SparseState({self._key(g, p, internal, 0, a): F(1)
                            for p in ports for a in range(60)}, 60*len(ports))

    def reflect(self, state):
        self._validate(state)
        blocks = defaultdict(dict)
        for (g, p, j, r, a), value in state.amplitudes:
            blocks[(g, p, j, r)][a] = value
        output = {}
        for block, values in blocks.items():
            twice_mean = sum(values.values(), F(0))/30
            for a in range(60):
                output[(*block, a)] = twice_mean-values.get(a, F(0))
        return SparseState(output, state.norm_denominator)

    def repair(self, state):
        """Q_g uses every central-link sign in the actual link-memory word."""
        self._validate(state)
        x = state.as_dict()
        output = dict(x)
        active = {(g, j, a) for (g, p, j, r, a), _ in state.amplitudes
                  if p in self.addresses[a]}
        for g, j, address in active:
            a, b = self.addresses[address]
            link = -1 if j >= 3 and (g >> (address//2) & 1) else 1
            keys = [(g, a, j, 0, address), (g, b, j, 0, address),
                    (g, a, j, 1, address), (g, b, j, 1, address)]
            signs = (1, -link, -1, -link)
            coefficient = sum((sign*x.get(key, F(0)) for key, sign in zip(keys, signs)), F(0))/2
            for key, sign in zip(keys, signs):
                output[key] = output.get(key, F(0))-sign*coefficient
        return SparseState(output, state.norm_denominator)

    def wave_step(self, state):
        return self.repair(self.reflect(state))

    def inverse_wave(self, state):
        return self.reflect(self.repair(state))

    def defect_control(self, state, edge):
        self._validate(state)
        if type(edge) is not int or edge not in self.law.interior:
            raise ValueError('controlled transport requires an interior edge')
        endpoints = self.law.edges[edge]
        return SparseState({(self._transports[(g, edge)] if p in endpoints else g, p, j, r, a): value
                            for (g, p, j, r, a), value in state.amplitudes}, state.norm_denominator)

    def controls(self, state, inverse=False):
        edges = reversed(self.law.interior) if inverse else self.law.interior
        for edge in edges:
            state = self.defect_control(state, edge)
        return state

    def step(self, state):
        return self.controls(self.wave_step(state))

    def inverse_step(self, state):
        return self.inverse_wave(self.controls(state, inverse=True))

    def preparation_filter(self, state):
        evolved = self.step(state)
        success = combine(state, evolved, F(2, 3), F(1, 3))
        failure = combine(state, evolved, s.sqrt(2)/3, -s.sqrt(2)/3)
        return (({'kind': 'joint_preparation', 'outcome': 'success'}, success),
                ({'kind': 'joint_preparation', 'outcome': 'failure'}, failure))

    def field_measurement(self, state):
        self._validate(state)
        return tuple(({'kind': 'link_word_readout', 'outcome': g},
                      SparseState({key: value for key, value in state.amplitudes if key[0] == g},
                                  state.norm_denominator)) for g in self.link_words)

    def region_face_measurement(self, state):
        """Five complete gauge-invariant defect-position outcomes in the region."""
        self._validate(state)
        return tuple(({'kind': 'regional_flux_readout', 'face_index': face,
                       'face': list(self.law.faces[face])},
                      SparseState({key: value for key, value in state.amplitudes
                                   if self.law.syndrome(key[0]) >> face & 1},
                                  state.norm_denominator)) for face in self.law.region)

    def wave_basis_measurement(self, state, port, internal=3, parity=0, address=0):
        self._validate(state)
        self._key(self.initial_word, port, internal, parity, address)
        coordinate = (port, internal, parity, address)
        return tuple(({'kind': 'wave_register_readout', 'coordinate': list(coordinate), 'outcome': yes},
                      SparseState({key: value for key, value in state.amplitudes
                                   if (key[1:] == coordinate) == yes}, state.norm_denominator))
                     for yes in (True, False))

    def wave_preparation_measurement(self, state, ports, internal=3):
        """Complete return to an owned uniform-address preparation, ignoring g."""
        self._validate(state)
        ports = tuple(ports)
        self.embedded_ports(self.initial_word, ports, internal)
        sums = defaultdict(lambda: F(0))
        for (g, p, j, r, a), value in state.amplitudes:
            if p in ports and j == internal and r == 0:
                sums[g] += value
        yes = SparseState({(g, p, internal, 0, a): value/(60*len(ports))
                           for g, value in sums.items() for p in ports for a in range(60)},
                          state.norm_denominator)
        no = combine(state, yes, 1, -1)
        return (({'kind': 'wave_preparation_return', 'outcome': True}, yes),
                ({'kind': 'wave_preparation_return', 'outcome': False}, no))
