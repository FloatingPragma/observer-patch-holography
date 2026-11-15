"""All-number CAR execution of the existing local Weyl/contact source words.

Field geometry and schedule resources do not allocate a global Fock matrix.
The explicit state evaluator has an independent, bounded host domain. It
retains all number sectors and the actual numerical local source matrices,
including vacuum scalars and tiny numerical off-sector residuals.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from ..interacting_field.lattice import SIGMA, real
from ..interacting_field import source
from ..common_source.chiral_candidate import gaussian_source as gaussian

MAX_MODES = 16
MAX_STATE_ENTRIES = 2**20


@dataclass(frozen=True)
class FieldGeometry:
    """Finite rectangular open box of the same two-component local action."""
    sides: tuple[int, int, int] = (2, 2, 1)
    speed: float = .125
    interaction: float = 1.

    def __post_init__(self):
        if type(self.sides) is not tuple or len(self.sides) != 3 or any(
                type(length) is not int or not 1 <= length <= 128 for length in self.sides):
            raise ValueError('three integer open-box side lengths in 1..128 required')
        if math.prod(self.sides) > 4096:
            raise ValueError('bounded geometry requires at most 4096 sites')
        object.__setattr__(self, 'speed', real(self.speed, 'coordinate speed', nonnegative=True))
        object.__setattr__(self, 'interaction', real(self.interaction, 'interaction', nonnegative=True))

    @property
    def sites(self):
        return tuple(product(*(range(length) for length in self.sides)))

    @property
    def volume(self):
        return math.prod(self.sides)

    @property
    def modes(self):
        return 2*self.volume

    @property
    def groups(self):
        indices = {site: index for index, site in enumerate(self.sites)}
        groups = [[] for _ in range(6)]
        for i, x in enumerate(self.sites):
            for axis in range(3):
                y = list(x); y[axis] += 1
                j = indices.get(tuple(y))
                if j is None:
                    continue
                for a in range(2):
                    for b in range(2):
                        z = 1j*self.speed*SIGMA[axis][a, b]/2
                        if z:
                            groups[2*axis+x[axis] % 2].append((2*i+a, 2*j+b, complex(z)))
        return tuple(tuple(group) for group in groups)

    @property
    def spatial_edges(self):
        return sum((self.sides[axis]-1)*math.prod(
            self.sides[j] for j in range(3) if j != axis) for axis in range(3))

    def single_h(self):
        rows, cols, values = [], [], []
        for group in self.groups:
            for p, q, z in group:
                rows.extend((p, q)); cols.extend((q, p)); values.extend((z, z.conjugate()))
        return sparse.coo_matrix((values, (rows, cols)),
                                 shape=(self.modes, self.modes), dtype=complex).tocsr()


class FullFock:
    """Bounded explicit complete Fock evaluator; no fixed-number projection."""
    def __init__(self, modes):
        if type(modes) is not int or not 2 <= modes <= MAX_MODES:
            raise ValueError(f'explicit full-Fock evaluator requires 2..{MAX_MODES} modes')
        self.modes, self.dimension = modes, 1 << modes
        self.addresses = np.arange(self.dimension, dtype=np.int64)

    def state(self, value):
        if np.ma.isMaskedArray(value):
            raise ValueError('masked state forbidden')
        state = np.asarray(value)
        if (state.ndim not in (1, 2) or state.shape[0] != self.dimension
                or (state.ndim == 2 and state.shape[1] == 0)
                or state.size > MAX_STATE_ENTRIES or state.dtype.kind not in 'ficu'
                or not np.isfinite(state).all()):
            raise ValueError('finite bounded full-Fock state/reference columns required')
        return np.asarray(state, dtype=complex)

    def _mode(self, mode):
        if type(mode) is not int or not 0 <= mode < self.modes:
            raise ValueError('mode outside the full Fock carrier')

    def annihilator(self, mode):
        self._mode(mode)
        inputs = self.addresses[(self.addresses & (1 << mode)) != 0]
        signs = np.array([(-1)**(int(i) & ((1 << mode)-1)).bit_count() for i in inputs])
        return sparse.coo_matrix((signs, (inputs ^ (1 << mode), inputs)),
                                 shape=(self.dimension, self.dimension), dtype=complex).tocsr()

    def number(self, mode):
        self._mode(mode)
        return sparse.diags(((self.addresses >> mode) & 1).astype(float), format='csr')

    def hamiltonian(self, geometry):
        if not isinstance(geometry, FieldGeometry) or geometry.modes != self.modes:
            raise ValueError('field geometry and Fock carrier differ')
        ladders = [self.annihilator(mode) for mode in range(self.modes)]
        result = sparse.csr_matrix((self.dimension, self.dimension), dtype=complex)
        for group in geometry.groups:
            for p, q, z in group:
                term = z*(ladders[p].conj().T@ladders[q])
                result += term+term.conj().T
        contact = sum(((self.addresses >> (2*x)) & 1)*
                      ((self.addresses >> (2*x+1)) & 1) for x in range(geometry.volume))
        return result+sparse.diags(geometry.interaction*contact, format='csr')

    def adjacent(self, state, left, gate):
        """Apply all sixteen local matrix entries, with no numerical projection."""
        self._mode(left)
        if left+1 >= self.modes:
            raise ValueError('adjacent gate requires a second mode')
        if np.ma.isMaskedArray(gate):
            raise ValueError('masked local gate forbidden')
        gate = np.asarray(gate)
        if gate.shape != (4, 4) or gate.dtype.kind not in 'ficu' or not np.isfinite(gate).all():
            raise ValueError('finite complete two-mode matrix required')
        out = self.state(state).copy()
        bit0, bit1 = 1 << left, 1 << (left+1)
        bases = self.addresses[(self.addresses & (bit0+bit1)) == 0]
        indices = np.array([bases, bases+bit0, bases+bit1, bases+bit0+bit1])
        out[indices] = np.einsum('ij,j...->i...', gate, out[indices])
        return out

    def routed(self, state, p, q, gate, *, source_routes=True):
        """Bring q next to p and undo every FSWAP after the complete local gate."""
        self._mode(p); self._mode(q)
        if p >= q or type(source_routes) is not bool:
            raise ValueError('ordered distinct modes and Boolean route choice required')
        swap = (gaussian.local_gate('fswap').fock_unitary if source_routes else
                np.array([[1,0,0,0], [0,0,1,0], [0,1,0,0], [0,0,0,-1]], complex))
        out = self.state(state).copy()
        route = tuple(range(q-1, p, -1))
        for left in route:
            out = self.adjacent(out, left, swap)
        out = self.adjacent(out, p, gate)
        for left in reversed(route):
            out = self.adjacent(out, left, swap)
        return out

    def exact_evolve(self, geometry, state, duration):
        duration = real(duration, 'Hamiltonian model duration')
        return expm_multiply(-1j*duration*self.hamiltonian(geometry), self.state(state))


def evolve(geometry, state, duration, steps, *, space=None):
    """Same symmetric hopping/contact source schedule on every number sector."""
    if not isinstance(geometry, FieldGeometry):
        raise ValueError('declared field geometry required')
    duration = real(duration, 'Hamiltonian model duration')
    if type(steps) is not int or steps < 1:
        raise ValueError('positive integer split count required')
    space = FullFock(geometry.modes) if space is None else space
    if not isinstance(space, FullFock) or space.modes != geometry.modes:
        raise ValueError('field geometry and Fock carrier differ')
    out = space.state(state).copy()
    half = duration/(2*steps)
    terms = tuple(term for group in geometry.groups for term in group)
    hops = {z: gaussian.local_gate('hopping', z, half).fock_unitary for _, _, z in terms}
    onsite, _ = source.contact(half, geometry.interaction)
    for _ in range(steps):
        for p, q, z in terms:
            out = space.routed(out, p, q, hops[z])
        for _ in range(2):
            for site in range(geometry.volume):
                out = space.adjacent(out, 2*site, onsite)
        for p, q, z in reversed(terms):
            out = space.routed(out, p, q, hops[z])
    return out


def number_pointer(space, state, mode):
    """Retain both source-operated occupation outcomes on the whole carrier."""
    if not isinstance(space, FullFock):
        raise ValueError('full Fock carrier required')
    space._mode(mode)
    kraus, _ = source.occupation()
    state = space.state(state)
    # A single-bit operator is applied literally as a tensor operation. The
    # intended number pointer is even; any numerical odd leakage is retained.
    empty = space.addresses[(space.addresses & (1 << mode)) == 0]
    indices = np.array([empty, empty+(1 << mode)])
    outputs = []
    for operator in kraus:
        out = np.zeros_like(state)
        out[indices] = np.einsum('ij,j...->i...', operator, state[indices])
        outputs.append(out)
    return tuple(outputs)


def resources(geometry, duration, steps, segments=1):
    """Same source count, independent of occupation number and reference size."""
    if not isinstance(geometry, FieldGeometry):
        raise ValueError('declared field geometry required')
    result = source.resources(geometry, duration, steps, segments)
    result['scope'] = 'all-number full-Fock evolution; preparation and pointers are separate'
    result.pop('number_sector_compression_residual_included_in_local_target_residuals')
    result['number_sector_projection_applied'] = False
    return result


def source_word(geometry, duration, steps, *, max_pulses=2**20):
    """Materialize a bounded actual source tape, including every routed word.

    This is optional: the general schedule and resource theorem do not need
    to materialize the tape or a global Fock operator.
    """
    if type(max_pulses) is not int or not 1 <= max_pulses <= 2**22:
        raise ValueError('bounded positive integer source-tape budget required')
    counts = resources(geometry, duration, steps)
    if counts['native_pulses'] > max_pulses:
        raise ValueError('source tape exceeds its declared host materialization budget')
    half = real(duration, 'Hamiltonian model duration')/(2*steps)
    swap = gaussian.local_gate('fswap').word
    terms = tuple(term for group in geometry.groups for term in group)
    hops = {z: gaussian.local_gate('hopping', z, half).word for _, _, z in terms}
    _, info = source.contact(half, geometry.interaction)
    tape = []
    def append(word, left):
        tape.extend({**event, 'targets': [left+target for target in event['targets']]}
                    for event in word)
    def hop(p, q, z):
        route = tuple(range(q-1, p, -1))
        for left in route:
            append(swap, left)
        append(hops[z], p)
        for left in reversed(route):
            append(swap, left)
    for _ in range(steps):
        for p, q, z in terms:
            hop(p, q, z)
        for _ in range(2):
            for site in range(geometry.volume):
                append(info['word'], 2*site)
        for p, q, z in reversed(terms):
            hop(p, q, z)
    if len(tape) != counts['native_pulses']:
        raise ArithmeticError('emitted all-number source tape differs from counted schedule')
    return tuple(tape)
