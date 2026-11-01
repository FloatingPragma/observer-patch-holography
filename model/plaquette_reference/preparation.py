"""A same-action Gibbs column ensemble, with every classical/export outcome.

No matrix is a physical primitive: polynomial columns determine angles for
the existing native amplitude compiler.  The exact ideal circuit and Taylor
error theorem are separate from the admitted floating source execution.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math

import numpy as np
from scipy import sparse

from ..field_stress.energy import GaussianInteger as GI
from ..field_symmetry.symmetry import operators
from ..full_fock_field.field import FieldGeometry
from ..observer.source_scalar import compiler
from ..reference_recovery.recovery import execute_complete
from ..reference_observation.observation import frozen


GEOMETRY = FieldGeometry((2, 2, 1))
MODES = 8
DIMENSION = 256
BETA = F(1, 32)
ORDER = 8
EXECUTION_ADMISSION = F(1, 50_000_000_000)


def error_theorem():
    """Rational bounds; no spectral or floating evaluation proves these.

    ||H-I|| <= B=3, q=beta*B/2.  T=exp[-beta(H-I)/2],
    P=sum_{k=0}^8 [-beta(H-I)/2]^k/k!.  The exponential tail is
    bounded by q^9/(9! (1-q)).  Normalized vectorizations differ by
    at most 2*r/(1-q), hence their density trace distance is <=4*r/(1-q).
    A CPTP regional restriction cannot increase that trace distance.
    """
    bound = F(3)
    q = BETA*bound/2
    remainder = q**(ORDER+1)/(math.factorial(ORDER+1)*(1-q))
    density_error = 4*remainder/(1-q)
    minimum_polynomial = (1-q-remainder)**2/(16*(1/(1-q)+remainder)**2)
    minimum_gibbs = (1-4*q)/16
    minimum = min(minimum_polynomial, minimum_gibbs)
    # The resolvent integral for log gives ||log A-log B|| <= ||A-B||/m.
    return {'H_centered_norm_bound': bound, 'q': q,
            'exponential_operator_error': remainder,
            'density_trace_error': density_error,
            'regional_minimum_weight': minimum,
            'regional_modular_operator_error': density_error/minimum,
            'full_minimum_weight': minimum/16,
            'full_modular_operator_error': 16*density_error/minimum}


@dataclass(frozen=True)
class Design:
    numerator: tuple
    denominator: int
    column_norm_squared: tuple
    total_norm_squared: int
    frame: np.ndarray
    hamiltonian: np.ndarray

    def probabilities(self):
        return np.array([x/self.total_norm_squared for x in self.column_norm_squared])

    def column(self, index):
        if type(index) is not int or not 0 <= index < DIMENSION:
            raise ValueError('one of all 256 complete occupation columns required')
        denominator = math.sqrt(self.column_norm_squared[index])
        return np.array([row[index]/denominator for row in self.numerator])

    def target_columns(self):
        real = np.array(self.numerator, dtype=float)/math.sqrt(self.total_norm_squared)
        return self.frame.conj().T@real

    def certificate(self):
        digest = hashlib.sha256()
        entries = 0
        for row in self.numerator:
            for value in row:
                digest.update((str(value)+';').encode())
                entries += value != 0
        return {'polynomial_denominator': str(self.denominator),
                'integer_numerator_sha256': digest.hexdigest(),
                'nonzero_entries': entries,
                'positive_column_norms': all(x > 0 for x in self.column_norm_squared),
                'total_integer_norm_squared': str(self.total_norm_squared),
                'probability_numerators': [str(x) for x in self.column_norm_squared]}


@lru_cache(maxsize=1)
def design():
    """Derive exact rational P8 from the unchanged full-Fock Weyl/contact H."""
    source = operators(GEOMETRY)
    identity = GI.from_exact(sparse.eye(DIMENSION))
    centered = source.hubbard_H-GI.from_exact(source.scale*sparse.eye(DIMENSION))
    if centered.imag.nnz or not (centered-centered.adjoint()).is_zero():
        raise ArithmeticError('the exact onsite frame must expose a real symmetric action')
    if not (source.frame@source.H@source.frame.adjoint()-source.hubbard_H).is_zero():
        raise ArithmeticError('the real frame is not the same original Weyl/contact action')
    row_bound = int(np.max(abs(centered.real).sum(axis=1)))
    if F(row_bound, source.scale) > error_theorem()['H_centered_norm_bound']:
        raise ArithmeticError('the exact full-Fock row norm exceeds the Taylor premise')
    step_denominator = 64*source.scale
    denominator = step_denominator**ORDER*math.factorial(ORDER)
    numerator = [[0]*DIMENSION for _ in range(DIMENSION)]
    power = identity
    for k in range(ORDER+1):
        factor = (-1)**k*step_denominator**(ORDER-k)*(math.factorial(ORDER)//math.factorial(k))
        entries = power.real.tocoo()
        for i, j, value in zip(entries.row, entries.col, entries.data, strict=True):
            numerator[int(i)][int(j)] += factor*int(value)
        if k < ORDER:
            power = power@centered
    norms = tuple(sum(numerator[i][j]**2 for i in range(DIMENSION))
                  for j in range(DIMENSION))
    if not all(norms):
        raise ArithmeticError('a complete column disappeared from the faithful reference')
    return Design(tuple(tuple(row) for row in numerator), denominator, norms, sum(norms),
                  frozen(source.frame.matrix().toarray()),
                  frozen(source.H.matrix(source.scale).toarray()))


def onsite_word(site):
    """Second-quantize R=(i sigma_x)^x(i sigma_y)^y(i sigma_z)^z.

    Labels are in little-endian up/down order. Every angle and basis change
    is compiled to the same six-level native source; even/empty occupations
    are included, with determinant-one action on the filled occupation.
    """
    if type(site) is not tuple or len(site) != 3 or any(type(x) is not int or x < 0 for x in site):
        raise ValueError('three nonnegative integer site coordinates required')
    blocks = ((('ZI', -math.pi/4), ('IZ', math.pi/4)),
              (('XY', math.pi/4), ('YX', -math.pi/4)),
              (('XX', math.pi/4), ('YY', math.pi/4)))
    word = []
    for coordinate, block in zip(reversed(site), blocks, strict=True):
        for _ in range(coordinate % 4):
            for label, angle in block:
                word += compiler.pauli_rotation(label, angle)
    return word


def inverse_frame_word():
    word = []
    for x, site in enumerate(GEOMETRY.sites):
        word += compiler.remap(onsite_word(site), [2*x, 2*x+1])
    return compiler.inverse(word)


def probability_word(data=None):
    data = design() if data is None else data
    return compiler.compile_real_state(np.sqrt(data.probabilities()))


def conditional_word(index, data=None):
    data = design() if data is None else data
    if type(index) is not int or not 0 <= index < DIMENSION:
        raise ValueError('one of all 256 retained source records required')
    # A reversible record-conditioned blanking of |j>, with its phase paid.
    word = [compiler.pulse('X', math.pi/2, [bit]) for bit in range(MODES) if index >> bit & 1]
    word += [compiler.pulse('I', -math.pi*index.bit_count()/2, [0])]
    return word+compiler.compile_real_state(data.column(index))+inverse_frame_word()


@dataclass(frozen=True)
class Readout:
    initial: np.ndarray
    columns: np.ndarray
    failure_weights: tuple
    events: int
    native_kraus: tuple

    def failure(self, event):
        """Reconstruct every actual rectangular abort vector at one read."""
        if type(event) is not int or not 0 <= event < self.events:
            raise ValueError('an existing readout event is required')
        branches = {0: self.initial.reshape(DIMENSION, 1)}
        current = 0
        for bit in range(MODES):
            following = {}
            for prefix, column in branches.items():
                if current == event:
                    arranged = np.moveaxis(column.reshape([2]*MODES+[1]),
                                           MODES-1-bit, 0).reshape(2, -1)
                    return tuple(frozen(k[2:, :2]@arranged) for k in self.native_kraus[2:])
                current += 1
                for value in (0, 1):
                    following[prefix | (value << bit)] = compiler._apply(
                        column, self.native_kraus[value][:2, :2], [bit], MODES)
            branches = following
        raise ArithmeticError('unreachable retained readout event')


def record_readout(state):
    """Execute eight complete local dephasing exports, retaining all records.

    A coherent load embeds the chosen bit in the two-dimensional proper
    code. The existing incoherent M6 export has two success outcomes and
    four failure outcomes. All 255 adaptive read events are evaluated; the
    four failures at each event are retained, including their exact zeros.
    """
    if np.ma.isMaskedArray(state):
        raise ValueError('masked source states are not admitted')
    raw = np.asarray(state)
    if raw.shape != (DIMENSION,) or raw.dtype.kind not in 'iufc' or not np.isfinite(raw).all():
        raise ValueError('one finite complete eight-mode source vector required')
    kraus = tuple(frozen(k) for k in compiler.transfer_kraus(2, coherent=False))
    if np.linalg.norm(sum(k.conj().T@k for k in kraus)-np.eye(6), 2) > 1e-14:
        raise ArithmeticError('readout export lost an outcome')
    branches = {0: np.asarray(raw, complex).reshape(DIMENSION, 1)}
    failures = []
    events = 0
    for bit in range(MODES):
        following = {}
        for prefix, column in branches.items():
            events += 1
            for value in (0, 1):
                gate = kraus[value][:2, :2]
                following[prefix | (value << bit)] = compiler._apply(column, gate, [bit], MODES)
            # Rectangular abort maps preserve the other seven spectator bits.
            arranged = column.reshape([2]*MODES+[1])
            arranged = np.moveaxis(arranged, MODES-1-bit, 0).reshape(2, -1)
            weights = tuple(float(np.vdot(k[2:, :2]@arranged, k[2:, :2]@arranged).real)
                            for k in kraus[2:])
            failures.append(weights)
        branches = following
    columns = np.column_stack([branches[j][:, 0] for j in range(DIMENSION)])
    total = np.vdot(columns, columns).real+math.fsum(map(sum, failures))
    if abs(total-np.vdot(raw, raw).real) > 2e-12:
        raise ArithmeticError('complete input readout probability was not retained')
    return Readout(frozen(raw), frozen(columns), tuple(failures), events, kraus)


@dataclass(frozen=True)
class Prepared:
    probability_run: object
    readout: Readout
    conditional_runs: tuple
    columns: np.ndarray
    failure_probability: float
    complete_output_error: float

    @property
    def field_density(self):
        """Unnormalized success field; its abort complement remains separate."""
        return self.columns@self.columns.conj().T


def prepare():
    """Execute every branch of the complete 256-record preparation instrument."""
    data = design()
    blank = np.zeros(DIMENSION); blank[0] = 1
    initial = execute_complete(probability_word(data), MODES, blank)
    readout = record_readout(initial.success[:, 0])
    runs = tuple(execute_complete(conditional_word(j, data), MODES, readout.columns[:, j])
                 for j in range(DIMENSION))
    columns = np.column_stack([run.success[:, 0] for run in runs])
    failures = (initial.failure_probability+math.fsum(map(sum, readout.failure_weights))
                +math.fsum(run.failure_probability for run in runs))
    error = math.sqrt(float(np.linalg.norm(columns-data.target_columns()))**2+failures)
    if not math.isfinite(error) or error > float(EXECUTION_ADMISSION):
        raise ArithmeticError('complete executed source ensemble exceeds its numerical admission')
    if abs(np.vdot(columns, columns).real+failures-1) > 2e-11:
        raise ArithmeticError('complete preparation failed to account for an input outcome')
    return Prepared(initial, readout, runs, frozen(columns), failures, error)
