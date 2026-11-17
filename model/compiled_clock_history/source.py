"""Emitted source apparatus, its complete 512-address chart and exact bounds.

Only the three union bits are addressed by the apparatus. Six high bits
label spectator configurations in their total-occupation shell. Apparatus
targets preserve local number, so this fixed chart implements native-local
operations on all 66 admitted columns. All 446 other addresses remain in the
executed maps. Proper-code loading, export, blanks and pointer measurement
are supplied services; the returned binary arrays have no normalization or
unitary correction applied.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import hashlib
import math
from pathlib import Path

import numpy as np

from ..common_source import exact_output as exact
from ..observer import native
from ..observer.source_scalar import compiler
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


UNION = (0, 1, 5)
ROOT = Path(__file__).resolve().parents[2]
SPECTATORS = tuple(p for p in range(12) if p not in UNION)
STATES = native.states(binary=True)
KINDS = ('preparation', 'apparatus', 'readout')
SOURCE_PATHS = (
    'model/__init__.py',
    'model/common_source/__init__.py', 'model/common_source/exact_output.py',
    'model/common_source/experiment.py', 'model/common_source/gauge_matter.py',
    'model/common_source/infinite_history.py', 'model/common_source/protocol.py',
    'model/common_source/records.py', 'model/common_source/scalar_coupling.py',
    'model/common_source/smooth_history.py', 'model/common_source/smooth_protocol.py',
    'model/compiled_clock_history/__init__.py', 'model/compiled_clock_history/source.py',
    'model/observer/__init__.py', 'model/observer/native.py',
    'model/observer/scalar_response/__init__.py', 'model/observer/scalar_response/experiment.py',
    'model/observer/source_scalar/__init__.py', 'model/observer/source_scalar/bounds.py',
    'model/observer/source_scalar/compiler.py', 'model/observer/source_scalar/cutoff.py',
    'model/observer/source_scalar/experiment.py', 'model/observer/source_scalar/field.py',
    'model/anharmonic_transfer/__init__.py', 'model/anharmonic_transfer/environment.py',
    'model/common_source/native_join.py', 'model/observer/p_response/__init__.py',
    'model/observer/p_response/experiment.py', 'model/primitive_observer/__init__.py',
    'model/relativistic_scalar/__init__.py', 'model/relativistic_scalar/lattice.py',
    'model/scalar_observations/__init__.py', 'model/scalar_observations/lattice.py',
    'model/scalar_observations/source.py',
    'model/scalar_observations/experimental_readout/__init__.py',
    'model/scalar_observations/experimental_readout/canonical_source.py',
    'model/scalar_observations/experimental_readout/executable_source.py',
    'code/source_selection_model/response.json', 'code/m1_source_realization/receipt.json',
    'model/requirements.txt',
)


@lru_cache(maxsize=1)
def _shells():
    return tuple(tuple(bits for bits in product((0, 1), repeat=9) if sum(bits) == total)
                 for total in range(3))


def encode(state):
    """Fixed binary-total-two source chart; no input-dependent routing."""
    if (type(state) is not tuple or len(state) != 12
            or any(type(v) is not int or v not in (0, 1) for v in state) or sum(state) != 2):
        raise ValueError('an original twelve-port binary total-two tuple is required')
    local = sum(state[p] << bit for bit, p in enumerate(UNION))
    spectators = tuple(state[p] for p in SPECTATORS)
    return local+8*_shells()[sum(spectators)].index(spectators)


def decode(code):
    """Decode an admitted column, refusing every padding address."""
    if type(code) is not int or not 0 <= code < 512:
        raise ValueError('an original integer address from zero through511 is required')
    local, label = code % 8, code//8
    total = 2-local.bit_count()
    if total < 0 or label >= len(_shells()[total]):
        raise ValueError('the address is outside the binary total-two source code')
    values = [0]*12
    for bit, p in enumerate(UNION):
        values[p] = (local >> bit) & 1
    for p, bit in zip(SPECTATORS, _shells()[total][label], strict=True):
        values[p] = bit
    return tuple(values)


def codes():
    """Admitted addresses in the existing canonical66-state order."""
    return tuple(encode(state) for state in STATES)


def _parameters(kind, gamma, reference):
    if type(kind) is not str or kind not in KINDS:
        raise ValueError('a declared preparation, apparatus or readout word is required')
    if type(gamma) is not int or gamma not in (0, 1) or type(reference) is not bool:
        raise ValueError('an original binary interaction switch and Boolean reference are required')


def _exchange(angle):
    # These Pauli products commute. Their combined action fixes00/11 and
    # rotates the ordered10/01 pair by [[cos,-sin],[sin,cos]].
    return (compiler.pauli_rotation('IXY', -angle/2)
            + compiler.pauli_rotation('IYX', angle/2))


def word(kind='apparatus', *, gamma=1, reference=False):
    """A fresh chronological tape using only existing source primitives."""
    _parameters(kind, gamma, reference)
    if kind == 'preparation':
        tape = [compiler.pulse('Y', -math.pi/2, [0]),
                compiler.pulse('Y', -math.pi/2, [1])]+_exchange(math.pi/4)
        return tape+compiler.cnot(2, 3) if reference else tape
    if kind == 'readout':
        return compiler.cnot(1, 3)
    # J(pi)=diag(-1,-1,-1,1); retain the scalar I(pi) for exact CZ.
    phase = [compiler.pulse('J', gamma*math.pi, [0, 1]),
             compiler.pulse('I', gamma*math.pi, [0])]
    return phase+_exchange(math.atan2(4, 3))


@lru_cache(maxsize=12)
def _executed(kind, gamma, reference):
    tape = word(kind, gamma=gamma, reference=reference)
    if kind == 'preparation':
        blank = np.zeros(16 if reference else 8, complex)
        blank[0] = 1
        result = compiler.execute(tape, 4 if reference else 3, blank)
    elif kind == 'readout':
        blank_pointer = np.zeros((16, 8), complex)
        blank_pointer[:8] = np.eye(8)
        result = compiler.execute(tape, 4, blank_pointer)
    else:
        result = compiler.execute(tape, 3)
    return result.shape, result.dtype.str, result.tobytes()


def _array(kind, gamma=1, reference=False):
    _parameters(kind, gamma, reference)  # Original controls are checked before cache lookup.
    shape, dtype, payload = _executed(kind, gamma, reference)
    return np.frombuffer(payload, dtype=np.dtype(dtype)).reshape(shape)


def prepared(reference=False):
    """Raw blank-input ket, with little-endian index local+8*reference.

    Compilation coordinate3 is a separate blank reference, not spectator
    chart bit3. The preparation starts from an apparatus-register blank;
    this is not a claim of native-repair reachability from a zero-load state.
    """
    return _array('preparation', reference=reference)


def preparation(reference=False):
    """Lift the raw prepared ket to512 work rows and one/two reference columns."""
    local = prepared(reference)
    result = np.zeros((512, 2 if reference else 1), complex)
    result[:8] = local.reshape(-1, 8).T
    return np.frombuffer(result.tobytes(), dtype=result.dtype).reshape(result.shape)


def local_apparatus(gamma=1):
    return _array('apparatus', gamma=gamma)


def apparatus(gamma=1):
    """All512 raw addresses, including the entire declared padding extension."""
    result = np.kron(np.eye(64, dtype=int), local_apparatus(gamma))
    return np.frombuffer(result.tobytes(), dtype=result.dtype).reshape(512, 512)


def readout_isometry():
    """All16x8 components; temporary bit3 is a pointer, never a spectator bit."""
    return _array('readout')


def local_readout():
    matrix = readout_isometry()
    return matrix[:8], matrix[8:]


def target(kind='apparatus', *, gamma=1, reference=False):
    """Independent sparse exact target: each real entry is a+b*sqrt(2).

    Keys have the original vector/matrix rank. Fresh Fractions never expose
    cached certificate state. This target is not consumed by the executor.
    """
    _parameters(kind, gamma, reference)
    if kind == 'preparation':
        return {(3,): (F(0), F(1, 2)), (5+8*int(reference),): (F(0), F(1, 2))}
    if kind == 'readout':
        return {(column+8*((column >> 1) & 1), column): (F(1), F(0))
                for column in range(8)}
    matrix = {}
    for column in range(8):
        phase = (-1)**(gamma*(column & 1)*((column >> 1) & 1))
        left, right = (column >> 1) & 1, (column >> 2) & 1
        if left == right:
            matrix[column, column] = F(phase), F(0)
        else:
            matrix[column, column] = F(3*phase, 5), F(0)
            matrix[column ^ 6, column] = F((4 if left else -4)*phase, 5), F(0)
    return matrix


def certify(kind='apparatus', matrix=None, *, gamma=1, reference=False):
    """Exact all-component binary-array comparison, including references.

    For preparation the input is the specified blank. Otherwise the ideal
    target is an isometry on all8 inputs, including padding. If epsilon bounds
    the Frobenius error then epsilon*(2+epsilon)/2 bounds half diamond error
    of the single map; pointer dephasing preserves this bound for the complete
    readout instrument. The raw map norm is at most(1+epsilon)^2. This certifies
    the supplied array's values, not its provenance or physical gate error.
    """
    _parameters(kind, gamma, reference)
    shape = ((16 if reference else 8,) if kind == 'preparation'
             else (16, 8) if kind == 'readout' else (8, 8))
    if matrix is None:
        matrix = _array(kind, gamma, reference)
    if type(matrix) is not np.ndarray:
        raise ValueError('an original base ndarray is required before any array dispatch')
    matrix = exact._array(matrix, shape)
    wanted = target(kind, gamma=gamma, reference=reference)
    mass, norm, radical, cross = F(0), F(0), F(0), F(0)
    for index in np.ndindex(shape):
        real, imag = exact._component(matrix[index])
        a, b = wanted.get(index, (F(0), F(0)))
        mass += real*real+imag*imag
        norm += a*a+2*b*b
        radical += 2*a*b
        cross += real*a
    cross_radical = sum((exact._component(matrix[index])[0]*b
                         for index, (_, b) in wanted.items()), F(0))
    expected_norm = 1 if kind == 'preparation' else 8
    if norm != expected_norm or radical:
        raise ArithmeticError('the declared exact isometry norm changed')
    cross_lo, cross_hi = exact._radical_sum({F(1): cross, F(2): cross_radical})
    squared_lo, squared_hi = exact._nonnegative(mass+norm-2*cross_hi, mass+norm-2*cross_lo)
    epsilon = exact._sqrt_interval(squared_hi)[1]
    bounds = {'operator_error_upper': epsilon,
            'frobenius_error_lower': exact._sqrt_interval(squared_lo)[0],
            'frobenius_error_upper': epsilon,
            'half_diamond_upper': epsilon*(2+epsilon)/2,
            'channel_norm_upper': (1+epsilon)**2,
            'reported_squared_frobenius_norm': mass,
            'exact_target_squared_frobenius_norm': F(expected_norm),
            'fixed_blank_preparation_only': kind == 'preparation',
            'arbitrary_input_and_reference': kind != 'preparation',
            'full_encoded_work_dimension': 512,
            'all_original_components_checked': int(matrix.size)}
    # Directed-root helpers cache mutable Fraction objects. Their cached
    # entries must never escape through this public certificate dictionary.
    return {key: F(value.numerator, value.denominator) if type(value) is F else value
            for key, value in bounds.items()}


def resources(kind='apparatus', *, gamma=1, reference=False):
    """Count every emitted request, including the gamma-zero phase pulses."""
    counts = compiler.resources(word(kind, gamma=gamma, reference=reference))
    counts['blank_work_registers'] = int(kind == 'preparation')
    counts['blank_reference_qubits'] = int(kind == 'preparation' and reference)
    counts['blank_pointer_qubits'] = int(kind == 'readout')
    counts['pointer_reads'] = int(kind == 'readout')
    return counts


def bindings():
    """Fresh source-file digests; callers must also bind their history executor."""
    return bound_pins({path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS})
