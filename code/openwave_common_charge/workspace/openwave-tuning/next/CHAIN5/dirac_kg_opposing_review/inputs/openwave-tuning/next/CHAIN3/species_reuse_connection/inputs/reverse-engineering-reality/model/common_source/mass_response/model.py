"""Explicit coherent continuations with coefficients from the native reference.

The reference-to-hopping/quartic assignment is a NEW architecture choice.
It is neither the unread native C/E channel nor an already identified physical
Hamiltonian. No measured mass or coupling is an input. The four-mode ring is
the executed cell; the Hubbard thermodynamic theorem is a separate law.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import math

import numpy as np
from scipy.linalg import eigh, expm

from .. import gauge_matter
from ...native_modular_reference import local_preparation
from ...observer import native
from ...observer.source_scalar import compiler


RING_EDGES = ((0, 1), (2, 3), (1, 2), (0, 3))


@dataclass(frozen=True)
class Parameters:
    weights: tuple[Fraction, Fraction, Fraction]
    intra: float
    inter: float
    interaction: float
    clock: float

    def record(self):
        return {
            "native_weights_per_configuration": list(map(str, self.weights)),
            "intra": self.intra, "inter": self.inter,
            "interaction": self.interaction, "clock_frequency": self.clock,
            "coefficient_map": {
                "intra": "log(w1/w0)", "inter": "log(w2/w1)",
                "interaction": "log(w1*w1/(w0*w2))", "clock": "log(w2/w0)",
            },
            "architecture_map_declared": True,
            "reference_is_exact_ideal_native_preparation_target": True,
            "finite_native_warmup_executed_in_this_history": False,
            "finite_reference_preparation_error_propagated": False,
            "empirical_mass_inputs": [], "empirical_coupling_inputs": [],
        }


@lru_cache(maxsize=1)
def parameters():
    """Recover the weights through the existing native preparation producer."""
    rows = tuple(zip(native.states(binary=True), local_preparation.reference(), strict=True))
    levels = []
    for q in range(3):
        values = {w for x, w in rows if sum(x[i] for i in local_preparation.CORE) == q}
        if len(values) != 1:
            raise ArithmeticError("the prepared reference lost its core-level symmetry")
        levels.append(values.pop())
    w0, w1, w2 = levels
    if not 0 < w0 < w1 < w2:
        raise ArithmeticError("positive ordered reference levels required")
    return Parameters(tuple(levels), math.log(float(w1/w0)), math.log(float(w2/w1)),
                      math.log(float(w1*w1/(w0*w2))), math.log(float(w2/w0)))


def _real(value, name, *, nonnegative=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, Fraction)):
        raise ValueError(f"finite real {name} required")
    value = float(value)
    if not math.isfinite(value) or (nonnegative and value < 0):
        raise ValueError(f"finite {'nonnegative ' if nonnegative else ''}{name} required")
    return value


@lru_cache(maxsize=1)
def annihilators():
    """Reuse the existing four-mode CAR occupation-sign implementation."""
    out = []
    for mode in range(4):
        matrix = np.zeros((16, 16))
        for occupation in range(16):
            image = gauge_matter.car_ladder(occupation, mode, False)
            if image is not None:
                sign, target = image
                matrix[target, occupation] = sign
        matrix.flags.writeable = False
        out.append(matrix)
    return tuple(out)


def number_diagonal():
    return np.array([mask.bit_count() for mask in range(16)], float)


def ring_terms(*, interaction=True, chemical_potential=0.0):
    """Independent occupation matrices; these are NOT used to execute words."""
    if type(interaction) is not bool:
        raise ValueError("Boolean interaction control required")
    mu = _real(chemical_potential, "chemical potential")
    p = parameters()
    ladders = annihilators()
    occupation = [a.T @ a for a in ladders]
    terms = {}
    for index, (i, j) in enumerate(RING_EDGES):
        strength = p.intra if index < 2 else p.inter
        terms[f"bond_{i}_{j}"] = -strength * (ladders[i].T @ ladders[j] + ladders[j].T @ ladders[i])
    identity = np.eye(16)
    terms["interaction"] = (p.interaction if interaction else 0.) * sum(
        (occupation[i] - identity/2) @ (occupation[j] - identity/2) for i, j in RING_EDGES)
    terms["chemical_potential"] = mu * np.diag(number_diagonal())
    return terms


def ring_hamiltonian(*, interaction=True, chemical_potential=0.0):
    return sum(ring_terms(interaction=interaction, chemical_potential=chemical_potential).values())


def sector_ground(hamiltonian, number):
    if type(number) is not int or not 0 <= number <= 4:
        raise ValueError("occupation sector 0 through 4 required")
    indices = [mask for mask in range(16) if mask.bit_count() == number]
    values, vectors = eigh(np.asarray(hamiltonian)[np.ix_(indices, indices)])
    vector = np.zeros(16, complex)
    vector[indices] = vectors[:, 0]
    vector *= np.exp(-1j*np.angle(vector[np.argmax(abs(vector))]))
    return float(values[0]), np.real_if_close(vector)


def charge_transfer_preparation(*, interaction=True):
    """Two identical cells, fixed TOTAL charge four, no odd-sector coherence.

    The two arms are (N_A,N_B)=(2,2) and (3,1). Both reservoirs evolve with
    the same law, so a common chemical potential cancels from their frequency.
    This is a declared known-state preparation, not a reader of unknown states.
    """
    h = ring_hamiltonian(interaction=interaction)
    e1, v1 = sector_ground(h, 1)
    e2, v2 = sector_ground(h, 2)
    e3, v3 = sector_ground(h, 3)
    reference_arm = np.kron(v2, v2)
    charged_arm = np.kron(v1, v3)  # cell A occupies the low four bits
    state = (reference_arm + charged_arm) / math.sqrt(2)
    return state, {
        "reference_arm_energy": 2*e2, "charged_arm_energy": e1+e3,
        "pair_gap": e1+e3-2*e2,
        "single_gap_under_particle_hole_symmetry": (e1+e3-2*e2)/2,
        "total_charge": 4,
        "reference_arm": "N_A=2,N_B=2", "charged_arm": "N_A=3,N_B=1",
        "state_preparation_is_an_architecture_service": True,
    }


def _bond_word(i, j, angle):
    i, j = sorted((i, j))
    words = []
    for axis in ("X", "Y"):
        label = ["I"] * 4
        label[i] = label[j] = axis
        for k in range(i+1, j):
            label[k] = "Z"
        words += compiler.pauli_rotation(label, angle)
    return words


def cell_word(duration, steps, *, interaction=True):
    """Local hopping and quartic Strang words; no spectral gap is an input."""
    duration = _real(duration, "model interval", nonnegative=True)
    if type(steps) is not int or not 1 <= steps <= 4096:
        raise ValueError("split count from 1 through 4096 required")
    if type(interaction) is not bool:
        raise ValueError("Boolean interaction control required")
    p = parameters()
    dt = duration / steps
    intra = sum((_bond_word(i, j, dt*p.intra/4) for i, j in RING_EDGES[:2]), [])
    inter = sum((_bond_word(i, j, dt*p.inter/4) for i, j in RING_EDGES[2:]), [])
    density = []
    if interaction:
        for i, j in RING_EDGES:
            label = ["I"] * 4
            label[i] = label[j] = "Z"
            density += compiler.pauli_rotation(label, -dt*p.interaction/4)
    step = intra + inter + density + inter + intra
    return step * steps


@lru_cache(maxsize=64)
def _source_cell_unitary(duration, steps, interaction):
    """Contract the actual local source word on ALL 16 input columns."""
    word = cell_word(duration, steps, interaction=interaction)
    raw = compiler.execute(word, 4)
    matrix = np.frombuffer(raw.tobytes(), dtype=raw.dtype).reshape(raw.shape)
    return matrix, compiler.resources(word)


def source_cell_unitary(duration, steps, interaction=True):
    """Immutable computed source matrix and detached resource metadata."""
    duration = _real(duration, "model interval", nonnegative=True)
    if type(steps) is not int or not 1 <= steps <= 4096 or type(interaction) is not bool:
        raise ValueError("positive bounded split count and Boolean interaction required")
    matrix, cost = _source_cell_unitary(duration, steps, interaction)
    return matrix, dict(cost)


def apply_subsystem(state, matrix, targets, qubits):
    """A tensor contraction of a computed matrix, preserving reference columns."""
    state = np.asarray(state, complex)
    if state.ndim != 2 or state.shape[0] != 1 << qubits:
        raise ValueError("work-by-reference matrix required")
    targets = tuple(targets)
    if len(set(targets)) != len(targets) or any(type(j) is not int or not 0 <= j < qubits for j in targets):
        raise ValueError("distinct supported subsystem bits required")
    matrix = np.asarray(matrix, complex)
    if matrix.shape != (1 << len(targets),) * 2 or not np.all(np.isfinite(matrix)):
        raise ValueError("finite subsystem matrix with matching dimensions required")
    # targets[0] is the subsystem's least significant bit.
    axes = [qubits-1-j for j in reversed(targets)]
    axes += [j for j in range(qubits) if j not in axes]
    axes += [qubits]
    arranged = state.reshape([2]*qubits + [state.shape[1]]).transpose(axes)
    output = (matrix @ arranged.reshape(len(matrix), -1)).reshape(arranged.shape)
    return output.transpose(np.argsort(axes)).reshape(state.shape)


def evolve(state, duration, steps, *, interaction=True, backend="source"):
    """The two cells AND both clocks advance in the same declared interval.

    The neutral clocks are programmed finite modular-reference clocks. No
    identification with laboratory seconds or an autonomous hardware clock is
    inferred. Caching matrix contractions does not reduce physical word costs.
    """
    duration = _real(duration, "model interval", nonnegative=True)
    if backend == "source":
        unitary, resource = source_cell_unitary(duration, steps, interaction)
    elif backend == "continuous":
        unitary = expm(-1j*duration*ring_hamiltonian(interaction=interaction))
        resource = {"native_pulses": 0, "proper_code_service_events": 0}
    else:
        raise ValueError("source or independent continuous backend required")
    output = apply_subsystem(state, unitary, range(4), 10)
    output = apply_subsystem(output, unitary, range(4, 8), 10)
    clock_word = [compiler.pulse("Z", -parameters().clock*duration/2, [bit]) for bit in (8, 9)]
    if backend == "source":
        output = compiler.execute(clock_word, 10, output)
    else:
        phase = parameters().clock*duration/2
        clock = np.diag([np.exp(-1j*phase), np.exp(1j*phase)])
        for bit in (8, 9):
            output = apply_subsystem(output, clock, [bit], 10)
    return output, {"native_pulses": 2*resource["native_pulses"] + len(clock_word),
                    "proper_code_service_events": 2*resource["proper_code_service_events"] + 2*len(clock_word)}


def split_operator_bound(duration, steps):
    """Exact rational ideal bound for BOTH cells (clock commutes exactly).

    The term norm sum is2a+2b+J=log(3125/64)<4. The exact positive Taylor
    partial sum through degree7 at4 already exceeds3125/64. The palindrome
    and exponential match to degree two, hence2*4^3*t^3/(3*n^2). This bound
    excludes source floating contraction, apparatus services and physical use.
    """
    if type(duration) not in (int, Fraction) or duration < 0 or type(steps) is not int or steps < 1:
        raise ValueError("exact nonnegative duration and positive split count required")
    w0, w1, w2 = parameters().weights
    ratio = (w2/w0)**2 * w1*w1/(w0*w2)
    exp_four_lower = sum((Fraction(4)**k/math.factorial(k) for k in range(8)), Fraction())
    if not 1 < ratio < exp_four_lower:
        raise ArithmeticError("exact source coefficient norm bound changed")
    return min(Fraction(2), 2 * 4**3 * Fraction(duration)**3 / (3*steps**2))


def ideal_clock_calibration():
    """Faithful reference on the q=0/q=2 pair, with modular gap log5."""
    w0, _, w2 = parameters().weights
    p0, p1 = w0/(w0+w2), w2/(w0+w2)
    return {"populations": [str(p0), str(p1)], "contrast": float(p0-p1),
            "generator": "K=constant*I+log(5)*Z/2",
            "parameter": "programmed modular parameter s, not physical seconds",
            "clock_services_supplied": True}
