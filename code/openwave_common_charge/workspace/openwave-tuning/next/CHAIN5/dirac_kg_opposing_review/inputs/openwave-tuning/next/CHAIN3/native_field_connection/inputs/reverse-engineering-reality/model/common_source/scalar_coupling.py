"""One finite register for scalar oscillators, charged CAR matter and a link.

The scalar action and the neutral scalar/matter coupling are declared model
choices. Every evolution and preparation routine is executed with the existing
twelve-port M6 pulse compiler. Target matrices are independent comparison
objects; they are never substituted for executed source words.
"""
from __future__ import annotations

from functools import lru_cache
from fractions import Fraction as F
import math
from numbers import Integral, Real

import numpy as np
from scipy.linalg import expm

from ..observer.source_scalar import compiler as c
from ..observer.source_scalar import field
from . import gauge_matter as gauge


COUPLING = 1 / 4


def _dimension(d):
    if isinstance(d, (bool, np.bool_)) or not isinstance(d, Integral) or d not in (2, 4):
        raise ValueError("the executed common model has onsite dimension 2 or 4")
    return int(d)


def _real(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be a finite real number")
    return value


def _parameters(d, steps, time, coupling):
    d = _dimension(d)
    if isinstance(steps, (bool, np.bool_)) or not isinstance(steps, Integral) or steps < 1:
        raise ValueError("a positive integer number of model-time steps required")
    time, coupling = _real(time, "model time"), _real(coupling, "coupling")
    if time < 0:
        raise ValueError("nonnegative model time required")
    return d, int(steps), time, coupling


def qubits(d):
    return 2 + 2 * (_dimension(d).bit_length() - 1)


def _spectral_word(values, basis, coefficient):
    word = c.compile_real_orthogonal(basis)
    return c.inverse(word) + c.compile_diagonal(coefficient * values) + word


def initialpreparation_word(d):
    """Blank gauge bits plus the actually compiled projected scalar vacuum."""
    d = _dimension(d)
    vacuum, _ = field.vacuum_vector(d)
    return c.remap(c.compile_real_state(vacuum), list(range(2, qubits(d))))


def _retained_resources(word):
    result = c.resources(word)
    # Eigenbasis and elimination roundoff can turn an analytically empty
    # angle into a tiny nonzero number on another LAPACK platform. Keep
    # every emitted operation and its pulse/service count, but do not pin
    # this numerical zero diagnostic in replayed evidence.
    result.pop("zero_angle_requests", None)
    return result


@lru_cache(maxsize=2)
def _prepared(d):
    word = initialpreparation_word(d)
    # Contract the scalar routine on every input column before applying the
    # existing guarded polar correction; do not consult a target matrix.
    unitary, guard = c.stabilize_executed_unitary(c.execute(word, qubits(d)))
    state = unitary[:, 0]
    state.flags.writeable = False
    metadata = {"resources": _retained_resources(word), "guards": {"preparation": guard},
                "register_qubits": qubits(d), "gauge_initial_bits": [0, 0],
                "scalar_preparation": field.vacuum_vector(d)[1],
                "physical_service_error_certified": False}
    return state, metadata


def preparation(d):
    """Return the native-source prepared vector and its resource/rounding data."""
    from copy import deepcopy
    state, metadata = _prepared(_dimension(d))
    return state.copy(), deepcopy(metadata)


def scalar_terms(d, coupling=COUPLING):
    """Independent matrices for onsite, bond and scalar/charge interaction."""
    d, coupling = _dimension(d), _real(coupling, "coupling")
    q = field.quadratures(d)[0]
    onsite = field.onsite_hamiltonian(d)
    identity = np.eye(d)
    # Bit zero denotes right occupation, so n_L0-n_R0 is Z_0.
    difference = np.diag([1., -1., 1., -1.])
    return {
        "onsite": np.kron(np.kron(identity, onsite) + np.kron(onsite, identity), np.eye(4)),
        "bond": -float(field.KAPPA) * np.kron(np.kron(q, q), np.eye(4)),
        "interaction": coupling * math.sqrt(2) * np.kron(np.kron(identity, q), difference),
    }


def scalar_subwords(d, dt, coupling=COUPLING):
    """Exact finite subterms, each compiled as native source pulses.

    The two onsite terms commute. Their half-step is followed by a half
    bond step and a full scalar/charge step; the caller mirrors the halves.
    Scalar sites occupy higher bits, and the charge word includes bit zero.
    """
    d, _, dt, coupling = _parameters(d, 1, dt, coupling)
    bits = d.bit_length() - 1
    q_values, q_basis = field.canonical_eigenbasis(field.quadratures(d)[0])
    h_values, h_basis = field.canonical_eigenbasis(field.onsite_hamiltonian(d))
    onsite = _spectral_word(h_values, h_basis, -dt / 2)
    q_word = c.compile_real_orthogonal(q_basis)
    both_q = q_word + c.remap(q_word, list(range(bits, 2 * bits)))
    bond = (c.inverse(both_q)
            + c.compile_diagonal(dt * float(field.KAPPA) * np.outer(q_values, q_values).ravel() / 2)
            + both_q)
    # Local bit zero is the gauge coordinate; scalar Q occupies higher bits.
    q_with_charge = c.remap(q_word, list(range(1, bits + 1)))
    phases = -dt * coupling * math.sqrt(2) * np.kron(q_values, np.array([1., -1.]))
    interaction = c.inverse(q_with_charge) + c.compile_diagonal(phases) + q_with_charge
    return {
        "onsite_half": c.remap(onsite, list(range(2, 2 + bits)))
                       + c.remap(onsite, list(range(2 + bits, 2 + 2 * bits))),
        "bond_half": c.remap(bond, list(range(2, 2 + 2 * bits))),
        "interaction": c.remap(interaction, [0, *range(2, 2 + bits)]),
    }


def step_word(d, dt, coupling=COUPLING):
    """Chronological symmetric word for one step of the common Hamiltonian."""
    d, _, dt, coupling = _parameters(d, 1, dt, coupling)
    terms = scalar_subwords(d, dt, coupling)
    gauge_half = gauge.source_step_word(dt / 2)
    return (gauge_half + terms["onsite_half"] + terms["bond_half"]
            + terms["interaction"] + terms["bond_half"]
            + terms["onsite_half"] + gauge_half)


def hamiltonian_terms(d, coupling=COUPLING):
    """Finite target terms in the same scalar tensor gauge basis."""
    d, coupling = _dimension(d), _real(coupling, "coupling")
    return {
        "hopping": np.kron(np.eye(d * d), gauge.numerical_matrix(gauge.hopping_hamiltonian())),
        "electric": np.kron(np.eye(d * d), gauge.numerical_matrix(gauge.electric_energy())),
        **scalar_terms(d, coupling),
    }


def ideal_hamiltonian(d, coupling=COUPLING):
    """H_hop + E² + Pi H_scalar Pi + g sqrt(2) Q0 (n_L0-n_R0)."""
    return sum(hamiltonian_terms(d, coupling).values())


def time_error_bound(d, steps, time=F(1), coupling=F(1, 4)):
    """Rational all-input bound for ideal finite product evolution.

    Both the palindrome and exp(-itH) match through degree two. Their
    third derivative norms are at most S^3, so integral Taylor remainder
    and unitary telescoping give ||U_n-exp(-itH)|| <= S^3 t^3/(3 n^2).
    Q_d^2 <= Pi Q^2 Pi <= 2N+1 supplies an explicit finite norm budget.
    This does not certify the floating contraction or an occupation limit.
    """
    d, steps, _, _ = _parameters(d, steps, 0, 0)
    if type(time) not in (int, F) or type(coupling) not in (int, F) or time < 0:
        raise ValueError("the certificate requires exact rational time and coupling")
    t, g = F(time), abs(F(coupling))
    root = field._root_upper(F(4 * d - 2))
    s = 3 + (1 + 2 * field.KAPPA) * (2 * d - 1) + g * root
    error = s**3 * t**3 / (3 * steps**2)
    return {
        "method": "palindromic degree-two agreement, integral Taylor remainder and unitary telescoping",
        "sum_generator_norms_upper": str(s),
        "unclipped_operator_error_upper": str(error),
        "ideal_split_operator_error_upper": str(min(F(2), error)),
        "ideal_channel_half_diamond_error_upper": str(min(F(1), error)),
        "model_time_exact": str(t), "coupling_exact": str(F(coupling)),
        "spatial_or_occupation_error_included": False,
        "floating_execution_or_physical_error_included": False,
    }


def independent_split_unitary(d, steps, time=1., coupling=COUPLING):
    """Compare using target-matrix exponentials, never source pulse matrices."""
    d, steps, time, coupling = _parameters(d, steps, time, coupling)
    dt = time / steps
    terms = hamiltonian_terms(d, coupling)
    hopping_quarter = expm(-.25j * dt * terms["hopping"])
    gauge_half = hopping_quarter @ expm(-.5j * dt * terms["electric"]) @ hopping_quarter
    onsite_half = expm(-.5j * dt * terms["onsite"])
    bond_half = expm(-.5j * dt * terms["bond"])
    step = (gauge_half @ onsite_half @ bond_half @ expm(-1j * dt * terms["interaction"])
            @ bond_half @ onsite_half @ gauge_half)
    return np.linalg.matrix_power(step, steps)


@lru_cache(maxsize=40)
def _source_unitary(d, steps, time, coupling):
    word = step_word(d, time / steps, coupling)
    raw = c.execute(word, qubits(d))
    step, guard = c.stabilize_executed_unitary(raw)
    result = np.linalg.matrix_power(step, steps)
    result.flags.writeable = False
    per_step = _retained_resources(word)
    additive = {"native_pulses", "single_code_pulses", "pair_code_pulses",
                "identity_phase_pulses", "proper_code_load_events",
                "proper_code_export_events", "proper_code_service_events",
                "integrated_processor_norm_upper"}
    resources = {key: value * steps if key in additive else value
                 for key, value in per_step.items()}
    metadata = {
        "resources": resources,
        "per_step_resources": per_step,
        "guards": {"executed_step": guard},
        "steps": steps, "model_time": time, "coupling": coupling,
        "register_qubits": qubits(d),
        "final_unitarity_residual": float(np.linalg.norm(result.conj().T @ result - np.eye(len(result)), 2)),
        "coupling_operator": "sqrt(2) Q0 (n_L0-n_R0)",
        "execution": "native source step contracted on every input column, then repeated",
        "spatial_or_occupation_continuum_claim": False,
        "physical_duration_or_service_error_certified": False,
        # The ideal certificate treats these actual finite binary inputs as
        # exact rationals. Angle evaluation and matrix roundoff are excluded.
        "finite_time_certificate": time_error_bound(d, steps, F(time), F(coupling)),
    }
    return result, metadata


def source_unitary(d, steps, time=1., coupling=COUPLING):
    """Full source-executed common-system unitary and counted word metadata."""
    from copy import deepcopy
    d, steps, time, coupling = _parameters(d, steps, time, coupling)
    result, metadata = _source_unitary(d, steps, time, coupling)
    return result.copy(), deepcopy(metadata)


def source_evolve(d, steps, time=1., coupling=COUPLING):
    """Prepare and evolve one common state, without an inserted probability."""
    initial, prep_metadata = preparation(d)
    evolution, metadata = source_unitary(d, steps, time, coupling)
    return evolution @ initial, {"preparation": prep_metadata, "evolution": metadata}
