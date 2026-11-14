"""Actual twelve-port source words for the fixed-band scalar observation code.

Every returned operation contracts emitted native source pulses.  Independent
targets only check those contractions.  The product schedule has a supplied
Hamiltonian model-time duration; native pulse counts do not calibrate clocks.
Blank pointers, code loading/export and routing retain their stated apparatus
scope.  Floating residuals do not certify analogue or physical service errors.
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
from math import pi
from numbers import Integral

import numpy as np
from scipy.linalg import expm

from ..observer.source_scalar import compiler as c
from ..relativistic_scalar.lattice import local_operators
from .lattice import GalerkinParameters, _mode, _real


COUNT_KEYS = (
    "native_pulses", "single_code_pulses", "pair_code_pulses",
    "identity_phase_pulses", "proper_code_load_events",
    "proper_code_export_events", "proper_code_service_events",
)


def _dimension(value):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
            or value not in (2, 4, 8)):
        raise ValueError("source occupation dimension must be two, four or eight")
    return int(value)


def _parameters(parameters):
    if parameters is None:
        return GalerkinParameters()
    if not isinstance(parameters, GalerkinParameters):
        raise ValueError("parameters must be a GalerkinParameters instance")
    return parameters


def _basis(matrix):
    values, vectors = np.linalg.eigh(matrix)
    # The first supported entry fixes signs without backend-dependent argmax ties.
    vectors[abs(vectors) < 1e-14] = 0
    for column in vectors.T:
        pivot = int(np.flatnonzero(abs(column) > 1e-12)[0])
        if column[pivot] < 0:
            column *= -1
    return values, vectors


def _counts(word):
    resources = c.resources(word)
    return {key: resources[key] for key in COUNT_KEYS}


def _execute(stages, qubits, target):
    word = [event for _, stage in stages for event in stage]
    raw = c.execute(word, qubits)
    residual = float(np.linalg.norm(raw-target, 2))
    if not np.isfinite(residual) or residual > 2e-10:
        raise ArithmeticError("full-phase source operation disagrees with the independent target")
    matrix, guard = c.stabilize_executed_unitary(raw)
    matrix.flags.writeable = False
    stage_counts = [{"stage": name, "resources": _counts(stage)} for name, stage in stages]
    resources = {key: sum(stage["resources"][key] for stage in stage_counts) for key in COUNT_KEYS}
    return matrix, {
        "resources": resources, "source_stages": stage_counts,
        "local_register_qubits": qubits, "full_phase_target_residual": residual,
        "numerical_stabilization": guard, "service_failure_outcomes_retained": True,
        "valid_code_failure_probability_in_exact_arithmetic": 0,
        "physical_service_error_certified": False,
    }


@lru_cache(maxsize=96)
def _onsite(dimension, mode, duration, parameters):
    operators = local_operators(dimension, parameters.omega)
    hamiltonian = (operators["p2"]/2
                   +parameters.frequency_squared(mode)*operators["q2"]/2
                   +parameters.onsite_quartic(mode)*operators["q4"])
    values, basis = _basis(hamiltonian)
    change = c.compile_real_orthogonal(basis)
    stages = [("onsite_inverse_basis", c.inverse(change)),
              ("onsite_eigenvalue_phase", c.compile_diagonal(-duration*values)),
              ("onsite_basis", change)]
    matrix, metadata = _execute(stages, dimension.bit_length()-1,
                                expm(-1j*duration*hamiltonian))
    metadata["mode"] = mode
    metadata["compressed_after_forming_oscillator_polynomials"] = True
    return matrix, metadata


def onsite(dimension, mode, duration, parameters=None):
    """Execute exp(-i duration h_mode) with its complete scalar phase."""
    matrix, metadata = _onsite(_dimension(dimension), _mode(mode),
                               _real(duration, "model duration"), _parameters(parameters))
    return matrix.copy(), deepcopy(metadata)


@lru_cache(maxsize=64)
def _pair(dleft, dright, coefficient, duration, omega):
    left = local_operators(dleft, omega)["q2"]
    right = local_operators(dright, omega)["q2"]
    left_values, left_basis = _basis(left)
    right_values, right_basis = _basis(right)
    left_bits, right_bits = dleft.bit_length()-1, dright.bit_length()-1
    left_change = c.compile_real_orthogonal(left_basis)
    right_change = c.remap(c.compile_real_orthogonal(right_basis),
                           list(range(left_bits, left_bits+right_bits)))
    phases = -duration*coefficient*np.kron(right_values, left_values)
    stages = [("pair_inverse_left_basis", c.inverse(left_change)),
              ("pair_inverse_right_basis", c.inverse(right_change)),
              ("pair_product_phase", c.compile_diagonal(phases)),
              ("pair_right_basis", right_change),
              ("pair_left_basis", left_change)]
    matrix, metadata = _execute(stages, left_bits+right_bits,
                                expm(-1j*duration*coefficient*np.kron(right, left)))
    metadata["generator"] = "coefficient*(P_left Q^2 P_left)*(P_right Q^2 P_right)"
    return matrix, metadata


def pair(dleft, dright, coefficient, duration, omega=2.):
    """Execute the positive q2*q2 interaction in separate local eigenbases.

    Index n_left+dleft*n_right orders the two local codes.  Squared full
    oscillator coordinates are compressed before this tensor product.
    """
    matrix, metadata = _pair(_dimension(dleft), _dimension(dright),
                             _real(coefficient, "pair coefficient", nonnegative=True),
                             _real(duration, "model duration"),
                             _real(omega, "basis frequency", positive=True))
    return matrix.copy(), deepcopy(metadata)


@lru_cache(maxsize=32)
def _instrument(dimension, eta, theta, omega):
    q = local_operators(dimension, omega)["q"]
    values, basis = _basis(q)
    angles = (eta*values+theta)/2
    cosine = (basis*np.cos(angles))@basis.T
    sine = (basis*np.sin(angles))@basis.T
    pointer = dimension.bit_length()-1
    change = c.compile_real_orthogonal(basis)
    # exp(i*pi*X/4) Z exp(-i*pi*X/4) = Y.  Tape order is chronological.
    # The resulting pointer rotation exp(-i*angle*Y) is [[cos,-sin],[sin,cos]].
    stages = [("ramsey_inverse_field_basis", c.inverse(change)),
              ("ramsey_pointer_inverse_axis", [c.pulse("X", -pi/4, [pointer])]),
              ("ramsey_conditional_phase", c.compile_diagonal(np.r_[-angles, angles])),
              ("ramsey_pointer_axis", [c.pulse("X", pi/4, [pointer])]),
              ("ramsey_field_basis", change)]
    target = np.block([[cosine, -sine], [sine, cosine]])
    matrix, metadata = _execute(stages, pointer+1, target)
    kraus = matrix[:dimension, :dimension], matrix[dimension:, :dimension]
    metadata.update({
        "instrument": "K0=cos((eta*Q_d+theta)/2); K1=sin((eta*Q_d+theta)/2)",
        "eta": eta, "theta": theta,
        "complete_pointer_outcomes_retained": True,
        "unrestricted_instrument_compression_claimed": False,
        "completeness_residual": float(np.linalg.norm(
            sum(k.conj().T@k for k in kraus)-np.eye(dimension), 2)),
    })
    return kraus, metadata


def instrument(dimension, eta=1., theta=pi/2, omega=2.):
    """Return both disturbed branches of the real cos/sin Ramsey dilation.

    Q_d is the compressed coordinate.  This declares the finite-code
    detector directly; f(Q_d) is not asserted equal to P f(Q) P.
    """
    kraus, metadata = _instrument(_dimension(dimension), _real(eta, "detector strength"),
                                  _real(theta, "detector phase"),
                                  _real(omega, "basis frequency", positive=True))
    return tuple(matrix.copy() for matrix in kraus), deepcopy(metadata)


@lru_cache(maxsize=32)
def _feedback(dimension, kick, omega):
    q = local_operators(dimension, omega)["q"]
    values, basis = _basis(q)
    change = c.compile_real_orthogonal(basis)
    stages = [("feedback_inverse_field_basis", c.inverse(change)),
              ("feedback_momentum_kick", c.compile_diagonal(kick*values)),
              ("feedback_field_basis", change)]
    matrix, metadata = _execute(stages, dimension.bit_length()-1, expm(1j*kick*q))
    metadata["generator"] = "exp(i*kick*Q_d)"
    metadata["kick"] = kick
    return matrix, metadata


def feedback(dimension, kick=.25, omega=2.):
    """Execute a record-selected momentum kick, with the signed phase retained."""
    matrix, metadata = _feedback(_dimension(dimension), _real(kick, "momentum kick"),
                                 _real(omega, "basis frequency", positive=True))
    return matrix.copy(), deepcopy(metadata)
