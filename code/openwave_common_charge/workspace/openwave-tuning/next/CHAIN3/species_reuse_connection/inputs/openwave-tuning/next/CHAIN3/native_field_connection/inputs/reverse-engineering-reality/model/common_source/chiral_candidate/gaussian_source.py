"""Native local Fock words, lifted through the finite Gaussian CAR identity.

Only four-column local Fock matrices are executed. Their completed Gaussian
action determines the many-mode Majorana transformation; no 2**M vector or
target BdG exponential is executed. Nonadjacent fermions are brought together
and restored with actual source FSWAP words, including their |11> minus sign.
Their entire routing cost is retained. The quadratic architecture, processor
code services, blank states, control precision and pointer readout are supplied.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
import math

import numpy as np

from ..mass_response import model as car_model
from ...observer.source_scalar import compiler


MAX_MODES = 512
TOLERANCE = 2e-10
_ADDITIVE_RESOURCES = (
    "native_pulses", "single_code_pulses", "pair_code_pulses", "zero_angle_requests",
    "identity_phase_pulses", "proper_code_load_events", "proper_code_export_events",
    "proper_code_service_events", "integrated_processor_norm_upper")


def _real(value, name):
    if isinstance(value, (bool, np.bool_)) or isinstance(value, (str, bytes)) or np.ndim(value) != 0:
        raise ValueError(f"finite real {name} required")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"finite real {name} required") from exc
    if not math.isfinite(result):
        raise ValueError(f"finite real {name} required")
    return result


def _complex(value, name):
    if isinstance(value, (bool, np.bool_)) or isinstance(value, (str, bytes)) or np.ndim(value) != 0:
        raise ValueError(f"finite complex {name} required")
    try:
        result = complex(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"finite complex {name} required") from exc
    if not math.isfinite(result.real) or not math.isfinite(result.imag) or not math.isfinite(abs(result)):
        raise ValueError(f"finite complex {name} required")
    return result


def _freeze(value):
    raw = np.asarray(value)
    return np.frombuffer(raw.tobytes(), dtype=raw.dtype).reshape(raw.shape)


def _matrix_record(value):
    value = np.asarray(value)
    return {"real": value.real.tolist(), "imag": value.imag.tolist()}


def _word_copy(word):
    return tuple({**event, "targets": list(event["targets"])} for event in word)


def _phase(number, angle):
    """exp(i angle*n_number), including its scalar phase."""
    return [compiler.pulse("I", angle/2, [number]),
            compiler.pulse("Z", -angle/2, [number])]


def _real_pair_word(kind, angle):
    # Native Pauli words implement exp(i theta P).
    return (compiler.pauli_rotation("XX", -angle/2)
            + compiler.pauli_rotation("YY", angle/2 if kind == "pairing" else -angle/2))


def local_word(kind, coefficient=0., duration=0.):
    """Chronological two-mode native source word for exp(-i duration H).

    hopping: H=z c0* c1+z* c1* c0;
    pairing: H=d c0* c1*+d* c1 c0;
    onsite: H=e n0 (mode1 is an identity spectator).
    fswap has fixed parameters0,0 and implements the fermion mode swap.
    """
    if type(kind) is not str or kind not in ("hopping", "pairing", "onsite", "fswap"):
        raise ValueError("hopping, pairing, onsite or fswap block required")
    z, time = _complex(coefficient, "coefficient"), _real(duration, "model duration")
    if kind == "fswap":
        if z != 0 or time != 0:
            raise ValueError("the FSWAP block has fixed zero request parameters")
        return (_real_pair_word("hopping", math.pi/2)
                + _phase(0, math.pi/2)+_phase(1, math.pi/2))
    if kind == "onsite":
        if z.imag != 0:
            raise ValueError("a real onsite coefficient required")
        return _phase(0, -time*z.real)
    angle = math.atan2(z.imag, z.real)
    # D=exp(i arg(z)*n0): H(z)=D H(abs(z)) D*. Execute D* first.
    return _phase(0, -angle)+_real_pair_word(kind, time*abs(z))+_phase(0, angle)


@lru_cache(maxsize=1)
def local_majoranas():
    # Existing CAR signs, restricted to the two low modes' complete Fock block.
    ladders = [a[:4, :4] for a in car_model.annihilators()[:2]]
    return tuple(_freeze(w) for a in ladders for w in (a+a.conj().T, -1j*(a-a.conj().T)))


@dataclass(frozen=True)
class LocalGate:
    kind: str
    coefficient: complex
    duration: float
    word: tuple[dict, ...]
    fock_unitary: np.ndarray
    orthogonal: np.ndarray
    resources: dict
    diagnostics: dict

    def record(self):
        return {"kind": self.kind, "coefficient": [self.coefficient.real, self.coefficient.imag],
                "duration": self.duration, "word": list(_word_copy(self.word)),
                "source_fock_unitary": _matrix_record(self.fock_unitary),
                "source_majorana_transform": self.orthogonal.tolist(),
                "resources": dict(self.resources), "diagnostics": dict(self.diagnostics)}


@lru_cache(maxsize=256)
def _local(kind, real, imag, duration):
    coefficient = complex(real, imag)
    word = local_word(kind, coefficient, duration)
    # Contract ALL four Fock input columns through the existing M6 source.
    unitary = compiler.execute(word, 2)
    majoranas = local_majoranas()
    evolved = [unitary.conj().T@w@unitary for w in majoranas]
    coefficients = np.array([[np.trace(w@image)/4 for w in majoranas] for image in evolved])
    orthogonal = coefficients.real
    linearity = max(float(np.max(abs(image-sum(value*w for value, w in zip(row, majoranas, strict=True)))))
                    for image, row in zip(evolved, orthogonal, strict=True))
    diagnostics = {"source_unitarity_residual": float(np.max(abs(unitary.conj().T@unitary-np.eye(4)))),
                   "majorana_linearity_residual": linearity,
                   "majorana_coefficient_imaginary_residual": float(np.max(abs(coefficients.imag))),
                   "majorana_orthogonality_residual": float(np.max(abs(orthogonal@orthogonal.T-np.eye(4)))),
                   "determinant": float(np.linalg.det(orthogonal)),
                   "contracted_fock_columns": 4, "polar_or_norm_correction_applied": False}
    if max(diagnostics[key] for key in ("source_unitarity_residual", "majorana_linearity_residual",
                                        "majorana_coefficient_imaginary_residual")) > TOLERANCE:
        raise ArithmeticError("the completed source block is not a resolved Gaussian CAR map")
    return LocalGate(kind, coefficient, duration, _word_copy(word), _freeze(unitary),
                     _freeze(orthogonal), compiler.resources(word), diagnostics)


def local_gate(kind, coefficient=0., duration=0.):
    """Detached word/resources and immutable actual source matrices."""
    z, time = _complex(coefficient, "coefficient"), _real(duration, "model duration")
    # Validate kind and the fixed onsite/FSWAP domains before consulting cache.
    local_word(kind, z, time)
    gate = _local(kind, z.real, z.imag, time)
    return LocalGate(gate.kind, gate.coefficient, gate.duration, _word_copy(gate.word),
                     gate.fock_unitary, gate.orthogonal, dict(gate.resources), dict(gate.diagnostics))


def vacuum_covariance(modes):
    """Covariance of declared blank fermion modes, not a native warmup proof."""
    if type(modes) is not int or not 1 <= modes <= 2*MAX_MODES:
        raise ValueError("a positive bounded number of modes required")
    gamma = np.zeros((2*modes, 2*modes))
    for mode in range(modes):
        gamma[2*mode, 2*mode+1], gamma[2*mode+1, 2*mode] = -1, 1
    return gamma


def validate_covariance(value):
    """Validate without skew symmetrizing, clipping or normalizing its state."""
    if np.ma.isMaskedArray(value):
        raise ValueError("masked covariance inadmissible")
    raw = np.asarray(value)
    if (raw.ndim != 2 or raw.shape[0] != raw.shape[1] or raw.shape[0] < 2
            or raw.shape[0] % 2 or raw.shape[0] > 4*MAX_MODES
            or raw.dtype.kind not in "iuf" or not np.isfinite(raw).all()):
        raise ValueError("finite real even-dimensional covariance required")
    gamma = raw.astype(float, copy=True)
    if np.max(abs(gamma+gamma.T)) > TOLERANCE:
        raise ValueError("antisymmetric Majorana covariance required")
    if np.linalg.eigvalsh(1j*gamma).max() > 1+TOLERANCE:
        raise ValueError("covariance violates the finite CAR state bound")
    return gamma


def occupations(covariance):
    """n_j=(1+Gamma[2j,2j+1])/2; no physical charge identification."""
    gamma = validate_covariance(covariance)
    return np.array([(1+gamma[2*j, 2*j+1])/2 for j in range(len(gamma)//2)])


def _quadratic(h, delta):
    if np.ma.isMaskedArray(h) or np.ma.isMaskedArray(delta):
        raise ValueError("masked quadratic coefficients inadmissible")
    raw_h, raw_d = np.asarray(h), np.asarray(delta)
    if (raw_h.ndim != 2 or raw_h.shape[0] != raw_h.shape[1] or not 1 <= len(raw_h) <= MAX_MODES
            or raw_d.shape != raw_h.shape or raw_h.dtype.kind not in "iufc"
            or raw_d.dtype.kind not in "iufc" or not np.isfinite(raw_h).all()
            or not np.isfinite(raw_d).all()):
        raise ValueError("finite matching quadratic coefficient matrices required")
    h, delta = raw_h.astype(complex), raw_d.astype(complex)
    hermitian, antisymmetric = float(np.max(abs(h-h.conj().T))), float(np.max(abs(delta+delta.T)))
    if hermitian > 2e-12 or antisymmetric > 2e-12:
        raise ValueError("Hermitian h and antisymmetric creation Delta required")
    terms = []
    for i in range(len(h)):
        if h[i, i].imag != 0:
            raise ValueError("onsite entries must be real")
        if h[i, i].real != 0:
            terms.append(("onsite", (i,), complex(h[i, i].real)))
    for kind, matrix in (("hopping", h), ("pairing", delta)):
        for i in range(len(h)):
            for j in range(i+1, len(h)):
                if matrix[i, j] != 0:
                    terms.append((kind, (i, j), complex(matrix[i, j])))
    return len(h), terms, {"input_h_Hermiticity_residual": hermitian,
                           "input_Delta_antisymmetry_residual": antisymmetric}


@dataclass(frozen=True)
class Evolution:
    orthogonal: np.ndarray
    covariance: np.ndarray | None
    resources: dict
    diagnostics: dict
    local_catalog: dict
    routed_block_sequence_sha256: str

    def record(self):
        return {"resources": dict(self.resources), "diagnostics": dict(self.diagnostics),
                "routed_block_sequence_sha256": self.routed_block_sequence_sha256,
                "local_source_catalog": self.local_catalog,
                "full_Fock_state_allocated": False,
                "target_BdG_exponential_used_by_source": False,
                "architecture_and_physical_services_derived": False}


def evolve(h, delta, duration, *, steps=1, covariance=None):
    """Symmetric local source evolution of H=c*hc+(Delta c*c*+h.c.)/2.

    Source matrices define U* w_a U=sum_b R_ab w_b and Gamma'=R Gamma R^T.
    Additional covariance modes are quantum references and remain untouched.
    A nonadjacent term moves its right mode left using FSWAP, applies its
    two-mode gate, then unroutes. Row updates compose those ACTUAL source
    Majorana maps, including the FSWAP numerical residual; there is no ideal
    permutation substitution or bosonic tensor embedding.
    """
    modes, terms, diagnostics = _quadratic(h, delta)
    time = _real(duration, "model duration")
    if type(steps) is not int or not 1 <= steps <= 4096:
        raise ValueError("a positive bounded integer split count required")
    gamma = None if covariance is None else validate_covariance(covariance)
    if gamma is not None and len(gamma) < 2*modes:
        raise ValueError("covariance omits a work mode")
    transform = np.eye(2*modes)
    resources = {key: 0 for key in _ADDITIVE_RESOURCES}
    resources.update({"completed_local_blocks": 0, "routing_fswap_calls": 0,
                      "onsite_calls": 0, "hopping_calls": 0, "pairing_calls": 0,
                      "retained_work_fermion_modes": modes,
                      "retained_reference_fermion_modes": 0 if gamma is None else len(gamma)//2-modes,
                      "maximum_executed_local_Fock_dimension": 4,
                      "source_mode_preparations": 0, "fresh_pointer_blanks": 0,
                      "physical_service_error_certified": False})
    catalog, catalog_keys = {}, {}
    digest = hashlib.sha256()
    max_local_gaussian_error = 0.

    def apply(gate, sites):
        nonlocal max_local_gaussian_error
        key = (gate.kind, gate.coefficient, gate.duration)
        if key not in catalog_keys:
            name = f"block_{len(catalog)}"
            catalog_keys[key] = name
            catalog[name] = gate.record()
        name = catalog_keys[key]
        digest.update(json.dumps([name, sites], separators=(",", ":")).encode())
        indices = [2*site+axis for site in sites for axis in (0, 1)]
        local = gate.orthogonal[:2, :2] if len(sites) == 1 else gate.orthogonal
        # Only four rows (two for onsite) are multiplied: O(M) per local block.
        transform[indices] = local @ transform[indices]
        for resource in _ADDITIVE_RESOURCES:
            resources[resource] += gate.resources[resource]
        resources["completed_local_blocks"] += 1
        resources["routing_fswap_calls" if gate.kind == "fswap" else gate.kind+"_calls"] += 1
        max_local_gaussian_error = max(max_local_gaussian_error, gate.diagnostics["majorana_linearity_residual"])

    dt = time/(2*steps)
    prepared = [(kind, sites, _local(kind, value.real, value.imag, dt)) for kind, sites, value in terms]
    swap = _local("fswap", 0., 0., 0.) if any(len(sites) == 2 and sites[1] > sites[0]+1
                                             for _, sites, _ in prepared) else None
    for _ in range(steps):
        for _, sites, gate in (*prepared, *reversed(prepared)):
            if len(sites) == 1:
                apply(gate, sites)
                continue
            left, right = sites
            route = [(j-1, j) for j in range(right, left+1, -1)]
            for pair in route:
                apply(swap, pair)
            apply(gate, (left, left+1))
            for pair in reversed(route):
                apply(swap, pair)
    output = None
    if gamma is not None:
        full = np.eye(len(gamma))
        full[:2*modes, :2*modes] = transform
        output = full @ gamma @ full.T
    diagnostics.update({"source_majorana_orthogonality_residual": float(np.max(abs(transform@transform.T-np.eye(2*modes)))),
                        "maximum_local_majorana_linearity_residual": max_local_gaussian_error,
                        "split_count": steps, "model_duration": time,
                        "quadratic_terms": len(terms), "local_source_catalog_size": len(catalog),
                        "all_local_Fock_columns_contracted": True,
                        "routing_includes_FSWAP_double_occupation_sign": True,
                        "reference_cross_covariances_retained": True,
                        "polar_or_covariance_normalization_applied": False,
                        "floating_execution_error_outward_certified": False,
                        "full_Fock_dimension_if_allocated": str(1 << modes)})
    if output is not None:
        diagnostics["output_covariance_antisymmetry_residual"] = float(np.max(abs(output+output.T)))
    return Evolution(_freeze(transform), None if output is None else _freeze(output), resources,
                     diagnostics, catalog, digest.hexdigest())


@dataclass(frozen=True)
class GaussianBranch:
    occupation: int
    probability: float
    weight: float
    covariance: np.ndarray | None


@dataclass(frozen=True)
class OccupationRead:
    branches: tuple[GaussianBranch, GaussianBranch]
    word: tuple[dict, ...]
    source_kraus: tuple[np.ndarray, np.ndarray]
    resources: dict
    diagnostics: dict


@lru_cache(maxsize=1)
def _occupation_pointer():
    # Bit0 is a field occupation; bit1 is a NEUTRAL record qubit, not a CAR mode.
    word = compiler.cnot(0, 1)
    matrix = compiler.execute(word, 2)
    kraus = (matrix[:2, :2], matrix[2:, :2])
    residual = max(float(np.max(abs(k-target))) for k, target in zip(
        kraus, (np.diag([1., 0.]), np.diag([0., 1.])), strict=True))
    if residual > TOLERANCE:
        raise ArithmeticError("source occupation pointer does not resolve its two projectors")
    return _word_copy(word), tuple(_freeze(k) for k in kraus), residual


def occupation_read(covariance, mode, *, weight=1.):
    """Both source-resolved occupation outcomes, with exact Gaussian conditioning.

    Gamma is a normalized Gaussian state descriptor; its probability weight
    is stored separately. All outcomes remain. Only an exactly zero-probability
    outcome has covarianceNone, since its conditional state is undefined.
    No clipping, success filtering or redistribution of weights occurs.
    Conditioning preserves the full unmeasured/reference covariance through
    Wick's finite-CAR identity. Strongly ill-conditioned roundoff remains a
    diagnostic limitation, not an invented state.
    """
    gamma = validate_covariance(covariance)
    if type(mode) is not int or not 0 <= mode < len(gamma)//2:
        raise ValueError("a retained fermion mode required")
    incoming = _real(weight, "branch weight")
    if incoming < 0:
        raise ValueError("nonnegative retained branch weight required")
    a, b = 2*mode, 2*mode+1
    raw = ((1-gamma[a, b])/2, (1+gamma[a, b])/2)
    if min(raw) < 0 or max(raw) > 1:
        raise ValueError("roundoff left the occupation probability domain; no clipping is applied")
    word, kraus, source_residual = _occupation_pointer()
    unmeasured = [i for i in range(len(gamma)) if i not in (a, b)]
    branches = []
    conditional_errors = []
    for occupied, probability in enumerate(raw):
        output = None
        if probability > 0:
            sign = 2*occupied-1
            output = gamma.copy()
            correction = (np.outer(gamma[:, a], gamma[b, :])-np.outer(gamma[:, b], gamma[a, :]))
            output[np.ix_(unmeasured, unmeasured)] += sign/(2*probability)*correction[np.ix_(unmeasured, unmeasured)]
            output[[a, b], :], output[:, [a, b]] = 0, 0
            output[a, b], output[b, a] = sign, -sign
            conditional_errors.append(max(0., float(np.linalg.eigvalsh(1j*output).max())-1))
            output = _freeze(output)
        branches.append(GaussianBranch(occupied, probability, incoming*probability, output))
    resources = {**compiler.resources(word), "fresh_pointer_blanks": 1,
                 "committed_classical_bits": 1, "source_mode_preparations": 0,
                 "physical_service_error_certified": False}
    diagnostics = {"source_projector_Kraus_residual": source_residual,
                   "maximum_conditional_CAR_bound_excess": max(conditional_errors, default=0.),
                   "all_outcomes_retained": True, "branch_weights_not_renormalized": True,
                   "neutral_record_qubit_is_not_a_fermion_mode": True,
                   "conditioning_uses_exact_Gaussian_Wick_identity": True,
                   "floating_and_physical_errors_outward_certified": False}
    return OccupationRead(tuple(branches), _word_copy(word), kraus, resources, diagnostics)
