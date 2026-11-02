"""Execute the two-site field and pointer using emitted twelve-port words.

The spectral decompositions below select a finite circuit. They never replace
an executed gate by its target matrix. Repeated subroutine matrices are cached
contractions of the source pulses, including their declared code services.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math

import numpy as np
from scipy.linalg import expm

from . import compiler as c, field
from .bounds import split_error_bound
from .cutoff import occupation_cutoff_bound, reference_probability_interval
from ..scalar_response.experiment import canonical, compare_receipt, exact_statistics, sample

TRIALS = 200_000
SEED = 20261012
DIMENSIONS = (4, 8, 16)
STEPS = (64, 128)
PINS = {
    c.RESPONSE: "226ab001894f357d3732d8a898fbae80e28b9e16ca41104a6f293190b441e053",
    c.CONTROL_RECEIPT: "cef14e5c14bdf02e14f21d6f470bcfc6c6f396dd301e8067ca6454194f25ab49",
    "model/observer/native.py": "dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d",
}


def sources():
    for path, expected in PINS.items():
        if hashlib.sha256((c.ROOT/path).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Source pin mismatch: {path}")
    return c.source_summary()


def encoded(state):
    return [[float(z.real), float(z.imag)] for z in state]


def _validate(d, steps, current=0):
    d = field.dimension(d)
    if d not in DIMENSIONS or type(steps) is not int or steps <= 0:
        raise ValueError("d=4, d=8 or d=16 and a positive integer step count required")
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError("a native current -1, 0 or 1 required")
    return d.bit_length()-1


def _spectral_word(values, basis_word, coefficient):
    return c.inverse(basis_word)+c.compile_diagonal(coefficient*values)+basis_word


@lru_cache(maxsize=8)
def program(d, steps):
    """Emit unique routines and the literal repeat recipe in little-endian order."""
    bits = _validate(d, steps)
    q_values, q_basis = field.canonical_eigenbasis(field.quadratures(d)[0])
    h_values, h_basis = field.canonical_eigenbasis(field.onsite_hamiltonian(d))
    q_word, h_word = c.compile_real_orthogonal(q_basis), c.compile_real_orthogonal(h_basis)
    target, _ = field.vacuum_vector(d)
    both_q = q_word+c.remap(q_word, list(range(bits, 2*bits)))
    # Pointer is the high bit of its local site+pointer routine. Prepare |+>,
    # couple exp(-i Z Q/2), then S^dag H readout maps Y=+1 to bit zero.
    hadamard = [c.pulse("H", math.pi/2, [bits]), c.pulse("I", -math.pi/2, [bits])]
    pointer = (hadamard+c.inverse(q_word)
               +c.compile_diagonal(np.concatenate((-.5*q_values, .5*q_values)))+q_word
               +[c.pulse("Z", math.pi/4, [bits]), c.pulse("I", -math.pi/4, [bits])]
               +hadamard)
    return {
        "qubits": 2*bits+1,
        "routines": {
            "prepare": c.compile_real_state(target),
            **{f"impulse_{j}": _spectral_word(q_values, q_word, j) for j in (-1, 0, 1)},
            "onsite_half": _spectral_word(h_values, h_word, -.5/steps),
            "bond": (c.inverse(both_q)
                     +c.compile_diagonal(float(field.KAPPA)*np.outer(q_values, q_values).ravel()/steps)
                     +both_q),
            "pointer": pointer,
        },
        "recipe": [
            {"routine": "prepare", "targets": list(range(2*bits))},
            {"routine": "impulse_j", "targets": list(range(bits))},
            {"repeat": steps, "body": [
                {"routine": "onsite_half", "targets": list(range(bits))},
                {"routine": "onsite_half", "targets": list(range(bits, 2*bits))},
                {"routine": "bond", "targets": list(range(2*bits))},
                {"routine": "onsite_half", "targets": list(range(bits))},
                {"routine": "onsite_half", "targets": list(range(bits, 2*bits))},
            ]},
            {"routine": "pointer", "targets": list(range(bits, 2*bits+1))},
            {"read": 2*bits, "outcomes": {"0": 1, "1": -1}},
        ],
        "initial_registers": "all blank |0>; pointer read is terminal",
        "services": "proper-code coherent load/export around EVERY primitive; complement failure outcomes retained",
    }


@lru_cache(maxsize=8)
def _executed(d, steps):
    bits = _validate(d, steps)
    words = program(d, steps)["routines"]
    diagnostics = {}

    def execute_routine(name, size):
        # Every exact pulse word is unitary. Roundoff from thousands of raw
        # matrix products can coherently accumulate over repeated routines.
        # Project ONLY the executed matrix, without consulting its target.
        # The guard rejects genuine nonunitarity and exposes the correction.
        unitary, diagnostics[name] = c.stabilize_executed_unitary(c.execute(words[name], size))
        return unitary

    prepared = execute_routine("prepare", 2*bits)[:, 0]
    half = execute_routine("onsite_half", bits)
    bond = execute_routine("bond", 2*bits)
    pointer = execute_routine("pointer", bits+1)
    impulses = {j: execute_routine(f"impulse_{j}", bits) for j in (-1, 0, 1)}
    half_both = np.kron(half, half)
    return prepared, half_both@bond@half_both, pointer, impulses, diagnostics


def _resource_counts(d, steps, current):
    words = program(d, steps)["routines"]
    uses = {"prepare": 1, f"impulse_{current}": 1, "onsite_half": 4*steps,
            "bond": steps, "pointer": 1}
    # The actual word is emitted on demand. The retained cross-platform receipt
    # binds stable schedule counts; exact floating-zero counts are word diagnostics.
    keys = ("native_pulses", "single_code_pulses", "pair_code_pulses",
            "identity_phase_pulses", "proper_code_load_events",
            "proper_code_export_events", "proper_code_service_events")
    resources = {name: c.resources(words[name]) for name in uses}
    return {**{key: sum(uses[name]*resources[name][key] for name in uses) for key in keys},
            "qubits": program(d, steps)["qubits"], "subroutine_uses": uses,
            "valid_code_failure_probability_in_exact_arithmetic": 0,
            "maximum_numerical_valid_code_leakage": max(
                row["maximum_numerical_valid_code_leakage"] for row in resources.values()),
            "initialization_storage_readout_physical_costs_included": False}


def source_history(d, steps, current, retain_program=False):
    """One actual source-compiled field instrument, before classical feedback."""
    _validate(d, steps, current)
    prepared, step, pointer, impulses, diagnostics = _executed(d, steps)
    state = np.kron(np.eye(d), impulses[current])@prepared
    for _ in range(steps):
        state = step@state
    # Pointer is initialized to zero by its blank register, never by a supplied
    # detector probability. This source-compiled unitary creates both branches.
    output = np.kron(pointer, np.eye(d))@np.concatenate((state, np.zeros(d*d)))
    branches = []
    for pointer_bit, outcome in ((0, 1), (1, -1)):
        branch = output[pointer_bit*d*d:(pointer_bit+1)*d*d]
        probability = float(np.vdot(branch, branch).real)
        branches.append({"outcome": outcome, "probability": probability,
                         "state": encoded(branch/math.sqrt(probability))})
    target, preparation = field.vacuum_vector(d)
    result = {"dimension": d, "steps": steps, "current": current, "model_time": "1",
              "preparation": preparation, "prepared_state": encoded(prepared),
              "finite_preparation_fidelity": float(abs(np.vdot(target, prepared))**2),
              "branches": branches, "resource_counts": _resource_counts(d, steps, current),
              "executed_routine_roundoff_diagnostics": diagnostics,
              "pre_read_state_norm": float(np.vdot(state, state).real),
              "post_read_total_probability": sum(row["probability"] for row in branches)}
    if retain_program:
        result["program"] = program(d, steps)
    return result


def _oracle(d, steps, current):
    """Independent target-matrix diagnostics, never used to produce source data."""
    q, _ = field.quadratures(d)
    onsite = field.onsite_hamiltonian(d)
    half = np.kron(np.eye(d), onsite)+np.kron(onsite, np.eye(d))
    u_half = expm(-.5j*half/steps)
    u_step = u_half@expm(1j*float(field.KAPPA)*np.kron(q, q)/steps)@u_half
    target, _ = field.vacuum_vector(d)
    forced = np.kron(np.eye(d), expm(1j*current*q))@target
    split_state = np.linalg.matrix_power(u_step, steps)@forced
    continuous = expm(-1j*field.physical_hamiltonian(d))@forced
    kraus = np.kron((expm(-.5j*q)-1j*expm(.5j*q))/2, np.eye(d))
    return {
        "split_plus_probability": float(np.linalg.norm(kraus@split_state)**2),
        "continuous_finite_plus_probability": float(np.linalg.norm(kraus@continuous)**2),
        "split_unitary_operator_error": float(np.linalg.norm(
            np.linalg.matrix_power(u_step, steps)-expm(-1j*field.physical_hamiltonian(d)), 2)),
    }


def run(trials=TRIALS, seed=SEED):
    source = sources()
    histories = exact_statistics()
    reference_interval = reference_probability_interval()
    checks = {}
    levels = []
    for d in DIMENSIONS:
        occupation = occupation_cutoff_bound(d)
        occupation_epsilon = F(occupation["instrument_trace_distance_upper"])
        rows = []
        for steps in STEPS:
            runs = [source_history(d, steps, current) for current in (-1, 0, 1)]
            probabilities = [item["branches"][0]["probability"] for item in runs]
            oracle = _oracle(d, steps, 1)
            bound = split_error_bound(d, steps)
            epsilon = F(bound["detector_probability_error_upper"])
            infinite_epsilon = min(F(1), occupation_epsilon+epsilon)
            ideal_lower = max(F(0), F(reference_interval["plus_probability_lower"])-infinite_epsilon)
            ideal_upper = min(F(1), F(reference_interval["plus_probability_upper"])+infinite_epsilon)
            response = probabilities[2]-.5
            prefix = f"d{d}_n{steps}"
            checks[prefix+"_instrument_normalized"] = all(abs(row["post_read_total_probability"]-1) < 2e-10 for row in runs)
            checks[prefix+"_blank_preparation_matches_target"] = all(abs(row["finite_preparation_fidelity"]-1) < 2e-10 for row in runs)
            checks[prefix+"_detector_parity"] = abs(probabilities[1]-.5) < 2e-10 and abs(probabilities[0]+probabilities[2]-1) < 2e-10
            checks[prefix+"_source_word_matches_independent_split"] = abs(probabilities[2]-oracle["split_plus_probability"]) < 2e-10
            checks[prefix+"_finite_exact_error_below_rational_bound"] = oracle["split_unitary_operator_error"] <= float(epsilon)
            checks[prefix+"_record_zero_current_has_full_service_schedule"] = len({row["resource_counts"]["native_pulses"] for row in runs}) == 1
            if d == 16:
                checks[prefix+"_ideal_probability_interval_excludes_fair_response"] = ideal_lower > F(1, 2)
            rows.append({"steps": steps, "response": response, "histories": runs,
                         "target_matrix_diagnostics": oracle,
                         "finite_split_certificate": bound,
                         "ideal_classical_history_tv_error_upper": str(F(5, 33)*epsilon),
                         "ideal_native_current_error_upper": str(F(5, 132)*epsilon),
                         "ideal_infinite_instrument_certificate": {
                             "reference_plus_probability_interval": reference_interval,
                             "instrument_trace_distance_upper": str(infinite_epsilon),
                             "detector_probability_error_upper": str(infinite_epsilon),
                             "ideal_finite_plus_probability_lower": str(ideal_lower),
                             "ideal_finite_plus_probability_upper": str(ideal_upper),
                             "ideal_nonzero_response_certified": ideal_lower > F(1, 2),
                             "classical_history_tv_error_upper_under_parity": str(F(5, 33)*infinite_epsilon),
                             "native_current_error_upper_under_parity": str(F(5, 132)*infinite_epsilon),
                             "composition": "continuous occupation cutoff plus finite Strang error by triangle inequality; the parity factors apply only to classical native records",
                             "scope": "ideal finite split instrument versus infinite two-oscillator reference, including terminal field states; no numerical execution, physical service or spatial continuum bound",
                         },
                         "forecast_second_native_current": 5*response/132})
        sampled = sample(rows[-1]["response"], trials, seed)
        for name, item in sampled["statistics"].items():
            checks[f"d{d}_{name}_within_six_standard_errors"] = abs(item["mean"]-item["forecast"]) <= 6*item["standard_error"]+1/trials
        checks[f"d{d}_finite_step_refinement_improves"] = rows[1]["target_matrix_diagnostics"]["split_unitary_operator_error"] < rows[0]["target_matrix_diagnostics"]["split_unitary_operator_error"]
        checks[f"d{d}_continuous_occupation_certificate_nontrivial"] = 0 < occupation_epsilon < 1
        levels.append({"dimension": d,
                       "continuous_instrument_occupation_certificate": occupation,
                       "step_refinement": rows, "sample": sampled})
    reference = field.ideal_prediction()
    differences = [abs(row["step_refinement"][-1]["histories"][-1]["branches"][0]["probability"]-reference["plus_probability"]) for row in levels]
    checks["occupation_refinement_improves_infinite_reference_diagnostic"] = all(
        later < earlier for earlier, later in zip(differences, differences[1:]))
    occupation_refinement = [
        {"dimension": d, "executed": d in DIMENSIONS,
         "execution_scope": ("source histories retained in levels" if d in DIMENSIONS
                             else "analytic occupation certificate only; no source program emitted or executed at this dimension"),
         "continuous_instrument_occupation_certificate": occupation_cutoff_bound(d)}
        for d in (4, 8, 16, 32, 64)
    ]
    occupation_errors = [F(row["continuous_instrument_occupation_certificate"]["instrument_trace_distance_upper"])
                         for row in occupation_refinement]
    checks["analytic_occupation_refinement_bound_decreases"] = all(
        later < earlier for earlier, later in zip(occupation_errors, occupation_errors[1:]))
    result = {
        "schema": "oph.observer.source_scalar.v1", "source_pins": PINS, "source": source,
        "scope": {
            "construction": "native repair current -> emitted source preparation, scalar evolution and quantum pointer -> committed bit -> native repair feedback",
            "target": "two sites with onsite occupation d=4,8,16, supplied mass=1 and spring kappa=9/32; not the earlier spatial-refinement family",
            "time": "declared scalar model coordinate t=1; not a source pulse clock, solver counter or calibrated physical clock",
            "preparation": "real entangled projected vacuum compiled from blank registers; target selected by supplied scalar action",
            "execution": "full M6 exponentials of twelve source controls; checked polar stabilization of source-executed routine matrices before cached repetition; no target-matrix substitution",
            "instrument": "both terminal conditional field states retained; only committed pointer bit enters feedback",
            "services": "coherent proper-code load/export, blank registers, storage, Born readout and analogue controls supplied",
            "rigorous_errors": "complete continuous occupation cutoff and ideal finite-Hamiltonian splitting bounded and composed by exact rational certificates, including both terminal field branches; d=32,64 occupation certificates are analytic refinements without source execution",
            "numeric_errors": "source compilation and full branch execution checked in floating arithmetic, not interval-certified",
            "open": "spatial continuum from this code, physical service budget, certified numerical execution error, emergent action/geometry, empirical physics validation",
        },
        "native_history_law": histories, "levels": levels,
        "occupation_refinement": occupation_refinement,
        "infinite_two_oscillator_reference": reference,
        "infinite_two_oscillator_reference_interval": reference_interval,
        "absolute_reference_differences_diagnostic_only": differences,
        "checks": checks, "all_pass": all(checks.values()),
    }
    if not result["all_pass"]:
        raise ValueError(f"Source scalar checks failed: {[name for name, ok in checks.items() if not ok]}")
    return result
