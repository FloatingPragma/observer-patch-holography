"""Compose deterministic native preparation with the complete phase readout.

Exact bounds reuse the native sweep certificate and a hash-bound complete
phase/work receipt. Numerical execution is a diagnostic, not its proof.
Warm-up records are ignored without conditioning. All subsequent heralds,
phase bits and unmeasured timeouts are retained.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

from ..native_equilibration import sweeps
from ..clocked_quantum_readout import construction as readout
from ..clocked_quantum_readout import experiment
from ..primitive_observer import runtime

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENT = "model/clocked_quantum_readout/results.json"
PARENT_SHA256 = "440ea8fce66ea5d16dd21788a36bc4a0fb92278a3d2b62651d67a897221d284d"
TOLERANCE = F(1, 2**24)


def preparation_certificate(tolerance=TOLERANCE):
    """Uniform exact bound, including one full-carrier absorption sweep.

    The all-law envelope admits fixed C/E and the native current-Gram law
    family under the controller contract in sweeps.py. It selects neither
    a source law nor the total-two sector. The fixed scan visits each seam.
    """
    report = sweeps.budget(tolerance, "all", full=True)
    # Even the vacuous tolerance1 needs a covering sweep to make the
    # binary readout domain available on every full-carrier input.
    if report["sweeps"] == 0:
        report = dict(sweeps.certificate(1, "all", full=True),
                      requested_tolerance=F(tolerance))
    return report


def parent_receipt():
    """Require the exact retained receipt and its source bindings."""
    raw = (ROOT/PARENT).read_bytes()
    if hashlib.sha256(raw).hexdigest() != PARENT_SHA256:
        raise ValueError("the complete phase/work parent receipt changed")
    data = json.loads(raw)
    if data.get("schema") != "oph.clocked_quantum_readout.v1" or data.get("all_pass") is not True:
        raise ValueError("a passing complete phase/work receipt is required")
    for path, digest in data["source_sha256"].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
            raise ValueError(f"phase/work source binding changed: {path}")
    return data


def native_outputs(rho, law, cap, *, reference_dimension=1):
    """Execute the complete ideal native suffix on an arbitrary binary input.

    Source index precedes the reference index. Subnormalization is retained.
    The input need not have the block-diagonal form used by the parent's
    exact rational oracle. This floating execution claims no rounding bound.
    """
    law = readout.clock.quantum.candidate(law)
    readout.clock._horizon(cap)
    state = sweeps.execute(rho, law, sweeps=0,
                          reference_dimension=reference_dimension, full=False)
    herald, pointer = readout.interface.herald_and_pointer(law)
    balanced = runtime.Instrument(kraus=((readout.interface.mixer(),),))
    initial = runtime._apply(herald, state, reference_dimension)
    outputs = {(0, 0, 0, -1): initial[0]}
    source = readout.clock.source(law).instrument
    for flag in (1, 2):
        waiting = runtime._apply(balanced, initial[flag], reference_dimension)[0]
        for n in range(1, cap+1):
            negative, waiting, positive = runtime._apply(source, waiting, reference_dimension)
            for current, stopped in ((-1, negative), (1, positive)):
                analyzed = runtime._apply(balanced, stopped, reference_dimension)[0]
                phases = runtime._apply(pointer, analyzed, reference_dimension)
                for bit in (0, 1):
                    outputs[flag, n, current, bit] = phases[bit if flag == 1 else 1-bit]
        outputs[flag, cap, 0, -1] = waiting
    return {record: outputs[record] for record in readout.records(cap)}


def execute(rho, law, cap, *, reference_dimension=1, tolerance=TOLERANCE):
    """Numerical native prefix and suffix on the full total-two input.

    The prefix ignores its tape, preserves the original reference marginal,
    and performs no successful-trial selection or reset. The final mapping
    retains every suffix record. No work Hamiltonian is executed by this
    native-only diagnostic; the complete work extension follows by the
    all-input contraction proof.
    """
    prep = preparation_certificate(tolerance)
    state = sweeps.execute(rho, law, prep["sweeps"],
                           reference_dimension=reference_dimension, full=True)
    binary = tuple(i for i, x in enumerate(sweeps.FULL_STATES) if max(x) == 1)
    if tuple(sweeps.FULL_STATES[i] for i in binary) != tuple(readout.clock.quantum.STATES):
        raise ArithmeticError("native preparation/readout basis mismatch")
    indices = (np.asarray(binary)[:, None]*reference_dimension
               +np.arange(reference_dimension)).ravel()
    complement = np.setdiff1d(np.arange(state.shape[0]), indices)
    if np.any(state[complement]) or np.any(state[:, complement]):
        raise ArithmeticError("the covering sweep did not absorb collision support")
    return native_outputs(state[np.ix_(indices, indices)], law, cap,
                          reference_dimension=reference_dimension)


def _interval(center, radius):
    return [str(center-radius), str(center+radius)]


def run():
    parent = parent_receipt()
    prep = preparation_certificate()
    epsilon = prep["trace_distance_upper"]
    implementation, _ = experiment.compiled_service()
    rows = []
    for row in parent["histories"]:
        law, cap = row["candidate"], row["cap"]
        exact = experiment.exact_native_statistics(law, cap)
        if str(exact["phase_match_and_success"]) != row["exact_native_phase_match_and_success"]:
            raise ArithmeticError("the exact native parent event changed")
        phase = exact["phase_match_and_success"]
        rows.append({"law": law, "cap": cap,
            "complete_joint_suffix_records": 128*cap+3,
            "ideal_phase_match_probability_interval": _interval(phase, epsilon),
            # Compiled binary matrices need not be exactly trace preserving.
            # The event-weight allowance therefore uses their full norm error.
            "source_mixer_phase_match_weight_interval": _interval(phase, epsilon+2*implementation),
            "complete_output_half_trace_error_from_parent_report_upper": str(
                epsilon+implementation+F(row["finite_ideal_service_output_error_upper"])),
            "maximum_native_attempts_including_prefix_herald_and_pointer":
                prep["native_attempts"]+cap+2,
            "exact_stationary_preparation_or_exact_herald_rate_assumed": False})
    difference = (experiment.exact_native_statistics("C", 64)["phase_match_and_success"]
                  -experiment.exact_native_statistics("E", 64)["phase_match_and_success"])
    files = ("model/native_equilibration/sweeps.py", "model/observer/native.py",
        "model/observer/quantum_mixing.py", "model/prepared_quantum_readout/certificate.py",
        "model/prepared_quantum_readout/PROOF.md")
    return {"schema": "oph.prepared_quantum_readout.v1", "all_pass": True,
        "parent": {"path": PARENT, "sha256": PARENT_SHA256},
        "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
        "preparation": prep,
        "source_mixer_two_call_half_norm_error_upper": implementation,
        "histories": rows,
        "cap64_source_mixer_C_minus_E_phase_match_weight_interval":
            _interval(difference, 2*epsilon+4*implementation),
        "scope": {"arbitrary_full_total_two_input_and_untouched_reference": True,
            "fixed_deterministic_seam_scan": True,
            "preparation_tape_ignored_without_selection": True,
            "continuation_cannot_access_preparation_tape_or_controller": True,
            "source_native_fair_coin_quantum_instrument_contract_supplied": True,
            "total_two_conservation_sector_selected_by_theorem": False,
            "balanced_coherent_service_work_action_work_preparation_and_model_clock_supplied": True,
            "complete_suffix_failure_signs_work_records_and_timeouts_retained": True,
            "exact_stationarity_after_finite_prefix": False,
            "Born_rule_derived_from_scalar_consensus": False,
            "physical_thermal_energy_temperature_or_EOS_derived": False,
            "unlimited_source_mixer_terminal_state_claimed": False}}


def _encode(value):
    if isinstance(value, F):
        return str(value)
    raise TypeError(type(value).__name__)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = run()
    data = json.dumps(receipt, default=_encode, indent=2, sort_keys=True)+"\n"
    path = HERE/"results.json"
    if args.write:
        path.write_text(data)
    if args.check and path.read_text() != data:
        raise SystemExit("prepared readout receipt differs from exact reproduction")
    print(json.dumps({"sweeps": receipt["preparation"]["sweeps"],
        "attempts": receipt["preparation"]["native_attempts"],
        "preparation_error_upper": str(receipt["preparation"]["trace_distance_upper"]),
        "complete_suffix_histories": len(receipt["histories"])}))


if __name__ == "__main__":
    main()
