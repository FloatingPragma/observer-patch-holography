"""Execute and certify all 66 transported-core clock words at parameter four.

This is the explicit successor service inventory, not a change to the native
repair law or the earlier four-word receipt. Every original 512 by 512 source
array is checked by rational phase enclosures. Proper-code controls, coherent
load/export and program storage retain their supplied-apparatus contracts.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np

from model.completion import observer
from model.completion.family import FamilyParameters
from model.native_modular_reference import common_caps, full_flow, history_clock, history_transport, local_preparation
from . import continuation
from .joins_contract import binding, check_receipt, source_bindings
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "clock_services_results.json"
PARAMETER = F(4)
OPERATOR_GATE = F(1, 10**9)
DIAGNOSTIC_CEILING = 5e-10
THEOREMS = ("Lean/Dynamics/PrivateInner.lean", "Lean/Dynamics/ChoiCPTP.lean")


def ratios(core):
    weights = continuation.reference(core)
    values = [F(1)]*512
    for address, weight in zip(full_flow.support()[0], weights, strict=True):
        values[address] = 60*weight
    return tuple(values)


def word(core):
    """The existing closed Gray-Walsh source grammar; no new primitive."""
    result = common_caps.diagonal_word(PARAMETER, ratios(core))
    if len(result) != 2042 or set(event["kind"] for event in result) != {"I", "Z", "H", "J"}:
        raise ArithmeticError("existing source grammar changed")
    return result


def certify(core, original_array):
    """Inspect every original component, including all 446 padding addresses."""
    bound = common_caps.certify_diagonal(PARAMETER, ratios(core), original_array)
    if bound > OPERATOR_GATE:
        raise ArithmeticError("clock word exceeds the certified operator gate")
    return bound


@lru_cache(maxsize=66)
def _verified(core):
    # Cache immutable original bytes, never a caller-mutable source array.
    continuation._core(core)
    original = full_flow.compiler.execute(word(core), 9)
    certify(core, original)
    return original.dtype.str, original.tobytes()


def matrix(core):
    """A fresh copy of the actually executed and checked source word."""
    dtype, payload = _verified(continuation._core(core))
    return np.frombuffer(payload, dtype=np.dtype(dtype)).reshape(512, 512).copy()


def export_outputs(original_array, density):
    """Success and the full reference-bearing abort; neither is normalized.

    The input is on H66 tensor R. The export traces out the escaped address,
    retaining its arbitrary reference state. This low-level map also permits
    injected arrays for independent failure controls; service() uses only
    the verified source array.
    """
    if (type(original_array) is not np.ndarray or original_array.shape != (512, 512)
            or original_array.dtype.kind != "c" or not np.all(np.isfinite(original_array))):
        raise ValueError("an original finite complex 512-square array required")
    if (type(density) is not np.ndarray or density.ndim != 2 or density.shape[0] != density.shape[1]
            or not density.shape[0] or density.shape[0] % 66
            or density.dtype.kind not in "iufc" or not np.all(np.isfinite(density))):
        raise ValueError("a finite H66 tensor reference matrix required")
    success, abort = history_clock._export(original_array, density)
    return {"success": success, "export_abort": abort}


def service(core, density):
    return export_outputs(matrix(core), density)


def after_preparation(blocks, core):
    """Consume the predecessor's complete flagged output; run only on success."""
    if type(blocks) is not dict or set(blocks) != set(local_preparation.FLAGS):
        raise ValueError("all predecessor preparation flags required")
    if any(type(value) is not np.ndarray for value in blocks.values()):
        raise ValueError("original reference-bearing preparation arrays required")
    shape = blocks["success"].shape
    if (len(shape) != 2 or shape[0] != shape[1] or not shape[0] or shape[0] % 66
            or any(value.shape != shape or value.dtype.kind not in "iufc"
                   or not np.all(np.isfinite(value)) for value in blocks.values())):
        raise ValueError("all preparation flags must share one H66 tensor reference space")
    result = {("preparation", flag): value.copy() for flag, value in blocks.items()
              if flag != "success"}
    result.update({("clock", flag): value for flag, value in service(core, blocks["success"]).items()})
    return result


def word_row(core):
    _verified(core)  # Do not certify an unexecuted word by its target alone.
    tape = word(core)
    encoded = full_flow.corner.canonical(tape).encode()
    kinds = Counter(event["kind"] for event in tape)
    return {"core": list(core), "word_sha256": hashlib.sha256(encoded).hexdigest(),
            "source_program_utf8_bytes": len(encoded), "native_pulses": len(tape),
            "pulse_kind_counts": dict(sorted(kinds.items())),
            "zero_angle_requests": sum(event["angle"] == 0 for event in tape),
            "proper_code_load_events": len(tape), "proper_code_export_events": len(tape),
            "original_array_entries_checked": 512**2,
            "operator_error_certified_below": str(OPERATOR_GATE),
            "padding_identity_extension_certified": True}


def error_budget(requests=24):
    if type(requests) is not int or not 0 <= requests <= 256:
        raise ValueError("integer request count in 0..256 required")
    gain = (1+OPERATOR_GATE)**2
    return {"clock_requests": requests, "operator_error_per_word": str(OPERATOR_GATE),
            "per_word_half_diamond_upper": str((gain-1)/2),
            "per_word_channel_norm_upper": str(gain),
            "per_word_native_export_abort_probability_upper": str(OPERATOR_GATE**2),
            "complete_interleaved_half_diamond_upper": str((gain**requests-1)/2),
            "native_instruments_in_comparison": "exact complete C/E instruments from continuation",
            "physical_service_error_budget_supplied": False}


def _probe():
    """An explicitly supplied entangled input, not a reference re-preparation."""
    configurations = ((0, 1), (0, 4))
    vector = np.zeros((66, 2), complex)
    for r, occupied in enumerate(configurations):
        i = continuation.BASIS.index(tuple(int(p in occupied) for p in range(12)))
        vector[i, r] = (1 if r == 0 else 1j)/np.sqrt(2)
    return vector


def _native_vector(vector, branch):
    result = np.zeros_like(vector)
    if branch.weight:
        for i, j in branch.pairs:
            result[j] = np.sqrt(float(branch.weight))*vector[i]
    return result


def _source_vector(core, vector):
    codes = np.array(full_flow.support()[0])
    output = matrix(core)[:, codes]@vector
    outside = np.array([i for i in range(512) if i not in set(codes)])
    escaped = output[outside]
    return output[codes], np.einsum("ar,as->rs", escaped, escaped.conj())


def _pure_distance(left, right):
    """Numerical trace-distance diagnostic, not the rigorous word bound."""
    # Bound ||aa* - bb*||_1, including unequal masses.
    a, b = left.ravel(), right.ravel()
    # Direct vector-distance majorant avoids cancellation when outputs agree.
    return (np.linalg.norm(a)+np.linalg.norm(b))*np.linalg.norm(a-b)/2


def complete_one_attempt(law):
    """Clock-before versus record-conditioned clock-after, with all outcomes.

    This uses actual source words and the exact native branch maps. It keeps
    every semantic zero slot and every export abort, including references.
    The source-word theorem supplies the rigorous arbitrary-input bound;
    floating vector propagation here is a separately labelled diagnostic.
    """
    if law not in ("C", "E"):
        raise ValueError("C or E required")
    params = FamilyParameters(c_weight=F(int(law == "C")))
    initial = _probe()
    before, before_abort = _source_vector(continuation.INITIAL_CORE, initial)
    residual, mass_before, mass_after = 0., float(np.trace(before_abort).real), 0.
    after_abort = np.zeros((2, 2), complex)
    positive, newly_covered = 0, set()
    returned_before = before_abort.copy()
    returned_after = np.zeros((2, 2), complex)
    menu = continuation.branches((), continuation.INITIAL_CORE, params)
    for branch in menu:
        left = _native_vector(before, branch)
        middle = _native_vector(initial, branch)
        if np.any(middle):
            positive += 1
            if branch.core not in continuation.compiled_cores():
                newly_covered.add(branch.core)
            right, abort = _source_vector(branch.core, middle)
        else:
            right, abort = middle, np.zeros((2, 2), complex)
        residual += _pure_distance(left, right)
        mass_before += float(np.vdot(left, left).real)
        mass_after += float(np.vdot(right, right).real+np.trace(abort).real)
        after_abort += abort
        returned_before += np.einsum("ar,as->rs", left, left.conj())
        returned_after += np.einsum("ar,as->rs", right, right.conj())+abort
    # Abort placement is kept as a distinct full flag in the comparison.
    residual += (float(np.trace(before_abort).real)+float(np.trace(after_abort).real))/2
    reference = np.einsum("ar,as->rs", initial, initial.conj())
    numerical = max(residual, abs(mass_before-1), abs(mass_after-1),
                    float(np.linalg.norm(returned_before-reference)),
                    float(np.linalg.norm(returned_after-reference)))
    if numerical > DIAGNOSTIC_CEILING or not newly_covered:
        raise ArithmeticError("complete source-clock continuation diagnostic failed")
    return {"law": law, "semantic_native_slots_retained": len(menu),
            "positive_probe_outcomes": positive, "new_cores_exercised": sorted(newly_covered),
            "clock_before_export_abort_retained": True,
            "clock_after_export_abort_retained_per_native_outcome": True,
            "arbitrary_reference_certificate": "two times the per-word half-diamond bound",
            "floating_output_residual_below_diagnostic_ceiling": str(F(1, 2_000_000_000)),
            "numerical_residual_is_the_rigorous_bound": False}


def exported_history(law):
    """Four actual repairs move one occupied port through newly compiled cores."""
    params = FamilyParameters(c_weight=F(int(law == "C")))
    history = (((0, 7), 0, 1), ((6, 7), 1, -1),
               ((2, 6), 1, -1), ((2, 3), 0, 1))
    core, memory, mapping, exports, weight = (0, 1), (), tuple((i, i) for i in range(66)), [], F(1)
    tape = []
    for seam, coin, current in history:
        # On a moving binary branch the E contrast equals the signed current.
        record = observer.Record("repair", "explore", seam, law, coin, current, current)
        branch = next(b for b in continuation.branches(memory, core, params) if b.record == record)
        step = dict(branch.pairs)
        mapping = tuple((i, step[j]) for i, j in mapping if j in step)
        if branch.exported is not None:
            exports.append(branch.exported)
        memory, core, weight = branch.memory, branch.core, weight*branch.weight
        tape.append(record)
    if tuple(exports)+memory != tuple(tape) or len(mapping) < 2 or core in continuation.compiled_cores():
        raise ArithmeticError("exported history failed to reach a new coherent corner")
    partial = np.zeros((66, 66), complex)
    for i, j in mapping:
        partial[j, i] = 1
    codes = np.array(full_flow.support()[0])
    before, after = (matrix(c)[np.ix_(codes, codes)] for c in ((0, 1), core))
    residual = float(np.linalg.norm(after@partial-partial@before, "fro"))
    if residual > DIAGNOSTIC_CEILING:
        raise ArithmeticError("new-core clock failed history covariance")
    return {"law": law, "complete_native_tape": [r.__dict__ for r in tape],
            "final_core": list(core), "branch_weight": str(weight),
            "history_corner_rank": len(mapping), "retained_records": len(memory),
            "exported_records": len(exports), "reference_reprepared": False,
            "clock_service_executed": True,
            "partial_isometry_intertwiner_operator_bound": str(2*OPERATOR_GATE),
            "floating_frobenius_residual_below_diagnostic_ceiling": str(F(1, 2_000_000_000))}


def run(progress=None):
    campaign = binding()
    rows = []
    for position, core in enumerate(continuation.CORES, 1):
        rows.append(word_row(core))
        if progress is not None and (position % 6 == 0 or position == 66):
            progress(position)
    old = continuation.compiled_cores()
    return {"schema": "oph.clock_defect_energy.clock_services.v1", "campaign": campaign,
            "source_sha256": bound_pins(source_bindings(set((
                "model/clock_defect_energy/clock_services.py", "model/clock_defect_energy/test_clock_services.py",
                "model/clock_defect_energy/continuation.py", "model/clock_defect_energy/joins_contract.py",
                "model/clock_defect_energy/contract.py", "tools/strict_json.py",
                "model/native_modular_reference/history_transport.py",
                "model/native_modular_reference/history_clock.py", "model/native_modular_reference/common_caps.py",
                "model/native_modular_reference/local_preparation.py",
                "model/native_modular_reference/full_flow.py", "model/observer/source_scalar/compiler.py",
                "code/source_selection_model/response.json", "code/m1_source_realization/receipt.json", *THEOREMS))
                | set(common_caps.SOURCE_PATHS) | set(history_transport.source_bindings())
                | {"model/clock_defect_energy/CLOCK_SERVICES.md"})),
            "inventory": {"modular_parameter": str(PARAMETER), "native_dimension": 66,
                          "source_register_dimension": 512, "predecessor_cores": sorted(old),
                          "source_words_executed_and_certified": len(rows),
                          "additional_words": len(rows)-len(old), "all_reached_cores_covered": True,
                          "other_parameters_require_fresh_source_certification": True},
            "words": rows,
            "resources": {"source_program_utf8_bytes_total": sum(r["source_program_utf8_bytes"] for r in rows),
                          "retained_programs": 66, "core_selector_states": 66,
                          "core_selector_bits": 7, "core_times_feedback_states": 132,
                          "quantum_address_qubits": 9, "source_pulses_per_request": 2042,
                          "load_export_service_calls_per_request": 4084,
                          "source_word_position_counter_bits": 11,
                          "source_pulses_for_24_clock_requests": 24*2042,
                          "load_export_service_calls_for_24_clock_requests": 24*4084,
                          "native_records_after_24_attempts": 24, "fifo_records": 1,
                          "exported_native_records_after_24_attempts": 23,
                          "source_interpreter_program_counter_counted_separately": True,
                          "bounded_total_environment_storage": False},
            "error_budget": error_budget(),
            "complete_one_attempt": {law: complete_one_attempt(law) for law in ("C", "E")},
            "exported_history": {law: exported_history(law) for law in ("C", "E")},
            "finite_preparation": history_transport.preparation_certificate(PARAMETER),
            "failure_contract": {"per_pulse_proper_code_exports": "all complement outcomes retained; exact invariant-code success on declared inputs",
                                 "terminal_native_export": "success plus reference-bearing padding abort",
                                 "exact_target_terminal_abort_probability": "0",
                                 "floating_output_mass_renormalized": False,
                                 "failed_preparation_flags": "unchanged predecessor failures bypass the conditional clock service"},
            "theorem_map": [
                {"path": THEOREMS[0], "theorems": ["finitePrivateStarAutomorphism_inner", "hamiltonianPropagator_mem_unitary"],
                 "model_premises": "Hilbert space C^66; faithful positive diagonal sigma_core; Hamiltonian H=-log(sigma_core); theorem time t=-s=-4 gives exp(-itH)=sigma_core^(-is), up to the explicitly compiled common native phase",
                 "application": "The specified modular action is unitary conjugation; Walsh source compilation implements that chosen action with the retained error.",
                 "does_not_supply": "source controls, a physical clock or a new formalized model instance"},
                {"path": THEOREMS[1], "theorems": ["isCompletelyPositive_of_kraus", "IsCPTP.comp", "IsCPTP.convexCombination"],
                 "model_premises": "complete native partial-injection Kraus branches; invariant-code exports; explicit success and complement export maps",
                 "application": "Finite Kraus amplification/composition justifies the arbitrary-reference ideal instrument; source error uses the separate exact numerical-array bound.",
                 "does_not_supply": "physical service error, free record erasure or equivalence of clock and field generators"}],
            "scope": {"new_native_primitive": False, "supplied_complex_controls": True,
                      "service_inventory_extended_under_successor_contract": True,
                      "predecessor_unavailable_service_receipt_rewritten": False,
                      "stationary_reference_substituted": False, "physical_clock_calibrated": False,
                      "integrated_cp01_accepted": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = run(lambda n: print(f"certified transported-core source words: {n}/66", flush=True))
    if args.write:
        OUTPUT.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    elif args.check:
        check_receipt(OUTPUT, report)
    print("all66 transported clock services: PASS")


if __name__ == "__main__":
    main()
