"""Exact real-operation obstruction for the constructed reference's modular flow.

This attaches the existing native real-operation boundary to unequal positive
weights in an actual stationary full-record block. It concerns the fixed
source-basis full matrix algebra of that block. It is not an exclusion of
every clock, encoding, carrier refinement or operational interpretation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import sympy as sp

from . import observer, observer_reference, stationary_reference
from .family import binding
from ..common_source.native_defect_decision import complete_record_maps
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def real_primitive_certificate():
    """Check the actual native support maps underlying the real-Kraus proof.

    Each native Kraus matrix is the checked integer partial-injection matrix
    times a real nonnegative scalar. Changing classical controller weights
    changes that scalar and the classical output address, not source phases.
    """
    counts = {}
    for law in ("C", "E"):
        rows = [row for row in complete_record_maps() if row["law"] == law]
        nonzero = 0
        addresses = 0
        for row in rows:
            pairs = row["pairs"]
            if any(type(i) is not int or type(j) is not int
                   or not 0 <= i < 78 or not 0 <= j < 78 for i, j in pairs):
                raise ArithmeticError("invalid native matrix address")
            if (len({i for i, _ in pairs}) != len(pairs)
                    or len({j for _, j in pairs}) != len(pairs)):
                raise ArithmeticError("native map is not a partial injection")
            nonzero += bool(pairs)
            addresses += len(pairs)
        counts[law] = {"complete_raw_branch_matrices": len(rows),
                       "nonzero_branch_matrices": nonzero,
                       "integer_matrix_entries_equal_one": addresses,
                       "unweighted_matrix_is_real": True,
                       "kraus_scalar": "1/sqrt(2), times real square roots of scheduler/law weights"}
    return {"C_E": counts,
            "holds_and_timeouts": "identity matrices; real",
            "classical_memory_operations": "classical transitions represented by real basis maps with nonnegative probability square roots",
            "closure": "products, classical sums, partial traces and positive-probability normalization preserve real density matrices",
            "ancilla_premise": "all supplied ancillary joint states are real in the declared basis; admitted diagonal independent ancillas suffice",
            "scope": "arbitrary finite admitted histories; no coherent phase service or change of quantum encoding"}


def select_witness(blocks):
    """Choose unequal positive weights in one actual stationary record block."""
    for memory, density in sorted(blocks.items(), key=lambda item: repr(item[0])):
        entries = []
        for (i, a, j, b), probability in sorted(density.items()):
            if (i, a) != (j, b) or a != 0 or probability <= 0:
                raise ArithmeticError("expected the constructed positive diagonal block reference")
            entries.append((i, probability))
        if not entries:
            continue
        low = min(entries, key=lambda item: (item[1], item[0]))
        high = max(entries, key=lambda item: (item[1], -item[0]))
        if high[1] > low[1]:
            return memory, entries, high, low
    raise ArithmeticError("no nontracial stationary block found")


def phase_control(ratio):
    """Exact two-dimensional modular witness; no numerical logarithm required."""
    if not isinstance(ratio, F) or ratio <= 1:
        raise ValueError("an exact rational ratio greater than one is required")
    q = sp.Rational(ratio.numerator, ratio.denominator)
    modular_parameter = sp.pi / (2 * sp.log(q))
    if sp.simplify(modular_parameter * sp.log(q)) != sp.pi / 2:
        raise ArithmeticError("modular phase cancellation failed")
    half = sp.Rational(1, 2)
    real_probe = sp.Matrix([[half, half], [half, half]])
    # The predual flow has relative phase exp(-i*s*log(q))=-i.
    modular_unitary_up_to_phase = sp.diag(-sp.I, 1)
    target = modular_unitary_up_to_phase * real_probe * modular_unitary_up_to_phase.adjoint()
    pauli_y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    if target != sp.Matrix([[half, -sp.I / 2], [sp.I / 2, half]]):
        raise ArithmeticError("incorrect modular orientation")
    if sp.trace(target * pauli_y) != 1 or target * target != target:
        raise ArithmeticError("target is not the positive-Y pure state")
    # A deliberately excluded complex phase gate attains the target.
    added_phase_gate = sp.diag(1, sp.I)
    if added_phase_gate * real_probe * added_phase_gate.adjoint() != target:
        raise ArithmeticError("complex-service countercontrol failed")
    closest_real_state = sp.eye(2) / 2
    difference = target - closest_real_state
    if difference * difference != sp.eye(2) / 4:
        raise ArithmeticError("sharp trace-distance countercontrol failed")
    return {"modular_parameter": "s_star=pi/(2*log(q)) > 0; dimensionless",
            "ratio_q": str(ratio),
            "predual_orientation": "rho -> sigma^(-i*s) rho sigma^(i*s)",
            "observable_orientation": "A -> sigma^(i*s) A sigma^(-i*s)",
            "relative_predual_phase_at_s_star": "-i",
            "probe": "(|i>+|j>)/sqrt(2), a real full-algebra test input",
            "target_density_on_ordered_pair": [["1/2", "-i/2"], ["i/2", "1/2"]],
            "target_pauli_y_expectation": "1",
            "every_normalized_real_output_pauli_y_expectation": "0",
            "trace_distance_lower_bound": "1/2",
            "sharp_real_state_control": "diag(1/2,1/2) has trace distance1/2 from target",
            "excluded_phase_service_control": "diag(1,i) maps the real probe to the target exactly; this complex gate is outside the campaign"}


def reference_witness(stationary_joint, law):
    c_weight = F(1) if law == "C" else F(0)
    blocks = observer_reference.lift(stationary_joint, c_weight)
    memory, entries, high, low = select_witness(blocks)
    block_mass = sum(p for _, p in entries)
    ratio = high[1] / low[1]
    encoded = [[i, str(p)] for i, p in entries]
    return {
        "law": law,
        "record": observer_reference.serialize_record(memory[0]),
        "positive_block_support_dimension": len(entries),
        "block_probability": str(block_mass),
        "block_exact_diagonal_sha256": hashlib.sha256(json.dumps(encoded, separators=(",", ":")).encode()).hexdigest(),
        "ordered_source_indices_i_j": [high[0], low[0]],
        "ordered_source_states_i_j": [list(observer.STATES[high[0]]), list(observer.STATES[low[0]])],
        "unnormalized_positive_weights": [str(high[1]), str(low[1])],
        "conditional_positive_weights": [str(high[1] / block_mass), str(low[1] / block_mass)],
        "central_record_weight_cancels_from_ratio": True,
        "phase_witness": phase_control(ratio),
    }


def report():
    quotient = stationary_reference.run()
    witnesses = {law: reference_witness(quotient["stationary_joint"], law) for law in ("C", "E")}
    files = (
        "model/completion/modular_access.py",
        "model/completion/observer.py",
        "model/completion/observer_reference.py",
        "model/completion/stationary_reference.py",
        "model/completion/adaptive_reference.py",
        "model/common_source/native_defect_decision.py",
        "model/native_interface/README.md",
    )
    return {
        "schema": "oph.completion.modular_access.v1",
        "gate": "CP-02",
        "campaign": binding(),
        "status": "FIXED_REFERENCE_ALGEBRA_IMPLEMENTATION_OBSTRUCTED",
        "full_CP02_gate_accepted": False,
        "whole_program_excluded": False,
        "source_sha256": bound_pins({path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in files}),
        "nearest_existing_result": "model/native_interface/README.md: the exact reachability boundary and real-effect complex-conjugation control",
        "addition": "Attach that real-operation invariant to unequal weights of the actual CP01 full-record stationary references and obtain an exact1/2 implementation-error witness for their nontrivial modular action on the fixed support algebra.",
        "primitive_certificate": real_primitive_certificate(),
        "reference_witnesses": witnesses,
        "proof": [
            "Every admitted native Kraus branch, hold and classical memory map is real in the fixed source/record basis. Classical probability weights multiply these maps by real scalars.",
            "A finite product maps a real joint density to a real subnormalized density. Summing classical branches, tracing ancillary registers, or dividing by a positive real branch probability preserves realness.",
            "In an actual stationary block choose positive unequal weights p_i>p_j. On its full support matrix algebra, the predual modular map at s_star=pi/(2log(p_i/p_j)) sends the real pair probe to the positive-Y pure state.",
            "Extend PauliY by zero outside that pair and record block. It has operator norm1. The target expectation is1, whereas every real normalized admitted output has expectation0, even with leakage or a changed classical record.",
            "Trace-norm duality therefore gives trace distance at least1/2 from the target. A positive-probability heralded native branch remains subject to the same bound. Norm limits cannot evade the closed real-output condition.",
        ],
        "boundaries": [
            "This excludes implementing the stated nontrivial modular action on the fixed source-basis full support algebra of these two references using the frozen operational grammar.",
            "The real pair probe tests an all-input map on that algebra; it is not claimed to be preparable by the admitted diagonal-input native grammar.",
            "PauliY is a mathematical trace-distance witness, not an executed or admitted native detector.",
            "The 1/2 bound compares normalized outputs; an unnormalized success branch of probability p contributes a p/2 bound when its success flag is retained.",
            "On the accessible commutative record algebra the modular automorphism flow is trivial; the result does not rule out its identity implementation.",
            "A new complex phase resource, a different quantum encoding, additional coherent registers, different reference/algebra, or carrier refinement changes the tested contract and is not excluded here.",
            "Supplied joint ancillary states must be real; individually real marginals would not suffice if they conceal complex joint correlations.",
            "A supplied non-real encoder or decoder is also outside the declared native implementation grammar.",
            "This is not a no-go for every modular-clock construction, a physical clock calibration, or the broader executable hypothesis.",
        ],
    }


def check_receipt(result, path):
    if json.loads(Path(path).read_text(encoding="utf-8")) != result:
        raise ValueError("modular-access receipt differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", nargs="?", type=Path, const=HERE / "modular_access_results.json")
    args = parser.parse_args()
    result = report()
    if args.write:
        (HERE / "modular_access_results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if args.check is not None:
        check_receipt(result, args.check)
    print(json.dumps({law: {"record": row["record"],
                                  "source_pair": row["ordered_source_indices_i_j"],
                                  "error_lower_bound": row["phase_witness"]["trace_distance_lower_bound"]}
                      for law, row in result["reference_witnesses"].items()}, indent=2))


if __name__ == "__main__":
    main()
