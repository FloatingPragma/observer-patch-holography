"""Construct the exact stationary reference of the adaptive source/bit kernel.

This is the same r=1, epsilon=1/2 controller as adaptive_reference.py, with
classical C/E records reduced to their exposed current sign. It prepares no
coherence and compiles no modular action. The stationary vector is computed
with exact rational linear algebra and certified by its exact residual and
the actual transition support.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import sympy as sp

from . import adaptive_reference
from .family import binding
from ..common_source.native_defect_decision import strongly_connected_components
from ..observer import native
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def transition_counts():
    """Complete classical projection, denominator120 per input row."""
    states = tuple((x, memory) for x in native.states(binary=True) for memory in (0, 1))
    index = {state: i for i, state in enumerate(states)}
    counts = [[0] * len(states) for _ in states]
    for i, state in enumerate(states):
        for output, mass in adaptive_reference.advance({state: F(1)}).items():
            integer = 120 * mass
            if integer.denominator != 1:
                raise ArithmeticError("declared rational kernel denominator changed")
            counts[i][index[output]] = integer.numerator
        if sum(counts[i]) != 120:
            raise ArithmeticError("joint transition row is not normalized")
    return states, counts


def run():
    states, counts = transition_counts()
    dimension = len(states)
    # Solve (K^T-I)pi=0 with one dependent equation replaced by sum pi=1.
    # Inverting over the rational domain also establishes uniqueness of this
    # normalized system; the exact residual below checks the replaced row.
    system = sp.Matrix(counts).T - 120 * sp.eye(dimension)
    system[dimension - 1, :] = sp.ones(1, dimension)
    rhs = sp.zeros(dimension, 1)
    rhs[dimension - 1] = 1
    solution = system.inv(method="DM") * rhs
    probabilities = [F(int(p.p), int(p.q)) for p in solution]
    if sum(probabilities) != 1 or min(probabilities) < 0:
        raise ArithmeticError("stationary solution is not a probability law")
    residual = [sum(probabilities[i] * counts[i][j] for i in range(dimension))
                - 120 * probabilities[j] for j in range(dimension)]
    if any(residual):
        raise ArithmeticError("stationary residual is nonzero")

    components = strongly_connected_components(counts)
    closed = [component for component in components
              if all(not counts[i][j] for i in component
                     for j in range(dimension) if j not in component)]
    support = tuple(i for i, p in enumerate(probabilities) if p > 0)
    if closed != [support]:
        raise ArithmeticError("positive stationary support is not the unique closed recurrent class")
    loops = [i for i in support if counts[i][i] > 0]
    if not loops:
        raise ArithmeticError("aperiodicity witness missing")
    joint = {state: p for state, p in zip(states, probabilities) if p}
    source = adaptive_reference.summarize(joint)
    if source["changed_source_entries"] == 0:
        raise ArithmeticError("adaptive stationary source unexpectedly became uniform")
    files = (
        "model/completion/stationary_reference.py",
        "model/completion/adaptive_reference.py",
        "model/common_source/native_defect_decision.py",
        "model/observer/native.py",
        "model/common_source/slow_geometry.py",
    )
    return {
        "schema": "oph.completion.stationary_reference.v1",
        "gate": "CP-02",
        "campaign": binding(),
        "canonical_premise_navigation": ["PR-07", "PR-23"],
        "source_sha256": bound_pins({p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in files}),
        "contract": {
            "kernel": "adaptive_reference.py with r=1, epsilon=1/2, m'=1[current>0]",
            "source_dimensions": 66,
            "classical_memory_states": 2,
            "joint_states": dimension,
            "tagged_record_laws": "all tagged C/E convex mixtures have this current-sign marginal",
            "time_parameter": "native attempt count",
            "initial_preparation": "uniform binary source and blank current-sign bit; supplied finite initial state",
            "active_recovery": False,
            "hold_or_halt_in_this_member": False,
        },
        "transition_denominator": 120,
        "stationary_joint": [
            {"state": list(x), "memory": memory, "probability": str(p)}
            for (x, memory), p in zip(states, probabilities)
        ],
        "stationary_source": source,
        "exact_certificate": {
            "normalized_system_unique_over_rationals": True,
            "stationarity_residual": "0",
            "normalization": "1",
            "joint_support_dimension": len(support),
            "joint_zero_probability_states": [
                {"state": list(states[i][0]), "memory": states[i][1]}
                for i, p in enumerate(probabilities) if not p
            ],
            "communicating_component_sizes": [len(c) for c in components],
            "closed_recurrent_components": 1,
            "closed_recurrent_size": len(support),
            "self_loop_witness": {"state": list(states[loops[0]][0]), "memory": states[loops[0]][1],
                                  "probability": str(F(counts[loops[0]][loops[0]], 120))},
            "source_marginal_support_dimension": sum(p > 0 for p in adaptive_reference.source_marginal(joint).values()),
        },
        "conclusions": [
            "The actual finite source/bit classical kernel has one stationary law, positive exactly on its unique129-state recurrent class.",
            "That class is aperiodic by the retained self-loop, so every initial distribution converges to this law in attempt count; no finite accuracy/mixing-time certificate is supplied.",
            "The stationary source marginal is nonuniform and positive on all66 binary basis states.",
            "Native branches preserve diagonal inputs, so this invariant classical law also gives an invariant diagonal density for the actual quantum/classical source-bit channel.",
        ],
        "boundaries": [
            "The joint density is not faithful on the full132-dimensional source-bit space; it is faithful on its129-dimensional support.",
            "Uniqueness is established for classical source/bit distributions, not all fixed states of the complete quantum channel.",
            "Support compression is an explicit mathematical algebra choice, not a proof of physical observer access to all support matrix units.",
            "The source marginal is not an autonomous closed source process: its stationary invariance belongs to the correlated joint source-bit state.",
            "On the accessible commutative record algebra the modular automorphism flow is trivial; nonuniform source weights alone supply no implemented noncommuting modular clock.",
            "No physical clock, energy normalization, G1--G6 realization or integrated CP-02 acceptance follows.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", nargs="?", type=Path,
                      const=HERE / "stationary_reference_results.json",
                      help="replay and check the retained receipt (or a supplied path)")
    args = parser.parse_args()
    result = run()
    if args.write:
        (HERE / "stationary_reference_results.json").write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.check is not None:
        check_receipt(result, args.check)
    print(json.dumps(result["exact_certificate"], indent=2))


def check_receipt(result, path):
    """Check the full exact state, custody and support rather than a summary."""
    if json.loads(Path(path).read_text(encoding="utf-8")) != result:
        raise ValueError("stationary-reference receipt differs from exact replay")


if __name__ == "__main__":
    main()
