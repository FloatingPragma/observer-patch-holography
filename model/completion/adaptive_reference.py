"""Exact countercontrol for transferring fixed-mixture stationarity to feedback.

The source is the existing total-two binary carrier. No new physical clock,
energy, quantum control or state is introduced. A finite classical bit selects
one of two preferred seams, mixed with uniform exploration. This demonstrates
why a stationary reference for each seam need not remain a stationary source
marginal after record-controlled scheduling correlates source and memory.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from ..observer import native
from ..common_source import slow_geometry
from .family import binding
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PREFERRED = ((0, 1), (0, 5))


def initial_joint():
    """Uniform binary source and declared blank one-bit classical memory."""
    states = native.states(binary=True)
    return {(x, 0): F(1, len(states)) for x in states}


def advance(joint, *, exploration=F(1, 2), repair=F(1), adaptive=True):
    """One exact attempt; identity holds retain the bit in this controller.

    On repair the bit becomes 1[current > 0]. Both C and E expose that
    current and have the same diagonal transition law. A tagged convex
    mixture consequently has the same source/bit marginal calculated here.
    This function does not replace their different complete quantum records.
    """
    if not isinstance(exploration, F) or not isinstance(repair, F):
        raise TypeError("exact Fraction probabilities required")
    if not 0 <= exploration <= 1 or not 0 <= repair <= 1:
        raise ValueError("probabilities must lie in [0,1]")
    if type(adaptive) is not bool:
        raise TypeError("explicit Boolean feedback selector required")
    states = set(native.states(binary=True))
    seams = tuple(slow_geometry.source().seams())
    out = defaultdict(F)
    for (x, memory), mass in joint.items():
        if x not in states or type(memory) is not int or memory not in (0, 1):
            raise ValueError("binary source and one-bit memory required")
        if not isinstance(mass, F) or mass < 0:
            raise ValueError("exact nonnegative input masses required")
        if not mass:
            continue
        out[x, memory] += mass * (1 - repair)
        preferred = PREFERRED[memory if adaptive else 0]
        for seam in seams:
            probability = exploration / len(seams)
            if seam == preferred:
                probability += 1 - exploration
            for coin in (0, 1):
                y, (_, current) = native.native(x, seam, coin)
                out[y, int(current > 0)] += mass * repair * probability / 2
    return {key: value for key, value in out.items() if value}


def source_marginal(joint):
    out = {x: F(0) for x in native.states(binary=True)}
    for (x, _), mass in joint.items():
        out[x] += mass
    return out


def summarize(joint):
    marginal = source_marginal(joint)
    uniform = F(1, len(marginal))
    deviations = {x: p - uniform for x, p in marginal.items()}
    return {
        "total_probability": str(sum(marginal.values())),
        "changed_source_entries": sum(d != 0 for d in deviations.values()),
        "maximum_absolute_deviation": str(max(abs(d) for d in deviations.values())),
        "total_variation_from_uniform": str(sum(abs(d) for d in deviations.values()) / 2),
        "source_probabilities": [
            {"state": list(x), "probability": str(p), "deviation_from_uniform": str(deviations[x])}
            for x, p in marginal.items()
        ],
    }


def run():
    blank = initial_joint()
    first = advance(blank)
    second = advance(first)
    fixed_second = advance(advance(blank, adaptive=False), adaptive=False)
    uniform_second = advance(advance(blank, exploration=F(1)), exploration=F(1))
    sources = (
        "model/completion/adaptive_reference.py",
        "model/observer/native.py",
        "model/common_source/slow_geometry.py",
        "evidence/native_geometric_source_20260925/spec.json",
    )
    return {
        "schema": "oph.completion.adaptive_reference.v1",
        "campaign": binding(),
        "gate": "CP-13",
        "canonical_premise_navigation": ["PR-07"],
        "source_sha256": bound_pins({p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in sources}),
        "contract": {
            "source": "existing 66-state binary shell of the total-two twelve-port carrier",
            "initial_state": "I_66/66 tensor blank classical bit m=0",
            "record_laws": "C, E, and every tagged convex mixture have this same diagonal marginal",
            "repair_probability": "1",
            "uniform_exploration_probability": "1/2",
            "preferred_seams_by_memory": [list(e) for e in PREFERRED],
            "memory_update_on_repair": "m'=1[current>0]",
            "hold_rule": "identity on source and memory; not used in this r=1 witness",
            "time_parameter": "native attempt count, not physical duration",
            "retention": "one current-sign bit; complete records may be retained separately without changing this marginal",
        },
        "first_attempt": summarize(first),
        "second_attempt": summarize(second),
        "fixed_policy_control": summarize(fixed_second),
        "uniform_schedule_control": summarize(uniform_second),
        "conclusion": "The initially uniform source marginal fails to remain uniform under this admissible record-feedback scheduler.",
        "scope": [
            "This refutes transfer of fixed state-independent mixture stationarity to arbitrary record-controlled source marginals.",
            "It does not refute mixture_stationary or its proof: the source/memory correlations violate that theorem's fixed-mixture identification.",
            "It does not exclude an invariant reference for the joint finite source/controller kernel or a nonuniform source marginal.",
            "It does not decide physical clock construction, defect protection, or the broader executable hypothesis.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true", help="retain this exact countercontrol beside the source")
    mode.add_argument("--check", nargs="?", type=Path,
                      const=HERE / "adaptive_reference_results.json",
                      help="replay and check the retained receipt (or a supplied path)")
    args = parser.parse_args()
    result = run()
    if args.write:
        target = HERE / "adaptive_reference_results.json"
        target.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.check is not None:
        check_receipt(result, args.check)
    print(json.dumps({"second_attempt": {k: v for k, v in result["second_attempt"].items()
                                               if k != "source_probabilities"},
                      "scope": result["scope"]}, indent=2))


def check_receipt(result, path):
    """Require the entire retained packet to equal the newly replayed result."""
    if json.loads(Path(path).read_text(encoding="utf-8")) != result:
        raise ValueError("adaptive-reference receipt differs from exact replay")


if __name__ == "__main__":
    main()
