"""An unread coherence table does not specify the committed instrument."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path

from tools import strict_json

from . import family as f


def remixed_native_c():
    """Change the measurement basis of two native record branches.

    The two chosen C branches have disjoint input supports, so the existing
    sparse single-output-per-column representation also represents their
    sum and difference. No source transition or unread channel changes.
    """
    native = f.INSTRUMENTS["C"]()
    left = native.pop(("coin_current", 0, 0))
    right = native.pop(("coin_current", 0, 1))
    if left.keys() & right.keys():
        raise AssertionError("the chosen native branches overlap on inputs")
    for sign in (1, -1):
        branch = {x: (y, a*f.ROOT_HALF) for x, (y, a) in left.items()}
        branch.update({x: (y, sign*a*f.ROOT_HALF) for x, (y, a) in right.items()})
        native[("record_superposition", sign)] = branch
    return native


def probability(branch, rho):
    output = f.apply_instrument({"outcome": branch}, rho)
    return sum((v for (x, y), v in output.items() if x == y), f.S())


def run():
    native = f.INSTRUMENTS["C"]()
    remixed = remixed_native_c()
    table = f.extract_table(native)
    x, y = f.representative((1, 0)), f.representative((0, 1))
    diagonal = {(x, x): f.S(1)}
    coherent = {(a, b): f.S(F(1, 2)) for a, b in product((x, y), repeat=2)}
    incoherent = {(a, a): f.S(F(1, 2)) for a in (x, y)}
    plus, minus = (remixed[("record_superposition", s)] for s in (1, -1))
    original = native[("coin_current", 0, 1)]
    same_channel = f.same_channel_everywhere(remixed, table)
    covariant, _ = f.difference_preserving(remixed)
    checks = {
        "native_complete": f.completeness(native),
        "remixed_complete": f.completeness(remixed),
        "same_unread_channel_on_all_matrix_units": same_channel,
        "same_unread_number_phase_covariance": covariant,
        "same_unread_classical_kernel": f.classical_kernel(native) == f.classical_kernel(remixed),
        "native_branches_have_definite_current": f.definite_currents(native),
        "remixed_branches_mix_currents": not f.definite_currents(remixed),
    }
    if not all(checks.values()):
        raise AssertionError("instrument/channel control failed")
    return {
        "schema": "oph.record_laws.instrument_contract.v1",
        "carrier_states": len(f.STATES),
        "matrix_units_checked": len(f.STATES)**2,
        "checks": checks,
        "witness": {
            "input_loads": [list(f.STATES[x]), list(f.STATES[y])],
            "basis_input_original_moving_record_probability": str(probability(original, diagonal)),
            "basis_input_remixed_plus_probability": str(probability(plus, diagonal)),
            "coherent_input_remixed_plus_probability": str(probability(plus, coherent)),
            "coherent_input_remixed_minus_probability": str(probability(minus, coherent)),
            "incoherent_input_remixed_plus_probability": str(probability(plus, incoherent)),
            "incoherent_input_remixed_minus_probability": str(probability(minus, incoherent)),
        },
        "scope": {
            "unread_channel_classification_preserved": True,
            "complete_instrument_selected_by_four_channel_requirements": False,
            "remixed_label_is_a_native_current": False,
            "physical_record_basis_control_implemented": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", type=Path)
    mode.add_argument("--write", type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        if strict_json.load(args.check) != result:
            raise ValueError("instrument-contract receipt differs")
        print("Unread-channel equality and distinct committed instruments replayed.")
        return
    payload = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n"
    if args.write:
        with args.write.open("x", encoding="utf-8") as stream:
            stream.write(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
