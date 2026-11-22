"""Replay the record-action measurements and compare the complete public receipt."""
import argparse
import json
import math
from pathlib import Path

from model.executable_model import cell, federation
from . import adjacency_bound, merged_signs, record_action, record_clock, record_readiness
from . import native_frequency, record_release, registered_forecast, sign_graph


RECEIPT = Path(__file__).with_name("record_measurements.json")
ARMS = ("primitive", "assembled", "off", "double", "untransported")


def circuit():
    return dict(local=record_action.run(), historical_contact=native_frequency.compare(), readiness=record_readiness.check(),
                donor_error=record_clock.parity_clock(), marker_error=record_clock.prepare_marker(),
                routed=record_clock.literal_lift(), shared=sign_graph.microscopic_path(),
                private=sign_graph.microscopic_path(True), masks=sign_graph.witness_masks(),
                projection={arm: sign_graph.full_projection(arm) for arm in ARMS[:2]},
                opposite_chart=sign_graph.full_projection("primitive", True),
                adjacency=[adjacency_bound.certificate(b) for b in (cell(), *(federation(i) for i in range(3)))])


def verify(actual, saved, path="receipt"):
    if type(actual) is not type(saved):
        raise ValueError(f"Different type at {path}")
    if isinstance(actual, dict):
        if actual.keys() != saved.keys():
            raise ValueError(f"Different fields at {path}")
        for key in actual:
            verify(actual[key], saved[key], f"{path}/{key}")
    elif isinstance(actual, list):
        if len(actual) != len(saved):
            raise ValueError(f"Different lengths at {path}")
        for i, (a, b) in enumerate(zip(actual, saved, strict=True)):
            verify(a, b, f"{path}/{i}")
    elif isinstance(actual, float):
        # Acceptance is computed from spectral lines and can be exactly zero
        # or many orders of magnitude below the detector's rounding tolerance.
        field = path.rsplit("/", 1)[-1]
        acceptance = field == "seventeen_read_acceptance"
        if acceptance and not (0 <= actual <= 1 and 0 <= saved <= 1):
            raise ValueError(f"Invalid probability at {path}")
        if field == "maximum_forbidden_read_probability" and not (
            0 <= actual < sign_graph.FREQUENCY_LEAKAGE_LIMIT and 0 <= saved < sign_graph.FREQUENCY_LEAKAGE_LIMIT
        ):
            raise ValueError(f"Clock leakage exceeds its registered limit at {path}")
        if not (math.isfinite(actual) and math.isfinite(saved)
                and math.isclose(actual, saved, abs_tol=0. if acceptance else 2e-11, rel_tol=2e-12)):
            raise ValueError(f"Different measurement at {path}")
    elif actual != saved:
        raise ValueError(f"Different result at {path}")


def replay(extended):
    packet = dict(circuit=circuit())
    if extended:
        packet["release"] = {f"{arm}/{n}": record_release.release(arm, n)
                             for n, arms in ((8, ARMS), (16, ARMS[:3])) for arm in arms}
        packet["externally_enrolled"] = [registered_forecast.compare(a) for a in ARMS[:3]]
        packet["merged_destinations"] = [merged_signs.compare(a) for a in ARMS[:2]]
    return json.loads(json.dumps(packet, allow_nan=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extended", action="store_true")
    parser.add_argument("--write", type=Path, help="create a new complete measurement file")
    args = parser.parse_args()
    if args.write and not args.extended:
        parser.error("--write requires the complete replay")
    measured = replay(args.extended)
    if args.write:
        with args.write.open("x", encoding="utf-8") as out:
            json.dump(measured, out, allow_nan=False, indent=2)
            out.write("\n")
    else:
        saved = json.loads(RECEIPT.read_text(encoding="utf-8"))
        if set(saved) != {"circuit", "release", "externally_enrolled", "merged_destinations"}:
            raise ValueError("Incomplete registered packet")
        verify(measured, saved if args.extended else {"circuit": saved["circuit"]})
        print("Verified all registered measurements." if args.extended else "Verified the circuit and block bounds.")


if __name__ == "__main__":
    main()
