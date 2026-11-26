"""Replay every enrolled simulated detector and compare with the Fourier forecast."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .gauge import flux_register
from .native import native_rows
from .walk import forecast, pointer, read_grating, rest_properties

HERE = Path(__file__).resolve().parent
PROTOCOL = json.loads((HERE / "protocol.json").read_text())
SOURCES = ("protocol.json", "walk.py", "native.py", "gauge.py", "__main__.py", "../law.py")


def digest():
    return {name: hashlib.sha256((HERE / name).read_bytes().replace(b"\r\n", b"\n")).hexdigest()
            for name in SOURCES}


# The November acquisition used the verifier in commit 8917041e6. These
# identities bind all deposited counts and the exact original source list.
ORIGINAL_VERIFIER = "861e63b442f4e0c58b145eb7dd3029660b8e31c479205784cf11147d4394aaa7"
DEPOSIT_DIGEST = "c9adc3bf650a1716dc6dd2bc7e567db0df1ed197ca71b30cd625ab5b243d578f"


def archived_acquisition(receipt):
    if set(receipt) != {"counts", "sources"}:
        raise ValueError("unexpected receipt fields")
    current = digest()
    if receipt["sources"] == current:
        return False
    original = current | {"__main__.py": ORIGINAL_VERIFIER}
    if receipt["sources"] != original:
        raise ValueError("source ancestry differs from the retained acquisition")
    encoded = json.dumps(receipt, separators=(",", ":"), sort_keys=True).encode()
    if hashlib.sha256(encoded).hexdigest() != DEPOSIT_DIGEST:
        raise ValueError("historical seed replay differs from the authenticated deposit")
    return True


def replay(retained=None, *, require_bitwise=False):
    historical = retained is not None and archived_acquisition(retained)
    differing_rows = 0
    largest_difference = 0
    p = PROTOCOL
    tol = p["tolerance"]
    rng = np.random.default_rng(p["seed"])
    trials = p["trial_count"]
    rows = {}
    max_residual = 0.
    max_score = 0.
    for period in p["grating_periods"]:
        # The smaller carrier also resolves the two transverse public axes.
        shape = (period, 3, 3) if period == min(p["grating_periods"]) else (period, 1, 1)
        for amplitude in p["grating_amplitudes"]:
            for strength in (amplitude, -amplitude):
                for charge in p["charges"]:
                    for spin in range(2):
                        for spectator in (False, True):
                            key = f"{period}/{strength}/{charge}/{spin}/{int(spectator)}"
                            args = (p["grating_private_length"], charge, strength, spin)
                            signal = read_grating(shape, *args, spectator=spectator)
                            expected = forecast(period, *args, spectator=spectator)
                            for channel, (actual, theory) in enumerate(zip(signal, (expected, 0j))):
                                record = f"{key}/{channel}"
                                delta = abs(actual - theory)
                                if not np.isfinite(delta) or delta > tol:
                                    raise ArithmeticError(f"wrong amplitude at {record}: {actual} vs {theory}")
                                max_residual = max(max_residual, float(delta))
                                counts = rng.multinomial(trials, pointer(actual)).tolist()
                                if retained is not None:
                                    observed = retained["counts"].get(record)
                                    if observed != counts:
                                        if not historical or require_bitwise:
                                            raise ValueError(f"seed replay differs at {record}")
                                        differing_rows += 1
                                        largest_difference = max(largest_difference, max(
                                            abs(a - b) for a, b in zip(observed, counts, strict=True)
                                        ))
                                    counts = observed
                                rows[record] = counts
                                prob = pointer(theory)
                                err = np.abs(np.array(counts) - trials * prob)
                                sd = np.sqrt(trials * prob * (1 - prob))
                                if sum(counts) != trials or np.any(err > 8 * sd + 1):
                                    raise ArithmeticError(f"count gate failed at {record}")
                                max_score = max(max_score, float(np.max(err / np.maximum(1., sd))))
    if retained is not None and rows.keys() != retained["counts"].keys():
        raise ValueError("unexpected records")
    flux_minus, bad_minus = flux_register(-1)
    flux_plus, bad_plus = flux_register(1)
    reflect = (-np.arange(flux_plus.shape[0])) % flux_plus.shape[0]
    inverse = flux_plus[np.ix_(reflect, reflect, reflect, np.arange(3))]
    if np.max(abs(inverse - flux_minus)) > tol:
        raise ArithmeticError("electric register failed charge reversal")
    native = native_rows(p["native_period"])
    summary = dict(settings=len(rows), preparations=trials * len(rows), max_count_sd=max_score,
                   max_amplitude_error=max_residual, native=native, gauss_failure=max(bad_minus, bad_plus),
                   charged_arrival=float(flux_plus[..., 1:].sum()),
                   lengths={str(n): rest_properties(n) for n in p["private_lengths"]})
    if retained is not None:
        summary["same_seed_replay"] = dict(differing_rows=differing_rows, max_bin_difference=largest_difference)
    return dict(sources=retained["sources"] if retained is not None else digest(), counts=rows), summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path, help="write an independent count file to an unused path")
    parser.add_argument("--require-bitwise", action="store_true",
                        help="also require the historic seeded draws to match on this platform")
    args = parser.parse_args()
    if args.require_bitwise and args.acquire is not None:
        parser.error("--require-bitwise checks the retained acquisition")
    if args.acquire is None:
        receipt = json.loads((HERE / "counts.json").read_text())
        _, summary = replay(receipt, require_bitwise=args.require_bitwise)
    else:
        data, summary = replay()
        with args.acquire.open("x", encoding="utf-8") as output:
            json.dump(data, output, separators=(",", ":"), sort_keys=True)
            output.write("\n")
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
