"""Invert the X meter using only calibration and launch counts."""

import argparse
from hashlib import sha256
import json
from math import ceil, exp, expm1, floor, log, log1p
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "results.json"
RECEIPT = HERE / "results.json"
ALPHA = 0.01


def check(condition, message):
    if not condition:
        raise ValueError(message)


def bernoulli_kl(p, q):
    # The bisections use interior q; endpoint observations take a closed form.
    return p * log1p((p - q) / q) + (1 - p) * log1p((q - p) / (1 - q))


def interval(k, n, log_penalty):
    """Invert each one-sided binomial Chernoff bound, rounding outwards."""
    check(0 <= k <= n and n > 0, "invalid eligible trial count")
    if not k:
        return 0.0, -expm1(-log_penalty / n)
    if k == n:
        return exp(-log_penalty / n), 1.0
    p = k / n

    def solve(outside, inside):
        for _ in range(46):
            mid = (outside + inside) / 2
            if n * bernoulli_kl(p, mid) >= log_penalty:
                outside = mid
            else:
                inside = mid
        return outside

    return solve(0.0, p), solve(1.0, p)


def charge_interval(observed, reference):
    # P(dark) is the convex mixture of three calibrated meter responses.
    r0lo, r0hi = reference[0]
    low = min(reference[-1][0], reference[1][0])
    high = max(reference[-1][1], reference[1][1])
    check(low > r0hi, "X meter cannot distinguish a charged standard from its host")
    lo = max(0.0, (observed[0] - r0hi) / (high - r0hi))
    hi = min(1.0, (observed[1] - r0lo) / (low - r0lo))
    check(lo <= hi, "charge calibration is incompatible with the launch counts")
    return lo, hi


def calibration(standards, shots, penalty):
    reference = {1: {}, 2: {}}
    for s in standards:
        counts = s["counts"]
        check(len(counts) == 9 and sum(counts) == shots, "invalid charge standard")
        reference[s["multiplier"]][s["region_charge"]] = interval(
            counts[2], shots, penalty
        )
    return reference


def original_cells(samples, standards, shots):
    # Three native launch cells, two dwell banks each. In each bank the
    # dark fraction among eligible trials and the refusal fraction in the
    # entire filtered band are separate binomial observations.
    selected = [
        s for s in samples if s[1] == "native_arrival" and s[0] == s[2] and s[3] == 0
    ]
    check(len(selected) == 3, "missing a native source cell")
    parameters = len(standards) + 4 * len(selected)
    penalty = log(2 * parameters / ALPHA)
    reference = calibration(standards, shots, penalty)
    rows = []
    for owner, _, _, _, readings in selected:
        bounds = []
        for s in readings:
            counts = s["counts"]
            eligible = sum(counts[4:13])
            refused = sum(counts[:4])
            _, fhi = interval(refused, refused + eligible, penalty)
            low, high = charge_interval(
                interval(counts[6], eligible, penalty), reference[s["multiplier"]]
            )
            bounds.append(((1 - fhi) * low, fhi + (1 - fhi) * high))
        lo, hi = max(b[0] for b in bounds), min(b[1] for b in bounds)
        check(lo <= hi, "conflicting whole-band charge estimates")
        rows.append(
            {
                "owner": owner,
                "region_charge_interval": [
                    floor(lo * 1e10) / 1e10,
                    ceil(hi * 1e10) / 1e10,
                ],
            }
        )
    return {
        "simultaneous_coverage": 1 - ALPHA,
        "binomial_parameters": parameters,
        "condition": "Original spectral filter passes, with every meter refusal included",
        "regions": rows,
    }


def read(raw):
    shots = raw["shots"]
    samples = []
    for launch in raw["launches"]:
        for mode, panels in launch["panels"].items():
            for profile, panel in zip(raw["profiles"], panels, strict=True):
                xs = [v for v in panel["complete_records"] if v["axis"] == "X"]
                check(
                    [v["multiplier"] for v in xs] == [1, 2], "missing X clock setting"
                )
                samples.append(
                    (launch["owner"], mode, profile["center"], profile["radius"], xs)
                )
    standards = [s for s in raw["charged_standards"] if s["axis"] == "X"]
    families = len(standards) + sum(len(s[-1]) for s in samples)
    check(families == 1086, "incomplete launch or calibration family")
    penalty = log(2 * families / ALPHA)
    reference = calibration(standards, shots, penalty)
    derived = []
    for owner, mode, center, radius, xs in samples:
        windows = []
        for s in xs:
            counts = s["counts"]
            check(
                len(counts) == 14 and sum(counts) == shots,
                "missing meter or refusal counts",
            )
            n = sum(counts[4:13])
            windows.append(
                charge_interval(
                    interval(counts[6], n, penalty), reference[s["multiplier"]]
                )
            )
        lo = max(w[0] for w in windows)
        hi = min(w[1] for w in windows)
        check(lo <= hi, "the two clock settings give incompatible charge intervals")
        derived.append(
            dict(
                owner=owner,
                control=mode,
                center=center,
                radius=radius,
                interval=[floor(lo * 1e10) / 1e10, ceil(hi * 1e10) / 1e10],
            )
        )
    return {
        "simultaneous_coverage": 1 - ALPHA,
        "binomial_parameters": families,
        "original_launches_per_count_bank": shots,
        "eligible_event": "Native source spectral filter passes and both meter sites are unit-filled",
        "conditional_charge_event": "Absolute signed load inside the specified region equals one",
        "estimator": "Two-sided Chernoff intervals on X dark fractions, inverted against the two charged standards and the host at each dwell. Both dwell bounds are intersected.",
        "meter_dark_calibration": reference,
        "regions": derived,
        "full_band_at_its_launch_cell": original_cells(samples, standards, shots),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    raw = SOURCE.read_bytes()
    result = read(json.loads(raw))
    result["counts_file_sha256"] = sha256(raw).hexdigest()
    encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.write:
        RECEIPT.write_text(encoded)
    else:
        check(
            encoded == RECEIPT.read_text(),
            "regional charge bounds differ from the deposit",
        )
        print("Verified simultaneous regional-charge intervals from 1,086 count banks.")


if __name__ == "__main__":
    main()
