"""A necessary common-time identity for two independently prepared gap clocks."""

import json
from pathlib import Path

import numpy as np
from model.experiments.count_intervals import clopper_pearson

HERE = Path(__file__).resolve().parent


def chebyshev_image(success, shots, tail, order):
    lo, hi = clopper_pearson(success, shots, tail)
    x, y = 2 * order * np.arcsin(np.sqrt([lo, hi]))
    lower, upper = np.minimum(np.cos(x), np.cos(y)), np.maximum(np.cos(x), np.cos(y))
    # Interior extrema occur at consecutive integer multiples of pi.
    for n in range(order + 1):
        crossed = (x <= n * np.pi) & (n * np.pi <= y)
        if n % 2:
            lower[crossed] = -1
        else:
            upper[crossed] = 1
    return np.column_stack([lower, upper])


def audit():
    with np.load(HERE / "counts.npz") as counts:
        summary = json.loads((HERE / "results.json").read_text())
        panels = summary["panels"]
        names = [key.removeprefix("k8_even_") for key in panels if key.startswith("k8_even_")]
        dates = counts["k8_even_blank"].shape[0]
        n = summary["shots"]
        tests = 2 * dates * len(names)
        result = {"simultaneous_coverage": .99, "probability_intervals": tests, "settings": {}}
        for name in names:
            images = []
            for arm, order in (("even", 33), ("odd", 32)):
                count = counts[f"k8_{arm}_{name}"][:, 2:].sum(axis=1)
                images.append(chebyshev_image(count, n, .01 / (2 * tests), order))
            left, right = images
            gap = np.maximum(left[:, 0] - right[:, 1], right[:, 0] - left[:, 1])
            peak = int(gap.argmax())
            result["settings"][name] = {
                "disjoint_appointments": int(sum(gap > 0)),
                "max_gap": float(max(0, gap[peak])), "appointment": peak if gap[peak] > 0 else None,
                "even_image_at_max_gap": left[peak].tolist() if gap[peak] > 0 else None,
                "odd_image_at_max_gap": right[peak].tolist() if gap[peak] > 0 else None,
            }
    return result


if __name__ == "__main__":
    import argparse
    from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = audit()
    target = HERE / "common_time.json"
    if args.write:
        target.write_text(json.dumps(result, indent=2) + "\n")
    else:
        compare(result, json.loads(target.read_text()), "common clock time")
    for key, row in result["settings"].items():
        print(key, row)
