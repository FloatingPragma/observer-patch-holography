"""Test finite weak-clock cadence in the family mass comparison."""

import argparse
import json
from pathlib import Path

import mpmath as mp
import numpy as np

SOURCE = (
    "https://github.com/FloatingPragma/observer-patch-holography/blob/"
    "3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/lepton_register_round/PROOF.md"
)


def moment(n, harmonic):
    p = mp.mpf(2) / 9
    return (1 + p * mp.expm1(mp.j * harmonic / n)) ** n


def levels(n, read):
    c1 = moment(n, 1)
    roots = [mp.expj(2 * mp.pi * j / 3) for j in range(3)]
    if read == "all_records":
        c2 = moment(n, 2)
        return sorted(
            2 + 2 * mp.sqrt(2) * (root * c1).real + (root**2 * c2).real
            for root in roots
        )
    if read == "phase_only":
        c1 /= abs(c1)
    elif read != "all_returns_joint":
        raise ValueError(read)
    return sorted((1 + mp.sqrt(2) * (root * c1).real) ** 2 for root in roots)


def explicit_check(n):
    """Sum binomial clock histories, averaging the three-path probabilities or amplitudes."""
    p = mp.mpf(2) / 9
    weights = [mp.binomial(n, k) * p**k * (1 - p) ** (n - k) for k in range(n + 1)]
    mean, second = [], []
    for j in range(3):
        a = [
            1 + mp.sqrt(2) * mp.cos(mp.mpf(k) / n + 2 * mp.pi * j / 3)
            for k in range(n + 1)
        ]
        mean.append(sum(w * v for w, v in zip(weights, a)) ** 2)
        second.append(sum(w * v**2 for w, v in zip(weights, a)))
    for measured, expected in (
        (sorted(mean), levels(n, "all_returns_joint")),
        (sorted(second), levels(n, "all_records")),
    ):
        if max(abs(a - b) for a, b in zip(measured, expected)) > mp.mpf("1e-45"):
            raise ValueError(
                "Finite clock history disagrees with the harmonic compression"
            )


def run():
    with mp.workdps(60):
        p = mp.mpf(2) / 9
        e = mp.mpf("0.51099895069")
        target = mp.mpf("206.7682827")
        uncertainty = mp.mpf("0.0000046")
        for n in (1, 10, 100, 302):
            explicit_check(n)
        results = {}
        rng = np.random.default_rng(178013)
        for read in ("phase_only", "all_returns_joint", "all_records"):

            def ratio(n):
                v = levels(n, read)
                return v[1] / v[0]

            def crossing(goal):
                upper = 2
                while ratio(upper) < goal and upper < 2**40:
                    upper *= 2
                if ratio(upper) < goal:
                    raise ValueError("No crossing in the search window")
                return mp.findroot(lambda n: ratio(n) - goal, (upper / 2, upper))

            root = crossing(target)
            fitted = int(mp.nint(root))
            tested = {}
            for n in sorted({1, 10, 100, 1000, 301, 302, 303, fitted}):
                v = levels(n, read)
                tested[str(n)] = {
                    "middle_ratio": float(v[1] / v[0]),
                    "tau_MeV_electron_anchor": float(e * v[2] / v[0]),
                }
            # Independent clock observable on the comparison with the fitted input supplied.
            c1 = moment(fitted, 1)
            dark = (1 - (c1 * mp.expj(-p)).real) / 2
            shots = 2**36
            results[read] = {
                "N_fitted_from_muon": fitted,
                "real_root_using_muon": float(root),
                "integer_window_near_fit_within_one_ratio_uncertainty": [
                    int(mp.ceil(crossing(target - uncertainty))),
                    int(mp.floor(crossing(target + uncertainty))),
                ],
                "at_N": tested,
                "clock_visibility_at_fit": float(abs(c1)),
                "clock_dark_probability_at_fit": float(dark),
                "clock_dark_counts_given_fitted_N": int(
                    rng.binomial(shots, float(dark))
                ),
            }
        return {
            "source": SOURCE,
            "separate_clock_shots": shots,
            "seed": 178013,
            "explicit_binomial_history_checks": [1, 10, 100, 302],
            "cadence_is_an_input": True,
            "finite_register_comparisons": results,
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x") as stream:
        json.dump(run(), stream, indent=2, allow_nan=False)
        stream.write("\n")
