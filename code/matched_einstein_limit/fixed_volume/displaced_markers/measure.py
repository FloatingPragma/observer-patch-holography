"""Read the leading area response from displaced marker cohorts with unknown local intensity."""
from itertools import product
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.stats import chi2


HERE = Path(__file__).resolve().parent
SEED = 632871


def shifted_count_mean(curvature, radius, density_quadratic=0):
    """Normalized mark count in a half-radius ball displaced by one radius."""
    with mp.workdps(40):
        d = mp.mpf(radius)
        s = d / 2
        k = mp.mpf(curvature)
        if k:
            freq = mp.sqrt(abs(k))
            sin = mp.sin if k > 0 else mp.sinh
            cos = mp.cos if k > 0 else mp.cosh

        def integrand(r):
            if k == 0:
                cap = (s**2 - (r - d)**2) / (2 * r * d)
            else:
                cap = 1 - (cos(freq * s) - cos(freq * r) * cos(freq * d)) / (
                    mp.sign(k) * sin(freq * r) * sin(freq * d))
            return r**2 * (1 + density_quadratic * r**2) * cap

        return float(3 * mp.quad(integrand, [d - s, d + s]) / (2 * s**3))


def decode(means, radius):
    central, half, *displaced = means
    if len(displaced) != 6:
        raise ValueError("Six displaced axis reads are required")
    density = (4 * half - central) / 3
    product_estimate = (18 * np.mean(displaced) + 22 * half - 40 * central) / radius**2
    return float(product_estimate / density)


def count_interval(count, budget, tail):
    # Two-sided Garwood limits for the independently repeated Poisson census.
    lower = chi2.ppf(tail, 2 * count) / 2 if count else 0
    upper = chi2.isf(tail, 2 * (count + 1)) / 2
    return lower / budget, upper / budget


def decoded_interval(counts, budget, radius, tail):
    bounds = [count_interval(int(x), budget, tail) for x in counts]
    central, half, *displaced = bounds
    # The two centered observations each enter the numerator and denominator.
    # Taking an interval quotient keeps their possible dependence.
    density = ((4 * half[0] - central[1]) / 3, (4 * half[1] - central[0]) / 3)
    if density[0] <= 0:
        raise ValueError("The data do not resolve a positive reference intensity")
    numerator = ((18 * np.mean([v[0] for v in displaced]) + 22 * half[0] - 40 * central[1]),
                 (18 * np.mean([v[1] for v in displaced]) + 22 * half[1] - 40 * central[0]))
    candidates = [n / (c * radius**2) for n, c in product(numerator, density)]
    return [float(min(candidates)), float(max(candidates))]


def measure():
    archive = json.loads((HERE.parent / "measurements.json").read_text())
    rng = np.random.default_rng(SEED)
    number_of_cases = sum(len(panel["cases"]) + 1 for panel in archive["panels"])
    tail = .01 / (2 * 8 * number_of_cases)
    rows = []
    for panel in archive["panels"]:
        radius = panel["flat_radius"]
        budget = int(2**17 / radius**6)
        cases = [(entry, 0.) for entry in panel["cases"]]
        flat = next(case for case in panel["cases"] if case["prepared_sectional_curvature"] == 0)
        cases.append((flat, .25))
        for case, quadratic in cases:
            curvature = case["prepared_sectional_curvature"]
            means = [1 + 3 * quadratic * s**2 / 5 for s in (radius, radius / 2)]
            means += [shifted_count_mean(curvature, radius, quadratic)] * 6
            counts = rng.poisson(np.array(means) * budget)
            estimate = decode(counts / budget, radius)
            area_contrast = case["inferred"]["volume_threshold"]["area_contrast"]
            # This independent, uniform proper-area census never enters the decoder.
            area_scalar = -15 * area_contrast / (2 * math.pi * radius**4)
            rows.append(dict(radius=radius, sectional_curvature=curvature,
                density_quadratic=quadratic, expected_normalized_counts=means,
                matched_area_count_scalar=area_scalar, area_contrast=area_contrast,
                normalized_count_budget=budget, counts=list(map(int, counts)),
                reconstructed_scalar=estimate, mean_decoder=decode(means, radius),
                mean_center_only=-40 * (means[0] - means[1]) / radius**2,
                counting_interval=decoded_interval(counts, budget, radius, tail)))
    return dict(seed=SEED, simultaneous_counting_probability=.99, rows=rows)


if __name__ == "__main__":
    print(json.dumps(measure(), indent=2, allow_nan=False))
