"""An independent chord ruler and a signed ordinary source on the same flat guide."""

import json
from math import isqrt, pi

import numpy as np

from model.candidate_laws.interval_echo.candidate import KAPPA
from ..collar_identifiability.refinement.admission import material_panel


def integrated_chords(points_per_radius):
    """Return integral(w)/radius^4 using both first exterior stations on each row."""
    n = points_per_radius
    total = 0
    for y in range(1 - n, n):
        for z in range(1 - n, n):
            bound = n * n - y * y - z * z
            if bound <= 0:
                continue
            # Interior x = 1-b,...,b-1. The first exits are the stations at ±b.
            b = isqrt(bound - 1) + 1
            total += (2 * b - 1) * b * b - (b - 1) * b * (2 * b - 1) // 3
    return total / (2 * n**5)


def run():
    bins = material_panel(1, 1, 8)["records"][-1]
    q = np.array(bins["q"])
    p = q + np.array(bins["d"])

    def gain(hit):
        return float(-(q - p) @ np.log(p + hit * (q - p)))

    energy = 6 * KAPPA
    gamma = energy / gain(0)
    rows = []
    for n in (4, 8, 16, 32, 64, 128):
        radius = 1 / (2048 * n)
        spacing = radius / n
        # One hole outside the sampled ball. A normalized uniform mode is exactly
        # stationary on the punctured guide when the independent anchor is vacant.
        unit_occupancy_weight = 2 * pi * integrated_chords(n) / (1 - spacing**3)
        d_hit = gamma * unit_occupancy_weight * radius**4
        hit = 0.5 * d_hit
        # Rounding each exit out by at most one spacing gives this uniform bound.
        peak = gamma * pi * radius * (1 + 1 / n) ** 2
        if peak > 1:
            raise ValueError("the local switch probability could exceed one")
        rows.append(
            dict(
                points_per_radius=n,
                radius=radius,
                spacing=spacing,
                normalized_kernel=unit_occupancy_weight,
                maximum_switch_bound=peak,
                normalized_entropy_derivative=gain(hit) * gamma * unit_occupancy_weight,
                normalized_area_derivative=0.0,
            )
        )
    return dict(
        energy_per_mobile=energy,
        ordinary_entropy_limit=energy * 8 * pi**2 / 15,
        gain=gain(0),
        switch_scale=gamma,
        panels=rows,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
