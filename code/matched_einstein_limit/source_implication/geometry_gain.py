"""Test an independently surveyed equal-volume shape against ordinary and pairing changes."""

import json
from itertools import combinations
from pathlib import Path

import numpy as np

from model.experiments.interval_echo_einstein_station.collar_pairing import (
    measure as pairing,
)
from model.experiments.interval_echo_einstein_station.collar_pairing.ranging import (
    measure as ranging,
)
from model.experiments.interval_echo_einstein_station.collar_pairing.ranging.pool import (
    pooled_ranges,
    distance_check,
)


def simplex_area_volume(ranges):
    """Fit clock gains and reconstruct total area and interior volume from the distances."""
    off = ~np.eye(4, dtype=bool)
    logs = np.zeros_like(ranges)
    logs[off] = np.log(ranges[off])
    gains = (logs - logs.T).mean(axis=1)
    corrected = logs - gains[:, None]
    sq = np.exp(corrected + corrected.T)
    np.fill_diagonal(sq, 0)

    def simplex(sites):
        d = sq[np.ix_(sites, sites)]
        gram = (d[0, 1:, None] + d[None, 0, 1:] - d[1:, 1:]) / 2
        return np.sqrt(np.linalg.det(gram)) / (6 if len(sites) == 4 else 2)

    volume = simplex(range(4))
    return sum(simplex(face) for face in combinations(range(4), 3)), volume


def equal_volume_area(ranges):
    area, volume = simplex_area_volume(ranges)
    return area / volume ** (2 / 3)


def area_gradient(ranges, difference):
    """Complex-step differentiation inside the positive, nondegenerate simplex domain."""
    return float(equal_volume_area(ranges + 1j * 1e-16 * difference).imag / 1e-16)


def returns():
    board = pairing.formation.federation(1)
    route = pairing.formation.region(board)
    joint, independent, _ = pairing.prepared_states(board, route, 1)
    joint /= np.trace(joint).real
    independent /= np.trace(independent).real
    results = []
    for odd in (False, True):
        h = ranging.generator(board, route, 1, 8, odd)
        (arm,) = pairing.material_read(h, route[1:3], {1: (joint, independent)})
        results.append(
            [
                arm[key]["appointments"][-1]["complete"]
                for key in ("formed", "independent")
            ]
        )
    even, odd = results
    return float(even[0] - even[1]), float(odd[1] - even[1])


def run():
    archive = json.loads(Path(ranging.__file__).with_name("results.json").read_text())
    banks = {
        (arm["probe"], arm["left_odd"]): arm["panels"]
        for arm in archive["runs"]
        if arm["contrast"] == 8
    }
    z_dark, z_ordinary = returns()
    panels = []
    for i in (0, 1):
        for horizon in (2, 4, 6):
            even = [banks[sign, False][i] for sign in (-1, 1)]
            odd = [banks[sign, True][i] for sign in (-1, 1)]
            baseline = pooled_ranges(even, "independent", horizon)
            domain = distance_check(baseline)
            if min(domain["centered_eigenvalues"]) <= 0:
                raise ArithmeticError(
                    "the projected distances fail the tetrahedral domain"
                )
            dark_change = pooled_ranges(even, "formed", horizon) - baseline
            ordinary_change = pooled_ranges(odd, "independent", horizon) - baseline
            a_dark = area_gradient(baseline, dark_change)
            a_ordinary = area_gradient(baseline, ordinary_change)
            forecast = z_ordinary * a_dark / a_ordinary
            panels.append(
                dict(
                    tick=even[0]["step"],
                    horizon=horizon,
                    area=equal_volume_area(baseline),
                    ordinary_area_derivative=a_ordinary,
                    pairing_area_derivative=a_dark,
                    ordinary_balanced_bit_forecast_for_pairing=forecast,
                    relative_transfer=forecast / z_dark,
                    determinant=z_ordinary * a_dark - z_dark * a_ordinary,
                    sign_of_ordinary_area_gain=int(
                        np.sign(z_ordinary * z_dark / a_ordinary)
                    ),
                )
            )
    return dict(
        pairing_completion=z_dark,
        ordinary_completion_per_exchange_energy=z_ordinary,
        panels=panels,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
