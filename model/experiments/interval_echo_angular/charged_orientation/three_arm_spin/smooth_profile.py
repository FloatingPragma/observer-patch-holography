"""Forecast held charge counters from a standing wave inferred on four stations."""

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial

from . import arm_wavelength as calibration
from . import measure

HERE = Path(__file__).resolve().parent
HELD = tuple(range(7, 15))
STARTS = 2**38
SEED = 6834157


def outward(low, high):
    return np.array([np.nextafter(low, -np.inf), np.nextafter(high, np.inf)])


def product(left, right):
    corners = np.outer(left, right)
    return outward(corners.min(), corners.max())


def plus(left, right):
    return outward(left[0] + right[0], left[1] + right[1])


def polynomial_range(polynomial, interval):
    value = np.zeros(2)
    for coefficient in reversed(polynomial.coef):
        value = plus(product(value, interval), [coefficient, coefficient])
    return value


def weights():
    # At indices 4, 5 and 6 the three interpolation weights form the identity.
    prior = [[Polynomial([x]) for x in row] for row in np.eye(3)]
    multiplier = Polynomial([1, 2])
    for _ in HELD:
        first, second, third = prior[-3:]
        prior.append(
            [
                a + multiplier * (c - b)
                for a, b, c in zip(first, second, third, strict=True)
            ]
        )
    return prior[3:]


def prediction(counts, starts):
    tail = 0.01 / (2 * (4 + len(HELD)))
    training = np.array([calibration.interval(v, starts, tail) for v in counts])
    numerator = calibration.contrast(training, (-1, 1, -1, 1))
    denominator = calibration.contrast(training, (0, -2, 2, 0))
    if denominator[0] <= 0:
        raise ValueError("Calibration does not resolve a positive difference.")
    cosine = product(numerator, outward(1 / denominator[1], 1 / denominator[0]))
    if cosine[0] <= -1 or cosine[1] >= 1:
        raise ValueError("Calibration does not isolate a standing spatial frequency.")
    partition = np.linspace(*cosine, num=65)
    projected = []
    for row in weights():
        cells = []
        for a, b in zip(partition[:-1], partition[1:], strict=True):
            result = np.zeros(2)
            for polynomial, prob in zip(row, training[1:], strict=True):
                result = plus(
                    result, product(polynomial_range(polynomial, [a, b]), prob)
                )
            cells.append(result)
        projected.append([max(0.0, min(x[0] for x in cells)), max(x[1] for x in cells)])
    flat_range = [training[:, 0].max(), training[:, 1].min()]
    third_difference = calibration.contrast(training, (-1, 3, -3, 1))
    return {
        "training_intervals": training.tolist(),
        "cosine_of_twice_phase": cosine.tolist(),
        "phase_per_step_aliases": [
            (np.arccos(cosine[::-1]) / 2).tolist(),
            (np.pi - np.arccos(cosine) / 2).tolist(),
        ],
        "predicted": projected,
        "flat_compatible_with_training": bool(flat_range[0] <= flat_range[1]),
        "third_difference_training": third_difference.tolist(),
        "quadratic_compatible_with_training": bool(
            third_difference[0] <= 0 <= third_difference[1]
        ),
    }


def held_counts():
    """Independently sample the finite-clock state on the registered sparse graph."""
    setup = json.loads((HERE / "results.json").read_text())
    clock = setup["clock"]
    energies, _, input_weights, lift, modes, _, _ = measure.build()
    phases = (energies - clock["center"]) * clock["tick"] / (2 * np.pi)
    delays = clock["number_delays"]
    exposure = np.exp(-1j * np.pi * (delays - 1) * phases)
    exposure *= np.sinc(delays * phases) / np.sinc(phases)
    state = lift @ (modes @ (input_weights * exposure))
    physical = measure.fixture()
    port_index = {p: j for j, p in enumerate(physical["ports"])}
    sites = [[leg[j] for leg in physical["legs"]] for j in HELD]
    density = np.sum(abs(state.reshape(-1, 64)) ** 2, axis=1) / 3
    probs = [sum(density[port_index[p]] for p in group) for group in sites]
    reads = np.random.default_rng(SEED).multinomial(STARTS, [*probs, 1 - sum(probs)])
    return sites, reads.tolist()


def run():
    archived = json.loads((HERE / "wavelength.json").read_text())
    assert archived["station_indices"] == [3, 4, 5, 6]
    counts = archived["original_start_counts"][:4]
    fit = prediction(counts, archived["launches"])
    sites, reads = held_counts()
    intervals = np.array(
        [calibration.interval(x, STARTS, 0.01 / 24) for x in reads[:-1]]
    )
    compared = []
    for bounds, observed in zip(fit["predicted"], intervals, strict=True):
        compared.append(bool(bounds[0] <= observed[1] and observed[0] <= bounds[1]))
    return {
        "training": {
            "source": "wavelength.json",
            "stations": archived["station_indices"],
            "original_starts": archived["launches"],
            "counts": counts,
        },
        "count_only_inference": fit,
        "held": {
            "stations": HELD,
            "physical_port_groups": sites,
            "original_starts": STARTS,
            "seed": SEED,
            "counts_with_other_final": reads,
            "intervals": intervals.tolist(),
            "forecast_overlaps": compared,
        },
        "simultaneous_confidence": 0.99,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    outcome = run()
    if args.write:
        (HERE / "smooth_profile.json").write_text(json.dumps(outcome, indent=2) + "\n")
    print(
        "Held counter intervals covered:",
        sum(outcome["held"]["forecast_overlaps"]),
        "/",
        len(HELD),
    )
