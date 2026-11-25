"""Measure a selected charged line's standing spatial phase using only local loads."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta as beta_distribution

from . import measure

HERE = Path(__file__).resolve().parent
# A four-station contrast chosen on the modeled profile, before the fresh counts.
STATIONS = (3, 4, 5, 6)
SHOTS = 2**38
SEED = 738193


def interval(count, n, tail):
    """A two-sided exact binomial interval, including boundary outcomes."""
    low = 0 if count == 0 else beta_distribution.ppf(tail, count, n - count + 1)
    high = 1 if count == n else beta_distribution.ppf(1 - tail, count + 1, n - count)
    return np.array([low, high])


def contrast(bounds, coefficients):
    low = sum(c * bounds[j, int(c < 0)] for j, c in enumerate(coefficients))
    high = sum(c * bounds[j, int(c > 0)] for j, c in enumerate(coefficients))
    return np.array([low, high])


def read_counts(counts, launches):
    """Infer cos(2k) up to its spatial alias; neither energy nor generator is an input."""
    bounds = np.array([interval(x, launches, 0.01 / 8) for x in counts])
    denominator = contrast(bounds, (0, -2, 2, 0))
    if denominator[0] <= 0:
        raise ValueError("The positive density contrast is unresolved.")
    numerator = contrast(bounds, (-1, 1, -1, 1))
    ratios = np.outer(numerator, 1 / denominator)
    cosine = np.array([ratios.min(), ratios.max()])
    if cosine[0] < -1 or cosine[1] > 1:
        raise ValueError("The counts do not resolve two distinct spatial aliases.")
    lower = np.arccos(cosine[::-1]) / 2
    return dict(
        cosine_of_double_phase_99=cosine.tolist(),
        phase_intervals_99=[lower.tolist(), (np.pi - lower[::-1]).tolist()],
        zero_phase_third_difference_99=contrast(bounds, (-1, 3, -3, 1)).tolist(),
    )


def run():
    records = json.loads((HERE / "results.json").read_text())
    energies, _, input_weights, lift, modes, native, _ = measure.build()
    clock = records["clock"]
    x = (energies - clock["center"]) * clock["tick"] / (2 * np.pi)
    number_delays = clock["number_delays"]
    response = np.exp(-1j * np.pi * (number_delays - 1) * x)
    response *= np.sinc(number_delays * x) / np.sinc(x)
    state = lift @ (modes @ (input_weights * response))
    fixture = measure.fixture()
    port_index = {p: j for j, p in enumerate(fixture["ports"])}
    station_ports = [[port_index[leg[j]] for leg in fixture["legs"]] for j in STATIONS]
    # Energy enters the preparation only. A detector reads the sharp charge,
    # summing over all 64 words of the six terminal head records.
    position_probability = (abs(state.reshape(-1, 64)) ** 2).sum(axis=1)
    probabilities = np.array([position_probability[s].sum() for s in station_ports]) / 3
    draws = np.random.default_rng(SEED).multinomial(
        SHOTS, [*probabilities, 1 - probabilities.sum()]
    )
    # Verify the uniform free-arm recurrence directly on all 64 incident columns.
    for leg in fixture["legs"]:
        for j in STATIONS:
            k = port_index[leg[j]]
            for word in range(64):
                col = native.getrow(64 * k + word)
                terms = dict(zip(col.indices, col.data, strict=True))
                assert terms == {
                    64 * k + word: 17,
                    64 * port_index[leg[j - 1]] + word: -0.5,
                    64 * port_index[leg[j + 1]] + word: -0.5,
                }
    return dict(
        launches=SHOTS,
        seed=SEED,
        local_port_groups=[
            [fixture["ports"][k] for k in group] for group in station_ports
        ],
        station_indices=list(STATIONS),
        original_start_counts=draws.tolist(),
        encoder_success=1 / 3,
        energy_selected_probability=float(np.vdot(state, state).real),
        inferred=read_counts(draws[:-1], SHOTS),
        selected_pole=float(energies[records["selected_standard_index"]]),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    actual = run()
    p = HERE / "wavelength.json"
    if args.write:
        p.write_text(json.dumps(actual, indent=2) + "\n")
        return
    saved = json.loads(p.read_text())
    for key in (
        "launches",
        "seed",
        "local_port_groups",
        "station_indices",
        "original_start_counts",
    ):
        if actual[key] != saved[key]:
            raise SystemExit(f"Wavelength receipt differs at {key}.")
    for key in ("encoder_success", "energy_selected_probability", "selected_pole"):
        np.testing.assert_allclose(actual[key], saved[key], rtol=1e-9)
    for key in actual["inferred"]:
        np.testing.assert_allclose(
            actual["inferred"][key], saved["inferred"][key], rtol=1e-9
        )


if __name__ == "__main__":
    main()
