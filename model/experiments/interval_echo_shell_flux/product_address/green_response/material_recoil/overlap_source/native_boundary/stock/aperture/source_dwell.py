"""Resolve the connected source with near-half-turn ordinary-pair interference."""

import argparse
from hashlib import sha256
import json
from math import ceil, log, pi, sqrt

import numpy as np

from . import record_receiver as receiver

contact = receiver.contact
source = receiver.source
HERE = source.HERE
DWELL = 2 * pi / 65


def forecast_mixed(calibration):
    energies = np.asarray(list(calibration.values()))
    assert energies.shape == (4, 2)
    z = np.mean(np.exp(-1j * DWELL * energies.mean(axis=1)))
    # Each displaced phase is within dwell times its energy error on the unit circle.
    radius = DWELL * np.mean(energies[:, 1] - energies[:, 0]) / 2
    distance = abs(1 - z)
    return np.array(
        [max(0.0, distance - radius) ** 2 / 2, (distance + radius) ** 2 / 2]
    )


def run():
    contact_run = json.loads((HERE / "pair_contact.json").read_bytes())
    sparse = json.loads((HERE / "record_receiver.json").read_bytes())
    mixed = forecast_mixed(contact_run["decoded"]["ordinary_energy_bands"])
    probabilities = receiver.sharp_bright(DWELL)
    expected = (1 + np.cos(pi / 65)) ** 2 / 8
    np.testing.assert_allclose(probabilities @ [1, -1, -1, 1], expected)
    assert mixed[0] <= expected <= mixed[1]
    _, words = source.classes()
    panel = [(len(w), source.counts(w)) for w in words]
    cases = {}
    rng = np.random.default_rng(512416)
    for label, menu in (
        ("accepted_2200", "none"),
        ("sparse_blind", "none"),
        ("sparse_third", "third_head"),
    ):
        fixed = label == "accepted_2200"
        weights = (
            [source.Q(k == 2) for k in range(5)]
            if fixed
            else source.input_weights(source.LEVELS[0])
        )
        a, b, _ = source.recover(source.combine(panel, weights, menu))
        denominator = sum(a) if fixed else 1
        if fixed:
            delta = receiver.probe.log_factor(contact_run["source"])
        else:
            delta = sparse["cases"][menu]["delta_from_charge"]
        forecast = receiver.probe.mix_bands(delta, mixed)
        sign = np.sign(forecast[0])
        d = min(sign * forecast)
        assert d > 0
        # Hoeffding on 2N independent 0/1 variables with weights +/-1/N.
        # At threshold d/2 both the null upper tail and alternative lower tail are <= exp(-N*d²/4).
        n = ceil(4 * log(100) / d**2)
        probs = [
            float(np.asarray(v, float) @ probabilities) / float(denominator)
            for v in (a, b)
        ]
        actual = probs[1] - probs[0]
        assert forecast[0] <= actual <= forecast[1]
        counts = [int(rng.binomial(n, p)) for p in probs]
        shift = (counts[1] - counts[0]) / n
        width = sqrt(log(200) / n)
        cases[label] = dict(
            forecast=forecast.tolist(),
            trial_basis="accepted" if fixed else "attempted",
            sufficient_per_arm=n,
            threshold=d / 2,
            sign=int(sign),
            bright_counts=counts,
            observed_shift=shift,
            independent_99_percent=[shift - width, shift + width],
            detected=bool(sign * shift > d / 2),
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "SOURCE_DWELL.md",
                HERE / "pair_contact.json",
                HERE / "record_receiver.json",
            )
        },
        dwell=DWELL,
        mixed_forecast=mixed.tolist(),
        mixed_exact=expected,
        cases=cases,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    path = HERE / "source_dwell.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), result
        )
    print(result["cases"])
