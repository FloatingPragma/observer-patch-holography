"""Resolve rare recoil using independent time-domain absorption and native receiver shots."""

import argparse
import json
import numpy as np

from . import source, transport
from .. import measure as apparatus
from ..closed_region.measure import check


def terminal(origin, alpha):
    live = np.zeros((10, 10, 10))
    live[origin] = 1
    deposit = np.zeros_like(live)
    for age in range(4096):
        deposit += alpha * live
        padded = np.pad(live, 1)
        next_live = 2 * live
        for axis in range(3):
            for shift in (1, -1):
                section = [slice(1, -1)] * 3
                section[axis] = slice(2, None) if shift == 1 else slice(None, -2)
                next_live = next_live + padded[tuple(section)]
        live = (1 - alpha) * next_live / 8
        if live.sum() < 1e-17:
            break
    assert live.sum() < 1e-17
    return deposit


def run():
    n = 10
    steps, conductance = apparatus.green.measured_steps(None)
    origins = transport.triangle(n, 0)
    alpha = 1 / 64
    modes, eig = transport.spectral(n, conductance[0])
    rng = np.random.default_rng(271130)
    trials = 2**41
    waves = np.zeros((n,) * 3, dtype=np.int64)
    forecast = np.zeros(waves.shape)
    source_counts = []
    for origin in origins:
        p = terminal(origin, alpha).ravel()
        mass = p.sum()
        pairs = p * p
        probs = np.r_[
            pairs, mass * mass - pairs.sum(), 2 * mass * (1 - mass), (1 - mass) ** 2
        ]
        sample = rng.multinomial(trials, probs)
        assert int(sample.sum()) == trials
        only = rng.binomial(sample[:-3], 1 / 4)
        high = rng.binomial(only, 1 / 2)
        emitted = rng.binomial(high, 9 / 512) + rng.binomial(only - high, 8 / 512)
        waves += emitted.reshape(waves.shape)
        source_counts.append(
            dict(
                origin=origin,
                coincident=int(sample[:-3].sum()),
                separate=int(sample[-3]),
                one_lost=int(sample[-2]),
                both_lost=int(sample[-1]),
                releases=int(emitted.sum()),
            )
        )
        theoretical = alpha * transport.green(origin, n, modes, eig, alpha)
        np.testing.assert_allclose(p.reshape(waves.shape), theoretical, atol=8e-17)
        forecast += float(source.run_case(8)[1]["sharp"]["gain"]) / 512 * theoretical**2
    release_count = int(waves.sum())
    expected = float(trials * forecast.sum())
    visits = np.zeros_like(waves)
    H = apparatus.receiver.parent.H
    for tick in range(H):
        if not waves.any():
            break
        visits += waves
        waves = transport.propagate(waves, steps, rng)
    assert not waves.any()
    predicted_visits = apparatus.green.poisson(forecast / forecast.sum(), conductance)[
        1:-1, 1:-1, 1:-1
    ]
    panels = []
    for probe in (1, 2):
        setting = apparatus.receiver.Setting("matched", (6, 5, 5, 0), probe, (2, 2, 2))
        arms = apparatus.receiver.count_panel(rng, visits, setting)
        observed = sum(sum(a["right"]) - sum(a["starts"]) / 2 for a in arms) / expected
        reward = apparatus.receiver.signed_readout(setting)
        predicted = float((reward * predicted_visits).sum())
        # Each original write yields either zero or one entire material flight.
        # For the fixed origin strata their summed second moments are bounded
        # by the release-weighted moment used here; use expected release count.
        width = apparatus.receiver.bound(reward, predicted_visits, 2, releases=expected)
        panels.append(
            dict(
                probe=probe,
                arms=arms,
                signed_per_expected_release=observed,
                predicted=predicted,
                simultaneous_99_half_width=width,
            )
        )
    return dict(
        trials_per_origin=trials,
        source_counts=source_counts,
        expected_releases=expected,
        emitted=release_count,
        appointments=int(visits.sum()),
        last_tick=tick,
        panels=panels,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    result = json.loads(json.dumps(run()))
    path = source.HERE / "precision.json"
    if p.parse_args().write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        check(json.loads(path.read_text()), result)
