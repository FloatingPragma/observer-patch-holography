"""Bound the measured thermal contribution across reporting rates and register sizes."""

import argparse
import json
from hashlib import sha256

import numpy as np

from . import sector_cadence
from . import sector_register_balls as balls
from . import sector_register_pairs as pairs

HERE = balls.HERE


def thermal_range():
    thermal = pairs.chart.source.thermal
    words, _, _, blocks, _ = thermal.apparatus()
    cal = json.loads((HERE / "sector_calorimetry.json").read_text())
    bands = pairs.chart.source.sector_regional.calibrate(
        cal["ordinary_clock_calibration"]
    )
    e, r = bands["exchanged"], bands["retained"]
    # Each sector has exactly the same normalized Gibbs operator at both
    # references. Their ratio is constant throughout the sector, so every
    # calibrated diagonal within the sector bounds this same scalar.
    intervals = np.array(
        [
            [
                np.log(e[g, 0] / r[g, 1]).max() - 1e-12,
                np.log(e[g, 1] / r[g, 0]).min() + 1e-12,
            ]
            for g in blocks
        ]
    )
    bath = thermal.bath()
    for g, (lo, hi) in zip(blocks, intervals, strict=True):
        ix = np.ix_(g, g)
        R, E = bath["retained"][ix], bath["exchanged"][ix]
        ratio = np.trace(E) / np.trace(R)
        assert lo <= np.log(ratio) <= hi
        np.testing.assert_allclose(E, ratio * R, atol=1e-13)
    span = float(intervals[:, 1].max() - intervals[:, 0].min())
    std = json.loads((HERE / "thermal_standard.json").read_text())
    cost = np.array(
        [0, std["single"]["energy_band"][0], std["double"]["energy_band"][0]]
    )
    minimum = float(cost[words[:, :4]].sum(1).min())
    return {
        "sector_count": len(blocks),
        "returned_words": len(words),
        "sector_log_likelihood_intervals": intervals.tolist(),
        "per_member_price_range_upper_nats": span,
        "minimum_ordinary_inventory_energy_per_member": minimum,
        "beta": thermal.BETA,
        "absolute_extra_to_ordinary_upper": span / (thermal.BETA * minimum),
    }


def conditional_history(sub, mats, flags, cat):
    clock = pairs.chart.source.experiment.ReturnClock(
        (mats[2] @ mats[1] @ mats[0]).T.tocsr(), sub
    )
    power, eta = sector_cadence.mixing(clock, len(flags))
    center = flags.mean(0)
    ticks = flags * clock.grid
    history, errors = [], []
    err = np.zeros(flags.shape[1])
    for _ in range(1024):
        ticks, step = clock.advance(ticks)
        err = err + step
        history.append(
            pairs.collective.sector_fronts_counts.group_readings(
                ticks / clock.grid - center, cat
            )
        )
        errors.append(err + 2e-12)
    series = np.array(history)
    errors = np.array(errors)
    # Under the uniform stationary law each centered column has zero mean.
    # Its sup norm is at most its oscillation, contracting by eta per power.
    tail = (np.ptp(ticks / clock.grid, axis=0) + 2 * err + 1e-10) * power / (1 - eta)
    pos, neg = (
        np.maximum(series, 0).reshape(len(series), -1),
        np.minimum(series, 0).reshape(len(series), -1),
    )
    return series, errors, tail, pos, neg, power, eta


def rate_response():
    (_, _, _, _, cfg, ids, _), refs, _, delta, eps = pairs.apparatus()
    size, bank = pairs.get_bank()
    counts, gamma = pairs.collective.sector_fronts_counts.moments(bank)
    cat = pairs.collective.peer_population.fibers.apparatus()[3]
    _, M = pairs.type_weights(cat, counts, gamma, size)
    M /= size - 1
    mabs = abs(M)
    signed = (cfg[:, None, :] + cfg[None, :, :]) @ np.array([1, 3, 9, 27])
    prices = {
        k: balls.prices(signed.ravel(), *refs[k]) for k in ("retained", "exchanged")
    }
    g, kr, lr, hr = prices["retained"]
    ge, ke, le, he = prices["exchanged"]
    assert np.array_equal(g, ge)
    K = (kr - ke)[g].reshape(signed.shape)
    E = np.maximum(kr - ke - (lr - he), (hr - le) - (kr - ke))[g].reshape(signed.shape)
    _, _, sub, _, mats = pairs.chart.source.sources((8,))
    series, errors, tail, pos, neg, power, eta = conditional_history(
        sub, mats, np.eye(len(cfg), dtype=int)[ids], cat
    )
    previous = json.loads((HERE / "sector_fronts_collective.json").read_text())
    rows = []
    ages = np.arange(len(series))
    for cell in previous["continuous_cadence"]["bins"]:
        a, b = cell["rate"]
        lo, hi = (1 - b) ** ages, (1 - a) ** ages
        width = hi @ errors + (1 - a) ** len(series) * tail + 1e-9
        lower = (lo @ pos + hi @ neg).reshape(series.shape[1:]) - width
        upper = (hi @ pos + lo @ neg).reshape(series.shape[1:]) + width
        x, dx = (lower + upper) / 2, (upper - lower) / 2
        d = x.T @ M @ x
        de = abs(x).T @ mabs @ dx + dx.T @ mabs @ (abs(x) + dx)
        val = np.sum(d * K)
        error = np.sum(abs(d) * E + de * (abs(K) + E)) + 1e-10
        scaled = np.multiply.outer([val - error, val + error], [a, b])
        rows.append(
            {
                "rate": [a, b],
                "nats_per_squared_rate": [val - error, val + error],
                "per_rate": [float(scaled.min()), float(scaled.max())],
            }
        )
    check = float(np.sum(delta["cross"] * K))
    check_err = float(
        np.sum(abs(delta["cross"]) * E) + eps["cross"] * (np.ptp(K) + E.max() * 2)
    )
    for cell in rows:
        if cell["rate"][0] <= 1 / 8 <= cell["rate"][1]:
            qlo, qhi = np.array(cell["nats_per_squared_rate"]) / 64
            assert qlo <= check - check_err and qhi >= check + check_err
    return {
        "bins": rows,
        "register_members": size,
        "members_in_reading": 2,
        "age_terms": len(series),
        "power": power,
        "oscillation_contraction": eta,
        "at_eighth_nats_per_pair": [check - check_err, check + check_err],
    }


def run():
    return {
        "sha256": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in (
                "SECTOR_REGISTER_BALLS.md",
                "thermal_standard.json",
                "sector_calorimetry.json",
                "sector_fronts_collective.json",
            )
        },
        "any_common_observation_of_product_references": thermal_range(),
        "fixed_volume_pooled_home_reading": rate_response(),
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    data = run()
    dst = HERE / "sector_register_source.json"
    if cli.parse_args().write:
        dst.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        pairs.chart.source.balls.radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(dst.read_text()), data
        )
    print(
        data["any_common_observation_of_product_references"][
            "absolute_extra_to_ordinary_upper"
        ],
        data["fixed_volume_pooled_home_reading"]["at_eighth_nats_per_pair"],
    )
