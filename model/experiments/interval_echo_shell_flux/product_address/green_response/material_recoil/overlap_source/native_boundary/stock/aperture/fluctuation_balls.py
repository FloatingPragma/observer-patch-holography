"""Signed home-load fluctuations on nested regions of the common material echo chart."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np

from . import closed_energy as energy
from . import thermal_radar as radar
from .measure import check

HERE = energy.HERE
SHOTS = 2**26


def run():
    _, w, cells, *_ = energy.pair.setup()
    board, home, seams, *_ = energy.returned(w)
    _, distance, _ = radar.chart(board, seams)
    d = w[:, :4] - 1
    np.testing.assert_array_equal(d.sum(1), 0)
    nonunit = d[:, 1] != 0
    np.testing.assert_array_equal(d[nonunit, 0], -d[nonunit, 1])
    radii = np.arange(distance.max() + 1)
    mask = (
        (distance[:, home, None] <= radii[None, None, :])
        * (1 << np.arange(4))[None, :, None]
    ).sum(1)
    assert np.all(mask[:, -1] == 15)
    delta = 4 * np.sqrt(np.log(2 * 32 / 0.01) / (2 * SHOTS))
    result = {}
    additive = {}
    totals = []
    bands = []
    for seed, (name, weights) in enumerate(
        (
            ("original", np.full(len(w), F(1), object)),
            ("repaired", energy.least_change(cells)),
        )
    ):
        weights = weights / sum(weights)
        assert all(weights @ d == 0)
        cov = d.T @ (weights[:, None] * d)
        assert all(cov[i, j] <= 0 for i in range(4) for j in range(i))
        draw = np.random.default_rng(843007 + seed).multinomial(
            SHOTS, np.asarray(weights, float)
        )
        exact = []
        bounds = []
        for bits in range(16):
            inside = [i for i in range(4) if bits >> i & 1]
            outside = [i for i in range(4) if i not in inside]
            signed = d[:, inside].sum(1)
            var = weights @ (signed**2)
            cut = -sum(cov[i, j] for i in inside for j in outside)
            assert var == cut
            local = sum(cov[i, i] for i in inside)
            internal = sum(
                (cov[i, j] for i in inside for j in inside if i < j), start=F()
            )
            assert var - 2 * internal == local
            variance = draw @ (signed**2) / SHOTS
            # Identically zero on the entire support for the empty and full cuts.
            lo, hi = (
                (0, 0)
                if not np.any(signed)
                else (max(0, variance - delta), variance + delta)
            )
            assert lo <= float(var) <= hi
            exact.append(str(var))
            bounds.append([lo, hi])
        centers = np.sqrt(np.array([float(F(x)) for x in exact]))[mask].sum(1)
        band = np.sqrt(np.array(bounds))[mask].sum(1)
        totals.append(centers)
        bands.append(band)
        additive[name] = dict(
            variance_units=str(sum(cov.diagonal())),
            rms_units=sum(np.sqrt(float(v)) for v in cov.diagonal()),
        )
        result[name] = dict(
            variances=exact,
            measured_variance_bands=bounds,
            minimum=float(centers.min()),
            maximum=float(centers.max()),
            middle_total=float(centers[home[1]]),
            middle_interval=band[home[1]].tolist(),
        )
    zero = np.flatnonzero(np.ptp(distance[:, home], axis=1) == 0)
    difference = totals[1] - totals[0]
    band = np.column_stack(
        (bands[1][:, 0] - bands[0][:, 1], bands[1][:, 1] - bands[0][:, 0])
    )
    restfile = HERE / "thermal_standard.json"
    rest = json.loads(restfile.read_bytes())
    mass_ratio = F(float(rest["single"]["energy_band"][1])) / F(
        float(rest["double"]["energy_band"][0])
    )
    factor = F(2, 9) * F(448, 2110) * mass_ratio
    dv = F(additive["repaired"]["variance_units"]) - F(
        additive["original"]["variance_units"]
    )
    ds = additive["repaired"]["rms_units"] - additive["original"]["rms_units"]
    additive.update(
        variance_increment=str(dv),
        rms_increment=ds,
        single_per_double_upper=float(mass_ratio),
        rms_increment_mass_ratio_upper=float(factor) * ds,
        all_possible_absolute_mass_ratio_upper=float(4 * factor),
    )
    return dict(
        additive=additive,
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                restfile,
                HERE / "FLUCTUATION_BALLS.md",
                HERE / "fluctuation_return.json",
                HERE / "closed_balls.json",
            )
        },
        middle_nonunit_opposite_partner_words=int(nonunit.sum()),
        home_ports=home,
        diameter=int(radii[-1]),
        shots_per_arm=SHOTS,
        equal_distance_centers=[[int(x), int(distance[x, home[0]])] for x in zero],
        arms=result,
        increment=dict(
            zero_by_identity=zero.tolist(),
            middle=float(difference[home[1]]),
            middle_interval=band[home[1]].tolist(),
            minimum=float(difference.min()),
            maximum=float(difference.max()),
            positive_resolved=int(sum(band[:, 0] > 0)),
            negative_resolved=int(sum(band[:, 1] < 0)),
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    dst = HERE / "fluctuation_balls.json"
    if cli.parse_args().write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        check(json.loads(dst.read_bytes()), result)
    print(result["increment"])
