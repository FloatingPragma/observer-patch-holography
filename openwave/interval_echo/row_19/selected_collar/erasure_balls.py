"""Compare the sign forced by data processing with fixed-proper-volume ball areas."""

import argparse
import json
from hashlib import sha256
from itertools import combinations

import numpy as np

from . import angular_erasure
from . import transmission_balance as measure

reg = measure.reg
iv = reg.iv
geom = reg.six.survey.geometry


def solid_angle_jets(square, site, others):
    """Van Oosterom–Strackee formula, differentiated on its positive branch."""
    a, b, c = (square[site, i] for i in others)
    g, h, k = (
        (square[site, i] + square[site, j] - square[i, j]) * 0.5
        for i, j in combinations(others, 2)
    )
    root = (a * b * c + 2 * g * h * k - a * k * k - b * h * h - c * g * g).sqrt()
    denominator = (
        a.sqrt() * b.sqrt() * c.sqrt() + g * c.sqrt() + h * b.sqrt() + k * a.sqrt()
    )
    assert denominator.v.a > 0
    ratio = root * denominator.unary(
        1 / denominator.v, -1 / denominator.v**2, 2 / denominator.v**3
    )
    return ratio.unary(
        2 * iv.atan2(ratio.v, 1),
        2 / (1 + ratio.v**2),
        -4 * ratio.v / (1 + ratio.v**2) ** 2,
    )


def run():
    iv.dps = 35
    paths = [
        measure.HERE / name
        for name in (
            "returned_head_lab.json",
            "product_pulse.json",
            "FIXED_VOLUME_ERASURE.md",
            "ANGULAR_ERASURE.md",
        )
    ]
    src, control = [json.loads(p.read_bytes()) for p in paths[:2]]
    tapes = {**src["tapes"], "product": control["tape"]}
    n = src["shots_per_setting"]
    incoming, local, orders, _, variations = measure.instrument(
        {k: t["incoming"] for k, t in tapes.items()}, n, False
    )
    # Half the spread of the input log-order bounds controls every channel's
    # transmitted score on a trace-zero ordinary redistribution.
    D, e = variations[0]
    norm = sum(np.linalg.norm(x, ord="nuc") for x in D) + 4 * sum(e) + 1e-14
    max_price = (
        norm * (max(v[1] for v in orders) - min(v[0] for v in orders)) / 2 + 1e-14
    )
    ordinary_any = reg.band(*(incoming[0] + local[0])) + reg.band(-max_price, max_price)
    # Marginals and parity agree exactly along the two joint directions.
    assert np.all((local[1:4, 0] <= 0) & (local[1:4, 1] >= 0))
    local[1:4] = 0
    outputs = np.array(
        [
            measure.instrument(
                {k: t["returned"][i] for k, t in tapes.items()}, n, True, orders
            )[0]
            for i in range(30)
        ]
    )
    raw = {k: r["time"] for k, r in tapes.items()}
    means, slopes = reg.timing_slopes(raw, n, measure.CONTRASTS)
    square = reg.six.squared_range_jets(means["reference"], slopes)

    def endpoints(key):
        arr = np.zeros((6, 6), dtype=object)
        for route, t in zip(reg.six.survey.ROUTES, means[key], strict=True):
            arr[route] = t
        return geom.observed_squared(arr)

    cube = endpoints("product")
    angles = {}
    finite = {}
    min_altitude = np.inf
    baseline = np.vectorize(
        lambda x: x.v if hasattr(x, "v") else iv.mpf(0), otypes=[object]
    )(square)
    for sites in combinations(range(6), 4):
        for i in sites:
            rest = tuple(j for j in sites if j != i)
            angles[sites, i] = solid_angle_jets(square, i, rest)
            for s in (baseline, cube):
                omega, vol = geom.solid_angle(s, i, rest)
                altitude = 3 * vol / geom.face_area(s, rest)
                min_altitude = min(min_altitude, float(altitude.a))
            finite[sites, i] = omega
    assert min_altitude > 1
    rows = []
    for inside in combinations(range(6), 3):
        if 0 not in inside or ((4 in inside) == (5 in inside)):
            continue
        cells = [s for s in combinations(range(6), 4) if len(set(s) & set(inside)) == 2]
        assert len(cells) == 9
        for i in range(6):
            omega = sum(angles[c, i] for c in cells if i in c)
            end = sum(finite[c, i] for c in cells if i in c)
            radius = iv.mpf("0.5")
            final_radius = radius * (omega.v / end) ** (iv.mpf(1) / 3)
            assert final_radius.b < min_altitude
            # V=omega*r^3/3. At fixed V, dA=(r^2/3)d(omega).
            dA = omega.g * (radius**2 / 3)
            change = omega.v * radius**2 * ((end / omega.v) ** (iv.mpf(1) / 3) - 1)
            selected = [j for j, r in enumerate(reg.six.survey.ROUTES) if r[0] == i]
            out = outputs[selected].mean(0)
            lost = incoming - out[:, ::-1]
            meters = [reg.band(*z) for z in (local + lost)]
            assert meters[3].a > 0
            gain = None if dA[0].a <= 0 <= dA[0].b else -meters[0] / dA[0]
            row = {
                "triangle": list(inside),
                "center_index": i,
                "center_port": src["ports"][i],
                "equilibrium_solid_angle": reg.six.bounds(omega.v),
                "proper_volume": reg.six.bounds(omega.v * radius**3 / 3),
                "final_equal_volume_radius": reg.six.bounds(final_radius),
                "finite_erasure_area_change": reg.six.bounds(change),
                "area_slopes": list(map(reg.six.bounds, dA)),
                "observer_total_scores": list(map(reg.six.bounds, meters[:4])),
                "absolute_lost_relative_entropy": (-lost[4, ::-1]).tolist(),
                "ordinary_area_gain": reg.six.bounds(gain)
                if gain is not None
                else None,
            }
            if gain is not None:
                row["residuals_at_ordinary_gain"] = [
                    reg.six.bounds(meters[j] + gain * dA[j]) for j in (1, 2, 3)
                ]
                gain_any = -ordinary_any / dA[0]
                row["any_receiver_ordinary_gain"] = reg.six.bounds(gain_any)
                # On erasure, q lies between zero and the INPUT Jeffreys
                # divergence. Any additional source invisible in the ordinary
                # calibration must close exactly this residual.
                row["additional_score_needed_any_receiver"] = reg.six.bounds(
                    -gain_any * dA[3] - reg.band(0, incoming[3, 1])
                )
                row["additional_score_needed_measured_receiver"] = reg.six.bounds(
                    -gain * dA[3] - meters[3]
                )
            row["moving_aperture"] = angular_erasure.compare(
                i,
                cells,
                angles,
                omega,
                dA,
                incoming,
                local,
                outputs,
                ordinary_any,
                orders,
            )
            rows.append(row)
    return {
        "inputs": {p.name: sha256(p.read_bytes()).hexdigest() for p in paths},
        "minimum_opposite_face_altitude_at_either_endpoint": min_altitude,
        "input_symmetric_divergence": incoming[3].tolist(),
        "balls": rows,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    data = run()
    dest = measure.HERE / "erasure_balls.json"
    if cli.parse_args().write:
        dest.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(dest.read_bytes()), data)
    print(
        len(data["balls"]), sum(row["area_slopes"][3][0] > 0 for row in data["balls"])
    )
