"""Read the modular energy of a station's actual hit/miss algebra on ordinary and overlap flights."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import distributed_source as dist

HERE, regions, iv, lab, reg = dist.HERE, dist.regions, dist.iv, dist.lab, dist.reg


def local_prices(readout):
    block, err = readout
    miss = 1 - np.trace(block, axis1=1, axis2=2).real.sum()
    dm = 4 * sum(err) + 1e-13
    e = 2 * max(err) + 1e-13
    floor = np.linalg.eigvalsh(block).min() - e
    assert floor > 0 and miss > dm
    K = -regions.c.station.log_blocks(block) + np.log(miss) * np.eye(4)
    energy = np.stack([np.zeros((4, 4)), np.eye(4)])
    base_odd = np.trace(block[1]).real
    odd = reg.band(base_odd - err[1] - 1e-13, base_odd + err[1] + 1e-13)
    return {
        "modular": (K, e / floor + dm / (miss - dm) + 1e-12),
        "parked": (energy, 0.0),
    }, odd


def same_sign(v):
    return bool(v.a > 0 or v.b < 0)


def run():
    iv.dps = 35
    source, clocks, ctrl = (
        dist.read(n)
        for n in (
            "distributed_source",
            "six_meter_clocks",
            "ordinary_service_distribution",
        )
    )
    ports = source["ports"]

    def band(x):
        return lab.six.bounds(x)

    # The held dense reference and all ordinary standards represent the same uniform dense background.
    priors = [
        local_prices(regions.at_port(source["target_counts"]["dense"], i))
        for i in range(6)
    ]
    dose = regions.density.DOSE
    means, _ = reg.timing_slopes(
        clocks["clock_counts"], clocks["starts"], (("held", "dense", 1),)
    )
    base = np.array(
        [
            (1 - dose) * u + dose * v
            for u, v in zip(means["reference"], means["dense"], strict=True)
        ],
        dtype=object,
    )
    routes = lab.six.survey.ROUTES
    raw = np.full((6, 6), iv.mpf(0), dtype=object)
    for route, m in zip(routes, base, strict=True):
        raw[route] = m
    square = regions.balls.geom.observed_squared(raw)
    gradients = regions.ball_area_rates(base, np.eye(30), tangent=True)
    timing = {}
    local = {}

    def score(hi, lo):
        return [
            regions.local_contrasts(
                priors[i][0], regions.at_port(hi, i), regions.at_port(lo, i)
            )
            for i in range(6)
        ]

    for key, arms in source["standard_counts"].items():
        t, h = key.split("/", 1)
        col = ctrl["windows"][t]["ordinary_columns"][h]
        timing[key] = [iv.mpf(col[f"clock/{a}/{b}"]) for a, b in routes]
        local[key] = score(arms["ordinary"], arms["reference"])
    histories = list(timing)
    timing["ordinary"] = [
        sum(z) / len(histories) for z in zip(*timing.values(), strict=True)
    ]
    local["ordinary"] = [
        {
            k: sum(local[h][i][k] for h in histories) / len(histories)
            for k in ("modular", "parked")
        }
        for i in range(6)
    ]
    timing["overlap"] = [
        dose * (a - b) for a, b in zip(means["held"], means["dense"], strict=True)
    ]
    local["overlap"] = score(
        source["target_counts"]["overlap"], source["target_counts"]["dense"]
    )
    rows = []
    for (tri, i), grad in gradients.items():
        cells = [t for t in regions.c.TETRA if i in t and len(set(t) & set(tri)) == 2]
        omega = sum(
            regions.balls.geom.solid_angle(square, i, tuple(j for j in t if j != i))[0]
            for t in cells
        )
        area0 = omega / 4
        E0 = dose * priors[i][1]
        readings = {}
        for k, direction in timing.items():
            area = sum(g * d for g, d in zip(grad, direction, strict=True))
            entropy, energy = local[k][i]["modular"], local[k][i]["parked"]
            half = iv.pi / 2
            readings[k] = {
                "area": area,
                "modular": entropy,
                "energy": energy,
                "fixed_radius_remainder": entropy - half * energy,
                "fixed_volume_remainder": entropy - half * (energy - E0 * area / area0),
            }
        allowed = [
            regions.admissible(band(readings[h]["modular"]), band(readings[h]["area"]))
            for h in histories
        ]
        gain = regions.overlap(allowed)
        o, y = readings["ordinary"], readings["overlap"]
        coefficient = -o["modular"] / o["area"] if same_sign(o["area"]) else None
        forecast = (
            -y["modular"] / coefficient
            if coefficient is not None and coefficient.a > 0
            else None
        )
        det = y["modular"] * o["area"] - o["modular"] * y["area"]
        rows.append(
            {
                "triangle": list(tri),
                "vertex": i,
                "port": ports[i],
                "base_area": band(area0),
                "baseline_parked_energy": band(E0),
                "readings": {
                    k: {n: band(z) for n, z in r.items()} for k, r in readings.items()
                },
                "individual_ordinary_gains": dict(zip(histories, allowed, strict=True)),
                "joint_ordinary_gain": gain,
                "uniform_ordinary_gain": None
                if coefficient is None
                else band(coefficient),
                "ordinary_normalized_overlap_area": None
                if forecast is None
                else band(forecast),
                "uniform_ordinary_overlap_determinant": band(det),
            }
        )
    return {
        "balls": rows,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in (
                "LOCAL_CAP_PROTOCOL.md",
                "distributed_source.json",
                "six_meter_clocks.json",
                "ordinary_service_distribution.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    path = HERE / "local_cap.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    rows = r["balls"]
    print(
        "ordinary compatible",
        sum(bool(z["joint_ordinary_gain"]) for z in rows),
        "positive uniform calibration",
        sum(z["ordinary_normalized_overlap_area"] is not None for z in rows),
        "held determinant excludes",
        sum(
            a * b > 0
            for a, b in (z["uniform_ordinary_overlap_determinant"] for z in rows)
        ),
    )
