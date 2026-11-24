"""Measure propagation separately from the static boundary spectrum."""

import argparse
import json

import numpy as np

from . import experiment as ex, schedule

M = 46
CENTRES = np.arange(0, M - 1 + 0.01, 1 / 8)
HOLDS = (0, 32, 64, 96, 128)
SHOTS = 2**22


def detectors(centre, direction):
    x = np.arange(M) - centre
    g = np.exp(-x * x / (4 * 3**2))
    v = g / np.linalg.norm(g)
    w = (x - np.sum(x * v * v)) * g
    w /= np.linalg.norm(w)
    direction = np.asarray(direction) / np.linalg.norm(direction)
    return np.column_stack((np.kron(v, direction), np.kron(w, direction)))


def run():
    j, h = ex.couplings(8)
    f = h / j
    dt = schedule.DELTA * (1 + f)
    adv = np.eye(2 * M)
    schedule.rotate_pairs(adv, np.arange(1, 2 * M - 1, 2), 2 * h * schedule.DELTA)
    schedule.rotate_pairs(adv, np.arange(0, 2 * M, 2), 2 * j * f * schedule.DELTA)
    oracle = {t: np.linalg.matrix_power(adv, t) for t in HOLDS}
    rng = np.random.default_rng(201128)
    scans = []
    bank = {}
    for direction in ((1, 1), (1, -1), (1, 0), (0, 1)):
        all_det = np.stack([detectors(x, direction).T for x in CENTRES])
        for start in (14, 31):
            init = detectors(start, direction)
            track = []
            for t in HOLDS:
                overlaps = all_det @ (oracle[t] @ init)
                contrast = np.linalg.det(overlaps)
                draws = rng.binomial(
                    SHOTS, np.round(np.clip((1 + contrast) / 2, 0, 1), 13)
                )
                bank[f"{direction[0]}_{direction[1]}_{start}_{t}"] = draws
                obs = 2 * draws / SHOTS - 1
                peak = int(np.argmax(obs))
                epsilon = np.sqrt(2 * np.log(2 * 40 * len(CENTRES) / 0.01) / SHOTS)
                contenders = CENTRES[obs >= obs[peak] - 2 * epsilon]
                track.append(
                    dict(
                        t=t,
                        centre=float(CENTRES[peak]),
                        contrast=float(obs[peak]),
                        peak_grid_bounds=[
                            float(min(contenders)),
                            float(max(contenders)),
                        ],
                    )
                )
            first, last = track[1], track[3]
            shift = last["centre"] - first["centre"]
            qualifies = (
                min(
                    first["centre"],
                    M - 1 - first["centre"],
                    last["centre"],
                    M - 1 - last["centre"],
                )
                >= 6
                and min(first["contrast"], last["contrast"]) > 0.75
                and abs(shift) >= 4
            )
            scans.append(
                dict(
                    orientation=direction,
                    start=start,
                    track=track,
                    qualified=qualifies,
                    sites_per_exchange_time=shift / (64 * dt),
                )
            )
    qualified = [s for s in scans if s["qualified"]]
    assert (
        len(qualified) == 2 and qualified[0]["start"] + qualified[1]["start"] == M - 1
    )
    bounds = []
    for s in qualified:
        x = s["track"][1]["peak_grid_bounds"]
        y = s["track"][3]["peak_grid_bounds"]
        ends = sorted(abs(z) / (64 * dt) for z in (y[0] - x[1], y[1] - x[0]))
        bounds.append(ends)
    speed = np.mean([abs(s["sites_per_exchange_time"]) for s in qualified])
    coefficient = np.pi * speed
    band = np.pi * np.mean(bounds, axis=0)
    return dict(
        shots_per_read=SHOTS,
        period=dt,
        centres=CENTRES.tolist(),
        scans=scans,
        velocity=speed,
        coefficient=coefficient,
        coefficient_grid_counting_bounds=band.tolist(),
    ), bank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out, bank = run()
    data = json.loads(json.dumps(ex.portable(out)))
    report, tape = ex.HERE / "flight.json", ex.HERE / "flight.npz"
    if args.write:
        report.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
        np.savez_compressed(tape, **bank)
    else:
        if data != json.loads(report.read_text()):
            raise ValueError("Flight summary differs")
        with np.load(tape) as stored:
            assert set(stored.files) == set(bank)
            for k, a in bank.items():
                np.testing.assert_array_equal(stored[k], a)
    print("Verified 40 spatial packet scans on the independently timed critical guide.")


if __name__ == "__main__":
    main()
