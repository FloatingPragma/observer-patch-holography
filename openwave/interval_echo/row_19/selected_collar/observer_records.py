"""Decode a collar comparison from counted clock, reception and phase-kick records only."""

import argparse
import json
from pathlib import Path

import numpy as np
from mpmath import iv
from scipy.linalg import expm

from matched_einstein_limit.source_implication.beam_response.counts import draw
from matched_einstein_limit.source_implication.directional_stress import phase_read
from matched_einstein_limit.source_implication.material_clocks.measure import (
    blank_clock,
)
from . import measure

HERE = Path(__file__).resolve().parent
N_CAMERA, N_CLOCK, N_FLUX = 2**72, 2**60, 2**62
OFFSET = 1 / 64
Q = -1
BEAM = measure.jacobi.beam
iv.dps = 40


def pixel_image(rows, shots=N_CAMERA, *, stock_count=4):
    # Union bound over all camera settings, launch choices and pixels.
    if len(rows) != 5:
        raise ValueError("Incomplete camera sender set")
    tail = iv.ln(iv.mpf(2 * 3 * 2 * stock_count * 5 * 5) / iv.mpf(".003"))
    pixels = np.empty((5, 5), dtype=object)
    for i, rec in enumerate(rows):
        if len(rec) != 6 or min(rec) < 0 or sum(rec) != shots:
            raise ValueError("Incomplete camera launch record")
        for j, n in enumerate(rec[:5]):
            r = iv.sqrt(2 * iv.mpf(n) * tail) + 2 * tail
            pixels[j, i] = (iv.mpf(n) + iv.mpf([-1, 1]) * r) / shots
    b = np.array(
        [[iv.mpf(int(t)) / BEAM.UNIT for t in row] for row in BEAM.PIXELS], dtype=object
    )
    mid = np.array([[float(x.mid) for x in row] for row in pixels])
    # Center first. Covariance is invariant under either translation.
    left = np.c_[
        np.ones(5, dtype=int), b - (mid.sum(axis=0) @ BEAM.BEARINGS / mid.sum())
    ]
    right = np.c_[
        np.ones(5, dtype=int), b - (mid.sum(axis=1) @ BEAM.BEARINGS / mid.sum())
    ]
    moments = right.T @ pixels @ left
    mass = moments[0, 0]
    return np.array(
        [
            [
                2
                * (moments[j + 1, i + 1] - moments[j + 1, 0] * moments[0, i + 1] / mass)
                / mass
                for i in range(2)
            ]
            for j in range(2)
        ],
        dtype=object,
    )


def focusing(images, clocks):
    dm, dp = clocks[1] - clocks[0], clocks[2] - clocks[1]
    assert dm.a > 0 and dp.a > 0
    left, mid, right = images
    second = np.array(
        [
            [
                2
                * ((right[j, i] - mid[j, i]) / dp - (mid[j, i] - left[j, i]) / dm)
                / (dp + dm)
                for i in range(2)
            ]
            for j in range(2)
        ],
        dtype=object,
    )
    det = mid[0, 0] * mid[1, 1] - mid[1, 0] * mid[0, 1]
    assert 0 not in det
    return (
        -(
            second[0, 0] * mid[1, 1]
            + second[1, 1] * mid[0, 0]
            - second[0, 1] * mid[1, 0]
            - second[1, 0] * mid[0, 1]
        )
        / det
    )


def recorded_clocks(counts, shots=N_CLOCK):
    # p=sin(t/2)^2 on the monotone branch pi < t < 2 pi.
    rad = iv.sqrt(iv.ln(iv.mpf(2 * 3) / iv.mpf(".003")) / (2 * shots))
    clock = []
    for n in counts:
        p = iv.mpf(n) / shots + iv.mpf([-1, 1]) * rad
        clock.append(2 * iv.pi - 2 * iv.atan2(iv.sqrt(p), iv.sqrt(1 - p)))
    return clock


def decode(record):
    flux_shots = record["shots_per_flux_setting"]
    camera_shots, clock_shots = record["shots_per_sender"], record["shots_per_clock"]
    if (camera_shots, clock_shots, flux_shots) not in (
        (N_CAMERA, N_CLOCK, N_FLUX),
        (N_CAMERA, N_CLOCK, 2**78),
        (2**78, 2**68, 2**78),
    ):
        raise ValueError("Unregistered camera, clock or phase-kick budget")
    if len(record["clock_plus_counts"]) != 3 or any(
        n < 0 or n > clock_shots for n in record["clock_plus_counts"]
    ):
        raise ValueError("Invalid clock count")
    clocks = recorded_clocks(record["clock_plus_counts"], clock_shots)
    fs, js = {}, {}
    nsettings = 19 * 4
    magnitude = 5376
    err = iv.mpf(magnitude) * iv.sqrt(
        2 * iv.ln(iv.mpf(2 * nsettings) / iv.mpf(".003")) / flux_shots
    )
    # Fourier differentiation of the nineteen phase-kick means, in their documented order.
    weights = [
        sum(k * iv.sin(2 * iv.pi * m * k / 19) for k in range(1, 10)) / 19
        for m in range(19)
    ]
    for name, data in record["stocks"].items():
        if len(data["flux_plus_counts"]) != 19 or any(
            n < 0 or n > flux_shots for n in data["flux_plus_counts"]
        ):
            raise ValueError("Invalid phase-kick count")
        f = {}
        for tag in ("forward", "reverse"):
            f[tag] = focusing(
                [pixel_image(rows, camera_shots) for rows in data[tag]], clocks
            )
        fs[name] = f["forward"] - f["reverse"]
        js[name] = sum(
            w * (magnitude * (2 * iv.mpf(n) / flux_shots - 1) + iv.mpf([-1, 1]) * err)
            for w, n in zip(weights, data["flux_plus_counts"], strict=True)
        )

    def contrast(values, key):
        return (values[key + "_plus"] - values[key + "_minus"]) / (
            2 * record["half_steps"][key]
        )

    ordinary_f, ordinary_j = contrast(fs, "ordinary"), contrast(js, "ordinary")
    assert 0 not in ordinary_j
    measured = contrast(fs, "collar")
    forecast = ordinary_f / ordinary_j * contrast(js, "collar")
    return dict(
        ordinary_focusing=str(ordinary_f),
        ordinary_current=str(ordinary_j),
        collar_current=str(contrast(js, "collar")),
        predicted_odd_focusing=str(forecast),
        measured_odd_focusing=str(measured),
        measured_minus_predicted=str(measured - forecast),
        distinct=bool((measured - forecast).a > 0 or (measured - forecast).b < 0),
        clocks=list(map(str, clocks)),
    )


def precision_draw(rng, shots, probabilities):
    # Each binomial or multinomial batch fits in a signed 64-bit integer.
    batch = 2**62
    if shots % batch:
        raise ValueError("Incomplete precision batch")
    if np.ndim(probabilities) == 0:
        return sum(map(int, rng.binomial(batch, probabilities, size=shots // batch)))
    rows = []
    for hits in probabilities.T:
        outcomes = np.r_[hits, 1 - hits.sum()]
        if outcomes.min() < 0:
            raise ValueError("Invalid joint reception probabilities")
        counts = rng.multinomial(batch, outcomes, size=shots // batch)
        rows.append(counts.sum(axis=0, dtype=object).tolist())
    return rows


def collect(selected=None, step=None, time=4, seed=190813, *, precision=False):
    ref, prior, *_ = BEAM.aperture.source_and_changes()
    if selected is None:
        selected, _, _, _, source = measure.select_source(ref)
        step = source["step"]
    steps = {"ordinary": 1 / 2048, "collar": step}
    states = {
        k + "_" + tag: ref + sign * steps[k] * d
        for k, d in [("ordinary", prior["ordinary"]), ("collar", selected)]
        for tag, sign in [("minus", -1), ("plus", 1)]
    }
    for x in states.values():
        assert np.min(np.linalg.eigvalsh(x)) > 0
        prep = measure.station.preparation(x)
        assert (
            min(prep["min_kept"], prep["min_dephased"]) > 0
            and prep["max_error"] < 1e-11
        )
    f = BEAM.front.optics.pairing.formation
    board = f.federation(1)
    route = f.region(board)
    rings = BEAM.screens(board, route)
    rng = np.random.default_rng(seed)
    camera_shots, clock_shots, flux_shots = (
        (2**78, 2**68, 2**78) if precision else (N_CAMERA, N_CLOCK, N_FLUX)
    )
    camera_draw = precision_draw if precision else draw
    stocks = {k: {"forward": [], "reverse": []} for k in states}
    clocks = []
    for t in (time - OFFSET, time, time + OFFSET):
        clock = expm(-1j * t * blank_clock(1, 8))
        clock_probability = abs(clock[1, 0]) ** 2
        clocks.append(
            precision_draw(rng, clock_shots, clock_probability)
            if precision
            else int(rng.binomial(clock_shots, clock_probability))
        )
        for tag, screen in [("forward", rings), ("reverse", rings[::-1])]:
            table = BEAM.effects(board, route, screen, Q, t)
            for k, x in states.items():
                stocks[k][tag].append(
                    camera_draw(rng, camera_shots, BEAM.probabilities(table, x))
                )
    _, masks = phase_read.sums_for_charge(board, route, rings, Q, time)
    for key, x in states.items():
        energy = np.einsum("btij,bji->t", np.array(masks["active"]), x).real
        probabilities = (1 + energy / 5376) / 2
        if precision:
            stocks[key]["flux_plus_counts"] = [
                precision_draw(rng, flux_shots, p) for p in probabilities
            ]
        else:
            stocks[key]["flux_plus_counts"] = rng.binomial(
                flux_shots, probabilities
            ).tolist()
    return dict(
        vacancy_charge=Q,
        seed=seed,
        shots_per_sender=camera_shots,
        shots_per_clock=clock_shots,
        shots_per_flux_setting=flux_shots,
        half_steps=steps,
        clock_plus_counts=clocks,
        stocks=stocks,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument(
        "--counts-only",
        action="store_true",
        help="decode archived counts without a simulator state",
    )
    args = parser.parse_args()
    path = HERE / "observer_records.json"
    data = json.loads(path.read_text()) if args.counts_only else collect()
    out = decode(data)
    if args.write:
        data["decoded"] = out
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        assert out == json.loads(path.read_text())["decoded"]
    print(json.dumps(out, indent=2))
