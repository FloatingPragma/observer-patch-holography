"""Physical spin gaps and a free single-dimer sign measurement on commissioned rings."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, svd

from ...row_30_collective import experiment as strand
from ...row_30_collective.rest.commission import prepare_ring

HERE = Path(__file__).resolve().parent
REST_RECEIPT = HERE.parents[1] / "submission/row_30_cautions/ring_repair.json"
METER = HERE.parents[1] / "row_30_collective/coupling_meter.json"
SHOTS = 65536
TIMES = np.arange(1024) * 0.08


def majorana(n, flip, coupling, transverse):
    shift = np.eye(n, k=-1)
    shift[0, -1] = -flip * (-1) ** n
    return svd(2 * (coupling * np.eye(n) - transverse * shift))


def autocorrelation(n, flip, coupling, transverse, times):
    # Physical X_0 is the dual bond joining B_0 and A_1.
    u, sigma, vt = majorana(n, flip, coupling, transverse)
    cs, sn = np.cos(np.outer(times, sigma)), np.sin(np.outer(times, sigma))
    a, b = u[1], vt[:, 0]
    return (cs @ (a * a)) * (cs @ (b * b)) + (sn @ (a * b)) ** 2


def check_physical_matrix(n, flip):
    # Independent physical orientation basis and exact physical flip projection.
    words = np.arange(2**n)
    z = 1 - 2 * ((words[:, None] >> np.arange(n)) & 1)
    h = np.diag(8.0 * (z * np.roll(z, 1, axis=1)).sum(axis=1))
    for i in range(n):
        h[words, words ^ (1 << i)] = -0.5
    ident = np.eye(len(words))
    orbit = (ident + flip * ident[::-1])[:, : len(words) // 2] / np.sqrt(2)
    energies, vectors = eigh(orbit.T @ h @ orbit)
    sigma = majorana(n, flip, 8, 0.5)[1]
    expected = [
        np.dot(np.asarray(occ) - 0.5, sigma)
        for occ in product((0, 1), repeat=n)
        if sum(occ) % 2 == n % 2
    ]
    np.testing.assert_allclose(energies, sorted(expected), atol=2e-11, rtol=0)
    sector_x = vectors.T @ (orbit.T @ orbit[words ^ 1]) @ vectors
    waits = [0.07, 0.37, 1.24, 2.15]
    direct = [
        np.sum(sector_x**2 * np.cos((energies[:, None] - energies) * t)) / len(energies)
        for t in waits
    ]
    np.testing.assert_allclose(
        direct, autocorrelation(n, flip, 8, 0.5, waits), atol=2e-12, rtol=0
    )


def record(signal, rng):
    # Terminal sign conditional on the first dimer's recorded X sign.
    if np.max(np.abs(signal)) > 1 + 1e-12:
        raise ArithmeticError("invalid sign probability")
    agree = np.round((1 + signal) / 2, 12)
    plus = rng.binomial(SHOTS, 0.5, size=len(agree))
    pp = rng.binomial(plus, agree)
    mm = rng.binomial(SHOTS - plus, agree)
    return np.column_stack([pp, plus - pp, SHOTS - plus - mm, mm])


def peaks(words, clocks):
    cos, sin = 2 * np.asarray(clocks) / SHOTS - 1
    phase = np.unwrap(np.arctan2(sin, cos))
    dt = np.polyfit(np.arange(len(phase)), phase, 1)[0]
    sign = (words[:, 0] + words[:, 3] - words[:, 1] - words[:, 2]) / SHOTS
    win = np.hanning(len(sign))
    power = 2 * np.abs(np.fft.rfft(win * (sign - sign.mean()), n=16384)) / win.sum()
    freqs = 2 * np.pi * np.fft.rfftfreq(16384, dt)

    def peak(low, high):
        band = np.flatnonzero((freqs >= low) & (freqs <= high))
        idx = band[power[band].argmax()]
        return [round(float(x), 8) for x in (freqs[idx], power[idx])]

    return dict(dt=round(float(dt), 10), rest=peak(14.9, 15.1), pair=peak(29.8, 34.2))


def run():
    strand.native_checks()
    for n in (5, 6, 7, 8):
        for flip in (-1, 1):
            check_physical_matrix(n, flip)
    lower, upper = json.loads(METER.read_text())["results"]["predicted_m"][
        "bootstrap95"
    ]
    if not 14.99 < lower < upper < 15.01:
        raise ValueError("independent meter interval left the enrolled exclusion")
    result = []
    for previous in json.loads(REST_RECEIPT.read_text())["rings"]:
        if not lower <= previous["readings"]["native"]["rest"] <= upper:
            raise ValueError("the two rest instruments disagree")
        seed = previous["seed"]
        ring, receipt = prepare_ring(1, seed)
        n = len(ring)
        if (
            receipt["circuit_sha256"] != previous["circuit_sha256"]
            or n != previous["dimers"]
        ):
            raise ValueError("different commissioned ring")
        flip = -((-1) ** n)
        radius = n * np.sqrt(513) / 31
        if not (n * n * 513 * 10000 < (1499 * 31) ** 2 and 32 - radius > 15.01):
            raise ValueError("rest interval intersects a possible beat")
        radius_grid = max(np.sum(abs(majorana(n, f, 8, 0.5)[1] - 16)) for f in (-1, 1))
        rng = np.random.default_rng([3403, seed])
        clock = [
            rng.binomial(SHOTS, np.round((1 + func(TIMES)) / 2, 12)).tolist()
            for func in (np.cos, np.sin)
        ]
        arms = {}
        for name, j, h in [
            ("native", 8, 0.5),
            ("exchange_only", 0, 0.5),
            ("planted", 3.75, 0),
        ]:
            counts = record(autocorrelation(n, flip, j, h, TIMES), rng)
            arms[name] = dict(counts=counts.tolist(), peaks=peaks(counts, clock))
        if (
            arms["native"]["peaks"]["rest"][1] >= 0.01
            or arms["native"]["peaks"]["pair"][1] <= 0.03
            or arms["planted"]["peaks"]["rest"][1] <= 0.35
            or np.asarray(arms["exchange_only"]["counts"])[:, 1:3].any()
        ):
            raise ValueError(
                f"failed detector gate: {[(a, d['peaks']) for a, d in arms.items()]}"
            )
        result.append(
            dict(
                dimers=n,
                seed=seed,
                circuit=receipt["circuit_sha256"],
                flip=flip,
                strict_bound=round(float(radius), 10),
                grid_sum=round(float(radius_grid), 10),
                shots=SHOTS,
                clock_counts=clock,
                arms=arms,
            )
        )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    output = parser.parse_args().output
    data = run()
    if output:
        with output.open("x") as f:
            json.dump(data, f, separators=(",", ":"))
            f.write("\n")
    elif data != json.loads((HERE / "records.json").read_text()):
        raise SystemExit("retained records disagree")
    for r in data:
        print(
            r["dimers"],
            r["strict_bound"],
            r["grid_sum"],
            {a: v["peaks"] for a, v in r["arms"].items()},
        )
