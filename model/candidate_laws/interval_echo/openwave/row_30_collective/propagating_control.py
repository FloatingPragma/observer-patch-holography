"""Broadband wall-pair read on matched preparations with and without a bulk gap."""

import argparse
import json
from pathlib import Path

import numpy as np

from model.executable_model import federation
from .experiment import DT, K, SAMPLES, SHOTS, dual, grow, transport, wedge_norm

HERE = Path(__file__).resolve().parent


def line(n, j, seed):
    # Identical sine launches, local Pauli readout and sampler in both arms.
    m = n + 1
    f, g = np.sqrt(2 / (m + 1)) * np.sin(
        np.pi * np.outer([1, 2], np.arange(1, m + 1)) / (m + 1)
    )
    u, e, _ = dual(n, j)
    fx, gx = u.T @ f, u.T @ g
    c = np.cos(DT * np.outer(np.arange(SAMPLES), e))
    signal = (c @ (fx**2)) * (c @ (gx**2)) - (c @ (fx * gx)) ** 2
    weight = wedge_norm(f, g) ** 2
    rng = np.random.default_rng(seed)
    nplus = rng.binomial(SHOTS, (1 + signal / weight) / 2)
    trace = weight * (2 * nplus / SHOTS - 1)
    window = np.hanning(SAMPLES)
    power = abs(np.fft.rfft(window * trace)) * 2 / window.sum()
    # Record the locations and amplitudes of all resolved positive local peaks.
    # The zero bin contains the preparation's time-independent return.
    (peaks,) = np.where(
        (power[1:-1] > power[:-2]) & (power[1:-1] > power[2:]) & (power[1:-1] > 0.02)
    )
    return [
        [round(float(i * 2 * np.pi / (SAMPLES * DT)), 5), round(float(power[i]), 5)]
        for i in peaks + 1
    ]


def run():
    result = []
    seed = 312143
    for level in range(3):
        strands = grow(federation(level), seed)
        route = max((r for closed, r in strands if not closed), key=len)
        assert len(route) % 2 == 1
        arms = []
        for j in [K, 0.5]:
            read = transport(route, j, seed)
            clicks = max(
                row["distant_max_count_difference"] for row in read["flights"][1:]
            )
            arms.append(
                dict(
                    flight_contrast=j,
                    peaks=line(len(route), j, seed + level),
                    max_remote_difference=clicks,
                    launches_per_arm=read["trials_per_arm"],
                )
            )
        result.append(dict(dimers=len(route), arms=arms))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    args = parser.parse_args()
    observed = run()
    path = HERE / "propagating_control.json"
    if args.record:
        with path.open("x") as f:
            json.dump(observed, f, indent=1)
            f.write("\n")
    elif json.loads(path.read_text()) != observed:
        raise SystemExit(
            "Broadband pair and remote wall controls differ from the receipt"
        )
    else:
        print(
            "Verified traveling critical and massive material on the same three odd strands."
        )
