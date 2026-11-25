"""Frequency and position records of charged one-, three- and five-defect cores."""

import argparse
import json
from math import pi
from pathlib import Path

import numpy as np

from model.experiments.interval_echo_gauge.flag_measure.measure import meter
from .sector import columns, core_word, eigenpair, line, orbits, prepare

HERE = Path(__file__).resolve().parent
N = 15
SHOTS = 2**26
FILTER_STEPS = 32768
POWERS = tuple(range(-7, 26))
# Delays are a fixed power-of-two ladder for every size and both laws.
FLIGHTS = (0, 4, 256, 1024, 4096, 65536, 1048576, 4194304, 16777216, 33554432)


def fit(banks, center):
    previous = 0.0
    for bank in banks:
        quadratures = []
        errors = []
        for bright, dark, *_ in bank["counts"]:
            accepted = bright + dark
            quadratures.append((bright - dark) / accepted)
            # Two bounded signs; union over all six sources, three modes, powers and axes.
            errors.append(
                np.sqrt(2 * np.log(2 * 6 * 3 * len(POWERS) * 2 / 0.01) / accepted)
            )
        z = complex(*quadratures)
        angle = np.angle(z)
        time = bank["time"]
        turns = np.rint((previous * time - angle) / (2 * pi))
        previous = (angle + 2 * pi * turns) / time
        # 4e-12 native frequency allowance for the numerical eigensolver/phases.
        radius = np.arcsin(np.hypot(*errors) / abs(z)) / time + 4e-12
        assert radius * time < 0.02
    return dict(energy=float(center + previous), radius_99=float(radius))


def filter_amplitude(offset, center, d):
    if d == 1:
        return 1.0
    x = offset / center
    return np.exp(-1j * FILTER_STEPS * x + FILTER_STEPS * np.log(np.cos(x)))


def clock_banks(rng, m, d, energy):
    center = energy["center"]
    delta = energy["offset"]
    w = energy["core_weight"]
    success = w * abs(filter_amplitude(delta, center, d)) ** 2
    rows = []
    for power in POWERS:
        t = float(2.0**power)
        theta = delta * t
        pairs = []
        for quadrature in (np.cos(theta), np.sin(theta)):
            # Choose the Fourier output before launch. Count every other output and every filter refusal.
            tagged = int(rng.binomial(SHOTS, 1 / N))
            admitted = int(rng.binomial(tagged, success))
            bright = int(rng.binomial(admitted, np.clip((1 + quadrature) / 2, 0, 1)))
            pairs.append([bright, admitted - bright, tagged - admitted, SHOTS - tagged])
        rows.append(dict(time=t, counts=pairs))
    return dict(mode=m, records=rows, fit=fit(rows, center))


def flight(rng, top, d, ring):
    lines = [top[min(i, N - i)] for i in range(N)]
    energies = np.array([row["offset"] for row in lines])
    w = np.array([row["core_weight"] for row in lines])
    filt = np.array([filter_amplitude(x, lines[0]["center"], d) for x in energies])
    success = np.mean(w * abs(filt) ** 2)
    source_owners = {p // 12 for p in ring[:d]}
    entirely_remote = [
        j
        for j in range(N)
        if all(ring[(j + i) % N] // 12 not in source_owners for i in range(d))
    ]
    rows = []
    for time in FLIGHTS:
        amp = np.fft.ifft(w * filt * np.exp(-1j * energies * time))
        probabilities = np.abs(amp) ** 2
        bins = np.r_[
            probabilities, max(0.0, success - sum(probabilities)), max(0.0, 1 - success)
        ]
        counts = rng.multinomial(SHOTS, bins / sum(bins))
        rows.append(
            dict(
                time=time,
                core_by_first_port=counts[:N].tolist(),
                other_loads=int(counts[-2]),
                source_refused=int(counts[-1]),
                all_outside_initial_owners=int(sum(counts[j] for j in entirely_remote)),
            )
        )
    return dict(
        original_owners=sorted(source_owners),
        disjoint_owner_bins=entirely_remote,
        records=rows,
    )


def local_recoil(rng, d, kappa, top):
    """An ordinary single-port vacancy read disturbs a rested core; resolve its outgoing translation."""
    center = top[0]["center"]
    _, v0, _, _ = eigenpair(N, d, kappa, 0)
    admitted = (
        top[0]["core_weight"]
        * abs(filter_amplitude(top[0]["offset"], center, d)) ** 2
        / N
    )
    words = np.array(orbits(N, d)[0])
    physical_tags = np.asarray(meter())[words]
    rows = []
    for m in (1, 2):
        _, vm, _, _ = eigenpair(N, d, kappa, m)
        shifts = np.arange(N)
        # (T^s w)[0] = w[-s].  The complete vacancy flag preserves all matching configurations.
        mask = (
            (physical_tags[:, -shifts % N] == 0) @ np.exp(-2j * pi * m * shifts / N) / N
        )
        projected = mask * v0
        mode = np.vdot(projected, projected).real
        hit = (d + 1) / (2 * N)
        band = (
            abs(np.vdot(vm, projected) * filter_amplitude(top[m]["offset"], center, d))
            ** 2
        )
        probabilities = np.array(
            [
                1 - admitted,
                admitted * (1 - hit),
                admitted * (hit - mode),
                admitted * (mode - band),
                admitted * band,
            ]
        )
        counts = rng.multinomial(
            SHOTS, np.maximum(probabilities, 0) / np.maximum(probabilities, 0).sum()
        )
        rows.append(
            dict(
                mode=m,
                counts=counts.tolist(),
                matched_no_probe_outgoing_mode=0,
                p_outgoing_bound_given_source=float(band),
                p_bound_given_hit_and_mode=float(band / mode),
            )
        )
    return rows


def source(rng, d, kappa, ring):
    top = [line(N, d, kappa, m) for m in range(N // 2 + 1)]
    center = top[0]["center"]
    bound = columns(N, d, kappa)[-1]
    w = core_word(N, d)
    cost = kappa * sum((a - b) ** 2 for a, b in zip(w, w[1:] + w[:1]))
    unequal = sum(a != b for a, b in zip(w, w[1:] + w[:1]))
    parity = int(rng.binomial(SHOTS, unequal / (2 * N)))
    assert cost + unequal / 2 == center
    if d > 1:
        assert bound == (4 * d - 4) * kappa + d + 2
    norm_tail = (
        float(np.cos((center - bound) / center) ** FILTER_STEPS) if d > 1 else 0.0
    )
    assert norm_tail < 1e-16
    banks = [clock_banks(rng, m, d, top[m]) for m in range(3)]
    energies = [p["fit"]["energy"] for p in banks]
    return dict(
        kappa=kappa,
        defects=d,
        orbit_dimension=len(orbits(N, d)[0]),
        sharp_random_bond_minus=parity,
        sharp_contrast=cost,
        sharp_energy_read=float(cost + N * parity / SHOTS),
        highest_separated_ceiling=bound if d > 1 else None,
        sharp_core_probability_floor=(center - bound) ** 2
        / ((center - bound) ** 2 + d + 1)
        if d > 1
        else 1.0,
        ideal_lines=top,
        filter_discarded_amplitude_bound=norm_tail,
        momentum_metrology=banks,
        recoil_after_rest=local_recoil(rng, d, kappa, top),
        measured_gaps=[e - energies[0] for e in energies[1:]],
        localized_flight=flight(rng, top, d, ring),
    )


def measure():
    preparations = [prepare(d) for d in (1, 3, 5)]
    ring = preparations[0]["route"]
    rng = np.random.default_rng(2301307)
    return dict(
        shots_per_setting=SHOTS,
        filter_steps=FILTER_STEPS,
        preparations=preparations,
        cases=[source(rng, d, k, ring) for k in (8, 4) for d in (1, 3, 5)],
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path)
    dest = p.parse_args().output
    record = measure()
    if dest:
        with dest.open("x") as f:
            json.dump(record, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        for row in record["cases"]:
            print(
                row["kappa"],
                row["defects"],
                row["momentum_metrology"][0]["fit"],
                row["measured_gaps"],
            )
