"""Momentum and free wall velocity on the calibrated 29- and 49-dimer strands."""

import argparse
import json
from pathlib import Path

import numpy as np
from model.experiments.count_intervals import clopper_pearson

from .. import experiment, packets

HERE = Path(__file__).resolve().parent
SHOTS = 2**26
MOMENTA = (0.10, 0.20, 0.35, 0.80, np.arccos(1 / 16), 1.9, 2.6)
EXPOSURES = np.array([0.0, 2.0, 4.0, 6.0])
# Fewer than 150,000 binomial marginals across all recorded settings.
TAIL = 0.01 / (2 * 150000)


def confidence(count, n=SHOTS):
    return np.stack(clopper_pearson(count, n, TAIL), axis=-1)


def phase_panel(f, rng):
    coherence = f[:-1].conj() * f[1:]
    total = abs(f[:-1]) ** 2 + abs(f[1:]) ** 2
    settings = []
    for quadrature in (2 * coherence.real, 2 * coherence.imag):
        probs = np.column_stack(
            ((total + quadrature) / 2, (total - quadrature) / 2, 1 - total)
        )
        settings.append(rng.multinomial(SHOTS, np.clip(probs, 0, 1)))
    return np.array(settings)


def angle_interval(center, xb, yb):
    offsets = [
        np.angle(np.exp(1j * (np.arctan2(y, x) - center))) for x in xb for y in yb
    ]
    if xb[0] <= 0 <= xb[1] and yb[0] <= 0 <= yb[1]:
        raise ValueError("Unresolved phase")
    return np.array([center + min(offsets), center + max(offsets)])


def momentum(rows):
    mean = (rows[:, :, 0] - rows[:, :, 1]).sum(axis=1) / (2 * SHOTS)
    p, q = confidence(rows[:, :, 0]), confidence(rows[:, :, 1])
    bands = 0.5 * (p - q[:, :, ::-1]).sum(axis=1)
    angle = np.arctan2(mean[1], mean[0])
    return float(angle), angle_interval(angle, *bands)


def clock(rows):
    hands = 2 * rows / SHOTS - 1
    times = np.unwrap(np.arctan2(hands[1], hands[0]))
    bounds = 2 * confidence(rows) - 1
    return times, np.array(
        [angle_interval(t, *bounds[:, i]) for i, t in enumerate(times)]
    )


def slope(times, tb, density):
    n = density.shape[-1]
    x = np.arange(n) - (n - 1) / 2
    diff = (density[0] - density[1]) @ x / (2 * SHOTS)
    a, b = confidence(density)
    per_site = 0.5 * (a - b[:, :, ::-1])
    moments = np.sum(np.sort(per_site * x[None, :, None], axis=2), axis=1)
    numer, denom = np.zeros(2), np.zeros(2)
    for j in range(1, len(times)):
        for i in range(j):
            dt = tb[j] - tb[i, ::-1]
            dm = moments[j] - moments[i, ::-1]
            assert dt[0] > 0
            prod = np.outer(dt, dm)
            numer += [prod.min(), prod.max()]
            denom += dt**2
    bounds = np.outer(numer, 1 / denom)
    return float(np.polyfit(times, diff, 1)[0]), [
        float(bounds.min()),
        float(bounds.max()),
    ]


def run():
    experiment.native_checks()
    original = json.loads((experiment.HERE / "results.json").read_text())["boards"]
    stock = [b for b in original if b["seed"] == packets.SEED and b["level"] in (1, 2)]
    assert [b["dimers"] for b in stock] == [29, 49]
    rng = np.random.default_rng(11823457)
    meters = json.loads((experiment.HERE / "coupling_meter.json").read_text())[
        "results"
    ]
    raw, summaries = {}, []
    for board in stock:
        n = board["dimers"]
        x = np.arange(n + 1)
        envelope = np.exp(-((x - n / 2) ** 2) / (4 * min(4.0, (n + 1) / 7) ** 2))
        for idx, k in enumerate(MOMENTA):
            pos = packets.synthesize(
                (envelope / np.linalg.norm(envelope)) * np.exp(1j * x * k)
            )
            neg = packets.synthesize(pos.conj())
            base = f"n{n}_p{idx}"
            wave = np.array([phase_panel(f, rng) for f in (pos, neg)])
            raw[base + "_wave"] = wave
            wp, bp = momentum(wave[0])
            wm, bm = momentum(wave[1])
            for j in (experiment.K, experiment.K / 2):
                tag = base + f"_j{int(j)}"
                timing = rng.binomial(
                    SHOTS, (1 + np.array([np.cos(EXPOSURES), np.sin(EXPOSURES)])) / 2
                )
                raw[tag + "_clock"] = timing
                times, tb = clock(timing)
                u, e, v = experiment.dual(n, j)
                prob = np.array(
                    [
                        [packets.occupations(f, u, e, v, t) for t in EXPOSURES]
                        for f in (pos, neg)
                    ]
                )
                assert prob.min() > -1e-12 and prob.max() < 1 + 1e-12
                density = rng.binomial(SHOTS, np.clip(prob, 0, 1))
                raw[tag + "_wall"] = density
                velocity, vb = slope(times, tb, density)
                h = meters["h"]["value"]
                strength = meters["J"]["value"] * j / experiment.K
                gap = 2 * abs(strength - h)
                c2 = 4 * strength * h
                k_read = float((wp - wm) / 2)
                kb = (bp - bm[::-1]) / 2
                grid = 2 * np.pi * np.fft.fftfreq(4096)
                packet_k = (
                    abs(np.fft.fft(envelope * np.exp(1j * x * k_read), 4096)) ** 2
                )
                theory = (
                    c2 * np.sin(grid) / np.sqrt(gap**2 + 2 * c2 * (1 - np.cos(grid)))
                )
                summaries.append(
                    dict(
                        dimers=n,
                        contrast=j,
                        kick_index=idx,
                        read_k=k_read,
                        read_k_band=kb.tolist(),
                        velocity=velocity,
                        velocity_band=vb,
                        envelope_band_velocity=float(
                            np.average(theory, weights=packet_k)
                        ),
                        continuum_speed=c2**0.5,
                        maximum_band_speed=2 * min(strength, h),
                        zero_k_response=c2 / gap,
                        finite_response=velocity / k_read,
                    )
                )
    return summaries, raw


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, help="Write the acquisition into a fresh directory"
    )
    directory = parser.parse_args().output
    result, raw = run()
    if directory:
        directory.mkdir(exist_ok=False)
        (directory / "results.json").write_text(json.dumps(result, indent=2) + "\n")
        np.savez_compressed(directory / "counts.npz", **raw)
    else:
        experiment.verify(result, json.loads((HERE / "results.json").read_text()))
        with np.load(HERE / "counts.npz") as archived:
            assert set(archived.files) == set(raw)
            for key, counts in raw.items():
                np.testing.assert_array_equal(counts, archived[key])
    for r in result:
        print(
            r["dimers"],
            r["contrast"],
            round(r["read_k"], 4),
            round(r["velocity"], 5),
            np.round(r["velocity_band"], 5),
            round(r["envelope_band_velocity"], 5),
        )
