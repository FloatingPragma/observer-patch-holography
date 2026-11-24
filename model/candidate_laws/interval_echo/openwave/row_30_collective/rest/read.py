"""Infer a periodic rest frequency from a saved ring detector count bank."""

import argparse
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks

PAIRS = ((0, 1), (0, 2), (1, 2), (0, 3), (1, 3))


def profile(n):
    f = np.cos(2 * np.pi * np.outer(np.arange(n), np.arange(4)) / n)
    return f / np.linalg.norm(f, axis=0)


def two_lines(trace, times):
    # Locate the two largest nonzero frequencies, then fit both (free phases,
    # amplitudes and offset) to the original trace. The fit has no mass input.
    t = times - times.mean()
    dt = t[1] - t[0]
    coarse = abs(
        np.fft.rfft(np.hanning(len(trace)) * (trace - trace.mean()), n=8 * len(t))
    )
    peaks = find_peaks(coarse, distance=10)[0]
    selected = sorted(peaks[np.argsort(coarse[peaks])[-2:]])
    w = np.asarray(selected) * 2 * np.pi / (8 * len(t) * dt)

    def design(xs):
        phase = t[:, None] * xs
        return np.column_stack([np.cos(phase), np.sin(phase), np.ones(len(t))])

    def residual(xs):
        a = design(xs)
        return a @ np.linalg.lstsq(a, trace, rcond=None)[0] - trace

    best = least_squares(residual, w, ftol=1e-11, xtol=1e-11).x
    return np.sort(best)


def decode(bank, tag, start=0, stop=None):
    n, shots = int(bank["dimers"]), int(bank["shots"])
    f = profile(n)

    def signs(name, reps):
        return 2 * bank[name][start:stop].astype(float) / reps - 1

    phi = np.unwrap(
        np.arctan2(
            signs("sin", int(bank["hand_shots"])), signs("cos", int(bank["hand_shots"]))
        )
    )
    sample = np.arange(len(phi))
    t = np.polyfit(sample, phi, 1)[0] * sample
    readings, planes = {}, {}
    for a, b in PAIRS:
        stem = f"{tag}_{a}{b}_"
        readings[a, b] = two_lines(signs(stem + "pair", shots), t)
        if stem + "scan0" not in bank:
            continue
        demod = np.hanning(len(t)) * np.exp(-1j * t * readings[a, b][1])
        scans = []
        for i, ref in enumerate((f[:, a], f[:, b])):
            data = signs(stem + "scan" + str(i), shots) * (np.abs(ref).sum() - abs(ref))
            obs = demod @ data
            scans.append((obs * np.exp(-1j * np.angle(obs[np.argmax(abs(obs))]))).real)
        planes[a, b] = np.linalg.qr(np.stack(scans, axis=1))[0]
    upper = {pair: ws[1] for pair, ws in readings.items()}
    rest = (upper[0, 1] + upper[0, 2] - upper[1, 2]) / 2
    omega = np.array([rest, *[upper[0, b] - rest for b in range(1, 4)]])
    out = dict(
        clock_step=float(t[1]),
        rest=float(rest),
        omega=omega.tolist(),
        pair_lines={f"{a},{b}": ws.tolist() for (a, b), ws in readings.items()},
        held_sum_error=float(abs(upper[1, 3] - omega[1] - omega[3])),
        max_difference_line_error=float(
            max(
                abs(ws[0] - abs(omega[a] - omega[b])) for (a, b), ws in readings.items()
            )
        ),
    )
    if planes:
        modes = []
        for a in range(3):
            p = sum(q @ q.T for pair, q in planes.items() if a in pair)
            w = np.linalg.eigh(p)[1][:, -1]
            modes.append(
                dict(
                    uniform_weight=float(np.sum(w) ** 2 / n),
                    spatial_gradient_norm=float(np.linalg.norm(w - np.roll(w, 1))),
                )
            )
        out["spatial"] = modes
    return out


if __name__ == "__main__":
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bank", type=Path)
    args = parser.parse_args()
    with np.load(args.bank, allow_pickle=False) as data:
        print(
            json.dumps(
                {tag: decode(data, tag) for tag in ("native", "changed")}, indent=2
            )
        )
