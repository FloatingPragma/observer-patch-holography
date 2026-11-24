"""Decode the representative observer bank using only NumPy and saved counts."""

import argparse
import json
from pathlib import Path

import numpy as np


def decode(path):
    with np.load(path) as bank:
        n, shots, clock_shots = (int(bank[k]) for k in ("n", "shots", "clock_shots"))

        def sign(k, exposure):
            return 2.0 * bank[k] / exposure - 1.0

        theta = np.unwrap(
            np.arctan2(sign("clock_sin", clock_shots), sign("clock_cos", clock_shots))
        )
        index = np.arange(len(theta))
        elapsed = np.polyfit(index, theta, 1)[0] * index
        window = np.hanning(len(theta))
        m = n + 1
        f = np.sqrt(2 / (m + 1)) * np.sin(
            np.pi * np.outer(np.arange(1, m + 1), np.arange(1, 4)) / (m + 1)
        )
        pair_sums, spaces = {}, {}
        for a, b in [(0, 1), (0, 2), (1, 2)]:
            tag = f"p{a + 1}_{b + 1}_"
            source = np.sum(
                abs(np.triu(np.outer(f[:, a], f[:, b]) - np.outer(f[:, b], f[:, a]), 1))
            )
            y = sign(tag + "sum", shots) * source**2
            spec = abs(np.fft.rfft(window * y))
            start = len(spec) // 2
            pos = int(start + np.argmax(spec[start:-1]))
            lo, hi, right = np.log(spec[pos - 1 : pos + 2])
            offset = np.clip((lo - right) / (2 * (lo - 2 * hi + right)), -0.5, 0.5)
            w = 2 * np.pi * (pos + offset) / (len(y) * elapsed[1])
            pair_sums[a, b] = w
            profiles = []
            for z, ref in enumerate((f[:, a], f[:, b])):
                response = (
                    source
                    * (sum(abs(ref)) - abs(ref))
                    * sign(tag + "scan" + str(z), shots)
                )
                profile = (window * np.exp(-1j * w * elapsed)) @ response
                profiles.append(
                    np.real(
                        profile
                        * np.exp(-1j * np.angle(profile[np.argmax(abs(profile))]))
                    )
                )
            spaces[a, b] = np.linalg.qr(np.asarray(profiles).T)[0]
        e0 = (pair_sums[0, 1] + pair_sums[0, 2] - pair_sums[1, 2]) / 2
        result = []
        for mode in range(3):
            projs = [q @ q.T for pair, q in spaces.items() if mode in pair]
            w = np.linalg.eigh(sum(projs))[1][:, -1]
            mid = w[1:-1]
            cos = np.dot(mid, w[:-2] + w[2:]) / (2 * np.dot(mid, mid))
            energy = e0 if mode == 0 else pair_sums[0, mode] - e0
            result.append(
                dict(
                    mode=mode + 1,
                    k=float(np.arccos(np.clip(cos, -1, 1))),
                    omega=float(energy),
                )
            )
        return dict(clock_step=elapsed[1], representative_modes=result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "path",
        type=Path,
        nargs="?",
        default=Path(__file__).with_name("example_counts.npz"),
    )
    print(json.dumps(decode(parser.parse_args().path), indent=2))
