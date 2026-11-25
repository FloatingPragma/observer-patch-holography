"""Finite two-port calibration with independent terminal material and counter records."""

import argparse
import json

import numpy as np
from scipy.linalg import expm

from . import measure as cur

SHOTS = 2**16
SEED = 122015


def run():
    rng = np.random.default_rng(SEED)
    exchange, potential = cur.native_law.local_term_parts()
    h = (exchange + cur.KAPPA * potential)[:9, :9]
    words = cur.native_law.SRC
    records = []
    for time in [0.0, np.pi / 2, np.pi]:
        evolve = expm(-1j * time * h)
        for cut in range(10):
            for a, b in cur.belt.ring(10).seams:
                tagged = abs(cur.crossing(a, b, cut))
                for src in [(1, 0), (0, 1), (0, 0)]:
                    p = np.zeros((len(words), 2))
                    initial = words.index(src)
                    for j, word in enumerate(words):
                        flip = tagged * ((word[0] - src[0]) % 2)
                        p[j, flip] = abs(evolve[j, initial]) ** 2
                    counts = rng.multinomial(SHOTS, p.ravel() / p.sum()).reshape(-1, 2)
                    records.append(
                        dict(
                            cut=cut,
                            edge=[a, b],
                            source=src,
                            wait_pi_units=time / np.pi,
                            occupied_bins=[
                                [int(j), int(bit), int(counts[j, bit])]
                                for j, bit in zip(*np.nonzero(counts), strict=True)
                            ],
                        )
                    )
    return dict(shots=SHOTS, seed=SEED, local_words=words, banks=records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = json.dumps(run(), separators=(",", ":"), allow_nan=False) + "\n"
    dest = cur.HERE / "calibration.json"
    if args.write:
        dest.write_text(receipt)
    elif dest.read_text() != receipt:
        raise ValueError("Contact calibration replay changed")
    print("1800 full two-port calibration banks replayed")


if __name__ == "__main__":
    main()
