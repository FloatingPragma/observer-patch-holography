"""Read the registered radial statistic with finite independent number measurements."""

import argparse
from hashlib import sha256
import json
from math import log, sqrt

import numpy as np

from . import measure as m

NAMES = (
    "target",
    "phase_feedback",
    "exchange_1_1",
    "balanced_sharp",
    "balanced_coherent",
    "balanced_mass",
)
SHOTS = 2**26
SEED = 230419


def acquire(data):
    rng = np.random.default_rng(SEED)
    _, _, dist, origins = m.fp.make_board()
    family = len(NAMES) * len(origins) * len(m.fp.TIMES)
    half_width = sqrt(4 * sum(dist**2) * log(2 * family / 0.01) / (2 * SHOTS))
    records = {}
    for name in NAMES:
        reference = data["standards"][name]
        rows = []
        for port, probe in zip(origins, reference["probes"], strict=True):
            # Four disjoint cohorts for the difference of two insertion contrasts.
            probs = np.array(
                [
                    probe["with_extra_occupations"][1],
                    reference["source_only_port_occupations"][1],
                    probe["with_extra_occupations"][0],
                    reference["source_only_port_occupations"][0],
                ]
            )[:, : len(m.fp.TIMES)]
            assert np.all((probs >= -1e-12) & (probs <= 1 + 1e-12))
            counts = rng.binomial(SHOTS, np.clip(probs, 0, 1))
            estimate = (counts[0] - counts[1] - counts[2] + counts[3]) @ dist / SHOTS
            rows.append(
                dict(port=port, counts=counts.tolist(), radial=estimate.tolist())
            )
        records[name] = rows
    return dict(
        protocol_sha256=sha256((m.HERE / "SHARP_READS.md").read_bytes()).hexdigest(),
        counts_order=["high+extra", "high", "low+extra", "low"],
        seed=SEED,
        shots=SHOTS,
        radial_comparisons=family,
        simultaneous_99_percent_half_width=half_width,
        records=records,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    path = m.HERE / "counts.json"
    data = acquire(json.loads((m.HERE / "results.json").read_text()))
    if args.write:
        path.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    else:
        assert data == json.loads(path.read_text())
    print("bound", data["simultaneous_99_percent_half_width"])
    for name in NAMES:
        print(name, data["records"][name][0]["radial"])
