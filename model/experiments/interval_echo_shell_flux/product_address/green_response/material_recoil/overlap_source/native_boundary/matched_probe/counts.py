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


def acquire(data, recorded=None):
    """Measure the radial contrast or decode an existing set of count banks."""
    rng = np.random.default_rng(SEED)
    _, _, dist, origins = m.fp.make_board()
    family = len(NAMES) * len(origins) * len(m.fp.TIMES)
    half_width = sqrt(4 * sum(dist**2) * log(2 * family / 0.01) / (2 * SHOTS))
    # Bernstein's bound checks all individual cohort counts against the law.
    log_tail = log(2 * 4 * len(dist) * family / 1e-6)
    records = {}
    for name in NAMES:
        reference = data["standards"][name]
        rows = []
        for i, (port, probe) in enumerate(zip(origins, reference["probes"], strict=True)):
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
            probs = np.clip(probs, 0, 1)
            if recorded is None:
                counts = rng.binomial(SHOTS, probs)
            else:
                counts = np.asarray(recorded["records"][name][i]["counts"])
            assert counts.shape == probs.shape
            assert np.issubdtype(counts.dtype, np.integer)
            assert np.all((0 <= counts) & (counts <= SHOTS))
            deviation = (
                np.sqrt(2 * SHOTS * probs * (1 - probs) * log_tail) + 2 * log_tail / 3
            )
            np.testing.assert_array_less(abs(counts - SHOTS * probs), deviation)
            estimate = (counts[0] - counts[1] - counts[2] + counts[3]) @ dist / SHOTS
            np.testing.assert_allclose(
                estimate, probe["radial"][: len(m.fp.TIMES)], atol=half_width, rtol=0
            )
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
    source = json.loads((m.HERE / "results.json").read_text())
    if args.write:
        data = acquire(source)
        path.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    else:
        archived = json.loads(path.read_text())
        data = acquire(source, recorded=archived)
        m.fp.check(archived, data)
    print("bound", data["simultaneous_99_percent_half_width"])
    for name in NAMES:
        print(name, data["records"][name][0]["radial"])
