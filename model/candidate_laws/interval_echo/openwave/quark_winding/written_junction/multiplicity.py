"""Exact lower bounds on the number of distinct stationary energies in a written block."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from . import backgrounds, reservoir


PRIME = 65521


def recurrence(sequence):
    """Shortest recurrence over the prime field, via Berlekamp–Massey."""
    size = len(sequence)
    c = np.zeros(size + 1, dtype=np.int64)
    c[0] = 1
    prior = c.copy()
    degree, shift, prior_delta = 0, 1, 1
    for k in range(size):
        delta = (
            int(sequence[k] + c[1 : degree + 1] @ sequence[k - degree : k][::-1])
            % PRIME
        )
        if not delta:
            shift += 1
            continue
        factor = delta * pow(prior_delta, -1, PRIME) % PRIME
        promote = 2 * degree <= k
        if promote:
            old = c.copy()
        c[shift:] = (c[shift:] - factor * prior[:-shift]) % PRIME
        if promote:
            degree, prior, prior_delta, shift = k + 1 - degree, old, delta, 1
        else:
            shift += 1
    return c[: degree + 1]


def certificate(h):
    m = (2 * h).tocsr()
    assert np.array_equal(m.data, m.data.round())
    m = m.astype(np.int64)
    assert (m != m.T).nnz == 0
    dimension = m.shape[0]
    rng = np.random.default_rng(12013005)
    u = rng.integers(1, PRIME, size=dimension, dtype=np.int64)
    v = u.copy()
    # Values in a dot product are reduced residues. Check the largest exact
    # dot product and sparse update both fit their integer accumulator.
    assert (PRIME - 1) ** 2 * dimension < np.iinfo(np.int64).max
    assert abs(m).sum(axis=1).max() * (PRIME - 1) < np.iinfo(np.int64).max
    seq = np.zeros(2 * dimension, dtype=np.int64)
    for i in range(len(seq)):
        seq[i] = int(u @ v) % PRIME
        v = (m @ v) % PRIME
    c = recurrence(seq)
    degree = len(c) - 1
    assert all(
        int(c @ seq[i - degree : i + 1][::-1]) % PRIME == 0
        for i in range(degree, len(seq))
    )
    return dict(
        dimension=dimension,
        distinct_energies_at_least=degree,
        maximum_multiplicity_at_most=dimension - degree + 1,
    )


def run(scans):
    result = {}
    for name, pair in scans["paired"].items():
        stats = []
        for row in pair["selected"]:
            length = row["length"]
            g = backgrounds.bits(row["gauss"], length)
            ref = scans["geometry"][str(length)]["reference"]
            h = backgrounds.generator(
                backgrounds.block(g, np.array(ref)), name == "extension"
            )
            s = dict(length=length, gauss=row["gauss"], ring=certificate(h))
            if row["n"] == length // 3:
                _, spokes, edges, faces, geo = reservoir.patch(
                    length, pair["additional_site"]
                )
                w, n, _, _ = reservoir.configurations(g, spokes)
                coupled = reservoir.generator(
                    w, n, edges, faces, np.array(geo["reference"]), name == "extension"
                )
                s["reservoir"] = certificate(coupled)
            stats.append(s)
        result[name] = stats
    first = scans["paired"]["extension"]["selected"][0]
    data = backgrounds.block(
        backgrounds.bits(first["gauss"], 12),
        np.array(scans["geometry"]["12"]["reference"]),
    )
    plain = backgrounds.generator(data, True)
    return dict(
        prime=PRIME,
        laws=result,
        three_identical_copies=certificate(sparse.block_diag([plain] * 3)),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    s = json.loads(Path(__file__).with_name("backgrounds.json").read_text())
    result = run(s)
    with args.output.open("x") as f:
        f.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
