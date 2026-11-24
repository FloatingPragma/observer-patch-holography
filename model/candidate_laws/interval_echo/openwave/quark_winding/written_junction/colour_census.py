"""Search every mobile twelve-port background for a possible exact colour multiplet."""

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy import sparse

from . import backgrounds, multiplicity


def modular_rank(a):
    a = a.copy() % multiplicity.PRIME
    row = 0
    for col in range(a.shape[1]):
        choices = np.flatnonzero(a[row:, col])
        if not len(choices):
            continue
        idx = int(choices[0]) + row
        a[[row, idx]] = a[[idx, row]]
        inv = pow(int(a[row, col]), -1, multiplicity.PRIME)
        pivot = (inv * a[row, col:]) % multiplicity.PRIME
        a[row + 1 :, col:] = (
            a[row + 1 :, col:] - a[row + 1 :, col, None] * pivot
        ) % multiplicity.PRIME
        row += 1
        if row == a.shape[0]:
            break
    return row


def capacity(h):
    cert = multiplicity.certificate(h)
    bound = cert["maximum_multiplicity_at_most"]
    if bound <= 2:
        return bound
    m = (2 * h).astype(np.int64)
    n = m.shape[0]
    rng = np.random.default_rng(12013005)
    blocks = []
    for k in range(1, min(bound, 5) + 1):
        v = rng.integers(1, multiplicity.PRIME, n)
        cols = []
        for _ in range(n):
            cols.append(v)
            v = (m @ v) % multiplicity.PRIME
        blocks.append(np.column_stack(cols))
        if k >= 2 and modular_rank(np.column_stack(blocks)) == n:
            return k
    return bound


def exact_lines(h):
    ev = np.linalg.eigvalsh(h.toarray())
    groups = np.split(ev, np.flatnonzero(np.diff(ev) > 1e-8) + 1)
    lines = []
    for a in groups:
        if len(a) < 3:
            continue
        val = round(2 * a[0])
        if abs(2 * a[0] - val) > 1e-8:
            lines.append(
                dict(approximate_energy=float(a[0]), numerical_multiplicity=len(a))
            )
            continue
        shifted = (2 * h - val * sparse.eye(h.shape[0])).astype(np.int64)
        matrix = sp.polys.matrices.DomainMatrix.from_Matrix(
            sp.SparseMatrix(shifted.toarray())
        ).to_field()
        lines.append(
            dict(twice_energy=val, exact_multiplicity=h.shape[0] - matrix.rank())
        )
    return lines


def run(scans):
    ref = np.array(scans["geometry"]["12"]["reference"])
    result = {}
    for ext, name in ((False, "native"), (True, "extension")):
        hist = Counter()
        possible = []
        for scan in scans["census"]:
            for row in scan["rows"]:
                g = backgrounds.bits(row[0], 12)
                data = backgrounds.block(g, ref)
                h = backgrounds.generator(data, ext)
                bound = capacity(h)
                hist[bound] += 1
                if bound >= 3:
                    possible.append(
                        dict(
                            gauss=row[0],
                            dimension=h.shape[0],
                            bound=bound,
                            tight_rows=row[2],
                            lines=exact_lines(h),
                        )
                    )
        result[name] = dict(
            multiplicity_upper_bound_histogram={
                str(k): int(v) for k, v in sorted(hist.items())
            },
            eligible_for_triplet=possible,
        )
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    s = json.loads(Path(__file__).with_name("backgrounds.json").read_text())
    data = run(s)
    with args.output.open("x") as f:
        f.write(json.dumps(data, indent=2, allow_nan=False) + "\n")
