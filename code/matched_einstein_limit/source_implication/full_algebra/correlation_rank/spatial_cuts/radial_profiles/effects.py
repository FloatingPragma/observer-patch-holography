"""Enclose the four-tick native radial and regional material transfer readings."""

import argparse
from fractions import Fraction
import json
from math import factorial, isqrt
from pathlib import Path

import numpy as np
import sympy as sp

from .. import measure as cuts

HERE = Path(__file__).parent
BITS = 45
BASE = 1 << BITS
LIMBS = 5
UNIT = BASE**LIMBS
DEGREE = 640
SCALE = 10**58


def coefficients():
    # Four time units on the integer 2H at norm at most 64: z=128.
    denom = factorial(DEGREE)
    assert Fraction(128 ** (DEGREE + 1), factorial(DEGREE + 1)) < Fraction(1, 10**70)
    result = []
    for k in range(DEGREE + 1):
        term = denom // factorial(k) * 64**k
        total = term
        for m in range(1, (DEGREE - k) // 2 + 1):
            term = term * 64**2 // (m * (m + k))
            total += -term if m % 2 else term
        c = (1 if k == 0 else 2) * total * UNIT // denom
        assert abs(c) <= 3 * UNIT
        result.append(c)
    return result


def columns(h, ends, starts, coefficient):
    """Evaluate the bounded Chebyshev recurrence with integer limbs, low first."""
    d = h.shape[0]
    row = int(abs(h).sum(axis=1).max())
    assert row < 64 and row * (BASE + 1) < 2**61
    eps = Fraction(3 * (1 + isqrt(d)) * DEGREE * (DEGREE + 1) * (DEGREE + 2), 6 * UNIT)
    eps += Fraction(3 * (DEGREE + 1) + 2, UNIT) + Fraction(1, 10**70)
    assert 4 * eps < Fraction(1, 10**56)
    # High limb is signed; lower limbs encode floor-division remainders.
    parts = np.zeros((LIMBS + 1, d, len(starts)), dtype=np.int64)
    parts[-1, starts] = np.eye(len(starts), dtype=np.int64)
    prev = parts.copy()
    out = np.zeros((2, len(ends), len(starts)), dtype=object)
    for k, c in enumerate(coefficient):
        val = parts[-1, ends].astype(object)
        for limb in reversed(parts[:-1]):
            val = val * BASE + limb[ends]
        out[k % 2] += (1 if k % 4 in (0, 3) else -1) * c * val
        if k == DEGREE:
            break
        assert np.abs(parts[-1]).max() < 8
        prod = np.asarray([h @ limb for limb in parts])
        for i in range(LIMBS):
            prod[i + 1] += prod[i] >> BITS
            prod[i] &= BASE - 1
        shift = 6 if k == 0 else 5
        step = prod >> shift
        step[:-1] += (prod[1:] % (1 << shift)) << (BITS - shift)
        if k:
            step -= prev
            for i in range(LIMBS):
                step[i + 1] += step[i] >> BITS
                step[i] &= BASE - 1
        prev, parts = parts, step
    return out // UNIT


def read_rows():
    native, reg = cuts.native, cuts.regional
    atlas = native.Atlas(1)
    board = atlas.board
    patch = native.chosen(atlas)
    route = reg.bridge.native.iso.pairing.formation.region(board)
    start = reg.bridge.beam.indices(sum(reg.bridge.beam.screens(board, route), []))
    edges = native.flux_on_edges(board, patch)
    # F_d is total outward transfer of cells at exactly graph distance d.
    shells = {d: {u for u in patch["cells"] if atlas.dist[0][u] == d} for d in range(6)}
    assert set.union(*shells.values()) == set(patch["cells"])
    shell_ops = [
        sum(op for (u, _), op in edges.items() if u in cells)
        for cells in shells.values()
    ]
    cells = sorted(patch["cells"])
    node_ops = [
        sum(op for (u, _), op in edges.items() if u == target) for target in cells
    ]
    bases = [reg.cert.prior.sym_matrix(b) for b in reg.cert.iso.basis(True)]
    co = coefficients()
    rows = []
    for charge in (-1, 1):
        channels = {}
        for odd in (False, True):
            h, ops = native.op_table(board, route, patch, charge, odd)
            assert (sum(shell_ops) - ops["boundary"][0]).count_nonzero() == 0
            native_ops = [(name, *ops[name]) for name in reg.NAMES]
            native_ops.extend((str(i), op, 4) for i, op in enumerate(shell_ops))
            native_ops.extend((f"cell-{u}", op, 4) for u, op in zip(cells, node_ops))
            idx = np.unique(
                np.concatenate([op.nonzero()[0] for _, op, _ in native_ops])
            )
            w = columns(h, idx, start, co)
            for name, op, unit in native_ops:
                ids = np.unique(op.nonzero()[0])
                sites = len(set(ids // 4))
                row_norm = int(abs(op).sum(axis=1).max()) / abs(unit)
                assert sites <= 960 and row_norm <= 200
                block = reg.region.compress(
                    w[:, np.searchsorted(idx, ids)], op[ids, :][:, ids], unit
                )
                # compress uses the parent fixed-point denominator. Restore the one used here.
                block *= sp.Rational(reg.bridge.native.prior.FRAC**2, UNIT**2)
                # Each exact compressed 4-by-4 port block errs by <1e-56.
                # For n ports the stacked error is sqrt(n)*1e-56. The exact
                # column is contractive. Pauli trace norm two. We use n in
                # place of sqrt(n), outward, also bounding the quadratic term.
                err = sp.Rational(5 * sites * max(1, int(np.ceil(row_norm))), 10**56)
                entries = channels.setdefault(name, [])
                for b in bases:
                    val = reg.cert.inner(block, b)
                    entries.append(
                        [
                            int(sp.floor((val - err) * SCALE)),
                            int(sp.ceiling((val + err) * SCALE)),
                        ]
                    )
        rows.extend(
            dict(charge=charge, name=name, row=entries)
            for name, entries in channels.items()
        )
    return dict(scale=SCALE, shells=[sorted(x) for x in shells.values()], rows=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = read_rows()
    if args.write:
        lines = [json.dumps(r, separators=(",", ":")) for r in data["rows"]]
        header = json.dumps({k: v for k, v in data.items() if k != "rows"})[:-1]
        (HERE / "effects.json").write_text(
            header + ', "rows": [\n' + ",\n".join(lines) + "\n]}\n"
        )
    else:
        assert data == json.loads((HERE / "effects.json").read_text())
    print(
        "Verified the separately enrolled native shell and regional four-tick readings."
    )
