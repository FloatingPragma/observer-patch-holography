"""Integer enclosures for one mobile load on arbitrary stationary empty writes."""

import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import sympy as sp

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation

HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "results.json"


def check(condition, message):
    if not condition:
        raise ArithmeticError(message)


def native_columns():
    """Exact conjugated generator on blank and vacuum-written seam factors."""
    law = candidate.native_law
    check(candidate.KAPPA == 8, "registered contrast changed")
    pattern = sp.Matrix(law._write_pattern().tolist()) / sp.sqrt(2)
    U = sp.eye(law.DIM)
    U[:9, :9] = sp.zeros(9)
    U[:9, 9:] = pattern.T
    U[9:, :9] = pattern
    U[9:, 9:] = sp.eye(90) - pattern * pattern.T
    r, d2 = law._label_ops()
    memory = (sp.eye(10) - sp.Matrix(r.tolist())) / 2 + 8 * sp.Matrix(d2.tolist())
    g = sp.diag(sp.zeros(9), sp.kronecker_product(memory, sp.eye(9)))
    direct = U * g * U.T
    blank = sp.eye(law.DIM)[:, :9]
    vacuum = law.SRC.index((0, 0))
    record = (U * blank[:, vacuum]).reshape(11, 9)[:, vacuum]
    check(record.dot(record) == 1 and record[0] == 0, "invalid empty-write factor")
    checks = []
    for a, b in ((0, 0), (1, 0), (0, 1), (2, 0), (0, 2)):
        i = law.SRC.index((a, b))
        j = law.SRC.index((b, a))
        load = sp.eye(9)[:, i]
        written = sp.kronecker_product(record, load)
        expected = (blank[:, i] - blank[:, j]) / 2 + 8 * (a - b) ** 2 * blank[:, i]
        check(
            direct * blank[:, i] == expected and direct * written == sp.zeros(99, 1),
            "full 99-state invariant-column identity failed",
        )
        checks.append([a, b])
    return checks


def bands(written_degrees, amount, max_degree):
    counts = Counter(int(d) for d in written_degrees)
    return [
        [s, counts[s], -8 * amount**2 * s, -8 * amount**2 * s + max_degree - s]
        for s in sorted(counts, reverse=True)
    ]


def audit_cell(board, edges, amount):
    cut = set(tuple(sorted(e)) for e in edges)
    degree = np.zeros(board.n, dtype=int)
    h = np.zeros((board.n, board.n))
    for a, b in board.seams:
        if (a, b) in cut:
            degree[[a, b]] += 1
        else:
            h[a, a] += 0.5
            h[b, b] += 0.5
            h[a, b] -= 0.5
            h[b, a] -= 0.5
    h -= np.diag(8 * amount**2 * degree)
    e = np.linalg.eigvalsh(h)
    enclosure = bands(degree, amount, int(max(board.degree)))
    offset = 0
    for _, multiplicity, lo, hi in enclosure:
        values = e[offset : offset + multiplicity]
        check(
            values.min() >= lo - 1e-11 and values.max() <= hi + 1e-11,
            "eigenvalue outside rank band",
        )
        offset += multiplicity
    check(
        offset == board.n
        and np.count_nonzero(e < -0.000001) == np.count_nonzero(degree),
        "negative index differs from the number of touched ports",
    )
    return enclosure


def report():
    local = native_columns()
    small = cell()
    fan = [e for e in small.seams if 0 in e]
    # Every subset of the five-seam fan includes isolated and intersecting writes.
    motifs = [list(c) for n in range(len(fan) + 1) for c in combinations(fan, n)]
    motifs.extend([list(small.seams), small.seams[::2], small.seams[1::3]])
    for motif in motifs:
        for a in (1, 2):
            audit_cell(small, motif, a)
    cases = []
    for level in (0, 1, 2):
        board = federation(level)
        check(max(board.degree) <= 6, "registered degree cap failed")
        for name, W in [
            ("blank", []),
            ("one", board.seams[:1]),
            ("star", [e for e in board.seams if 0 in e]),
            ("all_external", board.external),
            ("all_written", board.seams),
        ]:
            touches = np.bincount(np.asarray(W, dtype=int).ravel(), minlength=board.n)
            cases.append(
                dict(
                    cells=board.cells,
                    pattern=name,
                    writes=len(W),
                    negative_index=int(np.count_nonzero(touches)),
                    bands_by_amount={
                        str(a): bands(touches, a, int(max(board.degree)))
                        for a in (1, 2)
                    },
                )
            )
    return dict(
        native_columns=local,
        isolated_cell_patterns=len(motifs),
        full_cell_spectra_checked=2 * len(motifs),
        cases=cases,
        universal_bands_in_written_degree_order={
            str(a): [[lo, hi] for _, _, lo, hi in reversed(bands(range(7), a, 6))]
            for a in (1, 2)
        },
        reference="Blank-board uniform one-load energy 6*kappa*a^2. Whole-owner and extensive patterns use that same comparison.",
    )


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path)
    args = cli.parse_args()
    result = report()
    if args.output:
        with args.output.open("x") as f:
            f.write(json.dumps(result, indent=2) + "\n")
    else:
        check(
            result == json.loads(RECEIPT.read_text()),
            "core-pattern bounds disagree with deposit",
        )
        print(
            "Verified the empty-write band inequalities and counts for one-load motion."
        )


if __name__ == "__main__":
    main()
