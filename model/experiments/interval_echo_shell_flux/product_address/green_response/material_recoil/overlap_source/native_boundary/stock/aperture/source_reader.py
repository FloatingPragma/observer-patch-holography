"""Infer nested modular increments from independent reference and source counts."""

import argparse
from collections import Counter
from hashlib import sha256
import json

import networkx as nx
import numpy as np
from scipy.linalg import expm

from . import measure as native

HERE = native.HERE
MULT = np.array([1, 2, 2, 2, 2, 2])
LABELS = ("blank", "written_zero", "one_even", "one_odd", "two_even", "two_odd")
SETTINGS = ("blank", "outward", "folded", "folded_reversed")
TRIALS = 2**14
REFERENCE = 2**26


def spectrometer():
    labels = native.trail.measure.LAW.LABELS
    eye = np.eye(11)
    groups = [[eye[0]], [eye[1 + labels.index((c, 0))] for c in (0, 1)]]
    for d in (1, 2):
        for hand in (1, -1):
            groups.append(
                [
                    (
                        eye[1 + labels.index((c, d))]
                        + hand * eye[1 + labels.index((c, -d))]
                    )
                    / np.sqrt(2)
                    for c in (0, 1)
                ]
            )
    return np.array([np.array(g).T @ np.array(g) for g in groups])


def probabilities(k):
    effects = spectrometer()
    h = native.edge_operators(k)
    rho = expm(-native.BETA * h)
    rho /= rho.trace()
    ref = np.einsum("kij,ji->k", effects, rho)
    # Checking that the actual reference is scalar in each registered subspace
    # is essential if the count decoder is to recover its quantum logarithm.
    np.testing.assert_allclose(
        np.einsum("k,kij->ij", ref / MULT, effects), rho, atol=1e-13
    )
    u = expm(-4j * h)
    idx = native.trail.measure.LAW.LABELS
    stocks = {}
    for sign in (-1, 0, 1):
        v = u[:, 0 if sign == 0 else 1 + idx.index((0, sign))]
        p = np.einsum("i,kij,j->k", v.conj(), effects, v).real
        np.testing.assert_allclose(
            p, [1, 0, 0, 0, 0, 0] if not sign else [0, 0, 0.5, 0.5, 0, 0], atol=1e-13
        )
        p[abs(p) < 1e-14] = 0
        stocks[sign] = p / sum(p)
    return ref, stocks


def geometry():
    source = json.loads((HERE / "results.json").read_bytes())
    board = next(b for b in source["boards"] if b["owners"] == 80)
    carrier = native.trail.measure.observations.federation(1)
    graph = nx.Graph(carrier.seams)
    dist = nx.single_source_shortest_path_length(graph, 0)
    shells = Counter(dist[a] + dist[b] for a, b in carrier.seams)
    marks = {}
    for tag, layout in (("outward", "geodesic"), ("folded", "alternative")):
        path = board[layout]["path"]
        marks[tag] = Counter(
            dist[a] + dist[b] for a, b in zip(path[:-1], path[1:], strict=True)
        )
        assert sum(marks[tag].values()) == 32
    return shells, marks


def emit(k, seed, shells, marks):
    rng = np.random.default_rng(seed)
    ref, stocks = probabilities(k)
    counts = {}
    for name in SETTINGS:
        rows = []
        for shell, size in sorted(shells.items()):
            active = marks.get(
                name, marks["folded"] if name == "folded_reversed" else {}
            ).get(shell, 0)
            # The instrument is insensitive to the orientation of a sharp sign.
            sign = -1 if name == "folded_reversed" else 1
            row = rng.multinomial(TRIALS * active, stocks[sign])
            row += rng.multinomial(TRIALS * (size - active), stocks[0])
            rows.append([shell, size, row.tolist()])
        counts[name] = rows
    return dict(reference=rng.multinomial(REFERENCE, ref).tolist(), counts=counts)


def decode(bank, tail):
    """Uses the tape, spectral multiplicities and repeat counts only."""
    ref = np.array(bank["reference"]) / REFERENCE
    eps = np.sqrt(tail / (2 * REFERENCE))
    p, lo, hi = (x / MULT for x in (ref, ref - eps, ref + eps))
    assert min(lo) > 0
    gap = np.log(p[0] / p[2:4])
    lower = np.log(lo[0] / hi[2:4])
    upper = np.log(hi[0] / lo[2:4])
    output = {}
    for name, rows in bank["counts"].items():
        shells = []
        for shell, size, row in rows:
            assert sum(row) == size * TRIALS and sum(row[4:]) + row[1] == 0
            n = row[2] + row[3]
            if n == 0:
                shells.append([shell / 2, 0.0, 0.0, 0.0])
                continue
            assert n % TRIALS == 0
            written = n / TRIALS
            f = row[3] / n
            err = np.sqrt(tail / (2 * n))
            ends = [max(f - err, 0), min(f + err, 1)]
            mean = written * ((1 - f) * gap[0] + f * gap[1])
            low = written * min((1 - x) * lower[0] + x * lower[1] for x in ends)
            high = written * max((1 - x) * upper[0] + x * upper[1] for x in ends)
            shells.append([shell / 2, mean, low, high])
        table = np.array(shells)
        enclosed = [table[table[:, 0] <= r, 1:].sum(axis=0).tolist() for r in (12, 24)]
        # Bound the same two-radius contrast directly; each shell's uncertainty
        # enters once even if it was read in both balls.
        coeff = (table[:, 0] <= 24).astype(int) - 2 * (table[:, 0] <= 12)
        low = np.where(coeff >= 0, table[:, 2], table[:, 3]) @ coeff
        high = np.where(coeff >= 0, table[:, 3], table[:, 2]) @ coeff
        weights = np.array(
            [np.maximum(0, (r * r - table[:, 0] ** 2) / (2 * r)) for r in (12, 24)]
        )
        weighted = weights @ table[:, 1:]
        inversion = (48 * weights[1] - 24 * weights[0]) / (24**2 - 12**2)
        output[name] = dict(
            shells=shells,
            balls=enclosed,
            outer_minus_twice_inner=[low, high],
            diamond_weighted_balls=weighted.tolist(),
            two_radius_inversion=(inversion @ table[:, 1:]).tolist(),
        )
    return output


def run():
    sizes, marks = geometry()
    # Both laws, both labs, all reference labels and every possible nonempty bank.
    tail = np.log(2 * 4 * (len(LABELS) + len(SETTINGS) * len(sizes)) / 0.01)
    cases = {}
    for i, k in enumerate((8, 0)):
        for j, lab in enumerate(("A", "B")):
            seed = 94107 + 2 * i + j
            bank = emit(k, seed, sizes, marks)
            cases[f"{k}_{lab}"] = dict(
                seed=seed, records=bank, decoded=decode(bank, tail)
            )
    return dict(
        protocol_sha256=sha256((HERE / "SOURCE_TEST.md").read_bytes()).hexdigest(),
        native_sha256=sha256((HERE / "results.json").read_bytes()).hexdigest(),
        labels=LABELS,
        multiplicities=MULT.tolist(),
        shots_per_memory=TRIALS,
        reference_shots=REFERENCE,
        radii=[12, 24],
        cases=cases,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    out = json.loads(json.dumps(run(), allow_nan=False))
    path = HERE / "source_reader.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2) + "\n")
    else:
        native.check(json.loads(path.read_bytes()), out)
    for name, case in out["cases"].items():
        print(name, case["decoded"]["folded"]["balls"])
