"""Conditional equilibrium of a monitored path of registered material contacts."""

import argparse
import json
from fractions import Fraction
from itertools import combinations
from math import comb, log, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
N = 2**14
WORDS = tuple(sum(1 << j for j in bits) for bits in combinations(range(8), 4))
INDEX = {w: i for i, w in enumerate(WORDS)}
ALL = tuple(range(7))
CUT = (0, 2, 3, 4, 5, 6)
WAITS = (0, 1, 2, 5, 12, 32, 64)
RADIUS = sqrt(log(2 * 1024 / 0.01) / (2 * N))


def physical_path(level):
    board = federation(level)
    ends = np.asarray(board.seams)
    graph = coo_matrix((np.ones(len(ends)), ends.T), shape=(board.n, board.n)).tocsr()
    dist, pred = shortest_path(
        graph + graph.T, indices=0, unweighted=True, return_predecessors=True
    )
    for port in np.flatnonzero(dist >= 7):
        path = [int(port)]
        while path[-1] != 0:
            path.append(int(pred[path[-1]]))
        path = path[::-1]
        for j in range(len(path) - 7):
            block = path[j : j + 8]
            if block[0] // 12 != block[-1] // 12:
                return block
    raise ValueError("No inter-owner path")


def exchange_columns(disabled=False):
    e, c = native_law.local_term_parts()
    h = e + KAPPA * c if not disabled else KAPPA * c
    gate = expm(-1j * np.pi / 2 * h)[:, :9]
    if np.max(np.abs(gate[9:])) > 2e-13:
        raise ValueError("Blank seam leakage")
    return np.abs(gate[:9]) ** 2


def prepare_word(left):
    k = left.bit_count()
    return left | (((1 << (4 - k)) - 1) << 2)


def forward_kernel(edges, disabled=False):
    local = exchange_columns(disabled)
    k = np.eye(len(WORDS))
    for edge in edges:
        gate = np.zeros_like(k)
        for i, w in enumerate(WORDS):
            a, b = (w >> edge) & 1, (w >> (edge + 1)) & 1
            col = native_law.SRC.index((a, b))
            for ab in ((0, 0), (0, 1), (1, 0), (1, 1)):
                p = local[native_law.SRC.index(ab), col]
                if p < 1e-13:
                    continue
                v = (w & ~(3 << edge)) | ((ab[0] + 2 * ab[1]) << edge)
                gate[i, INDEX[v]] += p
        np.testing.assert_allclose(gate.sum(axis=1), 1, atol=1e-13, rtol=0)
        np.testing.assert_allclose(gate, gate.T, atol=1e-13, rtol=0)
        k = k @ gate
    return k


def sweep_maps(edges):
    table = np.empty((1 << len(edges), len(WORDS)), dtype=np.uint8)
    for tape in range(len(table)):
        for index, w in enumerate(WORDS):
            for j, e in enumerate(edges):
                pair = (w >> e) & 3
                if pair in (1, 2):
                    orient = 1 + ((tape >> j) & 1)
                    w = (w & ~(3 << e)) | (orient << e)
            table[tape, index] = INDEX[w]
    return table


def terminal_sample(rng, edges, inventory, bank):
    """Extend the common map one sweep into the past until its image is one word."""
    tape = sweep_maps(edges)
    allowed = np.array(
        [
            i
            for i, w in enumerate(WORDS)
            if len(edges) == 7 or (w & 3).bit_count() == inventory
        ]
    )
    indices = np.full(len(WORDS), -1, dtype=int)
    indices[allowed] = np.arange(len(allowed))
    reduced = indices[tape[:, allowed]]
    if np.min(reduced) < 0:
        raise ValueError("Sector not preserved")
    image = np.tile(np.arange(len(allowed), dtype=np.uint8), (bank, 1))
    counts = np.zeros(16, dtype=np.int64)
    depths = []
    while len(image):
        layers = rng.integers(len(tape), size=len(image))
        image = np.take_along_axis(image, reduced[layers], axis=1)
        done = (image == image[:, :1]).all(axis=1)
        terminal = np.asarray(WORDS)[allowed[image[done, 0]]] & 15
        counts += np.bincount(terminal, minlength=16)
        depths.append(int(done.sum()))
        image = image[~done]
    return counts, depths


def law(reference=False):
    output = []
    mass = (Fraction(3, 8), Fraction(1, 4), Fraction(3, 8))
    for word in range(16):
        k = (word & 3).bit_count()
        exposed_bath_count = (word >> 2).bit_count()
        c = comb(4, 4 - k - exposed_bath_count)
        p = Fraction(c, 70)
        if not reference:
            p *= mass[k] / Fraction(comb(2, k) * comb(6, 4 - k), 70)
        output.append(p)
    return output


def grouped(prob, mode):
    counts = {}
    for w, p in enumerate(prob):
        a = (w & 3).bit_count()
        b = (w >> 2).bit_count()
        key = {
            "pair_totals": (a, b),
            "bath_ports": (a, w >> 2),
            "left_window": (a, w & 7),
            "right_window": (a, w >> 1),
            "full": w,
        }[mode]
        counts[key] = counts.get(key, 0) + p
    return counts


def relative(p, q):
    return sum(float(v) * log(float(v / q[k])) for k, v in p.items() if v)


def run():
    reference, source = law(True), law()
    product = []
    for w in range(16):
        k = (w & 3).bit_count()
        b = (w >> 2).bit_count()
        left = [Fraction(3, 8), Fraction(1, 8), Fraction(1, 8), Fraction(3, 8)][w & 3]
        p = Fraction(4 - k, 6)
        product.append(left * p**b * (1 - p) ** (2 - b))
    modes = ["pair_totals", "bath_ports", "left_window", "right_window", "full"]
    witness = {}
    for view in modes:
        r, p, q = [grouped(x, view) for x in (reference, source, product)]
        witness[view] = {
            "source_divergence": relative(p, r),
            "foil_divergence": relative(q, r),
            "excess": relative(q, p),
        }
        np.testing.assert_allclose(
            relative(q, r) - relative(p, r), relative(q, p), atol=1e-14
        )
    records = []
    for level in range(3):
        rng = np.random.default_rng(12701 + level)
        path = physical_path(level)
        ref, hr = terminal_sample(rng, ALL, None, N)
        bits = rng.integers(2, size=(3, N))
        sectors = bits[0] + np.where(bits[2], bits[0], bits[1])
        active = np.zeros(16, dtype=int)
        depths = {}
        for k in range(3):
            counts, h = terminal_sample(rng, CUT, k, int((sectors == k).sum()))
            active += counts
            depths[str(k)] = h
        records.append(
            dict(
                path=path,
                reference=ref.tolist(),
                reference_depth=hr,
                active=active.tolist(),
                source_sectors=np.bincount(sectors, minlength=3).tolist(),
                active_depth=depths,
            )
        )
    waits = []
    x = np.zeros(len(WORDS))
    weights = [3 / 8, 1 / 8, 1 / 8, 3 / 8]
    for left, weight in enumerate(weights):
        x[INDEX[prepare_word(left)]] = weight
    rng = np.random.default_rng(718243)
    for name, edges, disabled in [
        ("native", CUT, False),
        ("opened", ALL, False),
        ("disabled", CUT, True),
    ]:
        kernel = forward_kernel(edges, disabled)
        for t in WAITS:
            z = x @ np.linalg.matrix_power(kernel, t)
            obs = np.bincount(np.asarray(WORDS) & 15, z, minlength=16)
            c = rng.multinomial(N, obs / obs.sum())
            waits.append(
                dict(mode=name, sweeps=t, counts=c.tolist(), expected=obs.tolist())
            )
    return dict(
        shots=N,
        joint_radius=RADIUS,
        reference=[str(q) for q in reference],
        source=[str(q) for q in source],
        foil=[str(q) for q in product],
        divergences=witness,
        records=records,
        forward=waits,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    ans = run()
    dest = HERE / "results.json"
    if parser.parse_args().write:
        dest.write_text(json.dumps(ans, indent=2, allow_nan=False) + "\n")
    elif dest.read_text() != json.dumps(ans, indent=2, allow_nan=False) + "\n":
        raise ValueError("Replay mismatch")
    print("6 coupled record banks and 21 independent forward banks")
