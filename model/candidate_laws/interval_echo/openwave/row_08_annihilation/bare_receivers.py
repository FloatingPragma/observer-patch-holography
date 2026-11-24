"""Free wall collisions with two initially uniform receiving seams."""

import argparse
import gzip
import json

import numpy as np
from scipy import sparse
from scipy.special import jv

from model.executable_model import federation
from . import experiment as ex, walls
from .verify import compare

N = 12
POWERS = 3 ** np.arange(N - 1, -1, -1)
SHOTS = 2**27
TIMES = (0, 0.5, 2, 8, 32)
CLASSES = ("ordered", "reversed", "two_walls", "four_walls", "other")
SEED = 81136


def geometry():
    ring, anchors, wings = walls.geometry()
    sites = ring + wings[0][:2] + wings[1][:2]
    edges = [(j, (j + 1) % 8) for j in range(8)]
    edges += [(anchors[0], 8), (8, 9), (anchors[1], 10), (10, 11)]
    board = federation(0)
    assert all(sites[q] in board.neighbours[sites[p]] for p, q in edges)
    return sites, edges, 6 - board.degree[sites]


def basis(mixing):
    codes = np.arange(3**N)
    words = (codes[:, None] // POWERS % 3).astype(np.int8)
    if mixing:
        use = words.sum(axis=1) == N
    else:
        use = ((words == 0).sum(axis=1) == 4) & ((words == 2).sum(axis=1) == 4)
    return codes[use], words[use]


def generator(codes, words, k, mixing):
    _, edges, missing = geometry()
    term = ex.term(k)
    local = term.copy()
    if mixing:
        for j in (2, 6):
            local[j, 4] = local[4, j] = 0.5 / np.sqrt(2)
    lookup = np.full(3**N, -1)
    lookup[codes] = np.arange(len(codes))
    diagonal = k * (words.astype(float) ** 2 @ missing)
    rows, cols, data = [], [], []
    for a, b in edges:
        incoming = 3 * words[:, a] + words[:, b]
        diagonal += local[incoming, incoming]
        for dest, src in np.argwhere(np.triu(abs(local) > 1e-12, 1)):
            keep = np.flatnonzero(incoming == src)
            new = (
                codes[keep]
                + (dest // 3 - src // 3) * POWERS[a]
                + (dest % 3 - src % 3) * POWERS[b]
            )
            target = lookup[new]
            assert np.all(target >= 0)
            rows.extend([target, keep])
            cols.extend([keep, target])
            data.extend([np.full(keep.size, local[dest, src])] * 2)
    length = np.arange(len(codes))
    h = sparse.coo_matrix(
        (
            np.concatenate([diagonal, *data]),
            (np.concatenate([length, *rows]), np.concatenate([length, *cols])),
        ),
        shape=(len(length), len(length)),
    ).tocsr()
    return h


def flights(h, x, times):
    # A Gershgorin enclosure makes the Chebyshev argument valid for every eigenvalue.
    diag = h.diagonal()
    radius = np.asarray(abs(h).sum(axis=1)).ravel() - abs(diag)
    low, high = np.min(diag - radius), np.max(diag + radius)
    mid, half = (high + low) / 2, (high - low) / 2
    scaled = (h - mid * sparse.eye(len(x), format="csr")) / half
    order = int(
        np.ceil(half * max(times) + 18 * (half * max(times) + 1) ** (1 / 3) + 40)
    )
    orders = np.arange(order + 1)
    weight = (
        2 * (-1j) ** (orders % 4) * jv(orders[None, :], np.array(times)[:, None] * half)
    )
    weight[:, 0] /= 2
    # Monitor the discarded coefficient scale and check the propagated norm and energy.
    assert np.max(abs(weight[:, -4:])) < 1e-19
    out = np.zeros((len(times), *x.shape), complex)
    p, q = x, scaled @ x
    for i, t in enumerate(times):
        out[i] = weight[i, 0] * p + weight[i, 1] * q
    for n in range(2, order + 1):
        nxt = 2 * (scaled @ q) - p
        p, q = q, nxt
        for i in np.flatnonzero(abs(weight[:, n]) > 1e-17):
            out[i] += weight[i, n] * q
    return out, order


def analyzers():
    table = np.zeros((9, 9))
    labels = []
    for a in range(3):
        for b in range(a, 3):
            for sign in (1,) if a == b else (1, -1):
                i = len(labels)
                labels.append((a, b, sign))
                table[3 * a + b, i] = 1 if a == b else 1 / np.sqrt(2)
                if a != b:
                    table[3 * b + a, i] = sign / np.sqrt(2)
    np.testing.assert_allclose(table.T @ table, np.eye(9), atol=1e-14)
    return labels, table


def centers():
    words = np.arange(3**8)[:, None] // POWERS[4:] % 3
    order = np.asarray(walls.loads(walls.reference(0)))
    classes = np.full(len(words), 4)
    valid = (words[:, ::2] + words[:, 1::2] == 2).all(axis=1) & (words != 1).all(axis=1)
    bits = words[:, ::2] == 2
    wall_count = (bits == np.roll(bits, -1, axis=1)).sum(axis=1)
    classes[valid & (wall_count == 2)] = 2
    classes[valid & (wall_count == 4)] = 3
    classes[(words == order).all(axis=1)] = 0
    classes[(words == 2 - order).all(axis=1)] = 1
    return classes


def read(codes, state, k, rng):
    full = np.zeros(3**N, complex)
    full[codes] = state
    full = full.reshape(-1, 9, 9)
    labels, table = analyzers()
    probabilities = abs(table.T @ full @ table) ** 2
    center = centers()
    bins = np.stack([probabilities[center == j].sum(axis=0) for j in range(5)])
    np.testing.assert_allclose(bins.sum(), 1, atol=3e-10, rtol=0)
    weights = bins.ravel().round(11)
    rank = weights.argsort(kind="stable")
    sorted_weights = weights[rank].copy()
    sorted_weights[-1] = 1 - sorted_weights[:-1].sum()
    counts = np.empty(weights.size, dtype=int)
    counts[rank] = rng.multinomial(SHOTS, sorted_weights)
    counts = counts.reshape(bins.shape)
    energy = np.array([k * (a - b) ** 2 + (parity == -1) for a, b, parity in labels])
    bright = energy > 0
    neutral = np.array([a + b == 2 for a, b, _ in labels])
    outside = []
    positions = abs(full) ** 2
    for i in (0, 1):
        port = np.arange(9) % 3 != 1
        outside.append(
            float(
                positions[:, port, :].sum() if i == 0 else positions[:, :, port].sum()
            )
        )
    return dict(
        probability=bins.round(11).tolist(),
        count=counts.tolist(),
        ordered=int(counts[:2].sum()),
        ordered_and_any_energy=int(
            counts[:2, bright, :].sum() + counts[:2, ~bright, :][:, :, bright].sum()
        ),
        both_neutral_bright=int(
            counts[:, bright & neutral, :][:, :, bright & neutral].sum()
        ),
        both_neutral_bright_probability=float(
            bins[:, bright & neutral, :][:, :, bright & neutral].sum()
        ),
        order_probability=float(bins[:2].sum()),
        outside_material_probability=outside,
        outside_material_counts=[int(rng.binomial(SHOTS, min(1, p))) for p in outside],
        mean_seam_energy=[
            float((bins.sum(axis=(0, 2)) @ energy)),
            float((bins.sum(axis=(0, 1)) @ energy)),
        ],
    )


def run():
    sites, edges, _ = geometry()
    rng = np.random.default_rng(SEED)
    cases = []
    for k, mixing in ((8, False), (0, False), (8, True)):
        times = (0.5,) if mixing else TIMES
        codes, words = basis(mixing)
        h = generator(codes, words, k, mixing)
        x = np.zeros((len(words), 2))
        for i in range(2):
            bits = tuple(
                z ^ int(i and j >= 2) for j, z in enumerate(walls.reference(0))
            )
            word = walls.loads(bits) + [1] * 4
            x[np.searchsorted(codes, word @ POWERS), i] = 1
        initial = np.sum(x * (h @ x), axis=0)
        states, order = flights(h, x, times)
        samples = []
        for time, both in zip(times, states, strict=True):
            np.testing.assert_allclose(
                np.sum(both.conj() * (h @ both), axis=0), initial, atol=3e-8, rtol=0
            )
            samples.append(
                dict(time=time, preparations=[read(codes, s, k, rng) for s in both.T])
            )
        cases.append(
            dict(
                contrast=k,
                pair_mixing=mixing,
                dimension=len(words),
                initial_energy=initial.tolist(),
                expansion_order=order,
                readings=samples,
            )
        )
    return dict(
        sites=sites,
        edges=edges,
        shots=SHOTS,
        seed=SEED,
        ring_classes=CLASSES,
        detector_labels=analyzers()[0],
        cases=cases,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    path = ex.HERE / "bare_receivers.json.gz"
    args = parser.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    if args.write:
        path.write_bytes(
            gzip.compress(
                (json.dumps(result, separators=(",", ":")) + "\n").encode(), mtime=0
            )
        )
    else:
        compare(
            result, json.loads(gzip.decompress(path.read_bytes())), "bare_receivers"
        )
        print(
            "Verified the collision, all remote energy and material records, and both alternative laws."
        )


if __name__ == "__main__":
    main()
