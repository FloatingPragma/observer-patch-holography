"""The first four time derivatives of the unparked near-minus-far material read."""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import factorial

import networkx as nx
import numpy as np

from . import clock_release as lab

HERE = lab.HERE


def apparatus(board, home, seams, h):
    # The k-th nested commutator with a two-port reading uses only edges
    # touching its (k-1)-hop neighborhood. This metric includes written bonds.
    distance = nx.multi_source_dijkstra_path_length(nx.Graph(board.seams), [12, 13])
    sites = sorted(v for v, d in distance.items() if d <= 4)
    ix = {v: i for i, v in enumerate(sites)}
    edges = [
        (ix[x], ix[y], seams.index((x, y)) if (x, y) in seams else -1)
        for x, y in board.seams
        if min(distance[x], distance[y]) < 4
    ]
    kept = [min(distance[x], distance[y]) < 4 for x, y in seams]
    ports = len(ix)
    padding = [32 * (6 - int(board.degree[v])) for v in sites]
    entries = []
    for col in range(99):
        changes = []
        for k in np.flatnonzero(h[:, col]):
            if k != col:
                head, p = divmod(int(k), 9)
                a, b = divmod(p, 3)
                changes.append((head, a, b, int(h[k, col])))
        entries.append(changes)
    diag = h.diagonal().tolist()
    near, far = ix[13], ix[12]

    def initial(word):
        state = [0] * (ports + 3)
        for v, n in zip(home, word[:4], strict=True):
            if v in ix:
                state[ix[v]] = int(n)
        for j, keep in enumerate(kept):
            if keep:
                state[ports + j] = int(word[4 + j])
        return tuple(state)

    @lru_cache(None)
    def column(state):
        energy = sum(k * n * n for k, n in zip(padding, state[:ports], strict=True))
        out = defaultdict(int)
        for x, y, j in edges:
            mark = 0 if j < 0 else state[ports + j]
            code = 9 * mark + 3 * state[x] + state[y]
            energy += diag[code]
            for new, a, b, value in entries[code]:
                dst = list(state)
                dst[x], dst[y] = a, b
                if j >= 0:
                    dst[ports + j] = new
                else:
                    assert new == 0
                out[tuple(dst)] += value
        out[state] += energy
        return tuple(out.items())

    def act(vec):
        out = defaultdict(int)
        for state, weight in vec.items():
            if weight:
                for k, value in column(state):
                    out[k] += weight * value
        return out

    def obs(vec):
        return {k: v * (k[near] - k[far]) for k, v in vec.items() if k[near] != k[far]}

    return initial, act, obs, near, far, len(ix), len(edges)


def dot(a, b):
    return sum(v * b.get(k, 0) for k, v in a.items())


def derivatives(vec, act, obs):
    n = sum(v * v for v in vec.values())
    first = act(vec)
    second = act(first)
    obs_first = obs(first)
    obs_init = act(obs(vec))
    d2 = 2 * (dot(first, obs_first) - dot(first, obs_init))
    d4 = (
        2 * dot(second, act(obs_init))
        - 8 * dot(second, act(obs_first))
        + 6 * dot(second, obs(second))
    )
    return F(d2, 4**2 * n), F(d4, 4**4 * n)


def run():
    _, w, cells, *_ = lab.lab.energy.pair.setup()
    board, home, marks, *_ = lab.lab.energy.returned(w)
    h = lab.lab.energy.pair.stock.kernels()[2]
    start, act, obs, near, far, ns, ne = apparatus(board, home, marks, h)
    shift = (lab.lab.energy.least_change(cells) - 1) / len(w)
    weights = defaultdict(F)
    for word, p in zip(w, shift, strict=True):
        weights[start(word)] += p
    witness = json.loads((HERE / "clock_release.json").read_bytes())["readers"]
    # Bound the full generator on *all* words by the sum of the native bond
    # norms (Gershgorin) and private maxima. No truncated-board flight is
    # used for the finite-time remainder.
    bond = max(np.abs(h).sum(1)) / 4
    Hbound = int(np.ceil(len(board.seams) * bond + 32 * sum(6 - board.degree)))
    rows = {}
    for t in (1, 2):
        values = [[F(), F()], [F(), F()]]
        for state, p in weights.items():
            if not p:
                continue
            ref = derivatives({state: 1}, act, obs)
            a, b = list(state), list(state)
            a[near] = b[far] = t
            for k, s in enumerate((1, -1)):
                response = derivatives({tuple(a): 1, tuple(b): s}, act, obs)
                for j in range(2):
                    values[k][j] += p * (response[j] - ref[j]) / t
        row = {}
        for parity, (d2, d4) in zip(("even", "odd"), values, strict=True):
            assert d2 == F(witness[str(t)][parity]["repair_increment"][0])
            # The difference of two preparations on probe and control has
            # trace norm at most four; ||n_near-n_far||<=2. Odd orders
            # vanish for the real native generator and real sharp preparations.
            first = d2 / 2 if d2 else d4 / 24
            order = 2 if d2 else 4
            remainder = F(8 * (2 * Hbound) ** (order + 2), t * factorial(order + 2))
            limit = (float(abs(first) / (2 * remainder))) ** 0.5 if first else 0
            safe = F(1, 10)
            while first and remainder * safe**2 > abs(first) / 2:
                safe /= 10
            row[parity] = dict(
                second=str(d2),
                fourth=str(d4),
                first_nonzero_order=order if first else None,
                uniform_certified_nonzero_time=str(safe) if first else None,
                half_margin_time_estimate=limit,
            )
        rows[str(t)] = row
    return dict(
        sha256={
            "clock_release.json": sha256(
                (HERE / "clock_release.json").read_bytes()
            ).hexdigest()
        },
        covered_sites=ns,
        covered_bonds=ne,
        full_H_norm_upper=Hbound,
        readers=rows,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "clock_release_time.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        lab.lab.check(json.loads(dst.read_bytes()), result)
    print(result["readers"])
