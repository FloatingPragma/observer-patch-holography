"""How long does the measured overlap-source contrast survive native write service?"""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from . import collective_source as source

tidal = source.tidal
HERE = tidal.HERE
GRID = 2**51


def setup():
    full, selected, cells, *_ = tidal.lab.energy.pair.setup()
    lookup = {tuple(w): i for i, w in enumerate(full)}
    original = np.full(len(full), Q(0), object)
    repaired = original.copy()
    initial = np.array([lookup[tuple(w)] for w in selected])
    original[initial] = Q(1, len(selected))
    repaired[initial] = tidal.lab.energy.least_change(cells) / len(selected)
    W, V, *_ = tidal.lab.energy.pair.stock.kernels()
    _, _, _, edges, *_ = tidal.lab.energy.pair.stock.carrier()
    service = []
    random = []
    for s, (a, b) in enumerate(edges):
        matrices = []
        for matrix in (W, V):
            dest = []
            start = []
            weights = []
            for i, word in enumerate(full):
                c = 9 * word[4 + s] + 3 * word[a] + word[b]
                for end in np.flatnonzero(matrix[:, c]):
                    w = word.copy()
                    w[4 + s], m = divmod(end, 9)
                    w[a], w[b] = divmod(m, 3)
                    dest.append(lookup[tuple(w)])
                    start.append(i)
                    weights.append(matrix[end, c])
            mat = csr_matrix(
                (weights, (dest, start)), shape=(len(full),) * 2, dtype=np.int64
            )
            assert (mat != mat.T).nnz == 0 and np.all(np.asarray(mat.sum(0)) == 64)
            matrices.append(mat)
        Wm, Vm = matrices
        service.append((Wm + Vm, Vm))
        random.append((Wm + Vm).tocsc())
    assert connected_components(sum(s[0] for s in service), directed=False)[0] == 1
    assert all(np.all(mat.diagonal() > 0) for mat, _ in service)
    # Fix the exact source under exchange-only service before comparing writes.
    for arm in (original, repaired):
        denom = math.lcm(*(x.denominator for x in arm))
        ticks = np.array([int(p * denom) for p in arm], np.int64)
        assert all(np.array_equal(v @ ticks, 64 * ticks) for _, v in service)
    return full, original, repaired, service, random


def gates(full):
    """Project each recorded sharp word onto a stationary all-port read
    after returning its heads and loads to the open 240-station federation.
    """
    table = tidal.fastest_relay()
    prob = {}
    columns = []
    for w in full:
        pop = np.bincount(w[:4], minlength=3)
        pop[0] += tidal.lab.V - 4
        pop = tuple(map(int, pop))
        p = []
        for h in w[4:]:
            label = tidal.lab.stationary_heads(pop)[h]
            key = (pop, label)
            if key not in prob:
                counts = tidal.lab.pre(pop, label)
                sums = table @ counts
                prob[key] = [Q(int(z), 32 * int(counts.sum())) for z in sums]
            p.append(prob[key])
        columns.append(p)
    return np.array(columns, object)


def draw(series, ticket_tables, initial, rng, observed):
    n = 2**17
    w = rng.choice(len(initial), n, p=np.array(initial, float))
    outputs = {}
    for step in range(33):
        if step in (0, 1, 4, 16, 32):
            # A local yes/no read after the open-board preparation.
            outputs[str(step)] = int((rng.random(n) < observed[w]).sum())
        if step == 32:
            break
        for j, gate in enumerate(series):
            r = rng.integers(128, size=n)
            # For every occupied starting column pick an actual output ticket.
            if step == 0:
                assert np.all(np.asarray(gate.sum(0)) == 128)
            lookup = ticket_tables[j]
            w = lookup[w, r]
    return outputs, n


def run():
    full, original, repaired, service, random = setup()
    mask = np.array(original > 0, bool)
    reading = gates(full)
    states = np.array([original, repaired], float).T
    read = np.array(reading, float)
    flags = mask.astype(float)
    killed = states.copy()
    raw_delta = repaired - original
    delta = np.array([round(x * GRID) for x in raw_delta], np.int64)
    assert int(max(abs(delta))) * 128 + 64 < 2**63
    # Every update is column- and row-stochastic. Rounding the signed
    # measure costs at most one half grid unit per coordinate and step.
    # Dobrushin contraction carries the L1 error to every later read.
    levels = (0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096)
    out = {}
    for sweep in range(max(levels) + 1):
        if sweep in levels:
            err = Q((1 + 3 * sweep) * len(full), 2 * GRID)
            tv = (Q(sum(map(abs, map(int, delta))), GRID) + err) / 2
            arm = states.T @ read[:, 2, :]
            gap = (states[:, 1] - states[:, 0]) @ read.reshape(len(states), 6)
            assert np.max(abs(gap)) <= float(tv) + 2e-13
            out[str(sweep)] = dict(
                total_variation_upper=str(tv),
                unconditioned_third_relay=arm.tolist(),
                repaired_minus_original_all_links=gap.reshape(3, 2).tolist(),
                eligible_share=(flags @ states).tolist(),
                continuously_eligible=killed.sum(0).tolist(),
                eligible_third_relay=(
                    states.T @ (read[:, 2, :] * flags[:, None])
                ).tolist(),
            )
        if sweep == max(levels):
            break
        for mat, _ in service:
            states = mat @ states / 128
            killed = (mat @ killed / 128) * flags[:, None]
            delta = (mat @ delta + 64) // 128
    np.testing.assert_allclose(states.sum(0), 1, rtol=1e-12)
    lookup = []
    for gate in random:
        lut = []
        for i in range(len(full)):
            a, b = gate.indptr[i : i + 2]
            lut.append(np.repeat(gate.indices[a:b], gate.data[a:b]))
        lookup.append(np.stack(lut))
    simulations = {}
    for j, p in enumerate((original, repaired)):
        bank, n = draw(
            random, lookup, p, np.random.default_rng(304251 + j), read[:, 2, 0]
        )
        eps = float(np.sqrt(np.log(2 * 10 / 0.01) / (2 * n)))
        for step, count in bank.items():
            assert abs(count / n - out[step]["unconditioned_third_relay"][j][0]) <= eps
        simulations[("original", "source")[j]] = dict(
            trials=n, counts=bank, joint99_error=eps
        )
    projection = reading[:, 2, 0]
    stationary = sum(projection, Q(0)) / len(full)
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "COLLECTIVE_TIDE.md", HERE / "closed_steering.json")
        },
        accessible_words=len(full),
        eligible_words=int(mask.sum()),
        unique_stationary_relay=str(stationary),
        native_readout={
            k: str((w @ projection))
            for k, w in (("original", original), ("source", repaired))
        },
        appointments_per_sweep=3,
        watched=out,
        direct_draws=simulations,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "collective_retention.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        assert json.loads(dst.read_bytes()) == result
    print(
        {k: float(Q(v["total_variation_upper"])) for k, v in result["watched"].items()}
    )
