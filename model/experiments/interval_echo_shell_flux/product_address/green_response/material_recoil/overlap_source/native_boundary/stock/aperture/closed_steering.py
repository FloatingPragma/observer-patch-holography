"""Reach prescribed re-pairings by native stochastic appointments on the same inventory."""

import argparse
from hashlib import sha256
import json
from math import lcm

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

from . import closed_energy as energy

pair = energy.pair
HERE = pair.HERE
TRIALS = 2**15


def actions(full):
    W, V, *_ = pair.stock.kernels()
    _, _, _, edges, *_ = pair.stock.carrier()
    lookup = {tuple(r): i for i, r in enumerate(full)}
    tickets = np.empty((len(full), 6, 64), np.int16)
    adj = np.zeros((len(full), len(full)), bool)
    for i, w in enumerate(full):
        for j, (a, b) in enumerate(edges):
            old = 9 * w[4 + j] + 3 * w[a] + w[b]
            for k, mat in enumerate((W, V)):
                ends = np.flatnonzero(mat[:, old])
                ids = []
                for col in ends:
                    out = w.copy()
                    out[4 + j], word = divmod(col, 9)
                    out[a], out[b] = divmod(word, 3)
                    ids.append(lookup[tuple(out)])
                tickets[i, 2 * j + k] = np.repeat(ids, mat[ends, old])
                adj[i, ids] = True
    assert np.array_equal(adj, adj.T)
    hops = shortest_path(csr_matrix(adj), directed=False, unweighted=True)
    assert np.isfinite(hops).all() and hops.max() < 127
    return tickets, hops.astype(np.int8), lookup


def run():
    full, w, cells, _, _, _, _, _ = pair.setup()
    tickets, distance, lookup = actions(full)
    z = energy.least_change(cells)
    accepted = np.array([lookup[tuple(x)] for x in w])
    local = np.full(len(full), -1, int)
    local[accepted] = np.arange(len(w))
    denom = lcm(*(x.denominator for x in z))
    thresholds = np.array([min(denom, int(x * denom)) for x in z])
    surplus = np.array([max(0, int((x - 1) * denom)) for x in z])
    rng = np.random.default_rng(63285107)
    initial = rng.integers(len(full), size=TRIALS)
    old = local[initial]
    ok = old >= 0
    move = ok & (rng.integers(denom, size=TRIALS) >= thresholds[old])
    current = initial[move].copy()
    target = accepted[
        np.searchsorted(
            surplus.cumsum(), rng.integers(surplus.sum(), size=sum(move)), side="right"
        )
    ]
    ages = np.zeros(len(target), int)
    rejected_steps = 0
    trace = sha256()
    counts = np.zeros(6, int)
    while (active := np.flatnonzero(current != target)).size:
        before = current[active]
        goal = target[active]
        d = distance[before, goal]
        routes = tickets[before]
        future = distance[routes, goal[:, None, None]]
        # Lexicographic choice: lowest mean distance among appointments
        # admitting a nearer outcome, then most nearer tickets and action ID.
        nearer = (future < d[:, None, None]).sum(axis=2)
        score = np.where(nearer > 0, 65 * future.sum(axis=2) - nearer, 10**8)
        act = score.argmin(axis=1)
        assert np.all(nearer[np.arange(len(active)), act] > 0)
        counts += np.bincount(act, minlength=6)
        ticket = rng.integers(64, size=len(active))
        after = routes[np.arange(len(active)), act, ticket]
        rejected_steps += int(sum(local[after] < 0))
        trace.update(
            np.column_stack((active, before, act, ticket, after))
            .astype("<i4")
            .tobytes()
        )
        current[active] = after
        ages[active] += 1
    end = initial.copy()
    end[move] = current
    assert np.array_equal(initial[~ok], end[~ok])
    assert np.all(full[initial, :4].sum(1) == full[end, :4].sum(1))
    # All three clocks and loads are read on the same identified run.
    t0 = np.bincount(cells[old[ok]], minlength=4)
    t1 = np.bincount(cells[local[end[ok]]], minlength=4)
    before = full[initial[ok], 1]
    after = full[end[ok], 1]
    d = (after != 1).astype(int) - (before != 1)
    dy = (cells[local[end[ok]]] == 0).astype(int) - (cells[old[ok]] == 0)
    assert np.array_equal(d, dy)

    def bond(n):
        return 2 * (n != 0) + 32 * n * n

    assert np.array_equal(bond(after) - bond(before), 31 * d + 65 * (after - before))
    assert np.array_equal(
        2 * ((after != 0).astype(int) - (before != 0)), (after - before) - d
    )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "CLOSED_STEERING.md", HERE / "closed_energy.json")
        },
        trials=TRIALS,
        source_rejects=int(sum(~ok)),
        accepted=int(sum(ok)),
        retained_without_move=int(sum(ok & ~move)),
        steered=int(sum(move)),
        mean_appointments_per_move=float(ages.mean()),
        longest_realized_path=int(ages.max()),
        visits_outside_accepted_set=rejected_steps,
        total_appointments=int(ages.sum()),
        appointments_by_action=counts.tolist(),
        graph_diameter=int(distance.max()),
        trace_sha256=trace.hexdigest(),
        original_cells=t0.tolist(),
        final_cells=t1.tolist(),
        middle_nonunit_increment=int(d.sum()),
        mean_load_change=float((after - before).mean()),
        bond_increment_accepted=float((bond(after) - bond(before)).mean() / 4),
        accepted_half_echo_change=float(
            ((after != 0).astype(int) - (before != 0)).mean()
        ),
        identity_violations=0,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "closed_steering.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        assert json.loads(path.read_bytes()) == out
    print({k: v for k, v in out.items() if k not in ("sha256", "trace_sha256")})
