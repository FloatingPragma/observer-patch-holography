"""Compare two independently prepared owner placements by bulk and boundary spectroscopy."""

import argparse
import json
from math import log, sqrt
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier

HERE = Path(__file__).resolve().parent
SHOTS = 2**42
law = candidate.native_law


def local_reader(kappa):
    exchange, contrast = law.local_term_parts()
    energies, vectors = np.linalg.eigh(exchange + kappa * contrast)
    levels = np.rint(energies).astype(int)
    np.testing.assert_allclose(energies, levels, rtol=0, atol=2e-13)
    memories = np.zeros((2, law.DIM // law.N_SRC))
    memories[0, 0] = 1
    empty = law.SRC.index((0, 0))
    memories[1] = law.echo()[:, empty].reshape(-1, law.N_SRC)[:, empty]
    # The 64 uniformly spaced exposure arms resolve these integer levels uniquely.
    count = 64
    phases = np.exp(-2j * np.pi * np.outer(np.arange(count), energies) / count)
    clock_bins = abs(np.fft.ifft(phases, axis=0)) ** 2
    np.testing.assert_allclose(clock_bins, np.eye(count)[:, levels], atol=2e-13)
    ans = np.zeros((2, 3, 3, count))
    for written, memory in enumerate(memories):
        for idx, (a, b) in enumerate(law.SRC):
            seed = np.kron(memory, np.eye(law.N_SRC)[:, idx])
            ans[written, a, b] = clock_bins @ abs(vectors.conj().T @ seed) ** 2
    return ans


def geometry(level):
    car = carrier.federation(level)
    edges = np.array(car.seams)
    owners = edges // 12
    external = owners[:, 0] != owners[:, 1]
    i, j = owners[external].T
    graph = csr_matrix((np.ones(len(i)), (i, j)), shape=(car.cells, car.cells))
    graph = graph + graph.T
    dist = shortest_path(graph, indices=0, directed=False, unweighted=True).astype(int)
    old = int(np.flatnonzero(dist == 1)[0])
    gateways = car.degree.reshape(-1, 12) == 6
    same = np.all(gateways == gateways[old], axis=1)

    def farthest(match):
        eligible = np.flatnonzero(match & (dist > 1))
        return int(eligible[np.argmax(dist[eligible])])

    return car, edges, external, dist, gateways, old, [farthest(same), farthest(~same)]


def relocate(edges, loads, memory, old, new):
    ports = np.arange(loads.size).reshape(-1, 12)
    before, after = ports[[old, new]]
    perm = ports.copy()
    perm[[old, new]] = perm[[new, old]]
    perm = perm.ravel()
    index = {tuple(edge): i for i, edge in enumerate(edges)}
    w = memory.copy()
    for group in (before, after):
        for a, b in carrier.CELL_SEAMS:
            i, j = group[[a, b]]
            w[index[i, j]] = memory[index[perm[i], perm[j]]]
    return loads[perm], w


def meter(rng, table, loads, memory, edges, mask, private, kappa):
    i, j = edges[mask].T
    # Aggregate exactly equivalent local preparations before drawing the whole cohort.
    codes = 9 * memory[mask] + 3 * loads[i] + loads[j]
    mass = np.bincount(codes, minlength=18) @ table.reshape(18, 64)
    private_levels = kappa * loads[private] ** 2
    mass += np.bincount(private_levels, minlength=64)
    n = len(i) + len(private_levels)
    if n == 0:
        return dict(terms=0, counts=[0] * 64)
    prob = mass / n
    # Rounding collapses machine-epsilon sidelobes; every declared output is archived.
    prob = np.round(prob, 14)
    if min(prob) < 0 or abs(prob.sum() - 1) > 1e-12:
        raise ArithmeticError("Invalid energy meter")
    return dict(terms=n, counts=rng.multinomial(SHOTS, prob / prob.sum()).tolist())


def read(bank, width):
    count = np.array(bank["counts"])
    if bank["terms"] == 0:
        return 0.0, 0.0
    if count.sum() != SHOTS:
        raise ValueError("Incomplete cohort")
    value = bank["terms"] * (count @ np.arange(64)) / SHOTS
    err = bank["terms"] * width * (sqrt(log(2 * 200 / 0.01) / (2 * SHOTS)) + 1e-12)
    return float(value), float(err)


def score(banks, kappa):
    work = {}
    width = 1 + 4 * kappa
    for source in ("whole", "exterior", "private"):
        readings = [read(b[source], width) for b in banks]
        work[source] = [
            readings[1][0] - readings[0][0],
            readings[1][1] + readings[0][1],
        ]
    wh, ex, pr = (np.array(work[k]) for k in ("whole", "exterior", "private"))
    return dict(
        work=work,
        boundary_residual=float(wh[0] - ex[0] - pr[0]),
        joint_allowance=float(wh[1] + ex[1] + pr[1]),
        zero_work_excluded=bool(abs(wh[0]) > wh[1]),
        exterior_only_excluded=bool(abs(wh[0] - ex[0]) > wh[1] + ex[1]),
    )


def acquire(sampler, tab, initial, target, edges, exterior, priv, marked, kappa):
    banks = []
    for state, collars in [initial, target]:
        banks.append(
            {
                "whole": meter(
                    sampler,
                    tab,
                    state,
                    collars,
                    edges,
                    np.ones(len(edges), bool),
                    priv,
                    kappa,
                ),
                "exterior": meter(
                    sampler, tab, state, collars, edges, exterior, [], kappa
                ),
                "private": meter(
                    sampler,
                    tab,
                    state,
                    collars,
                    edges,
                    np.zeros(len(edges), bool),
                    marked,
                    kappa,
                ),
            }
        )
    return banks


def run(private_control=False):
    reader = {k: local_reader(k) for k in (8, 0)}
    sampler = np.random.default_rng(433029 if private_control else 10022037)
    rows = []
    for level in range(3):
        car, edges, external, dist, gates, old, targets = geometry(level)
        prep = np.random.default_rng(221037 + level)
        loads = prep.integers(0, 3, car.n)
        memory = prep.integers(0, 2, len(edges))
        for cell, value in zip([old, *targets], [2, 0, 0], strict=True):
            loads[12 * cell : 12 * (cell + 1)] = value
        priv = np.flatnonzero(car.degree == 5)
        for new in targets[1:] if private_control else targets:
            changed = np.flatnonzero(gates[old] != gates[new])
            exterior = external & np.isin(edges // 12, [old, new]).any(axis=1)
            # Matched private slots cancel after permutation. Inspect only the unmatched slots.
            marked = np.intersect1d(
                priv, np.concatenate([12 * c + changed for c in (old, new)])
            )
            if private_control:
                for owner in (old, new):
                    loads[owner * 12 + changed] = 2 * gates[owner, changed]
            target = relocate(edges, loads, memory, old, new)
            for kappa, tab in reader.items():
                banks = acquire(
                    sampler,
                    tab,
                    (loads, memory),
                    target,
                    edges,
                    exterior,
                    priv,
                    marked,
                    kappa,
                )
                rows.append(
                    dict(
                        cells=car.cells,
                        old=old,
                        new=new,
                        distance=int(dist[new]),
                        kappa=kappa,
                        external_terms=int(exterior.sum()),
                        mismatched_labels=len(changed),
                        relocation_bound=int(
                            (1 + 4 * kappa) * exterior.sum() + 4 * kappa * len(changed)
                        ),
                        banks=banks,
                        **score(banks, kappa)
                    )
                )
    return dict(starts_per_nonempty_cohort=SHOTS, rows=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--private-control", action="store_true")
    args = parser.parse_args()
    dst = args.output
    out = run(args.private_control)
    if dst:
        dst.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        committed = "private.json" if args.private_control else "records.json"
        compare(out, json.loads((HERE / committed).read_bytes()), "endpoint work")
    for row in out["rows"]:
        print(
            row["cells"],
            row["distance"],
            row["kappa"],
            row["mismatched_labels"],
            row["work"]["whole"],
            row["boundary_residual"],
            row["joint_allowance"],
        )
