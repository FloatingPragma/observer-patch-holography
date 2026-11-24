"""Test all near-source shortest routes and count sharp and polarized arrival times."""

import argparse
import json

import networkx as nx
import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier

from .experiment import (
    HERE,
    frame,
    generator,
    path_and_gates,
    qualifies,
    spectral_groups,
)


TIMES = np.arange(129) / 2
SHOTS = 2**25
THRESHOLDS = (0.001, 0.005)
NREAD = 2 * 4 * 3 * len(TIMES)
EPS = float(np.sqrt(np.log(2 * NREAD / 0.01) / (2 * SHOTS)))


def routes(board):
    graph, gates = nx.Graph(), {}
    for a, b in board.external:
        u, v = a // 12, b // 12
        graph.add_edge(u, v)
        gates[u, v], gates[v, u] = a % 12, b % 12
    paths = sorted(
        tuple(p)
        for owner, distance in nx.single_source_shortest_path_length(
            graph, 0, cutoff=3
        ).items()
        if distance
        for p in nx.all_shortest_paths(graph, 0, owner)
    )
    return paths, gates


def crossing(p, threshold):
    # With a simultaneous uniform band, the first crossing must lie between
    # the first possible crossing and the first certain crossing on the grid.
    out = []
    for diff in (EPS, 0, -EPS):
        found = np.flatnonzero(p + diff >= threshold)
        out.append(None if len(found) == 0 else float(TIMES[found[0]]))
    return out


def flight(energy, vectors, board, rng, bank):
    path, gates = path_and_gates(board)
    source = np.column_stack([np.eye(12)[:, gates[0]], frame(gates[0])])
    spec = vectors[:12].T @ source
    batch = np.exp(-1j * np.outer(TIMES, energy))[:, :, None] * spec[None, :, :]
    measurements = {}
    for d, owner in enumerate(path[1:], 1):
        v = vectors[12 * owner : 12 * owner + 12]
        local = v @ batch
        prob = (abs(local) ** 2).sum(axis=1)
        counts = rng.binomial(SHOTS, prob.clip(0, 1).round(13))
        bank[f"n{board.cells}_d{d}_counts"] = counts
        measurements[str(d)] = {
            label: [crossing(counts[:, j] / SHOTS, q) for q in THRESHOLDS]
            for j, label in enumerate(("sharp", "x", "y", "z"))
        }
    return dict(path=path, readings=measurements)


def direction(energy, vectors, board):
    paths, gates = routes(board)
    groups = spectral_groups(energy)
    starts = [g[0] for g in groups]
    # Coherent erasure of M equally weighted quarter-tick emission times.
    # The geometric sum is stable at exact resonance and close eigenvalues.
    frequencies = energy[starts]
    detuning = (frequencies[:, None] - energy[None, :]) / 4
    filters = {}
    for arms in (256, 512):
        delta = (detuning + np.pi) % (2 * np.pi) - np.pi
        filters[arms] = (
            np.exp(0.5j * (arms - 1) * delta)
            * np.sinc(arms * delta / (2 * np.pi))
            / np.sinc(delta / (2 * np.pi))
        )
    prep = {}
    finite_passages = {arms: {} for arms in filters}
    passages, records, above, ideal = {}, [], set(), set()
    for p in paths:
        first, last = gates[0, p[1]], gates[p[-1], p[-2]]
        if first not in prep:
            prep[first] = vectors[:12].T @ frame(first)
        receive = frame(last).T @ vectors[12 * p[-1] : 12 * p[-1] + 12]
        product = receive.T[:, :, None] * prep[first][:, None, :]
        tr = np.add.reduceat(product, starts, axis=0)
        singular = np.linalg.svd(tr, compute_uv=False)
        transverse = np.linalg.svd(tr[:, :2, :2], compute_uv=False)
        lon = np.maximum(
            np.linalg.norm(tr[:, 2, :], axis=1), np.linalg.norm(tr[:, :, 2], axis=1)
        )
        data = np.column_stack([transverse, lon])
        accepted = np.flatnonzero(qualifies(data)).tolist()
        bright = np.flatnonzero(singular[:, 1] >= 0.001)
        above.update(map(int, bright))
        ideal.update(accepted)
        passages[p] = set(accepted)
        finite = {}
        for arms, filt in filters.items():
            avg = (filt @ product.reshape(len(energy), 9)).reshape(-1, 3, 3)
            sv = np.linalg.svd(avg[:, :2, :2], compute_uv=False)
            leak = np.maximum(
                np.linalg.norm(avg[:, 2, :], axis=1),
                np.linalg.norm(avg[:, :, 2], axis=1),
            )
            coords = np.column_stack([sv, leak])
            yes = np.flatnonzero(qualifies(coords))
            finite[str(arms)] = [[int(i), *coords[i].round(9).tolist()] for i in yes]
            finite_passages[arms][p] = set(yes)
        records.append(
            dict(
                path=list(p),
                port_axes=[first, last],
                pass_lines=accepted,
                finite_filter_pass=finite,
                possible_two_mode_lines=[
                    dict(
                        index=int(i),
                        frequency=round(float(energy[starts[i]]), 10),
                        multiplicity=len(groups[i]),
                        full=singular[i].round(9).tolist(),
                        transverse_longitudinal=data[i].round(9).tolist(),
                    )
                    for i in bright
                ],
            )
        )
    all_prefixes = [
        list(p)
        for p in paths
        if len(p) == 4 and bool(set.intersection(*(passages[p[:k]] for k in (2, 3, 4))))
    ]
    return dict(
        paths=records,
        pass_at_all_three=all_prefixes,
        finite_pass_at_all_three={
            str(arms): [
                list(p)
                for p in paths
                if len(p) == 4
                and bool(set.intersection(*(rays[p[:k]] for k in (2, 3, 4))))
            ]
            for arms, rays in finite_passages.items()
        },
        eligible_frequencies=[
            round(float(energy[starts[i]]), 9) for i in sorted(above)
        ],
        passing_frequencies=[round(float(energy[starts[i]]), 9) for i in sorted(ideal)],
    )


def run():
    rng = np.random.default_rng(290612)
    bank = {}
    out = []
    for level in (0, 1):
        board = carrier.federation(level)
        for material in ("unit", "vacancy"):
            e, v = eigh(generator(board, material, candidate.KAPPA), driver="evd")
            row = dict(
                cells=board.cells, material=material, survey=direction(e, v, board)
            )
            if material == "unit":
                row["flight"] = flight(e, v, board, rng, bank)
            out.append(row)
    return dict(
        boards=out,
        shots=SHOTS,
        simultaneous_99_probability_error=EPS,
        launches=["sharp", "x", "y", "z"],
        thresholds=THRESHOLDS,
        times=TIMES.tolist(),
    ), bank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    rows, bank = run()
    p = HERE / "survey.json"
    q = HERE / "flight.npz"
    if args.write:
        p.write_text(json.dumps(rows, indent=2, allow_nan=False) + "\n")
        np.savez_compressed(q, **bank)
    else:
        # A degeneracy's projector is invariant under changes of its numerical basis.
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(json.loads(json.dumps(rows)), json.loads(p.read_text()), "survey")
        with np.load(q) as saved:
            assert set(saved) == set(bank)
            for k, a in bank.items():
                np.testing.assert_array_equal(a, saved[k])
    print("Verified the complete three-step direction survey and held flight readings.")


if __name__ == "__main__":
    main()
