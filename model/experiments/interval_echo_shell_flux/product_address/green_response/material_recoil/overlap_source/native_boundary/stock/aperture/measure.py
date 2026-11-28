"""Native cumulative record energy on a written trail, before and after load extraction."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy.special import comb

from .. import trail
from ...population import measure as population
from ....closed_region.measure import check

HERE = Path(__file__).resolve().parent
BETA = float(np.log(8 / 7))


def edge_operators(k):
    law = trail.measure.LAW
    _, twice_g, twice_b, h4, _ = population.certificate(k)
    memory_ids = np.arange(11) * 9
    outside = np.setdiff1d(np.arange(99), memory_ids)
    assert not h4[np.ix_(outside, memory_ids)].any()
    vacuum = h4[np.ix_(memory_ids, memory_ids)] / 4
    reflected, d2 = law._label_ops()
    np.testing.assert_array_equal(
        4 * vacuum[1:, 1:], 2 * (np.eye(10) - reflected) + 4 * k * d2
    )
    assert not vacuum[0].any()
    # B_e on the ordinary reduction and the positive write-image credit C_e.
    # Sharp signs have memory energy k+1/2. Read both losses on the same joint density.
    tags = [c for c, (_, d) in enumerate(law.LABELS) if abs(d) == 1]
    pairs = [law.SRC.index(p) for p in ((0, 0), (1, 0), (0, 1))]
    ids = [9 * tag + p for tag in tags for p in pairs]
    credit = 2 * twice_g - h4[9:, 9:]
    loss = (credit + 2 * np.kron(np.eye(10, dtype=int), twice_b))[np.ix_(ids, ids)]
    empty = np.arange(0, 12, 3)
    assert not loss[empty].any()
    # Exact PSD of C_e from its Gram factor W B W*; additionally check its
    # spectrum and the blank term on this physical twelve-state restriction.
    assert np.linalg.eigvalsh(loss).min() > -1e-10
    assert np.linalg.eigvalsh(loss).max() <= 8 * k + 8 + 1e-10
    # The addressed blank detector seam exchanges 10 and 01 in time pi.
    one = [law.SRC.index(p) for p in ((1, 0), (0, 1))]
    hb = h4[np.ix_(one, one)] / 4
    np.testing.assert_allclose(
        abs(expm(-1j * np.pi * hb)) ** 2, [[0, 1], [1, 0]], atol=2e-13
    )
    return vacuum


def vacuum_readings(vacuum):
    labels = trail.measure.LAW.LABELS
    sharp = 1 + labels.index((0, 1))
    es, V = np.linalg.eigh(vacuum)
    rho = (V * np.exp(-BETA * es)) @ V.T / np.exp(-BETA * es).sum()
    modular = (V * (BETA * es + np.log(np.exp(-BETA * es).sum()))) @ V.T
    np.testing.assert_allclose(expm(-modular), rho, atol=2e-14)
    # Direct two-memory calculation of the reduced canonical logarithm.
    joint = np.kron(vacuum, np.eye(11)) + np.kron(np.eye(11), vacuum)
    thermal = expm(-BETA * joint)
    thermal /= np.trace(thermal)
    np.testing.assert_allclose(thermal, np.kron(rho, rho), atol=2e-14)
    np.testing.assert_allclose(
        thermal.reshape(11, 11, 11, 11).trace(axis1=1, axis2=3), rho
    )
    shots = 2**18
    rng = np.random.default_rng(191221)
    panels = []
    for time in (0.0, np.pi / 2, np.pi, 4.0, 16.0):
        psi = expm(-1j * time * vacuum)[:, sharp]
        probs = abs(V.T @ psi) ** 2
        np.testing.assert_allclose(
            [probs[np.isclose(es, k)].sum() for k in (8, 9)], [0.5, 0.5]
        )
        energy = (psi.conj() @ vacuum @ psi).real
        sign_block = np.diag(abs(psi) ** 2)
        dephased = np.trace(vacuum @ sign_block)
        np.testing.assert_allclose([energy, dephased], [8.5] * 2)
        panels.append(
            dict(
                time=time,
                sign_returns=float(abs(psi[sharp]) ** 2),
                mean_energy=float(energy),
                modular_increment=float(
                    (psi.conj() @ modular @ psi).real - modular[0, 0]
                ),
            )
        )
    sampled = []
    for n in (1, 14, 32, 64):
        p = comb(n, np.arange(n + 1)) / 2.0**n
        counts = rng.multinomial(shots, p)
        sampled.append(
            dict(
                records=n,
                native_energies=list(range(8 * n, 9 * n + 1)),
                counts=counts.tolist(),
                exact_modular_increment=BETA * 8.5 * n,
                sampled_modular_increment=BETA
                * (8 * n + counts @ np.arange(n + 1) / shots),
            )
        )
    return dict(
        beta=BETA,
        sector_dimension=11,
        partition_function=np.exp(-BETA * es).sum(),
        spectral_probabilities_per_head=[0.5, 0.5],
        trials_per_size=shots,
        vacuum_connected_memory_load_increment=0,
        panels=panels,
        sampled=sampled,
    )


def alternate(graph, length, rng):
    # Sample a path by a randomized depth-first search, saving its backtracks.
    # The supplied graph contains a path of this length from port zero.
    path = [0]

    def options(x):
        edges = sorted(graph[x])
        return iter(rng.permutation(edges).tolist())

    stack = [options(0)]
    removed = 0
    while len(path) <= length:
        for y in stack[-1]:
            if y not in path:
                path.append(y)
                stack.append(options(y))
                break
        else:
            path.pop()
            stack.pop()
            removed += 1
            if not stack:
                raise ValueError("No eligible route")
    return path, removed


def profile(board, route, terminal, dist, h4):
    edges = {edge: i for i, edge in enumerate(board.seams)}
    marks = np.zeros(len(edges), dtype=int)
    for u, v in zip(route, route[1:]):
        marks[edges[tuple(sorted((u, v)))]] = terminal[int(u > v)]
    idxs = [edges[tuple(sorted((u, v)))] for u, v in zip(route, route[1:])]
    filled, terms = trail.energy(board, route[-1], marks, h4)
    bare, base = trail.energy(board, route[-1], marks * 0, h4)
    excess = (terms - base)[idxs] / 4
    np.testing.assert_array_equal(excess, np.r_[[8.5] * (len(idxs) - 1), -4.25])
    midpoint = (dist[route[:-1]] + dist[route[1:]]) / 2
    panels = []
    for radius in range(int(dist.max()) + 1):
        inside = midpoint <= radius
        half = (
            (dist[route[:-1]] <= radius).astype(int)
            + (dist[route[1:]] <= radius).astype(int)
        ) / 2
        panels.append(
            dict(
                radius=radius,
                inside=int(inside.sum()),
                endpoint_share=half.sum(),
                closed_formation_excess=excess @ inside,
                half_shared_excess=excess @ half,
                evacuated_excess=8.5 * int(inside.sum()),
            )
        )
    return dict(
        path=route,
        terminal_distance=int(dist[route[-1]]),
        midpoints=midpoint.tolist(),
        board_energy=float(filled),
        matched_blank=float(bare),
        apertures=panels,
    )


def run():
    vac = edge_operators(8)
    control = edge_operators(0)
    labels = trail.measure.LAW.LABELS
    mark = 1 + labels.index((0, 1))
    np.testing.assert_allclose([vac[mark, mark], control[mark, mark]], [8.5, 0.5])
    boards = json.loads((trail.HERE / "trail.json").read_text())["boards"]
    W, V, H, _ = trail.measure.kernels()
    terminal = [trail.instrument(f, W, V, H)[1] for f in (True, False)]
    runs = []
    rng = np.random.default_rng(191222)
    for level, row in enumerate(boards):
        board = trail.measure.observations.federation(level)
        graph = nx.Graph(board.seams)
        dist = np.array(
            [
                d
                for _, d in sorted(
                    nx.single_source_shortest_path_length(graph, 0).items()
                )
            ]
        )
        route = row["path"]
        geodesic = profile(board, route, terminal, dist, H)
        np.testing.assert_array_equal(dist[route], np.arange(len(route)))
        n = len(route) - 1
        alternative, back = alternate(graph, n, rng)
        comparison = profile(board, alternative, terminal, dist, H)
        initial = float(Fraction(row["initial_energy"]))
        # In the complete one-load sector H_blank = 48I + L_graph/2.
        # L ≥ 0 and L ≤ 2*maxdegree; both are operator inequalities.
        assert board.degree.max() <= 6
        total = geodesic["board_energy"]
        dmin = 8.5 * n + 48 - total
        dmax = 8.5 * n + 48 + int(board.degree.max()) - total
        assert (dmin, dmax) == (10.25, 16.25)
        detector = initial  # one ordinary load at the equal-degree blank detector port
        runs.append(
            dict(
                owners=row["owners"],
                geodesic=geodesic,
                alternative=comparison,
                search_backtracks=back,
                initial=initial,
                whole_deficit_interval=[dmin, dmax],
                evacuated_board_energy=8.5 * n,
                detector_energy=detector,
                extraction_work=8.5 * n + detector - total,
            )
        )
    return dict(
        boards=runs,
        empty_board=vacuum_readings(vac),
        formation_sharp_gain_at_zero_stiffness=0.5,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    else:
        check(json.loads(path.read_text()), data)
