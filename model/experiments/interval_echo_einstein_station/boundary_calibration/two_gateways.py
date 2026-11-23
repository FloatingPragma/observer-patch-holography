"""Native coherent incidence across two gateways and reproducible routed joint readouts."""

import argparse
import hashlib
import json
from collections import deque
from itertools import combinations
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .measure import (
    HERE,
    SEED,
    SHOTS,
    TIMES,
    compile_wave,
    federation,
    inputs,
    operators,
    pair_index,
    product_ball,
    seed_packet,
)
from .phase_receiver import evaluate


def outside_route(board, left, right):
    parents = {left: None}
    queue = deque([left])
    while right not in parents:
        current = queue.popleft()
        for neighbor in sorted(board.neighbours[current]):
            if neighbor >= 12 and neighbor not in parents:
                parents[neighbor] = current
                queue.append(neighbor)
    route = [right]
    while route[-1] != left:
        route.append(parents[route[-1]])
    return list(reversed(route))


def local_carrier(board, sites):
    """Keep registered degree-reference phase marks on the isolated connected read graph."""
    index = {int(p): j for j, p in enumerate(sites)}
    neighbors = [[index[q] for q in board.neighbours[p] if q in index] for p in sites]
    seams = [(a, b) for a, row in enumerate(neighbors) for b in row if a < b]
    return SimpleNamespace(
        n=len(sites), degree=board.degree[sites], neighbours=neighbors, seams=seams
    )


def components(a, b):
    cross = a @ b.conj().T
    return ((a @ a.conj().T + b @ b.conj().T) / 2, (cross - cross.conj().T) / (2j))


def joint_blocks(branches, first, second, visible, route, board):
    n = len(visible)
    labels = np.full(board.n, -1, dtype=int)
    labels[visible] = np.arange(n)
    aa, bb = labels[first], labels[second]
    single = (aa >= 0) ^ (bb >= 0)
    local = np.maximum(aa, bb)[single]
    other = np.where(aa[single] >= 0, second[single], first[single])
    single_waves = np.zeros((2, n, board.n), complex)
    single_waves[:, local, other] = branches[single].T
    both = (aa >= 0) & (bb >= 0)
    pair_waves = np.zeros((2, n * (n - 1) // 2, 1), complex)
    pair_waves[:, pair_index(n, aa[both], bb[both]), 0] = branches[both].T
    empty = ~(single | both)
    blocks = [
        components(*branches[empty].T[:, None, :]),
        components(*single_waves),
        components(*pair_waves),
    ]
    # The trapped source cell and every other external site outside the route are traced out.
    single_waves[:, :, route] = 0.0
    return blocks, components(*single_waves)


def rail_effects(block):
    mean, slope = block
    n = len(mean)
    settings = []
    for a in range(n):
        settings.append(
            ([a], "port", [float(mean[a, a].real), float(slope[a, a].real)])
        )
    for a, b in combinations(range(n), 2):
        for name, z in (("real", 1.0), ("circular", 1j)):
            effect = [
                (m[a, a].real + m[b, b].real) / 2 + (z * m[a, b]).real for m in block
            ]
            settings.append(([a, b], name, effect))
    return settings


def measure(level):
    board = federation(level)
    dark, gates = inputs(board)
    first, second = product_ball(board, gates)
    kinetic, contrast, keys = operators(board, first, second)
    packet = seed_packet(keys, board.n, dark, gates)
    step = -0.25j * (kinetic + contrast)
    pairs = []
    for i, j in combinations(range(3), 2):
        a, b = int(gates[i]), int(gates[j])
        visible = np.sort(
            np.r_[
                np.arange(a - a % 12, a - a % 12 + 12),
                np.arange(b - b % 12, b - b % 12 + 12),
            ]
        )
        route = outside_route(board, a, b)
        auxiliary = [v for v in route if v not in visible]
        sites = np.sort(np.unique(np.r_[visible, route])).astype(int)
        reader = local_carrier(board, sites)
        probe = np.zeros(board.n, complex)
        probe[a], probe[b] = 1 / np.sqrt(2), 1j / np.sqrt(2)
        np.testing.assert_allclose(compile_wave(board, probe), probe, atol=1e-12)
        pairs.append((i, j, visible, auxiliary, sites, reader))
    approaches = ("incoherent", "coherent")
    records, traces = {k: [] for k in approaches}, {k: [] for k in approaches}
    for t in TIMES:
        packet = expm_multiply(step, packet, traceA=step.diagonal().sum())
        for i, j, visible, aux, sites, reader in pairs:
            a, b = packet[:, 2 * i : 2 * i + 2], packet[:, 2 * j : 2 * j + 2]
            contributions = [
                joint_blocks(v, first, second, visible, aux, board)
                for v in (a, b, (a + 1j * b) / np.sqrt(2))
            ]
            both = [
                ((x + y) / 2, (u + v) / 2)
                for (x, u), (y, v) in zip(contributions[0][0], contributions[1][0])
            ]
            routed = tuple(
                (x + y) / 2 for x, y in zip(contributions[0][1], contributions[1][1])
            )
            for k, (blocks, read_block) in zip(
                approaches, [(both, routed), contributions[2]], strict=True
            ):
                ledgers = dict(
                    gateways=[int(gates[i]), int(gates[j])],
                    exposure=t,
                    route_auxiliary_ports=aux,
                )
                np.testing.assert_allclose(
                    [sum(m.trace().real for m in part) for part in zip(*blocks)],
                    [1.0, 0.0],
                    atol=2e-13,
                )
                trace = 0.125 * sum(
                    np.abs(np.linalg.eigvalsh(v)).sum() for _, v in blocks
                )
                occupancy = 0.125 * sum(np.abs(np.diag(v)).sum() for _, v in blocks)
                traces[k].append(
                    dict(
                        **ledgers,
                        trace_distance=float(trace),
                        occupation_tv=float(occupancy),
                    )
                )
                for rails, phase, effect in rail_effects(read_block):
                    records[k].append(
                        dict(
                            **ledgers,
                            ports=[int(visible[v]) for v in rails],
                            analyzer=phase,
                            effect=effect,
                        )
                    )
    result = {}
    for k, approach in enumerate(approaches):
        score = evaluate(records[approach], traces[approach], level, 35 + k)
        winner = records[approach][score["selected_setting"]]
        selected = next(
            p
            for p in pairs
            if [int(gates[p[0]]), int(gates[p[1]])] == winner["gateways"]
        )
        sites, reader = selected[-2:]
        wave = np.zeros(len(sites), complex)
        indices = np.searchsorted(sites, winner["ports"])
        if len(indices) == 1:
            wave[indices] = 1.0
        else:
            wave[indices] = [
                1 / np.sqrt(2),
                (1j if winner["analyzer"] == "circular" else 1.0) / np.sqrt(2),
            ]
        np.testing.assert_allclose(compile_wave(reader, wave), wave, atol=1e-12)
        # Deterministic candidate order and the per-setting signed counts replay every choice.
        score["selected_reader"] = winner
        score["calibration_records"] = [
            dict(effect=r["effect"]) for r in score["calibration_records"]
        ]
        result[approach] = score
    return dict(cells=board.cells, studies=result)


def run():
    return dict(
        protocol_sha256=hashlib.sha256(
            (HERE / "TWO_GATEWAYS.md").read_bytes()
        ).hexdigest(),
        seed=SEED,
        shots_per_ensemble=SHOTS,
        boards=[measure(i) for i in range(3)],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = run()
    if args.output:
        with args.output.open("x") as stream:
            json.dump(data, stream, separators=(",", ":"), allow_nan=False)
            stream.write("\n")
    else:
        from model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.coarse_response import (
            check_replay,
        )

        check_replay(json.loads((HERE / "two_gateways.json").read_text()), data)
    for board in data["boards"]:
        for name, s in board["studies"].items():
            print(
                board["cells"],
                name,
                max(r["trace_distance"] for r in s["settings"]),
                max(r["occupation_tv"] for r in s["settings"]),
                s["held_difference"],
                s["held_standard_error"],
                s["selected_reader"],
            )
