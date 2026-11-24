"""Static source reversal in a native quantum polarization guide."""

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm, svdvals

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import federation

HERE = Path(__file__).resolve().parent
SHOTS = 2**23
TIMES = 0.05 * 2 ** np.arange(12)
SEED = 201124
BIN = (-2, -1, 0, 1, 2)


def guide(board):
    owners = nx.Graph()
    port_graph = nx.Graph(board.seams)
    gates = {}
    for a, b in board.external:
        i, j = a // 12, b // 12
        owners.add_edge(i, j)
        gates[i, j], gates[j, i] = a, b
    rings = nx.single_source_shortest_path_length(owners, 0)
    end = min(i for i, d in rings.items() if d == max(rings.values()))
    cells = nx.shortest_path(owners, 0, end)
    src = gates[cells[0], cells[1]]
    path = [min(set(board.neighbours[src]) & set(range(12))), src]
    for i in range(1, len(cells)):
        cell = cells[i]
        a = gates[cell, cells[i - 1]]
        path.append(a)
        if i == len(cells) - 1:
            path.append(min(p for p in board.neighbours[a] if p // 12 == cell))
            break
        b = gates[cell, cells[i + 1]]
        inside = [
            p
            for p in range(12 * cell, 12 * cell + 12)
            if board.degree[p] == 5 or p in (a, b)
        ]
        piece = min(
            (
                x
                for x in nx.all_simple_paths(
                    port_graph.subgraph(inside), a, b, cutoff=7
                )
                if len(x) % 2 == 0
            ),
            key=lambda x: (len(x), x),
        )
        path.extend(piece[1:])
    assert len(set(path)) == len(path)
    assert all(
        board.degree[a] == board.degree[b]
        for a, b in zip(path[1:-1:2], path[2:-1:2], strict=True)
    )
    distances = nx.single_source_shortest_path_length(port_graph, path[0])
    return cells, path, distances


def ranging(board, path, rng):
    """Route one amount-two token, keeping shutters and clock independent of spectroscopy."""
    exchange, contrast = law.local_term_parts()
    gate = expm(-1j * np.pi * (exchange + candidate.KAPPA * contrast))
    indices = [law.SRC.index(t) for t in ((2, 0), (0, 2))]
    pair = gate[np.ix_(indices, indices)]
    graph = nx.Graph(board.seams)
    result = []
    for n in range(2, (len(path) - 2) // 2 + 1):
        dest = path[2 * n + 1]
        route = nx.shortest_path(graph, path[0], dest)
        arms = []
        for omit in (None, len(route) // 2 - 1):
            psi = np.zeros(board.n, complex)
            psi[route[0]] = 1
            arrivals = []
            for k, (a, b) in enumerate(zip(route, route[1:])):
                if k != omit:
                    psi[[a, b]] = pair @ psi[[a, b]]
                arrivals.append(
                    int(rng.binomial(4096, round(float(abs(psi[dest]) ** 2), 12)))
                )
            arms.append(arrivals)
        result.append(dict(n=n, route=route, shots_per_wait=4096, open_and_held=arms))
    return result


def couplings(kappa):
    exchange, contrast = law.local_term_parts()
    local = exchange + kappa * contrast
    occupied, empty, rev = (law.SRC.index(t) for t in ((2, 0), (0, 0), (0, 2)))
    return float((local[occupied, occupied] - local[empty, empty]) / 2), float(
        -local[occupied, rev]
    )


def splitting(n, J, h):
    # The odd wall-parity sector starts with one fermion on the open dual chain.
    b = J * np.eye(n + 1) - h * np.eye(n + 1, k=-1)
    return float(2 * svdvals(b)[-1])


def phase_reader(energy, rng, times=TIMES):
    phase = times * energy
    p = (1 + np.array([np.cos(phase), np.sin(phase)])) / 2
    draws = rng.binomial(SHOTS, np.round(p.clip(0, 1), 13))
    z = 2 * draws / SHOTS - 1
    wrapped = np.arctan2(z[1], z[0])
    angle = wrapped[0]
    for current in wrapped[1:]:
        angle = current + 2 * np.pi * round((2 * angle - current) / (2 * np.pi))
    # A union bound over <= 30,000 binomial quadratures, transformed to phase.
    eps = np.sqrt(2 * np.log(2 * 30000 / 0.01) / SHOTS)
    angular_error = np.arcsin(np.sqrt(2) * eps)
    assert 3 * angular_error < np.pi
    return dict(
        counts_xy=draws.tolist(),
        energy=float(angle / times[-1]),
        simultaneous_error=float(angular_error / times[-1]),
    )


def sign_read(q):
    col = law.echo()[:, law.SRC.index((q, 1))].reshape(-1, law.N_SRC)
    weight = (col[1:] ** 2).sum(axis=1)
    return [
        int(
            round(
                SHOTS
                * sum(p for p, (_, delta) in zip(weight, law.LABELS) if delta == b)
            )
        )
        for b in BIN
    ]


def trial(n, kappa, rng):
    j, h = couplings(kappa)
    value = splitting(n, j, h)
    records = [phase_reader(s * value, rng) for s in (-1, 0, 1)]
    # Combining independent reversed sources removes common detector phase shifts.
    read = (records[2]["energy"] - records[0]["energy"]) / 2
    return dict(
        n=n,
        contrast=kappa,
        ordering=j,
        tunnelling=h,
        asymptote=2 * max(0, j - h),
        exact=value,
        source_minus_zero_plus=records,
        measured=read,
        simultaneous_error=records[0]["simultaneous_error"],
    )


def forecasts(rows, distance):
    training = [r for r in rows if r["n"] in (4, 6)]
    held = [r for r in rows if r["n"] >= 8]
    if not held or len(training) != 2:
        return None
    x, y = (
        np.array([1 / r[distance] for r in training]),
        np.array([r["measured"] for r in training]),
    )
    if x[0] == x[1]:
        return None
    errors = []
    for r in held:
        w = (1 / r[distance] - x[0]) / (x[1] - x[0])
        pred = (1 - w) * y[0] + w * y[1]
        err = abs(r["measured"] - pred)
        bound = (1 + abs(1 - w) + abs(w)) * r["simultaneous_error"]
        errors.append(
            dict(
                n=r["n"],
                predicted=pred,
                observed=r["measured"],
                error=err,
                bound=bound,
                accepted=bool(err < bound),
            )
        )
    return errors


def run():
    rng = np.random.default_rng(SEED)
    ruler_rng = np.random.default_rng(201129)
    boards = []
    for level in range(3):
        board = federation(level)
        owners, path, dist = guide(board)
        cases = []
        for k in (candidate.KAPPA, 0, 1 / 8):
            rows = []
            for n in range(2, (len(path) - 2) // 2 + 1, 2):
                row = trial(n, k, rng)
                row.update(guide_links=2 * n + 1, native_distance=dist[path[2 * n + 1]])
                rows.append(row)
            cases.append(
                dict(
                    contrast=k,
                    rows=rows,
                    guide_held=forecasts(rows, "guide_links"),
                    intact_held=forecasts(rows, "native_distance"),
                )
            )
        boards.append(
            dict(
                cells=board.cells,
                owners=owners,
                ports=path,
                studies=cases,
                shortest_route_reads=ranging(board, path, ruler_rng),
            )
        )
    return dict(
        seed=SEED,
        shots_per_quadrature=SHOTS,
        exposures=TIMES.tolist(),
        sign_bins=list(BIN),
        local_meter=[sign_read(s) for s in (0, 1, 2)],
        empty_medium=phase_reader(0, rng),
        boards=boards,
    )


def portable(d):
    if isinstance(d, list):
        return [portable(x) for x in d]
    if isinstance(d, dict):
        return {k: portable(v) for k, v in d.items()}
    if isinstance(d, (float, np.floating)):
        return round(float(d), 11)
    return d


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = portable(run())
    file = HERE / "records.json"
    if args.write:
        file.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    elif result != json.loads(file.read_text()):
        raise ValueError("Polarization guide records differ from the receipt")
    print(
        "Verified three federation routes, end-charge spectra and signed Ramsey reads."
    )


if __name__ == "__main__":
    main()
