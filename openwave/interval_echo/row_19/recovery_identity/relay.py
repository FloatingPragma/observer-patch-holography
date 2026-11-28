"""Transmit source differences to a third seam, then release fresh probes there."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from model.experiments.interval_echo_einstein_station.shared_geometry_probe.flight import (
    prepared_states,
)
from . import tidal

RANGING = tidal.front.optics.ranging
FORMATION = RANGING.formation
HERE = tidal.HERE
RECEIPT = HERE / "relay.json"
NAMES = (
    *tidal.STOCKS,
    "vacancy_formed",
    "vacancy_product",
    "excess_formed",
    "excess_product",
)
TIMES = (2, 4, 8)
DT = 0.25
SHOTS = 2**60
SEED = 191130
# Two stations, three times, eight stocks: three spin settings and two
# charges, each with two spatial launches at three appointments.
RAD = np.sqrt(np.log(2 * 2 * 3 * 8 * (3 + 2 * 2 * 3) / 0.01) / (2 * SHOTS))


def receiver_paths(board, route, radii):
    cells = np.array(board.external) // 12
    graph = sparse.coo_matrix((np.ones(len(cells)), cells.T), shape=(board.cells,) * 2)
    depth = shortest_path(
        graph.tocsr(),
        directed=False,
        unweighted=True,
        indices=[route[1] // 12, route[2] // 12],
    ).min(axis=0)
    paths = []
    for level in (1, 2):
        paths.append(
            min(
                (radii[a], a, b, c)
                for a in range(board.n)
                if depth[a // 12] == level
                for b in board.neighbours[a]
                if a // 12 == b // 12 and radii[b] == radii[a] + 1
                for c in board.neighbours[b]
                if b // 12 == c // 12 and radii[c] == radii[b] + 1
            )[1:]
        )
    return paths


def generator(board, contacts, q, cut_owner=None):
    """Replace the blank blocks at the enrolled seams by their native written blocks."""
    edges = [tuple(sorted(edge)) for edge, _ in contacts]
    blanks = np.array([e for e in board.seams if e not in edges])
    degree = np.bincount(blanks.ravel(), minlength=board.n)
    hopping = blanks
    if cut_owner is not None:
        hopping = blanks[
            (blanks[:, 0] // 12 == cut_owner) == (blanks[:, 1] // 12 == cut_owner)
        ]
    a, b = hopping.T
    adjacency = sparse.csr_matrix(
        (np.full(2 * len(a), -0.5), (np.r_[a, b], np.r_[b, a])),
        shape=(board.n, board.n),
    )
    diagonal = 8.5 * degree + 8 * (2 * q + q * q) * (6 - board.degree)
    nhead = len(contacts)
    h = sparse.kron(
        sparse.diags(diagonal) + adjacency, sparse.eye(2**nhead), format="csr"
    )
    for slot, (edge, parity) in enumerate(contacts):
        parts = RANGING.parity_blocks(q, 8, parity)
        i, j = sorted(edge)
        for u, dst in enumerate((i, j)):
            for v, src in enumerate((i, j)):
                block = parts[u + 1, v + 1] - (u == v) * parts[0, 0]
                head = np.kron(
                    np.eye(2**slot), np.kron(block, np.eye(2 ** (nhead - slot - 1)))
                )
                dyad = sparse.csr_matrix(
                    ([1], ([dst], [src])), shape=(board.n, board.n)
                )
                h += sparse.kron(dyad, head, format="csr")
    return h - h.diagonal().mean() * sparse.eye(h.shape[0], format="csr")


def sources(board, route):
    stocks, steps, interval, info = tidal.transfer.inputs()
    states = [stocks[k] for k in tidal.STOCKS]
    for s in (-1, 1):
        formed, _, _ = prepared_states(board, route, s)
        formed = formed / np.trace(formed)
        for state in (formed, np.kron(*FORMATION.marginals(formed))):
            states.append(np.stack((state, np.zeros_like(state))))
    return np.array(states), steps, (float(interval.a), float(interval.b)), info


def deliver(board, route, edge, stocks, cut=False):
    # Columns specify the sender coin basis, all with receiver coin zero.
    start = np.zeros((8 * board.n, 4))
    start[8 * route[0] + 2 * np.arange(4)] = np.eye(4)
    received = np.zeros((len(TIMES), len(stocks), 2, 2), complex)
    for b in (0, 1):
        contacts = ((route[:2], 1 - 2 * b), (route[2:], 1), (edge, 1))
        h = generator(board, contacts, 1, edge[0] // 12 if cut else None)
        evolved = expm_multiply(-1j * h, start, start=0, stop=8, num=5, traceA=0)
        for it, t in enumerate(TIMES):
            w = evolved[t // 2]
            np.testing.assert_allclose(w.conj().T @ w, np.eye(4), atol=4e-12, rtol=0)
            w = w.reshape(-1, 2, 4)
            # Trace over the courier and both sending heads, including all missed visits.
            channel = np.einsum("pca,pdb->cdab", w, w.conj(), optimize=True)
            received[it] += np.einsum("cdab,sab->scd", channel, stocks[:, b])
    np.testing.assert_allclose(
        received.trace(axis1=-2, axis2=-1), 1, atol=5e-12, rtol=0
    )
    if np.linalg.eigvalsh(received).min() < -1e-12:
        raise ValueError("Receiver state is not positive")
    return received


def position_reader(board, path, radii, q):
    h = generator(board, [(path[:2], 1)], q)
    launch = np.zeros((board.n, 2))
    for j in range(2):
        launch[list(path[j : j + 2]), j] = 2**-0.5
    psi = np.kron(launch, np.eye(2))
    acc = []
    for w in np.split(psi, 2, axis=1):
        v = h @ w
        weights = np.repeat(radii, 2)[:, None]
        second = w.T @ (weights * (h @ v))
        acc.append(2 * v.T @ (weights * v) - second - second.T)
    waves = expm_multiply(-1j * h, psi, start=0, stop=2 * DT, num=3, traceA=0)
    effects = []
    for wave in waves:
        w = wave.reshape(board.n, 2, 2, 2)
        shells = [
            np.einsum("pnri,pnrj->rij", w[radii == k].conj(), w[radii == k])
            for k in range(radii.max() + 1)
        ]
        effects.append(shells)
    return acc[1] - acc[0], np.array(effects)


def forecast(values, steps, ratio, endpoint_radius):
    ordinary, held = (
        (values[1] - values[0]) / (2 * steps["ordinary"]),
        (values[3] - values[2]) / (2 * steps["held"]),
    )
    ro, rh = endpoint_radius / steps["ordinary"], endpoint_radius / steps["held"]
    forecast = sorted(a * b for a in ratio for b in (ordinary - ro, ordinary + ro))
    resolved = abs(ordinary) > max(1e-11, ro)
    return dict(
        ordinary=float(ordinary),
        held=float(held),
        ordinary_radius=float(ro),
        held_radius=float(rh),
        ratio=float(held / ordinary) if resolved else None,
        predicted=[float(forecast[0]), float(forecast[-1])],
        compatible=bool(held + rh >= forecast[0] and held - rh <= forecast[-1])
        if resolved
        else None,
    )


def position_reads(states, reader, rng, steps, ratio):
    accel, effects = reader
    laws = np.einsum("tmrij,sji->srtm", effects, states).real
    np.testing.assert_allclose(laws.sum(axis=-1), 1, atol=5e-12, rtol=0)
    if laws.min() < -1e-15:
        raise ValueError("Negative shell probability")
    counts = np.empty(laws.shape, np.int64)
    for index in np.ndindex(laws.shape[:-1]):
        p = np.maximum(laws[index], 0)
        counts[index] = rng.multinomial(SHOTS, p / p.sum())
    if not np.all(counts.sum(axis=-1) == SHOTS):
        raise ValueError("Incomplete shell bank")
    shell = np.arange(laws.shape[-1])

    def difference(means):
        return np.diff(means, axis=1)[:, 0] @ np.array([1, -2, 1]) / DT**2

    exact = difference(laws @ shell)
    measured = difference((counts / SHOTS) @ shell)
    release = np.einsum("ij,sji->s", accel, states).real
    radius = 8 * shell[-1] * RAD / DT**2
    return dict(
        counts=counts.tolist(),
        acceleration_operator=accel.tolist(),
        release_prediction=forecast(release, steps, ratio, 0),
        finite_prediction=forecast(exact, steps, ratio, 0),
        sampled_prediction=forecast(measured, steps, ratio, radius),
        formed_release_contrasts=(release[[5, 7]] - release[[4, 6]]).tolist(),
        formed_population_contrasts=(exact[[5, 7]] - exact[[4, 6]]).tolist(),
        formed_sampled_contrasts=(measured[[5, 7]] - measured[[4, 6]]).tolist(),
        formed_sampled_radius=float(2 * radius),
    )


def run():
    board = FORMATION.federation(1)
    route = FORMATION.region(board)
    radii, _, _ = tidal.chart(board, route)
    stocks, steps, ratio, src = sources(board, route)
    rng = np.random.default_rng(SEED)
    stations = []
    for path in receiver_paths(board, route, radii):
        edge = path[:2]
        states = deliver(board, route, edge, stocks)
        cut = deliver(board, route, edge, stocks, cut=True)[-1]
        cut_error = float(np.max(abs(cut - np.diag([1, 0]))))
        np.testing.assert_allclose(
            cut, np.broadcast_to(np.diag([1, 0]), cut.shape), atol=5e-12, rtol=0
        )
        readers = {q: position_reader(board, path, radii, q) for q in (-1, 1)}
        appointments = {}
        for t, received in zip(TIMES, states, strict=True):
            bloch = np.einsum("kij,sji->sk", FORMATION.PAULI, received).real
            probabilities = (1 + bloch) / 2
            if probabilities.min() < -3e-12 or probabilities.max() > 1 + 3e-12:
                raise ValueError("Pauli probabilities outside support")
            plus = rng.binomial(SHOTS, probabilities.clip(0, 1))
            measured = 2 * (plus / SHOTS) - 1
            appointments[str(t)] = dict(
                bloch=bloch.tolist(),
                plus=plus.tolist(),
                formed_measured_pauli=(measured[[5, 7]] - measured[[4, 6]]).tolist(),
                probes={
                    str(q): position_reads(received, reader, rng, steps, ratio)
                    for q, reader in readers.items()
                },
            )
        stations.append(
            dict(
                path=list(path),
                radii=radii[list(path)].tolist(),
                cut_error=cut_error,
                appointments=appointments,
            )
        )
    return dict(
        protocol_sha256=sha256((HERE / "relay_protocol.md").read_bytes()).hexdigest(),
        source=src,
        stocks=list(NAMES),
        half_steps=steps,
        shots_per_setting=SHOTS,
        seed=SEED,
        pauli_difference_radius=float(4 * RAD),
        stations=stations,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with RECEIPT.open("x") as out:
            json.dump(result, out, separators=(",", ":"), allow_nan=False)
            out.write("\n")
    else:
        tidal.verify(result, json.loads(RECEIPT.read_text()))
    for station in result["stations"]:
        print(
            "path",
            station["path"],
            "radii",
            station["radii"],
            "cut",
            station["cut_error"],
        )
        for time, row in station["appointments"].items():
            print(
                time,
                "received",
                np.array(row["formed_measured_pauli"])[:, 2],
                "scalar",
                [
                    (
                        p["release_prediction"]["ratio"]
                        if abs(p["release_prediction"]["ordinary"])
                        * 2
                        * result["half_steps"]["ordinary"]
                        > 1e-11
                        else None,
                        p["sampled_prediction"]["compatible"],
                    )
                    for p in row["probes"].values()
                ],
            )
