"""Release one coherent probe on every gateway of each full radial cut."""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json

import networkx as nx
import numpy as np

from . import routed_contacts as routes

local = routes.local
HERE = local.HERE


def output_error(value, bound):
    error = bound + abs(F(float(value)) - value)
    return float(np.nextafter(float(error), np.inf))


def full_sector(board, adjacency, inside, cut, deposit, m, q):
    # Apply twice the native Hamiltonian to psi and X psi. Only the seams
    # touching a launched load can contribute to either vector.
    hpsi, hxpsi = defaultdict(int), defaultdict(int)
    for g in (x for edge in cut for x in edge):

        def canonical(k):
            return tuple(sorted(k)) if m == q else tuple(k)

        key = canonical([deposit, g] if m else [g])
        loads = dict(zip(key, (m, q) if m else (q,), strict=True))
        weight = sum(v // 12 in inside for v in key)
        diag = 16 * sum(
            (6 - int(board.degree[v])) * load * load for v, load in loads.items()
        )
        edges = {tuple(sorted((v, w))) for v in key for w in adjacency[v]}
        for v, w in edges:
            delta = loads.get(v, 0) - loads.get(w, 0)
            if not delta:
                continue
            diag += 1 + 16 * delta**2
            other = canonical(w if x == v else v if x == w else x for x in key)
            hpsi[other] -= 1
            hxpsi[other] -= weight
        hpsi[key] += diag
        hxpsi[key] += weight * diag
    total = sum(
        z * (sum(v // 12 in inside for v in k) * z - hxpsi[k]) for k, z in hpsi.items()
    )
    return F(total, 4 * len(cut))


def prepare_shell(tree, cut, q):
    # Split a real wave over calibrated blank contacts. Both output phase
    # settings use the same local reference as the clock recombination.
    root = 0
    targets = {v for edge in cut for v in edge}
    order = list(nx.bfs_tree(tree, root))
    children = dict(nx.bfs_successors(tree, root))
    weights = {}
    for v in reversed(order):
        weights[v] = int(v in targets) + sum(weights[w] for w in children.get(v, ()))
    assert weights[root] == len(targets)
    wave = np.zeros(len(order), complex)
    wave[root] = 1
    pulses = 0
    for v in order:
        remaining = weights[v]
        for w in children.get(v, ()):
            n = weights[w]
            if not n:
                continue
            sine = np.sqrt(n / remaining)
            cosine = np.sqrt((remaining - n) / remaining)
            duration = 2 * np.arctan2(sine, cosine)
            phase = np.exp(-1j * duration * (8 * q * q + 0.5))
            # Native single-load half-Laplacian block, followed by the two
            # phase settings that align transmitted and retained components.
            wave[v], wave[w] = phase * cosine * wave[v], phase * 1j * sine * wave[v]
            wave[v] *= phase.conjugate()
            wave[w] *= -1j * phase.conjugate()
            remaining -= n
            pulses += 1
    ideal = np.zeros(len(order))
    ideal[list(targets)] = 1 / np.sqrt(len(targets))
    err = float(np.linalg.norm(ideal - wave))
    assert err < 2e-13
    return {"partial_swaps": pulses, "amplitude_error": err}


def run():
    profiles = local.inputs("routed_contacts")["profiles"]
    base = set(local.survey.routes.codebook()[0][0])
    usable = set(local.survey.routes.cell().seams) - base
    out = []
    verified = 0
    for level in (0, 1, 2):
        board = local.survey.bulk.federation(level)
        steps, _, _ = routes.layout(board)
        marks = routes.receiver_codes(board, steps)
        gates = [list(np.flatnonzero(z != v)) for v, z in enumerate(steps)]

        def partner(v, g):
            return next(x + y - g for x, y in marks[v] if g in (x, y))

        where = [[partner(v, g) for g in ports] for v, ports in enumerate(gates)]
        graph = nx.Graph((a // 12, b // 12) for a, b in board.external)
        hops = nx.single_source_shortest_path_length(graph, 0)
        blank = defaultdict(set)
        sharp = defaultdict(set)
        flags = np.random.default_rng(99113).integers(2, size=board.cells)
        for a, b in board.seams:
            c = a // 12
            e = (a % 12, b % 12)
            if c == b // 12 and e in base:
                continue
            blank[a].add(b)
            blank[b].add(a)
            if c != b // 12 or not flags[c] or e not in marks[c]:
                sharp[a].add(b)
                sharp[b].add(a)
        tree = nx.Graph()
        tree.add_nodes_from(range(board.n))
        tree.add_edges_from(
            (a, b)
            for a in range(board.n)
            for b in sorted(blank[a])
            if a < b and (a // 12 != b // 12 or (a % 12, b % 12) not in marks[a // 12])
        )
        assert nx.is_connected(tree)
        for k in range(max(hops.values())):
            inside = {c for c in hops if hops[c] <= k}
            cut = [
                (a, b) if a // 12 in inside else (b, a)
                for a, b in board.external
                if (a // 12 in inside) != (b // 12 in inside)
            ]

            launch = [prepare_shell(tree, cut, q) for q in (1, 2)]

            def mean_drop(f):
                return sum(f(a) - f(b) for a, b in cut) / len(cut)

            if level == 0:
                # Every terminal cell, sharp mark patterns, and both ordinary
                # loads. Each check applies the full native two-load sector.
                for q in (1, 2):
                    baseline = full_sector(board, blank, inside, cut, 0, 0, q)
                    write = full_sector(board, sharp, inside, cut, 0, 0, q)
                    assert write - baseline == (F(1, 4) + 4 * q * q) * mean_drop(
                        lambda v: F(int(flags[v // 12]))
                    )
                    for c in range(board.cells):
                        s = int(12 * c + where[c][(k + c) % 3])
                        for m in (1, 2):
                            measured = full_sector(board, sharp, inside, cut, s, m, q)
                            weight = F(1 + (m == q), 4) + 8 * m * q
                            assert measured - write == weight * mean_drop(
                                lambda v: F(int(s in sharp[v]))
                            )
                            verified += 1
            for row in (r for r in profiles if r["stations"] == board.cells):
                D = list(map(F, row["diagonal_green"]))
                g = list(map(F, row["green_from_source"]))
                u = [x / (2 * y - 1) for x, y in zip(g, D, strict=True)]
                stop = F(1, 2 ** row["stop_bits"])
                pf = F(25, 36)

                def available(v):
                    c, j = divmod(v, 12)
                    count = sum(
                        tuple(sorted((int(j), int(p)))) in usable for p in where[c]
                    )
                    return stop * (count * g[c] - pf * D[c] * u[c]) / 3

                record = mean_drop(lambda v: u[v // 12]) / 144
                ordinary = mean_drop(available)
                eps = F(float(np.nextafter(row["green_rounding_bound"], np.inf)))
                uerr = 3 * eps / (2 * min(D) - 1)
                recorderr = 2 * uerr / 144
                ordinaryerr = 2 * stop * (eps + pf * (eps + max(D) * uerr) / 3)
                out.append(
                    {
                        "stations": board.cells,
                        "stop_bits": row["stop_bits"],
                        "shell": k,
                        "crossings": len(cut),
                        "launches": launch,
                        "record_gradient": float(record),
                        "ordinary_gradient": float(ordinary),
                        "record_error": output_error(record, recorderr),
                        "ordinary_error": output_error(ordinary, ordinaryerr),
                        "record_curvatures": [
                            float(17 * record / 4),
                            float(65 * record / 4),
                        ],
                        "ordinary_one_curvatures": [
                            float(17 * ordinary / 2),
                            float(65 * ordinary / 4),
                        ],
                        "ordinary_two_curvatures": [
                            float(65 * ordinary / 4),
                            float(65 * ordinary / 2),
                        ],
                    }
                )
    return {
        "native_two_load_comparisons": verified,
        "shells": out,
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in (
                "SHELL_RECEIVERS.md",
                "routed_contacts.json",
                "routed_contact_flight.json",
            )
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    z = run()
    path = HERE / "shell_receivers.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(z, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), z)
    rows = z["shells"]
    print(
        len(rows),
        z["native_two_load_comparisons"],
        min(r["record_gradient"] - r["record_error"] for r in rows),
        sum(r["ordinary_gradient"] > r["ordinary_error"] for r in rows),
    )
