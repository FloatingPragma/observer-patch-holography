"""Three resolved orbitals transferring a rail qubit between initially disjoint regions."""

import argparse
import json
from scipy.optimize import minimize_scalar

import networkx as nx
import numpy as np

from model.executable_model import federation
from model.experiments.interval_echo_angular.charged_orientation.paired_flight import (
    material as rail,
)
from model.experiments.interval_echo_angular.charged_orientation.paired_flight.measure import (
    logical_sources,
    read,
    spatial_read,
)
from .material import HERE, programmed_mark, synthesize
from ..conjugate_flight.measure import binomial_interval

N = 13
LEFT, RIGHT = np.arange(2, 5), np.arange(8, 11)
SITES = [route[j] for route in rail.RAILS for j in [*LEFT, *RIGHT]]
SHOTS = 2**24


def side_contacts():
    board = federation(0)
    both = set(sum(rail.RAILS, ()))
    graph = nx.Graph()
    for site in SITES:
        for q in sorted(board.neighbours[site]):
            if q not in both and q // 12 == site // 12 and board.degree[q] == 5:
                graph.add_edge(site, int(q))
    pairing = nx.bipartite.maximum_matching(graph, top_nodes=sorted(SITES))
    assert all(s in pairing for s in SITES)
    return sorted((p, pairing[p]) for p in SITES)


def store():
    # Confine the balancing charge to a separate dimer, away from every side port.
    b = federation(0)
    graph = nx.Graph(b.seams)
    dist = nx.multi_source_dijkstra_path_length(graph, set(sum(rail.RAILS, ())))
    a, c = next((a, c) for a, c in b.seams if min(dist[a], dist[c]) > 3)
    cuts = [e for e in b.seams if (a in e or c in e) and e != (a, c)]
    # Park the charge alternately at the two endpoints for each equal-load write.
    assert all((c if a in e else a) not in e for e in cuts)
    graph.remove_edges_from(cuts)
    assert set(nx.node_connected_component(graph, a)) == {a, c}
    return dict(dimer=[a, c], cut_count=len(cuts))


def beat(times, hits, shots):
    # Search the entire positive Nyquist interval using port counts alone.
    y = np.asarray(hits) / shots
    x = times / times[-1]

    def error(freq):
        f = 2 * np.pi * freq * x
        design = np.array([np.ones_like(x), np.cos(f), np.sin(f)]).T
        fit = np.linalg.lstsq(design, y, rcond=None)[0]
        return np.sum((design @ fit - y) ** 2)

    grid = np.arange(0.125, (len(times) - 1) / 2 + 0.001, 0.125)
    best = grid[np.argmin([error(f) for f in grid])]
    fit = minimize_scalar(
        error, bounds=(max(1e-4, best - 0.125), best + 0.125), method="bounded"
    )
    return 2 * np.pi * fit.x / times[-1]


def orbital_analyzer(columns, sign):
    """Compile the complete three-outcome orbital decoder with native adjacent contacts."""
    mat = columns.copy()
    gates = np.eye(3, dtype=complex)
    generator = rail.pair(8, sign)
    from scipy.linalg import expm

    def phase(site, angle):
        hold = (-angle % (2 * np.pi)) / 32
        gates[site] *= np.exp(-32j * hold)

    for j in range(2):
        for b in range(2, j, -1):
            a = b - 1
            theta = np.arctan2(mat[b, j], mat[a, j])
            c, s = np.cos(theta), np.sin(theta)
            mat[[a, b]] = np.array([[c, s], [-s, c]]) @ mat[[a, b]]
            t = (2 * theta) % (4 * np.pi)
            phase(b, -np.pi / 2)
            gates[[a, b]] = expm(-1j * t * generator) @ gates[[a, b]]
            phase(a, 8.5 * t)
            phase(b, 8.5 * t + np.pi / 2)
    for j in range(3):
        phase(j, -np.angle(mat[j, j]))
    np.testing.assert_allclose(gates @ columns, np.eye(3), atol=1e-12)
    return gates


def guide(sign, changed=False):
    h = rail.native_matrix(8, sign)
    for i, path in enumerate(rail.RAILS):
        for source, _ in side_contacts():
            if source in path and not (changed and i == 1 and source == path[9]):
                j = i * N + path.index(source)
                h[j, j] += programmed_mark(8)
    return h


def material_energy(v, h, start, rng):
    # An enrolled left orbital occupies three internal ports with two live
    # contacts apiece and one marked side head. Sample a uniform blank bond.
    # The port and its two silent neighbours read the decoded position.
    p = np.zeros(N, complex)
    p[start] = v
    site = rng.multinomial(SHOTS, abs(p) ** 2 / np.vdot(p, p).real)
    # Native exchange antisymmetry is the only term with moving load on the blank path.
    odd = np.sum(abs(p[:-1] - p[1:]) ** 2) / (2 * (N - 1))
    hits = int(rng.binomial(SHOTS, odd))
    ci = 48 + (N - 1) * binomial_interval(hits, SHOTS, 0.01 / 12)
    assert abs(np.vdot(p, h @ p).real - (48 + (N - 1) * odd)) < 2e-12
    return dict(
        ports=site.tolist(),
        odd=hits,
        measured=48 + (N - 1) * hits / SHOTS,
        interval=ci.tolist(),
    )


def run():
    rng = np.random.default_rng(601143)
    times = np.arange(257) * 32768
    cases = []
    for sign in (-1, 1):
        h = guide(sign)
        path = h[:N, :N]
        np.testing.assert_allclose(h, np.kron(np.eye(2), path), atol=1e-13)
        perturbed = guide(sign, changed=True)[N:, N:]
        eig, modes = np.linalg.eigh(path)
        alt, other = np.linalg.eigh(perturbed)
        local, u = np.linalg.eigh(path[np.ix_(LEFT, LEFT)])
        analyzer = orbital_analyzer(u[::-1, :], sign)
        for b in range(3):
            v = synthesize(u[:, b].astype(complex))
            src = np.zeros(N, complex)
            dst = src.copy()
            src[LEFT] = v
            dst[RIGHT] = analyzer[b].conjugate()
            weight = abs(modes.T @ src) ** 2
            chosen = abs(eig - local[b]) < 0.2
            assert sum(chosen) == 2
            coeff = (dst.conj() @ modes) * (modes.T @ src)
            amplitudes = np.exp(-1j * np.outer(times, eig)) @ coeff
            pilots = rng.binomial(2**14, np.clip(abs(amplitudes) ** 2, 0, 1))
            j = np.argmax(pilots)
            t = int(times[j])
            arrival = amplitudes[j] * np.eye(2)
            altered = (dst.conj() @ other) @ (np.exp(-1j * t * alt) * (other.T @ src))
            sources = logical_sources(8, sign)
            channel = read(arrival, sources, rng)
            # Blind three-outcome orbital read at the destination, incl. every miss.
            field = modes @ (np.exp(-1j * t * eig) * (modes.T @ src))
            orbital_prob = abs(analyzer @ field[RIGHT]) ** 2
            all_bins = rng.multinomial(
                SHOTS, [*orbital_prob, max(0, 1 - sum(orbital_prob))]
            )
            cases.append(
                dict(
                    sign=sign,
                    band=b,
                    launch_measurement=material_energy(v, path, LEFT, rng),
                    two_lines=eig[chosen].tolist(),
                    doublet_weight=float(sum(weight[chosen])),
                    initial_receiver_probability=float(abs(np.vdot(dst, src)) ** 2),
                    time=t,
                    pilot_counts=pilots.tolist(),
                    beat_from_pilot_counts=beat(times, pilots, 2**14),
                    remote_orbitals_and_miss=all_bins.tolist(),
                    channel=channel,
                    spatial_orders=spatial_read(arrival, sources, rng),
                    unmatched_head=read(
                        np.diag([amplitudes[j], altered]), sources, rng
                    ),
                )
            )
    return dict(
        shots_energy_and_orbital=SHOTS,
        shots_per_process_setting=2**22,
        side_contacts=side_contacts(),
        remote_store=store(),
        left=LEFT.tolist(),
        right=RIGHT.tolist(),
        pilot_times=times.tolist(),
        cases=cases,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    data = json.dumps(run(), indent=2) + "\n"
    if p.parse_args().write:
        (HERE / "paired.json").write_text(data)
    else:
        print(data, end="")
