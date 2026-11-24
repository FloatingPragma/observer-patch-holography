"""Three measured winding walls in a commensurate material reference."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import candidate, currents
from model.candidate_laws.contrast_memory_spherical_v1 import law
from ..written_junction import ring as written_ring
from . import belt, composite


def core_read(sites, length, ones, separate_sites=False):
    arc = composite.gap_arc(sites, length, ones, 3 - ones)
    if arc is None:
        return np.nan
    # Read local deviations of three-site sums from the single-load bulk.
    reads = []
    for a, b in zip(sites, (*sites[1:], sites[0])):
        gap = (b - a) % length
        if gap == 1:
            reads.append(((a, b), 2))
        if gap == 4:
            reads.append((((a + 2) % length,), -1))
    ports = [j for region, _ in reads for j in region]
    if len(ports) != len(set(ports)):
        return np.nan
    if separate_sites:
        footprints = [
            {(j + d) % length for j in region for d in (-1, 0, 1)}
            for region, _ in reads
        ]
        if sum(map(len, footprints)) != len(set.union(*footprints)):
            return np.nan
    occ = np.isin(np.arange(length), sites).astype(int)
    excess = np.roll(occ, 1) + occ + np.roll(occ, -1) - 1
    for region, target in reads:
        assert sum(excess[list(region)]) == target
    assert not excess[~np.isin(np.arange(length), ports)].any()
    return arc


def basis(length, ones):
    n = length // 3 + ones - 1
    sites = list(combinations(range(length), n))
    words = list(map(composite.mask, sites))
    ix = {w: i for i, w in enumerate(words)}
    moves = composite.hop_list(belt.ring(length).seams)
    rows = []
    cols = []
    for i, w in enumerate(words):
        dest = [ix[z] for z in composite.neighbours(w, moves)]
        rows.extend(dest)
        cols.extend([i] * len(dest))
    hop = sparse.csr_matrix(
        (np.full(len(rows), -0.5), (rows, cols)), shape=(len(ix), len(ix))
    )
    orbits = []
    seen = set()
    for w in words:
        if w in seen:
            continue
        orb = []
        while w not in seen:
            seen.add(w)
            orb.append(ix[w])
            w = ((w << 3) & ((1 << length) - 1)) | (w >> (length - 3))
        orbits.append(orb)
    arcs = np.array([core_read(s, length, ones) for s in sites])
    return sites, words, moves, hop, orbits, arcs


def momentum_basis(orbits, m, period, size):
    rows = []
    cols = []
    values = []
    reps = []
    for orb in orbits:
        if len(orb) * m % period:
            continue
        rows.extend(orb)
        cols.extend([len(reps)] * len(orb))
        reps.append(orb[0])
        values.extend(
            np.exp(2j * np.pi * m * np.arange(len(orb)) / period) / np.sqrt(len(orb))
        )
    return sparse.csr_matrix((values, (rows, cols)), shape=(size, len(reps))), np.array(
        reps
    )


def select(states, energies, tight, valid, best, m, limits):
    restricted = states[tight]
    probabilities = np.sum(abs(restricted) ** 2, axis=0)
    # The projector commutes with the three-site translation. Coherence between
    # momenta cannot improve the bound; within a repeated line it can matter.
    groups = np.split(
        np.arange(len(energies)), np.flatnonzero(np.diff(energies) > 1e-8) + 1
    )
    rotated = {}
    for g in groups:
        if len(g) > 1:
            r = restricted[:, g]
            vals, vec = np.linalg.eigh(r.conj().T @ r)
            probabilities[g] = -1
            probabilities[g[0]] = vals[-1]
            rotated[int(g[0])] = states[:, g] @ vec[:, -1]
    for tag, ok in [
        ("all", np.ones(len(energies), bool)),
        ("below", energies < limits[0] - 1e-7),
        ("above", energies > limits[1] + 1e-7),
    ]:
        ids = np.flatnonzero(ok)
        if not len(ids):
            continue
        i = ids[np.argmax(probabilities[ids])]
        if tag in best and probabilities[i] <= best[tag]["tight_probability"] + 1e-9:
            continue
        q = rotated.get(int(i), states[:, i])
        best[tag] = dict(
            energy=energies[i],
            momentum=m,
            tight_probability=probabilities[i],
            any_disjoint_wall_probability=float(np.sum(abs(q[valid]) ** 2)),
        )


def panel(length, ones, arrays, k, duty):
    sites, words, moves, hop, orbits, arcs = arrays
    # Four diagonal contact readings calibrate the actual reference coupling.
    a = np.diag(law.seam_term(k))

    def at(i, j):
        return a[law.SRC.index((i, j))]

    rung = at(1, 1) - at(0, 1) - at(1, 0) + at(0, 0)
    assert rung == -(2 * k + 1)
    potential = np.array(
        [
            composite.diagonal(w, moves, k) + duty * rung * sum(i % 3 != 0 for i in s)
            for w, s in zip(words, sites, strict=True)
        ]
    )
    h = hop + sparse.diags(potential)
    valid = np.isfinite(arcs)
    lo, hi = np.nanmin(arcs), np.nanmax(arcs)
    tight = arcs == lo
    far = arcs >= (lo + hi) / 2
    stretched = h[far][:, far].toarray()
    levels = np.linalg.eigvalsh(stretched)
    limits = [float(levels[0]), float(levels[-1])]
    primary_sites = np.array([core_read(s, length, ones, True) for s in sites])
    independent = np.isfinite(primary_sites)
    independent_best = {}
    best = {}
    ranks = 0
    period = length // 3
    for m in range(period // 2 + 1):
        tr, reps = momentum_basis(orbits, m, period, len(words))
        block = (tr.conj().T @ (h @ tr)).toarray()
        if np.abs(block.imag).max() < 1e-12:
            block = block.real
        e, v = eigh(block, driver="evd", check_finite=False)
        select(v, e, tight[reps], valid[reps], best, m, limits)
        if independent.any():
            select(
                v, e, independent[reps], independent[reps], independent_best, m, limits
            )
        ranks += len(e) * (1 if 2 * m % period == 0 else 2)
    assert ranks == len(words)
    return dict(
        length=length,
        material=len(sites[0]),
        positive_cores=ones,
        sharp_dimension=ranks,
        kappa=k,
        rung_duty=duty,
        measured_connected_rung=rung,
        arc_range=[lo, hi],
        disjoint_three_wall_words=int(sum(valid)),
        tight_words=int(sum(tight)),
        stretched_words=int(sum(far)),
        stretched_interval=limits,
        nonwall_sharp_energy_witnesses=[
            float(min(potential[~valid])),
            float(max(potential[~valid])),
        ],
        optimal_stationary_tight=best,
        disjoint_input_sites=dict(
            words=int(sum(independent)),
            tight_words=int(sum(independent & tight)),
            minimum_arc=float(np.min(primary_sites[independent]))
            if independent.any()
            else None,
            max_any_separation=independent_best["all"]["tight_probability"]
            if independent_best
            else 0,
        ),
    )


def tension(length, k, rng):
    # Common one-port terms cancel between fixed-count preparations.
    # The three-outcome seam analyzer reads the shifted native contact.
    ref = np.arange(length) % 3 != 0
    sites = np.arange(0, length, 3)
    shots = 2**18

    def measure(loc):
        n = np.isin(np.arange(length), loc)
        both = int(sum(n & ref))
        unequal = int(sum(n != ref))
        odd = int(rng.binomial(shots * unequal, 0.5))
        return dict(
            both_odd_other=[shots * both, odd, shots * (length - both) - odd],
            mean=unequal / 2 - 2 * k * both,
            measured=odd / shots - 2 * k * both,
        )

    baseline = measure(sites)
    rows = []
    for shift in (1, -2):
        for size in range(len(sites) + 1):
            moved = sites.copy()
            moved[:size] = (moved[:size] + shift) % length
            gaps = (np.roll(np.sort(moved), -1) - np.sort(moved)) % length
            r = measure(moved)
            rows.append(
                dict(
                    shift=shift,
                    shifted_loads=size,
                    anomalous_gaps=int(sum(gaps != 3)),
                    predicted=-(2 * k + 1) * size,
                    energy_change=r["mean"] - baseline["mean"],
                    measured_change=r["measured"] - baseline["measured"],
                    record=r["both_odd_other"],
                )
            )
    return dict(
        length=length, coupling=k, shots_per_rung=shots, baseline=baseline, rows=rows
    )


def run():
    board, cycles = written_ring.physical_cycles()
    g = nx.Graph(board.seams)
    out = []
    layouts = {}
    rng = np.random.default_rng(221213)
    tension_reads = []
    for length in (12, 15, 18):
        path = cycles[length]
        routes = [
            min(nx.all_shortest_paths(g, path[i], path[(i + 2) % length]))
            for i in range(length)
        ]
        # One distinct held reference port for every rung. Route a remote rung
        # along actual seams and undo the intervening swaps if needed.
        idle = sorted(set(range(board.n)) - set(path))
        refs = []
        rung_routes = []
        for j in path:
            r = min(idle, key=lambda a: (nx.shortest_path_length(g, j, a), a))
            idle.remove(r)
            refs.append(r)
            rung_routes.append(min(nx.all_shortest_paths(g, j, r)))
        layouts[str(length)] = dict(
            cycle=path,
            second_neighbour_routes=routes,
            reference_sites=refs,
            rung_routes=rung_routes,
        )
        tension_reads.extend(tension(length, k, rng) for k in (int(candidate.KAPPA), 0))
        if length == 18:
            continue
        for ones in (1, 2):
            a = basis(length, ones)
            for k, duty in (
                (int(candidate.KAPPA), 0),
                (int(candidate.KAPPA), 0.25),
                (int(candidate.KAPPA), 1),
                (0, 1),
            ):
                out.append(panel(length, ones, a, k, duty))
    return currents._portable(
        dict(geometry=layouts, panels=out, domain_energy=tension_reads)
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    data = run()
    with a.output.open("x") as f:
        f.write(json.dumps(data, indent=2, allow_nan=False) + "\n")
