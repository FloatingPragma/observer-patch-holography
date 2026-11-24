"""Energy-resolved three-wall measurements over all positions on an open federation."""

import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations
import json
from math import sqrt
from pathlib import Path

import networkx as nx
import numpy as np
import sympy as sp

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier
from ...row_13_walls.energy import is_wall
from . import belt, localized, orbital_singlet


@lru_cache
def local_data(n):
    seams = carrier.cell().seams
    words = sorted(sum(1 << p for p in v) for v in combinations(range(12), n))
    selected = {
        sum(1 << ports[j] for j in range(10) if w >> j & 1)
        for ports in belt.belts().values()
        for w in range(1024)
        if w.bit_count() == n and is_wall(w)
    }
    rows = {}
    for w in words:
        edges = [e for e in seams if ((w >> e[0]) ^ (w >> e[1])) & 1]
        stay = sum((w ^ (1 << u) ^ (1 << v)) in selected for u, v in edges)
        pairs = (5 * n - len(edges)) // 2
        rows[w] = (len(edges), stay, pairs)
    back = max(i[1] for w, i in rows.items() if w not in selected)
    return selected, rows, back


def piece(n, gates, k):
    selected, rows, _ = local_data(n)
    mask = sum(1 << j for j in gates)
    lo, hi, exits, trials = [], [], [], []
    for w, (cuts, stay, pairs) in rows.items():
        outgoing = (w & mask).bit_count()
        diagonal = (cuts + outgoing) / 2 - 2 * k * pairs
        trials.append((diagonal, w, pairs, cuts + outgoing))
        if w in selected:
            lo.append(diagonal - stay / 2)
            hi.append(diagonal + stay / 2)
            exits.append(cuts + outgoing - stay)
    return dict(lo=min(lo), hi=max(hi), out=max(exits), trial=min(trials))


def occupied_owners(board):
    graph = nx.Graph((a // 12, b // 12) for a, b in board.external)
    assert len(graph.edges) == len(board.external)
    assert max(sum(p in edge for edge in board.external) for p in range(board.n)) == 1
    max_edges = 3 if any(nx.triangles(graph).values()) else 2
    distance = dict(nx.all_pairs_shortest_path_length(graph))
    owners = [0]
    for _ in range(2):
        owners.append(
            max(
                (j for j in graph if j not in owners),
                key=lambda j: (min(distance[o][j] for o in owners), -j),
            )
        )
    distances = [distance[a][b] for a, b in combinations(owners, 2)]
    assert min(distances) > 1
    return max_edges, owners, distances


def band_bound(norm2, gap):
    return min(1.0, norm2 / gap**2) if gap > 0 else 1.0


def fraction(norm2, gap):
    return norm2 / (norm2 + gap**2) if gap > 0 else 1.0


def face_band(board, k):
    """Spectral separation for faces in the complete three-load configuration graph."""
    graph = [set(row) for row in board.neighbours]
    faces = {(a, b, c) for a, b in board.seams for c in graph[a] & graph[b] if b < c}
    back = Counter()
    deg = []
    for f in faces:
        w = set(f)
        exits = 0
        for p in f:
            for q in graph[p] - w:
                child = tuple(sorted(w - {p} | {q}))
                assert child not in faces
                back[child] += 1
                exits += 1
        deg.append(exits)
    norm2 = max(deg) * max(back.values()) / 4
    # K is positive. A nonface has zero, one or two occupied edges.
    ceiling = max(deg) / 2 - 6 * k
    floor = min(-2 * k * p for p in (0, 1, 2))
    separated = min(0, -2 * k)
    gap = floor - ceiling
    loss = band_bound(norm2, gap)
    return dict(
        faces=len(faces),
        degree_range=[min(deg), max(deg)],
        maximum_parents=max(back.values()),
        coupling=k,
        face_compression_maximum=ceiling,
        complement_floor=floor,
        spectral_gap_lower=max(0, gap),
        separated_pair_single_floor=separated,
        energy_below_separated_at_least=max(0, separated - ceiling),
        face_probability_in_low_band_lower=1 - loss if gap > 0 else 0,
        face_probability_after_arbitrary_wait_lower=max(0, 1 - 4 * loss)
        if gap > 0
        else 0,
    )


def contact_read(board, owners, trials, k, rng):
    occupied = {
        12 * c + p
        for c, row in zip(owners, trials, strict=True)
        for p in range(12)
        if row[1] >> p & 1
    }
    pairs = sum(a in occupied and b in occupied for a, b in board.seams)
    unequal = sum((a in occupied) != (b in occupied) for a, b in board.seams)
    assert (pairs, unequal) == (sum(r[2] for r in trials), sum(r[3] for r in trials))
    shots = 2**18
    odd = int(rng.binomial(unequal * shots, 0.5))
    counts = [pairs * shots, odd, len(board.seams) * shots - pairs * shots - odd]
    return dict(
        shots_per_seam=shots,
        total_both_odd_other=counts,
        energy=(odd - 2 * k * counts[0]) / shots,
    )


def energy_budget():
    h0, h1 = law.local_term_parts()
    rows = []
    for k in (int(candidate.KAPPA), 0):
        h = h0 + k * h1
        bound = 1 + 4 * k
        # Positivity and the largest eigenvalue also follow from the unitary
        # conjugation to label-reflection parity plus k times squared contrast.
        spectrum = np.linalg.eigvalsh(h)
        assert spectrum[0] >= -1e-12 and spectrum[-1] <= bound + 1e-12
        for b in range(3):
            flags = np.array(
                [
                    int(x != b) + int(y != b) + int(m != 0)
                    for m in range(11)
                    for x, y in law.SRC
                ]
            )
            assert np.linalg.norm(h[:, law.SRC.index((b, b))]) < 1e-12
            slack = np.linalg.eigvalsh(np.diag(bound * flags) - h)[0]
            assert slack >= -1e-12
            rows.append(
                dict(
                    kappa=k,
                    host=b,
                    largest_local_energy=spectrum[-1],
                    local_slack_minimum=slack,
                    lower_per_material=-k * b * b,
                    upper_per_material=6 * bound,
                    upper_per_mark=bound,
                )
            )
    return rows


def run():
    shapes = localized.owner_shapes()
    data = []
    faces = []
    rng = np.random.default_rng(131214)
    for level in range(3):
        board = carrier.federation(level)
        max_edges, owners, distances = occupied_owners(board)
        gates = [
            tuple(np.flatnonzero(board.degree[12 * c : 12 * c + 12] == 6))
            for c in owners
        ]
        kinds = [
            g
            for g, occurrences in shapes.items()
            if any(s == board.cells for s, _ in occurrences)
        ]
        # Every owner offers an exact one-load compact orbital in blank memory.
        # Choose the highest compact eigenvalue shared by enough separate owners.
        options = []
        for shape in kinds:
            detail, _, _ = orbital_singlet.compact_levels(shape)
            poly = sp.Poly(
                sp.sympify(detail["twice_kinetic_polynomial"].replace("lambda", "x"))
            )
            interval = poly.intervals(eps=sp.Rational(1, 10**10))[-1][0]
            lower = float(interval[0]) / 2
            count = next(v for s, v in shapes[shape] if s == board.cells)
            options.extend([lower] * count)
        for k in (int(candidate.KAPPA), 0, -int(candidate.KAPPA)):
            faces.append(dict(cells=board.cells, **face_band(board, k)))
            parts = {(n, g): piece(n, g, k) for n in (3, 4) for g in kinds}
            for loads in ((3, 3, 4), (3, 4, 4)):
                n = sum(loads)
                trials = [
                    piece(i, g, k)["trial"] for i, g in zip(loads, gates, strict=True)
                ]
                energy = sum(row[0] for row in trials)
                # Inter-owner double occupancy reduces K by one relative to the
                # sum of the two single-owner K diagonals. A W-to-W transfer
                # instead contributes one off-diagonal of magnitude 1/2.
                v = -(1 + 2 * k)
                low = sum(
                    min(parts[i, g]["lo"] for g in kinds) for i in loads
                ) + max_edges * min(0, v, -0.5)
                high = sum(
                    max(parts[i, g]["hi"] for g in kinds) for i in loads
                ) + max_edges * max(0, v, 0.5)
                forward = sum(max(parts[i, g]["out"] for g in kinds) for i in loads)
                back = max(local_data(i)[2] for i in loads) + max_edges
                norm2 = forward * back / 4
                upper_line = sum(sorted(options, reverse=True)[:n])
                data.append(
                    dict(
                        cells=board.cells,
                        coupling=k,
                        loads=loads,
                        owners=owners,
                        distances=distances,
                        separated_words=[row[1] for row in trials],
                        separated_pairs=sum(row[2] for row in trials),
                        separated_unequal=sum(row[3] for row in trials),
                        separated_mean=energy,
                        separated_read=contact_read(board, owners, trials, k, rng),
                        compact_single_fragment_line=upper_line,
                        wall_gershgorin=[low, high],
                        cut_degrees=[forward, back],
                        cut_norm_upper=sqrt(norm2),
                        individual_low_eigenstate_upper=fraction(norm2, low - energy),
                        individual_high_eigenstate_upper=fraction(
                            norm2, upper_line - high
                        ),
                        entire_lower_band_wall_upper=band_bound(norm2, low - energy),
                        entire_upper_band_wall_upper=band_bound(
                            norm2, upper_line - high
                        ),
                    )
                )
    return currents._portable(
        dict(
            wall_patterns={str(n): len(local_data(n)[0]) for n in (3, 4)},
            boards=data,
            finite_excitation_energy=energy_budget(),
            face_bands=faces,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    data = run()
    with args.output.open("x") as f:
        f.write(json.dumps(data, indent=2, allow_nan=False) + "\n")
