"""Two-position native work interferometers on entire graph cuts."""

from collections import deque

import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.contrast_memory_spherical_v1 import law

SOURCES = ('host', 'excess', 'vacancy', 'write02', 'write20', 'write11')
SHOTS = 2**22


def mount_geometry(board, seam):
    distance = shortest_path(board.adjacency(), directed=False, unweighted=True, indices=list(seam)).min(axis=0)
    cuts = {}
    for radius in (1, 2, 3):
        cuts[radius] = [(a, b) if distance[a] < distance[b] else (b, a)
                        for a, b in board.seams if {distance[a], distance[b]} == {radius, radius + 1}]
    outside = sorted((min((int(distance[p]), p) for p in range(12 * owner, 12 * owner + 12)), owner)
                     for owner in range(board.cells) if owner not in {a // 12 for a in seam})
    launches = [site for (_, site), _ in outside[:3]]
    return cuts, launches


def private_route(board, start, avoid):
    queue = deque([[start]])
    seen = {start, *avoid}
    while queue:
        path = queue.popleft()
        if board.degree[path[-1]] == 5:
            return path
        for other in board.neighbours[path[-1]]:
            if other not in seen:
                seen.add(other)
                queue.append([*path, other])
    raise ValueError('no spare native reference reachable')


def native_analyzers():
    analyzers = {}
    for sign in (1, -1):
        n = 1 + sign
        indices = [3 * n + 1, 3 + n]
        local = law.seam_term(8)[np.ix_(indices, indices)]
        split = expm(-1j * np.pi / 2 * local)
        wait = np.pi / 48 if sign == 1 else 3 * np.pi / 16
        private = np.diag(np.exp(-8j * wait * np.array([1, n * n])))
        phase = private[1, 1] / private[0, 0]
        np.testing.assert_allclose(phase, -1j, atol=2e-14)
        analyzers[sign] = (split[:, 0], split.conj().T, split.conj().T @ private)
    return analyzers


def material_costs(board, seam, source, local):
    """Energy differences from adding the test load at one port, with the source independent."""
    occupation = np.ones(board.n, int)
    occupation[seam[0]] += {'excess': 1, 'vacancy': -1}.get(source, 0)
    costs = {}
    for sign in (1, -1):
        change = np.zeros(board.n)
        for p in range(board.n):
            if p in seam:
                continue
            n = 1 + sign
            # Source and test configurations see identical central written memory at every boundary arm.
            change[p] = 8 * (6 - board.degree[p]) * (n*n - 1)
            for other in board.neighbours[p]:
                m = occupation[other]
                change[p] += local[3*n + m, 3*n + m] - local[3 + m, 3 + m]
        costs[sign] = change
    return costs


def work_counts(gap, time, analyzer, rng):
    calibrated, x_read, y_read = analyzer
    state = np.array([1, np.exp(-1j * time * gap)]) * calibrated
    result = []
    for read in (x_read, y_read):
        p = abs((read @ state)[0])**2
        plus = int(rng.binomial(SHOTS, min(1, max(0, p))))
        result.append(plus)
    return result


def amplitude(counts):
    return complex(*(2 * np.asarray(counts) / SHOTS - 1))


def run(board, seams):
    rng = np.random.default_rng(191101)
    analyzer = native_analyzers()
    local = law.seam_term(8)
    banks = []
    for seam in seams:
        cuts, _ = mount_geometry(board, seam)
        costs = {s: material_costs(board, seam, s, local) for s in SOURCES}
        for radius, edges in cuts.items():
            for time in (1/64, 1/32):
                reads = []
                for u, v in edges:
                    route = private_route(board, v, {u, *seam})
                    counts = []
                    for source in SOURCES:
                        counts.append([work_counts(costs[source][p][v] - costs[source][p][u], time, analyzer[p], rng)
                                       for p in (1, -1)])
                    reads.append(dict(edge=[u, v], receiving_owner=v // 12, phase_reference_route=route,
                                      plus_counts=counts))
                banks.append(dict(seam=list(seam), radius=radius, dwell=time, edges=reads))
    bank_count = sum(len(b['edges']) for b in banks) * len(SOURCES) * 4
    component = np.sqrt(2 * np.log(2 * bank_count / .01) / SHOTS)
    summaries = []
    for bank in banks:
        seen = np.zeros(len(SOURCES))
        for entry in bank['edges']:
            reference = [amplitude(c) for c in entry['plus_counts'][0]]
            for idx, source in enumerate(entry['plus_counts']):
                phases = [-np.angle(amplitude(c) / r) / bank['dwell'] for c, r in zip(source, reference)]
                seen[idx] += (phases[0] - phases[1]) / 2
        # Each of the two complex quadratures has modulus one before shot noise.
        radius = 2 * np.arcsin(np.sqrt(2) * component) * len(bank['edges']) / bank['dwell']
        summaries.append(dict(seam=bank['seam'], radius=bank['radius'], dwell=bank['dwell'],
                              cuts=len(bank['edges']), owners=sorted({row['receiving_owner'] for row in bank['edges']}),
                              flux=seen.tolist(), simultaneous_flux_radius=float(radius)))
    return dict(shots=SHOTS, seed=191101, source_order=SOURCES, test_sign_order=[1, -1],
                count_order='X-plus, Y-plus; each other exit is shots minus plus', banks=banks, summaries=summaries)
