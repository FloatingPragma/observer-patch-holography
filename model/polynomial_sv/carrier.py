"""The carrier: one cell (the icosahedral echosahedron) and the federations of 20, 80 and 320 cells.

A cell has 12 ports, 30 seams and 20 faces with the incidence of the icosahedron. Port `12 k + j` of a
federation is port `j` of cell `k`; every cell keeps all 30 of its seams with the labelling of `CELL_SEAMS`.
External seams join ports of different cells and form a partial matching: every cell has three external ports,
so a port has five or six seams. The federation arrays are in `carriers/FEDERATION_LEVEL{0,1,2}.npz`
(20, 80 and 320 cells).
"""
from __future__ import annotations

import collections
import hashlib
import json
from functools import lru_cache
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# The 30 seams of one cell, in the labelling every cell of the federations uses.
CELL_SEAMS = (
    (0, 1), (0, 5), (0, 7), (0, 10), (0, 11), (1, 5), (1, 7), (1, 8), (1, 9), (2, 3),
    (2, 4), (2, 6), (2, 10), (2, 11), (3, 4), (3, 6), (3, 8), (3, 9), (4, 5), (4, 9),
    (4, 11), (5, 9), (5, 11), (6, 7), (6, 8), (6, 10), (7, 8), (7, 10), (8, 9), (10, 11),
)

NAMES = {'cell': 'one cell', 0: '20 cells', 1: '80 cells', 2: '320 cells'}
# The emitter seam of each carrier: an internal seam (0, 1) or (2, 3) of the first cell; one end of it carries
# an external seam on every federation.
EMITTER_SEAM = {'cell': (0, 1), 0: (0, 1), 1: (2, 3), 2: (2, 3)}


class Carrier:
    """A graph of ports and seams with the derived tables the law needs."""

    def __init__(self, name: str, ports: int, seams):
        self.name = name
        self.n = int(ports)
        self.seams = sorted(tuple(sorted(map(int, e))) for e in seams)
        self.neighbours = [[] for _ in range(self.n)]
        for a, b in self.seams:
            self.neighbours[a].append(b)
            self.neighbours[b].append(a)
        for x in self.neighbours:
            x.sort()
        self.degree = np.array([len(x) for x in self.neighbours], dtype=np.int64)
        table = np.full((self.n, max(6, int(self.degree.max()))), -1, dtype=np.int64)
        for p, x in enumerate(self.neighbours):
            table[p, :len(x)] = x
        self.table = table
        self.external = [(a, b) for a, b in self.seams if a // 12 != b // 12]

    @property
    def cells(self) -> int:
        return self.n // 12

    def adjacency(self):
        from scipy import sparse
        a, b = np.array(self.seams).T
        return sparse.coo_matrix((np.ones(2 * len(a)), (np.r_[a, b], np.r_[b, a])), shape=(self.n, self.n)).tocsr()

    def laplacian(self):
        from scipy import sparse
        return (sparse.diags(self.degree.astype(float)) - self.adjacency()).tocsr()

    def faces(self):
        """Triangles of the carrier (all lie inside cells)."""
        nb = [set(x) for x in self.neighbours]
        return sorted((a, b, c) for a, b in self.seams for c in nb[a] & nb[b] if c > b)


def cell() -> Carrier:
    return Carrier('one cell', 12, CELL_SEAMS)


@lru_cache(maxsize=None)
def _federation_arrays(level: int):
    z = np.load(HERE / 'carriers' / f'FEDERATION_LEVEL{level}.npz')
    return int(len(z['degree'])), tuple(tuple(map(int, e)) for e in z['edges']), tuple(map(int, z['degree']))


def federation(level: int) -> Carrier:
    n, edges, degree = _federation_arrays(level)
    car = Carrier(NAMES[level], n, edges)
    if tuple(car.degree.tolist()) != degree:
        raise ValueError('stored degrees disagree with the seams')
    return car


def carrier(key) -> Carrier:
    return cell() if key == 'cell' else federation(key)


def icosahedron_from_coordinates():
    """The 12 vertices `(0, +-1, +-phi)` and cyclic shifts; seams join vertices at the shortest distance."""
    phi = (1 + 5 ** 0.5) / 2
    pts = []
    for s1 in (1, -1):
        for s2 in (phi, -phi):
            pts += [(0, s1, s2), (s1, s2, 0), (s2, 0, s1)]
    pts = np.array(pts)
    d = np.linalg.norm(pts[:, None] - pts[None], axis=2)
    seams = [(i, j) for i in range(12) for j in range(i + 1, 12) if abs(d[i, j] - 2) < 1e-9]
    return pts, seams


def isomorphism(seams_a, seams_b, n=12):
    """A vertex map taking the graph `seams_a` onto `seams_b`, found by backtracking (None if none exists)."""
    A = np.zeros((n, n), bool)
    B = np.zeros((n, n), bool)
    for a, b in seams_a:
        A[a, b] = A[b, a] = True
    for a, b in seams_b:
        B[a, b] = B[b, a] = True
    m = [-1] * n
    used = [False] * n

    def place(i):
        if i == n:
            return True
        for v in range(n):
            if used[v] or A[i].sum() != B[v].sum():
                continue
            if all(A[i, j] == B[v, m[j]] for j in range(i)):
                m[i], used[v] = v, True
                if place(i + 1):
                    return True
                m[i], used[v] = -1, False
        return False

    return list(m) if place(0) else None


def automorphisms(seams, n=12):
    """All vertex permutations that map the seam set onto itself."""
    adj = np.zeros((n, n), bool)
    for a, b in seams:
        adj[a, b] = adj[b, a] = True
    out = []
    m = [-1] * n
    used = [False] * n

    def place(i):
        if i == n:
            out.append(tuple(m))
            return
        for v in range(n):
            if not used[v] and all(adj[i, j] == adj[v, m[j]] for j in range(i)):
                m[i], used[v] = v, True
                place(i + 1)
                m[i], used[v] = -1, False

    place(0)
    return out


def cell_checks() -> dict:
    c = cell()
    faces = c.faces()
    pts, seams = icosahedron_from_coordinates()
    iso = isomorphism(c.seams, seams)
    A = c.adjacency().toarray()
    w = np.sort(np.linalg.eigvalsh(A))
    expected = np.sort([5.0] + [5 ** 0.5] * 3 + [-1.0] * 5 + [-(5 ** 0.5)] * 3)
    autos = automorphisms(c.seams)
    return {'ports': c.n, 'seams': len(c.seams), 'faces': len(faces),
            'degree_five_everywhere': bool((c.degree == 5).all()),
            'isomorphic_to_coordinate_icosahedron': iso is not None,
            'adjacency_spectrum_5_sqrt5x3_minus1x5_minus_sqrt5x3': bool(np.abs(w - expected).max() < 1e-12),
            'automorphisms': len(autos),
            'every_seam_in_two_faces': all(sum(1 for f in faces if a in f and b in f) == 2 for a, b in c.seams)}


def federation_checks(car: Carrier) -> dict:
    """The hypotheses of the binding theorem, counted on the carrier."""
    nb = [set(x) for x in car.neighbours]
    ext_ports = collections.Counter(p for e in car.external for p in e)
    cell0 = sorted(CELL_SEAMS)
    by_cell = collections.defaultdict(list)
    for a, b in car.seams:
        if a // 12 == b // 12:
            by_cell[a // 12].append((a % 12, b % 12))
    worst = 0
    for a in range(car.n):
        two = collections.Counter(y for x in nb[a] for y in nb[x] if y != a and y not in nb[a])
        worst = max(worst, max(two.values(), default=0))
    ends = collections.Counter((len(nb[a] & nb[b]), a // 12 == b // 12) for a, b in car.seams)
    per_cell = collections.Counter(p // 12 for p in ext_ports)
    return {'cells': car.cells, 'ports': car.n, 'seams': len(car.seams), 'external_seams': len(car.external),
            'degrees': {str(k): int(v) for k, v in sorted(collections.Counter(car.degree.tolist()).items())},
            'port_degrees_five_or_six': bool(set(car.degree.tolist()) <= {5, 6}),
            'every_cell_complete': len(by_cell) == car.cells and all(sorted(v) == cell0 for v in by_cell.values()),
            'external_seams_form_a_matching': max(ext_ports.values()) == 1,
            'external_ports_per_cell': sorted(set(per_cell.values())),
            'common_neighbours_of_internal_seam_ends': sorted({k[0] for k in ends if k[1]}),
            'common_neighbours_of_external_seam_ends': sorted({k[0] for k in ends if not k[1]}),
            'most_common_neighbours_of_a_nonadjacent_pair': int(worst),
            'faces': len(car.faces())}


def seam_digest(car: Carrier) -> str:
    return hashlib.sha256(json.dumps(car.seams).encode()).hexdigest()


def run(levels=(0, 1)) -> dict:
    out = {'cell': cell_checks()}
    for level in levels:
        car = federation(level)
        entry = federation_checks(car)
        entry['seam_list_sha256'] = seam_digest(car)
        a, b = EMITTER_SEAM[level]
        entry['emitter_seam'] = [a, b]
        entry['emitter_seam_end_has_external_seam'] = any(p in e for e in car.external for p in (a, b))
        out[f'level{level}'] = entry
    return out


def check(result: dict) -> list[str]:
    """Names of failed carrier conditions in a result of `run`."""
    bad = []
    c = result['cell']
    for key in ('degree_five_everywhere', 'isomorphic_to_coordinate_icosahedron',
                'adjacency_spectrum_5_sqrt5x3_minus1x5_minus_sqrt5x3', 'every_seam_in_two_faces'):
        if not c[key]:
            bad.append('cell.' + key)
    if (c['ports'], c['seams'], c['faces'], c['automorphisms']) != (12, 30, 20, 120):
        bad.append('cell.counts')
    for key, entry in result.items():
        if not key.startswith('level'):
            continue
        for flag in ('port_degrees_five_or_six', 'every_cell_complete', 'external_seams_form_a_matching',
                     'emitter_seam_end_has_external_seam'):
            if not entry[flag]:
                bad.append(f'{key}.{flag}')
        if entry['most_common_neighbours_of_a_nonadjacent_pair'] > 2 or entry['external_ports_per_cell'] != [3]:
            bad.append(f'{key}.local_hypotheses')
    return bad
