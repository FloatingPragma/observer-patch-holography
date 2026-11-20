"""The cage: a profile on a ring of four seams that the scalar law keeps.

With `C` the signed seam-by-wedge incidence (+1 at seam `ab`, -1 at seam `bc` of the wedge `(a, b, c)`) and `Q`
the unsigned port-by-seam incidence, the law of one object is

    D = [[m, c C], [c C^T, -m]],        m = 1,  c = 1/sqrt2,        K = c^2 C C^T .

Two facts in integers carry the construction:

    C C^T = diag(d_a + d_b) - Q^T Q        and        Q e = 0  for the profile of a ring,

so that `C C^T e = 10 e` when the four ports of the ring have five seams. With `e' = C^T e / sqrt(10)` the law
acts on the pair `(e, e')` as `[[m, sigma], [sigma, -m]]`, `sigma = sqrt(10) c`, and couples it to nothing else.

One cell is held in dense integer arrays. A federation is held in sparse integer arrays, and every statement
about it here is a product of integers.
"""
from __future__ import annotations

import itertools
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from ..carrier import EDGES, require

HERE = Path(__file__).resolve().parent
CARRIERS = HERE.parents[1] / 'contrast_memory_spherical_v1' / 'carriers'
RING = (0, 2, 4, 6)                          # the ring of the manual, in the labels of the manual's cell
REST = 1.0                                   # m
HOP = 1 / math.sqrt(2)                       # c
SIGMA = math.sqrt(10) * HOP                  # sqrt5
ENERGY = math.hypot(REST, SIGMA)             # sqrt6
STEP = 0.125                                 # one step of the observer's protocol, in units of 1/m
STRENGTH = 64.0                              # the converter of section 6.1
T_COHERENT = math.pi / (2 * math.sqrt(2 * ENERGY * (ENERGY - REST)))
T_SHARP = math.pi * REST / (2 * SIGMA * ENERGY)


def incidences(edges, ports, sparse=False):
    """Seams, wedges, neighbours, the seam index, and the integer incidences C and Q of a carrier."""
    edges = sorted({tuple(sorted((int(a), int(b)))) for a, b in edges})
    neighbours = [[] for _ in range(ports)]
    for a, b in edges:
        neighbours[a].append(b)
        neighbours[b].append(a)
    neighbours = [sorted(x) for x in neighbours]
    wedges = [(a, b, c) for b in range(ports) for a, c in itertools.combinations(neighbours[b], 2)]
    index = {e: i for i, e in enumerate(edges)}
    rows = [index[tuple(sorted(pair))] for a, b, c in wedges for pair in ((a, b), (b, c))]
    cols = [w for w in range(len(wedges)) for _ in (0, 1)]
    signs = [1, -1] * len(wedges)
    C = sp.csr_array((signs, (rows, cols)), shape=(len(edges), len(wedges)), dtype=np.int64)
    Q = sp.csr_array(([1] * (2 * len(edges)), ([p for e in edges for p in e], [i for i in range(len(edges)) for _ in (0, 1)])),
                     shape=(ports, len(edges)), dtype=np.int64)
    if not sparse:
        C, Q = C.toarray(), Q.toarray()
    return dict(edges=edges, wedges=wedges, neighbours=neighbours, index=index, C=C, Q=Q, ports=ports,
                degree=np.array([len(x) for x in neighbours], dtype=np.int64))


def cell():
    """One cell in the labels of the manual."""
    return incidences(EDGES, 12)


def federation(level):
    """The federation of 20, 80 or 320 cells (level 0, 1, 2), with sparse incidences."""
    require(level in (0, 1, 2), 'federation level must be 0, 1 or 2')
    with np.load(CARRIERS / f'FEDERATION_LEVEL{level}.npz', allow_pickle=False) as z:
        edges = np.asarray(z['edges'])
    return incidences(edges, int(edges.max()) + 1, sparse=True)


def seam_laplacian(car):
    """C C^T, the Laplacian of the graph of seams, in integers."""
    return car['C'] @ car['C'].T


def line_graph_identity(car):
    """C C^T = diag(d_a + d_b) - Q^T Q, in integers. Returns the largest difference, which is zero."""
    degree = car['degree']
    diagonal = np.array([degree[a] + degree[b] for a, b in car['edges']], dtype=np.int64)
    if sp.issparse(car['C']):
        difference = seam_laplacian(car) - (sp.diags_array(diagonal, dtype=np.int64) - car['Q'].T @ car['Q'])
        return int(abs(difference).max()) if difference.nnz else 0
    return int(np.abs(seam_laplacian(car) - (np.diag(diagonal) - car['Q'].T @ car['Q'])).max())


def ring_profile(car, ring):
    """The integer profile of a ring of four ports: +1, -1, +1, -1 on its four seams. The profile is half of it."""
    require(len(ring) == 4 and len(set(ring)) == 4, 'a ring has four different ports')
    out = np.zeros(len(car['edges']), dtype=np.int64)
    for k in range(4):
        seam = tuple(sorted((ring[k], ring[(k + 1) % 4])))
        require(seam in car['index'], 'the ports of a ring are joined in a cycle')
        out[car['index'][seam]] = (-1) ** k
    return out


def rings(car, ports=None):
    """All rings of four ports through the given ports, each once, as tuples that start at their smallest port."""
    nb = [set(x) for x in car['neighbours']]
    pool = range(car['ports']) if ports is None else ports
    seen = set()
    for a in pool:
        for b in nb[a]:
            for c in nb[b]:
                if c == a:
                    continue
                for d in nb[c]:
                    if d == a or d == b or a not in nb[d]:
                        continue
                    k = min(range(4), key=lambda i: (a, b, c, d)[i])
                    ring = ((a, b, c, d) * 2)[k:k + 4]
                    seen.add(ring if ring[1] < ring[3] else (ring[0], ring[3], ring[2], ring[1]))
    return sorted(seen)


def rank(rows):
    """Exact rank of an integer matrix, by elimination in fractions."""
    m = [[Fraction(int(x)) for x in row] for row in rows]
    r = 0
    for col in range(len(m[0]) if m else 0):
        pivot = next((i for i in range(r, len(m)) if m[i][col] != 0), None)
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        for i in range(r + 1, len(m)):
            if m[i][col] != 0:
                f = m[i][col] / m[r][col]
                m[i] = [x - f * y for x, y in zip(m[i], m[r])]
        r += 1
        if r == len(m):
            break
    return r


def small_supports(car, largest=4):
    """The sets of at most `largest` seams of one cell that carry a profile with Q e = 0.

    A profile on a set of seams exists exactly when the columns of Q on those seams are dependent. A set is
    listed when its columns are dependent and no column can be left out, which makes it the support of one
    profile.
    """
    Q = car['Q']
    out = {k: [] for k in range(1, largest + 1)}
    columns = [Q[:, i].tolist() for i in range(Q.shape[1])]
    for k in range(1, largest + 1):
        for subset in itertools.combinations(range(len(columns)), k):
            ports = [sum(columns[i][p] for i in subset) for p in range(car['ports'])]
            if any(x == 1 for x in ports):             # a port with one seam of the set cannot cancel
                continue
            if rank([columns[i] for i in subset]) == k - 1:
                out[k].append(tuple(car['edges'][i] for i in subset))
    return out


def pair(car, profile, rest=REST, hop=HOP):
    """The unit vectors e and e' of a profile, the block h, and the largest residual of D (e, e') = (e, e') h.

    The products use the incidences alone, so the same function serves one cell and a federation.
    """
    C = car['C']
    norm = math.sqrt(float(profile @ profile))
    image = C.T @ profile
    length = math.sqrt(float(image @ image))
    e, ep = profile / norm, image / length
    sigma = hop * length / norm
    h = np.array([[rest, sigma], [sigma, -rest]])
    # D (e, 0) = (m e, c C^T e) and D (0, e') = (c C e', -m e')
    residual = max(float(np.abs(hop * (C.T @ e) - sigma * ep).max()), float(np.abs(hop * (C @ ep) - sigma * e).max()))
    return e, ep, h, residual


def modes(car, profile, rest=REST, hop=HOP):
    """The pair (e, e') as two columns on the seams and wedges of one cell, the block h, and the residual."""
    E, W = car['C'].shape
    e, ep, h, residual = pair(car, profile, rest, hop)
    return np.column_stack([np.r_[e, np.zeros(W)], np.r_[np.zeros(E), ep]]), h, residual


def law(car, rest=REST, hop=HOP):
    """The law of one object on the seams and wedges of one cell."""
    require(not sp.issparse(car['C']), 'the dense law is built for one cell')
    E, W = car['C'].shape
    return np.block([[rest * np.eye(E), hop * car['C']], [hop * car['C'].T, -rest * np.eye(W)]])


def changed_law_term(car):
    """Seam hopping round every port: the term of the changed law, on the seams and wedges of one cell.

    At every port its five seams are joined in the cyclic order of the port's neighbours, each joined pair with
    the amplitude one. The term keeps the alphabet, the routing, the number of objects and the symmetries of the
    cell. It is added to the law in every instrument and in every wait of a run under the changed law.
    """
    E, W = car['C'].shape
    nb, index = car['neighbours'], car['index']
    J = np.zeros((E + W, E + W))
    for b in range(car['ports']):
        require(len(nb[b]) == 5, 'the changed law is defined on a cell')
        order = [nb[b][0]]
        while len(order) < 5:
            order.append([x for x in nb[order[-1]] if x in nb[b] and x not in order][0])
        for k in range(5):
            i = index[tuple(sorted((b, order[k])))]
            j = index[tuple(sorted((b, order[(k + 1) % 5])))]
            J[i, j] += 1
            J[j, i] += 1
    return J


def circulation(car, port=0):
    """The observer's reference at a port: three transits round it, at the level -m exactly."""
    a, c, f = car['neighbours'][port][:3]
    E = len(car['edges'])
    where = {w: i for i, w in enumerate(car['wedges'])}
    out = np.zeros(E + len(car['wedges']))
    for wedge, sign in (((a, port, c), 1), ((c, port, f), 1), ((a, port, f), -1)):
        out[E + where[wedge]] = sign / math.sqrt(3)
    return out


def exact_rings(car):
    """The rings whose four ports have five seams, each checked in integers: Q e = 0 and C C^T e = 10 e."""
    five = [p for p in range(car['ports']) if car['degree'][p] == 5]
    keep = set(five)
    laplacian = seam_laplacian(car)
    out = []
    for ring in rings(car, five):
        if all(p in keep for p in ring):
            e = ring_profile(car, ring)
            require(not np.any(car['Q'] @ e) and np.array_equal(laplacian @ e, 10 * e), 'a ring on ports with five seams lost its level')
            out.append(ring)
    return out


def leaking_ring(car):
    """A ring through a port with six seams, and the largest entry of C C^T e - 10 e on it, which is not zero."""
    six = next(p for p in range(car['ports']) if car['degree'][p] == 6)
    ring = next(r for r in rings(car, [six]) if six in r)
    e = ring_profile(car, ring)
    return ring, int(np.abs(seam_laplacian(car) @ e - 10 * e).max())


def flat_branch(car):
    """Counts of the flat branch of one cell: rank of Q, the profiles with Q e = 0, and the rank of the rings."""
    Q = car['Q']
    r = rank(Q.tolist())
    all_rings = rings(car)
    span = rank([ring_profile(car, ring).tolist() for ring in all_rings])
    return dict(rank_of_Q=r, profiles=Q.shape[1] - r, rings=len(all_rings), rank_of_the_rings=span)


def kept_profiles(car):
    """The profiles with Q e = 0 on seams between ports with five seams, counted exactly.

    Such a profile has C C^T e = 10 e by the incidence identity. The ports with five seams and the seams between
    them fall into connected pieces, and the count of a piece is its number of seams minus the rank of its
    incidence. Returns the count of every piece, in the order of its smallest port.
    """
    five = {p for p in range(car['ports']) if car['degree'][p] == 5}
    inner = [(a, b) for a, b in car['edges'] if a in five and b in five]
    parent = {p: p for p in five}

    def find(p):
        while parent[p] != p:
            parent[p] = parent[parent[p]]
            p = parent[p]
        return p
    for a, b in inner:
        parent[find(a)] = find(b)
    pieces = {}
    for a, b in inner:
        pieces.setdefault(find(a), []).append((a, b))
    out = []
    for root in sorted(pieces, key=lambda r: min(p for seam in pieces[r] for p in seam)):
        seams = pieces[root]
        ports = sorted({p for seam in seams for p in seam})
        where = {p: i for i, p in enumerate(ports)}
        rows = [[0] * len(seams) for _ in ports]
        for i, (a, b) in enumerate(seams):
            rows[where[a]][i] = rows[where[b]][i] = 1
        out.append(len(seams) - rank(rows))
    return out
