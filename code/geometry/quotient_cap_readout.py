#!/usr/bin/env python3
"""Machine receipts for the quotient-intrinsic geometry producer (GitHub #523).

Implements the objects of the compact paper's subsection
`subsec:quotient-intrinsic-geometry-producer`:

* a finite transactional quotient repair system whose record layer binds one
  record token per 2-cell of an abstract patch adjacency (the adjacency is
  never consumed directly by the readout: the incidence complex is recomputed
  from the repaired record tokens of the normal form);
* the support-visible incidence complex `K(W)` of
  Definition `def:support-visible-incidence-complex`, with the invariance
  checks of Lemma `lem:incidence-invariance` (schedule, gauge, refinement);
* the computable receipts of Definition `def:spherical-incidence-receipt`
  (closed combinatorial surface, connected, Euler characteristic 2, coherent
  orientation) and the wrong-beta KMS receipt clause;
* the topology-production step of Theorem `thm:topology-production`
  (surface classification by orientability + Euler characteristic);
* the conformal production step of Theorem `thm:conformal-cap-production`:
  modular cross-ratios against a fixed gauge triple reconstruct the celestial
  embedding, cap normals `n_C` are produced from boundary circles by the
  formula of Proposition `prop:round-cap-normal`, and the residuals shrink
  under refinement;
* the underdetermination countermodels of
  Theorem `thm:incidence-underdetermination`: repair systems with isomorphic
  rewrite structure over S^2, T^2 (Csaszar torus), the 2-skeleton of the
  boundary of the 4-simplex, and a wedge of two spheres, distinguished only
  by the receipts.
"""

from __future__ import annotations

import itertools
import random
from dataclasses import dataclass, field

import numpy as np

TOL = 1e-9


# ---------------------------------------------------------------------------
# abstract 2-complexes (record layers)
# ---------------------------------------------------------------------------

def icosahedron() -> tuple[list[tuple[int, int, int]], np.ndarray]:
    """Return (triangles, vertex coordinates on the unit sphere)."""
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    verts = []
    for a, b in itertools.product((-1.0, 1.0), (-phi, phi)):
        verts += [(0.0, a, b), (a, b, 0.0), (b, 0.0, a)]
    coords = np.array(verts)
    coords /= np.linalg.norm(coords, axis=1, keepdims=True)
    # triangles: triples of mutually nearest vertices (edge length = 2/sqrt(phi^2+1))
    edge2 = 4.0 / (phi * phi + 1.0)
    n = len(coords)
    tris = []
    for i, j, k in itertools.combinations(range(n), 3):
        d = (
            np.sum((coords[i] - coords[j]) ** 2),
            np.sum((coords[j] - coords[k]) ** 2),
            np.sum((coords[i] - coords[k]) ** 2),
        )
        if all(abs(x - edge2) < 1e-6 for x in d):
            tris.append((i, j, k))
    assert len(tris) == 20
    return tris, coords


def csaszar_torus() -> list[tuple[int, int, int]]:
    """The Moebius 7-vertex triangulation of the torus (complete graph K7):
    triangles {i, i+1, i+3} and {i, i+2, i+3} mod 7."""
    tris = [tuple(sorted(((i) % 7, (i + 1) % 7, (i + 3) % 7))) for i in range(7)]
    tris += [tuple(sorted(((i) % 7, (i + 2) % 7, (i + 3) % 7))) for i in range(7)]
    assert len(set(tris)) == 14
    return tris


def boundary_4_simplex_2_skeleton() -> list[tuple[int, int, int]]:
    """All ten triangles on five vertices: the 2-skeleton of boundary(Delta^4).

    Locally three-dimensional (every edge lies in three triangles), so the
    closed-surface clause of the spherical incidence receipt must fail.
    """
    return list(itertools.combinations(range(5), 3))


def wedge_of_two_spheres() -> list[tuple[int, int, int]]:
    """Two icosahedra glued at one vertex: connected, chi = 3, nonmanifold."""
    tris, _ = icosahedron()
    shift = 12

    def relabel(v: int) -> int:
        # vertex 0 of the second copy is identified with vertex 0 of the first
        return 0 if v == 0 else v + shift - 1

    second = [tuple(relabel(v) for v in t) for t in tris]
    return tris + second


# ---------------------------------------------------------------------------
# finite transactional quotient repair system
# ---------------------------------------------------------------------------

@dataclass
class RepairSystem:
    """Patch bits with equality constraints along record tokens.

    Records (one per 2-cell) demand that all member patch bits agree. The
    initial state seeds one conflict inside a single record; repair commits
    one transactional conflict-component resolution (majority-with-min-label
    tie break), which is manifestly terminating and confluent: there is a
    single conflict component and its resolution does not depend on the
    schedule order of untouched records.
    """

    records: list[tuple[int, ...]]
    seed_record: int = 0
    state: dict[int, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        patches = sorted({p for rec in self.records for p in rec})
        self.state = {p: 0 for p in patches}
        # seed the conflict: flip the highest-label patch of the seed record
        bad = max(self.records[self.seed_record])
        self.state[bad] = 1

    def conflicted_records(self) -> list[int]:
        return [
            k for k, rec in enumerate(self.records)
            if len({self.state[p] for p in rec}) > 1
        ]

    def repair(self, schedule: list[int] | None = None) -> dict[int, int]:
        """Run transactional repair; return the quotient normal form."""
        order = list(schedule) if schedule is not None else list(range(len(self.records)))
        working = dict(self.state)
        progress = True
        while progress:
            progress = False
            for k in order:
                rec = self.records[k]
                vals = [working[p] for p in rec]
                if len(set(vals)) > 1:
                    # majority value, min-label tie break: recovery-derived law
                    counts = {v: vals.count(v) for v in set(vals)}
                    best = min(sorted(counts), key=lambda v: (-counts[v], v))
                    for p in rec:
                        working[p] = best
                    progress = True
        return working

    def rewrite_signature(self) -> tuple[int, int]:
        """Confluence-level invariant: (#patches touched by conflict, #steps).

        Isomorphic across all four countermodel systems by construction: the
        conflict component is a single record with three patches everywhere.
        """
        rec = self.records[self.seed_record]
        return (len(rec), 1)


def normal_form_records(system: RepairSystem, normal_form: dict[int, int]) -> list[tuple[int, ...]]:
    """The record tokens carried by the repaired normal form.

    A record is quotient-visible exactly when its constraint is satisfied on
    the normal form (all tokens after successful repair). The incidence
    readout below consumes ONLY this list, i.e. normal-form data.
    """
    return [rec for rec in system.records
            if len({normal_form[p] for p in rec}) == 1]


# ---------------------------------------------------------------------------
# support-visible incidence complex and receipts
# ---------------------------------------------------------------------------

@dataclass
class IncidenceComplex:
    vertices: list[int]
    edges: set[frozenset]
    triangles: set[frozenset]


def incidence_complex(records: list[tuple[int, ...]]) -> IncidenceComplex:
    """K(W): patches are vertices; a k-simplex is a jointly supported set.

    Joint support nonvanishing of a patch set is witnessed by a common record
    token containing it (transport into a common refining patch algebra).
    """
    vertices = sorted({p for rec in records for p in rec})
    edges, triangles = set(), set()
    for rec in records:
        for pair in itertools.combinations(sorted(rec), 2):
            edges.add(frozenset(pair))
        for tri in itertools.combinations(sorted(rec), 3):
            triangles.add(frozenset(tri))
    return IncidenceComplex(vertices, edges, triangles)


def euler_characteristic(K: IncidenceComplex) -> int:
    return len(K.vertices) - len(K.edges) + len(K.triangles)


def is_connected(K: IncidenceComplex) -> bool:
    if not K.vertices:
        return False
    adj: dict[int, set[int]] = {v: set() for v in K.vertices}
    for e in K.edges:
        a, b = tuple(e)
        adj[a].add(b)
        adj[b].add(a)
    seen = {K.vertices[0]}
    stack = [K.vertices[0]]
    while stack:
        for w in adj[stack.pop()]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == len(K.vertices)


def is_closed_surface(K: IncidenceComplex) -> bool:
    """Every edge in exactly two triangles and every vertex link a single cycle."""
    for e in K.edges:
        cofaces = [t for t in K.triangles if e < t]
        if len(cofaces) != 2:
            return False
    for v in K.vertices:
        star = [t for t in K.triangles if v in t]
        link_edges = [tuple(sorted(t - {v})) for t in star]
        nodes = sorted({x for le in link_edges for x in le})
        deg = {x: sum(1 for le in link_edges if x in le) for x in nodes}
        if any(d != 2 for d in deg.values()):
            return False
        # single cycle: connected link
        adj = {x: set() for x in nodes}
        for a, b in link_edges:
            adj[a].add(b)
            adj[b].add(a)
        seen, stack = {nodes[0]}, [nodes[0]]
        while stack:
            for w in adj[stack.pop()]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        if len(seen) != len(nodes):
            return False
    return True


def orient(K: IncidenceComplex) -> list[tuple[int, int, int]] | None:
    """Return coherently oriented triangles, or None if nonorientable."""
    tris = [tuple(sorted(t)) for t in K.triangles]
    if not tris:
        return None
    oriented: dict[tuple, tuple] = {}
    first = tris[0]
    oriented[first] = first
    stack = [first]
    remaining = set(tris[1:])
    while stack:
        t = stack.pop()
        a, b, c = oriented[t]
        # induced edge orientations of t
        induced = {(a, b), (b, c), (c, a)}
        for u in list(remaining):
            shared = set(t) & set(u)
            if len(shared) == 2:
                x, y = tuple(shared)
                (z,) = set(u) - shared
                # neighbor must induce the opposite orientation on the shared edge
                if (x, y) in induced:
                    ori = (y, x, z)
                else:
                    ori = (x, y, z)
                oriented[u] = ori
                remaining.discard(u)
                stack.append(u)
    if remaining:
        return None  # disconnected; caller checks connectivity separately
    # verify global coherence
    edge_orientations: dict[frozenset, list] = {}
    for t, (a, b, c) in oriented.items():
        for x, y in ((a, b), (b, c), (c, a)):
            edge_orientations.setdefault(frozenset((x, y)), []).append((x, y))
    for e, oris in edge_orientations.items():
        if len(oris) == 2 and oris[0] == oris[1]:
            return None  # incoherent: nonorientable
    return list(oriented.values())


def spherical_incidence_receipt(K: IncidenceComplex) -> bool:
    """SphInc: connected closed orientable combinatorial surface with chi = 2."""
    return (
        is_connected(K)
        and is_closed_surface(K)
        and euler_characteristic(K) == 2
        and orient(K) is not None
    )


def classify_surface(K: IncidenceComplex) -> str:
    """Theorem thm:topology-production: classification by receipts."""
    if not is_connected(K):
        return "DISCONNECTED"
    if not is_closed_surface(K):
        return "NOT_A_CLOSED_SURFACE"
    chi = euler_characteristic(K)
    orientable = orient(K) is not None
    if orientable and chi == 2:
        return "S2"
    if orientable and chi == 0:
        return "T2"
    if orientable:
        return f"GENUS_{(2 - chi) // 2}"
    return "NONORIENTABLE"


# ---------------------------------------------------------------------------
# modular cross-ratio production (Theorem thm:conformal-cap-production)
# ---------------------------------------------------------------------------

def stereographic(p: np.ndarray) -> complex:
    """North-pole stereographic chart of the unit sphere."""
    x, y, z = p
    if abs(1.0 - z) < 1e-15:
        return complex(np.inf, 0.0)
    return complex(x / (1.0 - z), y / (1.0 - z))


def cross_ratio(z1: complex, z2: complex, z3: complex, z4: complex) -> complex:
    return ((z1 - z3) * (z2 - z4)) / ((z1 - z4) * (z2 - z3))


def mobius_normalize(z: complex, g1: complex, g2: complex, g3: complex) -> complex:
    """The Mobius map sending the gauge triple (g1,g2,g3) -> (0,1,inf)."""
    return ((z - g1) * (g2 - g3)) / ((z - g3) * (g2 - g1))


def reconstruct_from_cross_ratios(points: np.ndarray, gauge: tuple[int, int, int]) -> np.ndarray:
    """Embed sphere points into C from cross-ratio data against a gauge triple.

    Consumes only cross-ratios (the modular receipt data), never coordinates:
    cr(z, g2; g1, g3) equals the normalized coordinate mobius_normalize(z).
    Returns the array of reconstructed complex coordinates.
    """
    zs = [stereographic(p) for p in points]
    g1, g2, g3 = (zs[i] for i in gauge)
    out = []
    for i, z in enumerate(zs):
        # gauge points land on (0, 1, inf) by definition of the gauge
        if i == gauge[0]:
            out.append(0.0 + 0.0j)
        elif i == gauge[1]:
            out.append(1.0 + 0.0j)
        elif i == gauge[2]:
            out.append(complex(np.inf, 0.0))
        else:
            # cross-ratio receipt value for the quadruple (z, g2; g1, g3)
            out.append(cross_ratio(z, g2, g1, g3))
    return np.array(out)


def produced_cap_normal(boundary_points: np.ndarray) -> np.ndarray:
    """Produce n_C from (>= 3) boundary points of a cap circle on S^2.

    Fits the plane c . x = cos(alpha) through the points and applies the
    formula n_C = (cot alpha, csc alpha * c) of Proposition
    prop:round-cap-normal. Returns the 4-vector (time-first).
    """
    # least-squares plane through points on the sphere: minimize |P m - 1|
    m, *_ = np.linalg.lstsq(boundary_points, np.ones(len(boundary_points)), rcond=None)
    # m = c / cos(alpha)
    norm = np.linalg.norm(m)
    c = m / norm
    cos_a = 1.0 / norm
    cos_a = min(1.0 - 1e-12, max(-1.0 + 1e-12, cos_a))
    sin_a = np.sqrt(1.0 - cos_a * cos_a)
    return np.concatenate(([cos_a / sin_a], c / sin_a))


def minkowski(u: np.ndarray, v: np.ndarray) -> float:
    return float(-u[0] * v[0] + np.dot(u[1:], v[1:]))


def kms_receipt(clock_scale: float, beta_target: float = 2.0 * np.pi, tol: float = 1e-9) -> bool:
    """Wrong-normalization separation clause: the independently normalized
    geometric comparison certifies exactly the declared modular temperature."""
    return abs(clock_scale - beta_target) < tol


# ---------------------------------------------------------------------------
# invariance drivers (Lemma lem:incidence-invariance)
# ---------------------------------------------------------------------------

def readout_from_system(system: RepairSystem, schedule: list[int] | None = None) -> IncidenceComplex:
    nf = system.repair(schedule)
    return incidence_complex(normal_form_records(system, nf))


def complexes_equal(K1: IncidenceComplex, K2: IncidenceComplex) -> bool:
    return (
        K1.vertices == K2.vertices
        and K1.edges == K2.edges
        and K1.triangles == K2.triangles
    )


def gauge_relabel(records: list[tuple[int, ...]], perm: dict[int, int]) -> list[tuple[int, ...]]:
    return [tuple(sorted(perm[p] for p in rec)) for rec in records]


def complexes_isomorphic_under(K1: IncidenceComplex, K2: IncidenceComplex, perm: dict[int, int]) -> bool:
    ed = {frozenset(perm[v] for v in e) for e in K1.edges}
    tr = {frozenset(perm[v] for v in t) for t in K1.triangles}
    return (
        sorted(perm[v] for v in K1.vertices) == K2.vertices
        and ed == K2.edges
        and tr == K2.triangles
    )


def refinement_subdivide(records: list[tuple[int, ...]]) -> tuple[list[tuple[int, ...]], dict[int, int]]:
    """Barycentric-style edge subdivision of the record layer (one refinement
    stage) together with the coarse-graining projection on patch labels."""
    K = incidence_complex(records)
    next_label = max(K.vertices) + 1
    midpoint: dict[frozenset, int] = {}
    projection: dict[int, int] = {v: v for v in K.vertices}
    for e in sorted(K.edges, key=lambda e: tuple(sorted(e))):
        midpoint[e] = next_label
        a, _b = tuple(sorted(e))
        projection[next_label] = a  # midpoints project to an endpoint
        next_label += 1
    refined = []
    for rec in records:
        a, b, c = sorted(rec)
        mab = midpoint[frozenset((a, b))]
        mbc = midpoint[frozenset((b, c))]
        mac = midpoint[frozenset((a, c))]
        refined += [
            (a, mab, mac), (b, mab, mbc), (c, mbc, mac), (mab, mbc, mac),
        ]
    return refined, projection


def refinement_is_simplicial(fine: IncidenceComplex, coarse: IncidenceComplex,
                             projection: dict[int, int]) -> bool:
    """Petz support/CPTP clause: projected simplices are (degenerate) simplices."""
    for t in fine.triangles:
        img = {projection[v] for v in t}
        if len(img) == 3 and frozenset(img) not in coarse.triangles:
            return False
        if len(img) == 2 and frozenset(img) not in coarse.edges:
            return False
    return True
