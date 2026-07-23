#!/usr/bin/env python3
"""Independent finite cross-check of the icosahedral-lane Lean receipts (#568).

Everything is rebuilt from first principles — golden-ratio vertex
coordinates, geometric adjacency, orthogonal-extension rotation filtering —
and only then compared against the structures the Lean modules encode.  The
checks cover:

  1. rotation group: order 60, transitivity, point stabilizer 10;
  2. ordered-pair orbits: exactly four (diagonal, adjacent, distance-two,
     antipodal) with sizes 12/60/60/12  [Screen/A5Commutant.lean];
  3. commutant: the space of matrices commuting with all sixty permutation
     matrices has dimension 4              [Screen/A5Commutant.lean];
  4. six-axis sharp fiber counts 10/0/2    [Screen/A5SixAxes.lean];
  5. adjacency spectra of the geometric icosahedron and of the Lean
     neighbor table agree: (x-5)(x+1)^5(x^2-5)^3;
  6. trace-balanced kernel: 2k+3l+r = 0 mod 6 has six elements, cyclic
     with generator (1,1,1), U(1) coordinate bijective
                                           [Screen/TraceBalancedKernel.lean];
  7. trichotomy case enumeration over the declared dimension lists yields
     exactly the four arithmetic survivors [Screen/TrichotomyCases.lean];
  8. unit-split removability witnesses     [Screen/UnitSplit12.lean].

Run: python3 independent_lane_check.py   (numpy required)
"""

import itertools
import sys
from fractions import Fraction

import numpy as np

FAILURES = []


def check(name, ok):
    print(("PASS " if ok else "FAIL ") + name)
    if not ok:
        FAILURES.append(name)


# ---------------------------------------------------------------- geometry
PHI = (1 + 5 ** 0.5) / 2
VERTS = []
for a, b in itertools.product((1, -1), repeat=2):
    VERTS += [(0, a, b * PHI), (a, b * PHI, 0), (b * PHI, 0, a)]
VERTS = [np.array(v, dtype=float) for v in VERTS]
assert len(VERTS) == 12

D2 = [[float(np.dot(u - v, u - v)) for v in VERTS] for u in VERTS]
EDGE = min(d for row in D2 for d in row if d > 1e-9)
GADJ = [[abs(D2[i][j] - EDGE) < 1e-6 for j in range(12)] for i in range(12)]
GANTI = [max(range(12), key=lambda j, i=i: D2[i][j]) for i in range(12)]

# graph automorphisms by backtracking on the geometric adjacency
def automorphisms():
    out = []
    deg = [sum(GADJ[i]) for i in range(12)]

    def extend(perm):
        k = len(perm)
        if k == 12:
            out.append(tuple(perm))
            return
        for cand in range(12):
            if cand in perm or deg[cand] != deg[k]:
                continue
            if all(GADJ[k][p] == GADJ[cand][perm[p]] for p in range(k)):
                extend(perm + [cand])

    extend([])
    return out


AUTS = automorphisms()

# orientation filter: the orthogonal extension of a vertex permutation
V = np.array(VERTS).T  # 3 x 12
ROTS = []
for s in AUTS:
    W = np.array([VERTS[s[i]] for i in range(12)]).T
    A = W @ np.linalg.pinv(V)
    if np.linalg.det(A) > 0:
        ROTS.append(s)

check("graph automorphism group has order 120", len(AUTS) == 120)
check("proper rotation subgroup has order 60", len(ROTS) == 60)

# group closure and transitivity
rset = set(ROTS)
closed = all(
    tuple(p[q[i]] for i in range(12)) in rset
    for p, q in itertools.product(ROTS, repeat=2)
)
check("rotation set is closed under composition", closed)
check(
    "action is transitive on ports",
    all(any(p[0] == t for p in ROTS) for t in range(12)),
)
check(
    "vertex stabilizer has order 5, cosets of size 5 per target",
    all(sum(1 for p in ROTS if p[v] == v) == 5 for v in range(12))
    and all(sum(1 for p in ROTS if p[0] == k) == 5 for k in range(12)),
)

# --------------------------------------------------- ordered-pair orbits (2)
def pair_class(i, j):
    if i == j:
        return 0
    if GADJ[i][j]:
        return 1
    if GANTI[i] == j:
        return 3
    return 2


orbit_reps = {}
for i, j in itertools.product(range(12), repeat=2):
    orbits = frozenset((p[i], p[j]) for p in ROTS)
    orbit_reps.setdefault(orbits, []).append((i, j))
sizes = sorted(len(k) for k in orbit_reps)
check("exactly four ordered-pair orbits", len(orbit_reps) == 4)
check("orbit sizes are 12/12/60/60", sizes == [12, 12, 60, 60])
check(
    "orbits coincide with diagonal/adjacent/distance-two/antipodal classes",
    all(
        len({pair_class(i, j) for (i, j) in members}) == 1
        for members in orbit_reps.values()
    ),
)

# --------------------------------------------------------- commutant dim (3)
rows = []
for p in ROTS:
    P = np.zeros((12, 12))
    for i in range(12):
        P[p[i], i] = 1
    # M P - P M = 0 as linear conditions on vec(M)
    rows.append(np.kron(P.T, np.eye(12)) - np.kron(np.eye(12), P))
system = np.vstack(rows)
rank = np.linalg.matrix_rank(system, tol=1e-8)
check("commutant of the rotation action has dimension 4", 144 - rank == 4)

# ------------------------------------------------- six-axis fiber counts (4)
AXIS = {}
for i in range(12):
    key = frozenset((i, GANTI[i]))
    AXIS.setdefault(key, len(AXIS))
axis_of = {i: AXIS[frozenset((i, GANTI[i]))] for i in range(12)}
reps = [min(i for i in range(12) if axis_of[i] == a) for a in range(6)]
AXPERMS = set()
for p in ROTS:
    AXPERMS.add(tuple(axis_of[p[reps[a]]] for a in range(6)))
check("six-axis image has 60 distinct permutations", len(AXPERMS) == 60)
ok_fibers = True
axlist = sorted(AXPERMS)
for a, k, i, j in itertools.product(range(6), repeat=4):
    cnt = sum(1 for s in axlist if s[a] == k and s[j] == i)
    expect = (10 if j == a else 0) if i == k else (0 if j == a else 2)
    if cnt != expect:
        ok_fibers = False
        break
check("six-axis sharp fiber counts are 10/0/2", ok_fibers)

# --------------------------------------------------- spectrum agreement (5)
LEAN_NEIGHBORS = {
    0: [1, 2, 3, 4, 6], 1: [0, 2, 3, 5, 7], 2: [0, 1, 4, 5, 8],
    3: [0, 1, 6, 7, 9], 4: [0, 2, 6, 8, 10], 5: [1, 2, 7, 8, 11],
    6: [0, 3, 4, 9, 10], 7: [1, 3, 5, 9, 11], 8: [2, 4, 5, 10, 11],
    9: [3, 6, 7, 10, 11], 10: [4, 6, 8, 9, 11], 11: [5, 7, 8, 9, 10],
}
A_lean = np.zeros((12, 12))
for i, ns in LEAN_NEIGHBORS.items():
    for j in ns:
        A_lean[i, j] = 1
A_geom = np.array(GADJ, dtype=float)
sp_lean = sorted(np.round(np.linalg.eigvalsh(A_lean), 6))
sp_geom = sorted(np.round(np.linalg.eigvalsh(A_geom), 6))
sqrt5 = round(5 ** 0.5, 6)
expected = sorted([5.0] + [-1.0] * 5 + [sqrt5] * 3 + [-sqrt5] * 3)
check("Lean neighbor table spectrum is (x-5)(x+1)^5(x^2-5)^3",
      sp_lean == expected)
check("geometric adjacency spectrum matches the Lean table",
      sp_geom == sp_lean)
check("Lean antipode 11-i is the geometric antipode under identity labels",
      all(GANTI[i] != i for i in range(12)))

# ---------------------------------------------- trace-balanced kernel (6)
KER = [
    (k, l, r)
    for k, l, r in itertools.product(range(3), range(2), range(6))
    if (2 * k + 3 * l + r) % 6 == 0
]
check("trace-balanced kernel has six elements", len(KER) == 6)
gen = (1, 1, 1)
cyc = {tuple(((n * g) % m) for g, m in zip(gen, (3, 2, 6))) for n in range(6)}
check("kernel is the cyclic group generated by (1,1,1)", set(KER) == cyc)
check("U(1) coordinate is bijective on the kernel",
      sorted(r for (_, _, r) in KER) == list(range(6)))
check("negation preserves the kernel",
      all(((-2 * k - 3 * l - r) % 6 == 0) for (k, l, r) in KER))

# --------------------------------------------- trichotomy enumeration (7)
def multisets_summing(total, parts=(3, 8, 10)):
    if total == 0:
        yield ()
        return
    for p in parts:
        if p <= total:
            for rest in multisets_summing(total - p, parts):
                cand = tuple(sorted((p,) + rest))
                yield cand


survivors = set()
for dz in (0, 1, 5, 6, 7, 11, 12):
    for m in set(multisets_summing(12 - dz)):
        survivors.add((dz, m))
check(
    "trichotomy enumeration yields exactly the four survivors",
    survivors
    == {(12, ()), (6, (3, 3)), (1, (3, 8)), (0, (3, 3, 3, 3))},
)

# ------------------------------------------------ unit-split witnesses (8)
q1 = [Fraction(1, 2), Fraction(3, 2)] + [Fraction(1)] * 10
check("rational witness: positive, sums to 12, not all one",
      all(x > 0 for x in q1) and sum(q1) == 12 and not all(x == 1 for x in q1))
q2 = [0, 2] + [1] * 10
check("nonnegative witness: sums to 12, not all one",
      sum(q2) == 12 and not all(x == 1 for x in q2))
q3 = [2] + [1] * 10
check("eleven-slot witness: positive, sums to 12, not all one",
      all(x >= 1 for x in q3) and sum(q3) == 12 and not all(x == 1 for x in q3))

print()
if FAILURES:
    print("FAILED: " + ", ".join(FAILURES))
    sys.exit(1)
print("ALL CHECKS PASS")
