#!/usr/bin/env python3
"""Exact compact-real-locus classification for the issue #705 fourteen-parameter Jacobi variety.

Semantic inputs are the pinned stage-one Reynolds basis and the pinned stage-two
Jacobi reduction for issue #566. This producer

1. reproduces the complete Jacobi system from the stage-one basis and verifies
   the coefficient-row rank 38 and the exact 11+27 rowspace split against the
   pinned stage-two data;
2. rewrites the system in the pinned channel coordinates and computes the exact
   quadric spans of the system restricted to every channel subspace that the
   compact-type structure theory admits;
3. decides compactness through the Killing-form criterion, with every Killing
   matrix, signature, derived algebra, and center computed exactly over
   Q(sqrt(5));
4. emits a deterministic JSON certificate together with an independent verifier.

Classification statement certified here, conditional on the three named
textbook lemmas recorded in the certificate: the compact real locus of the
fourteen-parameter variety is the union of exactly three families, all inside
the derivation-free slice d_plus = d_minus = d_five = 0:

  P  the closed two-parameter plane spanned by t_pp_to_p and t_mm_to_m;
  F  the three-parameter cell t_pp_to_p = a, t_ff_to_p = b, t_pf_to_f =
     -sqrt(5)*a, t_mm_to_m = e with a nonzero and a*b < 0 (su(3)+so(3) type);
  G  the mirror cell t_mm_to_m = a, t_ff_to_m = b, t_mf_to_f = sqrt(5)*a,
     t_pp_to_p = e with a nonzero and a*b > 0.

The producer does not select a source, does not build port perturbations or
skew generators, does not construct recharting implementers or holonomy, and
does not close issue #705.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
STAGE1 = REPO / "code/a5_closure/issue_566_bracket_space_stage1"
STAGE2 = REPO / "code/a5_closure/issue_566_bracket_space_stage2"
STAGE1_BASIS = STAGE1 / "a5_alternating_bracket_reynolds_basis.json"
STAGE1_RECEIPT = STAGE1 / "a5_alternating_bracket_space_stage1.receipt.json"
STAGE2_SYSTEM = STAGE2 / "a5_jacobi_system_reduction.json"
STAGE2_RECEIPT = STAGE2 / "a5_jacobi_stage2.receipt.json"
CERTIFICATE_PATH = HERE / "b14_compact_locus.certificate.json"
VERIFIER_PATH = HERE / "verify.py"

PORTS = 12
PARAMS = 14
MONOMIALS = [(a, b) for a in range(PARAMS) for b in range(a, PARAMS)]
MIDX = {m: i for i, m in enumerate(MONOMIALS)}
SECTORS = ["fixed", "three_plus", "three_minus", "five"]
SECTOR_DIMS = {"fixed": 1, "three_plus": 3, "three_minus": 3, "five": 5}
CHANNEL_IDS = [
    "d_plus", "d_minus", "d_five", "t_pp_to_p", "t_mm_to_m", "t_ff_to_p",
    "t_ff_to_m", "t_pm_to_f", "t_pf_to_p", "t_pf_to_m", "t_pf_to_f",
    "t_mf_to_p", "t_mf_to_m", "t_mf_to_f",
]


class ClassifyError(RuntimeError):
    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise ClassifyError(code, message)


# ---------------------------------------------------------------- Q(sqrt(5))
class Q5:
    """Exact element a + b*sqrt(5) with rational a, b."""

    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = a if isinstance(a, Fraction) else Fraction(a)
        self.b = b if isinstance(b, Fraction) else Fraction(b)

    def __add__(self, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return Q5(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return Q5(-self.a, -self.b)

    def __sub__(self, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return Q5(self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        return Q5(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return Q5(self.a * o.a + 5 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def inv(self):
        n = self.a * self.a - 5 * self.b * self.b
        require(n != 0, "Q5_DIVISION", "division by zero in Q(sqrt(5))")
        return Q5(self.a / n, -self.b / n)

    def __truediv__(self, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return self * o.inv()

    def __bool__(self):
        return bool(self.a or self.b)

    def __eq__(self, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return self.a == o.a and self.b == o.b

    def __hash__(self):
        return hash((self.a, self.b))

    def sign(self) -> int:
        if not self:
            return 0
        if self.b == 0:
            return 1 if self.a > 0 else -1
        if self.a == 0:
            return 1 if self.b > 0 else -1
        if self.a * self.a > 5 * self.b * self.b:
            return 1 if self.a > 0 else -1
        return 1 if self.b > 0 else -1


Q5_ZERO = Q5()
Q5_ONE = Q5(1)
R5 = Q5(0, 1)


def q5_encode(v: Q5) -> list[int]:
    return [v.a.numerator, v.a.denominator, v.b.numerator, v.b.denominator]


def q5_decode(entry) -> Q5:
    return Q5(Fraction(entry[0], entry[1]), Fraction(entry[2], entry[3]))


def vec_encode(vec) -> list[list[int]]:
    return [[i, *q5_encode(v)] for i, v in enumerate(vec) if v]


def sym_encode(M) -> list[list[int]]:
    return [[i, j, *q5_encode(M[i][j])] for i in range(len(M)) for j in range(i, len(M)) if M[i][j]]


def rref_q5(rows):
    work = [list(r) for r in rows if any(r)]
    if not work:
        return [], []
    pivots = []
    pr = 0
    for col in range(len(work[0])):
        p = next((r for r in range(pr, len(work)) if work[r][col]), None)
        if p is None:
            continue
        work[pr], work[p] = work[p], work[pr]
        sc = work[pr][col].inv()
        work[pr] = [v * sc for v in work[pr]]
        for r in range(len(work)):
            if r != pr and work[r][col]:
                f = work[r][col]
                work[r] = [x - f * y for x, y in zip(work[r], work[pr])]
        pivots.append(col)
        pr += 1
        if pr == len(work):
            break
    return work[:pr], pivots


def rank_q5(rows) -> int:
    return len(rref_q5(rows)[0])


def in_rowspace(rref_rows, pivots, row) -> bool:
    residual = list(row)
    for base, col in zip(rref_rows, pivots):
        if residual[col]:
            f = residual[col]
            residual = [x - f * y for x, y in zip(residual, base)]
    return not any(residual)


def rref_fraction(rows):
    work = [list(r) for r in rows if any(r)]
    if not work:
        return [], []
    pivots = []
    pr = 0
    for col in range(len(work[0])):
        p = next((r for r in range(pr, len(work)) if work[r][col]), None)
        if p is None:
            continue
        work[pr], work[p] = work[p], work[pr]
        sc = work[pr][col]
        work[pr] = [v / sc for v in work[pr]]
        for r in range(len(work)):
            if r != pr and work[r][col]:
                f = work[r][col]
                work[r] = [x - f * y for x, y in zip(work[r], work[pr])]
        pivots.append(col)
        pr += 1
        if pr == len(work):
            break
    return work[:pr], pivots


def signature_q5(M):
    """Exact inertia (positive, negative, zero) of a symmetric Q5 matrix."""
    n = len(M)
    W = [row[:] for row in M]
    pos = neg = 0
    size = n
    while size:
        d = next((i for i in range(size) if W[i][i]), None)
        if d is None:
            found = None
            for i in range(size):
                for j in range(i + 1, size):
                    if W[i][j]:
                        found = (i, j)
                        break
                if found:
                    break
            if not found:
                break
            i, j = found
            for k in range(size):
                W[i][k] = W[i][k] + W[j][k]
            for k in range(size):
                W[k][i] = W[k][i] + W[k][j]
            d = i
        if d != 0:
            W[0], W[d] = W[d], W[0]
            for row in W:
                row[0], row[d] = row[d], row[0]
        piv = W[0][0]
        if piv.sign() > 0:
            pos += 1
        else:
            neg += 1
        Wn = []
        for i in range(1, size):
            f = W[i][0] / piv
            Wn.append([W[i][j] - f * W[0][j] for j in range(1, size)])
        W = Wn
        size -= 1
    return [pos, neg, n - pos - neg]


def canonical_bytes(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def object_sha256(value) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------------------------------------------------------- load pinned
def load_pinned():
    basis_raw = json.loads(STAGE1_BASIS.read_text(encoding="utf-8"))
    stage1_receipt = json.loads(STAGE1_RECEIPT.read_text(encoding="utf-8"))
    system = json.loads(STAGE2_SYSTEM.read_text(encoding="utf-8"))
    stage2_receipt = json.loads(STAGE2_RECEIPT.read_text(encoding="utf-8"))
    up = system["upstream"]
    require(up["basis_file_sha256"] == file_sha256(STAGE1_BASIS), "PIN_BASIS", "stage-one basis hash differs")
    require(up["receipt_file_sha256"] == file_sha256(STAGE1_RECEIPT), "PIN_RECEIPT", "stage-one receipt hash differs")
    require(up["basis_canonical_json_sha256"] == object_sha256(basis_raw), "PIN_BASIS_CANONICAL", "basis canonical hash differs")
    require(
        stage2_receipt["system_artifact"]["canonical_json_sha256"] == object_sha256(system),
        "PIN_SYSTEM",
        "stage-two system canonical hash differs",
    )
    receipt_copy = dict(stage2_receipt)
    stored = receipt_copy.pop("receipt_sha256", None)
    require(stored == object_sha256(receipt_copy), "PIN_STAGE2_SELF", "stage-two receipt self-hash differs")
    return basis_raw, stage1_receipt, system, stage2_receipt


def structure_constants(basis_raw):
    C = [[[[Fraction(0)] * PORTS for _ in range(PORTS)] for _ in range(PORTS)] for _ in range(PARAMS)]
    rows = basis_raw["basis"]
    require(len(rows) == PARAMS, "BASIS_COUNT", "expected fourteen basis rows")
    for a, row in enumerate(rows):
        for o, l, r, n, d in row["entries"]:
            require(l < r, "BASIS_ENTRY", "basis entry is not upper triangular")
            C[a][o][l][r] = Fraction(n, d)
            C[a][o][r][l] = -Fraction(n, d)
    return C


# ------------------------------------------------------- Jacobi reproduction
def parity(vals) -> int:
    return sum(vals[i] > vals[j] for i in range(len(vals)) for j in range(i + 1, len(vals))) % 2


def raw_jacobi_rows(C):
    pair_support: dict[tuple[int, int], list] = {}
    col_index: dict[tuple[int, int], list] = {}
    for a in range(PARAMS):
        for o in range(PORTS):
            for l in range(PORTS):
                for r in range(l + 1, PORTS):
                    v = C[a][o][l][r]
                    if v:
                        pair_support.setdefault((l, r), []).append((a, o, v))
    for b in range(PARAMS):
        for o in range(PORTS):
            for m in range(PORTS):
                for k in range(PORTS):
                    v = C[b][o][m][k]
                    if v:
                        col_index.setdefault((m, k), []).append((b, o, v))

    def pair_val(i, j):
        if i < j:
            return pair_support.get((i, j), [])
        return [(a, m, -v) for (a, m, v) in pair_support.get((j, i), [])]

    raw = {}
    for (i, j, k) in itertools.combinations(range(PORTS), 3):
        rows = [[Fraction(0)] * len(MONOMIALS) for _ in range(PORTS)]
        for (pi, pj, pk) in ((i, j, k), (j, k, i), (k, i, j)):
            for (a, m, v) in pair_val(pi, pj):
                for (b, o, w) in col_index.get((m, pk), []):
                    rows[o][MIDX[tuple(sorted((a, b)))]] += v * w
        for o in range(PORTS):
            require(any(rows[o]), "JACOBI_ZERO_ROW", "unexpected identically zero coordinate equation")
            raw[(o, i, j, k)] = tuple(rows[o])
    require(len(raw) == 2640, "JACOBI_COUNT", "wrong coordinate equation count")
    return raw


def signed_orbits(group):
    unassigned = {
        (o, i, j, k)
        for o in range(PORTS)
        for (i, j, k) in itertools.combinations(range(PORTS), 3)
    }
    orbits = []
    while unassigned:
        seed = min(unassigned)
        signs = {}
        for g in group:
            o, i, j, k = seed
            images = [g[i], g[j], g[k]]
            sign = -1 if parity(images) else 1
            image = (g[o], *sorted(images))
            require(signs.get(image, sign) == sign, "ORBIT_SIGN", "signed orbit has a negative stabilizer")
            signs[image] = sign
        unassigned -= set(signs)
        orbits.append((seed, signs))
    require(len(orbits) == 44 and all(len(s) == 60 for _, s in orbits), "ORBIT_COUNT", "expected forty-four size-sixty orbits")
    return orbits


def decode_pinned_equation(entries):
    row = [Fraction(0)] * len(MONOMIALS)
    for a, b, n, d in entries:
        row[MIDX[(a, b)]] = Fraction(n, d)
    return tuple(row)


# ------------------------------------------------- channel transform rebuild
def matmul(X, Y):
    n, m, p = len(X), len(Y), len(Y[0])
    Z = [[Q5_ZERO] * p for _ in range(n)]
    for i in range(n):
        for k in range(m):
            if X[i][k]:
                for j in range(p):
                    if Y[k][j]:
                        Z[i][j] = Z[i][j] + X[i][k] * Y[k][j]
    return Z


def identity_matrix(n):
    return [[Q5_ONE if i == j else Q5_ZERO for j in range(n)] for i in range(n)]


def pair_orbits(group):
    unassigned = {(i, j) for i in range(PORTS) for j in range(PORTS)}
    orbits = []
    while unassigned:
        seed = min(unassigned)
        orbit = {(g[seed[0]], g[seed[1]]) for g in group}
        unassigned -= orbit
        orbits.append(orbit)
    return orbits


def spectral_projectors(group):
    orbits = pair_orbits(group)
    require(sorted(len(o) for o in orbits) == [12, 12, 60, 60], "PAIR_ORBITS", "unexpected ordered-pair orbit sizes")
    candidates = sorted((o for o in orbits if len(o) == 60), key=lambda o: tuple(sorted(o)))
    orbital = candidates[0]
    A = [[Q5_ZERO] * PORTS for _ in range(PORTS)]
    for i, j in orbital:
        A[i][j] = Q5_ONE
    eigenvalues = {"fixed": Q5(5), "three_plus": R5, "three_minus": -R5, "five": Q5(-1)}
    projectors = {}
    for label, ev in eigenvalues.items():
        M = identity_matrix(PORTS)
        for label2, ev2 in eigenvalues.items():
            if label2 == label:
                continue
            factor = [[A[i][j] - (ev2 if i == j else Q5_ZERO) for j in range(PORTS)] for i in range(PORTS)]
            M = [[v * (ev - ev2).inv() for v in row] for row in matmul(M, factor)]
        projectors[label] = M
    total = [[sum((projectors[l][i][j] for l in SECTORS), Q5_ZERO) for j in range(PORTS)] for i in range(PORTS)]
    require(total == identity_matrix(PORTS), "PROJECTOR_SUM", "spectral projectors do not sum to identity")
    for label in SECTORS:
        P = projectors[label]
        require(matmul(P, P) == P, "PROJECTOR_IDEMPOTENT", f"{label} projector is not idempotent")
        trace = sum((P[i][i] for i in range(PORTS)), Q5_ZERO)
        require(trace == Q5(SECTOR_DIMS[label]), "PROJECTOR_TRACE", f"{label} projector trace differs")
    return orbital, projectors, len(orbits)


def derivation_matrices(C):
    result = []
    for a in range(PARAMS):
        M = [[Q5(sum((C[a][o][s][c] for s in range(PORTS)), Fraction(0))) for c in range(PORTS)] for o in range(PORTS)]
        result.append(M)
    return result


CHANNEL_SPECS = [
    ("t_pp_to_p", "three_plus", "three_plus", "three_plus", (1, 0, 0)),
    ("t_mm_to_m", "three_minus", "three_minus", "three_minus", (0, 1, 0)),
    ("t_ff_to_p", "five", "five", "three_plus", (-1, 0, 2)),
    ("t_ff_to_m", "five", "five", "three_minus", (0, -1, 2)),
    ("t_pm_to_f", "three_plus", "three_minus", "five", (1, 1, -1)),
    ("t_pf_to_p", "three_plus", "five", "three_plus", (0, 0, 1)),
    ("t_pf_to_m", "three_plus", "five", "three_minus", (1, -1, 1)),
    ("t_pf_to_f", "three_plus", "five", "five", (1, 0, 0)),
    ("t_mf_to_p", "three_minus", "five", "three_plus", (-1, 1, 1)),
    ("t_mf_to_m", "three_minus", "five", "three_minus", (0, 0, 1)),
    ("t_mf_to_f", "three_minus", "five", "five", (0, 1, 0)),
]


def channel_coordinate_form(vectors, Pl, Pr, Po, coordinate):
    output, li, ri = coordinate
    form = []
    for vector in vectors:
        value = Q5_ZERO
        for (so, p, q), coefficient in vector.items():
            of = Po[output][so]
            if not of:
                continue
            inf = Pl[p][li] * Pr[q][ri] - Pl[q][li] * Pr[p][ri]
            if inf:
                value = value + of * inf * coefficient
        form.append(value)
    return form


def channel_transform(vectors, C, group):
    orbital, projectors, orbit_count = spectral_projectors(group)
    derivations = derivation_matrices(C)
    eigenforms = {}
    for label in SECTORS:
        P = projectors[label]
        form = []
        for M in derivations:
            trace = Q5_ZERO
            for i in range(PORTS):
                for j in range(PORTS):
                    if P[i][j] and M[j][i]:
                        trace = trace + P[i][j] * M[j][i]
            form.append(trace / SECTOR_DIMS[label])
        for a, M in enumerate(derivations):
            left = matmul(M, P)
            right = [[form[a] * v for v in row] for row in P]
            require(left == right, "DERIVATION_SPECTRAL", f"fixed-line adjoint is not scalar on {label}")
        eigenforms[label] = form
    require(not any(eigenforms["fixed"]), "DERIVATION_FIXED", "fixed-line adjoint does not kill the fixed line")
    transform = [eigenforms["three_plus"], eigenforms["three_minus"], eigenforms["five"]]
    for _, left, right, output, _ in CHANNEL_SPECS:
        found = None
        for coordinate in itertools.product(range(PORTS), repeat=3):
            form = channel_coordinate_form(vectors, projectors[left], projectors[right], projectors[output], coordinate)
            if any(form):
                found = form
                break
        require(found is not None, "CHANNEL_EMPTY", "channel has no nonzero coordinate")
        transform.append(found)
    rr, _ = rref_q5(transform)
    require(len(rr) == PARAMS, "TRANSFORM_RANK", "channel transform is singular")
    aug = [list(row) + [Q5_ONE if i == j else Q5_ZERO for j in range(PARAMS)] for i, row in enumerate(transform)]
    solved, pivots = rref_q5(aug)
    require(pivots == list(range(PARAMS)), "TRANSFORM_INVERT", "channel transform inversion failed")
    inverse = [row[PARAMS:] for row in solved]
    require(matmul(transform, inverse) == identity_matrix(PARAMS), "TRANSFORM_INVERSE", "inverse is not exact")
    return orbital, projectors, orbit_count, derivations, transform, inverse


# ------------------------------------------------------------ restricted spans
def restricted_span_rows(reps, inverse, support):
    """Rewrite each representative x-quadric on the channel subspace spanned by
    `support` (list of channel indices); return rows over the restricted
    channel-pair monomial basis."""
    pairs = [(support[i], support[j]) for i in range(len(support)) for j in range(i, len(support))]
    pidx = {p: i for i, p in enumerate(pairs)}
    rows = []
    for row in reps:
        out = [Q5_ZERO] * len(pairs)
        for (alpha, beta), coefficient in zip(MONOMIALS, row):
            if not coefficient:
                continue
            for ci in support:
                va = inverse[alpha][ci]
                if not va:
                    continue
                for cj in support:
                    vb = inverse[beta][cj]
                    if vb:
                        out[pidx[tuple(sorted((ci, cj)))]] = out[pidx[tuple(sorted((ci, cj)))]] + coefficient * va * vb
        rows.append(out)
    return pairs, rows


# --------------------------------------------------------------- bracket tools
def bracket_tensor(C, inverse, channel_values):
    x = [sum((inverse[i][j] * v for j, v in channel_values.items()), Q5_ZERO) for i in range(PARAMS)]
    B = [[[Q5_ZERO] * PORTS for _ in range(PORTS)] for _ in range(PORTS)]
    for a in range(PARAMS):
        if not x[a]:
            continue
        for o in range(PORTS):
            for i in range(PORTS):
                row = C[a][o][i]
                for j in range(PORTS):
                    if row[j]:
                        B[o][i][j] = B[o][i][j] + x[a] * row[j]
    return x, B


def jacobi_defect(B):
    worst = Q5_ZERO
    for (i, j, k) in itertools.combinations(range(PORTS), 3):
        for o in range(PORTS):
            v = Q5_ZERO
            for m in range(PORTS):
                v = v + B[m][i][j] * B[o][m][k] + B[m][j][k] * B[o][m][i] + B[m][k][i] * B[o][m][j]
            if v:
                return v
    return worst


def killing_matrix(B):
    K = [[Q5_ZERO] * PORTS for _ in range(PORTS)]
    for i in range(PORTS):
        for j in range(i, PORTS):
            v = Q5_ZERO
            for m in range(PORTS):
                for o in range(PORTS):
                    if B[m][i][o] and B[o][j][m]:
                        v = v + B[m][i][o] * B[o][j][m]
            K[i][j] = v
            K[j][i] = v
    return K


def derived_span(B):
    rows = []
    for i in range(PORTS):
        for j in range(i + 1, PORTS):
            col = [B[o][i][j] for o in range(PORTS)]
            if any(col):
                rows.append(col)
    return rref_q5(rows)[0] if rows else []


def center_basis(B):
    rows = []
    for v in range(PORTS):
        flat = []
        for o in range(PORTS):
            for j in range(PORTS):
                flat.append(B[o][v][j])
        rows.append(flat)
    kernel = []
    work, piv = rref_q5([list(r) for r in transpose_matrix(rows)])
    pivset = set(piv)
    free = [c for c in range(PORTS) if c not in pivset]
    for fc in free:
        vec = [Q5_ZERO] * PORTS
        vec[fc] = Q5_ONE
        for row, pc in zip(work, piv):
            vec[pc] = -row[fc]
        kernel.append(vec)
    for vec in kernel:
        for o in range(PORTS):
            for j in range(PORTS):
                s = sum((vec[v] * B[o][v][j] for v in range(PORTS) if vec[v]), Q5_ZERO)
                require(not s, "CENTER_KERNEL", "kernel vector fails the centrality recheck")
    return kernel


def transpose_matrix(M):
    return [[M[i][j] for i in range(len(M))] for j in range(len(M[0]))]


def sector_membership(rref_rows, projectors):
    labels = []
    for label in SECTORS:
        P = projectors[label]
        for row in rref_rows:
            image = [sum((P[o][i] * row[i] for i in range(PORTS) if row[i]), Q5_ZERO) for o in range(PORTS)]
            if any(image):
                labels.append(label)
                break
    return labels


def restrict_form(K, basis_rows):
    return [
        [
            sum(
                (basis_rows[r][i] * K[i][j] * basis_rows[s][j] for i in range(PORTS) for j in range(PORTS) if basis_rows[r][i] and K[i][j] and basis_rows[s][j]),
                Q5_ZERO,
            )
            for s in range(len(basis_rows))
        ]
        for r in range(len(basis_rows))
    ]


# ------------------------------------------------------------------ families
def family_channel_values(family_id, params):
    a, b, e = (list(params) + [0, 0, 0])[:3]
    a, b, e = Q5(a) if not isinstance(a, Q5) else a, Q5(b) if not isinstance(b, Q5) else b, Q5(e) if not isinstance(e, Q5) else e
    if family_id == "PLANE":
        return {3: a, 4: b}
    if family_id == "SU3PF":
        return {3: a, 5: b, 10: Q5(0, -1) * a, 4: e}
    if family_id == "SU3MF":
        return {4: a, 6: b, 13: Q5(0, 1) * a, 3: e}
    if family_id == "DERIV":
        return {0: a, 1: b, 2: e}
    raise ClassifyError("FAMILY_ID", f"unknown family {family_id}")


FAMILY_PARAM_COUNT = {"PLANE": 2, "SU3PF": 3, "SU3MF": 3, "DERIV": 3}


def spanning_points(n):
    points = []
    for i in range(n):
        e = [0] * n
        e[i] = 1
        points.append(tuple(e))
    for i in range(n):
        for j in range(i + 1, n):
            e = [0] * n
            e[i] = e[j] = 1
            points.append(tuple(e))
    return points


def evaluate_rows_at(reps, x):
    values = []
    for row in reps:
        v = Q5_ZERO
        for (alpha, beta), coefficient in zip(MONOMIALS, row):
            if coefficient and x[alpha] and x[beta]:
                v = v + coefficient * x[alpha] * x[beta]
        values.append(v)
    return values


def killing_monomial_matrices(C, inverse, family_id):
    n = FAMILY_PARAM_COUNT[family_id]
    K_at = {}
    for point in spanning_points(n):
        _, B = bracket_tensor(C, inverse, family_channel_values(family_id, point))
        K_at[point] = killing_matrix(B)
    M = {}
    for i in range(n):
        e = [0] * n
        e[i] = 1
        M[(i, i)] = K_at[tuple(e)]
    for i in range(n):
        for j in range(i + 1, n):
            e = [0] * n
            e[i] = e[j] = 1
            Kij = K_at[tuple(e)]
            M[(i, j)] = [
                [Kij[r][c] - M[(i, i)][r][c] - M[(j, j)][r][c] for c in range(PORTS)]
                for r in range(PORTS)
            ]
    return {m: mat for m, mat in M.items() if any(any(row) for row in mat)}


# ----------------------------------------------------------------------- main
def main() -> dict:
    basis_raw, stage1_receipt, system, stage2_receipt = load_pinned()
    C = structure_constants(basis_raw)
    vectors = []
    for row in basis_raw["basis"]:
        vectors.append({(o, l, r): Fraction(n, d) for o, l, r, n, d in row["entries"]})
    group = [tuple(r) for r in stage1_receipt["proper_port_action"]["permutation_rows"]]
    require(len(group) == 60 and len(set(group)) == 60, "GROUP", "expected sixty distinct group rows")

    # (1) Jacobi reproduction and rank/split verification.
    raw = raw_jacobi_rows(C)
    orbits = signed_orbits(group)
    orbit_constancy = all(
        raw[coord] == tuple(s * v for v in raw[seed])
        for seed, signs in orbits
        for coord, s in signs.items()
    )
    require(orbit_constancy, "ORBIT_CONSTANCY", "coordinate equations are not sign-constant on orbits")
    reps = [raw[seed] for seed, _ in orbits]
    rep_rref, rep_pivots = rref_fraction(reps)
    require(len(rep_rref) == 38, "RANK38", "representative rank is not 38")
    jt = system["jacobi_tensor"]
    require(rep_pivots == jt["representative_rref_pivot_monomial_indices"], "PIVOTS", "pivot monomials differ from pinned")
    require(jt["representative_coefficient_rank_over_Q"] == 38, "PINNED_RANK", "pinned rank differs")
    pinned_by_coord = {tuple(p["target_coordinate"]): decode_pinned_equation(p["equation"]) for p in jt["representative_equations"]}
    matches = sum(1 for seed, _ in orbits if pinned_by_coord.get(seed) == raw[seed])
    require(matches == 44, "PINNED_REPS", "pinned representative equations differ")

    # fixed-line contraction: rank 11.
    contraction_rows = []
    for o in range(PORTS):
        for i, j in itertools.combinations(range(PORTS), 2):
            total = [Fraction(0)] * len(MONOMIALS)
            for s in range(PORTS):
                vals = sorted((s, i, j))
                if len(set(vals)) < 3:
                    continue
                sign = -1 if parity((s, i, j)) else 1
                row = raw[(o, *vals)]
                for idx, v in enumerate(row):
                    if v:
                        total[idx] += sign * v
            require(any(total), "CONTRACTION_ZERO", "unexpected zero contracted row")
            contraction_rows.append(total)
    require(len(contraction_rows) == 792, "CONTRACTION_COUNT", "wrong contracted row count")
    con_rref, con_pivots = rref_fraction(contraction_rows)
    require(len(con_rref) == 11, "RANK11", "contracted rank is not 11")
    fl = system["fixed_line_reduction"]
    require(fl["contracted_coefficient_rank_over_Q"] == 11, "PINNED_RANK11", "pinned contracted rank differs")
    pinned_contracted = []
    for entries in fl["contracted_integer_equations_in_x"]:
        row = [Fraction(0)] * len(MONOMIALS)
        for a, b, v in entries:
            row[MIDX[(a, b)]] = Fraction(v)
        pinned_contracted.append(row)
    require(len(pinned_contracted) == 11, "PINNED_CONTRACTED", "expected eleven pinned contracted rows")
    for row in pinned_contracted:
        residual = list(row)
        for base, col in zip(con_rref, con_pivots):
            if residual[col]:
                f = residual[col]
                residual = [x - f * y for x, y in zip(residual, base)]
        require(not any(residual), "PINNED_CONTRACTED_SPAN", "pinned contracted row leaves the contracted rowspace")
    require(len(rref_fraction(pinned_contracted)[0]) == 11, "PINNED_CONTRACTED_RANK", "pinned contracted rows are dependent")

    # channel transform reproduction.
    orbital, projectors, pair_orbit_count, derivations, transform, inverse = channel_transform(vectors, C, group)
    cdm = fl["channel_decomposition"]
    pinned_transform = [[Q5_ZERO] * PARAMS for _ in range(PARAMS)]
    for i, entries in enumerate(cdm["transform_rows"]):
        for j, n1, d1, n2, d2 in entries:
            pinned_transform[i][j] = Q5(Fraction(n1, d1), Fraction(n2, d2))
    pinned_inverse = [[Q5_ZERO] * PARAMS for _ in range(PARAMS)]
    for i, entries in enumerate(cdm["inverse_transform_rows"]):
        for j, n1, d1, n2, d2 in entries:
            pinned_inverse[i][j] = Q5(Fraction(n1, d1), Fraction(n2, d2))
    require(transform == pinned_transform, "TRANSFORM_MATCH", "recomputed transform differs from pinned")
    require(inverse == pinned_inverse, "INVERSE_MATCH", "recomputed inverse differs from pinned")
    pinned_orbital = [tuple(p) for p in cdm["spectral_decomposition"]["canonical_valency_five_orbital"]]
    require(sorted(orbital) == sorted(pinned_orbital), "ORBITAL_MATCH", "canonical orbital differs from pinned")

    # product equations and the 11+27 split over Q(sqrt(5)).
    product_rows = []
    for number, (_, _, _, _, weight) in enumerate(CHANNEL_SPECS):
        weight_form = [
            sum((Q5(w) * transform[t][col] for t, w in enumerate(weight) if w), Q5_ZERO)
            for col in range(PARAMS)
        ]
        channel_form = transform[3 + number]
        row = [Q5_ZERO] * len(MONOMIALS)
        for i, va in enumerate(weight_form):
            if not va:
                continue
            for j, vb in enumerate(channel_form):
                if vb:
                    row[MIDX[tuple(sorted((i, j)))]] = row[MIDX[tuple(sorted((i, j)))]] + va * vb
        product_rows.append(row)
    prod_rref, prod_pivots = rref_q5(product_rows)
    require(len(prod_rref) == 11, "PRODUCT_RANK", "product equations do not have rank eleven")
    con_q5 = [[Q5(v) for v in row] for row in con_rref]
    for row in con_q5:
        require(in_rowspace(prod_rref, prod_pivots, row), "SPLIT_PRODUCT", "contracted row leaves the product rowspace")
    con_q5_rref, con_q5_piv = rref_q5(con_q5)
    for row in prod_rref:
        require(in_rowspace(con_q5_rref, con_q5_piv, row), "SPLIT_PRODUCT_REV", "product row leaves the contracted rowspace")
    pinned_residual = []
    for entries in fl["residual_integer_equations_in_x"]:
        row = [Q5_ZERO] * len(MONOMIALS)
        for a, b, v in entries:
            row[MIDX[(a, b)]] = Q5(v)
        pinned_residual.append(row)
    require(len(pinned_residual) == 27, "RESIDUAL_COUNT", "expected twenty-seven pinned residual rows")
    require(rank_q5(pinned_residual) == 27, "RESIDUAL_RANK", "pinned residual rows are dependent")
    full_rref, full_pivots = rref_q5([[Q5(v) for v in row] for row in reps])
    require(len(full_rref) == 38, "FULL_RANK_Q5", "full rank over Q(sqrt(5)) is not 38")
    for row in pinned_residual:
        require(in_rowspace(full_rref, full_pivots, row), "RESIDUAL_IN_FULL", "residual row leaves the Jacobi rowspace")
    union_rank = rank_q5(prod_rref + pinned_residual)
    require(union_rank == 38, "SPLIT_UNION", "product plus residual rank is not 38")
    require(fl["full_rank_split"] == [11, 27, 38], "PINNED_SPLIT", "pinned split differs")

    # (2) structure lemmas.
    for a in range(PARAMS):
        for i in range(PORTS):
            for j in range(i + 1, PORTS):
                require(
                    sum((C[a][o][i][j] for o in range(PORTS)), Fraction(0)) == 0,
                    "NO_FIXED_OUTPUT",
                    "a bracket output has a fixed-sector component",
                )
    for a in range(PARAMS):
        claimed = [
            [
                transform[0][a] * projectors["three_plus"][i][j]
                + transform[1][a] * projectors["three_minus"][i][j]
                + transform[2][a] * projectors["five"][i][j]
                for j in range(PORTS)
            ]
            for i in range(PORTS)
        ]
        require(derivations[a] == claimed, "ADU_SPECTRAL", "ad_u spectral identity fails")
    for a in range(PARAMS):
        for b in range(PARAMS):
            gram = Q5_ZERO
            for i in range(PORTS):
                for j in range(PORTS):
                    if derivations[a][i][j] and derivations[b][j][i]:
                        gram = gram + derivations[a][i][j] * derivations[b][j][i]
            claimed = (
                Q5(3) * transform[0][a] * transform[0][b]
                + Q5(3) * transform[1][a] * transform[1][b]
                + Q5(5) * transform[2][a] * transform[2][b]
            )
            require(gram == claimed, "KILLING_FIXED_GRAM", "K(u,u) Gram identity fails")
    require(pair_orbit_count == 4, "COMMUTANT_DIM", "endomorphism algebra dimension is not four")

    # (3) case spans and decomposition facts.
    reps_q5 = [[Q5(v) for v in row] for row in reps]
    case_records = {}

    def span_of(support):
        pairs, rows = restricted_span_rows(reps_q5, inverse, support)
        rr, piv = rref_q5(rows)
        return pairs, rr, piv

    for case_id, support in [("A", [3]), ("B", [4]), ("C", [3, 4])]:
        pairs, rr, _ = span_of(support)
        require(len(rr) == 0, f"CASE_{case_id}", "restricted system is expected to vanish identically")
        case_records[case_id] = {"support": support, "span_dimension": 0, "pair_monomials": [list(p) for p in pairs]}

    def pair_row(pairs, terms):
        row = [Q5_ZERO] * len(pairs)
        for (i, j), v in terms.items():
            row[pairs.index(tuple(sorted((i, j))))] = row[pairs.index(tuple(sorted((i, j))))] + v
        return row

    # Case D: support {3,5,8,10}; claimed basis and decomposition facts.
    pairsD, rrD, pivD = span_of([3, 5, 8, 10])
    require(len(rrD) == 5, "CASE_D_DIM", "case D span dimension is not five")
    claimed_D = [
        pair_row(pairsD, {(3, 8): Q5_ONE}),
        pair_row(pairsD, {(5, 8): Q5_ONE}),
        pair_row(pairsD, {(8, 10): Q5_ONE}),
        pair_row(pairsD, {(10, 10): Q5_ONE, (3, 10): R5}),
        pair_row(pairsD, {(5, 10): R5, (3, 5): Q5(5), (8, 8): Q5(Fraction(-9, 10)) * (R5 + 1)}),
    ]
    for row in claimed_D:
        require(in_rowspace(rrD, pivD, row), "CASE_D_MEMBER", "claimed case D basis row leaves the span")
    require(rank_q5(claimed_D) == 5, "CASE_D_BASIS", "claimed case D basis is dependent")
    case_records["D"] = {
        "support": [3, 5, 8, 10],
        "span_dimension": 5,
        "pair_monomials": [list(p) for p in pairsD],
        "claimed_basis": [vec_encode(r) for r in claimed_D],
        "decomposition": [
            "the three product rows force c8=0 or c3=c5=c10=0; the last row then forces c8=0 outright",
            "with c8=0 the fourth row factors as c10*(c10+sqrt(5)*c3)",
            "on c10=-sqrt(5)*c3 the fifth row vanishes identically, giving the plane branch with c3, c5 free",
            "on c10=0 with c3 nonzero the fifth row forces c5=0, giving the c3 axis",
        ],
        "components": ["plane {c8=0, c10=-sqrt(5)*c3; c3,c5 free}", "line {c5=c8=c10=0; c3 free}"],
    }
    # verify the two stated decomposition constants
    fifth = claimed_D[4]
    # fifth with c8=0 equals c5*(sqrt5*c10 + 5*c3): at c10=-sqrt5*c3 the linear form is zero.
    lin_at_branch = R5 * (-R5) + Q5(5)
    require(not lin_at_branch, "CASE_D_BRANCH", "sqrt(5)*c10+5*c3 does not vanish on c10=-sqrt(5)*c3")
    require(Q5(Fraction(-9, 10)) * (R5 + 1), "CASE_D_C8SQ", "c8^2 coefficient vanishes unexpectedly")

    # Case E: support {4,6,12,13}; mirror.
    pairsE, rrE, pivE = span_of([4, 6, 12, 13])
    require(len(rrE) == 5, "CASE_E_DIM", "case E span dimension is not five")
    claimed_E = [
        pair_row(pairsE, {(4, 12): Q5_ONE}),
        pair_row(pairsE, {(6, 12): Q5_ONE}),
        pair_row(pairsE, {(12, 13): Q5_ONE}),
        pair_row(pairsE, {(13, 13): Q5_ONE, (4, 13): -R5}),
        pair_row(pairsE, {(6, 13): -R5, (4, 6): Q5(5), (12, 12): Q5(Fraction(9, 10)) * (R5 - 1)}),
    ]
    for row in claimed_E:
        require(in_rowspace(rrE, pivE, row), "CASE_E_MEMBER", "claimed case E basis row leaves the span")
    require(rank_q5(claimed_E) == 5, "CASE_E_BASIS", "claimed case E basis is dependent")
    lin_at_branch_E = -R5 * R5 + Q5(5)
    require(not lin_at_branch_E, "CASE_E_BRANCH", "-sqrt(5)*c13+5*c4 does not vanish on c13=sqrt(5)*c4")
    case_records["E"] = {
        "support": [4, 6, 12, 13],
        "span_dimension": 5,
        "pair_monomials": [list(p) for p in pairsE],
        "claimed_basis": [vec_encode(r) for r in claimed_E],
        "decomposition": [
            "the three product rows plus the fifth row force c12=0",
            "with c12=0 the fourth row factors as c13*(c13-sqrt(5)*c4)",
            "on c13=sqrt(5)*c4 the fifth row vanishes identically, giving the plane branch with c4, c6 free",
            "on c13=0 with c4 nonzero the fifth row forces c6=0, giving the c4 axis",
        ],
        "components": ["plane {c12=0, c13=sqrt(5)*c4; c4,c6 free}", "line {c6=c12=c13=0; c4 free}"],
    }

    # Cases F and G: coexistence; the span acquires no cross terms.
    pairsF, rrF, pivF = span_of([3, 4, 5, 8, 10])
    require(len(rrF) == 5, "CASE_F_DIM", "case F span dimension is not five")
    for row in rrF:
        for (i, j), v in zip(pairsF, row):
            if v:
                require(4 not in (i, j), "CASE_F_CROSS", "case F span contains a c4 monomial")
    pairsG, rrG, pivG = span_of([3, 4, 6, 12, 13])
    require(len(rrG) == 5, "CASE_G_DIM", "case G span dimension is not five")
    for row in rrG:
        for (i, j), v in zip(pairsG, row):
            if v:
                require(3 not in (i, j), "CASE_G_CROSS", "case G span contains a c3 monomial")
    case_records["F"] = {
        "support": [3, 4, 5, 8, 10],
        "span_dimension": 5,
        "no_monomial_involving": 4,
        "conclusion": "the case F variety is the c4 line times the case D variety",
    }
    case_records["G"] = {
        "support": [3, 4, 6, 12, 13],
        "span_dimension": 5,
        "no_monomial_involving": 3,
        "conclusion": "the case G variety is the c3 line times the case E variety",
    }

    # (4) per-channel structural facts for the family-level statements.
    sector_bases = {}
    for label in SECTORS:
        rows = [[projectors[label][i][j] for j in range(PORTS)] for i in range(PORTS)]
        sector_bases[label] = rref_q5(rows)[0]
    channel_facts = {}
    for k, (cid, left, right, output, _) in enumerate(CHANNEL_SPECS, start=3):
        _, B = bracket_tensor(C, inverse, {k: Q5_ONE})
        image = derived_span(B)
        image_sectors = sector_membership(image, projectors)
        fact = {
            "channel": CHANNEL_IDS[k],
            "image_dimension": len(image),
            "image_sectors": image_sectors,
            "image_is_full_output_sector": len(image) == SECTOR_DIMS[output] and image_sectors == [output],
        }
        if cid in ("t_pp_to_p", "t_mm_to_m", "t_ff_to_p", "t_ff_to_m"):
            domain = sector_bases[left]
            rows = []
            for vec in domain:
                flat = []
                for o in range(PORTS):
                    for j in range(PORTS):
                        flat.append(sum((vec[v] * B[o][v][j] for v in range(PORTS) if vec[v]), Q5_ZERO))
                rows.append(flat)
            fact["first_argument_kernel_dimension"] = SECTOR_DIMS[left] - rank_q5(rows)
            require(fact["first_argument_kernel_dimension"] == 0, "CHANNEL_KERNEL", f"{cid} has a first-argument kernel")
        channel_facts[CHANNEL_IDS[k]] = fact
    for cid in ("t_pp_to_p", "t_mm_to_m", "t_ff_to_p", "t_ff_to_m", "t_pf_to_f", "t_mf_to_f"):
        require(channel_facts[cid]["image_is_full_output_sector"], "CHANNEL_IMAGE", f"{cid} image is not the full output sector")

    # (5) families: Jacobi vanishing, Killing monomial matrices, verdicts.
    family_records = {}
    expected_monomials = {
        "PLANE": {(0, 0): "three_plus", (1, 1): "three_minus"},
        "SU3PF": {(0, 0): "three_plus", (0, 1): "five", (2, 2): "three_minus"},
        "SU3MF": {(0, 0): "three_minus", (0, 1): "five", (2, 2): "three_plus"},
        "DERIV": {(0, 0): "fixed", (1, 1): "fixed", (2, 2): "fixed"},
    }
    for family_id in ("PLANE", "SU3PF", "SU3MF", "DERIV"):
        n = FAMILY_PARAM_COUNT[family_id]
        for point in spanning_points(n):
            x, B = bracket_tensor(C, inverse, family_channel_values(family_id, point))
            values = evaluate_rows_at(reps_q5, x)
            require(not any(values), "FAMILY_JACOBI", f"{family_id} violates Jacobi at spanning point {point}")
        monomials = killing_monomial_matrices(C, inverse, family_id)
        require(set(monomials) == set(expected_monomials[family_id]), "FAMILY_MONOMIALS", f"{family_id} has unexpected Killing monomials")
        mono_records = {}
        for mono, M in monomials.items():
            supports = []
            for label in SECTORS:
                block = restrict_form(M, sector_bases[label])
                if any(any(row) for row in block):
                    supports.append((label, signature_q5(block)))
            require(len(supports) == 1 and supports[0][0] == expected_monomials[family_id][mono], "FAMILY_BLOCK", f"{family_id} monomial {mono} has unexpected sector support")
            label, sig = supports[0]
            mono_records["*".join(f"p{t}" for t in mono)] = {
                "monomial_parameter_indices": list(mono),
                "sector": label,
                "block_signature": sig,
                "matrix": sym_encode(M),
            }
        # cross-sector blocks vanish for every monomial matrix
        for mono, M in monomials.items():
            for ia, la in enumerate(SECTORS):
                for lb in SECTORS[ia + 1:]:
                    cross = [
                        [
                            sum((va[i] * M[i][j] * vb[j] for i in range(PORTS) for j in range(PORTS) if va[i] and M[i][j] and vb[j]), Q5_ZERO)
                            for vb in sector_bases[lb]
                        ]
                        for va in sector_bases[la]
                    ]
                    require(not any(any(row) for row in cross), "FAMILY_CROSS", f"{family_id} has a cross-sector Killing block")
        family_records[family_id] = {"killing_monomials": mono_records}

    # definiteness facts backing the verdicts
    def block_sig(family_id, mono):
        return family_records[family_id]["killing_monomials"]["*".join(f"p{t}" for t in mono)]["block_signature"]

    require(block_sig("PLANE", (0, 0)) == [0, 3, 0], "VERDICT_P1", "plane p-block is not negative definite")
    require(block_sig("PLANE", (1, 1)) == [0, 3, 0], "VERDICT_P2", "plane m-block is not negative definite")
    require(block_sig("SU3PF", (0, 0)) == [0, 3, 0], "VERDICT_F1", "SU3PF p-block is not negative definite")
    require(block_sig("SU3PF", (0, 1)) == [5, 0, 0], "VERDICT_F2", "SU3PF five-block is not positive definite")
    require(block_sig("SU3PF", (2, 2)) == [0, 3, 0], "VERDICT_F3", "SU3PF m-block is not negative definite")
    require(block_sig("SU3MF", (0, 0)) == [0, 3, 0], "VERDICT_G1", "SU3MF m-block is not negative definite")
    require(block_sig("SU3MF", (0, 1)) == [0, 5, 0], "VERDICT_G2", "SU3MF five-block is not negative definite")
    require(block_sig("SU3MF", (2, 2)) == [0, 3, 0], "VERDICT_G3", "SU3MF p-block is not negative definite")
    require(block_sig("DERIV", (0, 0)) == [1, 0, 0], "VERDICT_D1", "derivation a^2 block is not positive")
    require(block_sig("DERIV", (1, 1)) == [1, 0, 0], "VERDICT_D2", "derivation b^2 block is not positive")
    require(block_sig("DERIV", (2, 2)) == [1, 0, 0], "VERDICT_D3", "derivation e^2 block is not positive")

    # (6) sample points with the full compact-type criterion evaluated exactly.
    def analyze_sample(family_id, params):
        channel_values = family_channel_values(family_id, params)
        x, B = bracket_tensor(C, inverse, channel_values)
        defect = jacobi_defect(B)
        require(not defect, "SAMPLE_JACOBI", f"{family_id}{params} violates Jacobi")
        K = killing_matrix(B)
        derived = derived_span(B)
        center = center_basis(B)
        center_rref, center_piv = rref_q5(center) if center else ([], [])
        intersection = 0
        if derived and center_rref:
            combined = rank_q5(derived + center_rref)
            intersection = len(derived) + len(center_rref) - combined
        direct = intersection == 0 and len(derived) + len(center_rref) == PORTS
        negdef = True
        if derived:
            block = restrict_form(K, derived)
            sig = signature_q5(block)
            negdef = sig == [0, len(derived), 0]
        compact = direct and negdef
        return {
            "family": family_id,
            "parameters": [str(p) for p in params],
            "channel_values": {CHANNEL_IDS[k]: q5_encode(v) for k, v in sorted(channel_values.items())},
            "x_vector": vec_encode(x),
            "jacobi_zero": True,
            "killing_signature": signature_q5(K),
            "derived_dimension": len(derived),
            "derived_sectors": sector_membership(derived, projectors),
            "center_dimension": len(center_rref),
            "center_sectors": sector_membership(center_rref, projectors),
            "center_plus_derived_is_direct_and_full": direct,
            "killing_negative_definite_on_derived": negdef,
            "compact": compact,
        }

    samples = [
        ("PLANE", (0, 0), True), ("PLANE", (1, 0), True), ("PLANE", (1, 1), True),
        ("SU3PF", (1, -1, 0), True), ("SU3PF", (1, -1, 1), True),
        ("SU3PF", (1, 1, 0), False), ("SU3PF", (1, 0, 0), False), ("SU3PF", (0, 1, 0), False),
        ("SU3MF", (1, 1, 0), True), ("SU3MF", (1, 1, 1), True),
        ("SU3MF", (1, -1, 0), False), ("SU3MF", (1, 0, 0), False), ("SU3MF", (0, 1, 0), False),
        ("DERIV", (1, 0, 0), False), ("DERIV", (1, 1, 1), False),
    ]
    sample_records = []
    for family_id, params, expected in samples:
        record = analyze_sample(family_id, params)
        require(record["compact"] == expected, "SAMPLE_VERDICT", f"{family_id}{params} verdict differs from the derived expectation")
        sample_records.append(record)

    # (7) assemble the certificate.
    textbook_lemmas = {
        "TL1_compact_type_criterion": (
            "A finite-dimensional real Lie algebra is the Lie algebra of a compact group exactly when it is the "
            "vector-space direct sum of its center and its derived algebra and the Killing form is negative definite "
            "on the derived algebra. Consequence used for exclusion: the Killing form of a compact-type algebra is "
            "negative semidefinite."
        ),
        "TL2_compact_semisimple_dimensions": (
            "Every compact semisimple real Lie algebra is a direct sum of simple ideals; compact simple dimensions "
            "start 3, 8, 10, 14. There is no compact semisimple Lie algebra of dimension 5, and 8 forces su(3), "
            "3 forces so(3), while 11 splits only as 3+8. Automorphisms permute the simple ideals."
        ),
        "TL3_A5_simplicity": (
            "A5 is simple of order 60, so any action of A5 on at most four objects is trivial; every characteristic "
            "subspace of an equivariant bracket is therefore an A5-invariant subspace."
        ),
    }
    machine_premises = {
        "MP1_multiplicity_free": {
            "ordered_pair_orbit_count": 4,
            "statement": "dim End_A5(Q^12)=4, so the four sectors 1+3+3+5 are pairwise inequivalent irreducibles and every invariant subspace is a sum of sectors",
        },
        "MP2_no_fixed_sector_output": {
            "statement": "every basis bracket output sums to zero over ports, so the derived algebra of every bracket in the fourteen-parameter space avoids the fixed line",
        },
        "MP3_killing_on_fixed_generator": {
            "identity": "tr(ad_u ad_u) = 3*d_plus^2 + 3*d_minus^2 + 5*d_five^2 as an exact quadratic form on all fourteen parameters",
            "consequence": "a real point of compact type has d_plus=d_minus=d_five=0 by TL1; the fixed generator u is then central by MP4",
        },
        "MP4_ad_u_spectral": {
            "identity": "ad_u = d_plus*P_three_plus + d_minus*P_three_minus + d_five*P_five exactly",
        },
        "MP5_case_spans": "recorded per case in case_analysis",
        "MP6_families": "recorded per family in families and samples",
    }
    classification = {
        "statement": (
            "The compact real locus of the fourteen-parameter variety is contained in the slice "
            "d_plus=d_minus=d_five=0 and equals the union of the three families listed in compact_families. "
            "Completeness follows from TL1-TL3 with MP1-MP2: a compact-type point has derived algebra equal to a "
            "sum of sectors inside {three_plus, three_minus, five} carrying a compact semisimple bracket, the "
            "complementary sectors central; the admissible sector sums are {}, {p}, {m}, {p,m}, {p+f}, {m+f}, "
            "{p,m+f}, {m,p+f}; each case restricts the channels to one of the analyzed supports A-G, and the exact "
            "span decompositions plus the exact Killing blocks identify the compact points."
        ),
        "compact_families": [
            {
                "id": "P",
                "description": "closed plane: t_pp_to_p and t_mm_to_m free, all other channels zero",
                "dimension": 2,
                "compact_subset": "every point",
                "algebra": "so(3)+so(3)+R^6 generically, degenerating to so(3)+R^9 and R^12",
            },
            {
                "id": "F",
                "description": "t_pp_to_p=a, t_ff_to_p=b, t_pf_to_f=-sqrt(5)*a, t_mm_to_m=e, all other channels zero",
                "dimension": 3,
                "compact_subset": "a nonzero and a*b<0, e free",
                "algebra": "su(3)+so(3)+R generically (e nonzero); su(3)+R^4 at e=0",
                "sign_condition": "a*b < 0",
            },
            {
                "id": "G",
                "description": "t_mm_to_m=a, t_ff_to_m=b, t_mf_to_f=sqrt(5)*a, t_pp_to_p=e, all other channels zero",
                "dimension": 3,
                "compact_subset": "a nonzero and a*b>0, e free",
                "algebra": "su(3)+so(3)+R generically (e nonzero); su(3)+R^4 at e=0",
                "sign_condition": "a*b > 0",
            },
        ],
        "noncompact_strata_identified": [
            "sl(3,R)-type branch of family F (a nonzero, a*b>0): Killing signature (5,3,4)",
            "sl(3,R)-type branch of family G (a nonzero, a*b<0): Killing signature (5,3,4)",
            "so(3) acting on R^5: family F at b=0, a nonzero: degenerate Killing on the derived algebra",
            "two-step nilpotent axis: t_ff_to_p alone (and mirror): Killing identically zero, derived meets center",
            "derivation-only family: Killing signature (1,0,11) away from the origin",
        ],
    }
    calibration = {
        "family": "DERIV",
        "statement": (
            "The known three-parameter derivation-only family satisfies Jacobi identically; its Killing form is "
            "(3*a^2+3*b^2+5*e^2) on the fixed generator and zero elsewhere, so no point except the origin is of "
            "compact type. This executes the compact-type test on the known family as the calibration case."
        ),
    }
    certificate = {
        "schema": "oph.b14_compact_locus.certificate.v1",
        "issue": 705,
        "upstream": {
            "stage1_basis_path": str(STAGE1_BASIS.relative_to(REPO)),
            "stage1_basis_file_sha256": file_sha256(STAGE1_BASIS),
            "stage1_receipt_path": str(STAGE1_RECEIPT.relative_to(REPO)),
            "stage1_receipt_file_sha256": file_sha256(STAGE1_RECEIPT),
            "stage2_system_path": str(STAGE2_SYSTEM.relative_to(REPO)),
            "stage2_system_file_sha256": file_sha256(STAGE2_SYSTEM),
            "stage2_receipt_path": str(STAGE2_RECEIPT.relative_to(REPO)),
            "stage2_receipt_file_sha256": file_sha256(STAGE2_RECEIPT),
        },
        "reproduction": {
            "raw_jacobi_coordinate_count": 2640,
            "signed_orbit_count": 44,
            "orbit_sign_constancy": True,
            "representative_rank_over_Q": 38,
            "pivot_monomials_match_pinned": True,
            "representative_equations_match_pinned": 44,
            "contracted_row_count": 792,
            "contracted_rank_over_Q": 11,
            "pinned_contracted_rowspace_verified": True,
            "product_equation_rank_over_Qsqrt5": 11,
            "product_rowspace_equals_contracted": True,
            "residual_rank_over_Qsqrt5": 27,
            "product_plus_residual_rank": 38,
            "split_matches_pinned": [11, 27, 38],
            "channel_transform_matches_pinned": True,
        },
        "textbook_lemmas": textbook_lemmas,
        "machine_premises": machine_premises,
        "case_analysis": case_records,
        "channel_facts": channel_facts,
        "families": family_records,
        "samples": sample_records,
        "classification": classification,
        "calibration": calibration,
        "claim_boundary": (
            "This certificate classifies the compact real locus of the fourteen-parameter equivariant Jacobi "
            "variety, conditional on the three named textbook lemmas; every numerical statement is machine-checked "
            "exactly over Q(sqrt(5)). It does not enumerate the irreducible components of the noncompact part of "
            "the variety away from the analyzed channel supports, does not select a source or a preferred bracket, "
            "does not construct port perturbations, skew generators, recharting implementers, holonomy, or "
            "refinement intertwiners, and does not close issue #705."
        ),
        "implementation_pins": {
            "producer_sha256": file_sha256(Path(__file__).resolve()),
            "verifier_sha256": file_sha256(VERIFIER_PATH) if VERIFIER_PATH.exists() else "absent",
        },
    }
    certificate["certificate_sha256"] = object_sha256(certificate)
    CERTIFICATE_PATH.write_text(json.dumps(certificate, sort_keys=True, indent=1) + "\n", encoding="utf-8")
    return certificate


if __name__ == "__main__":
    result = main()
    print("certificate written:", CERTIFICATE_PATH)
    print("certificate_sha256:", result["certificate_sha256"])
    print("compact families:", [f["id"] for f in result["classification"]["compact_families"]])
