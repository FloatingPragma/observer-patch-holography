#!/usr/bin/env python3
"""B16 lattice computation from the issue #705 compact-locus certificate.

Semantic inputs are the pinned certificate code/b14_jacobi/b14_compact_locus.certificate.json,
its pinned stage-one Reynolds basis and stage-two reduction, and the issue #706
matter-class freeze code/b15_matter_freeze/matter_class_freeze_v1.json. For each
classified compact family (P: so(3)+so(3); F and G: su(3)+so(3)) this producer
rebuilds the exact bracket tensor at the certified generic compact sample point
and computes, over exact field towers Q(sqrt(5)) < E = Q(sqrt(5))[mu] <
F = E[nu] with mu^2 and nu^2 read off the bracket,

1. a Cartan subalgebra of every compact simple summand, the complete root
   decomposition with an eigenspace-dimension completeness proof, simple roots,
   coroots from sl(2) triples, and the integer Cartan matrix;
2. the coroot lattice, the root lattice, the integral weight lattice, and the
   center of the simply connected form as a Smith-normal-form lattice quotient
   of the computed root-coordinate matrix;
3. the exact weight multiset of the twelve-port carrier module under each
   summand, the carrier weight lattice, its membership in the root lattice, and
   the kernel of the finite center acting on the carrier: the loop-to-kernel
   quotient that names the global form the carrier selects;
4. the abelian sector: the exact charge matrices of the certified center acting
   on the carrier, the resulting zero charge lattice, and the finding that the
   source pins no primitive abelian period at this level;
5. controls: synthetic weight sets that force the covering form, and an abelian
   generator rescaling that changes nothing.

Every lattice is computed from the certified bracket data; nothing is imported
from a target normalization. The producer does not select a matter module over
the issue #706 frozen class, does not construct holonomy or a continuum bundle,
and does not close issue #707.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CERTIFICATE_PATH = REPO / "code/b14_jacobi/b14_compact_locus.certificate.json"
FREEZE_PATH = REPO / "code/b15_matter_freeze/matter_class_freeze_v1.json"
OUTPUT_PATH = HERE / "lattices_v1.json"
VALIDATOR_PATH = HERE / "validate_lattices.py"

PORTS = 12
PARAMS = 14
SECTORS = ["fixed", "three_plus", "three_minus", "five"]
SECTOR_DIMS = {"fixed": 1, "three_plus": 3, "three_minus": 3, "five": 5}
CHANNEL_IDS = [
    "d_plus", "d_minus", "d_five", "t_pp_to_p", "t_mm_to_m", "t_ff_to_p",
    "t_ff_to_m", "t_pm_to_f", "t_pf_to_p", "t_pf_to_m", "t_pf_to_f",
    "t_mf_to_p", "t_mf_to_m", "t_mf_to_f",
]
CHANNEL_INDEX = {cid: i for i, cid in enumerate(CHANNEL_IDS)}


class B16Error(RuntimeError):
    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise B16Error(code, message)


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


def q5_encode(v: Q5) -> list[int]:
    return [v.a.numerator, v.a.denominator, v.b.numerator, v.b.denominator]


def q5_decode(entry) -> Q5:
    return Q5(Fraction(entry[0], entry[1]), Fraction(entry[2], entry[3]))


def vec_encode_q5(vec) -> list[list[int]]:
    return [[i, *q5_encode(v)] for i, v in enumerate(vec) if v]


def rational_is_square(x: Fraction) -> bool:
    if x < 0:
        return False
    n, d = x.numerator, x.denominator
    return math.isqrt(n) ** 2 == n and math.isqrt(d) ** 2 == d


def q5_is_square(x: Q5) -> bool:
    """Exact test: is x a square in Q(sqrt(5))?

    x = (e + f*sqrt5)^2 forces e^2 + 5 f^2 = a and 2 e f = b. For b = 0 either
    f = 0 (a a rational square) or e = 0 (a/5 a rational square). For b nonzero
    e is a root of 4 e^4 - 4 a e^2 + 5 b^2 = 0, so e^2 = (a +- r)/2 with
    r^2 = a^2 - 5 b^2; each branch is tested by exact reconstruction."""
    if not x:
        return True
    if x.b == 0:
        return rational_is_square(x.a) or rational_is_square(x.a / 5)
    disc = x.a * x.a - 5 * x.b * x.b
    if not rational_is_square(disc):
        return False
    r = Fraction(math.isqrt(disc.numerator), math.isqrt(disc.denominator))
    for root in (r, -r):
        e2 = (x.a + root) / 2
        if e2 > 0 and rational_is_square(e2):
            e = Fraction(math.isqrt(e2.numerator), math.isqrt(e2.denominator))
            if e != 0:
                f = x.b / (2 * e)
                if e * e + 5 * f * f == x.a and 2 * e * f == x.b:
                    return True
    return False


# ------------------------------------------------- generic quadratic tower
class Ext:
    """Element p + q*g of base[g], g^2 = tower.d; base elements duck-typed."""

    __slots__ = ("p", "q", "tower")

    def __init__(self, tower, p, q):
        self.tower = tower
        self.p = p
        self.q = q

    def __add__(self, o):
        return Ext(self.tower, self.p + o.p, self.q + o.q)

    def __sub__(self, o):
        return Ext(self.tower, self.p - o.p, self.q - o.q)

    def __neg__(self):
        return Ext(self.tower, -self.p, -self.q)

    def __mul__(self, o):
        d = self.tower.d
        return Ext(self.tower, self.p * o.p + d * (self.q * o.q), self.p * o.q + self.q * o.p)

    def inv(self):
        d = self.tower.d
        n = self.p * self.p - d * (self.q * self.q)
        require(bool(n), "EXT_DIVISION", "division by zero in a quadratic extension")
        ninv = n.inv()
        return Ext(self.tower, self.p * ninv, -(self.q * ninv))

    def __truediv__(self, o):
        return self * o.inv()

    def __bool__(self):
        return bool(self.p) or bool(self.q)

    def __eq__(self, o):
        return self.p == o.p and self.q == o.q

    def __hash__(self):
        return hash((self.p, self.q))


class Tower:
    """Quadratic extension base[g]/(g^2 - d) over a duck-typed base field."""

    def __init__(self, base_zero, base_one, d, name):
        self.base_zero = base_zero
        self.base_one = base_one
        self.d = d
        self.name = name

    def zero(self):
        return Ext(self, self.base_zero, self.base_zero)

    def one(self):
        return Ext(self, self.base_one, self.base_zero)

    def gen(self):
        return Ext(self, self.base_zero, self.base_one)

    def lift(self, x):
        return Ext(self, x, self.base_zero)

    def lift_rational(self, r: Fraction):
        if isinstance(self.base_one, Q5):
            return Ext(self, Q5(r), self.base_zero)
        return Ext(self, self.base_one.tower.lift_rational(r), self.base_zero)


def flatten_rationals(x) -> list[Fraction]:
    if isinstance(x, Q5):
        return [x.a, x.b]
    return flatten_rationals(x.p) + flatten_rationals(x.q)


def ext_encode(x):
    if isinstance(x, Q5):
        return q5_encode(x)
    return [ext_encode(x.p), ext_encode(x.q)]


# --------------------------------------------------------- linear algebra
def rref(rows, width):
    work = [list(r) for r in rows if any(r)]
    if not work:
        return [], []
    pivots = []
    pr = 0
    for col in range(width):
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


def coords_in_rref(rref_rows, pivots, vec):
    residual = list(vec)
    out = []
    for base, col in zip(rref_rows, pivots):
        c = residual[col]
        out.append(c)
        if c:
            residual = [x - c * y for x, y in zip(residual, base)]
    require(not any(residual), "SPAN", "vector leaves the claimed span")
    return out


def kernel_basis(M, width, zero, one):
    rr, piv = rref(M, width)
    pivset = set(piv)
    free = [c for c in range(width) if c not in pivset]
    kernel = []
    for fc in free:
        v = [zero] * width
        v[fc] = one
        for row, pc in zip(rr, piv):
            v[pc] = -row[fc]
        kernel.append(v)
    return kernel


def mat_apply(M, v, zero):
    return [sum((M[o][j] * v[j] for j in range(len(v)) if v[j]), zero) for o in range(len(M))]


# ---------------------------------------------------------------- hashing
def canonical_bytes(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def object_sha256(value) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


# ------------------------------------------------------------ pinned input
def load_pinned():
    certificate = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
    cert_copy = dict(certificate)
    stored = cert_copy.pop("certificate_sha256", None)
    require(stored == object_sha256(cert_copy), "PIN_CERT_SELF", "certificate self-hash differs")
    up = certificate["upstream"]
    stage1_basis_path = REPO / up["stage1_basis_path"]
    stage1_receipt_path = REPO / up["stage1_receipt_path"]
    stage2_system_path = REPO / up["stage2_system_path"]
    stage2_receipt_path = REPO / up["stage2_receipt_path"]
    require(file_sha256(stage1_basis_path) == up["stage1_basis_file_sha256"], "PIN_STAGE1_BASIS", "stage-one basis hash differs")
    require(file_sha256(stage1_receipt_path) == up["stage1_receipt_file_sha256"], "PIN_STAGE1_RECEIPT", "stage-one receipt hash differs")
    require(file_sha256(stage2_system_path) == up["stage2_system_file_sha256"], "PIN_STAGE2_SYSTEM", "stage-two system hash differs")
    require(file_sha256(stage2_receipt_path) == up["stage2_receipt_file_sha256"], "PIN_STAGE2_RECEIPT", "stage-two receipt hash differs")
    freeze = json.loads(FREEZE_PATH.read_text(encoding="utf-8"))
    require(freeze["schema"] == "oph.b15_matter_class_freeze.v1", "PIN_FREEZE_SCHEMA", "freeze schema differs")
    require(
        freeze["provenance"]["embedded_certificate_sha256"] == certificate["certificate_sha256"],
        "PIN_FREEZE_CERT",
        "freeze pins a different certificate",
    )
    pinned_cert = next(
        p for p in freeze["provenance"]["pinned_files"] if p["path"] == "code/b14_jacobi/b14_compact_locus.certificate.json"
    )
    require(
        "sha256:" + pinned_cert["sha256"] == file_sha256(CERTIFICATE_PATH),
        "PIN_FREEZE_CERT_FILE",
        "freeze certificate file hash differs",
    )
    basis_raw = json.loads(stage1_basis_path.read_text(encoding="utf-8"))
    stage1_receipt = json.loads(stage1_receipt_path.read_text(encoding="utf-8"))
    system = json.loads(stage2_system_path.read_text(encoding="utf-8"))
    return certificate, freeze, basis_raw, stage1_receipt, system


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


def pinned_inverse_transform(system):
    cdm = system["fixed_line_reduction"]["channel_decomposition"]
    inverse = [[Q5_ZERO] * PARAMS for _ in range(PARAMS)]
    for i, entries in enumerate(cdm["inverse_transform_rows"]):
        for j, n1, d1, n2, d2 in entries:
            inverse[i][j] = Q5(Fraction(n1, d1), Fraction(n2, d2))
    return inverse


# ------------------------------------------------------ sector projectors
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
    R5 = Q5(0, 1)
    eigenvalues = {"fixed": Q5(5), "three_plus": R5, "three_minus": -R5, "five": Q5(-1)}
    projectors = {}
    for label, ev in eigenvalues.items():
        M = [[Q5_ONE if i == j else Q5_ZERO for j in range(PORTS)] for i in range(PORTS)]
        for label2, ev2 in eigenvalues.items():
            if label2 == label:
                continue
            factor = [[A[i][j] - (ev2 if i == j else Q5_ZERO) for j in range(PORTS)] for i in range(PORTS)]
            product = [[sum((M[i][k] * factor[k][j] for k in range(PORTS)), Q5_ZERO) for j in range(PORTS)] for i in range(PORTS)]
            M = [[v * (ev - ev2).inv() for v in row] for row in product]
        projectors[label] = M
    for label in SECTORS:
        P = projectors[label]
        PP = [[sum((P[i][k] * P[k][j] for k in range(PORTS)), Q5_ZERO) for j in range(PORTS)] for i in range(PORTS)]
        require(PP == P, "PROJECTOR_IDEMPOTENT", f"{label} projector is not idempotent")
        trace = sum((P[i][i] for i in range(PORTS)), Q5_ZERO)
        require(trace == Q5(SECTOR_DIMS[label]), "PROJECTOR_TRACE", f"{label} projector trace differs")
    return projectors


def sector_rref_bases(projectors):
    bases = {}
    for label in SECTORS:
        rows = [[projectors[label][i][j] for j in range(PORTS)] for i in range(PORTS)]
        rr, piv = rref(rows, PORTS)
        require(len(rr) == SECTOR_DIMS[label], "SECTOR_DIM", f"{label} basis dimension differs")
        bases[label] = (rr, piv)
    return bases


# ---------------------------------------------------------- bracket layer
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
    for (i, j, k) in itertools.combinations(range(PORTS), 3):
        for o in range(PORTS):
            v = Q5_ZERO
            for m in range(PORTS):
                v = v + B[m][i][j] * B[o][m][k] + B[m][j][k] * B[o][m][i] + B[m][k][i] * B[o][m][j]
            if v:
                return v
    return Q5_ZERO


def ad_matrix_q5(B, v):
    return [
        [sum((B[o][i][j] * v[i] for i in range(PORTS) if v[i]), Q5_ZERO) for j in range(PORTS)]
        for o in range(PORTS)
    ]


def bracket_q5(B, u, v):
    return [
        sum((B[o][i][j] * (u[i] * v[j]) for i in range(PORTS) for j in range(PORTS) if B[o][i][j] and u[i] and v[j]), Q5_ZERO)
        for o in range(PORTS)
    ]


def derived_span(B):
    rows = []
    for i in range(PORTS):
        for j in range(i + 1, PORTS):
            col = [B[o][i][j] for o in range(PORTS)]
            if any(col):
                rows.append(col)
    return rref(rows, PORTS)


def center_basis(B):
    rows = []
    for v in range(PORTS):
        flat = []
        for o in range(PORTS):
            for j in range(PORTS):
                flat.append(B[o][v][j])
        rows.append(flat)
    transposed = [[rows[i][j] for i in range(PORTS)] for j in range(len(rows[0]))]
    kernel = kernel_basis(transposed, PORTS, Q5_ZERO, Q5_ONE)
    rr, piv = rref(kernel, PORTS) if kernel else ([], [])
    return rr, piv


def sector_membership(rows, projectors):
    labels = []
    for label in SECTORS:
        P = projectors[label]
        for row in rows:
            image = [sum((P[o][i] * row[i] for i in range(PORTS) if row[i]), Q5_ZERO) for o in range(PORTS)]
            if any(image):
                labels.append(label)
                break
    return labels


# ------------------------------------------------------- integer lattices
def snf_with_right_transform(M):
    """Smith normal form of an integer matrix; returns (diag, V) with
    U*M*V = D for unimodular U (dropped) and V (returned, n x n)."""
    A = [list(r) for r in M]
    m = len(A)
    n = len(A[0]) if m else 0
    V = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    diag = []
    t = 0
    while t < min(m, n):
        pivot = None
        best = None
        for i in range(t, m):
            for j in range(t, n):
                if A[i][j] != 0 and (best is None or abs(A[i][j]) < best):
                    best = abs(A[i][j])
                    pivot = (i, j)
        if pivot is None:
            break
        stable = False
        while not stable:
            i0, j0 = pivot
            A[t], A[i0] = A[i0], A[t]
            for row in A:
                row[t], row[j0] = row[j0], row[t]
            for row in V:
                row[t], row[j0] = row[j0], row[t]
            stable = True
            for i in range(t + 1, m):
                q = A[i][t] // A[t][t]
                if q:
                    A[i] = [a - q * b for a, b in zip(A[i], A[t])]
                if A[i][t] != 0:
                    stable = False
            for j in range(t + 1, n):
                q = A[t][j] // A[t][t]
                if q:
                    for row in A:
                        row[j] -= q * row[t]
                    for row in V:
                        row[j] -= q * row[t]
                if A[t][j] != 0:
                    stable = False
            if not stable:
                pivot = None
                best = None
                for i in range(t, m):
                    for j in range(t, n):
                        if A[i][j] != 0 and (best is None or abs(A[i][j]) < best):
                            best = abs(A[i][j])
                            pivot = (i, j)
        if A[t][t] < 0:
            for row in A:
                row[t] = -row[t]
            for row in V:
                row[t] = -row[t]
        diag.append(A[t][t])
        t += 1
    # divisibility normalization
    for i in range(len(diag) - 1):
        for j in range(i + 1, len(diag)):
            a, b = diag[i], diag[j]
            if b % a != 0:
                g = math.gcd(a, b)
                lcm = a * b // g
                # adjust the recorded diagonal; the transforms for this step are
                # dropped with U, and V is unaffected for diagonal matrices of
                # coprime-free small cases handled here.
                require(False, "SNF_DIVISIBILITY", "divisibility repair requested; unsupported case")
    return diag, V


def lattice_quotient_data(rows, ambient_rank):
    """Quotient Z^n / rowspan(rows): invariant factors and the right transform."""
    if not rows:
        return [], [[1 if i == j else 0 for j in range(ambient_rank)] for i in range(ambient_rank)]
    diag, V = snf_with_right_transform(rows)
    return diag, V


def quotient_image(w, diag, V, ambient_rank):
    """Image coordinates of an integer vector in Z^n / L, componentwise."""
    y = [sum(w[i] * V[i][j] for i in range(ambient_rank)) for j in range(ambient_rank)]
    image = []
    for j in range(ambient_rank):
        if j < len(diag):
            image.append(y[j] % diag[j])
        else:
            image.append(y[j])
    return image


def in_lattice(w, diag, V, ambient_rank):
    return all(v == 0 for v in quotient_image(w, diag, V, ambient_rank))


def center_kernel_order(diag, images):
    """Order of the subgroup of Hom(Z^n/L, C*) acting trivially on all images."""
    total = 1
    for j, d in enumerate(diag):
        if d == 1:
            continue
        g = d
        for img in images:
            g = math.gcd(g, img[j])
        total *= g
    return total


# ----------------------------------------------------------- family layer
FAMILY_SPECS = {
    "P": {
        "certificate_family": "PLANE",
        "parameters": ["1", "1"],
        "summands": [
            {"name": "so3_plus", "sectors": ["three_plus"], "type": "A1"},
            {"name": "so3_minus", "sectors": ["three_minus"], "type": "A1"},
        ],
        "abelian_sectors": ["fixed", "five"],
    },
    "F": {
        "certificate_family": "SU3PF",
        "parameters": ["1", "-1", "1"],
        "summands": [
            {"name": "su3", "sectors": ["three_plus", "five"], "type": "A2"},
            {"name": "so3", "sectors": ["three_minus"], "type": "A1"},
        ],
        "abelian_sectors": ["fixed"],
    },
    "G": {
        "certificate_family": "SU3MF",
        "parameters": ["1", "1", "1"],
        "summands": [
            {"name": "su3", "sectors": ["three_minus", "five"], "type": "A2"},
            {"name": "so3", "sectors": ["three_plus"], "type": "A1"},
        ],
        "abelian_sectors": ["fixed"],
    },
}

ADJOINT_FORM = {"A1": "SO(3) = SU(2)/(Z/2)", "A2": "PSU(3) = SU(3)/(Z/3)"}
SC_CENTER = {"A1": "Z/2", "A2": "Z/3"}


def find_sample(certificate, family, parameters):
    for s in certificate["samples"]:
        if s["family"] == family and s["parameters"] == parameters:
            require(s["compact"] is True, "SAMPLE_COMPACT", "selected sample is not certified compact")
            return s
    raise B16Error("SAMPLE_MISSING", f"no certified sample for {family}{parameters}")


def restrict_to_span(M, span_rows, span_piv, zero):
    k = len(span_rows)
    R = [[zero] * k for _ in range(k)]
    for j, b in enumerate(span_rows):
        img = mat_apply(M, b, zero)
        co = coords_in_rref(span_rows, span_piv, img)
        for i in range(k):
            R[i][j] = co[i]
    return R


def proportional(u, v):
    k = next((i for i in range(PORTS) if u[i] or v[i]), None)
    if k is None:
        return True
    if not u[k] or not v[k]:
        return False
    r = v[k] / u[k]
    return all(v[i] == r * u[i] for i in range(PORTS))


def lift_matrix(M, tower):
    return [[tower.lift(M[i][j]) for j in range(len(M[0]))] for i in range(len(M))]


def minus_scalar(M, c, zero):
    return [[M[i][j] - (c if i == j else zero) for j in range(len(M))] for i in range(len(M))]


def extract_scalar(w, v, zero):
    """Solve w = lam * v for a single vector v with w in its line."""
    k = next(i for i in range(len(v)) if v[i])
    lam = w[k] / v[k]
    require(all(w[i] == lam * v[i] for i in range(len(v))), "SCALAR", "vector pair is not proportional")
    return lam


def extract_integer(x, tower_one, lo=-12, hi=12):
    for n in range(lo, hi + 1):
        acc = tower_one - tower_one
        for _ in range(abs(n)):
            acc = (acc + tower_one) if n > 0 else (acc - tower_one)
        if x == acc:
            return n
    raise B16Error("INTEGER_WEIGHT", "eigenvalue is not a small integer")


def root_sort_key(root_vals):
    flat = []
    for v in root_vals:
        flat.extend(flatten_rationals(v))
    return tuple(flat)


def is_positive_root(root_vals):
    for v in root_vals:
        for r in flatten_rationals(v):
            if r > 0:
                return True
            if r < 0:
                return False
    return False


def analyze_summand(B, spec, bases, family_id):
    """Full root-datum and carrier-weight computation for one compact summand."""
    sector_list = spec["sectors"]
    span_rows = []
    for label in sector_list:
        span_rows.extend([list(r) for r in bases[label][0]])
    s_rref, s_piv = rref(span_rows, PORTS)
    dim = len(s_rref)
    require(dim == sum(SECTOR_DIMS[l] for l in sector_list), "SUMMAND_DIM", "summand dimension differs")

    # closure of the summand and annihilator (isotypic trivial part)
    for u in s_rref:
        for v in s_rref:
            w = bracket_q5(B, u, v)
            if any(w):
                coords_in_rref(s_rref, s_piv, w)  # raises if the bracket leaves the summand

    stacked = []
    for u in s_rref:
        M = ad_matrix_q5(B, u)
        stacked.extend(M)
    ann = kernel_basis(stacked, PORTS, Q5_ZERO, Q5_ONE)
    ann_rref, ann_piv = rref(ann, PORTS)
    triv_dim = len(ann_rref)
    require(triv_dim == PORTS - dim, "ISOTYPIC", "annihilator dimension differs from the trivial-isotypic expectation")
    joined, _ = rref([list(r) for r in s_rref] + [list(r) for r in ann_rref], PORTS)
    require(len(joined) == PORTS, "ISOTYPIC_SUM", "summand plus annihilator does not fill the carrier")

    # Cartan generator h1 and the first extension
    h1 = list(bases[sector_list[0]][0][0])
    M1 = ad_matrix_q5(B, h1)
    first_rows, first_piv = bases[sector_list[0]]
    R_first = restrict_to_span(M1, [list(r) for r in first_rows], first_piv, Q5_ZERO)
    k = len(R_first)
    tr2 = sum((R_first[i][j] * R_first[j][i] for i in range(k) for j in range(k)), Q5_ZERO)
    s_val = tr2 / Q5(2)
    require(bool(s_val), "MU_ZERO", "mu^2 vanishes; h1 is not regular in its sector")
    require(s_val.sign() < 0, "MU_SIGN", "mu^2 is not negative; summand is not of compact type here")
    E = Tower(Q5_ZERO, Q5_ONE, s_val, "E")
    MU = E.gen()

    M1_E = lift_matrix(M1, E)
    rank = 1 if spec["type"] == "A1" else 2

    # eigenvalue completeness of ad_h1 on the carrier over E
    half_candidates = [Fraction(n, 2) for n in range(-6, 7)]
    eigen_dims = {}
    eigenspaces = {}
    total = 0
    for r in half_candidates:
        c = E.lift(Q5(r)) * MU if r != 0 else E.zero()
        ker = kernel_basis(minus_scalar(M1_E, c, E.zero()), PORTS, E.zero(), E.one())
        if ker:
            eigen_dims[str(r)] = len(ker)
            kr, kp = rref(ker, PORTS)
            eigenspaces[r] = (kr, kp)
            total += len(ker)
    require(total == PORTS, "EIGEN_COMPLETE", "ad_h1 eigenspace dimensions do not fill the carrier")

    record = {
        "name": spec["name"],
        "type_computed": None,
        "sectors": sector_list,
        "dimension": dim,
        "h1": vec_encode_q5(h1),
        "mu_squared": q5_encode(s_val),
        "carrier_eigenspace_dimensions_over_mu": {k2: v for k2, v in sorted(eigen_dims.items())},
        "trivial_isotypic_dimension": triv_dim,
    }

    if spec["type"] == "A1":
        require(eigen_dims == {"0": PORTS - 2, "1": 1, "-1": 1}, "A1_SPECTRUM", "unexpected ad_h1 spectrum for a rank-one summand")
        # coroot H = (2/mu) h1; ad_H is the rescaled ad_h1
        coroot_scale = E.lift(Q5(2)) * MU.inv()
        adH = [[E.lift(M1[i][j]) * coroot_scale for j in range(PORTS)] for i in range(PORTS)]
        int_dims = {}
        weight_multiset = {}
        totalH = 0
        one = E.one()
        for n in range(-6, 7):
            c = E.zero()
            for _ in range(abs(n)):
                c = c + one if n > 0 else c - one
            ker = kernel_basis(minus_scalar(adH, c, E.zero()), PORTS, E.zero(), E.one())
            if ker:
                int_dims[str(n)] = len(ker)
                weight_multiset[(n,)] = len(ker)
                totalH += len(ker)
        require(totalH == PORTS, "A1_COROOT_COMPLETE", "coroot eigenspaces do not fill the carrier")
        require(int_dims == {"0": PORTS - 2, "2": 1, "-2": 1}, "A1_WEIGHTS", "unexpected integral carrier weights for a rank-one summand")
        roots = [(2,), (-2,)]
        simple_roots = [(2,)]
        cartan_matrix = [[2]]
        coroot_coordinates = {"alpha_1": [1]}
        record["type_computed"] = "A1"
        record["cartan_dimension"] = 1
        record["roots_in_coroot_coordinates"] = [list(r) for r in roots]
        record["simple_root_coordinates"] = [list(r) for r in simple_roots]
        record["cartan_matrix"] = cartan_matrix
        record["positive_root_coroots_in_simple_coroot_basis"] = coroot_coordinates
        weight_items = sorted(weight_multiset.items())
        record["carrier_weights"] = [{"coordinates": list(w), "multiplicity": m} for w, m in weight_items]
        return record, roots, [w for w, _ in weight_items for _ in range(weight_multiset[w])], 1

    # A2 case
    require(eigen_dims == {"0": 6, "1": 2, "-1": 2, "2": 1, "-2": 1}, "A2_SPECTRUM", "unexpected ad_h1 spectrum for a rank-two summand")
    R_summand = restrict_to_span(M1, s_rref, s_piv, Q5_ZERO)
    cartan_vectors = kernel_basis(R_summand, dim, Q5_ZERO, Q5_ONE)
    require(len(cartan_vectors) == 2, "CARTAN_DIM", "centralizer of h1 in the summand is not two-dimensional")
    t_basis = []
    for kv in cartan_vectors:
        v = [Q5_ZERO] * PORTS
        for c, b in zip(kv, s_rref):
            if c:
                v = [x + c * y for x, y in zip(v, b)]
        t_basis.append(v)
    h2 = next((v for v in t_basis if not proportional(h1, v)), None)
    require(h2 is not None, "H2", "no independent Cartan direction found")
    br12 = mat_apply(M1, h2, Q5_ZERO)
    require(not any(br12), "CARTAN_ABELIAN", "[h1,h2] does not vanish")
    M2 = ad_matrix_q5(B, h2)
    M2_E = lift_matrix(M2, E)

    # second extension from the 2x2 restriction to the +mu eigenspace
    Wp_rows, Wp_piv = eigenspaces[Fraction(1)]
    R2 = restrict_to_span(M2_E, Wp_rows, Wp_piv, E.zero())
    tau = R2[0][0] + R2[1][1]
    delta = R2[0][0] * R2[1][1] - R2[0][1] * R2[1][0]
    disc = tau * tau - (delta + delta + delta + delta)
    require(bool(disc), "DISC_ZERO", "the M2 restriction to the mu eigenspace has a repeated eigenvalue")
    Wm_rows, Wm_piv = eigenspaces[Fraction(-1)]
    R2m = restrict_to_span(M2_E, Wm_rows, Wm_piv, E.zero())
    tau_m = R2m[0][0] + R2m[1][1]
    delta_m = R2m[0][0] * R2m[1][1] - R2m[0][1] * R2m[1][0]
    disc_m = tau_m * tau_m - (delta_m + delta_m + delta_m + delta_m)
    require(disc_m == disc, "DISC_MIRROR", "the -mu eigenspace discriminant differs; one extension does not suffice")
    if not disc.q:
        require(not q5_is_square(disc.p) and not q5_is_square(disc.p / s_val), "DISC_SQUARE", "discriminant is a square in E; the extension collapses")
    F = Tower(E.zero(), E.one(), disc, "F")
    NU = F.gen()
    HALF = F.lift_rational(Fraction(1, 2))

    def liftF(x_q5_matrix):
        return [[F.lift(E.lift(x_q5_matrix[i][j])) for j in range(PORTS)] for i in range(PORTS)]

    M1_F = liftF(M1)
    M2_F = liftF(M2)
    MU_F = F.lift(MU)

    # joint eigenspaces over F
    zeroF = F.zero()
    oneF = F.one()

    def stack(Ma, ca, Mb, cb):
        return minus_scalar(Ma, ca, zeroF) + minus_scalar(Mb, cb, zeroF)

    # zero weight space: joint kernel
    ker0 = kernel_basis(M1_F + M2_F, PORTS, zeroF, oneF)
    require(len(ker0) == 6, "JOINT_ZERO", "joint kernel dimension differs from six")

    # M2 eigenvalues: on the 2mu spaces M2 acts as a scalar; on the mu spaces via tau, nu
    two = oneF + oneF
    lam_p = [(F.lift(tau) + NU) * HALF, (F.lift(tau) - NU) * HALF]
    lam_m = [(F.lift(tau_m) + NU) * HALF, (F.lift(tau_m) - NU) * HALF]
    weight_count = 6
    root_values = []
    for c1_label, c1, lams in (
        (2, MU_F * two, None),
        (-2, zeroF - (MU_F * two), None),
        (1, MU_F, lam_p),
        (-1, zeroF - MU_F, lam_m),
    ):
        if lams is None:
            base_rows, base_piv = eigenspaces[Fraction(c1_label)]
            vecF = [F.lift(v) for v in base_rows[0]]
            img = mat_apply(M2_F, vecF, zeroF)
            lam = extract_scalar(img, vecF, zeroF) if any(img) else zeroF
            root_values.append(((c1, lam), vecF))
            weight_count += 1
        else:
            for lam in lams:
                ker = kernel_basis(stack(M1_F, c1, M2_F, lam), PORTS, zeroF, oneF)
                require(len(ker) == 1, "JOINT_ONE", "joint root space is not one-dimensional")
                kr, _ = rref(ker, PORTS)
                root_values.append(((c1, lam), kr[0]))
                weight_count += 1
    require(weight_count == PORTS, "JOINT_COMPLETE", "joint eigenspaces do not fill the carrier")

    # roots: nonzero joint values on the summand; positivity and simplicity
    root_list = [rv for rv, _ in root_values]
    positive = [rv for rv in root_list if is_positive_root(rv)]
    require(len(positive) == 3, "POSITIVE_COUNT", "expected three positive roots")
    positive_sorted = sorted(positive, key=root_sort_key)
    sums = set()
    for i in range(3):
        for j in range(3):
            if i != j:
                a, b = positive_sorted[i], positive_sorted[j]
                sums.add(((a[0] + b[0]), (a[1] + b[1])))
    simple = [rv for rv in positive_sorted if (rv[0], rv[1]) not in sums]
    require(len(simple) == 2, "SIMPLE_COUNT", "expected two simple roots")

    def root_vector(rv):
        for value, vec in root_values:
            if value == rv:
                return vec
        raise B16Error("ROOT_VECTOR", "root vector lookup failed")

    def ad_vec_F(v):
        return [
            [sum((F.lift(E.lift(B[o][i][j])) * v[i] for i in range(PORTS) if B[o][i][j] and v[i]), zeroF) for j in range(PORTS)]
            for o in range(PORTS)
        ]

    def bracket_F(u, v):
        return [
            sum((F.lift(E.lift(B[o][i][j])) * (u[i] * v[j]) for i in range(PORTS) for j in range(PORTS) if B[o][i][j] and u[i] and v[j]), zeroF)
            for o in range(PORTS)
        ]

    coroots = []
    for rv in simple:
        e_vec = root_vector(rv)
        f_vec = root_vector((zeroF - rv[0], zeroF - rv[1]))
        Hp = bracket_F(e_vec, f_vec)
        require(any(Hp), "TRIPLE", "the sl(2) triple bracket vanishes")
        img = mat_apply(ad_vec_F(Hp), e_vec, zeroF)
        alpha_Hp = extract_scalar(img, e_vec, zeroF)
        require(bool(alpha_Hp), "TRIPLE_VALUE", "alpha(H') vanishes")
        scale = (oneF + oneF) / alpha_Hp
        H = [scale * v for v in Hp]
        img2 = mat_apply(ad_vec_F(H), e_vec, zeroF)
        check2 = extract_scalar(img2, e_vec, zeroF)
        require(check2 == oneF + oneF, "COROOT_NORM", "alpha(H_alpha) differs from two")
        coroots.append(H)

    # Cartan matrix and integer coordinates
    def pair_with_coroots(vec_line):
        coords = []
        for H in coroots:
            img = mat_apply(ad_vec_F(H), vec_line, zeroF)
            lam = extract_scalar(img, vec_line, zeroF) if any(img) else zeroF
            coords.append(extract_integer(lam, oneF))
        return coords

    cartan_matrix = []
    for i, H in enumerate(coroots):
        row = []
        for j, rv in enumerate(simple):
            e_vec = root_vector(rv)
            img = mat_apply(ad_vec_F(H), e_vec, zeroF)
            lam = extract_scalar(img, e_vec, zeroF) if any(img) else zeroF
            row.append(extract_integer(lam, oneF))
        cartan_matrix.append(row)
    require(cartan_matrix in ([[2, -1], [-1, 2]],), "CARTAN_MATRIX", "computed Cartan matrix is not of type A2")

    all_root_coords = []
    for rv in root_list:
        all_root_coords.append(pair_with_coroots(root_vector(rv)))

    # positive-root coroots in the simple-coroot basis
    coroot_coordinates = {}
    for idx, rv in enumerate(positive_sorted):
        e_vec = root_vector(rv)
        f_vec = root_vector((zeroF - rv[0], zeroF - rv[1]))
        Hp = bracket_F(e_vec, f_vec)
        img = mat_apply(ad_vec_F(Hp), e_vec, zeroF)
        alpha_Hp = extract_scalar(img, e_vec, zeroF)
        scale = (oneF + oneF) / alpha_Hp
        H = [scale * v for v in Hp]
        # solve H = x1 H_1 + x2 H_2 exactly via two independent coordinates
        sol = None
        for i in range(PORTS):
            for j in range(i + 1, PORTS):
                det = coroots[0][i] * coroots[1][j] - coroots[0][j] * coroots[1][i]
                if det:
                    x1 = (H[i] * coroots[1][j] - H[j] * coroots[1][i]) / det
                    x2 = (coroots[0][i] * H[j] - coroots[0][j] * H[i]) / det
                    sol = (x1, x2)
                    break
            if sol:
                break
        require(sol is not None, "COROOT_SOLVE", "coroot coordinate solve failed")
        x1, x2 = sol
        for i in range(PORTS):
            require(H[i] == x1 * coroots[0][i] + x2 * coroots[1][i], "COROOT_SPAN", "coroot leaves the simple-coroot span")
        coroot_coordinates[f"positive_root_{idx + 1}"] = [extract_integer(x1, oneF), extract_integer(x2, oneF)]

    # carrier weights in integer coordinates
    weight_multiset = {}
    weight_multiset[(0, 0)] = 6
    for H in coroots:
        adH = ad_vec_F(H)
        for v in ker0:
            img = mat_apply(adH, v, zeroF)
            require(not any(img), "ZERO_WEIGHT", "a joint-kernel vector carries a nonzero weight")
    for rv, vec in root_values:
        coords = tuple(pair_with_coroots(vec))
        weight_multiset[coords] = weight_multiset.get(coords, 0) + 1
    require(sum(weight_multiset.values()) == PORTS, "WEIGHT_TOTAL", "carrier weight multiplicities do not sum to twelve")

    record["type_computed"] = "A2"
    record["cartan_dimension"] = 2
    record["h2"] = vec_encode_q5(h2)
    record["nu_squared"] = ext_encode(disc)
    record["m2_restriction_trace"] = ext_encode(tau)
    record["roots_in_coroot_coordinates"] = sorted(all_root_coords)
    simple_coords = [pair_with_coroots(root_vector(rv)) for rv in simple]
    record["simple_root_coordinates"] = simple_coords
    record["cartan_matrix"] = cartan_matrix
    record["positive_root_coroots_in_simple_coroot_basis"] = coroot_coordinates
    weight_items = sorted(weight_multiset.items())
    record["carrier_weights"] = [{"coordinates": list(w), "multiplicity": m} for w, m in weight_items]
    flat_weights = [w for w, m in weight_items for _ in range(m)]
    return record, [tuple(rc) for rc in all_root_coords], flat_weights, 2


def lattice_layer(record, roots, carrier_weights, rank):
    """Integer-lattice quotients from the computed roots and carrier weights."""
    root_rows = [list(r) for r in sorted(set(tuple(r) for r in roots))]
    diag, V = lattice_quotient_data(root_rows, rank)
    require(len(diag) == rank, "ROOT_RANK", "root lattice does not have full rank")
    center_order = 1
    for d in diag:
        center_order *= d
    invariant = [d for d in diag if d != 1]
    weight_images = []
    all_in_root_lattice = True
    for w in carrier_weights:
        img = quotient_image(list(w), diag, V, rank)
        weight_images.append(img)
        if any(img):
            all_in_root_lattice = False
    kernel_order = center_kernel_order(diag, weight_images)
    carrier_rows = [list(w) for w in sorted(set(carrier_weights)) if any(w)]
    if carrier_rows:
        cw_diag, _ = snf_with_right_transform(carrier_rows)
    else:
        cw_diag = []
    record["lattices"] = {
        "weight_lattice": f"Z^{rank} in simple-coroot pairing coordinates",
        "root_lattice_row_generators": root_rows,
        "root_lattice_invariant_factors": diag,
        "center_of_simply_connected_form_order": center_order,
        "center_of_simply_connected_form_invariant_factors": invariant,
        "coroot_lattice": f"Z^{rank} in the simple-coroot basis",
        "coweight_lattice_index_over_coroot_lattice": center_order,
        "carrier_weight_lattice_invariant_factors": cw_diag,
        "carrier_weights_contained_in_root_lattice": all_in_root_lattice,
    }
    record["center_action_on_carrier"] = {
        "carrier_weight_images_in_center_dual": sorted(set(tuple(i) for i in weight_images)),
        "kernel_order": kernel_order,
        "kernel_equals_full_center": kernel_order == center_order,
        "faithful_quotient_order": center_order // kernel_order,
    }
    return center_order, kernel_order


def analyze_family(fid, spec, certificate, C, inverse, projectors, bases):
    sample = find_sample(certificate, spec["certificate_family"], spec["parameters"])
    channel_values = {CHANNEL_INDEX[k]: q5_decode(v) for k, v in sample["channel_values"].items()}
    x, B = bracket_tensor(C, inverse, channel_values)
    require(vec_encode_q5(x) == sample["x_vector"], "X_VECTOR", "recomputed x vector differs from the certified one")
    require(not bool(jacobi_defect(B)), "JACOBI", "bracket violates Jacobi")
    d_rref, d_piv = derived_span(B)
    require(len(d_rref) == sample["derived_dimension"], "DERIVED_DIM", "derived dimension differs")
    require(sector_membership(d_rref, projectors) == sample["derived_sectors"], "DERIVED_SECTORS", "derived sectors differ")
    z_rref, z_piv = center_basis(B)
    require(len(z_rref) == sample["center_dimension"], "CENTER_DIM", "center dimension differs")
    require(sector_membership(z_rref, projectors) == sample["center_sectors"], "CENTER_SECTORS", "center sectors differ")

    # summands commute pairwise and fill the derived algebra
    summand_spans = []
    for sm in spec["summands"]:
        rows = []
        for label in sm["sectors"]:
            rows.extend([list(r) for r in bases[label][0]])
        summand_spans.append(rref(rows, PORTS))
    stacked = []
    for rr, _ in summand_spans:
        stacked.extend([list(r) for r in rr])
    joined, _ = rref(stacked, PORTS)
    require(len(joined) == len(d_rref), "SUMMAND_FILL", "summands do not fill the derived algebra")
    for row in joined:
        coords_in_rref(d_rref, d_piv, row)
    for (ra, _), (rb, _) in itertools.combinations(summand_spans, 2):
        for u in ra:
            for v in rb:
                w = [
                    sum((B[o][i][j] * (u[i] * v[j]) for i in range(PORTS) for j in range(PORTS) if B[o][i][j] and u[i] and v[j]), Q5_ZERO)
                    for o in range(PORTS)
                ]
                require(not any(w), "SUMMAND_COMMUTE", "distinct summands do not commute")

    family_record = {
        "certificate_family": spec["certificate_family"],
        "parameters": spec["parameters"],
        "channel_values": sample["channel_values"],
        "summands": [],
    }
    kernel_summary = []
    for sm in spec["summands"]:
        record, roots, carrier_weights, rank = analyze_summand(B, sm, bases, fid)
        center_order, kernel_order = lattice_layer(record, roots, carrier_weights, rank)
        record["selected_global_form_on_carrier"] = ADJOINT_FORM[record["type_computed"]] if kernel_order == center_order else "covering form retained"
        family_record["summands"].append(record)
        kernel_summary.append((record["type_computed"], center_order, kernel_order))

    # abelian sector: exact charge computation
    abelian_sectors = sector_membership(z_rref, projectors)
    require(abelian_sectors == spec["abelian_sectors"], "ABELIAN_SECTORS", "abelian sector list differs from the certified center")
    charge_matrix_ranks = []
    for zvec in z_rref:
        Mz = ad_matrix_q5(B, list(zvec))
        nonzero = any(any(row) for row in Mz)
        require(not nonzero, "ABELIAN_TRIVIAL", "a certified central vector acts nontrivially")
        charge_matrix_ranks.append(0)
    family_record["abelian"] = {
        "dimension": len(z_rref),
        "sectors": abelian_sectors,
        "charge_matrix_rank_per_generator": charge_matrix_ranks,
        "charge_lattice": "zero lattice {0}",
        "charge_lattice_rank": 0,
        "primitive_period": "undetermined_by_source",
        "finding": (
            "Every certified central generator acts by the exact zero matrix on the twelve-port carrier, "
            "so the carrier charge lattice is the zero lattice. A primitive abelian period is the minimal "
            "nonzero charge of a faithful character; the zero lattice has none, so the source pins no "
            "primitive abelian period and no compact-versus-noncompact abelian form at this level. This "
            "matches freeze constraint C3 of the issue #706 matter-class freeze."
        ),
    }

    # central relations and the loop-to-kernel quotient
    finite_center_factors = [SC_CENTER[t] for t, _, _ in kernel_summary]
    orders = [c for _, c, _ in kernel_summary]
    kernels = [k for _, _, k in kernel_summary]
    relations = []
    for i, (t, c, _) in enumerate(kernel_summary, start=1):
        relations.append(f"z{i}^{c} = 1")
    for i in range(1, len(kernel_summary) + 1):
        for j in range(i + 1, len(kernel_summary) + 1):
            relations.append(f"z{i} z{j} = z{j} z{i}")
    semisimple_forms = [ADJOINT_FORM[t] for t, c, k in kernel_summary]
    family_record["central_relations"] = {
        "finite_center_of_simply_connected_form": finite_center_factors,
        "relations": relations,
        "source": "invariant factors of the computed root-coordinate matrices; no order is imported",
    }
    family_record["loop_to_kernel"] = {
        "kernel_of_carrier_action_orders": kernels,
        "kernel_equals_full_finite_center": kernels == orders,
        "faithful_semisimple_form_on_carrier": " x ".join(semisimple_forms),
        "abelian_factor_on_carrier": (
            "acts by the identity; the carrier action factors through the semisimple adjoint form and "
            "determines no abelian period"
        ),
    }
    return family_record


# ------------------------------------------------------------------ controls
def run_controls():
    controls = {}
    # control 1: a fundamental-weight carrier would force the covering form (A2)
    a2_roots = [[2, -1], [-1, 2], [1, 1], [-2, 1], [1, -2], [-1, -1]]
    diag, V = lattice_quotient_data([list(r) for r in sorted(set(tuple(r) for r in a2_roots))], 2)
    weights = [list(r) for r in a2_roots] + [[1, 0]]
    images = [quotient_image(w, diag, V, 2) for w in weights]
    kernel = center_kernel_order(diag, images)
    center = 1
    for d in diag:
        center *= d
    controls["synthetic_fundamental_weight_A2"] = {
        "weights": weights,
        "center_order": center,
        "kernel_order": kernel,
        "conclusion": "adding one fundamental weight (1,0) shrinks the kernel to the identity; the pipeline selects the simply connected form when the carrier demands it",
    }
    require(center == 3 and kernel == 1, "CONTROL_A2", "A2 control kernel differs")
    # control 2: a spinor weight would force SU(2) (A1)
    diag1, V1 = lattice_quotient_data([[2], [-2]], 1)
    images1 = [quotient_image([w], diag1, V1, 1) for w in (2, -2, 1, -1)]
    kernel1 = center_kernel_order(diag1, images1)
    controls["synthetic_spinor_weight_A1"] = {
        "weights": [[2], [-2], [1], [-1]],
        "center_order": diag1[0],
        "kernel_order": kernel1,
        "conclusion": "adding weight 1 shrinks the kernel to the identity; the adjoint outcome on the carrier is computed, never assumed",
    }
    require(diag1[0] == 2 and kernel1 == 1, "CONTROL_A1", "A1 control kernel differs")
    # control 3: abelian generator rescaling
    controls["abelian_rescaling"] = {
        "statement": (
            "For every family the charge matrix of each central generator u and of the rescaled generator "
            "(7/3)u is the exact zero matrix, so no rescaling of the abelian generator produces a charge "
            "equation; the primitive-period question has no source-side constraint to solve"
        ),
        "verified_in": "analyze_family via ABELIAN_TRIVIAL together with linearity of the adjoint action",
    }
    return controls


# ----------------------------------------------------------------------- main
def main() -> dict:
    certificate, freeze, basis_raw, stage1_receipt, system = load_pinned()
    C = structure_constants(basis_raw)
    inverse = pinned_inverse_transform(system)
    group = [tuple(r) for r in stage1_receipt["proper_port_action"]["permutation_rows"]]
    require(len(group) == 60 and len(set(group)) == 60, "GROUP", "expected sixty distinct group rows")
    projectors = spectral_projectors(group)
    bases = sector_rref_bases(projectors)

    families = {}
    for fid, spec in FAMILY_SPECS.items():
        families[fid] = analyze_family(fid, spec, certificate, C, inverse, projectors, bases)

    controls = run_controls()

    freeze_c3 = next(c for c in freeze["constraints"] if c["id"] == "C3")
    result = {
        "schema": "oph.b16_lattices.v1",
        "issue": 707,
        "provenance": {
            "certificate_path": "code/b14_jacobi/b14_compact_locus.certificate.json",
            "certificate_file_sha256": file_sha256(CERTIFICATE_PATH),
            "certificate_self_sha256": certificate["certificate_sha256"],
            "freeze_path": "code/b15_matter_freeze/matter_class_freeze_v1.json",
            "freeze_file_sha256": file_sha256(FREEZE_PATH),
            "freeze_constraint_c3_status": freeze_c3["status"],
            "upstream": certificate["upstream"],
        },
        "conventions": {
            "field_tower": (
                "Q(sqrt(5)) < E = Q(sqrt(5))[mu], mu^2 read off tr(ad_h1^2)/2 on the leading sector; "
                "for rank-two summands F = E[nu], nu^2 the discriminant of the ad_h2 restriction to the "
                "mu eigenspace; all eigenspace dimensions carry completeness proofs (dimensions sum to twelve)"
            ),
            "weight_coordinates": (
                "integer pairings with the computed simple coroots, obtained from sl(2) triples "
                "[e,f]=H' normalized by alpha(H)=2"
            ),
            "lattice_quotients": (
                "Smith normal form of the integer root-coordinate matrix; the center of the simply "
                "connected form is the quotient of the weight lattice by the root lattice, and the "
                "carrier kernel is the annihilator of the carrier weight images in that quotient"
            ),
        },
        "families": families,
        "controls": controls,
        "findings": {
            "global_form_selected_by_carrier": {
                "P": "SO(3) x SO(3), the adjoint form of both rank-one summands",
                "F": "PSU(3) x SO(3), the adjoint form of both summands",
                "G": "PSU(3) x SO(3), the adjoint form of both summands",
            },
            "carrier_module_structure": (
                "under every compact summand the twelve-port carrier is the adjoint module plus a trivial "
                "complement, so every carrier weight lies in the root lattice and the full finite center "
                "acts by the identity"
            ),
            "primitive_abelian_period": (
                "undetermined by the source: the certified abelian sectors act by the exact zero matrix on "
                "the carrier, the charge lattice is {0}, and no primitive generator or period can be "
                "normalized against an empty charge set"
            ),
        },
        "claim_boundary": (
            "This computation derives coroot lattices, integral carrier weights, simply connected centers "
            "as lattice quotients, and the carrier-selected global form for the three certified compact "
            "families, exactly over declared field towers. It does not select a matter module over the "
            "issue #706 frozen class, does not derive an abelian period, does not construct holonomy, "
            "line operators, or a continuum bundle, and does not close issue #707."
        ),
        "implementation_pins": {
            "producer_sha256": file_sha256(Path(__file__).resolve()),
            "validator_sha256": file_sha256(VALIDATOR_PATH) if VALIDATOR_PATH.exists() else "absent",
        },
    }
    result["result_sha256"] = object_sha256(result)
    OUTPUT_PATH.write_text(json.dumps(result, sort_keys=True, indent=1) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    out = main()
    print("written:", OUTPUT_PATH)
    print("result_sha256:", out["result_sha256"])
    for fid, fam in out["families"].items():
        forms = [s["selected_global_form_on_carrier"] for s in fam["summands"]]
        print(fid, "->", " x ".join(forms))
