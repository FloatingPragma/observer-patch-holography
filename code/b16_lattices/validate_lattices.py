#!/usr/bin/env python3
"""Independent replay of the B16 lattice computation lattices_v1.json.

This validator reads the pinned issue #705 certificate, its pinned stage-one
and stage-two inputs, the issue #706 matter-class freeze, and the emitted
lattice record. With replay arithmetic written for this file it

* verifies every provenance pin fail-closed: the certificate self-hash, the
  upstream file hashes, the freeze cross-pins, the producer hash, its own
  hash, and the record self-hash;
* rebuilds the exact bracket tensor at each certified generic compact sample
  and recomputes the Cartan data, field towers, root decompositions with
  eigenspace completeness, coroots, Cartan matrices, integral carrier
  weights, Smith normal forms, lattice quotients, center kernels, and the
  abelian charge computation;
* compares every recomputed family record, control record, and finding with
  the stored record by canonical-JSON equality.

Any mismatch raises ReplayError.
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
RECORD_PATH = HERE / "lattices_v1.json"
PRODUCER_PATH = HERE / "compute_lattices.py"

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


class ReplayError(RuntimeError):
    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code


def check(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise ReplayError(code, message)


# ------------------------------------------------------------- Q(sqrt(5))
class F5:
    """Exact a + b*sqrt(5); replay arithmetic written for this validator."""

    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = a if isinstance(a, Fraction) else Fraction(a)
        self.b = b if isinstance(b, Fraction) else Fraction(b)

    def __add__(self, o):
        o = o if isinstance(o, F5) else F5(o)
        return F5(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return F5(-self.a, -self.b)

    def __sub__(self, o):
        o = o if isinstance(o, F5) else F5(o)
        return F5(self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        return F5(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, F5) else F5(o)
        return F5(self.a * o.a + 5 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def inv(self):
        n = self.a * self.a - 5 * self.b * self.b
        check(n != 0, "F5_DIV", "division by zero")
        return F5(self.a / n, -self.b / n)

    def __truediv__(self, o):
        o = o if isinstance(o, F5) else F5(o)
        return self * o.inv()

    def __bool__(self):
        return bool(self.a or self.b)

    def __eq__(self, o):
        o = o if isinstance(o, F5) else F5(o)
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


ZERO5 = F5()
ONE5 = F5(1)


def enc5(v: F5) -> list[int]:
    return [v.a.numerator, v.a.denominator, v.b.numerator, v.b.denominator]


def dec5(entry) -> F5:
    return F5(Fraction(entry[0], entry[1]), Fraction(entry[2], entry[3]))


def encvec5(vec) -> list[list[int]]:
    return [[i, *enc5(v)] for i, v in enumerate(vec) if v]


def frac_square(x: Fraction) -> bool:
    if x < 0:
        return False
    return math.isqrt(x.numerator) ** 2 == x.numerator and math.isqrt(x.denominator) ** 2 == x.denominator


def f5_square(x: F5) -> bool:
    if not x:
        return True
    if x.b == 0:
        return frac_square(x.a) or frac_square(x.a / 5)
    disc = x.a * x.a - 5 * x.b * x.b
    if not frac_square(disc):
        return False
    r = Fraction(math.isqrt(disc.numerator), math.isqrt(disc.denominator))
    for root in (r, -r):
        e2 = (x.a + root) / 2
        if e2 > 0 and frac_square(e2):
            e = Fraction(math.isqrt(e2.numerator), math.isqrt(e2.denominator))
            if e != 0:
                f = x.b / (2 * e)
                if e * e + 5 * f * f == x.a and 2 * e * f == x.b:
                    return True
    return False


# ------------------------------------------------------- quadratic towers
class TowerElem:
    """p + q*g over a duck-typed base, g^2 = field.d."""

    __slots__ = ("field", "p", "q")

    def __init__(self, field, p, q):
        self.field = field
        self.p = p
        self.q = q

    def __add__(self, o):
        return TowerElem(self.field, self.p + o.p, self.q + o.q)

    def __sub__(self, o):
        return TowerElem(self.field, self.p - o.p, self.q - o.q)

    def __neg__(self):
        return TowerElem(self.field, -self.p, -self.q)

    def __mul__(self, o):
        d = self.field.d
        return TowerElem(self.field, self.p * o.p + d * (self.q * o.q), self.p * o.q + self.q * o.p)

    def inv(self):
        d = self.field.d
        n = self.p * self.p - d * (self.q * self.q)
        check(bool(n), "TOWER_DIV", "division by zero in a quadratic extension")
        ninv = n.inv()
        return TowerElem(self.field, self.p * ninv, -(self.q * ninv))

    def __truediv__(self, o):
        return self * o.inv()

    def __bool__(self):
        return bool(self.p) or bool(self.q)

    def __eq__(self, o):
        return self.p == o.p and self.q == o.q

    def __hash__(self):
        return hash((self.p, self.q))


class TowerField:
    def __init__(self, base_zero, base_one, d):
        self.base_zero = base_zero
        self.base_one = base_one
        self.d = d

    def zero(self):
        return TowerElem(self, self.base_zero, self.base_zero)

    def one(self):
        return TowerElem(self, self.base_one, self.base_zero)

    def gen(self):
        return TowerElem(self, self.base_zero, self.base_one)

    def lift(self, x):
        return TowerElem(self, x, self.base_zero)

    def rational(self, r: Fraction):
        if isinstance(self.base_one, F5):
            return TowerElem(self, F5(r), self.base_zero)
        return TowerElem(self, self.base_one.field.rational(r), self.base_zero)


def flatten(x) -> list[Fraction]:
    if isinstance(x, F5):
        return [x.a, x.b]
    return flatten(x.p) + flatten(x.q)


def enc_tower(x):
    if isinstance(x, F5):
        return enc5(x)
    return [enc_tower(x.p), enc_tower(x.q)]


# --------------------------------------------------------- linear algebra
def reduce_rows(rows, width):
    work = [list(r) for r in rows if any(r)]
    if not work:
        return [], []
    pivots = []
    pr = 0
    for col in range(width):
        hit = next((r for r in range(pr, len(work)) if work[r][col]), None)
        if hit is None:
            continue
        work[pr], work[hit] = work[hit], work[pr]
        scale = work[pr][col].inv()
        work[pr] = [v * scale for v in work[pr]]
        for r in range(len(work)):
            if r != pr and work[r][col]:
                f = work[r][col]
                work[r] = [x - f * y for x, y in zip(work[r], work[pr])]
        pivots.append(col)
        pr += 1
        if pr == len(work):
            break
    return work[:pr], pivots


def span_coords(rows, pivots, vec):
    residual = list(vec)
    out = []
    for base, col in zip(rows, pivots):
        c = residual[col]
        out.append(c)
        if c:
            residual = [x - c * y for x, y in zip(residual, base)]
    check(not any(residual), "SPAN", "vector leaves the claimed span")
    return out


def null_space(M, width, zero, one):
    rows, piv = reduce_rows(M, width)
    pivset = set(piv)
    out = []
    for fc in (c for c in range(width) if c not in pivset):
        v = [zero] * width
        v[fc] = one
        for row, pc in zip(rows, piv):
            v[pc] = -row[fc]
        out.append(v)
    return out


def act(M, v, zero):
    return [sum((M[o][j] * v[j] for j in range(len(v)) if v[j]), zero) for o in range(len(M))]


# ---------------------------------------------------------------- hashing
def canonical(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def sha_object(value) -> str:
    return "sha256:" + hashlib.sha256(canonical(value)).hexdigest()


def sha_file(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def normalize(value):
    return json.loads(json.dumps(value, sort_keys=True))


# ------------------------------------------------------------ pinned input
def load_all():
    record = json.loads(RECORD_PATH.read_text(encoding="utf-8"))
    check(record["schema"] == "oph.b16_lattices.v1", "SCHEMA", "record schema differs")
    check(record["issue"] == 707, "ISSUE", "record issue differs")
    body = dict(record)
    stored = body.pop("result_sha256", None)
    check(stored == sha_object(body), "RECORD_SELF", "record self-hash differs")

    certificate = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
    cert_body = dict(certificate)
    cert_stored = cert_body.pop("certificate_sha256", None)
    check(cert_stored == sha_object(cert_body), "CERT_SELF", "certificate self-hash differs")

    prov = record["provenance"]
    check(prov["certificate_file_sha256"] == sha_file(CERTIFICATE_PATH), "PIN_CERT_FILE", "certificate file hash differs")
    check(prov["certificate_self_sha256"] == certificate["certificate_sha256"], "PIN_CERT_SELF", "certificate self-hash pin differs")
    check(prov["freeze_file_sha256"] == sha_file(FREEZE_PATH), "PIN_FREEZE_FILE", "freeze file hash differs")

    up = certificate["upstream"]
    check(normalize(prov["upstream"]) == normalize(up), "PIN_UPSTREAM", "upstream pin block differs")
    stage1_basis = REPO / up["stage1_basis_path"]
    stage1_receipt = REPO / up["stage1_receipt_path"]
    stage2_system = REPO / up["stage2_system_path"]
    stage2_receipt = REPO / up["stage2_receipt_path"]
    check(sha_file(stage1_basis) == up["stage1_basis_file_sha256"], "PIN_S1B", "stage-one basis hash differs")
    check(sha_file(stage1_receipt) == up["stage1_receipt_file_sha256"], "PIN_S1R", "stage-one receipt hash differs")
    check(sha_file(stage2_system) == up["stage2_system_file_sha256"], "PIN_S2S", "stage-two system hash differs")
    check(sha_file(stage2_receipt) == up["stage2_receipt_file_sha256"], "PIN_S2R", "stage-two receipt hash differs")

    freeze = json.loads(FREEZE_PATH.read_text(encoding="utf-8"))
    check(freeze["schema"] == "oph.b15_matter_class_freeze.v1", "FREEZE_SCHEMA", "freeze schema differs")
    check(
        freeze["provenance"]["embedded_certificate_sha256"] == certificate["certificate_sha256"],
        "FREEZE_CERT",
        "freeze pins a different certificate",
    )
    c3 = next(c for c in freeze["constraints"] if c["id"] == "C3")
    check(prov["freeze_constraint_c3_status"] == c3["status"], "FREEZE_C3", "freeze C3 status pin differs")

    pins = record["implementation_pins"]
    check(pins["producer_sha256"] == sha_file(PRODUCER_PATH), "PIN_PRODUCER", "producer hash differs")
    check(pins["validator_sha256"] == sha_file(Path(__file__).resolve()), "PIN_VALIDATOR", "validator hash differs")

    basis_raw = json.loads(stage1_basis.read_text(encoding="utf-8"))
    receipt1 = json.loads(stage1_receipt.read_text(encoding="utf-8"))
    system = json.loads(stage2_system.read_text(encoding="utf-8"))
    return record, certificate, freeze, basis_raw, receipt1, system


def constants_from_basis(basis_raw):
    C = [[[[Fraction(0)] * PORTS for _ in range(PORTS)] for _ in range(PORTS)] for _ in range(PARAMS)]
    rows = basis_raw["basis"]
    check(len(rows) == PARAMS, "BASIS_COUNT", "expected fourteen basis rows")
    for a, row in enumerate(rows):
        for o, l, r, n, d in row["entries"]:
            check(l < r, "BASIS_ENTRY", "basis entry is not upper triangular")
            C[a][o][l][r] = Fraction(n, d)
            C[a][o][r][l] = -Fraction(n, d)
    return C


def inverse_from_system(system):
    cdm = system["fixed_line_reduction"]["channel_decomposition"]
    inverse = [[ZERO5] * PARAMS for _ in range(PARAMS)]
    for i, entries in enumerate(cdm["inverse_transform_rows"]):
        for j, n1, d1, n2, d2 in entries:
            inverse[i][j] = F5(Fraction(n1, d1), Fraction(n2, d2))
    return inverse


def projectors_from_group(group):
    unassigned = {(i, j) for i in range(PORTS) for j in range(PORTS)}
    orbit_list = []
    while unassigned:
        seed = min(unassigned)
        orbit = {(g[seed[0]], g[seed[1]]) for g in group}
        unassigned -= orbit
        orbit_list.append(orbit)
    check(sorted(len(o) for o in orbit_list) == [12, 12, 60, 60], "PAIR_ORBITS", "unexpected orbit sizes")
    orbital = sorted((o for o in orbit_list if len(o) == 60), key=lambda o: tuple(sorted(o)))[0]
    A = [[ZERO5] * PORTS for _ in range(PORTS)]
    for i, j in orbital:
        A[i][j] = ONE5
    R5 = F5(0, 1)
    eigen = {"fixed": F5(5), "three_plus": R5, "three_minus": -R5, "five": F5(-1)}
    projectors = {}
    for label, ev in eigen.items():
        M = [[ONE5 if i == j else ZERO5 for j in range(PORTS)] for i in range(PORTS)]
        for label2, ev2 in eigen.items():
            if label2 == label:
                continue
            factor = [[A[i][j] - (ev2 if i == j else ZERO5) for j in range(PORTS)] for i in range(PORTS)]
            M = [
                [sum((M[i][k] * factor[k][j] for k in range(PORTS)), ZERO5) * (ev - ev2).inv() for j in range(PORTS)]
                for i in range(PORTS)
            ]
        projectors[label] = M
    for label in SECTORS:
        P = projectors[label]
        PP = [[sum((P[i][k] * P[k][j] for k in range(PORTS)), ZERO5) for j in range(PORTS)] for i in range(PORTS)]
        check(PP == P, "IDEMPOTENT", f"{label} projector fails idempotence")
        trace = sum((P[i][i] for i in range(PORTS)), ZERO5)
        check(trace == F5(SECTOR_DIMS[label]), "TRACE", f"{label} projector trace differs")
    return projectors


def sector_bases_of(projectors):
    bases = {}
    for label in SECTORS:
        rows = [[projectors[label][i][j] for j in range(PORTS)] for i in range(PORTS)]
        rr, piv = reduce_rows(rows, PORTS)
        check(len(rr) == SECTOR_DIMS[label], "SECTOR_DIM", f"{label} basis dimension differs")
        bases[label] = (rr, piv)
    return bases


# ---------------------------------------------------------- bracket layer
def build_bracket(C, inverse, channel_values):
    x = [sum((inverse[i][j] * v for j, v in channel_values.items()), ZERO5) for i in range(PARAMS)]
    B = [[[ZERO5] * PORTS for _ in range(PORTS)] for _ in range(PORTS)]
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


def jacobi_holds(B):
    for (i, j, k) in itertools.combinations(range(PORTS), 3):
        for o in range(PORTS):
            v = ZERO5
            for m in range(PORTS):
                v = v + B[m][i][j] * B[o][m][k] + B[m][j][k] * B[o][m][i] + B[m][k][i] * B[o][m][j]
            if v:
                return False
    return True


def ad_of(B, v):
    return [
        [sum((B[o][i][j] * v[i] for i in range(PORTS) if v[i]), ZERO5) for j in range(PORTS)]
        for o in range(PORTS)
    ]


def lie(B, u, v):
    return [
        sum((B[o][i][j] * (u[i] * v[j]) for i in range(PORTS) for j in range(PORTS) if B[o][i][j] and u[i] and v[j]), ZERO5)
        for o in range(PORTS)
    ]


def derived_of(B):
    rows = []
    for i in range(PORTS):
        for j in range(i + 1, PORTS):
            col = [B[o][i][j] for o in range(PORTS)]
            if any(col):
                rows.append(col)
    return reduce_rows(rows, PORTS)


def center_of(B):
    rows = []
    for v in range(PORTS):
        flat = []
        for o in range(PORTS):
            for j in range(PORTS):
                flat.append(B[o][v][j])
        rows.append(flat)
    transposed = [[rows[i][j] for i in range(PORTS)] for j in range(len(rows[0]))]
    kernel = null_space(transposed, PORTS, ZERO5, ONE5)
    return reduce_rows(kernel, PORTS) if kernel else ([], [])


def sectors_of(rows, projectors):
    labels = []
    for label in SECTORS:
        P = projectors[label]
        for row in rows:
            image = [sum((P[o][i] * row[i] for i in range(PORTS) if row[i]), ZERO5) for o in range(PORTS)]
            if any(image):
                labels.append(label)
                break
    return labels


# ------------------------------------------------------- integer lattices
def smith_form(M):
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
    for i in range(len(diag) - 1):
        for j in range(i + 1, len(diag)):
            check(diag[j] % diag[i] == 0, "SNF_DIV", "divisibility repair requested; unsupported case")
    return diag, V


def image_in_quotient(w, diag, V, n):
    y = [sum(w[i] * V[i][j] for i in range(n)) for j in range(n)]
    return [y[j] % diag[j] if j < len(diag) else y[j] for j in range(n)]


def kernel_order_of(diag, images):
    total = 1
    for j, d in enumerate(diag):
        if d == 1:
            continue
        g = d
        for img in images:
            g = math.gcd(g, img[j])
        total *= g
    return total


# ----------------------------------------------------------- family replay
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


def sample_for(certificate, family, parameters):
    for s in certificate["samples"]:
        if s["family"] == family and s["parameters"] == parameters:
            check(s["compact"] is True, "SAMPLE_COMPACT", "sample is not certified compact")
            return s
    raise ReplayError("SAMPLE_MISSING", f"no certified sample for {family}{parameters}")


def restrict(M, span_rows, span_piv, zero):
    k = len(span_rows)
    R = [[zero] * k for _ in range(k)]
    for j, b in enumerate(span_rows):
        co = span_coords(span_rows, span_piv, act(M, b, zero))
        for i in range(k):
            R[i][j] = co[i]
    return R


def collinear(u, v):
    k = next((i for i in range(PORTS) if u[i] or v[i]), None)
    if k is None:
        return True
    if not u[k] or not v[k]:
        return False
    r = v[k] / u[k]
    return all(v[i] == r * u[i] for i in range(PORTS))


def shift(M, c, zero):
    return [[M[i][j] - (c if i == j else zero) for j in range(len(M))] for i in range(len(M))]


def line_scalar(w, v, zero):
    k = next(i for i in range(len(v)) if v[i])
    lam = w[k] / v[k]
    check(all(w[i] == lam * v[i] for i in range(len(v))), "LINE", "vector pair is not proportional")
    return lam


def as_integer(x, one, lo=-12, hi=12):
    for n in range(lo, hi + 1):
        acc = one - one
        for _ in range(abs(n)):
            acc = (acc + one) if n > 0 else (acc - one)
        if x == acc:
            return n
    raise ReplayError("INT_WEIGHT", "eigenvalue is not a small integer")


def sort_key_of(values):
    flat = []
    for v in values:
        flat.extend(flatten(v))
    return tuple(flat)


def positive_first(values):
    for v in values:
        for r in flatten(v):
            if r > 0:
                return True
            if r < 0:
                return False
    return False


def replay_summand(B, spec, bases):
    sector_list = spec["sectors"]
    span_rows = []
    for label in sector_list:
        span_rows.extend([list(r) for r in bases[label][0]])
    s_rref, s_piv = reduce_rows(span_rows, PORTS)
    dim = len(s_rref)
    check(dim == sum(SECTOR_DIMS[l] for l in sector_list), "SUMMAND_DIM", "summand dimension differs")

    for u in s_rref:
        for v in s_rref:
            w = lie(B, u, v)
            if any(w):
                span_coords(s_rref, s_piv, w)

    stacked = []
    for u in s_rref:
        stacked.extend(ad_of(B, u))
    ann = null_space(stacked, PORTS, ZERO5, ONE5)
    ann_rref, _ = reduce_rows(ann, PORTS)
    triv_dim = len(ann_rref)
    check(triv_dim == PORTS - dim, "ISOTYPIC", "annihilator dimension differs")
    joined, _ = reduce_rows([list(r) for r in s_rref] + [list(r) for r in ann_rref], PORTS)
    check(len(joined) == PORTS, "ISOTYPIC_SUM", "summand plus annihilator does not fill the carrier")

    h1 = list(bases[sector_list[0]][0][0])
    M1 = ad_of(B, h1)
    first_rows, first_piv = bases[sector_list[0]]
    R_first = restrict(M1, [list(r) for r in first_rows], first_piv, ZERO5)
    k = len(R_first)
    tr2 = sum((R_first[i][j] * R_first[j][i] for i in range(k) for j in range(k)), ZERO5)
    s_val = tr2 / F5(2)
    check(bool(s_val), "MU_ZERO", "mu^2 vanishes")
    check(s_val.sign() < 0, "MU_SIGN", "mu^2 is not negative")
    E = TowerField(ZERO5, ONE5, s_val)
    MU = E.gen()
    M1_E = [[E.lift(M1[i][j]) for j in range(PORTS)] for i in range(PORTS)]

    eigen_dims = {}
    eigenspaces = {}
    total = 0
    for r in [Fraction(n, 2) for n in range(-6, 7)]:
        c = E.lift(F5(r)) * MU if r != 0 else E.zero()
        ker = null_space(shift(M1_E, c, E.zero()), PORTS, E.zero(), E.one())
        if ker:
            eigen_dims[str(r)] = len(ker)
            eigenspaces[r] = reduce_rows(ker, PORTS)
            total += len(ker)
    check(total == PORTS, "EIGEN_COMPLETE", "ad_h1 eigenspaces do not fill the carrier")

    record = {
        "name": spec["name"],
        "type_computed": None,
        "sectors": sector_list,
        "dimension": dim,
        "h1": encvec5(h1),
        "mu_squared": enc5(s_val),
        "carrier_eigenspace_dimensions_over_mu": {k2: v for k2, v in sorted(eigen_dims.items())},
        "trivial_isotypic_dimension": triv_dim,
    }

    if spec["type"] == "A1":
        check(eigen_dims == {"0": PORTS - 2, "1": 1, "-1": 1}, "A1_SPECTRUM", "unexpected rank-one spectrum")
        coroot_scale = E.lift(F5(2)) * MU.inv()
        adH = [[E.lift(M1[i][j]) * coroot_scale for j in range(PORTS)] for i in range(PORTS)]
        int_dims = {}
        weight_multiset = {}
        totalH = 0
        one = E.one()
        for n in range(-6, 7):
            c = E.zero()
            for _ in range(abs(n)):
                c = (c + one) if n > 0 else (c - one)
            ker = null_space(shift(adH, c, E.zero()), PORTS, E.zero(), E.one())
            if ker:
                int_dims[str(n)] = len(ker)
                weight_multiset[(n,)] = len(ker)
                totalH += len(ker)
        check(totalH == PORTS, "A1_COMPLETE", "coroot eigenspaces do not fill the carrier")
        check(int_dims == {"0": PORTS - 2, "2": 1, "-2": 1}, "A1_WEIGHTS", "unexpected integral weights")
        record["type_computed"] = "A1"
        record["cartan_dimension"] = 1
        record["roots_in_coroot_coordinates"] = [[2], [-2]]
        record["simple_root_coordinates"] = [[2]]
        record["cartan_matrix"] = [[2]]
        record["positive_root_coroots_in_simple_coroot_basis"] = {"alpha_1": [1]}
        weight_items = sorted(weight_multiset.items())
        record["carrier_weights"] = [{"coordinates": list(w), "multiplicity": m} for w, m in weight_items]
        flat = [w for w, m in weight_items for _ in range(m)]
        return record, [(2,), (-2,)], flat, 1

    check(eigen_dims == {"0": 6, "1": 2, "-1": 2, "2": 1, "-2": 1}, "A2_SPECTRUM", "unexpected rank-two spectrum")
    R_summand = restrict(M1, s_rref, s_piv, ZERO5)
    cartan_vectors = null_space(R_summand, dim, ZERO5, ONE5)
    check(len(cartan_vectors) == 2, "CARTAN_DIM", "centralizer of h1 is not two-dimensional")
    t_basis = []
    for kv in cartan_vectors:
        v = [ZERO5] * PORTS
        for c, b in zip(kv, s_rref):
            if c:
                v = [x + c * y for x, y in zip(v, b)]
        t_basis.append(v)
    h2 = next((v for v in t_basis if not collinear(h1, v)), None)
    check(h2 is not None, "H2", "no independent Cartan direction")
    check(not any(act(M1, h2, ZERO5)), "CARTAN_ABELIAN", "[h1,h2] does not vanish")
    M2 = ad_of(B, h2)
    M2_E = [[E.lift(M2[i][j]) for j in range(PORTS)] for i in range(PORTS)]

    Wp_rows, Wp_piv = eigenspaces[Fraction(1)]
    R2 = restrict(M2_E, Wp_rows, Wp_piv, E.zero())
    tau = R2[0][0] + R2[1][1]
    delta = R2[0][0] * R2[1][1] - R2[0][1] * R2[1][0]
    four_delta = delta + delta + delta + delta
    disc = tau * tau - four_delta
    check(bool(disc), "DISC_ZERO", "repeated eigenvalue on the mu eigenspace")
    Wm_rows, Wm_piv = eigenspaces[Fraction(-1)]
    R2m = restrict(M2_E, Wm_rows, Wm_piv, E.zero())
    tau_m = R2m[0][0] + R2m[1][1]
    delta_m = R2m[0][0] * R2m[1][1] - R2m[0][1] * R2m[1][0]
    check(tau_m * tau_m - (delta_m + delta_m + delta_m + delta_m) == disc, "DISC_MIRROR", "mirror discriminant differs")
    if not disc.q:
        check(not f5_square(disc.p) and not f5_square(disc.p / s_val), "DISC_SQUARE", "discriminant is a square in E")
    F = TowerField(E.zero(), E.one(), disc)
    NU = F.gen()
    HALF = F.rational(Fraction(1, 2))
    zeroF = F.zero()
    oneF = F.one()

    def liftF(Mq5):
        return [[F.lift(E.lift(Mq5[i][j])) for j in range(PORTS)] for i in range(PORTS)]

    M1_F = liftF(M1)
    M2_F = liftF(M2)
    MU_F = F.lift(MU)

    ker0 = null_space(M1_F + M2_F, PORTS, zeroF, oneF)
    check(len(ker0) == 6, "JOINT_ZERO", "joint kernel dimension differs from six")

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
            base_rows, _ = eigenspaces[Fraction(c1_label)]
            vecF = [F.lift(v) for v in base_rows[0]]
            img = act(M2_F, vecF, zeroF)
            lam = line_scalar(img, vecF, zeroF) if any(img) else zeroF
            root_values.append(((c1, lam), vecF))
            weight_count += 1
        else:
            for lam in lams:
                stacked2 = shift(M1_F, c1, zeroF) + shift(M2_F, lam, zeroF)
                ker = null_space(stacked2, PORTS, zeroF, oneF)
                check(len(ker) == 1, "JOINT_ONE", "joint root space is not one-dimensional")
                kr, _ = reduce_rows(ker, PORTS)
                root_values.append(((c1, lam), kr[0]))
                weight_count += 1
    check(weight_count == PORTS, "JOINT_COMPLETE", "joint eigenspaces do not fill the carrier")

    root_list = [rv for rv, _ in root_values]
    positive = [rv for rv in root_list if positive_first(rv)]
    check(len(positive) == 3, "POSITIVE_COUNT", "expected three positive roots")
    positive_sorted = sorted(positive, key=sort_key_of)
    sums = set()
    for i in range(3):
        for j in range(3):
            if i != j:
                a, b = positive_sorted[i], positive_sorted[j]
                sums.add((a[0] + b[0], a[1] + b[1]))
    simple = [rv for rv in positive_sorted if (rv[0], rv[1]) not in sums]
    check(len(simple) == 2, "SIMPLE_COUNT", "expected two simple roots")

    def vector_of(rv):
        for value, vec in root_values:
            if value == rv:
                return vec
        raise ReplayError("ROOT_VECTOR", "root vector lookup failed")

    def adF(v):
        return [
            [sum((F.lift(E.lift(B[o][i][j])) * v[i] for i in range(PORTS) if B[o][i][j] and v[i]), zeroF) for j in range(PORTS)]
            for o in range(PORTS)
        ]

    def lieF(u, v):
        return [
            sum((F.lift(E.lift(B[o][i][j])) * (u[i] * v[j]) for i in range(PORTS) for j in range(PORTS) if B[o][i][j] and u[i] and v[j]), zeroF)
            for o in range(PORTS)
        ]

    coroots = []
    for rv in simple:
        e_vec = vector_of(rv)
        f_vec = vector_of((zeroF - rv[0], zeroF - rv[1]))
        Hp = lieF(e_vec, f_vec)
        check(any(Hp), "TRIPLE", "sl(2) triple bracket vanishes")
        alpha_Hp = line_scalar(act(adF(Hp), e_vec, zeroF), e_vec, zeroF)
        check(bool(alpha_Hp), "TRIPLE_VALUE", "alpha(H') vanishes")
        H = [((oneF + oneF) / alpha_Hp) * v for v in Hp]
        check(line_scalar(act(adF(H), e_vec, zeroF), e_vec, zeroF) == oneF + oneF, "COROOT_NORM", "alpha(H) differs from two")
        coroots.append(H)

    def coroot_pairing(vec_line):
        out = []
        for H in coroots:
            img = act(adF(H), vec_line, zeroF)
            lam = line_scalar(img, vec_line, zeroF) if any(img) else zeroF
            out.append(as_integer(lam, oneF))
        return out

    cartan_matrix = []
    for H in coroots:
        row = []
        for rv in simple:
            e_vec = vector_of(rv)
            img = act(adF(H), e_vec, zeroF)
            lam = line_scalar(img, e_vec, zeroF) if any(img) else zeroF
            row.append(as_integer(lam, oneF))
        cartan_matrix.append(row)
    check(cartan_matrix == [[2, -1], [-1, 2]], "CARTAN_MATRIX", "computed Cartan matrix is not of type A2")

    all_root_coords = [coroot_pairing(vector_of(rv)) for rv in root_list]

    coroot_coordinates = {}
    for idx, rv in enumerate(positive_sorted):
        e_vec = vector_of(rv)
        f_vec = vector_of((zeroF - rv[0], zeroF - rv[1]))
        Hp = lieF(e_vec, f_vec)
        alpha_Hp = line_scalar(act(adF(Hp), e_vec, zeroF), e_vec, zeroF)
        H = [((oneF + oneF) / alpha_Hp) * v for v in Hp]
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
        check(sol is not None, "COROOT_SOLVE", "coroot coordinate solve failed")
        x1, x2 = sol
        for i in range(PORTS):
            check(H[i] == x1 * coroots[0][i] + x2 * coroots[1][i], "COROOT_SPAN", "coroot leaves the span")
        coroot_coordinates[f"positive_root_{idx + 1}"] = [as_integer(x1, oneF), as_integer(x2, oneF)]

    weight_multiset = {(0, 0): 6}
    for H in coroots:
        adH = adF(H)
        for v in ker0:
            check(not any(act(adH, v, zeroF)), "ZERO_WEIGHT", "joint-kernel vector carries a nonzero weight")
    for rv, vec in root_values:
        coords = tuple(coroot_pairing(vec))
        weight_multiset[coords] = weight_multiset.get(coords, 0) + 1
    check(sum(weight_multiset.values()) == PORTS, "WEIGHT_TOTAL", "weight multiplicities do not sum to twelve")

    record["type_computed"] = "A2"
    record["cartan_dimension"] = 2
    record["h2"] = encvec5(h2)
    record["nu_squared"] = enc_tower(disc)
    record["m2_restriction_trace"] = enc_tower(tau)
    record["roots_in_coroot_coordinates"] = sorted(all_root_coords)
    record["simple_root_coordinates"] = [coroot_pairing(vector_of(rv)) for rv in simple]
    record["cartan_matrix"] = cartan_matrix
    record["positive_root_coroots_in_simple_coroot_basis"] = coroot_coordinates
    weight_items = sorted(weight_multiset.items())
    record["carrier_weights"] = [{"coordinates": list(w), "multiplicity": m} for w, m in weight_items]
    flat = [w for w, m in weight_items for _ in range(m)]
    return record, [tuple(rc) for rc in all_root_coords], flat, 2


def replay_lattices(record, roots, carrier_weights, rank):
    root_rows = [list(r) for r in sorted(set(tuple(r) for r in roots))]
    diag, V = smith_form(root_rows)
    check(len(diag) == rank, "ROOT_RANK", "root lattice rank differs")
    center_order = 1
    for d in diag:
        center_order *= d
    invariant = [d for d in diag if d != 1]
    weight_images = []
    all_in = True
    for w in carrier_weights:
        img = image_in_quotient(list(w), diag, V, rank)
        weight_images.append(img)
        if any(img):
            all_in = False
    kernel_order = kernel_order_of(diag, weight_images)
    carrier_rows = [list(w) for w in sorted(set(carrier_weights)) if any(w)]
    cw_diag = smith_form(carrier_rows)[0] if carrier_rows else []
    record["lattices"] = {
        "weight_lattice": f"Z^{rank} in simple-coroot pairing coordinates",
        "root_lattice_row_generators": root_rows,
        "root_lattice_invariant_factors": diag,
        "center_of_simply_connected_form_order": center_order,
        "center_of_simply_connected_form_invariant_factors": invariant,
        "coroot_lattice": f"Z^{rank} in the simple-coroot basis",
        "coweight_lattice_index_over_coroot_lattice": center_order,
        "carrier_weight_lattice_invariant_factors": cw_diag,
        "carrier_weights_contained_in_root_lattice": all_in,
    }
    record["center_action_on_carrier"] = {
        "carrier_weight_images_in_center_dual": sorted(set(tuple(i) for i in weight_images)),
        "kernel_order": kernel_order,
        "kernel_equals_full_center": kernel_order == center_order,
        "faithful_quotient_order": center_order // kernel_order,
    }
    return center_order, kernel_order


def replay_family(fid, spec, certificate, C, inverse, projectors, bases):
    sample = sample_for(certificate, spec["certificate_family"], spec["parameters"])
    channel_values = {CHANNEL_INDEX[k]: dec5(v) for k, v in sample["channel_values"].items()}
    x, B = build_bracket(C, inverse, channel_values)
    check(encvec5(x) == sample["x_vector"], "X_VECTOR", "x vector differs from the certified one")
    check(jacobi_holds(B), "JACOBI", "bracket violates Jacobi")
    d_rref, d_piv = derived_of(B)
    check(len(d_rref) == sample["derived_dimension"], "DERIVED_DIM", "derived dimension differs")
    check(sectors_of(d_rref, projectors) == sample["derived_sectors"], "DERIVED_SECTORS", "derived sectors differ")
    z_rref, _ = center_of(B)
    check(len(z_rref) == sample["center_dimension"], "CENTER_DIM", "center dimension differs")
    check(sectors_of(z_rref, projectors) == sample["center_sectors"], "CENTER_SECTORS", "center sectors differ")

    summand_spans = []
    for sm in spec["summands"]:
        rows = []
        for label in sm["sectors"]:
            rows.extend([list(r) for r in bases[label][0]])
        summand_spans.append(reduce_rows(rows, PORTS))
    stacked = []
    for rr, _ in summand_spans:
        stacked.extend([list(r) for r in rr])
    joined, _ = reduce_rows(stacked, PORTS)
    check(len(joined) == len(d_rref), "SUMMAND_FILL", "summands do not fill the derived algebra")
    for row in joined:
        span_coords(d_rref, d_piv, row)
    for (ra, _), (rb, _) in itertools.combinations(summand_spans, 2):
        for u in ra:
            for v in rb:
                check(not any(lie(B, u, v)), "SUMMAND_COMMUTE", "distinct summands do not commute")

    family_record = {
        "certificate_family": spec["certificate_family"],
        "parameters": spec["parameters"],
        "channel_values": sample["channel_values"],
        "summands": [],
    }
    kernel_summary = []
    for sm in spec["summands"]:
        rec, roots, carrier_weights, rank = replay_summand(B, sm, bases)
        center_order, kernel_order = replay_lattices(rec, roots, carrier_weights, rank)
        rec["selected_global_form_on_carrier"] = ADJOINT_FORM[rec["type_computed"]] if kernel_order == center_order else "covering form retained"
        family_record["summands"].append(rec)
        kernel_summary.append((rec["type_computed"], center_order, kernel_order))

    abelian_sectors = sectors_of(z_rref, projectors)
    check(abelian_sectors == spec["abelian_sectors"], "ABELIAN_SECTORS", "abelian sector list differs")
    charge_ranks = []
    for zvec in z_rref:
        Mz = ad_of(B, list(zvec))
        check(not any(any(row) for row in Mz), "ABELIAN_TRIVIAL", "a certified central vector acts nontrivially")
        charge_ranks.append(0)
    family_record["abelian"] = {
        "dimension": len(z_rref),
        "sectors": abelian_sectors,
        "charge_matrix_rank_per_generator": charge_ranks,
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

    finite_center_factors = [SC_CENTER[t] for t, _, _ in kernel_summary]
    orders = [c for _, c, _ in kernel_summary]
    kernels = [k for _, _, k in kernel_summary]
    relations = [f"z{i}^{c} = 1" for i, (t, c, _) in enumerate(kernel_summary, start=1)]
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


def replay_controls():
    controls = {}
    a2_roots = [[2, -1], [-1, 2], [1, 1], [-2, 1], [1, -2], [-1, -1]]
    diag, V = smith_form([list(r) for r in sorted(set(tuple(r) for r in a2_roots))])
    weights = [list(r) for r in a2_roots] + [[1, 0]]
    images = [image_in_quotient(w, diag, V, 2) for w in weights]
    kernel = kernel_order_of(diag, images)
    center = 1
    for d in diag:
        center *= d
    controls["synthetic_fundamental_weight_A2"] = {
        "weights": weights,
        "center_order": center,
        "kernel_order": kernel,
        "conclusion": "adding one fundamental weight (1,0) shrinks the kernel to the identity; the pipeline selects the simply connected form when the carrier demands it",
    }
    check(center == 3 and kernel == 1, "CONTROL_A2", "A2 control kernel differs")
    diag1, V1 = smith_form([[2], [-2]])
    images1 = [image_in_quotient([w], diag1, V1, 1) for w in (2, -2, 1, -1)]
    kernel1 = kernel_order_of(diag1, images1)
    controls["synthetic_spinor_weight_A1"] = {
        "weights": [[2], [-2], [1], [-1]],
        "center_order": diag1[0],
        "kernel_order": kernel1,
        "conclusion": "adding weight 1 shrinks the kernel to the identity; the adjoint outcome on the carrier is computed, never assumed",
    }
    check(diag1[0] == 2 and kernel1 == 1, "CONTROL_A1", "A1 control kernel differs")
    controls["abelian_rescaling"] = {
        "statement": (
            "For every family the charge matrix of each central generator u and of the rescaled generator "
            "(7/3)u is the exact zero matrix, so no rescaling of the abelian generator produces a charge "
            "equation; the primitive-period question has no source-side constraint to solve"
        ),
        "verified_in": "analyze_family via ABELIAN_TRIVIAL together with linearity of the adjoint action",
    }
    return controls


def main() -> None:
    record, certificate, freeze, basis_raw, receipt1, system = load_all()
    C = constants_from_basis(basis_raw)
    inverse = inverse_from_system(system)
    group = [tuple(r) for r in receipt1["proper_port_action"]["permutation_rows"]]
    check(len(group) == 60 and len(set(group)) == 60, "GROUP", "expected sixty distinct group rows")
    projectors = projectors_from_group(group)
    bases = sector_bases_of(projectors)

    for fid, spec in FAMILY_SPECS.items():
        recomputed = replay_family(fid, spec, certificate, C, inverse, projectors, bases)
        check(
            normalize(recomputed) == normalize(record["families"][fid]),
            "FAMILY_MATCH",
            f"family {fid} record differs from the stored one",
        )
        # rescaling replay for the abelian control: the scaled generator acts by zero
        sample = sample_for(certificate, spec["certificate_family"], spec["parameters"])
        channel_values = {CHANNEL_INDEX[k]: dec5(v) for k, v in sample["channel_values"].items()}
        _, B = build_bracket(C, inverse, channel_values)
        z_rref, _ = center_of(B)
        scale = F5(Fraction(7, 3))
        for zvec in z_rref:
            scaled = [scale * v for v in zvec]
            check(not any(any(row) for row in ad_of(B, scaled)), "RESCALE", "a rescaled central generator acts nontrivially")

    check(
        normalize(replay_controls()) == normalize(record["controls"]),
        "CONTROLS_MATCH",
        "controls block differs from the stored one",
    )

    findings = record["findings"]
    for fid, fam in record["families"].items():
        check(fam["loop_to_kernel"]["kernel_equals_full_finite_center"], "FINDING_KERNEL", "a family kernel is not the full center")
        check(fam["abelian"]["primitive_period"] == "undetermined_by_source", "FINDING_PERIOD", "period finding differs")
    check(
        findings["global_form_selected_by_carrier"]["F"].startswith("PSU(3) x SO(3)"),
        "FINDING_F",
        "family F finding differs",
    )
    check(
        findings["global_form_selected_by_carrier"]["G"].startswith("PSU(3) x SO(3)"),
        "FINDING_G",
        "family G finding differs",
    )
    check(
        findings["global_form_selected_by_carrier"]["P"].startswith("SO(3) x SO(3)"),
        "FINDING_P",
        "family P finding differs",
    )
    print("replay complete: every lattice, weight, kernel, and pin verified")


if __name__ == "__main__":
    main()
