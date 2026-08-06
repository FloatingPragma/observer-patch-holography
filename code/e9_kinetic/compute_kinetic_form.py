#!/usr/bin/env python3
"""E9 sector 1: the gauge kinetic sector from the certified brackets (issue #716).

Semantic inputs are the pinned issue #705 compact-locus certificate
code/b14_jacobi/b14_compact_locus.certificate.json (with its pinned stage-one
Reynolds basis and stage-two reduction) and the pinned issue #707 lattice record
code/b16_lattices/lattices_v1.json. For each certified compact family (P:
so(3)+so(3); F and G: su(3)+so(3)) this producer rebuilds the exact bracket
tensor at the certified generic compact sample and computes, over Q(sqrt(5)),

1. the Killing form as the exact 12x12 matrix in the certified port basis, its
   exact decomposition as a combination of the four A5 sector projectors with
   one eigenvalue per sector, and the exact identity radical(K) = certified
   center;
2. the restriction of the Killing form to the gauge directions (the derived
   algebra) in the certified factor-concatenated basis, its exact negative
   definiteness, and the induced positive-definite gauge metric -K;
3. the Yang-Mills-shaped kinetic quadratic form F -> -(1/4) K(F, F): the exact
   kinetic matrix, its exact per-factor block split (the su(3) and so(3) blocks
   are Killing-orthogonal, verified entry by entry), and the exact statement
   that the space of ad-invariant symmetric forms on the gauge directions has
   dimension equal to the number of simple factors, so the source pins one
   coupling ray per simple factor and no scale;
4. coroot normalization receipts tying the Killing matrices to the issue #707
   coroot data: K(h1, h1) equals 12 mu^2 (A2) and 2 mu^2 (A1) exactly, and the
   carrier-weight second-moment Gram on the simple coroots equals 6 times the
   A2 Cartan matrix and 4 times the A1 Cartan matrix, all computed;
5. the structural comparison receipt: the split kinetic form has the Standard
   Model kinetic shape (one coupling per factor, no cross terms) as a
   structural statement, together with the structural delta at carrier level
   (two coupling rays su(3), so(3) against three Standard Model couplings; the
   certified center acts by zero on the carrier, so the source pins no abelian
   coupling ray). No coupling values and no measured data enter this record.

The producer fails closed on every pin and every exact identity. It selects no
coupling scale and does not close issue #716.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CERTIFICATE_PATH = REPO / "code/b14_jacobi/b14_compact_locus.certificate.json"
LATTICES_PATH = REPO / "code/b16_lattices/lattices_v1.json"
OUTPUT_PATH = HERE / "kinetic_form_v1.json"
VALIDATOR_PATH = HERE / "validate_kinetic_form.py"

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


class E9Error(RuntimeError):
    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise E9Error(code, message)


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


def vec_encode(vec) -> list[list[int]]:
    return [[i, *q5_encode(v)] for i, v in enumerate(vec) if v]


def vec_decode(entries, length) -> list[Q5]:
    out = [Q5_ZERO] * length
    for e in entries:
        out[e[0]] = q5_decode(e[1:])
    return out


def sym_encode(M) -> list[list[int]]:
    n = len(M)
    return [[i, j, *q5_encode(M[i][j])] for i in range(n) for j in range(i, n) if M[i][j]]


# --------------------------------------------------------- linear algebra
def rref(rows, width):
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


def coords_in_rref(rref_rows, pivots, vec, code="SPAN"):
    residual = list(vec)
    out = []
    for base, col in zip(rref_rows, pivots):
        c = residual[col]
        out.append(c)
        if c:
            residual = [x - c * y for x, y in zip(residual, base)]
    require(not any(residual), code, "vector leaves the claimed span")
    return out


def kernel_basis(M, width, zero, one):
    rr, piv = rref(M, width)
    pivset = set(piv)
    out = []
    for fc in (c for c in range(width) if c not in pivset):
        v = [zero] * width
        v[fc] = one
        for row, pc in zip(rr, piv):
            v[pc] = -row[fc]
        out.append(v)
    return out


def inertia(M):
    n = len(M)
    W = [row[:] for row in M]
    pos = neg = 0
    size = n
    while size:
        d = next((i for i in range(size) if W[i][i]), None)
        if d is None:
            hit = None
            for i in range(size):
                for j in range(i + 1, size):
                    if W[i][j]:
                        hit = (i, j)
                        break
                if hit:
                    break
            if hit is None:
                break
            i, j = hit
            for k in range(size):
                W[i][k] = W[i][k] + W[j][k]
            for k in range(size):
                W[k][i] = W[k][i] + W[k][j]
            d = i
        if d:
            W[0], W[d] = W[d], W[0]
            for row in W:
                row[0], row[d] = row[d], row[0]
        piv = W[0][0]
        if piv.sign() > 0:
            pos += 1
        else:
            neg += 1
        W = [[W[i][j] - (W[i][0] / piv) * W[0][j] for j in range(1, size)] for i in range(1, size)]
        size -= 1
    return [pos, neg, n - pos - neg]


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
    cert_body = dict(certificate)
    stored = cert_body.pop("certificate_sha256", None)
    require(stored == object_sha256(cert_body), "PIN_CERT_SELF", "certificate self-hash differs")
    up = certificate["upstream"]
    stage1_basis_path = REPO / up["stage1_basis_path"]
    stage1_receipt_path = REPO / up["stage1_receipt_path"]
    stage2_system_path = REPO / up["stage2_system_path"]
    stage2_receipt_path = REPO / up["stage2_receipt_path"]
    require(file_sha256(stage1_basis_path) == up["stage1_basis_file_sha256"], "PIN_STAGE1_BASIS", "stage-one basis hash differs")
    require(file_sha256(stage1_receipt_path) == up["stage1_receipt_file_sha256"], "PIN_STAGE1_RECEIPT", "stage-one receipt hash differs")
    require(file_sha256(stage2_system_path) == up["stage2_system_file_sha256"], "PIN_STAGE2_SYSTEM", "stage-two system hash differs")
    require(file_sha256(stage2_receipt_path) == up["stage2_receipt_file_sha256"], "PIN_STAGE2_RECEIPT", "stage-two receipt hash differs")
    lattices = json.loads(LATTICES_PATH.read_text(encoding="utf-8"))
    require(lattices["schema"] == "oph.b16_lattices.v1", "PIN_LATTICE_SCHEMA", "lattice record schema differs")
    lat_body = dict(lattices)
    lat_stored = lat_body.pop("result_sha256", None)
    require(lat_stored == object_sha256(lat_body), "PIN_LATTICE_SELF", "lattice record self-hash differs")
    prov = lattices["provenance"]
    require(prov["certificate_file_sha256"] == file_sha256(CERTIFICATE_PATH), "PIN_LATTICE_CERT_FILE", "lattice record pins a different certificate file")
    require(prov["certificate_self_sha256"] == certificate["certificate_sha256"], "PIN_LATTICE_CERT_SELF", "lattice record pins a different certificate body")
    basis_raw = json.loads(stage1_basis_path.read_text(encoding="utf-8"))
    stage1_receipt = json.loads(stage1_receipt_path.read_text(encoding="utf-8"))
    system = json.loads(stage2_system_path.read_text(encoding="utf-8"))
    return certificate, lattices, basis_raw, stage1_receipt, system


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
def spectral_projectors(group):
    unassigned = {(i, j) for i in range(PORTS) for j in range(PORTS)}
    orbits = []
    while unassigned:
        seed = min(unassigned)
        orbit = {(g[seed[0]], g[seed[1]]) for g in group}
        unassigned -= orbit
        orbits.append(orbit)
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


def jacobi_holds(B):
    for (i, j, k) in itertools.combinations(range(PORTS), 3):
        for o in range(PORTS):
            v = Q5_ZERO
            for m in range(PORTS):
                v = v + B[m][i][j] * B[o][m][k] + B[m][j][k] * B[o][m][i] + B[m][k][i] * B[o][m][j]
            if v:
                return False
    return True


def lie(B, u, v):
    return [
        sum((B[o][i][j] * (u[i] * v[j]) for i in range(PORTS) for j in range(PORTS) if B[o][i][j] and u[i] and v[j]), Q5_ZERO)
        for o in range(PORTS)
    ]


def ad_matrix(B, v):
    return [
        [sum((B[o][i][j] * v[i] for i in range(PORTS) if v[i]), Q5_ZERO) for j in range(PORTS)]
        for o in range(PORTS)
    ]


def killing_matrix(B):
    K = [[Q5_ZERO] * PORTS for _ in range(PORTS)]
    for i in range(PORTS):
        for j in range(i, PORTS):
            v = sum(
                (B[m][i][o] * B[o][j][m] for m in range(PORTS) for o in range(PORTS) if B[m][i][o] and B[o][j][m]),
                Q5_ZERO,
            )
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
    return rref(rows, PORTS)


def center_rref(B):
    rows = []
    for v in range(PORTS):
        flat = []
        for o in range(PORTS):
            for j in range(PORTS):
                flat.append(B[o][v][j])
        rows.append(flat)
    transposed = [[rows[i][j] for i in range(PORTS)] for j in range(len(rows[0]))]
    kernel = kernel_basis(transposed, PORTS, Q5_ZERO, Q5_ONE)
    return rref(kernel, PORTS) if kernel else ([], [])


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


def form_value(K, u, v):
    return sum(
        (u[i] * K[i][j] * v[j] for i in range(PORTS) for j in range(PORTS) if u[i] and K[i][j] and v[j]),
        Q5_ZERO,
    )


def gram_on_rows(K, rows):
    n = len(rows)
    G = [[Q5_ZERO] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            v = form_value(K, rows[i], rows[j])
            G[i][j] = v
            G[j][i] = v
    return G


# ----------------------------------------------------------- family layer
FAMILY_SPECS = {
    "P": {
        "certificate_family": "PLANE",
        "parameters": ["1", "1"],
        "factors": [
            {"name": "so3_plus", "sectors": ["three_plus"], "type": "A1"},
            {"name": "so3_minus", "sectors": ["three_minus"], "type": "A1"},
        ],
    },
    "F": {
        "certificate_family": "SU3PF",
        "parameters": ["1", "-1", "1"],
        "factors": [
            {"name": "su3", "sectors": ["three_plus", "five"], "type": "A2"},
            {"name": "so3", "sectors": ["three_minus"], "type": "A1"},
        ],
    },
    "G": {
        "certificate_family": "SU3MF",
        "parameters": ["1", "1", "1"],
        "factors": [
            {"name": "su3", "sectors": ["three_minus", "five"], "type": "A2"},
            {"name": "so3", "sectors": ["three_plus"], "type": "A1"},
        ],
    },
}

CARTAN_RATIO = {"A2": 12, "A1": 2}


def find_sample(certificate, family, parameters):
    for s in certificate["samples"]:
        if s["family"] == family and s["parameters"] == parameters:
            require(s["compact"] is True, "SAMPLE_COMPACT", "selected sample is not certified compact")
            return s
    raise E9Error("SAMPLE_MISSING", f"no certified sample for {family}{parameters}")


def coroot_receipts(K, B, projectors, factor_spec, summand):
    """Exact receipts tying the Killing matrix to the pinned coroot data."""
    h1 = vec_decode(summand["h1"], PORTS)
    mu2 = q5_decode(summand["mu_squared"])
    k_h1 = form_value(K, h1, h1)
    ratio = CARTAN_RATIO[factor_spec["type"]]
    require(k_h1 == Q5(ratio) * mu2, "COROOT_KH1", "K(h1, h1) differs from the spectral multiple of mu^2")
    # independent mu^2 receipt: tr(ad_h1^2 P_leading) = 2 mu^2 on the leading sector
    adh = ad_matrix(B, h1)
    P = projectors[factor_spec["sectors"][0]]
    tr = Q5_ZERO
    for i in range(PORTS):
        for j in range(PORTS):
            for k in range(PORTS):
                if P[i][j] and adh[j][k] and adh[k][i]:
                    tr = tr + P[i][j] * adh[j][k] * adh[k][i]
    require(tr == Q5(2) * mu2, "COROOT_MU2", "leading-sector trace of ad_h1^2 differs from 2 mu^2")
    receipts = {
        "killing_on_h1": q5_encode(k_h1),
        "mu_squared_from_lattice_record": q5_encode(mu2),
        "killing_on_h1_over_mu_squared": ratio,
        "leading_sector_trace_equals_two_mu_squared": True,
    }
    if factor_spec["type"] == "A2":
        h2 = vec_decode(summand["h2"], PORTS)
        require(not any(lie(B, h1, h2)), "COROOT_CARTAN", "[h1, h2] does not vanish")
        receipts["killing_on_h1_h2"] = q5_encode(form_value(K, h1, h2))
        receipts["killing_on_h2"] = q5_encode(form_value(K, h2, h2))
    # carrier-weight second moment on the simple coroots, from pinned integer data
    rank = summand["cartan_dimension"]
    gram = [[0] * rank for _ in range(rank)]
    for w in summand["carrier_weights"]:
        coords = w["coordinates"]
        mult = w["multiplicity"]
        for i in range(rank):
            for j in range(rank):
                gram[i][j] += mult * coords[i] * coords[j]
    cartan = summand["cartan_matrix"]
    require(cartan[0][0] != 0 and gram[0][0] % cartan[0][0] == 0, "GRAM_RATIO", "Gram ratio is not integral")
    gram_ratio = gram[0][0] // cartan[0][0]
    require(gram_ratio > 0, "GRAM_POSITIVE", "Gram ratio is not positive")
    require(
        all(gram[i][j] == gram_ratio * cartan[i][j] for i in range(rank) for j in range(rank)),
        "GRAM_CARTAN",
        "carrier-weight Gram is not the computed multiple of the Cartan matrix",
    )
    receipts["carrier_weight_gram_on_simple_coroots"] = gram
    receipts["gram_equals_ratio_times_cartan_matrix"] = True
    receipts["gram_ratio"] = gram_ratio
    return receipts


def invariant_form_space(B, derived_rows, derived_piv, projectors, factor_specs):
    """Dimension of ad-invariant symmetric forms on the derived algebra, and the
    invariance of the per-factor Killing blocks, all exact."""
    n = len(derived_rows)
    # ad matrices of the derived basis, in derived coordinates
    adm = []
    for bx in derived_rows:
        cols = []
        for by in derived_rows:
            w = lie(B, bx, by)
            cols.append(coords_in_rref(derived_rows, derived_piv, w, "DERIVED_CLOSED"))
        adm.append(cols)  # adm[x][j][i]: coefficient of e_i in [b_x, b_j]
    pair_index = {}
    count = 0
    for i in range(n):
        for j in range(i, n):
            pair_index[(i, j)] = count
            count += 1
    rows = []
    for Mx in adm:
        for u in range(n):
            for v in range(u, n):
                row = [Q5_ZERO] * count
                for k in range(n):
                    c = Mx[u][k]
                    if c:
                        a, b = min(k, v), max(k, v)
                        row[pair_index[(a, b)]] = row[pair_index[(a, b)]] + c
                    c2 = Mx[v][k]
                    if c2:
                        a, b = min(u, k), max(u, k)
                        row[pair_index[(a, b)]] = row[pair_index[(a, b)]] + c2
                if any(row):
                    rows.append(row)
    ech, _ = rref(rows, count)
    dimension = count - len(ech)
    # invariance of each per-factor Killing block, checked in derived coordinates:
    # with Ax[i][j] the matrix of ad(b_x), the condition is Ax^T Qb + Qb Ax = 0
    K = killing_matrix(B)
    for spec in factor_specs:
        Pf = [[sum((projectors[label][i][j] for label in spec["sectors"]), Q5_ZERO) for j in range(PORTS)] for i in range(PORTS)]
        projected = [
            [sum((Pf[o][i] * vec[i] for i in range(PORTS) if vec[i]), Q5_ZERO) for o in range(PORTS)]
            for vec in derived_rows
        ]
        Qb = [[form_value(K, projected[u], projected[v]) for v in range(n)] for u in range(n)]
        for Mx in adm:
            Ax = [[Mx[j][i] for j in range(n)] for i in range(n)]
            for u in range(n):
                for v in range(n):
                    val = sum((Ax[k][u] * Qb[k][v] + Qb[u][k] * Ax[k][v] for k in range(n) if Ax[k][u] or Qb[u][k]), Q5_ZERO)
                    require(not val, "BLOCK_INVARIANT", "a per-factor Killing block fails ad-invariance")
    return dimension


def build_family(fid, spec, certificate, lattices, C, inverse, projectors, bases):
    sample = find_sample(certificate, spec["certificate_family"], spec["parameters"])
    lat_family = lattices["families"][fid]
    require(lat_family["certificate_family"] == spec["certificate_family"], "LAT_FAMILY", "lattice record family differs")
    require(lat_family["parameters"] == spec["parameters"], "LAT_PARAMS", "lattice record parameters differ")
    require(lat_family["channel_values"] == sample["channel_values"], "LAT_CHANNELS", "lattice record channels differ")
    channel_values = {CHANNEL_INDEX[k]: q5_decode(v) for k, v in sample["channel_values"].items()}
    x, B = bracket_tensor(C, inverse, channel_values)
    require(vec_encode(x) == sample["x_vector"], "X_VECTOR", "recomputed x vector differs from the certified one")
    require(jacobi_holds(B), "JACOBI", "bracket violates Jacobi")
    d_rref, d_piv = derived_span(B)
    require(len(d_rref) == sample["derived_dimension"], "DERIVED_DIM", "derived dimension differs")
    require(sector_membership(d_rref, projectors) == sample["derived_sectors"], "DERIVED_SECTORS", "derived sectors differ")
    z_rref, _ = center_rref(B)
    require(len(z_rref) == sample["center_dimension"], "CENTER_DIM", "center dimension differs")
    require(sector_membership(z_rref, projectors) == sample["center_sectors"], "CENTER_SECTORS", "center sectors differ")

    K = killing_matrix(B)

    # sector eigenvalues: K equals the exact projector combination
    eigenvalues = {}
    for label in SECTORS:
        tr = Q5_ZERO
        for i in range(PORTS):
            for k in range(PORTS):
                if K[i][k] and projectors[label][k][i]:
                    tr = tr + K[i][k] * projectors[label][k][i]
        eigenvalues[label] = tr / Q5(SECTOR_DIMS[label])
    for i in range(PORTS):
        for j in range(PORTS):
            combo = sum((eigenvalues[label] * projectors[label][i][j] for label in SECTORS), Q5_ZERO)
            require(K[i][j] == combo, "KILLING_PROJECTOR_FORM", "Killing matrix leaves the projector span")

    # radical of K equals the certified center
    radical = kernel_basis([list(row) for row in K], PORTS, Q5_ZERO, Q5_ONE)
    rad_rref, _ = rref(radical, PORTS) if radical else ([], [])
    require(len(rad_rref) == len(z_rref), "RADICAL_DIM", "radical dimension differs from the center dimension")
    require(rad_rref == z_rref, "RADICAL_CENTER", "radical of the Killing form differs from the certified center")

    # per-factor bases: union of the sector reduced-row-echelon bases
    factor_rows = []
    factor_records = []
    for factor_spec, summand in zip(spec["factors"], lat_family["summands"]):
        require(summand["name"] == factor_spec["name"], "LAT_SUMMAND_NAME", "lattice summand name differs")
        require(summand["sectors"] == factor_spec["sectors"], "LAT_SUMMAND_SECTORS", "lattice summand sectors differ")
        require(summand["type_computed"] == factor_spec["type"], "LAT_SUMMAND_TYPE", "lattice summand type differs")
        rows = []
        for label in factor_spec["sectors"]:
            rows.extend([list(r) for r in bases[label][0]])
        dim = len(rows)
        require(dim == summand["dimension"], "FACTOR_DIM", "factor dimension differs from the lattice record")
        # subalgebra closure inside the factor
        f_rref, f_piv = rref([list(r) for r in rows], PORTS)
        require(len(f_rref) == dim, "FACTOR_RANK", "factor basis is dependent")
        for u in rows:
            for v in rows:
                w = lie(B, u, v)
                if any(w):
                    coords_in_rref(f_rref, f_piv, w, "FACTOR_CLOSED")
        # membership in the derived algebra
        for u in rows:
            coords_in_rref(d_rref, d_piv, u, "FACTOR_IN_DERIVED")
        Kf = gram_on_rows(K, rows)
        inert = inertia([row[:] for row in Kf])
        require(inert == [0, dim, 0], "FACTOR_NEGDEF", "factor Killing block is not negative definite")
        minus = [[-Kf[i][j] for j in range(dim)] for i in range(dim)]
        require(inertia([row[:] for row in minus]) == [dim, 0, 0], "FACTOR_POSDEF", "factor gauge metric is not positive definite")
        kinetic = [[Q5(Fraction(-1, 4)) * Kf[i][j] for j in range(dim)] for i in range(dim)]
        sector_eigs = {label: q5_encode(eigenvalues[label]) for label in factor_spec["sectors"]}
        factor_records.append({
            "name": factor_spec["name"],
            "sectors": factor_spec["sectors"],
            "dimension": dim,
            "type_from_lattice_record": factor_spec["type"],
            "cartan_matrix_from_lattice_record": summand["cartan_matrix"],
            "basis_port_vectors": [vec_encode(r) for r in rows],
            "killing_matrix_upper_entries": sym_encode(Kf),
            "killing_inertia": [0, dim, 0],
            "killing_negative_definite": True,
            "gauge_metric_inertia": [dim, 0, 0],
            "sector_eigenvalues_on_factor": sector_eigs,
            "kinetic_block_upper_entries": sym_encode(kinetic),
            "coupling_ray": (
                "lambda times minus the Killing block, lambda > 0; the source pins the ray and pins no scale"
            ),
            "coroot_receipts": coroot_receipts(K, B, projectors, factor_spec, summand),
        })
        factor_rows.append(rows)

    # pairwise commuting and exact Killing orthogonality across factors
    ortho_pairs = []
    for (ia, ra), (ib, rb) in itertools.combinations(enumerate(factor_rows), 2):
        for u in ra:
            for v in rb:
                require(not any(lie(B, u, v)), "FACTOR_COMMUTE", "distinct factors do not commute")
                require(not form_value(K, u, v), "FACTOR_ORTHO", "distinct factors are not Killing orthogonal")
        ortho_pairs.append({
            "left": spec["factors"][ia]["name"],
            "right": spec["factors"][ib]["name"],
            "cross_block_zero": True,
        })
    for rows in factor_rows:
        for u in rows:
            for z in z_rref:
                require(not form_value(K, u, list(z)), "CENTER_ORTHO", "a factor direction pairs with the center")

    # derived-basis Killing and kinetic matrices on the factor concatenation
    stacked = [r for rows in factor_rows for r in rows]
    s_rref, _ = rref([list(r) for r in stacked], PORTS)
    require(len(s_rref) == len(d_rref), "CONCAT_RANK", "factor concatenation does not span the derived algebra")
    ranges = []
    start = 0
    for rows in factor_rows:
        ranges.append([start, start + len(rows)])
        start += len(rows)
    K_der = gram_on_rows(K, stacked)
    n = len(stacked)
    for i in range(n):
        for j in range(n):
            inside = any(a <= i < b and a <= j < b for a, b in ranges)
            if not inside:
                require(not K_der[i][j], "OFF_BLOCK", "off-block Killing entry is nonzero")
    require(inertia([row[:] for row in K_der]) == [0, n, 0], "DERIVED_NEGDEF", "derived Killing block is not negative definite")
    minus_der = [[-K_der[i][j] for j in range(n)] for i in range(n)]
    require(inertia([row[:] for row in minus_der]) == [n, 0, 0], "DERIVED_POSDEF", "gauge metric is not positive definite")
    kinetic_der = [[Q5(Fraction(-1, 4)) * K_der[i][j] for j in range(n)] for i in range(n)]

    inv_dim = invariant_form_space(B, d_rref, d_piv, projectors, spec["factors"])
    require(inv_dim == len(spec["factors"]), "INVARIANT_DIM", "invariant form space dimension differs from the factor count")

    return {
        "certificate_family": spec["certificate_family"],
        "parameters": spec["parameters"],
        "channel_values": sample["channel_values"],
        "killing_port_basis": {
            "matrix_upper_entries": sym_encode(K),
            "sector_eigenvalues": {label: q5_encode(eigenvalues[label]) for label in SECTORS},
            "equals_projector_combination": True,
            "radical_dimension": len(rad_rref),
            "radical_equals_certified_center": True,
        },
        "factors": factor_records,
        "factor_orthogonality": {
            "cross_blocks_zero": ortho_pairs,
            "factors_orthogonal_to_center": True,
        },
        "derived_algebra": {
            "dimension": n,
            "sectors": sample["derived_sectors"],
            "basis": "per-factor sector bases concatenated in the recorded factor order",
            "block_index_ranges": ranges,
            "kinetic_form_matrix_upper_entries": sym_encode(kinetic_der),
            "off_block_entries_zero": True,
            "killing_inertia": [0, n, 0],
            "gauge_metric_inertia": [n, 0, 0],
            "kinetic_form_positive_definite": True,
        },
        "invariant_form_space": {
            "dimension_on_gauge_directions": inv_dim,
            "equals_simple_factor_count": True,
            "factor_killing_blocks_are_invariant": True,
            "consequence": (
                "every ad-invariant symmetric bilinear form on the gauge directions is a per-factor "
                "multiple of the Killing blocks, so the kinetic term splits with exactly one coupling "
                "ray per simple factor and zero cross terms"
            ),
        },
    }


def shape_comparison(families, lattices):
    per_family = {}
    for fid, fam in families.items():
        types = [f["type_from_lattice_record"] for f in fam["factors"]]
        names = [f["name"] for f in fam["factors"]]
        abelian = lattices["families"][fid]["abelian"]
        delta = {
            "coupling_rays_here": len(fam["factors"]),
            "standard_model_couplings": 3,
            "absent_at_carrier_level": (
                ["su(3) coupling ray", "u(1) abelian coupling ray"] if "A2" not in types else ["u(1) abelian coupling ray"]
            ),
            "abelian_reason_from_source": (
                "the certified central generators act by the exact zero matrix on the carrier and the "
                "charge lattice is the zero lattice (pinned issue #707 record), so the source pins no "
                "abelian kinetic ray at carrier level"
            ),
            "abelian_charge_lattice_rank": abelian["charge_lattice_rank"],
        }
        if "A2" in types:
            delta["su2_identification"] = (
                "the A1 factor carries the su(2) Lie algebra; the carrier-selected global form is "
                "SO(3) = SU(2)/(Z/2) per the pinned issue #707 record"
            )
        per_family[fid] = {
            "factor_names": names,
            "factor_types": types,
            "coupling_ray_count": len(fam["factors"]),
            "block_diagonal_no_cross_terms": True,
            "one_coupling_ray_per_factor": True,
            "matches_standard_model_kinetic_shape": True,
            "structural_delta": delta,
        }
    return {
        "standard_model_reference_shape": {
            "factors": ["su(3)", "su(2)", "u(1)"],
            "coupling_count": 3,
            "kinetic_structure": "block diagonal, one coupling per factor, zero cross terms",
            "role": "shape reference only; no coupling values and no measured data enter this record",
        },
        "per_family": per_family,
        "finding": (
            "Families F and G realize the two-coupling structure su(3) + so(3): the Standard Model "
            "kinetic shape with one coupling per simple factor and zero cross terms, without an abelian "
            "factor at carrier level. Family P realizes so(3) + so(3). The delta against the three "
            "Standard Model couplings is the absent abelian ray (and for P the absent su(3) ray)."
        ),
    }


SCALE_PINNING = {
    "pinned_by_source": (
        "per simple factor, the ray of ad-invariant positive quadratic forms "
        "{lambda * (minus the Killing block) : lambda > 0}; the invariant-form space on the gauge "
        "directions has dimension exactly equal to the simple-factor count, so the block split and "
        "the rays are forced"
    ),
    "distinguished_representative": (
        "the ambient Killing form selects the representative tuple (minus each Killing block) inside "
        "the product of rays; whether the physical normalization coincides with this tuple is part of "
        "the residual"
    ),
    "unpinned_by_source": "one positive scale per simple factor (two per family)",
    "what_would_pin_the_scales": (
        "a source-selected normalization: a source-derived functional that selects one form out of "
        "each ray, equivalently one positive number per simple factor"
    ),
    "residual": {"id": "E9_SECTOR1_COUPLING_SCALES", "status": "open"},
}


# ----------------------------------------------------------------------- main
def main() -> dict:
    certificate, lattices, basis_raw, stage1_receipt, system = load_pinned()
    C = structure_constants(basis_raw)
    inverse = pinned_inverse_transform(system)
    group = [tuple(r) for r in stage1_receipt["proper_port_action"]["permutation_rows"]]
    require(len(group) == 60 and len(set(group)) == 60, "GROUP", "expected sixty distinct group rows")
    projectors = spectral_projectors(group)
    bases = sector_rref_bases(projectors)

    families = {}
    for fid, spec in FAMILY_SPECS.items():
        families[fid] = build_family(fid, spec, certificate, lattices, C, inverse, projectors, bases)

    result = {
        "schema": "oph.e9_kinetic.v1",
        "issue": 716,
        "sector": 1,
        "provenance": {
            "certificate_path": "code/b14_jacobi/b14_compact_locus.certificate.json",
            "certificate_file_sha256": file_sha256(CERTIFICATE_PATH),
            "certificate_self_sha256": certificate["certificate_sha256"],
            "lattices_path": "code/b16_lattices/lattices_v1.json",
            "lattices_file_sha256": file_sha256(LATTICES_PATH),
            "lattices_result_sha256": lattices["result_sha256"],
            "upstream": certificate["upstream"],
        },
        "conventions": {
            "field": "Q(sqrt(5)); scalars encode as [n1, d1, n2, d2] meaning n1/d1 + (n2/d2)*sqrt(5)",
            "certified_basis": (
                "the twelve-port coordinate basis of the carrier; sector subspaces carry the reduced "
                "row echelon bases of the A5 sector projectors in the pinned port order"
            ),
            "matrix_encoding": (
                "symmetric matrices list the nonzero upper-triangle entries as [i, j, n1, d1, n2, d2]"
            ),
            "vector_encoding": "port vectors list the nonzero entries as [i, n1, d1, n2, d2]",
            "kinetic_shape": (
                "F -> -(1/4) K(F, F) on the gauge directions; every recorded kinetic matrix equals "
                "minus one quarter of the corresponding Killing matrix"
            ),
            "inertia": "signatures are recorded as [positive, negative, zero] counts of an exact symmetric reduction",
        },
        "families": families,
        "shape_comparison": shape_comparison(families, lattices),
        "scale_pinning": SCALE_PINNING,
        "claim_boundary": (
            "This record computes, for the three certified compact families at their certified generic "
            "samples, the exact Killing form in the certified port basis, its sector-projector "
            "decomposition, its restriction to the gauge directions, the per-factor split of the "
            "Yang-Mills-shaped kinetic quadratic form, and the coupling-ray receipts, all over "
            "Q(sqrt(5)). The Standard Model appears only as a kinetic shape reference; no coupling "
            "values and no measured data enter. The record selects no coupling scale, constructs no "
            "holonomy, spacetime action, or continuum bundle, and does not close issue #716."
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
        names = [f["name"] for f in fam["factors"]]
        dims = [f["dimension"] for f in fam["factors"]]
        print(fid, "->", " + ".join(f"{n}[{d}]" for n, d in zip(names, dims)),
              "| gauge metric positive definite on", fam["derived_algebra"]["dimension"], "directions")
