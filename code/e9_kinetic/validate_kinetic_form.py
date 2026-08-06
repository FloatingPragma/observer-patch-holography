#!/usr/bin/env python3
"""Independent replay of the E9 sector-1 kinetic record kinetic_form_v1.json.

This validator reads the pinned issue #705 certificate, its pinned stage-one
and stage-two inputs, the pinned issue #707 lattice record, and the emitted
kinetic record. With replay arithmetic written for this file it

* verifies every provenance pin fail-closed: the certificate self-hash, the
  upstream file hashes, the lattice-record self-hash and its cross pins, the
  producer hash, its own hash, and the record self-hash;
* rebuilds the exact bracket tensor at each certified generic compact sample
  and recomputes the Killing matrices, sector eigenvalues, radical, factor
  split, orthogonality, definiteness, kinetic matrices, coroot receipts, and
  the invariant-form dimension;
* checks the stored record twice: targeted exact checks on the stored data
  (kinetic scaling, inertia, off-block zeros, projector form, delta content)
  and canonical-JSON equality of every recomputed block with the stored one.

Any mismatch raises ReplayError. The flag --tamper-audit reruns the check
suite against six tampered in-memory variants of the record and requires a
rejection for each, printing the rejecting code.
"""

from __future__ import annotations

import copy
import hashlib
import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CERTIFICATE_PATH = REPO / "code/b14_jacobi/b14_compact_locus.certificate.json"
LATTICES_PATH = REPO / "code/b16_lattices/lattices_v1.json"
RECORD_PATH = HERE / "kinetic_form_v1.json"
PRODUCER_PATH = HERE / "compute_kinetic_form.py"

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


ZERO = F5()
ONE = F5(1)


def enc5(v: F5) -> list[int]:
    return [v.a.numerator, v.a.denominator, v.b.numerator, v.b.denominator]


def dec5(entry) -> F5:
    return F5(Fraction(entry[0], entry[1]), Fraction(entry[2], entry[3]))


def encvec(vec) -> list[list[int]]:
    return [[i, *enc5(v)] for i, v in enumerate(vec) if v]


def decvec(entries, length) -> list[F5]:
    out = [ZERO] * length
    for e in entries:
        out[e[0]] = dec5(e[1:])
    return out


def enc_sym(M) -> list[list[int]]:
    n = len(M)
    return [[i, j, *enc5(M[i][j])] for i in range(n) for j in range(i, n) if M[i][j]]


def dec_sym(entries, n):
    M = [[ZERO] * n for _ in range(n)]
    for e in entries:
        i, j = e[0], e[1]
        v = dec5(e[2:])
        M[i][j] = v
        M[j][i] = v
    return M


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


def span_coords(rows, pivots, vec, code="SPAN"):
    residual = list(vec)
    out = []
    for base, col in zip(rows, pivots):
        c = residual[col]
        out.append(c)
        if c:
            residual = [x - c * y for x, y in zip(residual, base)]
    check(not any(residual), code, "vector leaves the claimed span")
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


def sym_inertia(M):
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
def canonical(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def sha_object(value) -> str:
    return "sha256:" + hashlib.sha256(canonical(value)).hexdigest()


def sha_file(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def normalize(value):
    return json.loads(json.dumps(value, sort_keys=True))


# ------------------------------------------------------------ pinned input
def load_env():
    certificate = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
    cert_body = dict(certificate)
    cert_stored = cert_body.pop("certificate_sha256", None)
    check(cert_stored == sha_object(cert_body), "CERT_SELF", "certificate self-hash differs")
    up = certificate["upstream"]
    stage1_basis_path = REPO / up["stage1_basis_path"]
    stage1_receipt_path = REPO / up["stage1_receipt_path"]
    stage2_system_path = REPO / up["stage2_system_path"]
    stage2_receipt_path = REPO / up["stage2_receipt_path"]
    check(sha_file(stage1_basis_path) == up["stage1_basis_file_sha256"], "PIN_STAGE1_BASIS", "stage-one basis hash differs")
    check(sha_file(stage1_receipt_path) == up["stage1_receipt_file_sha256"], "PIN_STAGE1_RECEIPT", "stage-one receipt hash differs")
    check(sha_file(stage2_system_path) == up["stage2_system_file_sha256"], "PIN_STAGE2_SYSTEM", "stage-two system hash differs")
    check(sha_file(stage2_receipt_path) == up["stage2_receipt_file_sha256"], "PIN_STAGE2_RECEIPT", "stage-two receipt hash differs")
    lattices = json.loads(LATTICES_PATH.read_text(encoding="utf-8"))
    check(lattices["schema"] == "oph.b16_lattices.v1", "LATTICE_SCHEMA", "lattice record schema differs")
    lat_body = dict(lattices)
    lat_stored = lat_body.pop("result_sha256", None)
    check(lat_stored == sha_object(lat_body), "LATTICE_SELF", "lattice record self-hash differs")
    prov = lattices["provenance"]
    check(prov["certificate_file_sha256"] == sha_file(CERTIFICATE_PATH), "LATTICE_CERT_FILE", "lattice record pins a different certificate file")
    check(prov["certificate_self_sha256"] == certificate["certificate_sha256"], "LATTICE_CERT_SELF", "lattice record pins a different certificate body")
    basis_raw = json.loads(stage1_basis_path.read_text(encoding="utf-8"))
    receipt1 = json.loads(stage1_receipt_path.read_text(encoding="utf-8"))
    system = json.loads(stage2_system_path.read_text(encoding="utf-8"))
    return {
        "certificate": certificate,
        "lattices": lattices,
        "basis_raw": basis_raw,
        "receipt1": receipt1,
        "system": system,
        "certificate_file_sha256": sha_file(CERTIFICATE_PATH),
        "lattices_file_sha256": sha_file(LATTICES_PATH),
    }


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
    inverse = [[ZERO] * PARAMS for _ in range(PARAMS)]
    for i, entries in enumerate(cdm["inverse_transform_rows"]):
        for j, n1, d1, n2, d2 in entries:
            inverse[i][j] = F5(Fraction(n1, d1), Fraction(n2, d2))
    return inverse


def projectors_from_group(group):
    unassigned = {(i, j) for i in range(PORTS) for j in range(PORTS)}
    orbits = []
    while unassigned:
        seed = min(unassigned)
        orbit = {(g[seed[0]], g[seed[1]]) for g in group}
        unassigned -= orbit
        orbits.append(orbit)
    check(sorted(len(o) for o in orbits) == [12, 12, 60, 60], "PAIR_ORBITS", "unexpected ordered-pair orbit sizes")
    candidates = sorted((o for o in orbits if len(o) == 60), key=lambda o: tuple(sorted(o)))
    orbital = candidates[0]
    A = [[ZERO] * PORTS for _ in range(PORTS)]
    for i, j in orbital:
        A[i][j] = ONE
    R5 = F5(0, 1)
    eigenvalues = {"fixed": F5(5), "three_plus": R5, "three_minus": -R5, "five": F5(-1)}
    projectors = {}
    for label, ev in eigenvalues.items():
        M = [[ONE if i == j else ZERO for j in range(PORTS)] for i in range(PORTS)]
        for label2, ev2 in eigenvalues.items():
            if label2 == label:
                continue
            factor = [[A[i][j] - (ev2 if i == j else ZERO) for j in range(PORTS)] for i in range(PORTS)]
            product = [[sum((M[i][k] * factor[k][j] for k in range(PORTS)), ZERO) for j in range(PORTS)] for i in range(PORTS)]
            M = [[v * (ev - ev2).inv() for v in row] for row in product]
        projectors[label] = M
    for label in SECTORS:
        P = projectors[label]
        PP = [[sum((P[i][k] * P[k][j] for k in range(PORTS)), ZERO) for j in range(PORTS)] for i in range(PORTS)]
        check(PP == P, "PROJECTOR_IDEMPOTENT", f"{label} projector is not idempotent")
        trace = sum((P[i][i] for i in range(PORTS)), ZERO)
        check(trace == F5(SECTOR_DIMS[label]), "PROJECTOR_TRACE", f"{label} projector trace differs")
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
    x = [sum((inverse[i][j] * v for j, v in channel_values.items()), ZERO) for i in range(PARAMS)]
    B = [[[ZERO] * PORTS for _ in range(PORTS)] for _ in range(PORTS)]
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


def jacobi_ok(B):
    for (i, j, k) in itertools.combinations(range(PORTS), 3):
        for o in range(PORTS):
            v = ZERO
            for m in range(PORTS):
                v = v + B[m][i][j] * B[o][m][k] + B[m][j][k] * B[o][m][i] + B[m][k][i] * B[o][m][j]
            if v:
                return False
    return True


def lie_of(B, u, v):
    return [
        sum((B[o][i][j] * (u[i] * v[j]) for i in range(PORTS) for j in range(PORTS) if B[o][i][j] and u[i] and v[j]), ZERO)
        for o in range(PORTS)
    ]


def ad_of(B, v):
    return [
        [sum((B[o][i][j] * v[i] for i in range(PORTS) if v[i]), ZERO) for j in range(PORTS)]
        for o in range(PORTS)
    ]


def killing_of(B):
    K = [[ZERO] * PORTS for _ in range(PORTS)]
    for i in range(PORTS):
        for j in range(i, PORTS):
            v = sum(
                (B[m][i][o] * B[o][j][m] for m in range(PORTS) for o in range(PORTS) if B[m][i][o] and B[o][j][m]),
                ZERO,
            )
            K[i][j] = v
            K[j][i] = v
    return K


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
    kernel = null_space(transposed, PORTS, ZERO, ONE)
    return reduce_rows(kernel, PORTS) if kernel else ([], [])


def sectors_of(rows, projectors):
    labels = []
    for label in SECTORS:
        P = projectors[label]
        for row in rows:
            image = [sum((P[o][i] * row[i] for i in range(PORTS) if row[i]), ZERO) for o in range(PORTS)]
            if any(image):
                labels.append(label)
                break
    return labels


def form_val(K, u, v):
    return sum(
        (u[i] * K[i][j] * v[j] for i in range(PORTS) for j in range(PORTS) if u[i] and K[i][j] and v[j]),
        ZERO,
    )


def gram_rows(K, rows):
    n = len(rows)
    G = [[ZERO] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            v = form_val(K, rows[i], rows[j])
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


def sample_for(certificate, family, parameters):
    for s in certificate["samples"]:
        if s["family"] == family and s["parameters"] == parameters:
            check(s["compact"] is True, "SAMPLE_COMPACT", "selected sample is not certified compact")
            return s
    raise ReplayError("SAMPLE_MISSING", f"no certified sample for {family}{parameters}")


def replay_coroot_receipts(K, B, projectors, factor_spec, summand):
    h1 = decvec(summand["h1"], PORTS)
    mu2 = dec5(summand["mu_squared"])
    k_h1 = form_val(K, h1, h1)
    ratio = CARTAN_RATIO[factor_spec["type"]]
    check(k_h1 == F5(ratio) * mu2, "COROOT_KH1", "K(h1, h1) differs from the spectral multiple of mu^2")
    adh = ad_of(B, h1)
    P = projectors[factor_spec["sectors"][0]]
    tr = ZERO
    for i in range(PORTS):
        for j in range(PORTS):
            for k in range(PORTS):
                if P[i][j] and adh[j][k] and adh[k][i]:
                    tr = tr + P[i][j] * adh[j][k] * adh[k][i]
    check(tr == F5(2) * mu2, "COROOT_MU2", "leading-sector trace of ad_h1^2 differs from 2 mu^2")
    receipts = {
        "killing_on_h1": enc5(k_h1),
        "mu_squared_from_lattice_record": enc5(mu2),
        "killing_on_h1_over_mu_squared": ratio,
        "leading_sector_trace_equals_two_mu_squared": True,
    }
    if factor_spec["type"] == "A2":
        h2 = decvec(summand["h2"], PORTS)
        check(not any(lie_of(B, h1, h2)), "COROOT_CARTAN", "[h1, h2] does not vanish")
        receipts["killing_on_h1_h2"] = enc5(form_val(K, h1, h2))
        receipts["killing_on_h2"] = enc5(form_val(K, h2, h2))
    rank = summand["cartan_dimension"]
    gram = [[0] * rank for _ in range(rank)]
    for w in summand["carrier_weights"]:
        coords = w["coordinates"]
        mult = w["multiplicity"]
        for i in range(rank):
            for j in range(rank):
                gram[i][j] += mult * coords[i] * coords[j]
    cartan = summand["cartan_matrix"]
    check(cartan[0][0] != 0 and gram[0][0] % cartan[0][0] == 0, "GRAM_RATIO", "Gram ratio is not integral")
    gram_ratio = gram[0][0] // cartan[0][0]
    check(gram_ratio > 0, "GRAM_POSITIVE", "Gram ratio is not positive")
    check(
        all(gram[i][j] == gram_ratio * cartan[i][j] for i in range(rank) for j in range(rank)),
        "GRAM_CARTAN",
        "carrier-weight Gram is not the computed multiple of the Cartan matrix",
    )
    receipts["carrier_weight_gram_on_simple_coroots"] = gram
    receipts["gram_equals_ratio_times_cartan_matrix"] = True
    receipts["gram_ratio"] = gram_ratio
    return receipts


def replay_invariant_dimension(B, derived_rows, derived_piv, projectors, factor_specs):
    n = len(derived_rows)
    adm = []
    for bx in derived_rows:
        cols = []
        for by in derived_rows:
            w = lie_of(B, bx, by)
            cols.append(span_coords(derived_rows, derived_piv, w, "DERIVED_CLOSED"))
        adm.append(cols)
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
                row = [ZERO] * count
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
    ech, _ = reduce_rows(rows, count)
    dimension = count - len(ech)
    K = killing_of(B)
    for spec in factor_specs:
        Pf = [[sum((projectors[label][i][j] for label in spec["sectors"]), ZERO) for j in range(PORTS)] for i in range(PORTS)]
        projected = [
            [sum((Pf[o][i] * vec[i] for i in range(PORTS) if vec[i]), ZERO) for o in range(PORTS)]
            for vec in derived_rows
        ]
        Qb = [[form_val(K, projected[u], projected[v]) for v in range(n)] for u in range(n)]
        for Mx in adm:
            Ax = [[Mx[j][i] for j in range(n)] for i in range(n)]
            for u in range(n):
                for v in range(n):
                    val = sum((Ax[k][u] * Qb[k][v] + Qb[u][k] * Ax[k][v] for k in range(n) if Ax[k][u] or Qb[u][k]), ZERO)
                    check(not val, "BLOCK_INVARIANT", "a per-factor Killing block fails ad-invariance")
    return dimension


def replay_family(fid, spec, env, C, inverse, projectors, bases):
    certificate = env["certificate"]
    lattices = env["lattices"]
    sample = sample_for(certificate, spec["certificate_family"], spec["parameters"])
    lat_family = lattices["families"][fid]
    check(lat_family["certificate_family"] == spec["certificate_family"], "LAT_FAMILY", "lattice record family differs")
    check(lat_family["parameters"] == spec["parameters"], "LAT_PARAMS", "lattice record parameters differ")
    check(lat_family["channel_values"] == sample["channel_values"], "LAT_CHANNELS", "lattice record channels differ")
    channel_values = {CHANNEL_INDEX[k]: dec5(v) for k, v in sample["channel_values"].items()}
    x, B = build_bracket(C, inverse, channel_values)
    check(encvec(x) == sample["x_vector"], "X_VECTOR", "x vector differs from the certified one")
    check(jacobi_ok(B), "JACOBI", "bracket violates Jacobi")
    d_rref, d_piv = derived_of(B)
    check(len(d_rref) == sample["derived_dimension"], "DERIVED_DIM", "derived dimension differs")
    check(sectors_of(d_rref, projectors) == sample["derived_sectors"], "DERIVED_SECTORS", "derived sectors differ")
    z_rref, _ = center_of(B)
    check(len(z_rref) == sample["center_dimension"], "CENTER_DIM", "center dimension differs")
    check(sectors_of(z_rref, projectors) == sample["center_sectors"], "CENTER_SECTORS", "center sectors differ")

    K = killing_of(B)
    eigenvalues = {}
    for label in SECTORS:
        tr = ZERO
        for i in range(PORTS):
            for k in range(PORTS):
                if K[i][k] and projectors[label][k][i]:
                    tr = tr + K[i][k] * projectors[label][k][i]
        eigenvalues[label] = tr / F5(SECTOR_DIMS[label])
    for i in range(PORTS):
        for j in range(PORTS):
            combo = sum((eigenvalues[label] * projectors[label][i][j] for label in SECTORS), ZERO)
            check(K[i][j] == combo, "KILLING_PROJECTOR_FORM", "Killing matrix leaves the projector span")

    radical = null_space([list(row) for row in K], PORTS, ZERO, ONE)
    rad_rref, _ = reduce_rows(radical, PORTS) if radical else ([], [])
    check(len(rad_rref) == len(z_rref), "RADICAL_DIM", "radical dimension differs from the center dimension")
    check(rad_rref == z_rref, "RADICAL_CENTER", "radical of the Killing form differs from the certified center")

    factor_rows = []
    factor_records = []
    for factor_spec, summand in zip(spec["factors"], lat_family["summands"]):
        check(summand["name"] == factor_spec["name"], "LAT_SUMMAND_NAME", "lattice summand name differs")
        check(summand["sectors"] == factor_spec["sectors"], "LAT_SUMMAND_SECTORS", "lattice summand sectors differ")
        check(summand["type_computed"] == factor_spec["type"], "LAT_SUMMAND_TYPE", "lattice summand type differs")
        rows = []
        for label in factor_spec["sectors"]:
            rows.extend([list(r) for r in bases[label][0]])
        dim = len(rows)
        check(dim == summand["dimension"], "FACTOR_DIM", "factor dimension differs from the lattice record")
        f_rref, f_piv = reduce_rows([list(r) for r in rows], PORTS)
        check(len(f_rref) == dim, "FACTOR_RANK", "factor basis is dependent")
        for u in rows:
            for v in rows:
                w = lie_of(B, u, v)
                if any(w):
                    span_coords(f_rref, f_piv, w, "FACTOR_CLOSED")
        for u in rows:
            span_coords(d_rref, d_piv, u, "FACTOR_IN_DERIVED")
        Kf = gram_rows(K, rows)
        check(sym_inertia([row[:] for row in Kf]) == [0, dim, 0], "FACTOR_NEGDEF", "factor Killing block is not negative definite")
        minus = [[-Kf[i][j] for j in range(dim)] for i in range(dim)]
        check(sym_inertia([row[:] for row in minus]) == [dim, 0, 0], "FACTOR_POSDEF", "factor gauge metric is not positive definite")
        kinetic = [[F5(Fraction(-1, 4)) * Kf[i][j] for j in range(dim)] for i in range(dim)]
        sector_eigs = {label: enc5(eigenvalues[label]) for label in factor_spec["sectors"]}
        factor_records.append({
            "name": factor_spec["name"],
            "sectors": factor_spec["sectors"],
            "dimension": dim,
            "type_from_lattice_record": factor_spec["type"],
            "cartan_matrix_from_lattice_record": summand["cartan_matrix"],
            "basis_port_vectors": [encvec(r) for r in rows],
            "killing_matrix_upper_entries": enc_sym(Kf),
            "killing_inertia": [0, dim, 0],
            "killing_negative_definite": True,
            "gauge_metric_inertia": [dim, 0, 0],
            "sector_eigenvalues_on_factor": sector_eigs,
            "kinetic_block_upper_entries": enc_sym(kinetic),
            "coupling_ray": (
                "lambda times minus the Killing block, lambda > 0; the source pins the ray and pins no scale"
            ),
            "coroot_receipts": replay_coroot_receipts(K, B, projectors, factor_spec, summand),
        })
        factor_rows.append(rows)

    ortho_pairs = []
    for (ia, ra), (ib, rb) in itertools.combinations(enumerate(factor_rows), 2):
        for u in ra:
            for v in rb:
                check(not any(lie_of(B, u, v)), "FACTOR_COMMUTE", "distinct factors do not commute")
                check(not form_val(K, u, v), "FACTOR_ORTHO", "distinct factors are not Killing orthogonal")
        ortho_pairs.append({
            "left": spec["factors"][ia]["name"],
            "right": spec["factors"][ib]["name"],
            "cross_block_zero": True,
        })
    for rows in factor_rows:
        for u in rows:
            for z in z_rref:
                check(not form_val(K, u, list(z)), "CENTER_ORTHO", "a factor direction pairs with the center")

    stacked = [r for rows in factor_rows for r in rows]
    s_rref, _ = reduce_rows([list(r) for r in stacked], PORTS)
    check(len(s_rref) == len(d_rref), "CONCAT_RANK", "factor concatenation does not span the derived algebra")
    ranges = []
    start = 0
    for rows in factor_rows:
        ranges.append([start, start + len(rows)])
        start += len(rows)
    K_der = gram_rows(K, stacked)
    n = len(stacked)
    for i in range(n):
        for j in range(n):
            inside = any(a <= i < b and a <= j < b for a, b in ranges)
            if not inside:
                check(not K_der[i][j], "OFF_BLOCK", "off-block Killing entry is nonzero")
    check(sym_inertia([row[:] for row in K_der]) == [0, n, 0], "DERIVED_NEGDEF", "derived Killing block is not negative definite")
    minus_der = [[-K_der[i][j] for j in range(n)] for i in range(n)]
    check(sym_inertia([row[:] for row in minus_der]) == [n, 0, 0], "DERIVED_POSDEF", "gauge metric is not positive definite")
    kinetic_der = [[F5(Fraction(-1, 4)) * K_der[i][j] for j in range(n)] for i in range(n)]

    inv_dim = replay_invariant_dimension(B, d_rref, d_piv, projectors, spec["factors"])
    check(inv_dim == len(spec["factors"]), "INVARIANT_DIM", "invariant form space dimension differs from the factor count")

    return {
        "certificate_family": spec["certificate_family"],
        "parameters": spec["parameters"],
        "channel_values": sample["channel_values"],
        "killing_port_basis": {
            "matrix_upper_entries": enc_sym(K),
            "sector_eigenvalues": {label: enc5(eigenvalues[label]) for label in SECTORS},
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
            "kinetic_form_matrix_upper_entries": enc_sym(kinetic_der),
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


def replay_shape(families, lattices):
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


EXPECTED_SCALE_PINNING = {
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

EXPECTED_CONVENTIONS = {
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
}

EXPECTED_CLAIM = (
    "This record computes, for the three certified compact families at their certified generic "
    "samples, the exact Killing form in the certified port basis, its sector-projector "
    "decomposition, its restriction to the gauge directions, the per-factor split of the "
    "Yang-Mills-shaped kinetic quadratic form, and the coupling-ray receipts, all over "
    "Q(sqrt(5)). The Standard Model appears only as a kinetic shape reference; no coupling "
    "values and no measured data enter. The record selects no coupling scale, constructs no "
    "holonomy, spacetime action, or continuum bundle, and does not close issue #716."
)


# ------------------------------------------------------- stored-data checks
def stored_checks(record, projectors):
    """Targeted exact checks on the stored record, independent of the replay
    equality: kinetic scaling, inertia, off-block zeros, projector form, and
    delta content."""
    quarter = F5(Fraction(-1, 4))
    for fid, fam in record["families"].items():
        for factor in fam["factors"]:
            dim = factor["dimension"]
            Kf = dec_sym(factor["killing_matrix_upper_entries"], dim)
            Qf = dec_sym(factor["kinetic_block_upper_entries"], dim)
            for i in range(dim):
                for j in range(dim):
                    check(Qf[i][j] == quarter * Kf[i][j], "KINETIC_SCALING", f"{fid} {factor['name']} kinetic block is not -(1/4) times the Killing block")
            check(sym_inertia([row[:] for row in Kf]) == factor["killing_inertia"], "STORED_INERTIA", f"{fid} {factor['name']} stored Killing inertia differs")
            check(factor["killing_inertia"] == [0, dim, 0], "STORED_NEGDEF", f"{fid} {factor['name']} stored inertia is not negative definite")
            minus = [[-Kf[i][j] for j in range(dim)] for i in range(dim)]
            check(sym_inertia([row[:] for row in minus]) == factor["gauge_metric_inertia"], "STORED_GAUGE", f"{fid} {factor['name']} stored gauge inertia differs")
            check(factor["gauge_metric_inertia"] == [dim, 0, 0], "STORED_POSDEF", f"{fid} {factor['name']} stored gauge inertia is not positive definite")
        der = fam["derived_algebra"]
        n = der["dimension"]
        Qd = dec_sym(der["kinetic_form_matrix_upper_entries"], n)
        Kd = [[F5(-4) * Qd[i][j] for j in range(n)] for i in range(n)]
        ranges = der["block_index_ranges"]
        for i in range(n):
            for j in range(n):
                inside = any(a <= i < b and a <= j < b for a, b in ranges)
                if not inside:
                    check(not Kd[i][j], "STORED_OFFBLOCK", f"{fid} stored derived matrix has a nonzero off-block entry")
        check(sym_inertia([row[:] for row in Kd]) == der["killing_inertia"], "STORED_DER_INERTIA", f"{fid} stored derived Killing inertia differs")
        check(der["killing_inertia"] == [0, n, 0], "STORED_DER_NEGDEF", f"{fid} stored derived inertia is not negative definite")
        minus_der = [[-Kd[i][j] for j in range(n)] for i in range(n)]
        check(sym_inertia([row[:] for row in minus_der]) == der["gauge_metric_inertia"], "STORED_DER_GAUGE", f"{fid} stored derived gauge inertia differs")
        check(der["gauge_metric_inertia"] == [n, 0, 0], "STORED_DER_POSDEF", f"{fid} stored derived gauge inertia is not positive definite")
        kp = fam["killing_port_basis"]
        K = dec_sym(kp["matrix_upper_entries"], PORTS)
        eigen = {label: dec5(kp["sector_eigenvalues"][label]) for label in SECTORS}
        for i in range(PORTS):
            for j in range(PORTS):
                combo = sum((eigen[label] * projectors[label][i][j] for label in SECTORS), ZERO)
                check(K[i][j] == combo, "STORED_PROJ_FORM", f"{fid} stored Killing matrix leaves the stored projector combination")
    for fid, shape in record["shape_comparison"]["per_family"].items():
        fam = record["families"][fid]
        check(shape["coupling_ray_count"] == len(fam["factors"]), "STORED_RAY_COUNT", f"{fid} coupling ray count differs from the factor count")
        delta = shape["structural_delta"]
        check(delta["coupling_rays_here"] == len(fam["factors"]), "STORED_DELTA_COUNT", f"{fid} delta ray count differs")
        check("u(1) abelian coupling ray" in delta["absent_at_carrier_level"], "STORED_DELTA", f"{fid} delta does not record the absent abelian ray")


def run_checks(record, env, expected_families, expected_shape, projectors):
    check(record.get("schema") == "oph.e9_kinetic.v1", "SCHEMA", "record schema differs")
    check(record.get("issue") == 716, "ISSUE", "record issue differs")
    check(record.get("sector") == 1, "SECTOR", "record sector differs")
    body = dict(record)
    stored = body.pop("result_sha256", None)
    check(stored == sha_object(body), "RECORD_SELF", "record self-hash differs")

    prov = record["provenance"]
    check(prov["certificate_file_sha256"] == env["certificate_file_sha256"], "PROV_CERT_FILE", "record pins a different certificate file")
    check(prov["certificate_self_sha256"] == env["certificate"]["certificate_sha256"], "PROV_CERT_SELF", "record pins a different certificate body")
    check(prov["lattices_file_sha256"] == env["lattices_file_sha256"], "PROV_LATTICE_FILE", "record pins a different lattice file")
    check(prov["lattices_result_sha256"] == env["lattices"]["result_sha256"], "PROV_LATTICE_SELF", "record pins a different lattice body")
    check(normalize(prov["upstream"]) == normalize(env["certificate"]["upstream"]), "PROV_UPSTREAM", "record upstream block differs")
    pins = record["implementation_pins"]
    check(pins["producer_sha256"] == sha_file(PRODUCER_PATH), "IMPL_PRODUCER", "producer hash differs")
    check(pins["validator_sha256"] == sha_file(Path(__file__).resolve()), "IMPL_VALIDATOR", "validator hash differs")

    stored_checks(record, projectors)

    for fid in FAMILY_SPECS:
        check(
            normalize(record["families"][fid]) == normalize(expected_families[fid]),
            f"FAMILY_MATCH_{fid}",
            f"family {fid} record differs from the replay",
        )
    check(normalize(record["shape_comparison"]) == normalize(expected_shape), "SHAPE_MATCH", "shape comparison block differs from the replay")
    check(normalize(record["scale_pinning"]) == normalize(EXPECTED_SCALE_PINNING), "SCALE_MATCH", "scale pinning block differs from the replay")
    check(normalize(record["conventions"]) == normalize(EXPECTED_CONVENTIONS), "CONVENTIONS_MATCH", "conventions block differs from the replay")
    check(record["claim_boundary"] == EXPECTED_CLAIM, "CLAIM_MATCH", "claim boundary differs from the replay")


# ------------------------------------------------------------ tamper audit
def rehash(record):
    body = dict(record)
    body.pop("result_sha256", None)
    record["result_sha256"] = sha_object(body)
    return record


def tamper_variants(record):
    variants = []

    t = copy.deepcopy(record)
    entry = t["families"]["F"]["factors"][0]["killing_matrix_upper_entries"][0]
    entry[2] = -entry[2]
    variants.append(("killing_entry_sign_flip", rehash(t)))

    t = copy.deepcopy(record)
    t["families"]["F"]["derived_algebra"]["kinetic_form_matrix_upper_entries"].append([0, 8, 1, 1, 0, 1])
    variants.append(("off_block_entry_injected", rehash(t)))

    t = copy.deepcopy(record)
    t["shape_comparison"]["per_family"]["F"]["structural_delta"]["absent_at_carrier_level"] = []
    variants.append(("structural_delta_forged", rehash(t)))

    t = copy.deepcopy(record)
    t["result_sha256"] = "sha256:" + "0" * 64
    variants.append(("self_hash_tampered", t))

    t = copy.deepcopy(record)
    t["provenance"]["certificate_file_sha256"] = "sha256:" + "f" * 64
    variants.append(("provenance_pin_tampered", rehash(t)))

    t = copy.deepcopy(record)
    t["families"]["P"]["derived_algebra"]["killing_inertia"] = [1, 5, 0]
    variants.append(("inertia_claim_forged", rehash(t)))

    return variants


def main() -> None:
    env = load_env()
    C = constants_from_basis(env["basis_raw"])
    inverse = inverse_from_system(env["system"])
    group = [tuple(r) for r in env["receipt1"]["proper_port_action"]["permutation_rows"]]
    check(len(group) == 60 and len(set(group)) == 60, "GROUP", "expected sixty distinct group rows")
    projectors = projectors_from_group(group)
    bases = sector_bases_of(projectors)

    expected_families = {}
    for fid, spec in FAMILY_SPECS.items():
        expected_families[fid] = replay_family(fid, spec, env, C, inverse, projectors, bases)
    expected_shape = replay_shape(expected_families, env["lattices"])

    record = json.loads(RECORD_PATH.read_text(encoding="utf-8"))
    run_checks(record, env, expected_families, expected_shape, projectors)
    print("replay complete: every Killing matrix, split, receipt, and pin verified")

    if "--tamper-audit" in sys.argv[1:]:
        rejected = 0
        for name, variant in tamper_variants(record):
            try:
                run_checks(variant, env, expected_families, expected_shape, projectors)
            except ReplayError as err:
                rejected += 1
                print(f"tamper rejected: {name} -> {err.code}")
            else:
                print(f"tamper ACCEPTED: {name}")
                raise SystemExit(1)
        check(rejected >= 3, "TAMPER_COUNT", "fewer than three tamper rejections")
        print(f"tamper audit complete: {rejected} of {rejected} variants rejected")


if __name__ == "__main__":
    main()
