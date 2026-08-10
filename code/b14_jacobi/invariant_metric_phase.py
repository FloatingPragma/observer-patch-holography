#!/usr/bin/env python3
"""Exact carrier-induced invariant-metric phase diagram for issue B14 / #705.

The pinned twelve-port carrier splits under the proper port action into four
pairwise inequivalent real irreducible sectors 1 + 3 + 3' + 5 (multiplicity
free).  Every port-action-invariant inner product on the carrier is therefore a
positive sector-scale metric g = alpha*P1 + beta*P3 + gamma*P3' + delta*P5, and
every such metric functorially induces a Hilbert--Schmidt metric on the
structure-constant space Lambda^2 V* (x) V.  This producer computes the exact
squared distances from the pinned oriented-face bracket to the classified
compact families P, F, G as Laurent polynomials in the sector scales and
certifies the resulting phase diagram:

  1. COMPLETENESS   the commutant of the proper port action has Q-dimension 4
                    and is spanned by the four symmetric spectral projectors,
                    so the sector-scale cone is the complete family of
                    invariant carrier inner products;
  2. DIAGONALITY    the fourteen channel lines are pairwise orthogonal under
                    every induced metric (all sector-block cross pairings
                    vanish identically);
  3. CLOSED FORMS   d_P^2, d_F^2, d_G^2 are explicit three-term Laurent
                    polynomials in (beta, gamma, delta); alpha never appears;
  4. P-EXCLUSION    d_P^2 - d_G^2 and d_P^2 - d_F^2 are coefficient-positive,
                    so P is strictly excluded for every invariant metric;
  5. BALANCED SLICE for every sector-balanced metric (beta = gamma) the family
                    G is strictly nearest;
  6. PHASE BOX      for every metric with beta/delta in [1/50, 6] the family G
                    is strictly nearest for all gamma and delta;
  7. F REGION       the F family becomes nearest at the exact witness
                    (beta, gamma, delta) = (8, 1, 1); the tie surface is the
                    explicit curve h(gamma/delta) = k(beta/delta);
  8. GALOIS DUALITY the F and G coefficient tables are the sqrt(5) |-> -sqrt(5)
                    conjugates of one another with the sector swap
                    beta <-> gamma;
  9. CONTROL        a channel-diagonal invariant metric on the bracket space
                    that is NOT carrier-induced (it breaks the forced equality
                    of the normalized t_pp_to_p and t_pf_to_f weights) reverses
                    the balanced-point selection, so carrier-inducedness is
                    load-bearing.

BOUNDARY: the comparison ranges over the compact locus classified by the
pinned b14_compact_locus certificate, conditional on its three named textbook
lemmas.  The nearest-point repair rule and the restriction to carrier-induced
metrics are declared choices of the discriminator, not source selections.  The
sector-balance property is a property of a metric, not a derived source law.
No physical current, laboratory gauge field, or source bracket selection is
claimed.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
MANIFEST = REPO / "code/a5_closure/manifests/echosahedral_federation_reference.json"
BASIS = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json"
RECEIPT = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json"
COMPACT = HERE / "b14_compact_locus.certificate.json"
SELECTOR = HERE / "oriented_face_bracket_selector.certificate.json"
OUTPUT = HERE / "invariant_metric_phase.certificate.json"
VERIFIER = HERE / "verify_invariant_metric_phase.py"

N = 12
M = 14
SECTORS = ["fixed", "three_plus", "three_minus", "five"]
SECTOR_DIMS = {"fixed": 1, "three_plus": 3, "three_minus": 3, "five": 5}
CHANNELS = [
    "d_plus", "d_minus", "d_five", "t_pp_to_p", "t_mm_to_m",
    "t_ff_to_p", "t_ff_to_m", "t_pm_to_f", "t_pf_to_p",
    "t_pf_to_m", "t_pf_to_f", "t_mf_to_p", "t_mf_to_m",
    "t_mf_to_f",
]
CH = {name: i for i, name in enumerate(CHANNELS)}
FACE_CHANNELS = ["t_pp_to_p", "t_mm_to_m", "t_ff_to_p", "t_ff_to_m",
                 "t_pf_to_f", "t_mf_to_f"]


class CertificateError(RuntimeError):
    pass


def require(value: bool, message: str) -> None:
    if not value:
        raise CertificateError(message)


class Q5:
    """Exact element a + b sqrt(5)."""

    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = a if isinstance(a, Fraction) else Fraction(a)
        self.b = b if isinstance(b, Fraction) else Fraction(b)

    def __add__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return Q5(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q5(-self.a, -self.b)

    def __sub__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return Q5(self.a * other.a + 5 * self.b * other.b,
                  self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def reciprocal(self):
        norm = self.a * self.a - 5 * self.b * self.b
        require(norm != 0, "division by zero in Q(sqrt(5))")
        return Q5(self.a / norm, -self.b / norm)

    def __truediv__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return self * other.reciprocal()

    def __eq__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return self.a == other.a and self.b == other.b

    def __hash__(self):
        return hash((self.a, self.b))

    def __bool__(self):
        return bool(self.a or self.b)

    def conj(self):
        return Q5(self.a, -self.b)

    def __repr__(self):
        return f"Q5({self.a},{self.b})"


ZERO = Q5()
ONE = Q5(1)
ROOT5 = Q5(0, 1)


def sign(x: Q5) -> int:
    if not x:
        return 0
    if x.a == 0:
        return 1 if x.b > 0 else -1
    if x.b == 0:
        return 1 if x.a > 0 else -1
    if (x.a > 0) == (x.b > 0):
        return 1 if x.a > 0 else -1
    comparison = x.a * x.a - 5 * x.b * x.b
    require(comparison != 0, "unexpected zero in Q(sqrt(5)) sign test")
    if x.a > 0:
        return 1 if comparison > 0 else -1
    return 1 if comparison < 0 else -1


def enc(x: Q5) -> list[int]:
    return [x.a.numerator, x.a.denominator, x.b.numerator, x.b.denominator]


def canonical_bytes(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def object_hash(value) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_hash(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


# --------------------------------------------------------------- linear algebra
def matmul(X, Y):
    n, m, p = len(X), len(Y), len(Y[0])
    Z = [[ZERO] * p for _ in range(n)]
    for i in range(n):
        for k in range(m):
            if X[i][k]:
                xik = X[i][k]
                for j in range(p):
                    if Y[k][j]:
                        Z[i][j] = Z[i][j] + xik * Y[k][j]
    return Z


def identity(n):
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def transpose(X):
    return [list(row) for row in zip(*X)]


def rational_rank(rows):
    """Rank over Q of a matrix given as list of Fraction rows."""
    rows = [list(r) for r in rows if any(r)]
    rank = 0
    cols = len(rows[0]) if rows else 0
    pivot_col = 0
    while rows and pivot_col < cols:
        pivot_row = next((i for i, r in enumerate(rows) if r[pivot_col]), None)
        if pivot_row is None:
            pivot_col += 1
            continue
        rows[0], rows[pivot_row] = rows[pivot_row], rows[0]
        head = rows.pop(0)
        inv = 1 / head[pivot_col]
        head = [x * inv for x in head]
        rows = [[x - r[pivot_col] * h for x, h in zip(r, head)]
                for r in rows]
        rows = [r for r in rows if any(r)]
        rank += 1
        pivot_col += 1
    return rank


# ------------------------------------------------------------- pinned loading
def load_inputs():
    manifest = json.loads(MANIFEST.read_text())
    basis_raw = json.loads(BASIS.read_text())
    receipt = json.loads(RECEIPT.read_text())
    compact = json.loads(COMPACT.read_text())
    selector = json.loads(SELECTOR.read_text())
    return manifest, basis_raw, receipt, compact, selector


def group_rows(receipt):
    group = [tuple(row) for row in receipt["proper_port_action"]["permutation_rows"]]
    require(len(group) == 60 and len(set(group)) == 60, "expected sixty distinct group rows")
    for row in group:
        require(sorted(row) == list(range(N)), "group row is not a permutation")
    return group


def pair_orbits(group):
    unassigned = {(i, j) for i in range(N) for j in range(N)}
    orbits = []
    while unassigned:
        seed = min(unassigned)
        orbit = {(g[seed[0]], g[seed[1]]) for g in group}
        unassigned -= orbit
        orbits.append(orbit)
    return orbits


def spectral_projectors(group):
    orbits = pair_orbits(group)
    require(sorted(len(o) for o in orbits) == [12, 12, 60, 60], "unexpected ordered-pair orbit sizes")
    candidates = sorted((o for o in orbits if len(o) == 60), key=lambda o: tuple(sorted(o)))
    orbital = candidates[0]
    require(all((j, i) in orbital for (i, j) in orbital), "orbital adjacency is not symmetric")
    A = [[ZERO] * N for _ in range(N)]
    for i, j in orbital:
        A[i][j] = ONE
    eigenvalues = {"fixed": Q5(5), "three_plus": ROOT5, "three_minus": -ROOT5, "five": Q5(-1)}
    projectors = {}
    for label, ev in eigenvalues.items():
        Mx = identity(N)
        for label2, ev2 in eigenvalues.items():
            if label2 == label:
                continue
            factor = [[A[i][j] - (ev2 if i == j else ZERO) for j in range(N)] for i in range(N)]
            scale = (ev - ev2).reciprocal()
            Mx = [[v * scale for v in row] for row in matmul(Mx, factor)]
        projectors[label] = Mx
    total = [[sum((projectors[l][i][j] for l in SECTORS), ZERO) for j in range(N)] for i in range(N)]
    require(total == identity(N), "spectral projectors do not sum to identity")
    for label in SECTORS:
        P = projectors[label]
        require(matmul(P, P) == P, f"{label} projector is not idempotent")
        require(transpose(P) == P, f"{label} projector is not symmetric")
        trace = sum((P[i][i] for i in range(N)), ZERO)
        require(trace == Q5(SECTOR_DIMS[label]), f"{label} projector trace differs")
    for a in SECTORS:
        for b in SECTORS:
            if a < b:
                prod = matmul(projectors[a], projectors[b])
                require(all(not v for row in prod for v in row),
                        f"projectors {a} and {b} are not orthogonal")
    return orbital, projectors


def commutant_dimension(group):
    """Q-dimension of {X : X rho(g) = rho(g) X for all g}, computed exactly.

    The constraint for a permutation g is X[i][g(j)] = X[g^-1(i)][j]; we build
    the rational solution space dimension via the orbit count of the induced
    action on matrix positions (Burnside-free direct method: positions (i,j)
    and (g(i), g(j)) must carry equal entries, so the dimension equals the
    number of orbits of the diagonal action on ordered pairs)."""
    orbits = pair_orbits(group)
    return len(orbits)


def structure_constants(basis_raw):
    require([row["basis_id"] for row in basis_raw["basis"]] == [f"R{i:02d}" for i in range(M)],
            "unexpected Reynolds basis order")
    C = [[[[Fraction(0)] * N for _ in range(N)] for _ in range(N)] for _ in range(M)]
    for a, row in enumerate(basis_raw["basis"]):
        for out, left, right, numerator, denominator in row["entries"]:
            require(left < right, "noncanonical Reynolds entry")
            value = Fraction(numerator, denominator)
            C[a][out][left][right] = value
            C[a][out][right][left] = -value
    return C


CHANNEL_SPECS = [
    ("t_pp_to_p", "three_plus", "three_plus", "three_plus"),
    ("t_mm_to_m", "three_minus", "three_minus", "three_minus"),
    ("t_ff_to_p", "five", "five", "three_plus"),
    ("t_ff_to_m", "five", "five", "three_minus"),
    ("t_pm_to_f", "three_plus", "three_minus", "five"),
    ("t_pf_to_p", "three_plus", "five", "three_plus"),
    ("t_pf_to_m", "three_plus", "five", "three_minus"),
    ("t_pf_to_f", "three_plus", "five", "five"),
    ("t_mf_to_p", "three_minus", "five", "three_plus"),
    ("t_mf_to_m", "three_minus", "five", "three_minus"),
    ("t_mf_to_f", "three_minus", "five", "five"),
]


def channel_transform(vectors, C, projectors):
    """The pinned channel transform convention of the compact-locus producer:
    three derivation eigenform rows followed by eleven unnormalized
    first-nonzero coordinate forms in CHANNEL_SPECS order."""
    import itertools

    derivations = []
    for a in range(M):
        Mx = [[Q5(sum((C[a][o][s][c] for s in range(N)), Fraction(0)))
               for c in range(N)] for o in range(N)]
        derivations.append(Mx)
    eigenforms = {}
    for label in SECTORS:
        P = projectors[label]
        form = []
        for Mx in derivations:
            trace = ZERO
            for i in range(N):
                for j in range(N):
                    if P[i][j] and Mx[j][i]:
                        trace = trace + P[i][j] * Mx[j][i]
            form.append(trace / Q5(SECTOR_DIMS[label]))
        for a, Mx in enumerate(derivations):
            left = matmul(Mx, P)
            right = [[form[a] * v for v in row] for row in P]
            require(left == right, f"fixed-line adjoint is not scalar on {label}")
        eigenforms[label] = form
    require(not any(eigenforms["fixed"]), "fixed-line adjoint does not kill the fixed line")

    def coordinate_form(Pl, Pr, Po, coordinate):
        output, li, ri = coordinate
        form = []
        for vector in vectors:
            value = ZERO
            for (so, p, q), coefficient in vector.items():
                of = Po[output][so]
                if not of:
                    continue
                inf = Pl[p][li] * Pr[q][ri] - Pl[q][li] * Pr[p][ri]
                if inf:
                    value = value + of * inf * Q5(coefficient)
            form.append(value)
        return form

    transform = [eigenforms["three_plus"], eigenforms["three_minus"], eigenforms["five"]]
    for _, left, right, output in CHANNEL_SPECS:
        found = None
        for coordinate in itertools.product(range(N), repeat=3):
            form = coordinate_form(projectors[left], projectors[right],
                                   projectors[output], coordinate)
            if any(form):
                found = form
                break
        require(found is not None, "channel has no nonzero coordinate")
        transform.append(found)
    # invert the 14x14 transform over Q(sqrt5)
    aug = [list(row) + [ONE if i == j else ZERO for j in range(M)]
           for i, row in enumerate(transform)]
    for col in range(M):
        pivot = next((r for r in range(col, M) if aug[r][col]), None)
        require(pivot is not None, "channel transform is singular")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        inv = aug[col][col].reciprocal()
        aug[col] = [inv * x for x in aug[col]]
        for r in range(M):
            if r != col and aug[r][col]:
                f = aug[r][col]
                aug[r] = [x - f * y for x, y in zip(aug[r], aug[col])]
    inverse = [row[M:] for row in aug]
    require(matmul(transform, inverse) == identity(M), "inverse is not exact")
    return transform, inverse


# ------------------------------------------------------------ metric machinery
def q5_tensor(C, coeffs):
    """Assemble sum_a coeffs[a] * R_a as a full Q5 tensor T[o][i][j]."""
    T = [[[ZERO] * N for _ in range(N)] for _ in range(N)]
    for a, c in enumerate(coeffs):
        if not c:
            continue
        for o in range(N):
            for i in range(N):
                row = C[a][o][i]
                for j in range(N):
                    if row[j]:
                        T[o][i][j] = T[o][i][j] + c * row[j]
    return T


def sector_blocks(projectors, T, S):
    """M[(so,su,sv)] = sum over full index range of
    T^o_ij S^p_kl P_so[o][p] P_su[i][k] P_sv[j][l]; the induced pairing under
    sector scales lam is (1/2) sum lam_so lam_su^-1 lam_sv^-1 M[...]."""
    blocks = {}
    P = projectors
    for so in SECTORS:
        Pso = P[so]
        A = [[[ZERO] * N for _ in range(N)] for _ in range(N)]  # A[p][i][j]
        for o in range(N):
            for i in range(N):
                for j in range(N):
                    v = T[o][i][j]
                    if v:
                        for p in range(N):
                            if Pso[o][p]:
                                A[p][i][j] = A[p][i][j] + v * Pso[o][p]
        for su in SECTORS:
            Psu = P[su]
            B = [[[ZERO] * N for _ in range(N)] for _ in range(N)]  # B[p][k][j]
            for p in range(N):
                for i in range(N):
                    for j in range(N):
                        v = A[p][i][j]
                        if v:
                            for k in range(N):
                                if Psu[i][k]:
                                    B[p][k][j] = B[p][k][j] + v * Psu[i][k]
            for sv in SECTORS:
                Psv = P[sv]
                total = ZERO
                for p in range(N):
                    for k in range(N):
                        for j in range(N):
                            v = B[p][k][j]
                            if v:
                                for l in range(N):
                                    if Psv[j][l] and S[p][k][l]:
                                        total = total + v * Psv[j][l] * S[p][k][l]
                if total:
                    blocks[(so, su, sv)] = total
    return blocks


def pairing_at(blocks, scales):
    """Evaluate the induced pairing at concrete sector scales (Fractions)."""
    total = ZERO
    for (so, su, sv), value in blocks.items():
        factor = Fraction(scales[so]) / (Fraction(scales[su]) * Fraction(scales[sv]))
        total = total + value * Q5(factor)
    return total * Q5(Fraction(1, 2))


# --------------------------------------------------------------------- laurent
# A Laurent table maps monomial keys to Q5 coefficients.  Keys are strings
# over the multiplicative basis {beta/delta^2, 1/beta, gamma/delta^2, 1/gamma}.
MONO_KEYS = ["beta_over_delta_sq", "inv_beta", "gamma_over_delta_sq", "inv_gamma"]


def laurent_eval(table, beta, gamma, delta):
    b, g, d = Fraction(beta), Fraction(gamma), Fraction(delta)
    values = {
        "beta_over_delta_sq": b / (d * d),
        "inv_beta": 1 / b,
        "gamma_over_delta_sq": g / (d * d),
        "inv_gamma": 1 / g,
    }
    total = ZERO
    for key, coeff in table.items():
        total = total + coeff * Q5(values[key])
    return total


def laurent_sub(x, y):
    out = dict(x)
    for key, coeff in y.items():
        out[key] = out.get(key, ZERO) - coeff
        if not out[key]:
            del out[key]
    return out


def enc_table(table):
    return {key: enc(value) for key, value in sorted(table.items())}


# ------------------------------------------------------------------------ main
def main() -> dict:
    manifest, basis_raw, receipt, compact, selector = load_inputs()
    group = group_rows(receipt)
    orbital, projectors = spectral_projectors(group)
    C = structure_constants(basis_raw)
    vectors = [{(o, l, r): Fraction(num, den) for o, l, r, num, den in row["entries"]}
               for row in basis_raw["basis"]]

    # (1) completeness: the commutant dimension equals the sector count.
    commutant_dim = commutant_dimension(group)
    require(commutant_dim == 4, "commutant dimension is not four")
    # the four symmetric orthogonal projectors with nonzero traces are
    # independent commutant elements, hence a basis.
    for label in SECTORS:
        P = projectors[label]
        for g in group:
            image = [[P[g[i]][g[j]] for j in range(N)] for i in range(N)]
            require(image == P, f"{label} projector is not invariant under the port action")

    # (2) channel machinery and diagonality.
    transform, inverse = channel_transform(vectors, C, projectors)
    channel_tensors = [q5_tensor(C, [inverse[a][k] for a in range(M)]) for k in range(M)]
    self_blocks = [sector_blocks(projectors, channel_tensors[k], channel_tensors[k])
                   for k in range(M)]
    cross_zero_pairs = 0
    for k1 in range(M):
        for k2 in range(k1 + 1, M):
            blocks = sector_blocks(projectors, channel_tensors[k1], channel_tensors[k2])
            require(not blocks, f"channels {CHANNELS[k1]} and {CHANNELS[k2]} are not orthogonal")
            cross_zero_pairs += 1
    require(cross_zero_pairs == M * (M - 1) // 2, "cross-pair count mismatch")

    # (3) face bracket: rebuild from the pinned manifest faces and locate in
    # channel coordinates; check against the pinned selector certificate.
    ports = manifest["carrier"]["ports"]
    require(ports == [f"p{i:02d}" for i in range(N)], "unexpected source port order")
    index = {port: i for i, port in enumerate(ports)}
    faces = [tuple(index[p] for p in face) for face in manifest["carrier"]["oriented_faces"]]
    require(len(faces) == 20 and len(set(faces)) == 20, "expected twenty distinct oriented faces")
    face_tensor = [[[ZERO] * N for _ in range(N)] for _ in range(N)]
    for a, b, c in faces:
        for out, left, right in ((c, a, b), (a, b, c), (b, c, a)):
            face_tensor[out][left][right] = face_tensor[out][left][right] + ONE
            face_tensor[out][right][left] = face_tensor[out][right][left] - ONE
    # channel coordinates of the face bracket via the transform on x-coords.
    # First express the face tensor in Reynolds x-coordinates through the
    # pinned selector certificate value B_face = 60 R13, verified directly.
    sixty_r13 = q5_tensor(C, [Q5(60) if a == 13 else ZERO for a in range(M)])
    require(face_tensor == sixty_r13, "face bracket is not 60 R13")
    x_face = [ZERO] * M
    x_face[13] = Q5(60)
    face_channels = {CHANNELS[k]: sum((transform[k][a] * x_face[a] for a in range(M)), ZERO)
                     for k in range(M)}
    pinned_channels = {name: Q5(Fraction(v[0], v[1]), Fraction(v[2], v[3]))
                       for name, v in selector["channel_coordinates"].items()}
    require(face_channels == pinned_channels,
            "face channel coordinates differ from the pinned selector certificate")
    for name in CHANNELS:
        expected = ROOT5 * Q5(Fraction(1, 20))
        value = face_channels[name]
        if name in ("t_pp_to_p", "t_ff_to_p", "t_pf_to_f"):
            require(value == expected, f"face channel {name} differs")
        elif name in ("t_mm_to_m", "t_ff_to_m", "t_mf_to_f"):
            require(value == -expected, f"face channel {name} differs")
        else:
            require(not value, f"face channel {name} should vanish")

    # (4) per-channel monomial structure for the six face channels.
    # Each block key (so,su,sv) must correspond to one of the four monomials.
    def monomial_of(block_key):
        so, su, sv = block_key
        pair = tuple(sorted((su, sv)))
        if so == "three_plus" and pair == ("three_plus", "three_plus"):
            return "inv_beta"
        if so == "three_minus" and pair == ("three_minus", "three_minus"):
            return "inv_gamma"
        if so == "three_plus" and pair == ("five", "five"):
            return "beta_over_delta_sq"
        if so == "three_minus" and pair == ("five", "five"):
            return "gamma_over_delta_sq"
        if so == "five" and pair == ("five", "three_plus"):
            return "inv_beta"
        if so == "five" and pair == ("five", "three_minus"):
            return "inv_gamma"
        return None

    channel_norm_tables = {}
    for name in FACE_CHANNELS:
        k = CH[name]
        table = {}
        for key, value in self_blocks[k].items():
            mono = monomial_of(key)
            require(mono is not None, f"unexpected block {key} in channel {name}")
            table[mono] = table.get(mono, ZERO) + value * Q5(Fraction(1, 2))
        require(len(table) == 1, f"channel {name} spans more than one monomial")
        channel_norm_tables[name] = table

    expected_norms = {
        "t_pp_to_p": ("inv_beta", Q5(600, -120)),
        "t_mm_to_m": ("inv_gamma", Q5(600, 120)),
        "t_ff_to_p": ("beta_over_delta_sq", Q5(600, -120)),
        "t_ff_to_m": ("gamma_over_delta_sq", Q5(600, 120)),
        "t_pf_to_f": ("inv_beta", Q5(1200, -240)),
        "t_mf_to_f": ("inv_gamma", Q5(1200, 240)),
    }
    for name, (mono, value) in expected_norms.items():
        table = channel_norm_tables[name]
        require(list(table.keys()) == [mono] and table[mono] == value,
                f"channel {name} norm table differs from the expected monomial")

    # (5) closed-form distances.  Channel orthogonality decouples the least
    # squares per channel.  Free family parameters match their target channels
    # exactly; the coupled F/G parameter minimizes a one-dimensional quadratic
    # whose two weights share one monomial, so the minimizer is metric
    # independent.
    face_sq = {name: face_channels[name] * face_channels[name] for name in FACE_CHANNELS}

    def contribution(name):
        mono, value = expected_norms[name]
        return mono, face_sq[name] * value

    def add(table, mono, value):
        table[mono] = table.get(mono, ZERO) + value
        if not table[mono]:
            del table[mono]

    # face squared norm table
    norm_face = {}
    for name in FACE_CHANNELS:
        mono, value = contribution(name)
        add(norm_face, mono, value)

    # P covers t_pp_to_p and t_mm_to_m exactly.
    dP = {}
    for name in ("t_ff_to_p", "t_ff_to_m", "t_pf_to_f", "t_mf_to_f"):
        mono, value = contribution(name)
        add(dP, mono, value)

    # G covers t_pp_to_p (e) and t_ff_to_m (b) exactly; the coupled parameter a
    # fits t_mm_to_m (residual a + sqrt5/20) and t_mf_to_f (residual
    # sqrt5 (a + 1/20)); both carry monomial 1/gamma.
    #   K_G = min_a  n_mmm (a + sqrt5/20)^2 + 5 n_mff (a + 1/20)^2.
    n_mmm = expected_norms["t_mm_to_m"][1]
    n_mff = expected_norms["t_mf_to_f"][1]
    wA = n_mmm
    wB = Q5(5) * n_mff
    pA = -ROOT5 * Q5(Fraction(1, 20))
    pB = -Q5(Fraction(1, 20))
    a_G = (wA * pA + wB * pB) / (wA + wB)
    K_G = wA * (a_G - pA) * (a_G - pA) + wB * (a_G - pB) * (a_G - pB)
    dG = {}
    for name in ("t_ff_to_p", "t_pf_to_f"):
        mono, value = contribution(name)
        add(dG, mono, value)
    add(dG, "inv_gamma", K_G)
    b_G = face_channels["t_ff_to_m"]
    require(sign(a_G) < 0 and sign(b_G) < 0 and sign(a_G * b_G) > 0,
            "G minimizer does not satisfy the open compact-stratum sign a*b > 0")

    # F covers t_mm_to_m (e) and t_ff_to_p (b) exactly; the coupled parameter a
    # fits t_pp_to_p (residual a - sqrt5/20) and t_pf_to_f (family value
    # -sqrt5 a against target sqrt5/20, residual sqrt5 (a + 1/20)); both carry
    # monomial 1/beta.
    n_ppp = expected_norms["t_pp_to_p"][1]
    n_pff = expected_norms["t_pf_to_f"][1]
    wA_F = n_ppp
    wB_F = Q5(5) * n_pff
    pA_F = ROOT5 * Q5(Fraction(1, 20))
    pB_F = -Q5(Fraction(1, 20))
    a_F = (wA_F * pA_F + wB_F * pB_F) / (wA_F + wB_F)
    K_F = wA_F * (a_F - pA_F) * (a_F - pA_F) + wB_F * (a_F - pB_F) * (a_F - pB_F)
    dF = {}
    for name in ("t_ff_to_m", "t_mf_to_f"):
        mono, value = contribution(name)
        add(dF, mono, value)
    add(dF, "inv_beta", K_F)
    b_F = face_channels["t_ff_to_p"]
    require(sign(a_F) < 0 and sign(b_F) > 0 and sign(a_F * b_F) < 0,
            "F minimizer does not satisfy the open compact-stratum sign a*b < 0")

    require(K_G == Q5(Fraction(60, 11), Fraction(-12, 11)), "K_G differs from (60-12sqrt5)/11")
    require(K_F == Q5(Fraction(60, 11), Fraction(12, 11)), "K_F differs from (60+12sqrt5)/11")
    require(K_F == K_G.conj(), "K_F is not the Galois conjugate of K_G")

    # (6) reference-point calibration against the pinned selector L2 values.
    ref = {"fixed": 1, "three_plus": 1, "three_minus": 1, "five": 1}
    dP_ref = laurent_eval(dP, 1, 1, 1)
    dG_ref = laurent_eval(dG, 1, 1, 1)
    dF_ref = laurent_eval(dF, 1, 1, 1)
    require(dP_ref == Q5(45), "reference d_P^2 differs from the pinned 45")
    require(dG_ref == Q5(Fraction(615, 22), Fraction(-123, 22)),
            "reference d_G^2 differs from the pinned (615-123sqrt5)/22")
    require(dF_ref == Q5(Fraction(615, 22), Fraction(123, 22)),
            "reference d_F^2 differs from the pinned (615+123sqrt5)/22")

    # (7) independent brute-force check at random rational scale points: solve
    # the exact normal equations over the full family spans in the ambient
    # tensor space with the induced metric, and compare with the Laurent
    # tables.
    def family_span_tensors(family):
        if family == "P":
            gens = [{"t_pp_to_p": ONE}, {"t_mm_to_m": ONE}]
        elif family == "F":
            gens = [{"t_pp_to_p": ONE, "t_pf_to_f": -ROOT5},
                    {"t_ff_to_p": ONE}, {"t_mm_to_m": ONE}]
        else:
            gens = [{"t_mm_to_m": ONE, "t_mf_to_f": ROOT5},
                    {"t_ff_to_m": ONE}, {"t_pp_to_p": ONE}]
        tensors = []
        for gen in gens:
            coeffs = [sum((inverse[a][CH[name]] * value for name, value in gen.items()), ZERO)
                      for a in range(M)]
            tensors.append(q5_tensor(C, coeffs))
        return tensors

    def solve_normal(gram, rhs):
        n = len(rhs)
        aug = [list(row) + [rhs[i]] for i, row in enumerate(gram)]
        for col in range(n):
            pivot = next((r for r in range(col, n) if aug[r][col]), None)
            require(pivot is not None, "singular family Gram matrix")
            aug[col], aug[pivot] = aug[pivot], aug[col]
            inv = aug[col][col].reciprocal()
            aug[col] = [inv * x for x in aug[col]]
            for r in range(n):
                if r != col and aug[r][col]:
                    f = aug[r][col]
                    aug[r] = [x - f * y for x, y in zip(aug[r], aug[col])]
        return [aug[i][-1] for i in range(n)]

    sample_points = [(1, 1, 1), (2, 3, 5), (7, 2, 3), (1, 9, 4),
                     (Fraction(1, 3), Fraction(5, 2), 1), (11, 1, 6)]
    family_tensors = {family: family_span_tensors(family) for family in ("P", "F", "G")}
    face_blocks_with = {
        family: [sector_blocks(projectors, face_tensor, t) for t in family_tensors[family]]
        for family in family_tensors
    }
    span_blocks = {
        family: [[sector_blocks(projectors, ti, tj) for tj in family_tensors[family]]
                 for ti in family_tensors[family]]
        for family in family_tensors
    }
    face_self_blocks = sector_blocks(projectors, face_tensor, face_tensor)
    samples = []
    for beta, gamma, delta in sample_points:
        scales = {"fixed": 1, "three_plus": beta, "three_minus": gamma, "five": delta}
        norm_face_val = pairing_at(face_self_blocks, scales)
        point = {"scales": {"beta": [Fraction(beta).numerator, Fraction(beta).denominator],
                            "gamma": [Fraction(gamma).numerator, Fraction(gamma).denominator],
                            "delta": [Fraction(delta).numerator, Fraction(delta).denominator]},
                 "families": {}}
        for family, table in (("P", dP), ("F", dF), ("G", dG)):
            gram = [[pairing_at(span_blocks[family][i][j], scales)
                     for j in range(len(family_tensors[family]))]
                    for i in range(len(family_tensors[family]))]
            rhs = [pairing_at(face_blocks_with[family][i], scales)
                   for i in range(len(family_tensors[family]))]
            coeffs = solve_normal(gram, rhs)
            projected = sum((coeffs[i] * rhs[i] for i in range(len(rhs))), ZERO)
            brute = norm_face_val - projected
            closed = laurent_eval(table, beta, gamma, delta)
            require(brute == closed,
                    f"brute-force distance to {family} differs at {(beta, gamma, delta)}")
            point["families"][family] = enc(closed)
        samples.append(point)

    # (8) phase-diagram inequalities.
    # P-exclusion: coefficient-positive gap tables.
    gap_PG = laurent_sub(dP, dG)
    gap_PF = laurent_sub(dP, dF)
    require(set(gap_PG) == {"gamma_over_delta_sq", "inv_gamma"}, "unexpected P-G gap support")
    require(set(gap_PF) == {"beta_over_delta_sq", "inv_beta"}, "unexpected P-F gap support")
    for table in (gap_PG, gap_PF):
        for coeff in table.values():
            require(sign(coeff) > 0, "P-exclusion gap has a nonpositive coefficient")
    require(gap_PG["gamma_over_delta_sq"] == Q5(Fraction(15, 2), Fraction(3, 2)), "P-G quadratic term")
    require(gap_PG["inv_gamma"] == Q5(Fraction(105, 11), Fraction(45, 11)), "P-G reciprocal term")
    require(gap_PF["beta_over_delta_sq"] == Q5(Fraction(15, 2), Fraction(-3, 2)), "P-F quadratic term")
    require(gap_PF["inv_beta"] == Q5(Fraction(105, 11), Fraction(-45, 11)), "P-F reciprocal term")

    # balanced slice: dF - dG at gamma = beta has positive coefficients.
    # dF - dG = aF1 gamma/d^2 + C_plus/gamma - aG1 beta/d^2 - C_minus/beta.
    gap_FG = laurent_sub(dF, dG)
    aF1 = gap_FG["gamma_over_delta_sq"]
    aG1 = -gap_FG["beta_over_delta_sq"]
    C_plus = gap_FG["inv_gamma"]
    C_minus = -gap_FG["inv_beta"]
    require(aF1 == Q5(Fraction(15, 2), Fraction(3, 2)), "aF1 differs")
    require(aG1 == Q5(Fraction(15, 2), Fraction(-3, 2)), "aG1 differs")
    require(C_plus == Q5(Fraction(105, 11), Fraction(45, 11)), "C_plus differs")
    require(C_minus == Q5(Fraction(105, 11), Fraction(-45, 11)), "C_minus differs")
    require(aF1 == aG1.conj() and C_plus == C_minus.conj(),
            "phase-boundary constants are not Galois conjugate pairs")
    balanced_linear = aF1 - aG1
    balanced_reciprocal = C_plus - C_minus
    require(balanced_linear == Q5(0, 3), "balanced slice linear coefficient is not 3 sqrt5")
    require(balanced_reciprocal == Q5(0, Fraction(90, 11)),
            "balanced slice reciprocal coefficient is not (90/11) sqrt5")

    # phase box: for beta/delta in [1/50, 6], k(u) = aG1 u + C_minus/u < 28,
    # certified by endpoint signs of the upward parabola aG1 u^2 - 28 u + C_minus,
    # and h(v) = aF1 v + C_plus/v > 28 for all v > 0, certified by the
    # discriminant sign 28^2 - 4 aF1 C_plus < 0.
    K_box = Q5(28)
    end_low = aG1 * Q5(Fraction(1, 2500)) - K_box * Q5(Fraction(1, 50)) + C_minus
    end_high = aG1 * Q5(36) - K_box * Q5(6) + C_minus
    discriminant = K_box * K_box - Q5(4) * aF1 * C_plus
    require(sign(end_low) < 0, "phase box fails at beta/delta = 1/50")
    require(sign(end_high) < 0, "phase box fails at beta/delta = 6")
    require(sign(discriminant) < 0, "phase box discriminant is not negative")

    # F-region witness at (8,1,1) and G witness at reference.
    dF_w = laurent_eval(dF, 8, 1, 1)
    dG_w = laurent_eval(dG, 8, 1, 1)
    require(sign(dG_w - dF_w) > 0, "F witness fails at (8,1,1)")
    require(sign(dF_ref - dG_ref) > 0, "G witness fails at the reference point")

    # (9) Galois duality of the closed forms: swapping beta <-> gamma and
    # conjugating coefficients maps dG to dF.
    swap = {"beta_over_delta_sq": "gamma_over_delta_sq",
            "gamma_over_delta_sq": "beta_over_delta_sq",
            "inv_beta": "inv_gamma", "inv_gamma": "inv_beta"}
    dG_swapped_conj = {swap[key]: value.conj() for key, value in dG.items()}
    require(dG_swapped_conj == dF, "Galois duality fails between dG and dF")

    # (10) non-induced control: scale the t_pf_to_f channel weight by 6 at the
    # balanced reference point.  Carrier-induced metrics force the normalized
    # t_pp_to_p and t_pf_to_f weights to share the monomial 1/beta, so this
    # diagonal reweighting is invariant but not carrier-induced; it reverses
    # the balanced-point selection.
    s = Q5(6)
    # G's t_pf_to_f cost scales by six; F refits its coupled parameter with
    # the scaled weight.
    dG_control = dG_ref + (s - ONE) * face_sq["t_pf_to_f"] * n_pff
    wB_c = s * Q5(5) * n_pff
    a_Fc = (wA_F * pA_F + wB_c * pB_F) / (wA_F + wB_c)
    K_Fc = wA_F * (a_Fc - pA_F) * (a_Fc - pA_F) + wB_c * (a_Fc - pB_F) * (a_Fc - pB_F)
    dF_control = laurent_eval({key: value for key, value in dF.items() if key != "inv_beta"}, 1, 1, 1) + K_Fc
    require(sign(dG_control - dF_control) > 0,
            "non-induced control does not reverse the balanced-point selection")

    certificate = {
        "schema": "b14-invariant-metric-phase/1",
        "issue": "https://github.com/FloatingPragma/observer-patch-holography/issues/705",
        "upstream": {
            "manifest_sha256": file_hash(MANIFEST),
            "reynolds_basis_sha256": file_hash(BASIS),
            "stage1_receipt_sha256": file_hash(RECEIPT),
            "compact_locus_sha256": file_hash(COMPACT),
            "selector_sha256": file_hash(SELECTOR),
        },
        "completeness": {
            "commutant_dimension_over_Q": commutant_dim,
            "sector_dimensions": SECTOR_DIMS,
            "projectors_symmetric_orthogonal_idempotent_invariant": True,
            "statement": (
                "the ordered-pair orbit count of the proper port action is four, so the "
                "commutant is four dimensional and is spanned by the four symmetric "
                "invariant spectral projectors; the commutation system is rational, so "
                "its solution dimension is the same over the reals, and every invariant "
                "inner product on the carrier is a positive sector-scale metric"
            ),
        },
        "diagonality": {
            "cross_channel_pairs_checked": cross_zero_pairs,
            "all_cross_pairings_vanish": True,
        },
        "face_channel_coordinates": {name: enc(face_channels[name]) for name in CHANNELS},
        "channel_norm_monomials": {
            name: {"monomial": mono, "norm": enc(value)}
            for name, (mono, value) in expected_norms.items()
        },
        "closed_forms": {
            "monomial_basis": MONO_KEYS,
            "norm_face": enc_table(norm_face),
            "d_P_squared": enc_table(dP),
            "d_F_squared": enc_table(dF),
            "d_G_squared": enc_table(dG),
            "coupled_minimizers": {
                "G": {"a": enc(a_G), "b": enc(b_G), "K": enc(K_G),
                      "open_stratum_sign": "a*b > 0 holds strictly"},
                "F": {"a": enc(a_F), "b": enc(b_F), "K": enc(K_F),
                      "open_stratum_sign": "a*b < 0 holds strictly"},
            },
            "alpha_absent": True,
        },
        "reference_point": {
            "d_P_squared": enc(dP_ref),
            "d_G_squared": enc(dG_ref),
            "d_F_squared": enc(dF_ref),
            "matches_pinned_selector_l2": True,
        },
        "samples": samples,
        "phase_diagram": {
            "P_exclusion": {
                "gap_P_minus_G": enc_table(gap_PG),
                "gap_P_minus_F": enc_table(gap_PF),
                "statement": "both gaps are coefficient-positive Laurent polynomials, so P is strictly excluded for every invariant carrier metric",
            },
            "balanced_slice": {
                "gap_F_minus_G_at_balanced": {
                    "beta_over_delta_sq": enc(balanced_linear),
                    "inv_beta": enc(balanced_reciprocal),
                },
                "statement": "at gamma = beta the F-G gap is 3 sqrt5 beta/delta^2 + (90 sqrt5/11)/beta > 0, so every sector-balanced invariant metric selects G",
            },
            "phase_box": {
                "K": enc(K_box),
                "interval": [[1, 50], [6, 1]],
                "parabola_endpoint_low": enc(end_low),
                "parabola_endpoint_high": enc(end_high),
                "h_discriminant": enc(discriminant),
                "statement": "for beta/delta in [1/50, 6] the k-branch stays below 28 while the h-branch exceeds 28 for every gamma, so G is selected on the whole box for all gamma and delta",
            },
            "F_region_witness": {
                "point": {"beta": [8, 1], "gamma": [1, 1], "delta": [1, 1]},
                "d_F_squared": enc(dF_w),
                "d_G_squared": enc(dG_w),
                "statement": "at beta = 8 delta the F family is strictly nearer, so the sector-asymmetry threshold is real",
            },
            "tie_surface": "h(gamma/delta) = k(beta/delta) with h(v) = aF1 v + C_plus/v and k(u) = aG1 u + C_minus/u",
            "constants": {
                "aF1": enc(aF1), "aG1": enc(aG1),
                "C_plus": enc(C_plus), "C_minus": enc(C_minus),
                "galois_conjugate_pairs": True,
            },
        },
        "galois_duality": {
            "statement": "conjugating sqrt5 and swapping beta <-> gamma maps d_G to d_F exactly",
            "verified": True,
        },
        "non_induced_control": {
            "channel_scaled": "t_pf_to_f",
            "scale": [6, 1],
            "d_G_control": enc(dG_control),
            "d_F_control": enc(dF_control),
            "statement": (
                "scaling the t_pf_to_f channel weight by six at the balanced reference "
                "point reverses the selection to F; carrier-induced metrics force the "
                "t_pp_to_p and t_pf_to_f weights to share the monomial 1/beta, so this "
                "control metric is invariant on the bracket space but not carrier-induced, "
                "and carrier-inducedness is load-bearing"
            ),
        },
        "claim_boundary": (
            "This certificate proves the exact phase diagram of the nearest-classified-"
            "compact-family rule over the complete four-parameter cone of invariant "
            "carrier inner products: P is excluded everywhere, every sector-balanced "
            "metric and every metric with beta/delta in [1/50, 6] selects G, and F "
            "requires an explicit sector-asymmetric (chiral) metric past the exact tie "
            "surface.  The comparison ranges over the compact locus classified by the "
            "pinned certificate, conditional on its three named textbook lemmas.  The "
            "nearest-point repair rule and the carrier-induced metric class are declared "
            "discriminator choices; no norm, bracket, current, or physical gauge "
            "structure is source-selected, and no laboratory identification is made."
        ),
        "implementation_pins": {
            "producer_sha256": file_hash(Path(__file__).resolve()),
            "verifier_sha256": file_hash(VERIFIER) if VERIFIER.exists() else "absent",
        },
    }
    certificate["certificate_sha256"] = object_hash(certificate)
    return certificate


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true",
                        help="recompute and compare against the committed certificate")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    cli = parser.parse_args()
    result = main()
    payload = json.dumps(result, sort_keys=True, indent=1) + "\n"
    if cli.check:
        committed = OUTPUT.read_text(encoding="utf-8")
        if committed != payload:
            raise CertificateError("recomputed certificate differs from the committed bytes")
        print("committed certificate is current")
    else:
        cli.output.write_text(payload, encoding="utf-8")
        print("certificate written:", cli.output)
        print("certificate_sha256:", result["certificate_sha256"])
        print("phase box:", result["phase_diagram"]["phase_box"]["statement"])
