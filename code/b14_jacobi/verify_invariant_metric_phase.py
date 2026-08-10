#!/usr/bin/env python3
"""Independent exact replay of the B14 invariant-metric phase certificate.

No code is imported from the producer.  This verifier separately rebuilds the
sector projectors, computes the commutant dimension by solving the actual
linear commutation system for a generating pair of the port action, rebuilds
the face and family tensors, evaluates the induced metric directly through the
raw carrier metric matrices g and g^{-1} at each certificate sample point, and
replays every phase-diagram sign fact from the serialized Laurent tables.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT = HERE / "invariant_metric_phase.certificate.json"
PRODUCER = HERE / "invariant_metric_phase.py"
MANIFEST = REPO / "code/a5_closure/manifests/echosahedral_federation_reference.json"
BASIS = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json"
RECEIPT = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json"
COMPACT = HERE / "b14_compact_locus.certificate.json"
SELECTOR = HERE / "oriented_face_bracket_selector.certificate.json"

N = 12
M = 14
SECTORS = ["fixed", "three_plus", "three_minus", "five"]
SECTOR_DIMS = {"fixed": 1, "three_plus": 3, "three_minus": 3, "five": 5}
MONO_KEYS = ["beta_over_delta_sq", "inv_beta", "gamma_over_delta_sq", "inv_gamma"]
FACE_CHANNELS = ["t_pp_to_p", "t_mm_to_m", "t_ff_to_p", "t_ff_to_m",
                 "t_pf_to_f", "t_mf_to_f"]


class VerificationError(RuntimeError):
    pass


def check(value: bool, message: str) -> None:
    if not value:
        raise VerificationError(message)


class F5:
    """Exact r + s sqrt(5), independent implementation."""

    __slots__ = ("r", "s")

    def __init__(self, r=0, s=0):
        self.r = r if isinstance(r, Fraction) else Fraction(r)
        self.s = s if isinstance(s, Fraction) else Fraction(s)

    def __add__(self, other):
        other = other if isinstance(other, F5) else F5(other)
        return F5(self.r + other.r, self.s + other.s)

    __radd__ = __add__

    def __neg__(self):
        return F5(-self.r, -self.s)

    def __sub__(self, other):
        other = other if isinstance(other, F5) else F5(other)
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        other = other if isinstance(other, F5) else F5(other)
        return F5(self.r * other.r + 5 * self.s * other.s,
                  self.r * other.s + self.s * other.r)

    __rmul__ = __mul__

    def inverse(self):
        norm = self.r * self.r - 5 * self.s * self.s
        check(norm != 0, "division by zero")
        return F5(self.r / norm, -self.s / norm)

    def __truediv__(self, other):
        other = other if isinstance(other, F5) else F5(other)
        return self * other.inverse()

    def __eq__(self, other):
        other = other if isinstance(other, F5) else F5(other)
        return self.r == other.r and self.s == other.s

    def __bool__(self):
        return bool(self.r or self.s)

    def conj(self):
        return F5(self.r, -self.s)


Z = F5()
ONE = F5(1)
RT5 = F5(0, 1)


def sgn(x: F5) -> int:
    if not x:
        return 0
    if x.r == 0:
        return 1 if x.s > 0 else -1
    if x.s == 0:
        return 1 if x.r > 0 else -1
    if (x.r > 0) == (x.s > 0):
        return 1 if x.r > 0 else -1
    comparison = x.r * x.r - 5 * x.s * x.s
    check(comparison != 0, "unexpected zero in sign test")
    if x.r > 0:
        return 1 if comparison > 0 else -1
    return 1 if comparison < 0 else -1


def dec(v) -> F5:
    return F5(Fraction(v[0], v[1]), Fraction(v[2], v[3]))


def file_hash(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_bytes(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def object_hash(value) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def matmul(X, Y):
    n, m, p = len(X), len(Y), len(Y[0])
    out = [[Z] * p for _ in range(n)]
    for i in range(n):
        for k in range(m):
            if X[i][k]:
                xik = X[i][k]
                for j in range(p):
                    if Y[k][j]:
                        out[i][j] = out[i][j] + xik * Y[k][j]
    return out


def identity(n):
    return [[ONE if i == j else Z for j in range(n)] for i in range(n)]


def build_projectors(group):
    pairs = {(i, j) for i in range(N) for j in range(N)}
    orbits = []
    while pairs:
        seed = min(pairs)
        orbit = {(g[seed[0]], g[seed[1]]) for g in group}
        pairs -= orbit
        orbits.append(orbit)
    check(sorted(len(o) for o in orbits) == [12, 12, 60, 60], "unexpected pair orbits")
    orbital = sorted((o for o in orbits if len(o) == 60), key=lambda o: tuple(sorted(o)))[0]
    A = [[Z] * N for _ in range(N)]
    for i, j in orbital:
        A[i][j] = ONE
    eigen = {"fixed": F5(5), "three_plus": RT5, "three_minus": -RT5, "five": F5(-1)}
    projectors = {}
    for label, ev in eigen.items():
        Mx = identity(N)
        for other, ev2 in eigen.items():
            if other == label:
                continue
            factor = [[A[i][j] - (ev2 if i == j else Z) for j in range(N)] for i in range(N)]
            scale = (ev - ev2).inverse()
            Mx = [[scale * v for v in row] for row in matmul(Mx, factor)]
        projectors[label] = Mx
    total = [[sum((projectors[l][i][j] for l in SECTORS), Z) for j in range(N)] for i in range(N)]
    check(total == identity(N), "projector partition of unity fails")
    for label in SECTORS:
        P = projectors[label]
        check(matmul(P, P) == P, "projector idempotence fails")
        check([list(r) for r in zip(*[[v for v in row] for row in P])] == P, "projector symmetry fails")
        tr = sum((P[i][i] for i in range(N)), Z)
        check(tr == F5(SECTOR_DIMS[label]), "projector trace fails")
    return orbits, projectors


def group_generators(group):
    """A generating subset of the port action, certified by exact closure."""
    generators = []
    generated = {tuple(range(N))}
    for candidate in group:
        if candidate in generated:
            continue
        generators.append(candidate)
        closure = set(generated)
        closure.add(candidate)
        queue = list(closure)
        while queue:
            item = queue.pop()
            for other in list(closure):
                for prod in (tuple(item[other[i]] for i in range(N)),
                             tuple(other[item[i]] for i in range(N))):
                    if prod not in closure:
                        closure.add(prod)
                        queue.append(prod)
        generated = closure
        if len(generated) == 60:
            break
    check(len(generated) == 60, "generators do not generate the port action")
    return generators


def commutant_nullity(group):
    """Q-dimension of the commutant, solved as the conjugation-invariance
    linear system X[g(i)][g(j)] = X[i][j] for a certified generating subset.
    For permutation matrices this is exactly X rho(g) = rho(g) X, and
    invariance under generators extends to the generated group because
    position conjugation is a group action."""
    generators = group_generators(group)
    rows = []
    for g in generators:
        for i in range(N):
            for j in range(N):
                row = [Fraction(0)] * (N * N)
                row[g[i] * N + g[j]] += 1
                row[i * N + j] -= 1
                if any(row):
                    rows.append(row)
    # rank over Q
    rank = 0
    cols = N * N
    pivot_col = 0
    work = [r for r in rows if any(r)]
    while work and pivot_col < cols:
        pivot = next((k for k, r in enumerate(work) if r[pivot_col]), None)
        if pivot is None:
            pivot_col += 1
            continue
        work[0], work[pivot] = work[pivot], work[0]
        head = work.pop(0)
        inv = 1 / head[pivot_col]
        head = [x * inv for x in head]
        work = [[x - r[pivot_col] * h for x, h in zip(r, head)] for r in work]
        work = [r for r in work if any(r)]
        rank += 1
        pivot_col += 1
    return N * N - rank


def sparse_entries(T):
    out = []
    for o in range(N):
        for i in range(N):
            for j in range(N):
                if T[o][i][j]:
                    out.append((o, i, j, T[o][i][j]))
    return out


def metric_image(T_entries, Gm, Gi):
    """T3[p][k][l] = sum_{o,i,j} T^o_ij Gm[o][p] Gi[i][k] Gi[j][l], staged."""
    T1 = {}
    for o, i, j, tv in T_entries:
        for l in range(N):
            g = Gi[j][l]
            if g:
                key = (o, i, l)
                T1[key] = T1.get(key, Z) + tv * g
    T2 = {}
    for (o, i, l), tv in T1.items():
        if not tv:
            continue
        for k in range(N):
            g = Gi[i][k]
            if g:
                key = (o, k, l)
                T2[key] = T2.get(key, Z) + tv * g
    T3 = {}
    for (o, k, l), tv in T2.items():
        if not tv:
            continue
        for p in range(N):
            g = Gm[o][p]
            if g:
                key = (p, k, l)
                T3[key] = T3.get(key, Z) + tv * g
    return T3


def induced_pairing_image(T3, S_entries):
    """(1/2) sum T3[p][k][l] S^p_kl with T3 precomputed by metric_image."""
    total = Z
    for p, k, l, sv in S_entries:
        value = T3.get((p, k, l))
        if value:
            total = total + value * sv
    return total * F5(Fraction(1, 2))


def laurent_eval(table, beta, gamma, delta):
    b, g, d = Fraction(beta), Fraction(gamma), Fraction(delta)
    values = {
        "beta_over_delta_sq": b / (d * d),
        "inv_beta": 1 / b,
        "gamma_over_delta_sq": g / (d * d),
        "inv_gamma": 1 / g,
    }
    total = Z
    for key, coeff in table.items():
        total = total + coeff * F5(values[key])
    return total


_FACTS_CACHE: dict = {}


def independent_facts() -> dict:
    """Input-derived exact facts, computed once per process."""
    if _FACTS_CACHE:
        return _FACTS_CACHE

    manifest = json.loads(MANIFEST.read_text())
    basis_raw = json.loads(BASIS.read_text())
    receipt = json.loads(RECEIPT.read_text())
    group = [tuple(row) for row in receipt["proper_port_action"]["permutation_rows"]]
    check(len(group) == 60 and len(set(group)) == 60, "group rows")

    orbits, projectors = build_projectors(group)

    # completeness: independent commutant nullity via the linear system, plus
    # explicit membership of the four independent symmetric projectors, so the
    # projectors span the commutant and every invariant inner product is a
    # positive sector-scale metric.
    nullity = commutant_nullity(group)
    check(nullity == 4, "commutant nullity is not four")
    check(len(orbits) == 4, "pair-orbit count is not four")
    for label in SECTORS:
        P = projectors[label]
        for g in group_generators(group):
            conjugated = [[P[g[i]][g[j]] for j in range(N)] for i in range(N)]
            check(conjugated == P, f"projector {label} is not in the commutant")

    # Reynolds tensors.
    check([row["basis_id"] for row in basis_raw["basis"]] == [f"R{i:02d}" for i in range(M)],
          "basis order")
    C = [[[[Z] * N for _ in range(N)] for _ in range(N)] for _ in range(M)]
    for a, row in enumerate(basis_raw["basis"]):
        for out, left, right, num, den in row["entries"]:
            value = F5(Fraction(num, den))
            C[a][out][left][right] = value
            C[a][out][right][left] = -value

    # face tensor from the manifest and the pinned identity with 60 R13.
    ports = manifest["carrier"]["ports"]
    index = {port: i for i, port in enumerate(ports)}
    faces = [tuple(index[p] for p in face) for face in manifest["carrier"]["oriented_faces"]]
    check(len(faces) == 20, "face count")
    face = [[[Z] * N for _ in range(N)] for _ in range(N)]
    for a, b, c in faces:
        for out, left, right in ((c, a, b), (a, b, c), (b, c, a)):
            face[out][left][right] = face[out][left][right] + ONE
            face[out][right][left] = face[out][right][left] - ONE
    sixty_r13 = [[[F5(60) * C[13][o][i][j] for j in range(N)] for i in range(N)] for o in range(N)]
    check(face == sixty_r13, "face tensor differs from 60 R13")

    # channel transform rebuilt the same way classify.py pins it.
    vectors = [{(o, l, r): Fraction(num, den) for o, l, r, num, den in row["entries"]}
               for row in basis_raw["basis"]]
    import itertools

    derivations = []
    for a in range(M):
        Mx = [[F5(sum((Fraction(C[a][o][s][c].r) for s in range(N)), Fraction(0)))
               for c in range(N)] for o in range(N)]
        derivations.append(Mx)
    eigenforms = {}
    for label in SECTORS:
        P = projectors[label]
        form = []
        for Mx in derivations:
            tr = Z
            for i in range(N):
                for j in range(N):
                    if P[i][j] and Mx[j][i]:
                        tr = tr + P[i][j] * Mx[j][i]
            form.append(tr / F5(SECTOR_DIMS[label]))
        eigenforms[label] = form
    specs = [
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
    channel_names = ["d_plus", "d_minus", "d_five"] + [s[0] for s in specs]
    transform = [eigenforms["three_plus"], eigenforms["three_minus"], eigenforms["five"]]
    for _, left, right, output in specs:
        found = None
        for coordinate in itertools.product(range(N), repeat=3):
            out, li, ri = coordinate
            form = []
            for vector in vectors:
                value = Z
                for (so, p, q), coefficient in vector.items():
                    of = projectors[output][out][so]
                    if not of:
                        continue
                    inf = (projectors[left][p][li] * projectors[right][q][ri]
                           - projectors[left][q][li] * projectors[right][p][ri])
                    if inf:
                        value = value + of * inf * F5(coefficient)
                form.append(value)
            if any(form):
                found = form
                break
        check(found is not None, "empty channel")
        transform.append(found)
    aug = [list(row) + [ONE if i == j else Z for j in range(M)]
           for i, row in enumerate(transform)]
    for col in range(M):
        pivot = next((r for r in range(col, M) if aug[r][col]), None)
        check(pivot is not None, "transform singular")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        inv = aug[col][col].inverse()
        aug[col] = [inv * x for x in aug[col]]
        for r in range(M):
            if r != col and aug[r][col]:
                f = aug[r][col]
                aug[r] = [x - f * y for x, y in zip(aug[r], aug[col])]
    inverse = [row[M:] for row in aug]

    chan = {name: k for k, name in enumerate(channel_names)}

    def tensor_from_channels(channel_values):
        coeffs = [sum((inverse[a][chan[name]] * value for name, value in channel_values.items()), Z)
                  for a in range(M)]
        T = [[[Z] * N for _ in range(N)] for _ in range(N)]
        for a, cval in enumerate(coeffs):
            if not cval:
                continue
            for o in range(N):
                for i in range(N):
                    for j in range(N):
                        if C[a][o][i][j]:
                            T[o][i][j] = T[o][i][j] + cval * C[a][o][i][j]
        return T

    spans = {
        "P": [{"t_pp_to_p": ONE}, {"t_mm_to_m": ONE}],
        "F": [{"t_pp_to_p": ONE, "t_pf_to_f": -RT5}, {"t_ff_to_p": ONE}, {"t_mm_to_m": ONE}],
        "G": [{"t_mm_to_m": ONE, "t_mf_to_f": RT5}, {"t_ff_to_m": ONE}, {"t_pp_to_p": ONE}],
    }
    span_tensors = {family: [sparse_entries(tensor_from_channels(gen)) for gen in gens]
                    for family, gens in spans.items()}
    face_entries = sparse_entries(face)

    distance_cache: dict = {}

    def replay_distance(family: str, beta: Fraction, gamma: Fraction, delta: Fraction) -> F5:
        """Exact least-squares distance through the raw carrier metric."""
        key = (family, beta, gamma, delta)
        if key in distance_cache:
            return distance_cache[key]
        lam = {"fixed": Fraction(1), "three_plus": beta, "three_minus": gamma, "five": delta}
        Gm = [[sum((F5(lam[s]) * projectors[s][i][j] for s in SECTORS), Z)
               for j in range(N)] for i in range(N)]
        Gi = [[sum((F5(1 / lam[s]) * projectors[s][i][j] for s in SECTORS), Z)
               for j in range(N)] for i in range(N)]
        check(matmul(Gm, Gi) == identity(N), "metric inverse fails")
        face_image = metric_image(face_entries, Gm, Gi)
        norm_face = induced_pairing_image(face_image, face_entries)
        gens = span_tensors[family]
        k = len(gens)
        images = [metric_image(gen, Gm, Gi) for gen in gens]
        gram = [[induced_pairing_image(images[i], gens[j]) for j in range(k)] for i in range(k)]
        rhs = [induced_pairing_image(face_image, gens[i]) for i in range(k)]
        aug2 = [list(gram[i]) + [rhs[i]] for i in range(k)]
        for col in range(k):
            pivot = next((r for r in range(col, k) if aug2[r][col]), None)
            check(pivot is not None, "singular Gram")
            aug2[col], aug2[pivot] = aug2[pivot], aug2[col]
            inv = aug2[col][col].inverse()
            aug2[col] = [inv * x for x in aug2[col]]
            for r in range(k):
                if r != col and aug2[r][col]:
                    f = aug2[r][col]
                    aug2[r] = [x - f * y for x, y in zip(aug2[r], aug2[col])]
        coeffs = [aug2[i][-1] for i in range(k)]
        projected = sum((coeffs[i] * rhs[i] for i in range(k)), Z)
        distance = norm_face - projected
        distance_cache[key] = distance
        return distance

    def channel_norm_at(name: str, beta: Fraction, gamma: Fraction, delta: Fraction) -> F5:
        """Squared norm of the pinned channel basis tensor under the metric."""
        entries = sparse_entries(tensor_from_channels({name: ONE}))
        lam = {"fixed": Fraction(1), "three_plus": beta, "three_minus": gamma, "five": delta}
        Gm = [[sum((F5(lam[s]) * projectors[s][i][j] for s in SECTORS), Z)
               for j in range(N)] for i in range(N)]
        Gi = [[sum((F5(1 / lam[s]) * projectors[s][i][j] for s in SECTORS), Z)
               for j in range(N)] for i in range(N)]
        return induced_pairing_image(metric_image(entries, Gm, Gi), entries)

    _FACTS_CACHE.update({
        "nullity": nullity,
        "orbit_count": len(orbits),
        "replay_distance": replay_distance,
        "channel_norm_at": channel_norm_at,
    })
    return _FACTS_CACHE


def verify_certificate(certificate_path: Path) -> dict:
    cert = json.loads(Path(certificate_path).read_text())
    check(cert["schema"] == "b14-invariant-metric-phase/1", "unexpected schema")

    # hash pinning
    up = cert["upstream"]
    check(up["manifest_sha256"] == file_hash(MANIFEST), "manifest hash drift")
    check(up["reynolds_basis_sha256"] == file_hash(BASIS), "basis hash drift")
    check(up["stage1_receipt_sha256"] == file_hash(RECEIPT), "receipt hash drift")
    check(up["compact_locus_sha256"] == file_hash(COMPACT), "compact locus hash drift")
    check(up["selector_sha256"] == file_hash(SELECTOR), "selector hash drift")
    body = {key: value for key, value in cert.items() if key != "certificate_sha256"}
    check(object_hash(body) == cert["certificate_sha256"], "certificate hash mismatch")

    facts = independent_facts()
    check(cert["completeness"]["commutant_dimension_over_Q"] == facts["nullity"] == 4,
          "certificate completeness field")

    # face channel coordinates against the pinned selector certificate.
    sel = json.loads(SELECTOR.read_text())
    check(set(cert["face_channel_coordinates"]) == set(sel["channel_coordinates"]),
          "face channel key set differs from the pinned selector")
    for name, value in cert["face_channel_coordinates"].items():
        check(dec(value) == dec(sel["channel_coordinates"][name]),
              f"face channel {name} differs from the pinned selector")

    # closed-form tables.
    dP = {key: dec(value) for key, value in cert["closed_forms"]["d_P_squared"].items()}
    dF = {key: dec(value) for key, value in cert["closed_forms"]["d_F_squared"].items()}
    dG = {key: dec(value) for key, value in cert["closed_forms"]["d_G_squared"].items()}
    for table in (dP, dF, dG):
        check(set(table) <= set(MONO_KEYS), "unexpected monomial key")

    # sample-point replay: every certificate sample row gets both the
    # closed-form check and the full independent contraction replay.
    check(len(cert["samples"]) >= 4, "too few sample points")
    for point in cert["samples"]:
        beta = Fraction(*point["scales"]["beta"])
        gamma = Fraction(*point["scales"]["gamma"])
        delta = Fraction(*point["scales"]["delta"])
        for family, table in (("P", dP), ("F", dF), ("G", dG)):
            value = dec(point["families"][family])
            check(laurent_eval(table, beta, gamma, delta) == value,
                  f"sample row disagrees with the closed form for {family}")
            check(facts["replay_distance"](family, beta, gamma, delta) == value,
                  f"distance replay differs for {family} at {(beta, gamma, delta)}")

    # reference-point pinned values.
    check(laurent_eval(dP, 1, 1, 1) == F5(45), "reference d_P")
    check(laurent_eval(dG, 1, 1, 1) == F5(Fraction(615, 22), Fraction(-123, 22)), "reference d_G")
    check(laurent_eval(dF, 1, 1, 1) == F5(Fraction(615, 22), Fraction(123, 22)), "reference d_F")

    # phase-diagram sign replay.
    pd = cert["phase_diagram"]
    gap_PG = {key: dec(value) for key, value in pd["P_exclusion"]["gap_P_minus_G"].items()}
    gap_PF = {key: dec(value) for key, value in pd["P_exclusion"]["gap_P_minus_F"].items()}
    for key in MONO_KEYS:
        left = dP.get(key, Z) - dG.get(key, Z)
        check(left == gap_PG.get(key, Z), "P-G gap table mismatch")
        left = dP.get(key, Z) - dF.get(key, Z)
        check(left == gap_PF.get(key, Z), "P-F gap table mismatch")
    for table in (gap_PG, gap_PF):
        check(table and all(sgn(v) > 0 for v in table.values()), "P-exclusion positivity fails")

    aF1 = dec(pd["constants"]["aF1"])
    aG1 = dec(pd["constants"]["aG1"])
    C_plus = dec(pd["constants"]["C_plus"])
    C_minus = dec(pd["constants"]["C_minus"])
    check(aF1 == dF.get("gamma_over_delta_sq", Z) - dG.get("gamma_over_delta_sq", Z), "aF1")
    check(aG1 == dG.get("beta_over_delta_sq", Z) - dF.get("beta_over_delta_sq", Z), "aG1")
    check(C_plus == dF.get("inv_gamma", Z) - dG.get("inv_gamma", Z), "C_plus")
    check(C_minus == dG.get("inv_beta", Z) - dF.get("inv_beta", Z), "C_minus")
    check(aF1 == aG1.conj() and C_plus == C_minus.conj(), "Galois pairs fail")
    check(aF1 - aG1 == F5(0, 3), "balanced linear coefficient")
    check(C_plus - C_minus == F5(0, Fraction(90, 11)), "balanced reciprocal coefficient")
    check(sgn(aF1) > 0 and sgn(aG1) > 0 and sgn(C_plus) > 0 and sgn(C_minus) > 0,
          "phase constants must be positive")

    # phase box replay.
    box = pd["phase_box"]
    K = dec(box["K"])
    check(K == F5(28), "box constant")
    low = aG1 * F5(Fraction(1, 2500)) - K * F5(Fraction(1, 50)) + C_minus
    high = aG1 * F5(36) - K * F5(6) + C_minus
    disc = K * K - F5(4) * aF1 * C_plus
    check(low == dec(box["parabola_endpoint_low"]), "box low endpoint")
    check(high == dec(box["parabola_endpoint_high"]), "box high endpoint")
    check(disc == dec(box["h_discriminant"]), "box discriminant")
    check(sgn(low) < 0 and sgn(high) < 0 and sgn(disc) < 0, "box sign facts fail")

    # F witness and Galois duality.
    wit = pd["F_region_witness"]
    dFw = laurent_eval(dF, Fraction(*wit["point"]["beta"]),
                       Fraction(*wit["point"]["gamma"]), Fraction(*wit["point"]["delta"]))
    dGw = laurent_eval(dG, Fraction(*wit["point"]["beta"]),
                       Fraction(*wit["point"]["gamma"]), Fraction(*wit["point"]["delta"]))
    check(dFw == dec(wit["d_F_squared"]) and dGw == dec(wit["d_G_squared"]), "witness values")
    check(sgn(dGw - dFw) > 0, "F witness ordering fails")
    swap = {"beta_over_delta_sq": "gamma_over_delta_sq",
            "gamma_over_delta_sq": "beta_over_delta_sq",
            "inv_beta": "inv_gamma", "inv_gamma": "inv_beta"}
    check({swap[key]: value.conj() for key, value in dG.items()} == dF, "Galois duality fails")

    # channel norm monomials against the independent metric replay.
    mono_points = {
        "beta_over_delta_sq": None, "inv_beta": None,
        "gamma_over_delta_sq": None, "inv_gamma": None,
    }
    probe_points = [(Fraction(1), Fraction(1), Fraction(1)),
                    (Fraction(2), Fraction(1), Fraction(1)),
                    (Fraction(1), Fraction(3), Fraction(1)),
                    (Fraction(1), Fraction(1), Fraction(5))]
    for name, entry in cert["channel_norm_monomials"].items():
        mono = entry["monomial"]
        check(mono in mono_points, "unknown channel monomial")
        norm = dec(entry["norm"])
        for beta, gamma, delta in probe_points:
            values = {
                "beta_over_delta_sq": beta / (delta * delta),
                "inv_beta": 1 / beta,
                "gamma_over_delta_sq": gamma / (delta * delta),
                "inv_gamma": 1 / gamma,
            }
            expected = norm * F5(values[mono])
            check(facts["channel_norm_at"](name, beta, gamma, delta) == expected,
                  f"channel norm replay fails for {name} at {(beta, gamma, delta)}")

    # non-induced control replay.
    ctrl = cert["non_induced_control"]
    check(ctrl["channel_scaled"] == "t_pf_to_f" and ctrl["scale"] == [6, 1], "control shape")
    dGc = dec(ctrl["d_G_control"])
    dFc = dec(ctrl["d_F_control"])
    face_sq_pff = dec(cert["face_channel_coordinates"]["t_pf_to_f"])
    face_sq_pff = face_sq_pff * face_sq_pff
    n_pff = dec(cert["channel_norm_monomials"]["t_pf_to_f"]["norm"])
    n_ppp = dec(cert["channel_norm_monomials"]["t_pp_to_p"]["norm"])
    check(laurent_eval(dG, 1, 1, 1) + F5(5) * face_sq_pff * n_pff == dGc, "control d_G")
    wA = n_ppp
    wB = F5(6) * F5(5) * n_pff
    pA = RT5 * F5(Fraction(1, 20))
    pB = -F5(Fraction(1, 20))
    a_star = (wA * pA + wB * pB) / (wA + wB)
    K_c = wA * (a_star - pA) * (a_star - pA) + wB * (a_star - pB) * (a_star - pB)
    base = {key: value for key, value in dF.items() if key != "inv_beta"}
    check(laurent_eval(base, 1, 1, 1) + K_c == dFc, "control d_F")
    check(sgn(dGc - dFc) > 0, "control does not reverse the selection")

    return {
        "commutant_dimension": facts["nullity"],
        "samples_replayed": len(cert["samples"]),
        "P_excluded_everywhere": True,
        "balanced_slice_selects": "G",
        "phase_box": "beta/delta in [1/50, 6] selects G for all gamma and delta",
        "F_witness": "(8, 1, 1)",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, default=DEFAULT)
    args = parser.parse_args()
    summary = verify_certificate(args.certificate)
    print("invariant-metric phase certificate verified")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
