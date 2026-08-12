#!/usr/bin/env python3
"""Independent exact replay of the B14 port-dual metric-selector certificate.

No code is imported from the producer, and the two computations take
different routes wherever the pinned data allows: the faces are read
from the manifest's ``oriented_faces`` table and cross-checked against
an edge-triangle reconstruction (the producer builds them from edges
alone); the sector scales are extracted through trace ratios (the
producer matches scaled projectors entrywise); and the squared distances
are evaluated from the phase certificate's serialized Laurent
coefficient tables (the producer evaluates its own closed forms).  The
spectral projectors are rebuilt by the same interpolation convention as
the producer, so they are a re-typed common route rather than an
independent one.  Every certificate field is recomputed and compared,
the three controls are replayed, the pinned input hashes and
implementation pins are required to match, the top-level key set is
closed, and the load-bearing prose (the declared-rules block and the
boundary paragraph) is bound to exact expected content together with a
canonical self-hash of the whole body.

Usage:

    python3 code/b14_jacobi/verify_port_dual_metric_selector.py [certificate]
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT = HERE / "port_dual_metric_selector.certificate.json"
MANIFEST = REPO / "code/a5_closure/manifests/echosahedral_federation_reference.json"
RECEIPT = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json"
SELECTOR = HERE / "oriented_face_bracket_selector.certificate.json"
PHASE = HERE / "invariant_metric_phase.certificate.json"

N = 12
SECTORS = ["fixed", "three_plus", "three_minus", "five"]
SECTOR_DIMS = {"fixed": 1, "three_plus": 3, "three_minus": 3, "five": 5}

EXPECTED_TOP_LEVEL_KEYS = {
    "schema", "issue", "pinned_inputs", "declared_rules",
    "port_dual_weight", "sector_scales", "sector_balanced",
    "phase_box_member", "squared_distances_at_measure_point", "selection",
    "gaps", "scaling_receipt", "controls", "lean_companions", "boundary",
    "implementation_pins", "certificate_sha256",
}

EXPECTED_DECLARED_RULES = {
    "face_weighting": "equal normalized weights on the twenty faces",
    "vertex_share": "1/3",
    "measure_to_metric": (
        "the normalized port-dual measure is read as the diagonal "
        "carrier inner product it induces on port functions"
    ),
    "repair_rule": "nearest classified compact family",
}

EXPECTED_BOUNDARY = (
    "The equal-face-weight rule, barycentric one-third share, "
    "diagonal measure-to-metric rule, and nearest-point repair rule "
    "are declared, and the compact locus is conditional on its three "
    "named textbook lemmas; the certificate upgrades the committed "
    "metric-conditional selection to a measure-conditional selection "
    "whose measure is pinned by a committed Lean theorem. No "
    "physical current, coupling, laboratory gauge field, or source "
    "bracket selection is claimed."
)

EXPECTED_ISSUE = (
    "https://github.com/muellerberndt/reverse-engineering-reality/issues/705"
)


class VerificationError(RuntimeError):
    pass


def check(value: bool, message: str) -> None:
    if not value:
        raise VerificationError(message)


class F5:
    """Exact r + s*sqrt5, independent implementation."""

    __slots__ = ("r", "s")

    def __init__(self, r=0, s=0):
        self.r = Fraction(r)
        self.s = Fraction(s)

    @staticmethod
    def of(pair) -> "F5":
        return F5(Fraction(pair[0]), Fraction(pair[1]))

    @staticmethod
    def of4(quad) -> "F5":
        return F5(Fraction(quad[0], quad[1]), Fraction(quad[2], quad[3]))

    def __add__(self, o):
        o = o if isinstance(o, F5) else F5(o)
        return F5(self.r + o.r, self.s + o.s)

    __radd__ = __add__

    def __sub__(self, o):
        o = o if isinstance(o, F5) else F5(o)
        return F5(self.r - o.r, self.s - o.s)

    def __mul__(self, o):
        o = o if isinstance(o, F5) else F5(o)
        return F5(self.r * o.r + 5 * self.s * o.s,
                  self.r * o.s + self.s * o.r)

    __rmul__ = __mul__

    def inv(self):
        norm = self.r * self.r - 5 * self.s * self.s
        check(norm != 0, "inverse of zero in F5")
        return F5(self.r / norm, -self.s / norm)

    def __eq__(self, o):
        o = o if isinstance(o, F5) else F5(o)
        return self.r == o.r and self.s == o.s

    def __hash__(self):
        return hash((self.r, self.s))

    def positive(self) -> bool:
        r, s = self.r, self.s
        if r == 0 and s == 0:
            return False
        if r >= 0 and s >= 0:
            return True
        if r <= 0 and s <= 0:
            return False
        if r > 0:
            return r * r > 5 * s * s
        return 5 * s * s > r * r


def enc(x: F5):
    return [str(x.r), str(x.s)]


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_bytes(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True).encode("ascii")


def object_hash(value) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def load_faces():
    carrier = json.loads(MANIFEST.read_text())["carrier"]
    ports = carrier["ports"]
    check(len(ports) == N, "expected twelve ports")
    index = {name: i for i, name in enumerate(ports)}
    declared = {
        frozenset(index[name] for name in face)
        for face in carrier["oriented_faces"]
    }
    check(len(declared) == 20, "expected twenty declared faces")
    check(all(len(f) == 3 for f in declared), "declared face is not a triangle")
    edges = {frozenset((index[a], index[b])) for a, b in carrier["edges"]}
    check(len(edges) == 30, "expected thirty edges")
    triangles = set()
    for i in range(N):
        for j in range(i + 1, N):
            if frozenset((i, j)) not in edges:
                continue
            for k in range(j + 1, N):
                if (frozenset((i, k)) in edges
                        and frozenset((j, k)) in edges):
                    triangles.add(frozenset((i, j, k)))
    check(triangles == declared,
          "declared oriented faces disagree with the edge triangles")
    return sorted(tuple(sorted(f)) for f in declared)


def load_group():
    receipt = json.loads(RECEIPT.read_text())
    group = [tuple(row) for row in receipt["proper_port_action"]["permutation_rows"]]
    check(len(group) == 60 and len(set(group)) == 60, "expected sixty rows")
    for row in group:
        check(sorted(row) == list(range(N)), "row is not a permutation")
    return group


def projectors_from(group):
    unassigned = {(i, j) for i in range(N) for j in range(N)}
    orbits = []
    while unassigned:
        seed = min(unassigned)
        orbit = {(g[seed[0]], g[seed[1]]) for g in group}
        unassigned -= orbit
        orbits.append(orbit)
    check(sorted(len(o) for o in orbits) == [12, 12, 60, 60],
          "unexpected orbit sizes")
    orbital = sorted((o for o in orbits if len(o) == 60),
                     key=lambda o: tuple(sorted(o)))[0]
    A = [[F5(0)] * N for _ in range(N)]
    for i, j in orbital:
        A[i][j] = F5(1)

    def matmul(X, Y):
        return [
            [sum((X[i][k] * Y[k][j] for k in range(N)), F5(0))
             for j in range(N)]
            for i in range(N)
        ]

    eigenvalues = {"fixed": F5(5), "three_plus": F5(0, 1),
                   "three_minus": F5(0, -1), "five": F5(-1)}
    projectors = {}
    for label, ev in eigenvalues.items():
        Mx = [[F5(1) if i == j else F5(0) for j in range(N)] for i in range(N)]
        for label2, ev2 in eigenvalues.items():
            if label2 == label:
                continue
            factor = [[A[i][j] - (ev2 if i == j else F5(0))
                       for j in range(N)] for i in range(N)]
            scale = (ev - ev2).inv()
            Mx = [[v * scale for v in row] for row in matmul(Mx, factor)]
        projectors[label] = Mx
    for label in SECTORS:
        P = projectors[label]
        check(matmul(P, P) == P, f"{label} projector is not idempotent")
        check(sum((P[i][i] for i in range(N)), F5(0)) == F5(SECTOR_DIMS[label]),
              f"{label} projector trace differs")
    return projectors, matmul


def trace_scales(metric, projectors, matmul):
    """Scales via trace ratios, then a full reconstruction equality."""
    scales = {}
    for label in SECTORS:
        P = projectors[label]
        prod = matmul(metric, P)
        tr = sum((prod[i][i] for i in range(N)), F5(0))
        scales[label] = tr * F5(SECTOR_DIMS[label]).inv()
    combo = [
        [sum((scales[l] * projectors[l][i][j] for l in SECTORS), F5(0))
         for j in range(N)]
        for i in range(N)
    ]
    if combo != metric:
        return None
    if not all(scales[l].positive() for l in SECTORS):
        return None
    return scales


def laurent_distances(beta: F5, gamma: F5, delta: F5):
    """Evaluate the phase certificate's serialized Laurent tables."""
    forms = json.loads(PHASE.read_text())["closed_forms"]
    values = {}
    monomials = {
        "beta_over_delta_sq": beta * (delta * delta).inv(),
        "gamma_over_delta_sq": gamma * (delta * delta).inv(),
        "inv_beta": beta.inv(),
        "inv_gamma": gamma.inv(),
    }
    for family, key in (("P", "d_P_squared"), ("F", "d_F_squared"),
                        ("G", "d_G_squared")):
        total = F5(0)
        for mono, coeff in forms[key].items():
            total = total + F5.of4(coeff) * monomials[mono]
        values[family] = total
    return values


def main() -> None:
    cert_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT
    cert = json.loads(cert_path.read_text())
    check(cert.get("schema") == "oph.b14_port_dual_metric_selector.v1",
          "unexpected schema")
    check(set(cert) == EXPECTED_TOP_LEVEL_KEYS,
          "unexpected top-level certificate key set")
    body = {k: v for k, v in cert.items() if k != "certificate_sha256"}
    check(object_hash(body) == cert["certificate_sha256"],
          "certificate self-hash mismatch")
    check(cert["issue"] == EXPECTED_ISSUE, "issue reference mismatch")
    check(cert["declared_rules"] == EXPECTED_DECLARED_RULES,
          "declared-rules block deviates from the pinned content")
    check(cert["boundary"] == EXPECTED_BOUNDARY,
          "boundary paragraph deviates from the pinned content")
    imp = cert["implementation_pins"]
    producer = HERE / "port_dual_metric_selector.py"
    check(imp["producer_sha256"] == file_hash(producer),
          "producer implementation pin mismatch")
    check(imp["verifier_sha256"] == file_hash(Path(__file__).resolve()),
          "verifier implementation pin mismatch")
    pins = cert["pinned_inputs"]
    expected = {
        "echosahedral_federation_reference.json": MANIFEST,
        "a5_alternating_bracket_space_stage1.receipt.json": RECEIPT,
        "oriented_face_bracket_selector.certificate.json": SELECTOR,
        "invariant_metric_phase.certificate.json": PHASE,
    }
    check(set(pins) == set(expected), "pinned input set mismatch")
    for name, path in expected.items():
        check(pins[name] == file_hash(path), f"pinned hash mismatch: {name}")

    faces = load_faces()
    counts = [sum(1 for f in faces if p in f) for p in range(N)]
    check(counts == [5] * N, "expected five incident faces per port")
    weights = [Fraction(0)] * N
    for face in faces:
        for p in face:
            weights[p] += Fraction(1, 3) * Fraction(1, len(faces))
    check(all(w == Fraction(1, 12) for w in weights),
          "port-dual measure is not uniform 1/12")
    check(cert["port_dual_weight"] == "1/12", "certificate weight field")

    metric = [[F5(weights[i]) if i == j else F5(0) for j in range(N)]
              for i in range(N)]
    group = load_group()
    for g in group:
        moved = [[metric[g[i]][g[j]] for j in range(N)] for i in range(N)]
        check(moved == metric, "metric not invariant under the port action")

    projectors, matmul = projectors_from(group)
    scales = trace_scales(metric, projectors, matmul)
    check(scales is not None, "metric has no positive sector-scale coordinates")
    for label in SECTORS:
        check(enc(scales[label]) == cert["sector_scales"][label],
              f"sector scale mismatch: {label}")
        check(scales[label] == F5(Fraction(1, 12)),
              f"sector scale is not 1/12: {label}")
    beta, gamma, delta = (scales["three_plus"], scales["three_minus"],
                          scales["five"])
    check(beta == gamma, "not sector-balanced")
    check((F5(50) * beta - delta).positive() or F5(50) * beta == delta,
          "phase box lower bound fails")
    check((F5(6) * delta - beta).positive() or F5(6) * delta == beta,
          "phase box upper bound fails")
    check(cert["sector_balanced"] is True and cert["phase_box_member"] is True,
          "certificate cone flags")

    distances = laurent_distances(beta, gamma, delta)
    for family in ("P", "F", "G"):
        check(enc(distances[family])
              == cert["squared_distances_at_measure_point"][family],
              f"distance mismatch: {family}")
    gap_FG = distances["F"] - distances["G"]
    gap_PG = distances["P"] - distances["G"]
    check(gap_FG.positive() and gap_PG.positive(),
          "G is not strictly nearest at the measure point")
    check(enc(gap_FG) == cert["gaps"]["F_minus_G"], "F gap mismatch")
    check(enc(gap_PG) == cert["gaps"]["P_minus_G"], "P gap mismatch")
    check(cert["selection"] == "G", "certificate selection field")

    fams = json.loads(SELECTOR.read_text())[
        "orthogonal_compact_locus_discriminator"]["families"]
    c = F5(Fraction(1, 12))
    for family in ("P", "F", "G"):
        an, ad, bn, bd = fams[family]["squared_distance"]
        check(c * distances[family] == F5(Fraction(an, ad), Fraction(bn, bd)),
              f"scaling receipt fails: {family}")

    check(trace_scales(projectors["five"], projectors, matmul) is None,
          "degenerate projector control unexpectedly passed")
    skewed = [Fraction(0)] * N
    for k, face in enumerate(faces):
        fw = Fraction(2 if k == 0 else 1, len(faces) + 1)
        for p in face:
            skewed[p] += Fraction(1, 3) * fw
    skewed_metric = [[F5(skewed[i]) if i == j else F5(0) for j in range(N)]
                     for i in range(N)]
    check(trace_scales(skewed_metric, projectors, matmul) is None,
          "skewed-measure control unexpectedly has sector scales")
    witness = laurent_distances(F5(8), F5(1), F5(1))
    check((witness["G"] - witness["F"]).positive(),
          "F-witness control fails at (8,1,1)")
    controls = cert["controls"]
    check(controls == {"degenerate_projector_rejected": True,
                       "skewed_face_weighting_rejected": True,
                       "f_witness_nonvacuous": True},
          "certificate controls block mismatch")

    print("VERIFICATION PASSED:", cert_path.name)
    print("  uniform 1/12 port-dual measure -> balanced cone point "
          "(1/12, 1/12, 1/12, 1/12) -> unique nearest family G, "
          "replayed through independent routes")


if __name__ == "__main__":
    main()
