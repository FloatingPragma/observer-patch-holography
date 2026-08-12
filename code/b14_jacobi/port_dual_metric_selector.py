#!/usr/bin/env python3
"""Measure-conditional metric selection for issue B14 / #705.

The committed dual-measure module
``Lean/Screen/PrimitivePortDualMeasure.lean`` proves that on the pinned
twelve-port incidence table, equal normalized weights on the twenty
reconstructed triangular faces and the declared barycentric one-third
vertex share give port-dual weight exactly ``1/12`` at every port.  The
committed phase module ``Lean/Screen/OrientedFaceInvariantMetric.lean``
proves that every sector-balanced invariant carrier metric selects the
compact family ``G`` as unique nearest to the pinned oriented-face
bracket.  This producer certifies the bridge between the two committed
layers on the pinned data:

  1. UNIFORMITY      the port-dual measure recomputed from the pinned
                     manifest is exactly ``1/12`` at every port (twenty
                     faces from the pinned edge table, five incident
                     faces per port, barycentric share ``1/3``);
  2. METRIC          under the declared diagonal measure-to-metric rule
                     the measure induces the carrier inner product
                     ``diag(w) = (1/12) I``, invariant under all sixty
                     pinned proper port permutations;
  3. CONE POINT      that metric lies in the positive sector-scale cone
                     with coordinates ``(alpha, beta, gamma, delta) =
                     (1/12, 1/12, 1/12, 1/12)``, verified against the
                     four exact spectral projectors: it is
                     sector-balanced (``beta = gamma``) and lies in the
                     committed phase box (``delta <= 50 beta``,
                     ``beta <= 6 delta``);
  4. SELECTION       the exact squared distances at that point are
                     ``d_G^2 = (3690 - 738 sqrt5)/11 <
                     d_F^2 = (3690 + 738 sqrt5)/11 < d_P^2 = 540``, so
                     the pinned measure point selects ``G`` with exact
                     positive gaps, and ``c * d^2(c,c,c) = d^2(1,1,1)``
                     reproduces the pinned selector reference distances
                     at ``c = 1/12``;
  5. CONTROLS        the degenerate five-sector spectral projector is
                     rejected by the cone-membership requirement; a
                     non-equal face weighting produces a non-uniform
                     port measure whose diagonal metric has no
                     sector-scale coordinates (the equal-weight premise
                     is load-bearing); the committed F-witness
                     ``(8, 1, 1)`` keeps the phase diagram nonvacuous,
                     so landing on the balanced locus carries content.

BOUNDARY: the equal-face-weight rule and the barycentric one-third share
are the declared symmetric weighting premises of the committed
dual-measure module, and the diagonal measure-to-metric rule and the
nearest-point repair rule are declared discriminator choices; none of
the four is derived from the source.  The comparison ranges over the
compact locus classified by the pinned ``b14_compact_locus``
certificate, conditional on its three named textbook lemmas.  The
result upgrades the committed metric-conditional selection to a
measure-conditional selection with the measure pinned by a committed
Lean theorem.  No physical current, coupling, laboratory gauge field,
or source bracket selection is claimed.

Companion Lean corollary: ``Lean/Screen/PortDualMetricSelection.lean``.

Usage (from the repository root):

    python3 code/b14_jacobi/port_dual_metric_selector.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
MANIFEST = REPO / "code/a5_closure/manifests/echosahedral_federation_reference.json"
RECEIPT = REPO / "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json"
SELECTOR = HERE / "oriented_face_bracket_selector.certificate.json"
PHASE = HERE / "invariant_metric_phase.certificate.json"
OUTPUT = HERE / "port_dual_metric_selector.certificate.json"
VERIFIER = HERE / "verify_port_dual_metric_selector.py"

N = 12
SECTORS = ["fixed", "three_plus", "three_minus", "five"]
SECTOR_DIMS = {"fixed": 1, "three_plus": 3, "three_minus": 3, "five": 5}


def require(value: bool, message: str) -> None:
    if not value:
        raise SystemExit(f"FAILED: {message}")


class Q5:
    """Exact arithmetic in Q(sqrt5): a + b*sqrt5 with rational a, b."""

    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = Fraction(a)
        self.b = Fraction(b)

    def __add__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return Q5(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q5(-self.a, -self.b)

    def __sub__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return Q5(self.a - other.a, self.b - other.b)

    def __rsub__(self, other):
        return Q5(other) - self

    def __mul__(self, other):
        other = other if isinstance(other, Q5) else Q5(other)
        return Q5(self.a * other.a + 5 * self.b * other.b,
                  self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def reciprocal(self):
        norm = self.a * self.a - 5 * self.b * self.b
        require(norm != 0, "division by zero in Q5")
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
        return self.a != 0 or self.b != 0

    def __repr__(self):
        return f"Q5({self.a},{self.b})"


ZERO = Q5(0)
ONE = Q5(1)
ROOT5 = Q5(0, 1)


def sign(x: Q5) -> int:
    """Exact sign of a + b*sqrt5."""
    a, b = x.a, x.b
    if a == 0 and b == 0:
        return 0
    if a >= 0 and b >= 0:
        return 1
    if a <= 0 and b <= 0:
        return -1
    lhs, rhs = a * a, 5 * b * b
    if a > 0:
        return 1 if lhs > rhs else -1
    return 1 if rhs > lhs else -1


def enc(x: Q5):
    return [str(x.a), str(x.b)]


def matmul(X, Y):
    n = len(X)
    return [
        [sum((X[i][k] * Y[k][j] for k in range(n)), ZERO) for j in range(n)]
        for i in range(n)
    ]


def identity(n):
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_bytes(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True).encode("ascii")


def object_hash(value) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def load_carrier():
    manifest = json.loads(MANIFEST.read_text())
    carrier = manifest["carrier"]
    ports = carrier["ports"]
    require(len(ports) == N, "expected twelve ports")
    index = {name: i for i, name in enumerate(ports)}
    edges = {
        frozenset((index[a], index[b])) for a, b in carrier["edges"]
    }
    require(len(edges) == 30, "expected thirty edges")
    return index, edges


def reconstruct_faces(edges):
    """Triangles of the pinned edge table, one canonical triple each."""
    faces = []
    for i in range(N):
        for j in range(i + 1, N):
            if frozenset((i, j)) not in edges:
                continue
            for k in range(j + 1, N):
                if (frozenset((i, k)) in edges
                        and frozenset((j, k)) in edges):
                    faces.append((i, j, k))
    return faces


def port_dual_weights(faces):
    """Equal normalized face weights, barycentric one-third share."""
    face_weight = Fraction(1, len(faces))
    share = Fraction(1, 3)
    weights = [Fraction(0)] * N
    for face in faces:
        for p in face:
            weights[p] += share * face_weight
    return weights


def group_rows():
    receipt = json.loads(RECEIPT.read_text())
    group = [tuple(row) for row in receipt["proper_port_action"]["permutation_rows"]]
    require(len(group) == 60 and len(set(group)) == 60,
            "expected sixty distinct group rows")
    for row in group:
        require(sorted(row) == list(range(N)), "group row is not a permutation")
    return group


def spectral_projectors(group):
    """The four sector projectors from the lexicographically minimal
    60-orbit of the pair action (the same convention as the committed
    phase producer, so sector labels agree across the two layers)."""
    unassigned = {(i, j) for i in range(N) for j in range(N)}
    orbits = []
    while unassigned:
        seed = min(unassigned)
        orbit = {(g[seed[0]], g[seed[1]]) for g in group}
        unassigned -= orbit
        orbits.append(orbit)
    require(sorted(len(o) for o in orbits) == [12, 12, 60, 60],
            "unexpected ordered-pair orbit sizes")
    orbital = sorted((o for o in orbits if len(o) == 60),
                     key=lambda o: tuple(sorted(o)))[0]
    require(all((j, i) in orbital for (i, j) in orbital),
            "orbital adjacency is not symmetric")
    A = [[ZERO] * N for _ in range(N)]
    for i, j in orbital:
        A[i][j] = ONE
    eigenvalues = {"fixed": Q5(5), "three_plus": ROOT5,
                   "three_minus": -ROOT5, "five": Q5(-1)}
    projectors = {}
    for label, ev in eigenvalues.items():
        Mx = identity(N)
        for label2, ev2 in eigenvalues.items():
            if label2 == label:
                continue
            factor = [[A[i][j] - (ev2 if i == j else ZERO)
                       for j in range(N)] for i in range(N)]
            scale = (ev - ev2).reciprocal()
            Mx = [[v * scale for v in row] for row in matmul(Mx, factor)]
        projectors[label] = Mx
    total = [[sum((projectors[l][i][j] for l in SECTORS), ZERO)
              for j in range(N)] for i in range(N)]
    require(total == identity(N), "spectral projectors do not sum to identity")
    for label in SECTORS:
        P = projectors[label]
        require(matmul(P, P) == P, f"{label} projector is not idempotent")
        trace = sum((P[i][i] for i in range(N)), ZERO)
        require(trace == Q5(SECTOR_DIMS[label]), f"{label} projector trace differs")
    return projectors


def sector_scales(metric, projectors):
    """Coordinates of a diagonal metric in the sector-scale cone.

    Returns the four exact scales, or None when the metric is not a
    positive combination of the four projectors."""
    scales = {}
    for label in SECTORS:
        P = projectors[label]
        prod = matmul(metric, P)
        candidate = None
        for i in range(N):
            for j in range(N):
                if P[i][j]:
                    candidate = prod[i][j] / P[i][j]
                    break
            if candidate is not None:
                break
        scaled = [[candidate * v for v in row] for row in P]
        if prod != scaled:
            return None
        if sign(candidate) <= 0:
            return None
        scales[label] = candidate
    combo = [[sum((scales[l] * projectors[l][i][j] for l in SECTORS), ZERO)
              for j in range(N)] for i in range(N)]
    if combo != metric:
        return None
    return scales


def d_squared(beta: Q5, gamma: Q5, delta: Q5):
    """The three committed Laurent closed forms."""
    d2 = delta * delta
    dP = ((15 - 3 * ROOT5) / beta + (15 + 3 * ROOT5) / gamma
          + ((15 - 3 * ROOT5) / 2) * beta / d2
          + ((15 + 3 * ROOT5) / 2) * gamma / d2)
    dF = (((15 + 3 * ROOT5) / 2) * gamma / d2 + (15 + 3 * ROOT5) / gamma
          + (60 + 12 * ROOT5) / (11 * beta))
    dG = (((15 - 3 * ROOT5) / 2) * beta / d2 + (15 - 3 * ROOT5) / beta
          + (60 - 12 * ROOT5) / (11 * gamma))
    return {"P": dP, "F": dF, "G": dG}


def build_certificate() -> dict:
    for path in (MANIFEST, RECEIPT, SELECTOR, PHASE):
        require(path.is_file(), f"missing pinned input {path.name}")

    index, edges = load_carrier()
    faces = reconstruct_faces(edges)
    require(len(faces) == 20, "expected twenty triangular faces")
    incident_counts = [sum(1 for f in faces if p in f) for p in range(N)]
    require(incident_counts == [5] * N, "expected five incident faces per port")

    weights = port_dual_weights(faces)
    require(all(w == Fraction(1, 12) for w in weights),
            "port-dual measure is not uniform 1/12")
    require(sum(weights) == 1, "port-dual measure is not normalized")

    # The declared diagonal measure-to-metric rule.
    metric = [[Q5(weights[i]) if i == j else ZERO for j in range(N)]
              for i in range(N)]

    group = group_rows()
    for g in group:
        moved = [[metric[g[i]][g[j]] for j in range(N)] for i in range(N)]
        require(moved == metric, "metric is not invariant under the port action")

    projectors = spectral_projectors(group)
    scales = sector_scales(metric, projectors)
    require(scales is not None, "metric has no positive sector-scale coordinates")
    require(all(scales[l] == Q5(Fraction(1, 12)) for l in SECTORS),
            "sector scales are not uniformly 1/12")
    beta, gamma, delta = scales["three_plus"], scales["three_minus"], scales["five"]
    require(beta == gamma, "metric is not sector-balanced")
    require(sign(Q5(50) * beta - delta) >= 0 and sign(Q5(6) * delta - beta) >= 0,
            "metric lies outside the committed phase box")

    distances = d_squared(beta, gamma, delta)
    require(distances["G"] == Q5(Fraction(3690, 11), Fraction(-738, 11)),
            "d_G^2 at the measure point differs from the pinned value")
    require(distances["F"] == Q5(Fraction(3690, 11), Fraction(738, 11)),
            "d_F^2 at the measure point differs from the pinned value")
    require(distances["P"] == Q5(540),
            "d_P^2 at the measure point differs from the pinned value")
    gap_FG = distances["F"] - distances["G"]
    gap_PG = distances["P"] - distances["G"]
    require(sign(gap_FG) > 0 and sign(gap_PG) > 0,
            "G is not the strict unique nearest family at the measure point")

    # Scaling receipt against the pinned selector reference distances.
    selector = json.loads(SELECTOR.read_text())
    fams = selector["orthogonal_compact_locus_discriminator"]["families"]
    ref = {}
    for family, entry in fams.items():
        an, ad, bn, bd = entry["squared_distance"]
        ref[family] = Q5(Fraction(an, ad), Fraction(bn, bd))
    c = Q5(Fraction(1, 12))
    for family in ("P", "F", "G"):
        require(c * distances[family] == ref[family],
                f"scaling receipt fails for {family}")

    # Control 1: a degenerate candidate is rejected by cone membership.
    require(sector_scales(projectors["five"], projectors) is None,
            "degenerate projector control unexpectedly passed")

    # Control 2: a non-equal face weighting has no sector-scale coordinates.
    skewed = [Fraction(0)] * N
    total = Fraction(0)
    for k, face in enumerate(faces):
        fw = Fraction(2 if k == 0 else 1, len(faces) + 1)
        total += fw
        for p in face:
            skewed[p] += Fraction(1, 3) * fw
    require(total == 1, "skewed control weighting is not normalized")
    require(len(set(skewed)) > 1, "skewed control weighting is uniform")
    skewed_metric = [[Q5(skewed[i]) if i == j else ZERO for j in range(N)]
                     for i in range(N)]
    require(sector_scales(skewed_metric, projectors) is None,
            "skewed-measure control unexpectedly has sector scales")

    # Control 3: the committed F-witness keeps the phase diagram nonvacuous.
    witness = d_squared(Q5(8), ONE, ONE)
    require(sign(witness["G"] - witness["F"]) > 0,
            "F-witness control fails: F is not nearer at (8,1,1)")

    certificate = {
        "schema": "oph.b14_port_dual_metric_selector.v1",
        "issue": "https://github.com/muellerberndt/reverse-engineering-reality/issues/705",
        "pinned_inputs": {
            "echosahedral_federation_reference.json": file_hash(MANIFEST),
            "a5_alternating_bracket_space_stage1.receipt.json": file_hash(RECEIPT),
            "oriented_face_bracket_selector.certificate.json": file_hash(SELECTOR),
            "invariant_metric_phase.certificate.json": file_hash(PHASE),
        },
        "declared_rules": {
            "face_weighting": "equal normalized weights on the twenty faces",
            "vertex_share": "1/3",
            "measure_to_metric": (
                "the normalized port-dual measure is read as the diagonal "
                "carrier inner product it induces on port functions"
            ),
            "repair_rule": "nearest classified compact family",
        },
        "port_dual_weight": "1/12",
        "sector_scales": {label: enc(scales[label]) for label in SECTORS},
        "sector_balanced": True,
        "phase_box_member": True,
        "squared_distances_at_measure_point": {
            family: enc(distances[family]) for family in ("P", "F", "G")
        },
        "selection": "G",
        "gaps": {"F_minus_G": enc(gap_FG), "P_minus_G": enc(gap_PG)},
        "scaling_receipt": (
            "c * d^2(c,c,c) = d^2(1,1,1) verified at c = 1/12 against the "
            "pinned selector reference distances for P, F, G"
        ),
        "controls": {
            "degenerate_projector_rejected": True,
            "skewed_face_weighting_rejected": True,
            "f_witness_nonvacuous": True,
        },
        "lean_companions": [
            "Lean/Screen/PrimitivePortDualMeasure.lean",
            "Lean/Screen/OrientedFaceInvariantMetric.lean",
            "Lean/Screen/PortDualMetricSelection.lean",
        ],
        "boundary": (
            "The equal-face-weight rule, barycentric one-third share, "
            "diagonal measure-to-metric rule, and nearest-point repair rule "
            "are declared, and the compact locus is conditional on its three "
            "named textbook lemmas; the certificate upgrades the committed "
            "metric-conditional selection to a measure-conditional selection "
            "whose measure is pinned by a committed Lean theorem. No "
            "physical current, coupling, laboratory gauge field, or source "
            "bracket selection is claimed."
        ),
        "implementation_pins": {
            "producer_sha256": file_hash(Path(__file__).resolve()),
            "verifier_sha256": file_hash(VERIFIER),
        },
    }
    certificate["certificate_sha256"] = object_hash(certificate)
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true",
                        help="require the committed certificate to match "
                             "a fresh rebuild instead of overwriting it")
    args = parser.parse_args()
    certificate = build_certificate()
    text = json.dumps(certificate, indent=1, sort_keys=True) + "\n"
    if args.check:
        require(args.output.exists(), f"missing certificate {args.output}")
        require(json.loads(args.output.read_text()) == certificate,
                "committed certificate is stale")
        print(f"CHECK PASSED: {args.output.name} matches a fresh rebuild")
    else:
        # write_bytes keeps LF endings on every platform, so fresh
        # production stays byte-identical to the committed certificate.
        args.output.write_bytes(text.encode("utf-8"))
        print(f"wrote {args.output} ({len(text.encode())} bytes)")
    print(f"certificate_sha256: {certificate['certificate_sha256']}")
    print("selection: G at sector scales (1/12, 1/12, 1/12, 1/12); "
          "d_G^2 = (3690 - 738*sqrt5)/11, d_F^2 = (3690 + 738*sqrt5)/11, "
          "d_P^2 = 540")


if __name__ == "__main__":
    main()
