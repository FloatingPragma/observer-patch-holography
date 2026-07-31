#!/usr/bin/env python3
"""Exact certificate for GitHub issue #566: the physical port-current algebra.

The input is a port-current response manifest.  It declares only:

* the certified twelve-port echosahedral carrier manifest (by path and hash);
* four signed A5-equivariant response coefficients, as exact rationals;
* the typing split between reversible response automorphisms (current
  sources) and irreversible strict-descent repairs (never current sources).

From that packet plus the carrier source data the verifier derives, rather
than assumes:

* an oriented realization of the canonical rank-three port frame, unique up
  to the proper action, with exact vertex vectors in Q(sqrt5);
* an injective port-to-generator map K : P_12 -> u(H) with twelve-dimensional
  image on a faithful charged response space H;
* exact skew-adjointness, commutator closure, compact type, the
  one-dimensional central u(1), derived dimension eleven, and adjoint rank
  eleven;
* A5 covariance, the icosahedral intertwiner, and innerness of the induced
  A5 action through sixty exact rotation normal-form witnesses;
* refinement naturality along the declared carrier tower;
* the four-dimensional equivariant response moduli (four signed band
  coefficients are exactly the declared algebraic freedom).

Abelian-record and rank-deficient response models fail the conditional
algebraic gate with typed error codes.  No Standard Model representation, particle
assignment, measured coupling, or gauge target is accepted in a source
manifest.  Every arithmetic decision is exact in Q(sqrt5); no floating point
appears in a proof step.
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402

SCHEMA = "oph.port_current_response_manifest.v5"
RECEIPT_SCHEMA = "oph.port_current_inner_receipt.v5"
NEGATIVE_SCHEMA = "oph.port_current_inner_negative_controls.v5"
ARTIFACT_SCHEMA = "oph.charged_response_semantic_artifact.v3"

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json
compose = e565.compose
inverse = e565.inverse


# ---------------------------------------------------------------------------
# Exact arithmetic in Q(sqrt5)
# ---------------------------------------------------------------------------


class F5:
    """An element a + b*sqrt(5) of Q(sqrt5) with exact Fraction coefficients."""

    __slots__ = ("a", "b")

    def __init__(self, a: Any = 0, b: Any = 0) -> None:
        self.a = Fraction(a)
        self.b = Fraction(b)

    def __add__(self, other: "F5") -> "F5":
        return F5(self.a + other.a, self.b + other.b)

    def __sub__(self, other: "F5") -> "F5":
        return F5(self.a - other.a, self.b - other.b)

    def __neg__(self) -> "F5":
        return F5(-self.a, -self.b)

    def __mul__(self, other: "F5") -> "F5":
        return F5(
            self.a * other.a + 5 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    def inv(self) -> "F5":
        norm = self.a * self.a - 5 * self.b * self.b
        if norm == 0:
            raise ZeroDivisionError("zero element of Q(sqrt5)")
        return F5(self.a / norm, -self.b / norm)

    def __truediv__(self, other: "F5") -> "F5":
        return self * other.inv()

    def conj(self) -> "F5":
        """The Galois conjugate sqrt5 -> -sqrt5."""

        return F5(self.a, -self.b)

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def is_positive(self) -> bool:
        if self.b == 0:
            return self.a > 0
        if self.a == 0:
            return self.b > 0
        if self.a > 0 and self.b > 0:
            return True
        if self.a < 0 and self.b < 0:
            return False
        if self.a > 0:  # b < 0
            return self.a * self.a > 5 * self.b * self.b
        return 5 * self.b * self.b > self.a * self.a  # a < 0, b > 0

    def __eq__(self, other: object) -> bool:
        return isinstance(other, F5) and self.a == other.a and self.b == other.b

    def __hash__(self) -> int:
        return hash((self.a, self.b))

    def __repr__(self) -> str:
        return f"F5({self.a},{self.b})"

    def text(self) -> str:
        if self.b == 0:
            return str(self.a)
        if self.a == 0:
            return f"{self.b}*sqrt(5)"
        return f"{self.a} + {self.b}*sqrt(5)"


ZERO = F5(0)
ONE = F5(1)
PHI = F5(Fraction(1, 2), Fraction(1, 2))  # golden ratio (1+sqrt5)/2
VERTEX_NORM = F5(2) + PHI  # squared vertex norm 2+phi


class C5:
    """An element re + i*im with re, im in Q(sqrt5)."""

    __slots__ = ("re", "im")

    def __init__(self, re: F5 = ZERO, im: F5 = ZERO) -> None:
        self.re = re
        self.im = im

    def __add__(self, other: "C5") -> "C5":
        return C5(self.re + other.re, self.im + other.im)

    def __sub__(self, other: "C5") -> "C5":
        return C5(self.re - other.re, self.im - other.im)

    def __neg__(self) -> "C5":
        return C5(-self.re, -self.im)

    def __mul__(self, other: "C5") -> "C5":
        return C5(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        )

    def conj(self) -> "C5":
        return C5(self.re, -self.im)

    def is_zero(self) -> bool:
        return self.re.is_zero() and self.im.is_zero()

    def __eq__(self, other: object) -> bool:
        return isinstance(other, C5) and self.re == other.re and self.im == other.im

    def text(self) -> str:
        if self.im.is_zero():
            return self.re.text()
        if self.re.is_zero():
            return f"i*({self.im.text()})"
        return f"({self.re.text()}) + i*({self.im.text()})"


Vec3 = tuple[F5, F5, F5]
RMat = list[list[F5]]
CMat = list[list[C5]]


def rzeros(n: int, m: int) -> RMat:
    return [[ZERO for _ in range(m)] for _ in range(n)]


def czeros(n: int) -> CMat:
    return [[C5() for _ in range(n)] for _ in range(n)]


def rmul(x: RMat, y: RMat) -> RMat:
    rows, inner, cols = len(x), len(y), len(y[0])
    out = rzeros(rows, cols)
    for i in range(rows):
        for k in range(inner):
            xv = x[i][k]
            if xv.is_zero():
                continue
            for j in range(cols):
                out[i][j] = out[i][j] + xv * y[k][j]
    return out


def rtranspose(x: RMat) -> RMat:
    return [[x[i][j] for i in range(len(x))] for j in range(len(x[0]))]


def cmul(x: CMat, y: CMat) -> CMat:
    n = len(x)
    out = czeros(n)
    for i in range(n):
        for k in range(n):
            xv = x[i][k]
            if xv.is_zero():
                continue
            for j in range(n):
                out[i][j] = out[i][j] + xv * y[k][j]
    return out


def cadd(x: CMat, y: CMat) -> CMat:
    return [[x[i][j] + y[i][j] for j in range(len(x))] for i in range(len(x))]


def csub(x: CMat, y: CMat) -> CMat:
    return [[x[i][j] - y[i][j] for j in range(len(x))] for i in range(len(x))]


def cdagger(x: CMat) -> CMat:
    n = len(x)
    return [[x[j][i].conj() for j in range(n)] for i in range(n)]


def commutator(x: CMat, y: CMat) -> CMat:
    return csub(cmul(x, y), cmul(y, x))


def ctrace(x: CMat) -> C5:
    total = C5()
    for i in range(len(x)):
        total = total + x[i][i]
    return total


def c_is_zero(x: CMat) -> bool:
    return all(entry.is_zero() for row in x for entry in row)


def rref(matrix: RMat) -> tuple[RMat, list[int]]:
    """Exact reduced row echelon form over Q(sqrt5); returns pivot columns."""

    m = [row[:] for row in matrix]
    rows = len(m)
    cols = len(m[0]) if rows else 0
    pivots: list[int] = []
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if not m[i][c].is_zero()), None)
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        scale = m[r][c].inv()
        m[r] = [entry * scale for entry in m[r]]
        for i in range(rows):
            if i != r and not m[i][c].is_zero():
                factor = m[i][c]
                m[i] = [m[i][j] - factor * m[r][j] for j in range(cols)]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return m, pivots


def nullspace(matrix: RMat) -> list[list[F5]]:
    """Exact right-nullspace basis of a matrix over Q(sqrt5)."""

    reduced, pivots = rref(matrix)
    cols = len(matrix[0]) if matrix else 0
    free = [c for c in range(cols) if c not in pivots]
    basis: list[list[F5]] = []
    for f in free:
        vec = [ZERO for _ in range(cols)]
        vec[f] = ONE
        for row_index, p in enumerate(pivots):
            vec[p] = -reduced[row_index][f]
        basis.append(vec)
    return basis


def rank(matrix: RMat) -> int:
    return len(rref(matrix)[1])


def det3(m: RMat) -> F5:
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


# ---------------------------------------------------------------------------
# Coordinate model of the icosahedral frame (all entries in Q(sqrt5))
# ---------------------------------------------------------------------------


def standard_vertices() -> list[Vec3]:
    """The twelve unnormalized icosahedron vertices: cyclic perms of (0,±1,±phi)."""

    verts: list[Vec3] = []
    for s1 in (ONE, -ONE):
        for s2 in (PHI, -PHI):
            verts.append((ZERO, s1, s2))
    for s1 in (ONE, -ONE):
        for s2 in (PHI, -PHI):
            verts.append((s1, s2, ZERO))
    for s1 in (ONE, -ONE):
        for s2 in (PHI, -PHI):
            verts.append((s2, ZERO, s1))
    return verts


def dot(u: Vec3, v: Vec3) -> F5:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def vertex_distance(u: Vec3, v: Vec3) -> int:
    value = dot(u, v)
    if value == VERTEX_NORM:
        return 0
    if value == PHI:
        return 1
    if value == -PHI:
        return 2
    if value == -VERTEX_NORM:
        return 3
    raise CertificateError("VERTEX_MODEL", "unexpected vertex inner product")


def hat(x: Vec3) -> RMat:
    return [
        [ZERO, -x[2], x[1]],
        [x[2], ZERO, -x[0]],
        [-x[1], x[0], ZERO],
    ]


def apply3(m: RMat, x: Vec3) -> Vec3:
    return (
        m[0][0] * x[0] + m[0][1] * x[1] + m[0][2] * x[2],
        m[1][0] * x[0] + m[1][1] * x[1] + m[1][2] * x[2],
        m[2][0] * x[0] + m[2][1] * x[1] + m[2][2] * x[2],
    )


def outer(u: Vec3, v: Vec3) -> RMat:
    return [[u[i] * v[j] for j in range(3)] for i in range(3)]


def to_cmat(real: RMat, imag: RMat | None = None) -> CMat:
    n = len(real)
    return [
        [
            C5(real[i][j], imag[i][j] if imag is not None else ZERO)
            for j in range(n)
        ]
        for i in range(n)
    ]


# ---------------------------------------------------------------------------
# Manifest validation
# ---------------------------------------------------------------------------


def parse_rational(value: Any, code: str) -> Fraction:
    try:
        return Fraction(str(value))
    except (ValueError, ZeroDivisionError) as exc:
        raise CertificateError(code, f"cannot parse exact rational {value!r}") from exc


def validate_manifest(
    manifest: Mapping[str, Any],
    base_dir: Path | None = None,
    *,
    allow_control_models: bool = False,
) -> dict[str, Any]:
    e565.enforce_source_firewall(manifest)
    require(manifest.get("schema") == SCHEMA, "SCHEMA", f"expected {SCHEMA}")

    response = manifest.get("reversible_response_automorphisms")
    require(isinstance(response, Mapping), "RESPONSE_TYPING", "reversible_response_automorphisms is missing")
    require(response.get("reversible") is True, "RESPONSE_TYPING", "response automorphisms must be typed reversible")
    require(response.get("defines_currents") is True, "RESPONSE_TYPING", "response automorphisms must be the declared current source")

    # The semantic artifact binds incidence, inverse-port response constraints,
    # the oriented frame, and refinement maps.  It does not select a Lie
    # bracket.  The matrix current construction is therefore a declared
    # algebraic fixture until an ordered-response/overlap-holonomy producer
    # supplies the missing source bridge.
    contract = manifest.get("response_declaration_contract")
    require(isinstance(contract, Mapping), "RESPONSE_TYPING", "response_declaration_contract is missing")
    require(
        contract.get("distinct_from_register_relabeling") is True,
        "REGISTER_RELABELING_CONFLATION",
        "the response fields must be typed distinct from register relabeling",
    )
    require(
        "measurement_artifact" not in contract,
        "RESPONSE_ARTIFACT",
        "measurement data enters only through the reviewed semantic artifact schema",
    )

    repairs = manifest.get("strict_descent_repairs")
    require(isinstance(repairs, Mapping), "REPAIR_TYPING", "strict_descent_repairs is missing")
    require(repairs.get("irreversible") is True, "REPAIR_TYPING", "strict-descent repairs must be typed irreversible")
    require(
        repairs.get("defines_currents") is False,
        "REPAIR_RESPONSE_CONFLATION",
        "irreversible strict-descent repairs may not be declared as current sources",
    )
    ledger = repairs.get("ledger")
    require(isinstance(ledger, list), "REPAIR_TYPING", "repair ledger must be a list")
    for row in ledger:
        require(isinstance(row, Mapping), "REPAIR_TYPING", "each repair record must be an object")
        require(row.get("reversible") is False, "REPAIR_TYPING", "repair records must be irreversible")
        require(
            row.get("defines_currents") is False,
            "REPAIR_RESPONSE_CONFLATION",
            f"repair record {row.get('repair_id')} is conflated with the current source",
        )

    artifact_ref = manifest.get("semantic_response_artifact")
    model = manifest.get("construction_model")
    if model is not None:
        require(
            model
            in ("charged_double_triplet", "abelian_record", "symmetric_record_control"),
            "RESPONSE_MODEL",
            "unknown control construction model",
        )
        if allow_control_models:
            artifact_ref = None
            response_status = "declared_branch_premise"
        else:
            require(
                model == "charged_double_triplet",
                "RESPONSE_MODEL",
                "the production algebraic fixture must be charged_double_triplet",
            )
            require(
                isinstance(artifact_ref, Mapping),
                "ARTIFACT_REFERENCE",
                "the production fixture must bind the response-constraint artifact",
            )
            require(
                contract.get("status")
                == "declared_current_fixture_with_source_bound_response_constraints",
                "RESPONSE_TYPING",
                "the production fixture must be typed as a declared current "
                "fixture with source-bound response constraints",
            )
            response_status = (
                "declared_current_fixture_with_source_bound_response_constraints"
            )
    else:
        require(
            allow_control_models,
            "CONSTRUCTION_MODEL_STRING",
            "a production manifest must name its algebraic current fixture; "
            "the response artifact does not select it",
        )
        model = "charged_double_triplet"
        artifact_ref = None
        response_status = "declared_branch_premise"

    scales_raw = manifest.get("response_band_scales")
    require(isinstance(scales_raw, Mapping), "RESPONSE_SCALES", "response_band_scales is missing")
    band_names = ("unit_band", "quintet_band", "frame_band", "kernel_band")
    require(set(scales_raw) == set(band_names), "RESPONSE_SCALES", f"scales must name exactly the bands {band_names}")
    scales = {name: parse_rational(scales_raw[name], "RESPONSE_SCALES") for name in band_names}

    axis_scales_raw = manifest.get("even_quintet_axis_scales", ["1"] * 6)
    require(
        isinstance(axis_scales_raw, list) and len(axis_scales_raw) == 6,
        "AXIS_SCALES",
        "even_quintet_axis_scales must list six exact rationals",
    )
    axis_scales = [parse_rational(x, "AXIS_SCALES") for x in axis_scales_raw]

    odd_axis_signs_raw = manifest.get("odd_axis_signs", [1] * 6)
    require(
        isinstance(odd_axis_signs_raw, list)
        and len(odd_axis_signs_raw) == 6
        and all(s in (1, -1) for s in odd_axis_signs_raw),
        "ODD_AXIS_SIGNS",
        "odd_axis_signs must list six values in {+1,-1}",
    )

    return {
        "model": model,
        "scales": scales,
        "axis_scales": axis_scales,
        "odd_axis_signs": [int(s) for s in odd_axis_signs_raw],
        "repair_ledger_rows": len(ledger),
        "response_status": response_status,
        "artifact_ref": artifact_ref,
    }


_CARRIER_CACHE: dict[str, tuple[Any, dict[str, Any], list[Any], list[Any]]] = {}


def load_carrier(manifest: Mapping[str, Any], base_dir: Path) -> tuple[Any, dict[str, Any], list[Any], Mapping[str, Any]]:
    path_raw = manifest.get("carrier_manifest_path")
    require(isinstance(path_raw, str), "CARRIER_REFERENCE", "carrier_manifest_path is missing")
    path = Path(path_raw)
    if not path.is_absolute():
        path = base_dir / path
    carrier_manifest = load_json(path)
    digest = sha256_json(carrier_manifest)
    declared = manifest.get("carrier_manifest_sha256")
    require(declared == digest, "CARRIER_HASH", "carrier manifest hash does not match the declared pin")
    if digest not in _CARRIER_CACHE:
        carrier = e565.validate_carrier(carrier_manifest)
        group_row, plus, minus = e565.group_certificate(carrier)
        _CARRIER_CACHE[digest] = (carrier, group_row, plus, minus)
    carrier, group_row, plus, _ = _CARRIER_CACHE[digest]
    return carrier, group_row, plus, carrier_manifest


# ---------------------------------------------------------------------------
# Semantic response artifact binding (issue #599)
# ---------------------------------------------------------------------------


def parse_f5_text(value: Any, code: str) -> F5:
    """Parse an exact Q(sqrt5) string: 'a', 'b*sqrt(5)', or 'a + b*sqrt(5)'."""

    require(isinstance(value, str) and value.strip(), code, f"expected an exact Q(sqrt5) string, got {value!r}")
    text = value.strip()
    try:
        if "sqrt(5)" not in text:
            return F5(Fraction(text))
        if "+" in text:
            left, right = text.split("+", 1)
            rational = Fraction(left.strip())
        else:
            left, right = "0", text
            rational = Fraction(0)
        radical = right.strip()
        require(radical.endswith("*sqrt(5)"), code, f"malformed radical part in {value!r}")
        coefficient = Fraction(radical[: -len("*sqrt(5)")])
        return F5(rational, coefficient)
    except (ValueError, ZeroDivisionError) as exc:
        raise CertificateError(code, f"cannot parse exact Q(sqrt5) value {value!r}") from exc


def artifact_self_hash(artifact: Mapping[str, Any]) -> str:
    body = {key: value for key, value in artifact.items() if key != "artifact_sha256"}
    return "sha256:" + sha256_json(body)


def load_semantic_artifact(
    artifact_ref: Mapping[str, Any],
    base_dir: Path,
    allow_inline_artifact_for_tests: bool,
) -> dict[str, Any]:
    if "value" in artifact_ref:
        require(
            allow_inline_artifact_for_tests,
            "ARTIFACT_REFERENCE",
            "inline artifacts are test-only; production binds a hash-pinned path",
        )
        artifact = artifact_ref["value"]
    else:
        path_raw = artifact_ref.get("path")
        require(isinstance(path_raw, str), "ARTIFACT_REFERENCE", "semantic artifact path is missing")
        path = Path(path_raw)
        if not path.is_absolute():
            path = base_dir / path
        artifact = load_json(path)
        declared = artifact_ref.get("artifact_sha256")
        require(
            isinstance(declared, str) and declared == artifact.get("artifact_sha256"),
            "ARTIFACT_HASH",
            "the manifest pin does not match the artifact self-hash field",
        )
    require(isinstance(artifact, Mapping), "ARTIFACT_SCHEMA", "the semantic artifact must be an object")
    require(
        artifact.get("schema") == ARTIFACT_SCHEMA,
        "ARTIFACT_SCHEMA",
        f"expected semantic artifact schema {ARTIFACT_SCHEMA}",
    )
    require(artifact.get("issue") == 599, "ARTIFACT_SCHEMA", "the semantic artifact must name issue 599")
    require(
        artifact.get("artifact_sha256") == artifact_self_hash(artifact),
        "ARTIFACT_HASH",
        "the artifact self-hash does not match its content",
    )
    return dict(artifact)


ARTIFACT_BAND_ORDER = ("unit_band", "quintet_band", "frame_band", "kernel_band")
ARTIFACT_SECTOR_DIMENSIONS = {
    "unit_band": 1,
    "quintet_band": 5,
    "frame_band": 3,
    "kernel_band": 3,
}
ARTIFACT_CHANNEL_VALUES = {
    "unit_band": F5(5),
    "quintet_band": F5(-1),
    "frame_band": F5(0, 1),
    "kernel_band": F5(0, -1),
}


def _adjacency_matrix(carrier: Any) -> RMat:
    return [
        [ONE if carrier.distances[i][j] == 1 else ZERO for j in range(12)]
        for i in range(12)
    ]


def _recompute_isotypic_channels(carrier: Any) -> dict[str, RMat]:
    """Exact spectral projectors of the carrier adjacency, recomputed here.

    This is the paper-side recomputation of the artifact's central claim: the
    source incidence presents exactly the 1 + 3 + 3' + 5 sector structure
    with Galois-paired triplet channels.
    """

    adjacency = _adjacency_matrix(carrier)
    identity: RMat = [[ONE if i == j else ZERO for j in range(12)] for i in range(12)]
    projectors: dict[str, RMat] = {}
    for band in ARTIFACT_BAND_ORDER:
        eigenvalue = ARTIFACT_CHANNEL_VALUES[band]
        product = identity
        for other_band in ARTIFACT_BAND_ORDER:
            if other_band == band:
                continue
            other = ARTIFACT_CHANNEL_VALUES[other_band]
            shifted = [
                [adjacency[i][j] - (other if i == j else ZERO) for j in range(12)]
                for i in range(12)
            ]
            product = rmul(product, shifted)
            scale = (eigenvalue - other).inv()
            product = [[scale * entry for entry in row] for row in product]
        projectors[band] = product

    residual = [[entry for entry in row] for row in identity]
    for band, projector in projectors.items():
        square = rmul(projector, projector)
        require(
            all(
                (square[i][j] - projector[i][j]).is_zero()
                for i in range(12)
                for j in range(12)
            ),
            "ARTIFACT_SECTORS",
            f"the recomputed {band} projector is not idempotent; the carrier "
            "does not present the charged sector structure",
        )
        image = rmul(_adjacency_matrix(carrier), projector)
        require(
            all(
                (
                    image[i][j]
                    - ARTIFACT_CHANNEL_VALUES[band] * projector[i][j]
                ).is_zero()
                for i in range(12)
                for j in range(12)
            ),
            "ARTIFACT_SECTORS",
            f"the recomputed {band} projector fails its eigenrelation",
        )
        trace = ZERO
        for i in range(12):
            trace = trace + projector[i][i]
        require(
            trace == F5(ARTIFACT_SECTOR_DIMENSIONS[band]),
            "ARTIFACT_SECTORS",
            f"the recomputed {band} sector dimension is not "
            f"{ARTIFACT_SECTOR_DIMENSIONS[band]}",
        )
        residual = [
            [residual[i][j] - projector[i][j] for j in range(12)] for i in range(12)
        ]
    require(
        all(residual[i][j].is_zero() for i in range(12) for j in range(12)),
        "ARTIFACT_SECTORS",
        "the recomputed sector projectors do not resolve the identity",
    )
    frame_projector = projectors["frame_band"]
    kernel_projector = projectors["kernel_band"]
    require(
        all(
            (frame_projector[i][j].conj() - kernel_projector[i][j]).is_zero()
            for i in range(12)
            for j in range(12)
        ),
        "ARTIFACT_SECTORS",
        "Galois conjugation does not pair the recomputed triplet sectors",
    )
    return projectors


def bind_semantic_artifact(
    artifact: Mapping[str, Any],
    carrier: Any,
    carrier_manifest: Mapping[str, Any],
    verts: Sequence[Vec3],
    matched: Sequence[tuple[int, ...]],
    proper_actions: Sequence[Sequence[int]],
    params: Mapping[str, Any],
) -> dict[str, Any]:
    """Verify the semantic artifact against paper-side recomputation.

    Returns the artifact-selected oriented frame assignment, the physical
    refinement port maps, and the binding report for the receipt. Every
    exact claim in the artifact is recomputed here before use.
    """

    binding = artifact.get("carrier_binding")
    require(isinstance(binding, Mapping), "ARTIFACT_CARRIER", "carrier_binding block is missing")
    require(
        binding.get("carrier_manifest_sha256") == sha256_json(carrier_manifest),
        "ARTIFACT_CARRIER",
        "the artifact does not bind the pinned carrier manifest",
    )
    require(
        list(binding.get("port_order", [])) == list(carrier.ports),
        "ARTIFACT_CARRIER",
        "the artifact port order does not match the carrier",
    )
    antipode_map = binding.get("antipode")
    require(isinstance(antipode_map, Mapping), "ARTIFACT_CARRIER", "antipode map missing")
    for position, port in enumerate(carrier.ports):
        require(
            antipode_map.get(port) == carrier.ports[carrier.antipode[position]],
            "ARTIFACT_CARRIER",
            "the artifact antipode map does not match the carrier",
        )

    basis = artifact.get("response_basis")
    require(isinstance(basis, Mapping), "ARTIFACT_SECTORS", "response_basis block is missing")
    require(
        list(basis.get("sector_order", [])) == list(ARTIFACT_BAND_ORDER),
        "ARTIFACT_SECTORS",
        "the artifact sector order does not name the four response bands",
    )
    dimensions = basis.get("sector_dimensions")
    require(
        isinstance(dimensions, Mapping)
        and {name: dimensions.get(name) for name in ARTIFACT_BAND_ORDER}
        == ARTIFACT_SECTOR_DIMENSIONS,
        "ARTIFACT_SECTORS",
        "the artifact sector dimensions are not (1, 5, 3, 3)",
    )
    channels = basis.get("adjacency_channel_values")
    require(isinstance(channels, Mapping), "ARTIFACT_SECTORS", "channel values missing")
    projectors = _recompute_isotypic_channels(carrier)
    for band in ARTIFACT_BAND_ORDER:
        parsed = parse_f5_text(channels.get(band), "ARTIFACT_SECTORS")
        require(
            parsed == ARTIFACT_CHANNEL_VALUES[band],
            "ARTIFACT_SECTORS",
            f"the artifact {band} channel value does not match the recomputed spectrum",
        )
    pairing = basis.get("galois_pairing")
    require(
        isinstance(pairing, Mapping)
        and pairing.get("frame_and_kernel_swapped_by_conjugation") is True
        and pairing.get("unit_and_quintet_galois_stable") is True,
        "ARTIFACT_SECTORS",
        "the artifact does not record the recomputed Galois pairing",
    )

    orientation = artifact.get("orientation_convention")
    require(isinstance(orientation, Mapping), "ARTIFACT_ORIENTATION", "orientation convention missing")
    require(
        "charge conjugation" in str(orientation.get("overall_u1_charge_sign", "")),
        "ARTIFACT_ORIENTATION",
        "the artifact must expose the common U(1) response sign as conventional",
    )
    require(
        ARTIFACT_CHANNEL_VALUES["frame_band"].is_positive()
        and not ARTIFACT_CHANNEL_VALUES["kernel_band"].is_positive(),
        "ARTIFACT_ORIENTATION",
        "orientation invariant failed",
    )

    frame_map = artifact.get("port_vertex_frame")
    require(isinstance(frame_map, Mapping), "ARTIFACT_FRAME", "port_vertex_frame is missing")
    vertex_index = {tuple(v): position for position, v in enumerate(verts)}
    psi: list[int] = []
    for position, port in enumerate(carrier.ports):
        coordinates = frame_map.get(port)
        require(
            isinstance(coordinates, list) and len(coordinates) == 3,
            "ARTIFACT_FRAME",
            f"port {port} lacks exact frame coordinates",
        )
        vector = tuple(parse_f5_text(text, "ARTIFACT_FRAME") for text in coordinates)
        slot = vertex_index.get(vector)
        require(
            slot is not None,
            "ARTIFACT_FRAME",
            f"port {port} is not assigned a standard icosahedron vertex",
        )
        psi.append(slot)
    require(
        tuple(psi) in set(matched),
        "ARTIFACT_FRAME",
        "the artifact frame assignment is not an oriented incidence realization",
    )

    structural = artifact.get("structural_audits")
    require(isinstance(structural, Mapping), "ARTIFACT_CHANNELS", "structural_audits missing")
    frame_audit = structural.get("frame_normalization")
    require(isinstance(frame_audit, Mapping), "ARTIFACT_CHANNELS", "frame normalization audit missing")
    gram = rzeros(3, 3)
    for position in range(12):
        vertex = verts[psi[position]]
        for i in range(3):
            for j in range(3):
                gram[i][j] = gram[i][j] + vertex[i] * vertex[j]
    tight = F5(10, 2)
    require(
        all(
            (gram[i][j] - (tight if i == j else ZERO)).is_zero()
            for i in range(3)
            for j in range(3)
        ),
        "ARTIFACT_CHANNELS",
        "the recomputed tight-frame constant is not (10 + 2*sqrt(5))",
    )
    require(
        parse_f5_text(frame_audit.get("unit_channel_constant"), "ARTIFACT_CHANNELS") == tight
        and frame_audit.get("response_sign_not_inferred_here") is True,
        "ARTIFACT_CHANNELS",
        "the artifact frame normalization does not match the recomputation",
    )
    rotation = structural.get("rotation_automorphisms")
    require(isinstance(rotation, Mapping), "ARTIFACT_CHANNELS", "rotation automorphism audit missing")
    require(
        rotation.get("incidence_automorphism_group_order") == 120
        and rotation.get("commuting_incidence_automorphisms") == len(matched)
        and rotation.get("response_sign_not_inferred_here") is True,
        "ARTIFACT_CHANNELS",
        "the artifact rotation audit does not match the recomputed group data",
    )

    source_response = artifact.get("source_response")
    require(isinstance(source_response, Mapping), "ARTIFACT_RESPONSE", "source_response missing")
    protocol = source_response.get("impulse_readback_protocol")
    require(
        isinstance(protocol, Mapping)
        and protocol.get("status") == "source_bound_operational_producer"
        and protocol.get("input") == "delta impulse at each unlabeled carrier port"
        and protocol.get("unique_solution_rank") == 4
        and protocol.get("unique_farthest_port_per_source") is True
        and protocol.get("nearer_shells_cancelled") is True
        and protocol.get("target_labels_used") is False
        and protocol.get("downstream_labels_used") is False,
        "ARTIFACT_RESPONSE",
        "the artifact does not expose the target-blind impulse/readback producer",
    )
    antipode = tuple(int(value) for value in carrier.antipode)
    require(
        source_response.get("operator") == "negative_graph_antipode_involution"
        and source_response.get("source")
        == "target_blind_maximal_distance_impulse_readback"
        and tuple(source_response.get("antipode_port_map", ())) == antipode
        and source_response.get("commutes_with_propagation_generator") is True
        and source_response.get("self_adjoint_unitary_involution") is True
        and source_response.get("impulse_readback_response_executed") is True
        and source_response.get("physical_perturb_readback_source_bound") is True,
        "ARTIFACT_RESPONSE",
        "the artifact does not bind the impulse/readback-derived negative antipode response",
    )
    identity_permutation = tuple(range(12))
    full_actions = {
        tuple(int(value) for value in action)
        for action in proper_actions
    }
    full_actions |= {
        tuple(antipode[int(action[index])] for index in range(12))
        for action in proper_actions
    }
    require(
        len(full_actions) == 120,
        "ARTIFACT_RESPONSE",
        "the proper actions and antipode do not generate 120 incidence automorphisms",
    )
    central_involutions = [
        action
        for action in full_actions
        if all(action[action[index]] == index for index in range(12))
        and all(
            action[other[index]] == other[action[index]]
            for other in full_actions
            for index in range(12)
        )
    ]
    require(
        set(central_involutions) == {identity_permutation, antipode}
        and source_response.get("unique_nonidentity_central_involution") is True
        and source_response.get("central_involution_count_including_identity") == 2,
        "ARTIFACT_RESPONSE",
        "the antipode is not the unique nonidentity central involution",
    )
    response = rzeros(12, 12)
    antipode_matrix = rzeros(12, 12)
    for index, partner in enumerate(antipode):
        response[index][partner] = -ONE
        antipode_matrix[index][partner] = ONE
    adjacency = _adjacency_matrix(carrier)
    adjacency_squared = rmul(adjacency, adjacency)
    adjacency_cubed = rmul(adjacency_squared, adjacency)
    powers = [
        [[ONE if i == j else ZERO for j in range(12)] for i in range(12)],
        adjacency,
        adjacency_squared,
        adjacency_cubed,
    ]
    diameter = max(
        carrier.distances[i][j] for i in range(12) for j in range(12)
    )
    require(diameter == 3, "ARTIFACT_RESPONSE", "carrier diameter is not three")
    target = [
        [ONE if carrier.distances[i][j] == diameter else ZERO for j in range(12)]
        for i in range(12)
    ]
    require(
        all(sum(target[i][j] == ONE for j in range(12)) == 1 for i in range(12)),
        "ARTIFACT_RESPONSE",
        "an impulse source lacks a unique maximal-distance readback port",
    )
    augmented = [
        [powers[k][i][j] for k in range(4)] + [target[i][j]]
        for i in range(12)
        for j in range(12)
    ]
    reduced, pivots = rref(augmented)
    require(
        pivots == [0, 1, 2, 3],
        "ARTIFACT_RESPONSE",
        "the common maximal-distance impulse filter is not uniquely solvable",
    )
    filter_coefficients = [reduced[index][4] for index in range(4)]
    require(
        filter_coefficients
        == [ONE, F5(Fraction(-1, 2)), F5(Fraction(-2, 5)), F5(Fraction(1, 10))]
        and [
            parse_f5_text(value, "ARTIFACT_RESPONSE")
            for value in protocol.get("homogeneous_filter_coefficients", [])
        ]
        == filter_coefficients,
        "ARTIFACT_RESPONSE",
        "the artifact impulse filter does not match the independent exact solve",
    )
    solved_target = rzeros(12, 12)
    for coefficient, power in zip(filter_coefficients, powers, strict=True):
        solved_target = [
            [
                solved_target[i][j] + coefficient * power[i][j]
                for j in range(12)
            ]
            for i in range(12)
        ]
    require(
        solved_target == target == antipode_matrix,
        "ARTIFACT_RESPONSE",
        "the solved impulse filter does not isolate the carrier antipode",
    )
    tenth = F5(Fraction(1, 10))
    polynomial_antipode = [
        [
            tenth
            * (
                adjacency_cubed[i][j]
                - F5(4) * adjacency_squared[i][j]
                - F5(5) * adjacency[i][j]
                + (F5(10) if i == j else ZERO)
            )
            for j in range(12)
        ]
        for i in range(12)
    ]
    require(
        polynomial_antipode == antipode_matrix
        and source_response.get("antipode_polynomial_identity")
        == "J = (A^3 - 4*A^2 - 5*A + 10*I)/10",
        "ARTIFACT_RESPONSE",
        "the exact adjacency-polynomial identity for the antipode failed",
    )
    response_signs: dict[str, int] = {}
    for band, projector in projectors.items():
        image = rmul(response, projector)
        if image == projector:
            response_signs[band] = 1
        elif image == [[-entry for entry in row] for row in projector]:
            response_signs[band] = -1
        else:
            raise CertificateError(
                "ARTIFACT_RESPONSE",
                f"the {band} sector is not an eigenspace of the antipode response",
            )
    artifact_signs = source_response.get("sector_eigenvalues")
    require(
        isinstance(artifact_signs, Mapping)
        and {
            band: artifact_signs.get(band)
            for band in ARTIFACT_BAND_ORDER
        }
        == response_signs,
        "ARTIFACT_RESPONSE",
        "the artifact response signs do not match the exact antipode eigenspaces",
    )

    derived = artifact.get("derived")
    require(isinstance(derived, Mapping), "ARTIFACT_SCALES", "derived block missing")
    lift_status = derived.get("current_lift_status")
    require(
        "construction" not in derived
        and "construction_provenance" not in derived
        and isinstance(lift_status, Mapping)
        and lift_status.get("source_selected") is False
        and lift_status.get("commutator_reconstructed_from_ordered_response")
        is False
        and lift_status.get("overlap_holonomy_internality_certified") is False
        and lift_status.get("charged_double_triplet_forced") is False,
        "ARTIFACT_CURRENT_BOUNDARY",
        "the response artifact must fail closed at the unconstructed current lift",
    )
    derived_scales = derived.get("response_band_scales")
    require(isinstance(derived_scales, Mapping), "ARTIFACT_SCALES", "derived scales missing")
    for band in ARTIFACT_BAND_ORDER:
        value = parse_rational(derived_scales.get(band), "ARTIFACT_SCALES")
        require(value != 0, "ARTIFACT_SCALES", f"the derived {band} scale is zero")
        sign = 1 if value > 0 else -1
        require(
            sign == response_signs[band],
            "ARTIFACT_ORIENTATION",
            f"the derived {band} sign does not match the impulse/readback-derived response sign",
        )
        manifest_scale = params["scales"][band]
        require(
            manifest_scale == value,
            "COEFFICIENT_MISMATCH",
            f"the manifest {band} coefficient does not equal the artifact-derived value",
        )

    refinement = artifact.get("physical_refinement_maps")
    require(isinstance(refinement, Mapping), "ARTIFACT_REFINEMENT", "physical refinement maps missing")
    rows = refinement.get("port_persistence_maps")
    require(isinstance(rows, list) and rows, "ARTIFACT_REFINEMENT", "no physical refinement maps recorded")
    port_maps: list[dict[str, Any]] = []
    for row in rows:
        require(isinstance(row, Mapping), "ARTIFACT_REFINEMENT", "malformed refinement row")
        port_map = row.get("port_map")
        require(
            isinstance(port_map, list) and sorted(port_map) == list(range(12)),
            "ARTIFACT_REFINEMENT",
            "a physical refinement map is not a port permutation",
        )
        for i in range(12):
            for j in range(12):
                require(
                    carrier.distances[i][j]
                    == carrier.distances[port_map[i]][port_map[j]],
                    "ARTIFACT_REFINEMENT",
                    "a physical refinement map does not preserve the incidence",
                )
        body = {key: value for key, value in row.items() if key != "map_hash"}
        require(
            row.get("map_hash") == "sha256:" + sha256_json(body),
            "ARTIFACT_REFINEMENT",
            "a physical refinement map hash does not recompute",
        )
        port_maps.append(
            {
                "source_level": row.get("source_level"),
                "target_level": row.get("target_level"),
                "port_map": [int(v) for v in port_map],
                "origin": row.get("origin"),
            }
        )

    provenance = artifact.get("provenance")
    require(isinstance(provenance, Mapping), "ARTIFACT_SCHEMA", "provenance block missing")
    runtime = provenance.get("runtime_binding")
    require(
        isinstance(runtime, Mapping)
        and runtime.get("spectrum_multiplicities") == [1, 3, 3, 5]
        and runtime.get("equivariance_receipt") is True
        and runtime.get("charged_response_operator_receipt") is True,
        "ARTIFACT_SCHEMA",
        "the artifact runtime binding does not certify the propagated dynamics",
    )

    return {
        "psi": tuple(psi),
        "port_maps": port_maps,
        "report": {
            "artifact_sha256": artifact.get("artifact_sha256"),
            "producer": provenance.get("producer"),
            "dynamics_sha256": runtime.get("dynamics_sha256"),
            "carrier_manifest_sha256": binding.get("carrier_manifest_sha256"),
            "sector_structure_recomputed": True,
            "galois_pairing_recomputed": True,
            "tight_frame_constant_recomputed": "10 + 2*sqrt(5)",
            "frame_assignment_source": "artifact port_vertex_frame",
            "current_lift_source_selected": lift_status.get("source_selected"),
            "commutator_reconstructed_from_ordered_response": lift_status.get(
                "commutator_reconstructed_from_ordered_response"
            ),
            "overlap_holonomy_internality_certified": lift_status.get(
                "overlap_holonomy_internality_certified"
            ),
            "charged_double_triplet_forced": lift_status.get(
                "charged_double_triplet_forced"
            ),
            "response_operator": source_response.get("operator"),
            "response_source": source_response.get("source"),
            "impulse_readback_status": protocol.get("status"),
            "impulse_readback_filter_coefficients_recomputed": [
                value.text() for value in filter_coefficients
            ],
            "physical_perturb_readback_source_bound": True,
            "antipode_polynomial_identity_recomputed": True,
            "unique_nonidentity_central_involution_recomputed": True,
            "response_sector_eigenvalues_recomputed": response_signs,
            "overall_u1_charge_sign": orientation.get("overall_u1_charge_sign"),
            "derived_response_band_scales": {
                band: str(derived_scales.get(band)) for band in ARTIFACT_BAND_ORDER
            },
            "normalization_rule": derived.get("normalization_rule"),
            "physical_refinement_map_count": len(port_maps),
        },
    }


# ---------------------------------------------------------------------------
# Frame realization: ports -> exact vertex vectors, unique up to A5
# ---------------------------------------------------------------------------


def orientation_matched_assignments(carrier: Any, verts: Sequence[Vec3]) -> list[tuple[int, ...]]:
    coord_distances = [
        [vertex_distance(u, v) for v in verts] for u in verts
    ]
    isometries = e565.enumerate_distance_isometries(carrier.distances, coord_distances)
    require(len(isometries) == 120, "FRAME_ISOMETRY_COUNT", f"expected 120 distance isometries, got {len(isometries)}")
    matched: list[tuple[int, ...]] = []
    for psi in isometries:
        signs = set()
        for a, b, c in carrier.faces:
            m = [list(verts[psi[a]]), list(verts[psi[b]]), list(verts[psi[c]])]
            d = det3(m)
            require(not d.is_zero(), "FRAME_FACE_DEGENERATE", "mapped face is degenerate")
            signs.add(d.is_positive())
            if len(signs) > 1:
                break
        if signs == {True}:
            matched.append(psi)
    require(len(matched) == 60, "FRAME_ORIENTATION_COUNT", f"expected 60 orientation-matched realizations, got {len(matched)}")
    return matched


class FrameRealization:
    """An oriented exact realization of the port frame and its band data."""

    def __init__(self, carrier: Any, psi: tuple[int, ...], verts: Sequence[Vec3]) -> None:
        self.carrier = carrier
        self.psi = psi
        self.vertex: list[Vec3] = [verts[psi[p]] for p in range(12)]
        for p in range(12):
            anti = carrier.antipode[p]
            require(
                all((self.vertex[p][k] + self.vertex[anti][k]).is_zero() for k in range(3)),
                "FRAME_ANTIPODE",
                "realized antipodes are not opposite vectors",
            )
        self.axis_reps: list[int] = [p for p in range(12) if p < carrier.antipode[p]]
        require(len(self.axis_reps) == 6, "FRAME_AXES", "expected six axis representatives")
        self.axis_vectors: list[Vec3] = [self.vertex[p] for p in self.axis_reps]
        self._base_triple: tuple[int, int, int] | None = None
        self._rotation_cache: dict[tuple[int, ...], RMat] = {}

    def even_odd(self, field: Sequence[Fraction]) -> tuple[list[F5], list[F5]]:
        """Split a port field into even and odd axis coordinates."""

        even: list[F5] = []
        odd: list[F5] = []
        half = F5(Fraction(1, 2))
        for rep in self.axis_reps:
            anti = self.carrier.antipode[rep]
            even.append(half * (F5(field[rep]) + F5(field[anti])))
            odd.append(half * (F5(field[rep]) - F5(field[anti])))
        return even, odd

    def frame_map(self, odd: Sequence[F5]) -> Vec3:
        """U d = sum_i d_i u_i over the axis vectors."""

        total = [ZERO, ZERO, ZERO]
        for coef, u in zip(odd, self.axis_vectors, strict=True):
            for k in range(3):
                total[k] = total[k] + coef * u[k]
        return (total[0], total[1], total[2])

    def galois_frame_map(self, odd: Sequence[F5]) -> Vec3:
        """sigma(U) d, the Galois-twisted frame map (kills the frame band)."""

        total = [ZERO, ZERO, ZERO]
        for coef, u in zip(odd, self.axis_vectors, strict=True):
            for k in range(3):
                total[k] = total[k] + coef * u[k].conj()
        return (total[0], total[1], total[2])

    def rotation_of(self, permutation: Sequence[int]) -> RMat:
        """The unique matrix R with R v_p = v_{g(p)}; verified orthogonal, det 1."""

        key = tuple(permutation)
        cached = self._rotation_cache.get(key)
        if cached is not None:
            return cached
        if self._base_triple is None:
            base = None
            for i in range(12):
                for j in range(i + 1, 12):
                    for k in range(j + 1, 12):
                        m = [list(self.vertex[i]), list(self.vertex[j]), list(self.vertex[k])]
                        if not det3(rtranspose(m)).is_zero():
                            base = (i, j, k)
                            break
                    if base:
                        break
                if base:
                    break
            require(base is not None, "FRAME_SPAN", "vertex vectors do not span")
            self._base_triple = base
        i, j, k = self._base_triple
        # Solve R [v_i v_j v_k] = [v_gi v_gj v_gk] exactly through the
        # transposed system [S^T | T^T] -> [I | R^T].
        source_rows = [list(self.vertex[i]), list(self.vertex[j]), list(self.vertex[k])]
        target_rows = [
            list(self.vertex[permutation[i]]),
            list(self.vertex[permutation[j]]),
            list(self.vertex[permutation[k]]),
        ]
        augmented = [source_rows[r][:] + target_rows[r][:] for r in range(3)]
        reduced, pivots = rref(augmented)
        require(pivots == [0, 1, 2], "FRAME_SOLVE", "vertex triple is not invertible")
        r_transpose = [[reduced[r][3 + c] for c in range(3)] for r in range(3)]
        rotation = rtranspose(r_transpose)
        for p in range(12):
            image = apply3(rotation, self.vertex[p])
            expected = self.vertex[permutation[p]]
            require(
                all((image[t] - expected[t]).is_zero() for t in range(3)),
                "FRAME_ROTATION",
                "solved rotation does not transport every vertex",
            )
        product = rmul(rtranspose(rotation), rotation)
        identity = [[ONE if a == b else ZERO for b in range(3)] for a in range(3)]
        require(product == identity, "FRAME_ROTATION_ORTHOGONAL", "implementer is not orthogonal")
        require(det3(rotation) == ONE, "FRAME_ROTATION_PROPER", "implementer is not proper")
        self._rotation_cache[key] = rotation
        return rotation


# ---------------------------------------------------------------------------
# Response models: port-to-generator maps and implementers
# ---------------------------------------------------------------------------


BASIS_FIELDS: list[list[Fraction]] = [
    [Fraction(1) if q == p else Fraction(0) for q in range(12)] for p in range(12)
]


class ChargedDoubleTripletModel:
    """K(f) = (frame-band + i(unit+quintet), sign * kernel-band) on C^3 (+) C^3."""

    def __init__(self, frame: FrameRealization, params: Mapping[str, Any]) -> None:
        self.frame = frame
        scales = params["scales"]
        self.lam_unit = F5(scales["unit_band"])
        self.lam_quintet = F5(scales["quintet_band"])
        self.lam_frame = F5(scales["frame_band"])
        self.lam_kernel = F5(scales["kernel_band"])
        self.axis_scales = [F5(x) for x in params["axis_scales"]]
        self.odd_axis_signs = [F5(s) for s in params["odd_axis_signs"]]
        self.blocks = (3, 3)

    def generator(self, field: Sequence[Fraction]) -> tuple[CMat, CMat]:
        even, odd = self.frame.even_odd(field)
        odd = [s * d for s, d in zip(self.odd_axis_signs, odd, strict=True)]
        sixth = F5(Fraction(1, 6))
        mean = (even[0] + even[1] + even[2] + even[3] + even[4] + even[5]) * sixth
        centered = [b - mean for b in even]

        sym = rzeros(3, 3)
        for w, b0, u in zip(self.axis_scales, centered, self.frame.axis_vectors, strict=True):
            block = outer(u, u)
            coef = self.lam_quintet * w * b0
            for i in range(3):
                for j in range(3):
                    sym[i][j] = sym[i][j] + coef * block[i][j]
        central = self.lam_unit * mean
        for i in range(3):
            sym[i][i] = sym[i][i] + central

        skew_e = hat(self.frame.frame_map(odd))
        real_e = [[self.lam_frame * skew_e[i][j] for j in range(3)] for i in range(3)]
        block_e = to_cmat(real_e, sym)

        skew_w = hat(self.frame.galois_frame_map(odd))
        real_w = [[self.lam_kernel * skew_w[i][j] for j in range(3)] for i in range(3)]
        block_w = to_cmat(real_w)
        return block_e, block_w

    def implementer(self, permutation: Sequence[int]) -> tuple[CMat, CMat]:
        rotation = self.frame.rotation_of(permutation)
        conjugate = [[rotation[i][j].conj() for j in range(3)] for i in range(3)]
        return to_cmat(rotation), to_cmat(conjugate)


class AbelianRecordModel:
    """Coefficient-record control: K(f) = i*diag(f) on C^12; records commute."""

    def __init__(self, frame: FrameRealization, params: Mapping[str, Any]) -> None:
        self.frame = frame
        self.blocks = (12,)

    def generator(self, field: Sequence[Fraction]) -> tuple[CMat, ...]:
        block = czeros(12)
        for p in range(12):
            block[p][p] = C5(ZERO, F5(field[p]))
        return (block,)

    def implementer(self, permutation: Sequence[int]) -> tuple[CMat, ...]:
        block = czeros(12)
        for p in range(12):
            block[permutation[p]][p] = C5(ONE, ZERO)
        return (block,)


class SymmetricRecordControl(ChargedDoubleTripletModel):
    """Control that drops the i on the even response, breaking skew-adjointness."""

    def generator(self, field: Sequence[Fraction]) -> tuple[CMat, CMat]:
        block_e, block_w = super().generator(field)
        broken = [[C5(entry.re + entry.im, ZERO) for entry in row] for row in block_e]
        return broken, block_w


MODELS: dict[str, type] = {
    "charged_double_triplet": ChargedDoubleTripletModel,
    "abelian_record": AbelianRecordModel,
    "symmetric_record_control": SymmetricRecordControl,
}


def flatten(blocks: Sequence[CMat]) -> list[F5]:
    out: list[F5] = []
    for block in blocks:
        for row in block:
            for entry in row:
                out.append(entry.re)
                out.append(entry.im)
    return out


# ---------------------------------------------------------------------------
# Certificate checks
# ---------------------------------------------------------------------------


def check_skew_adjoint(generators: Sequence[Sequence[CMat]]) -> None:
    for blocks in generators:
        for block in blocks:
            n = len(block)
            for i in range(n):
                for j in range(n):
                    require(
                        (block[j][i].conj() + block[i][j]).is_zero(),
                        "SKEW_ADJOINTNESS_BROKEN",
                        "a response generator is not skew-adjoint",
                    )


def check_covariance(
    model: Any,
    plus: Sequence[Sequence[int]],
) -> dict[str, Any]:
    implementers = {g: model.implementer(g) for g in plus}
    checked = 0
    for g, pi in implementers.items():
        pi_dagger = tuple(cdagger(block) for block in pi)
        ginv = inverse(tuple(g))
        for field in BASIS_FIELDS:
            # (g . f)(p) = f(g^{-1}(p)).
            moved = [field[ginv[p]] for p in range(12)]
            left = model.generator(moved)
            right_source = model.generator(field)
            for block_index in range(len(pi)):
                conjugated = cmul(cmul(pi[block_index], right_source[block_index]), pi_dagger[block_index])
                require(
                    c_is_zero(csub(left[block_index], conjugated)),
                    "COVARIANCE_BROKEN",
                    "K(g.f) != Pi(g) K(f) Pi(g)* for a group element and basis field",
                )
            checked += 1
    return {"pairs_checked": checked, "implementers": implementers}


def check_homomorphism(
    implementers: Mapping[Sequence[int], Sequence[CMat]],
    plus: Sequence[Sequence[int]],
) -> int:
    table = {tuple(g): pi for g, pi in implementers.items()}
    checked = 0
    for g in plus:
        for h in plus:
            gh = compose(tuple(g), tuple(h))
            product = tuple(
                cmul(bg, bh) for bg, bh in zip(table[tuple(g)], table[tuple(h)], strict=True)
            )
            expected = table[gh]
            for a, b in zip(product, expected, strict=True):
                require(
                    c_is_zero(csub(a, b)),
                    "IMPLEMENTER_HOMOMORPHISM",
                    "implementers do not compose as the group",
                )
            checked += 1
    distinct = {
        tuple(entry.re for block in pi for row in block for entry in row)
        + tuple(entry.im for block in pi for row in block for entry in row)
        for pi in implementers.values()
    }
    require(len(distinct) == len(implementers), "IMPLEMENTER_FAITHFUL", "implementers are not faithful")
    return checked


def solve_in_span(basis_flat: Sequence[Sequence[F5]], target: Sequence[F5]) -> list[F5]:
    """Solve sum_k x_k basis[k] = target exactly; fail closed if unsolvable."""

    rows = len(basis_flat)
    stacked = [list(basis_flat[k]) for k in range(rows)]
    reduced, pivots = rref(stacked)
    require(len(pivots) == rows, "BASIS_DEPENDENT", "generator basis is linearly dependent")
    # Solve via the pivot columns: restrict to a square invertible system.
    square = [[basis_flat[k][c] for k in range(rows)] for c in pivots]
    rhs = [target[c] for c in pivots]
    augmented = [square[r][:] + [rhs[r]] for r in range(rows)]
    solved, spivots = rref(augmented)
    require(spivots == list(range(rows)), "CLOSURE_SOLVE", "pivot system is singular")
    x = [solved[r][rows] for r in range(rows)]
    # Verify on all coordinates, not only pivots.
    for c in range(len(target)):
        acc = ZERO
        for k in range(rows):
            acc = acc + x[k] * basis_flat[k][c]
        require(
            (acc - target[c]).is_zero(),
            "COMMUTATOR_NOT_CLOSED",
            "a commutator leaves the span of the current generators",
        )
    return x


def rotation_normal_form(rotation: RMat) -> dict[str, Any]:
    """Exact Rodrigues normal form: proves rotation = exp(theta * hat(axis))."""

    identity = [[ONE if a == b else ZERO for b in range(3)] for a in range(3)]
    difference = [[rotation[i][j] - identity[i][j] for j in range(3)] for i in range(3)]
    if all(entry.is_zero() for row in difference for entry in row):
        return {"identity": True, "cosine": "1", "axis": ["0", "0", "0"]}
    kernel = nullspace(difference)
    require(len(kernel) == 1, "INNERNESS_AXIS", "rotation axis is not one-dimensional")
    axis = (kernel[0][0], kernel[0][1], kernel[0][2])
    eta = dot(axis, axis)
    require(not eta.is_zero(), "INNERNESS_AXIS", "axis vector is null")
    trace = rotation[0][0] + rotation[1][1] + rotation[2][2]
    cosine = (trace - ONE) / F5(2)
    hat_axis = hat(axis)
    hat_sq = rmul(hat_axis, hat_axis)
    coef = (ONE - cosine) / eta
    residual = [
        [difference[i][j] - coef * hat_sq[i][j] for j in range(3)]
        for i in range(3)
    ]
    # residual must equal t * hat(axis) with t^2 * eta = 1 - cosine^2.
    t_value: F5 | None = None
    for i in range(3):
        for j in range(3):
            h = hat_axis[i][j]
            r = residual[i][j]
            if h.is_zero():
                require(r.is_zero(), "INNERNESS_RODRIGUES", "residual is not proportional to the axis generator")
            else:
                candidate = r / h
                if t_value is None:
                    t_value = candidate
                else:
                    require(
                        (candidate - t_value).is_zero(),
                        "INNERNESS_RODRIGUES",
                        "residual proportionality is inconsistent",
                    )
    require(t_value is not None, "INNERNESS_RODRIGUES", "rotation residual vanished unexpectedly")
    require(
        (t_value * t_value * eta - (ONE - cosine * cosine)).is_zero(),
        "INNERNESS_RODRIGUES",
        "sine consistency t^2 |n|^2 = 1 - cos^2 fails",
    )
    return {
        "identity": False,
        "cosine": cosine.text(),
        "axis": [component.text() for component in axis],
    }


def band_projectors(frame: FrameRealization) -> dict[str, RMat]:
    """Exact port-space projectors onto the four isotypic response bands."""

    carrier = frame.carrier
    n = 12
    antipode_perm = rzeros(n, n)
    for p in range(n):
        antipode_perm[p][carrier.antipode[p]] = ONE
    identity = [[ONE if a == b else ZERO for b in range(n)] for a in range(n)]
    half = F5(Fraction(1, 2))
    even = [[(identity[i][j] + antipode_perm[i][j]) * half for j in range(n)] for i in range(n)]
    odd = [[(identity[i][j] - antipode_perm[i][j]) * half for j in range(n)] for i in range(n)]
    twelfth = F5(Fraction(1, 12))
    unit = [[twelfth for _ in range(n)] for _ in range(n)]
    quintet = [[even[i][j] - unit[i][j] for j in range(n)] for i in range(n)]

    # Frame-band projector: V^T V / (4 * vertex_norm) in port coordinates,
    # where V has columns v_p (signed vertex vectors) and V V^T = 4N * I_3.
    scale = (F5(4) * VERTEX_NORM).inv()
    gram = [[dot(frame.vertex[i], frame.vertex[j]) * scale for j in range(n)] for i in range(n)]
    frame_band = rmul(odd, gram)
    kernel_band = [[odd[i][j] - frame_band[i][j] for j in range(n)] for i in range(n)]

    for name, projector in (("unit", unit), ("quintet", quintet), ("frame", frame_band), ("kernel", kernel_band)):
        square = rmul(projector, projector)
        require(
            all((square[i][j] - projector[i][j]).is_zero() for i in range(n) for j in range(n)),
            "BAND_PROJECTOR",
            f"{name} band projector is not idempotent",
        )
    return {"unit_band": unit, "quintet_band": quintet, "frame_band": frame_band, "kernel_band": kernel_band}


def certificate_payload(
    manifest: Mapping[str, Any],
    base_dir: Path | None = None,
    *,
    allow_control_models: bool = False,
    allow_inline_artifact_for_tests: bool = False,
) -> dict[str, Any]:
    base = base_dir or MODULE_DIR
    params = validate_manifest(manifest, base, allow_control_models=allow_control_models)
    carrier, group_row, plus, carrier_manifest = load_carrier(manifest, base)

    verts = standard_vertices()
    matched = orientation_matched_assignments(carrier, verts)
    artifact_binding: dict[str, Any] | None = None
    if params["artifact_ref"] is not None:
        artifact = load_semantic_artifact(
            params["artifact_ref"], base, allow_inline_artifact_for_tests
        )
        artifact_binding = bind_semantic_artifact(
            artifact, carrier, carrier_manifest, verts, matched, plus, params
        )
        psi = artifact_binding["psi"]
    else:
        psi = matched[0]
    frame = FrameRealization(carrier, psi, verts)

    model_cls = MODELS[params["model"]]
    model = model_cls(frame, params)

    generators = [model.generator(field) for field in BASIS_FIELDS]
    check_skew_adjoint(generators)

    covariance = check_covariance(model, plus)
    homomorphism_pairs = check_homomorphism(covariance["implementers"], plus)

    flat = [flatten(blocks) for blocks in generators]
    image_rank = rank([row[:] for row in flat])
    require(
        image_rank == 12,
        "IMAGE_RANK_DEFICIENT",
        f"the port-to-generator map must be injective with twelve-dimensional image, got rank {image_rank}",
    )

    # Exact block identification: the even-block projections span the full
    # nine-dimensional skew-Hermitian algebra u(3), the kernel-block
    # projections are real and span the three-dimensional so(3).
    block_dimensions_verified = False
    if isinstance(model, ChargedDoubleTripletModel):
        for blocks in generators:
            for row in blocks[1]:
                for entry in row:
                    require(
                        entry.im.is_zero(),
                        "KERNEL_BLOCK_NOT_REAL",
                        "the kernel block must be a real skew matrix",
                    )
        even_rank = rank([row[:18] for row in flat])
        kernel_rank = rank([row[18:] for row in flat])
        require(
            even_rank == 9 and kernel_rank == 3,
            "BLOCK_DIMENSIONS",
            f"expected block dimensions (9, 3), got ({even_rank}, {kernel_rank})",
        )
        block_dimensions_verified = True

    structure_constants: dict[str, list[str]] = {}
    commutator_flats: list[list[F5]] = []
    for i in range(12):
        for j in range(i + 1, 12):
            bracket = tuple(
                commutator(generators[i][b], generators[j][b]) for b in range(len(generators[i]))
            )
            target = flatten(bracket)
            coeffs = solve_in_span(flat, target)
            commutator_flats.append(target)
            key = f"[{i},{j}]"
            structure_constants[key] = [c.text() for c in coeffs]

    derived_dimension = rank([row[:] for row in commutator_flats]) if commutator_flats else 0
    require(
        derived_dimension == 11,
        "CENTER_NOT_ONE_DIMENSIONAL",
        f"the conditional current construction needs derived dimension eleven; records that commute give {derived_dimension}",
    )

    # Exact type identification of the derived algebra: commutators are
    # traceless, so an eight-dimensional even derived block inside u(3) is
    # exactly su(3), and a three-dimensional kernel derived block is so(3).
    derived_even_rank = rank([row[:18] for row in commutator_flats])
    derived_kernel_rank = rank([row[18:] for row in commutator_flats])
    require(
        derived_even_rank == 8 and derived_kernel_rank == 3,
        "DERIVED_TYPE",
        f"expected derived block dimensions (8, 3), got ({derived_even_rank}, {derived_kernel_rank})",
    )

    # Center: x with [K(x), K_j] = 0 for all j, from the exact commutators.
    center_rows: list[list[F5]] = []
    bracket_table: dict[tuple[int, int], list[F5]] = {}
    index = 0
    for i in range(12):
        for j in range(i + 1, 12):
            bracket_table[(i, j)] = commutator_flats[index]
            index += 1
    ambient = len(flat[0])
    for j in range(12):
        for c in range(ambient):
            row = [ZERO] * 12
            for i in range(12):
                if i == j:
                    continue
                if i < j:
                    row[i] = bracket_table[(i, j)][c]
                else:
                    row[i] = -bracket_table[(j, i)][c]
            center_rows.append(row)
    center_basis = nullspace(center_rows)
    require(
        len(center_basis) == 1,
        "CENTER_NOT_ONE_DIMENSIONAL",
        f"expected a one-dimensional central u(1), got dimension {len(center_basis)}",
    )
    constant_direction = [ONE] * 12
    center_vector = center_basis[0]
    pivot = next(k for k in range(12) if not center_vector[k].is_zero())
    normalized_center = [entry / center_vector[pivot] for entry in center_vector]
    require(
        all((normalized_center[k] - ONE).is_zero() for k in range(12)),
        "CENTER_NOT_CONSTANT_LINE",
        "the central u(1) is not the constant even port line",
    )
    adjoint_rank = 12 - len(center_basis)
    require(adjoint_rank == 11, "ADJOINT_RANK", "adjoint rank must be eleven")

    # Central charge on the response space: K(constant field).
    constant_generator = model.generator([Fraction(1)] * 12)
    central_nonzero = any(not c_is_zero(block) for block in constant_generator)
    require(central_nonzero, "CHARGE_DEAD", "the central generator acts as zero; the response space is not charged")

    # Hilbert-Schmidt pullback: B_pq = -Re tr(K_p K_q), band decomposition.
    def hs_pairing(x: Sequence[CMat], y: Sequence[CMat]) -> F5:
        total = ZERO
        for bx, by in zip(x, y, strict=True):
            total = total + ctrace(cmul(bx, by)).re
        return -total

    hs_gram = [[hs_pairing(generators[p], generators[q]) for q in range(12)] for p in range(12)]
    for p in range(12):
        for q in range(12):
            require((hs_gram[p][q] - hs_gram[q][p]).is_zero(), "HS_SYMMETRY", "pullback form is not symmetric")
    # Positive definiteness by exact pivots of symmetric elimination.
    work = [row[:] for row in hs_gram]
    pivot_texts: list[str] = []
    for step in range(12):
        pivot_value = work[step][step]
        require(pivot_value.is_positive(), "COMPACT_TYPE", "pullback form is not positive definite")
        pivot_texts.append(pivot_value.text())
        for i in range(step + 1, 12):
            factor = work[i][step] / pivot_value
            for j in range(step, 12):
                work[i][j] = work[i][j] - factor * work[step][j]

    projectors = band_projectors(frame)
    band_coefficients: dict[str, str] = {}
    reconstruction = rzeros(12, 12)
    for name, projector in projectors.items():
        trace_bp = ZERO
        trace_p = ZERO
        product = rmul(hs_gram, projector)
        for d in range(12):
            trace_bp = trace_bp + product[d][d]
            trace_p = trace_p + projector[d][d]
        coefficient = trace_bp / trace_p
        band_coefficients[name] = coefficient.text()
        for i in range(12):
            for j in range(12):
                reconstruction[i][j] = reconstruction[i][j] + coefficient * projector[i][j]
    require(
        all((reconstruction[i][j] - hs_gram[i][j]).is_zero() for i in range(12) for j in range(12)),
        "BAND_DECOMPOSITION",
        "the pullback form is not band-scalar",
    )
    for g in plus:
        for p in range(12):
            for q in range(12):
                require(
                    (hs_gram[g[p]][g[q]] - hs_gram[p][q]).is_zero(),
                    "HS_INVARIANCE",
                    "pullback form is not A5-invariant",
                )

    # Assignment independence: an alternative orientation-matched realization
    # must reproduce the same exact band coefficients.
    alternative_frame = FrameRealization(carrier, matched[1], verts)
    alternative_model = model_cls(alternative_frame, params)
    alternative_generators = [alternative_model.generator(field) for field in BASIS_FIELDS]
    alternative_gram = [
        [hs_pairing(alternative_generators[p], alternative_generators[q]) for q in range(12)]
        for p in range(12)
    ]
    alternative_projectors = band_projectors(alternative_frame)
    for name, projector in alternative_projectors.items():
        trace_bp = ZERO
        trace_p = ZERO
        product = rmul(alternative_gram, projector)
        for d in range(12):
            trace_bp = trace_bp + product[d][d]
            trace_p = trace_p + projector[d][d]
        require(
            ((trace_bp / trace_p).text()) == band_coefficients[name],
            "ASSIGNMENT_INDEPENDENCE",
            f"band coefficient {name} depends on the frame realization choice",
        )

    # Innerness: every implementer is exp of an element of the current image.
    require(
        isinstance(model, ChargedDoubleTripletModel),
        "INNERNESS_MODEL",
        "innerness witnesses are defined for the charged double-triplet response",
    )
    frame_scale = model.lam_frame
    kernel_scale = model.lam_kernel
    require(
        not frame_scale.is_zero() and not kernel_scale.is_zero(),
        "INNERNESS_BLOCK_SKEW",
        "block skew-adjoint pairs are not in the current image",
    )
    # The odd map d -> (U d, sigma(U) d) must be a linear isomorphism onto R^3+R^3.
    odd_matrix: list[list[F5]] = []
    for a in range(6):
        odd_field = [Fraction(0)] * 12
        odd_field[frame.axis_reps[a]] = Fraction(1, 2)
        odd_field[carrier.antipode[frame.axis_reps[a]]] = Fraction(-1, 2)
        _, odd_coords = frame.even_odd(odd_field)
        u_image = frame.frame_map(odd_coords)
        w_image = frame.galois_frame_map(odd_coords)
        odd_matrix.append(list(u_image) + list(w_image))
    require(rank(odd_matrix) == 6, "INNERNESS_BLOCK_SKEW", "odd response bands do not span both skew blocks")

    innerness_rows = []
    order_of = e565.permutation_order
    for g in plus:
        rotation = frame.rotation_of(g)
        galois_rotation = [[rotation[i][j].conj() for j in range(3)] for i in range(3)]
        row_e = rotation_normal_form(rotation)
        row_w = rotation_normal_form(galois_rotation)
        innerness_rows.append(
            {
                "element_order": order_of(tuple(g)),
                "even_block": row_e,
                "kernel_block": row_w,
                "exp_witness": "Pi(g) = exp(diag(theta_E hat(n_E), theta_W hat(n_W))) with both blocks in the current image",
            }
        )
    require(len(innerness_rows) == 60, "INNERNESS_COUNT", "expected sixty innerness witnesses")

    # Register-relabeling no-go: the current lift provably cannot arise from
    # register relabeling.  Three exact facts, all basis-independent:
    # (1) the port action is a relabeling representation (non-negative
    #     integer character), and relabeled records generate only the abelian
    #     algebra that fails the gate (finite control);
    # (2) both charged-sector characters are irrational on every order-five
    #     element, so neither sector carries a signed register-relabeling
    #     action in any register basis;
    # (3) both sector characters have exact norm one (absolutely
    #     irreducible), and the element orders are exactly {1,2,3,5}, so no
    #     order-20 subgroup and hence no index-three subgroup exists: a
    #     three-dimensional irreducible monomial realization (relabeling with
    #     phases) is impossible.
    element_orders = sorted({order_of(tuple(g)) for g in plus})
    require(
        element_orders == [1, 2, 3, 5],
        "RELABELING_NO_GO",
        f"expected element orders [1, 2, 3, 5], got {element_orders}",
    )
    order_five_irrational = 0
    even_norm = ZERO
    kernel_norm = ZERO
    for g in plus:
        rotation = frame.rotation_of(g)
        even_trace = rotation[0][0] + rotation[1][1] + rotation[2][2]
        kernel_trace = even_trace.conj()
        even_norm = even_norm + even_trace * even_trace
        kernel_norm = kernel_norm + kernel_trace * kernel_trace
        if order_of(tuple(g)) == 5:
            require(
                even_trace.b != 0 and kernel_trace.b != 0,
                "RELABELING_NO_GO",
                "an order-five implementer has an integer character value",
            )
            order_five_irrational += 1
    require(order_five_irrational == 24, "RELABELING_NO_GO", "expected 24 order-five elements")
    group_order = F5(60)
    require(
        ((even_norm / group_order) - ONE).is_zero(),
        "RELABELING_NO_GO",
        "the even charged sector is not absolutely irreducible",
    )
    require(
        ((kernel_norm / group_order) - ONE).is_zero(),
        "RELABELING_NO_GO",
        "the kernel charged sector is not absolutely irreducible",
    )

    # Refinement naturality: every declared tower map is intertwined by K.
    refinement_row = e565.validate_refinement(carrier_manifest, carrier, plus, e565.gram_matrix(carrier))
    tower = carrier_manifest["refinement_tower"]
    naturality_maps = []
    for item in tower["maps"]:
        permutation = e565.parse_port_permutation(item["port_map"], carrier)
        pi = model.implementer(permutation)
        pi_dagger = tuple(cdagger(block) for block in pi)
        ginv = inverse(permutation)
        for field in BASIS_FIELDS:
            moved = [field[ginv[p]] for p in range(12)]
            left = model.generator(moved)
            right = model.generator(field)
            for b in range(len(pi)):
                conjugated = cmul(cmul(pi[b], right[b]), pi_dagger[b])
                require(
                    c_is_zero(csub(left[b], conjugated)),
                    "REFINEMENT_NATURALITY",
                    "a refinement map is not intertwined by the current lift",
                )
        naturality_maps.append({"source": item["source"], "target": item["target"], "intertwined": True})

    # Physical refinement maps: every artifact-bound persistence map is a
    # port permutation preserving the incidence and intertwined by K.
    physical_naturality_maps = []
    if artifact_binding is not None:
        for row in artifact_binding["port_maps"]:
            permutation = tuple(row["port_map"])
            pi = model.implementer(permutation)
            pi_dagger = tuple(cdagger(block) for block in pi)
            ginv = inverse(permutation)
            for field in BASIS_FIELDS:
                moved = [field[ginv[p]] for p in range(12)]
                left = model.generator(moved)
                right = model.generator(field)
                for b in range(len(pi)):
                    conjugated = cmul(cmul(pi[b], right[b]), pi_dagger[b])
                    require(
                        c_is_zero(csub(left[b], conjugated)),
                        "ARTIFACT_REFINEMENT",
                        "a physical refinement map is not intertwined by the derived current lift",
                    )
            physical_naturality_maps.append(
                {
                    "source_level": row["source_level"],
                    "target_level": row["target_level"],
                    "origin": row["origin"],
                    "intertwined": True,
                }
            )

    # Equivariant response moduli: dim Hom_A5(P12, g) via the permutation rank.
    fixed_point_squares = sum(sum(1 for p in range(12) if g[p] == p) ** 2 for g in plus)
    require(fixed_point_squares % 60 == 0, "MODULI_ARITHMETIC", "Burnside sum is not divisible by the group order")
    moduli_dimension = fixed_point_squares // 60
    require(
        moduli_dimension == 4,
        "MODULI_DIMENSION",
        f"expected a four-dimensional equivariant response moduli, got {moduli_dimension}",
    )

    # Schur rigidity: the kernel band admits no equivariant image in the even
    # block, so the block allocation is forced.  Character arithmetic:
    # multiplicity of the kernel band inside the even-block operator module.
    kernel_character: list[F5] = []
    even_block_character: list[F5] = []
    for g in plus:
        rotation = frame.rotation_of(g)
        galois_rotation = [[rotation[i][j].conj() for j in range(3)] for i in range(3)]
        chi_3 = rotation[0][0] + rotation[1][1] + rotation[2][2]
        chi_3p = galois_rotation[0][0] + galois_rotation[1][1] + galois_rotation[2][2]
        kernel_character.append(chi_3p)
        # u(3) block as a real A5 module: 1 + quintet + vector = chi_3^2.
        even_block_character.append(chi_3 * chi_3)
    sixty = F5(60)
    pairing_total = ZERO
    for chi_a, chi_b in zip(kernel_character, even_block_character, strict=True):
        pairing_total = pairing_total + chi_a * chi_b
    schur_multiplicity = pairing_total / sixty
    require(
        schur_multiplicity.is_zero(),
        "SCHUR_RIGIDITY",
        "the kernel band unexpectedly embeds into the even block",
    )

    gate = {
        "injective_twelve_dimensional_image": True,
        "skew_adjoint": True,
        "commutator_closed": True,
        "compact_type_positive_definite_invariant_form": True,
        "central_u1_kernel_dimension_one": True,
        "derived_dimension_eleven": True,
        "adjoint_rank_eleven": True,
        "faithful_charged_response_space": True,
        "a5_covariant_icosahedral_intertwiner": True,
        "induced_a5_action_inner": True,
        "refinement_natural": True,
        "repairs_distinct_from_responses": True,
    }
    require(all(gate.values()), "GATE", "conditional algebraic gate did not pass")

    constraints_bound = artifact_binding is not None
    current_source_bound = bool(
        constraints_bound
        and artifact_binding["report"].get("current_lift_source_selected") is True
        and artifact_binding["report"].get(
            "commutator_reconstructed_from_ordered_response"
        )
        is True
        and artifact_binding["report"].get(
            "overlap_holonomy_internality_certified"
        )
        is True
    )
    physical_gate = {
        "response_model_source_bound": current_source_bound,
        "response_coefficients_source_bound": constraints_bound,
        "target_blind_impulse_readback_recomputed": constraints_bound,
        "response_constraint_refinement_source_bound": constraints_bound,
        "ordered_response_commutator_reconstructed": current_source_bound,
        "a2_overlap_holonomy_fullness_source_bound": current_source_bound,
        "same_current_internal_implementers_source_bound": current_source_bound,
        "physical_refinement_intertwining_source_bound": current_source_bound
        and bool(physical_naturality_maps),
        "passed": current_source_bound and bool(physical_naturality_maps),
    }

    return {
        "schema": RECEIPT_SCHEMA,
        "issue": 566,
        "manifest_sha256": sha256_json(manifest),
        "carrier_manifest_sha256": sha256_json(carrier_manifest),
        "source_firewall": {
            "forbidden_dependency_hits": [],
            "uses_only": [
                "certified twelve-port carrier packet",
                "declared charged-double-triplet algebraic current fixture",
                "reversible response automorphism typing",
                "hash-pinned response-constraint artifact deriving R=-J and its "
                "four exact sector signs from carrier incidence",
                "irreversible strict-descent repair ledger (excluded from currents)",
            ],
        },
        "source_definedness": {
            "domain": "real port fields on the twelve primitive central atoms of the certified carrier",
            "operators": (
                "K is the declared charged-double-triplet fixture, evaluated on "
                "the artifact-bound oriented frame and exact response signs"
            ),
            "inner_product": "standard Hermitian pairing on the charged double-triplet response space C^3 (+) C^3",
            "response_pairing": "Hilbert-Schmidt pullback -Re tr(K(f)K(f')) with exact band coefficients",
            "refinement_maps": (
                "the source-bound carrier persistence maps are intertwined by "
                "the declared current fixture; this verifies a conditional "
                "intertwining identity, not a source selection of K"
            ),
            "carrier_and_refinement_provenance": "derived from the hash-pinned certified carrier packet",
            "response_construction_status": params["response_status"],
            "response_data_provenance": (
                "the response signs and oriented carrier constraints are derived "
                "from the unique nonidentity central antipode involution; the "
                "charged-double-triplet matrix lift and its Lie bracket are declared"
            ),
            "carrier_and_response_constraint_maps_source_bound": constraints_bound,
            "response_model_declared_as_branch_premise": True,
            "physical_response_source_bound": current_source_bound,
            "algebraic_construction_verified": True,
        },
        "frame_realization": {
            "coordinate_model": "twelve unnormalized icosahedron vertices, cyclic permutations of (0, +/-1, +/-phi)",
            "arithmetic_field": "Q(sqrt5), no floating point in any proof decision",
            "orientation_matched_assignments": 60,
            "equivalence": "realizations form one proper orbit; exact band coefficients verified equal for an alternative orientation-matched realization",
            "canonical_assignment": list(psi),
            "assignment_source": (
                "semantic artifact port_vertex_frame"
                if constraints_bound
                else "first orientation-matched realization"
            ),
            "axis_representatives": [carrier.ports[p] for p in frame.axis_reps],
        },
        "port_to_generator_map": {
            "model": params["model"],
            "construction_provenance": (
                "declared_charged_double_triplet_fixture_constrained_by_source_response"
                if constraints_bound
                else "declared_construction_model"
            ),
            "signed_response_band_coefficients": {k: str(v) for k, v in params["scales"].items()},
            "injective": True,
            "image_real_dimension": image_rank,
            "skew_adjoint": True,
            "block_dimensions": {"even_block_u3": 9, "kernel_block_so3": 3},
            "block_dimensions_verified": block_dimensions_verified,
            "kernel_block_real_verified": True,
            "compact_lie_type": "u(3) (+) so(3) = u(1) (+) su(3) (+) su(2)",
            "band_realization": {
                "unit_band": "i * (scale) * identity on the even block: the central u(1)",
                "quintet_band": "i * traceless symmetric axis quadrupoles on the even block",
                "frame_band": "real cross-product generators hat(U d) on the even block",
                "kernel_band": "Galois-twisted cross-product generators hat(sigma(U) d) on the kernel block",
            },
        },
        "closure": {
            "commutator_closed": True,
            "structure_constants_field": "Q(sqrt5)",
            "structure_constants": structure_constants,
            "derived_dimension": derived_dimension,
            "derived_block_dimensions": {"even_block_su3": derived_even_rank, "kernel_block_so3": derived_kernel_rank},
            "derived_type_identification": "commutators are traceless, so the eight-dimensional even derived block is exactly su(3) and the three-dimensional kernel derived block is so(3)",
            "center_dimension": len(center_basis),
            "center_is_constant_even_port_line": True,
            "adjoint_rank": adjoint_rank,
            "central_u1_kernel": "the adjoint action kills exactly the constant even port line",
        },
        "compactness": {
            "pullback_form": "-Re tr(K(f) K(f')) on port fields",
            "positive_definite": True,
            "elimination_pivots": pivot_texts,
            "hilbert_schmidt_pullback_band_coefficients": band_coefficients,
            "a5_invariant": True,
            "conclusion": "the current algebra carries an invariant positive-definite form, hence is of compact type",
        },
        "charged_faithfulness": {
            "matrix_faithful": True,
            "central_charge": "the central generator acts as i*(unit scale) on the even sector and 0 on the kernel sector",
            "charged": True,
            "equivariant_charge_rigidity": "an equivariant lift either charges the response space through the scalar central character or degenerates to a rank-deficient map",
        },
        "icosahedral_intertwiner": {
            "covariance_checks": covariance["pairs_checked"],
            "implementer_homomorphism_pairs": homomorphism_pairs,
            "implementers_faithful": True,
            "kernel_band_schur_rigidity": {
                "multiplicity_of_kernel_band_in_even_block": "0",
                "conclusion": "the kernel band cannot act on the even sector; the block allocation is forced",
            },
        },
        "inner_action": {
            "block_skew_pairs_in_image": True,
            "witness_count": len(innerness_rows),
            "witnesses": innerness_rows,
            "conclusion": "every implementer is exp of an element of the current image, so the induced A5 action lies in Int(g)",
        },
        "response_versus_register_relabeling": {
            "port_action_is_register_relabeling": "the port action is a permutation representation; relabeled records generate only the abelian algebra, which fails the gate (finite control abelian_record_model)",
            "element_orders": element_orders,
            "order_five_elements_with_irrational_sector_characters": order_five_irrational,
            "sector_character_norms": {"even_block": "1", "kernel_block": "1"},
            "no_index_three_subgroup": "element orders {1,2,3,5} exclude every group of order 20, so the derived action has no index-three subgroup",
            "no_go": "both charged sectors are absolutely irreducible with irrational order-five characters and no index-three subgroup exists, so neither sector admits a signed or phased register-relabeling (monomial) realization in any register basis",
            "conclusion": "the current lift factors through genuinely non-relabeling response sectors; register relabeling provably cannot generate these currents",
        },
        "refinement": {
            "carrier_tower": refinement_row,
            "naturality": naturality_maps,
            "physical_naturality": physical_naturality_maps,
            "natural": True,
        },
        "response_moduli": {
            "equivariant_lift_dimension": moduli_dimension,
            "burnside_rank_check": "sum of squared fixed-port counts over A5 equals 240 = 4 * 60",
            "source_data_status": (
                "response_signs_determined_current_lift_open"
                if constraints_bound
                else "open"
            ),
            "band_coefficient_provenance": (
                {
                    "unit_band": "negative eigenspace of the conventional representative R=-J",
                    "quintet_band": "negative eigenspace of the conventional representative R=-J",
                    "frame_band": "positive eigenspace of the conventional representative R=-J",
                    "kernel_band": "positive eigenspace of the conventional representative R=-J",
                }
                if constraints_bound
                else {
                    "unit_band": "signed unit-band response coefficient (open)",
                    "quintet_band": "signed quintet-band response coefficient (open)",
                    "frame_band": "signed frame-band response coefficient (open)",
                    "kernel_band": "signed kernel-band response coefficient (open)",
                }
            ),
        },
        "repair_response_distinction": {
            "reversible_response_automorphisms_define_currents": True,
            "responses_closed_under_inverse": True,
            "strict_descent_repairs_typed_irreversible": True,
            "repair_ledger_rows": params["repair_ledger_rows"],
            "disjoint": True,
        },
        "classification_vs_realization": {
            "coefficient_layer": "the compact-Lie trichotomy classifies coefficient algebras; coefficient records can commute",
            "conditional_current_layer": "given the declared charged-double-triplet construction, this receipt constructs current operators with nonabelian closure and an inner A5 action",
            "separating_witness": "the abelian record model matches the coefficient module yet fails the conditional algebraic gate with CENTER_NOT_ONE_DIMENSIONAL",
            "distinguished": True,
        },
        "conditional_algebraic_gate": {**gate, "passed": True},
        "physical_source_gate": physical_gate,
        "semantic_response_binding": (
            artifact_binding["report"] if constraints_bound else None
        ),
        "lean_cross_check": {
            "module": "Lean/Screen/A5IncidenceResponse.lean",
            "declarations": [
                "OPH.A5IncidenceResponse.distance_three_partner_unique",
                "OPH.A5IncidenceResponse.antipode_polynomial",
            ],
            "scope": (
                "independently checks unique graph-distance-three readback and "
                "the integral identity 10J=A^3-4A^2-5A+10I; the operational "
                "impulse/readback contract and laboratory attachment remain "
                "outside these Lean incidence theorems"
            ),
            "standard_axioms_only": True,
        },
        "derivation_chain": [
            {
                "step": 1,
                "premise": (
                    "response manifest with firewall, repair/response typing, a "
                    "declared matrix-current fixture, and the hash-pinned "
                    "response-constraint artifact"
                ),
                "uses": ["schema validation", "forbidden-token firewall", "reversible/irreversible typing split"],
                "source_artifact": "manifests/port_current_response_reference.json",
                "conclusion": (
                    "the conditional algebraic packet is admissible: the "
                    "charged-double-triplet construction is declared, its exact "
                    "response constraints are bound, and repairs are excluded "
                    "from currents"
                ),
            },
            {
                "step": "1a",
                "premise": "semantic response artifact derived from the finite carrier incidence and bound to the runtime dynamics",
                "uses": [
                    "artifact self-hash and manifest pin",
                    "paper-side recomputation of the isotypic projectors, Galois pairing, unique central antipode involution, tight-frame constant, and rotation group data",
                    "exact R=-J eigenspace and coefficient equality between the manifest and the artifact",
                    "independent Lean checks of the unique distance-three partner and 10J adjacency-polynomial identity",
                ],
                "source_artifact": (
                    "manifests/charged_response_semantic_artifact.json"
                    if constraints_bound
                    else "absent: declared control lane"
                ),
                "conclusion": (
                    "the four signed response constraints, oriented frame, and "
                    "carrier persistence maps are determined by finite carrier "
                    "structure; no current representation or Lie bracket is selected"
                    if constraints_bound
                    else "no semantic response-constraint artifact is bound"
                ),
            },
            {
                "step": 2,
                "premise": "hash-pinned certified carrier manifest",
                "uses": ["echosahedral_selector_certificate.validate_carrier", "group_certificate"],
                "source_artifact": "manifests/echosahedral_federation_reference.json",
                "conclusion": "unit split, antipode, faithful proper A5 action (60 permutations), distances, oriented faces, refinement tower re-derived",
            },
            {
                "step": 3,
                "premise": "exact coordinate model: cyclic permutations of (0, +/-1, +/-phi) over Q(sqrt5)",
                "uses": ["distance-isometry enumeration", "exact face determinant signs"],
                "source_artifact": "standard_vertices()",
                "conclusion": "120 distance isometries, exactly 60 orientation-matched: one proper orbit of oriented frame realizations",
            },
            {
                "step": 4,
                "premise": "antipode band split and the frame map U d = (1/2) sum_p f_p v_p",
                "uses": ["Galois automorphism sigma(sqrt5) = -sqrt5 applied entrywise"],
                "source_artifact": "FrameRealization.frame_map / galois_frame_map",
                "conclusion": "the four isotypic response bands are separated inside Q(sqrt5); sigma(U) kills the frame band and isolates the kernel band",
            },
            {
                "step": 5,
                "premise": (
                    "the declared charged-double-triplet fixture with exact band "
                    "scales constrained by the response artifact"
                ),
                "uses": ["skew-adjointness check", "exact rank over Q(sqrt5)"],
                "source_artifact": "ChargedDoubleTripletModel.generator",
                "conclusion": "injective port-to-generator map, image real dimension 12 with verified block dimensions (9, 3) and real kernel block",
            },
            {
                "step": 6,
                "premise": "implementers solved exactly from the frame transport R_g v_p = v_{g(p)}",
                "uses": ["720 covariance identities", "3600 homomorphism products", "faithfulness"],
                "source_artifact": "FrameRealization.rotation_of",
                "conclusion": "K(g.f) = Pi(g) K(f) Pi(g)* for the full derived icosahedral action",
            },
            {
                "step": 7,
                "premise": "all 66 basis brackets solved in the span",
                "uses": ["exact structure constants", "centralizer null space", "derived block ranks (8, 3)"],
                "source_artifact": "solve_in_span",
                "conclusion": "commutator-closed algebra with one-dimensional center on the constant even port line, derived dimension 11, adjoint rank 11; derived type su(3) (+) so(3)",
            },
            {
                "step": 8,
                "premise": "Hilbert-Schmidt pullback of the response pairing",
                "uses": ["twelve positive elimination pivots", "band-scalar reconstruction", "A5 invariance", "assignment independence"],
                "source_artifact": "certificate_payload",
                "conclusion": "invariant positive-definite form: compact type, with exact band coefficients",
            },
            {
                "step": 9,
                "premise": "odd bands jointly surject onto both skew blocks (exact rank six)",
                "uses": ["60 exact Rodrigues rotation normal forms in Q(sqrt5)"],
                "source_artifact": "rotation_normal_form",
                "conclusion": "every implementer is exp of an element of the current image, so the induced A5 action lies in Int(g)",
            },
            {
                "step": 10,
                "premise": "declared refinement tower maps",
                "uses": ["per-map intertwining checks", "carrier tower cocycle"],
                "source_artifact": "echosahedral_selector_certificate.validate_refinement",
                "conclusion": "the current construction is natural along the declared refinement tower",
            },
            {
                "step": 11,
                "premise": "Burnside count of the port action and exact character arithmetic",
                "uses": ["sum of squared fixed-port counts = 240", "kernel-band multiplicity 0 in the even block"],
                "source_artifact": "certificate_payload",
                "conclusion": "equivariant lifts form exactly a four-dimensional family of signed band coefficients; the block allocation is forced",
            },
            {
                "step": 12,
                "premise": "exact character arithmetic and element orders of the derived action",
                "uses": ["irrational order-five sector characters", "character norms one", "element orders {1,2,3,5}"],
                "source_artifact": "certificate_payload",
                "conclusion": "register-relabeling no-go: neither charged sector admits a monomial realization, so the currents cannot come from register relabeling",
            },
            {
                "step": 13,
                "premise": "gate aggregation and finite countermodels",
                "uses": ["typed negative controls"],
                "source_artifact": "negative_controls/issue_566_negative_controls.json",
                "conclusion": (
                    "the conditional algebraic gate passes on the declared "
                    "fixture and fails on every algebraic countermodel; the "
                    "physical current-source gate remains false"
                ),
            },
        ],
        "factor_origins": {
            "band_coefficients_1/4_3+sqrt5_5+sqrt5_5-sqrt5": "traces of the pullback form against the exact band projectors at unit scales; the frame/kernel pair is Galois-conjugate because the kernel band is realized through the Galois-twisted frame intertwiner",
            "dimensions_12_9_3": "exact ranks of the flattened generators and their block projections over Q(sqrt5)",
            "derived_11_center_1_adjoint_rank_11": "rank of the bracket span, nullity of the exact centralizer system, and their difference",
            "counts_120_60": "distance isometries onto the coordinate model and the orientation-matched proper subset",
            "checks_720_3600": "60 automorphisms times 12 basis fields; 60 times 60 implementer products",
            "moduli_4": "Burnside sum 240 of squared fixed-port counts divided by the group order 60",
            "order_five_cosines_(-1+-sqrt5)/4": "traces of the order-five rotation implementers, Galois-paired across the two blocks",
        },
        "branch_scope": {
            "branch": (
                "declared charged-double-triplet current fixture on the "
                "source-bound echosahedral response-constraint branch"
            ),
            "carrier": "the certified quotient-visible twelve-port carrier lineage of the pinned reference manifest",
            "response_data": (
                "the charged-double-triplet matrix lift is a declared algebraic "
                "fixture; the four relative A5-equivariant response signs and "
                "carrier maps are source-bound constraints"
            ),
            "not_claimed": (
                "no source selection of the current representation or Lie "
                "bracket, no statement about arbitrary OPH carriers, and no "
                "identification with the physical Standard Model gauge group"
            ),
        },
        "acceptance_criteria_status": {
            "current_operators_and_physical_refinement_source_defined": current_source_bound,
            "closure_compactness_rank_faithfulness_icosahedral_intertwiner_proved": True,
            "abelian_record_and_rank_deficient_models_fail_physical_current_gate": True,
            "coefficient_classification_distinguished_from_physical_current_realization": True,
            "no_measured_coupling_particle_assignment_or_standard_model_current_input": True,
        },
        "issue_closure_condition": {
            "produced_locally": (
                "the conditional full-rank compact skew-adjoint "
                "commutator-closed algebraic lift for the declared fixture, "
                "with inner A5 action and exact covariance over Q(sqrt5)"
            ),
            "response_field_provenance": params["response_status"],
            "conditional_algebraic_gate_passed": True,
            "physical_source_realization_gate_passed": physical_gate["passed"],
            "met_locally": physical_gate["passed"],
            "remaining_producer": (
                "a source producer must reconstruct twelve current generators "
                "and their commutator from ordered response histories, realize "
                "every proper carrier recharting by same-current closed overlap "
                "holonomy, and verify refinement intertwining"
            ),
        },
        "dependency_acyclicity_note": {
            "upstream": [
                "manifests/echosahedral_federation_reference.json and its receipt (strictly upstream carrier packet)",
                "the exact Q(sqrt5) coordinate model and the declared response manifest",
            ],
            "downstream": [
                "a5_screen_sm_closure.py and exterior_sm_completion.py reference this closure in their gate ledgers; this receipt does not consume them",
            ],
            "summary": "the proof-level dependency graph is acyclic: carrier packet -> response packet -> current receipt -> ledger references",
        },
        "verifier_command": "python3 code/a5_closure/port_current_inner_certificate.py verify --manifest code/a5_closure/manifests/port_current_response_reference.json --receipt code/a5_closure/receipts/port_current_inner_reference.receipt.json",
        "claim_boundary": {
            "proves": (
                "the conditional exact port-current algebra for the declared "
                "charged-double-triplet response construction"
            ),
            "status": "proved_conditional_on_declared_response_representation",
            "does_not_close": [
                "PORT-CURRENT-INNER as a physical source-bound receipt",
                "reconstruction of the current generators and bracket from ordered physical response",
                "A2 overlap-holonomy fullness and same-current internal implementers",
                "physical refinement intertwining beyond the declared algebraic tower maps",
                "block determinant balance and PORT-SPIN-LIFT",
                "physical Z6 deck/line descent (AXIS-CENTER-DESCENT)",
                "matter attachment, family structure, and exterior package selection",
                "continuum Yang-Mills quantum field theory, couplings, masses, or any measured number",
            ],
        },
    }


# ---------------------------------------------------------------------------
# Negative controls
# ---------------------------------------------------------------------------


def _control_lane(manifest: Mapping[str, Any], model: str) -> dict[str, Any]:
    """A declared-lane mutant exercising the algebraic countermodels."""

    mutant = copy.deepcopy(manifest)
    mutant["construction_model"] = model
    mutant.pop("semantic_response_artifact", None)
    return mutant


def _mutated_artifact(
    artifact: Mapping[str, Any], mutate: Callable[[dict[str, Any]], None], *, rehash: bool = True
) -> dict[str, Any]:
    mutated = copy.deepcopy(dict(artifact))
    mutate(mutated)
    if rehash:
        mutated["artifact_sha256"] = artifact_self_hash(mutated)
    return mutated


def _with_inline_artifact(
    manifest: Mapping[str, Any], artifact: Mapping[str, Any]
) -> dict[str, Any]:
    mutant = copy.deepcopy(manifest)
    mutant["semantic_response_artifact"] = {"value": artifact}
    return mutant


def negative_control_cases(
    manifest: Mapping[str, Any], base_dir: Path
) -> list[tuple[str, dict[str, Any], str]]:
    cases: list[tuple[str, dict[str, Any], str]] = []

    cases.append(
        ("abelian_record_model", _control_lane(manifest, "abelian_record"), "CENTER_NOT_ONE_DIMENSIONAL")
    )

    rank_deficient = _control_lane(manifest, "charged_double_triplet")
    rank_deficient["response_band_scales"]["kernel_band"] = "0"
    cases.append(("rank_deficient_kernel_band", rank_deficient, "IMAGE_RANK_DEFICIENT"))

    dead_center = _control_lane(manifest, "charged_double_triplet")
    dead_center["response_band_scales"]["unit_band"] = "0"
    cases.append(("rank_deficient_dead_center", dead_center, "IMAGE_RANK_DEFICIENT"))

    non_equivariant = _control_lane(manifest, "charged_double_triplet")
    non_equivariant["even_quintet_axis_scales"] = ["2", "1", "1", "1", "1", "1"]
    cases.append(("non_equivariant_axis_response", non_equivariant, "COVARIANCE_BROKEN"))

    uncommon_sign = _control_lane(manifest, "charged_double_triplet")
    uncommon_sign["odd_axis_signs"] = [1, 1, 1, 1, 1, -1]
    cases.append(("odd_axis_sign_not_common", uncommon_sign, "COVARIANCE_BROKEN"))

    cases.append(
        ("symmetric_record_pairing", _control_lane(manifest, "symmetric_record_control"), "SKEW_ADJOINTNESS_BROKEN")
    )

    conflated = copy.deepcopy(manifest)
    conflated["strict_descent_repairs"]["ledger"][0]["defines_currents"] = True
    cases.append(("repair_conflated_with_response", conflated, "REPAIR_RESPONSE_CONFLATION"))

    relabeling = copy.deepcopy(manifest)
    relabeling["response_declaration_contract"]["distinct_from_register_relabeling"] = False
    cases.append(("register_relabeling_conflated_as_response", relabeling, "REGISTER_RELABELING_CONFLATION"))

    unsupported_measurement = copy.deepcopy(manifest)
    unsupported_measurement["response_declaration_contract"]["measurement_artifact"] = {
        "path": "unrelated.json",
        "sha256": "0" * 64,
    }
    cases.append(("unsupported_measurement_upgrade", unsupported_measurement, "RESPONSE_ARTIFACT"))

    forbidden = copy.deepcopy(manifest)
    forbidden["downstream_hint"] = {"measured_coupling_target": "alpha_inverse"}
    cases.append(("inject_downstream_target", forbidden, "FORBIDDEN_DEPENDENCY"))

    production_model_string = copy.deepcopy(manifest)
    production_model_string.pop("construction_model", None)
    cases.append(
        ("production_missing_construction_fixture", production_model_string, "CONSTRUCTION_MODEL_STRING")
    )
    # This control, like all artifact mutations, runs through production validation.

    artifact_ref = manifest.get("semantic_response_artifact")
    if isinstance(artifact_ref, Mapping) and "path" in artifact_ref:
        artifact = load_semantic_artifact(artifact_ref, base_dir, False)

        unrelated = copy.deepcopy(manifest)
        unrelated["semantic_response_artifact"] = {
            "value": {"schema": "unrelated.receipt.v1", "rows": []}
        }
        cases.append(("artifact_unrelated_json", unrelated, "ARTIFACT_SCHEMA"))

        tampered = _mutated_artifact(
            artifact,
            lambda a: a["derived"].update(construction="abelian_record"),
            rehash=False,
        )
        cases.append(
            ("artifact_self_hash_tamper", _with_inline_artifact(manifest, tampered), "ARTIFACT_HASH")
        )

        wrong_carrier = _mutated_artifact(
            artifact,
            lambda a: a["carrier_binding"].update(carrier_manifest_sha256="0" * 64),
        )
        cases.append(
            ("artifact_wrong_carrier", _with_inline_artifact(manifest, wrong_carrier), "ARTIFACT_CARRIER")
        )

        missing_sector = _mutated_artifact(
            artifact,
            lambda a: a["response_basis"]["sector_dimensions"].pop("kernel_band"),
        )
        cases.append(
            ("artifact_missing_sector", _with_inline_artifact(manifest, missing_sector), "ARTIFACT_SECTORS")
        )

        def _swap_channels(a: dict[str, Any]) -> None:
            channels = a["response_basis"]["adjacency_channel_values"]
            channels["frame_band"], channels["kernel_band"] = (
                channels["kernel_band"],
                channels["frame_band"],
            )

        swapped = _mutated_artifact(artifact, _swap_channels)
        cases.append(
            ("artifact_orientation_swap", _with_inline_artifact(manifest, swapped), "ARTIFACT_SECTORS")
        )

        mismatched = _mutated_artifact(
            artifact,
            lambda a: a["derived"]["response_band_scales"].update(frame_band="7"),
        )
        cases.append(
            ("artifact_coefficient_mismatch", _with_inline_artifact(manifest, mismatched), "COEFFICIENT_MISMATCH")
        )

        false_source_selection = _mutated_artifact(
            artifact,
            lambda a: a["derived"]["current_lift_status"].update(
                source_selected=True
            ),
        )
        cases.append(
            (
                "artifact_forged_current_source_selection",
                _with_inline_artifact(manifest, false_source_selection),
                "ARTIFACT_CURRENT_BOUNDARY",
            )
        )

        def _doctor_map(a: dict[str, Any]) -> None:
            row = a["physical_refinement_maps"]["port_persistence_maps"][0]
            doctored_map = list(range(12))
            doctored_map[0], doctored_map[1] = doctored_map[1], doctored_map[0]
            row["port_map"] = doctored_map
            body = {key: value for key, value in row.items() if key != "map_hash"}
            row["map_hash"] = "sha256:" + sha256_json(body)

        doctored = _mutated_artifact(artifact, _doctor_map)
        cases.append(
            ("artifact_doctored_refinement_map", _with_inline_artifact(manifest, doctored), "ARTIFACT_REFINEMENT")
        )

        def _flip_sign(a: dict[str, Any]) -> None:
            a["derived"]["response_band_scales"]["unit_band"] = "1"

        sign_flip = _mutated_artifact(artifact, _flip_sign)
        cases.append(
            ("artifact_derived_sign_flip", _with_inline_artifact(manifest, sign_flip), "ARTIFACT_ORIENTATION")
        )

        old_schema = _mutated_artifact(
            artifact,
            lambda a: a.update(schema="oph.charged_response_semantic_artifact.v2"),
        )
        cases.append(
            ("artifact_v2_schema_rejected", _with_inline_artifact(manifest, old_schema), "ARTIFACT_SCHEMA")
        )

        def _doctor_antipode(a: dict[str, Any]) -> None:
            antipode = a["source_response"]["antipode_port_map"]
            antipode[0] = 0

        doctored_antipode = _mutated_artifact(artifact, _doctor_antipode)
        cases.append(
            (
                "artifact_doctored_source_antipode",
                _with_inline_artifact(manifest, doctored_antipode),
                "ARTIFACT_RESPONSE",
            )
        )

        source_sign_flip = _mutated_artifact(
            artifact,
            lambda a: a["source_response"]["sector_eigenvalues"].update(unit_band=1),
        )
        cases.append(
            (
                "artifact_source_eigenvalue_flip",
                _with_inline_artifact(manifest, source_sign_flip),
                "ARTIFACT_RESPONSE",
            )
        )

        protocol_filter_tamper = _mutated_artifact(
            artifact,
            lambda a: a["source_response"]["impulse_readback_protocol"].update(
                homogeneous_filter_coefficients=["1", "0", "0", "0"]
            ),
        )
        cases.append(
            (
                "artifact_impulse_filter_tamper",
                _with_inline_artifact(manifest, protocol_filter_tamper),
                "ARTIFACT_RESPONSE",
            )
        )

    return cases


CONTROL_LANE_CASES = {
    "abelian_record_model",
    "rank_deficient_kernel_band",
    "rank_deficient_dead_center",
    "non_equivariant_axis_response",
    "odd_axis_sign_not_common",
    "symmetric_record_pairing",
}


def negative_control_payload(manifest: Mapping[str, Any], base_dir: Path | None = None) -> dict[str, Any]:
    results: list[dict[str, Any]] = []
    for name, mutant, expected_code in negative_control_cases(manifest, base_dir or MODULE_DIR):
        actual_code = "ACCEPTED"
        try:
            certificate_payload(
                mutant,
                base_dir,
                allow_control_models=name in CONTROL_LANE_CASES,
                allow_inline_artifact_for_tests=name.startswith("artifact_"),
            )
        except CertificateError as exc:
            actual_code = exc.code
        require(
            actual_code == expected_code,
            "NEGATIVE_CONTROL_FAILED",
            f"{name}: expected {expected_code}, got {actual_code}",
        )
        results.append({"name": name, "expected_error": expected_code, "actual_error": actual_code, "passed": True})
    return {
        "schema": NEGATIVE_SCHEMA,
        "issue": 566,
        "manifest_sha256": sha256_json(manifest),
        "finite_controls": results,
        "countermodel_witnesses": {
            "abelian_record": {
                "model": "K(f) = i*diag(f) on C^12 with permutation implementers",
                "matches_coefficient_module": True,
                "derived_dimension": 0,
                "center_dimension": 12,
                "induced_action_on_abelian_algebra": "nontrivial, hence not inner: Int of an abelian algebra is trivial",
                "conclusion": "coefficient records that commute never pass the conditional algebraic gate",
            },
            "rank_deficient": {
                "kernel_band_scale_zero_image_dimension": 9,
                "unit_band_scale_zero_image_dimension": 11,
                "conclusion": "degenerate response coefficients cannot produce the twelve-dimensional algebraic current construction",
            },
            "equivariance": {
                "per_axis_rescaling": "breaks K(g.f) = Pi(g) K(f) Pi(g)* on any element moving the rescaled axis",
                "per_axis_sign_flip": "breaks equivariance and covariance",
            },
            "typing": {
                "repair_conflation": "an irreversible strict-descent repair declared as a current source fails closed",
                "relabeling_conflation": "a response contract not typed distinct from register relabeling fails closed",
                "unsupported_measurement_upgrade": "measurement data outside the reviewed semantic artifact schema is rejected",
                "firewall": "a measured-coupling target in the source manifest fails closed",
            },
            "artifact": {
                "construction_fixture": "a production manifest without an explicit charged-double-triplet fixture fails closed; the response artifact does not select a matrix current",
                "hash_binding": "self-hash tamper, wrong carrier pin, and unrelated JSON fail closed",
                "structure_binding": "all inline artifact mutations run through production validation; old schemas, missing sectors, swapped channels, forged source selection, source-antipode or eigensign tampering, coefficient mismatches, and doctored refinement maps fail closed",
            },
        },
    }


def verify_receipt(manifest: Mapping[str, Any], receipt: Mapping[str, Any], base_dir: Path | None = None) -> None:
    expected = certificate_payload(manifest, base_dir)
    require(receipt == expected, "RECEIPT_MISMATCH", "receipt is stale, malformed, or tampered")


def default_paths() -> tuple[Path, Path, Path]:
    return (
        MODULE_DIR / "manifests" / "port_current_response_reference.json",
        MODULE_DIR / "receipts" / "port_current_inner_reference.receipt.json",
        MODULE_DIR / "negative_controls" / "issue_566_negative_controls.json",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    certify = sub.add_parser("certify", help="create the deterministic exact receipt")
    certify.add_argument("--manifest", type=Path, required=True)
    certify.add_argument("--output", type=Path, required=True)
    verify = sub.add_parser("verify", help="recompute and compare a receipt")
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--receipt", type=Path, required=True)
    negative = sub.add_parser("negative-controls", help="run and write the finite countermodel bundle")
    negative.add_argument("--manifest", type=Path, required=True)
    negative.add_argument("--output", type=Path, required=True)
    all_cmd = sub.add_parser("all", help="regenerate receipt and negative controls at repository-default paths")
    all_cmd.add_argument("--manifest", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "certify":
        manifest = load_json(args.manifest)
        receipt = certificate_payload(manifest, args.manifest.resolve().parent.parent)
        write_json(args.output, receipt)
        print(json.dumps({"status": "PASS", "receipt": str(args.output), "sha256": sha256_json(receipt)}, indent=2))
    elif args.command == "verify":
        manifest = load_json(args.manifest)
        receipt = load_json(args.receipt)
        verify_receipt(manifest, receipt, args.manifest.resolve().parent.parent)
        print(json.dumps({"status": "PASS", "receipt": str(args.receipt)}, indent=2))
    elif args.command == "negative-controls":
        manifest = load_json(args.manifest)
        payload = negative_control_payload(manifest, args.manifest.resolve().parent.parent)
        write_json(args.output, payload)
        print(json.dumps({"status": "PASS", "negative_controls": str(args.output)}, indent=2))
    else:
        default_manifest, default_receipt, default_negative = default_paths()
        manifest_path = args.manifest or default_manifest
        manifest = load_json(manifest_path)
        write_json(default_receipt, certificate_payload(manifest))
        write_json(default_negative, negative_control_payload(manifest))
        print(
            json.dumps(
                {"status": "PASS", "receipt": str(default_receipt), "negative_controls": str(default_negative)},
                indent=2,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
