#!/usr/bin/env python3
"""Exact certificate for GitHub issue #566: the physical port-current algebra.

The input is a port-current response manifest.  It declares only:

* the certified twelve-port echosahedral carrier manifest (by path and hash);
* four A5-equivariant reversible-response band scales and one common
  odd-response sign, as exact rationals;
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
* the four-dimensional equivariant response moduli (four band scales and one
  common odd sign are exactly the open source data).

Abelian-record and rank-deficient response models fail the physical-current
gate with typed error codes.  No Standard Model representation, particle
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

SCHEMA = "oph.port_current_response_manifest.v1"
RECEIPT_SCHEMA = "oph.port_current_inner_receipt.v1"
NEGATIVE_SCHEMA = "oph.port_current_inner_negative_controls.v1"

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


def validate_manifest(manifest: Mapping[str, Any]) -> dict[str, Any]:
    e565.enforce_source_firewall(manifest)
    require(manifest.get("schema") == SCHEMA, "SCHEMA", f"expected {SCHEMA}")

    response = manifest.get("reversible_response_automorphisms")
    require(isinstance(response, Mapping), "RESPONSE_TYPING", "reversible_response_automorphisms is missing")
    require(response.get("reversible") is True, "RESPONSE_TYPING", "response automorphisms must be typed reversible")
    require(response.get("defines_currents") is True, "RESPONSE_TYPING", "response automorphisms must be the declared current source")

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

    model = manifest.get("response_model")
    require(
        model in ("charged_double_triplet", "abelian_record", "symmetric_record_control"),
        "RESPONSE_MODEL",
        f"unknown response model {model!r}",
    )

    scales_raw = manifest.get("response_band_scales")
    require(isinstance(scales_raw, Mapping), "RESPONSE_SCALES", "response_band_scales is missing")
    band_names = ("unit_band", "quintet_band", "frame_band", "kernel_band")
    require(set(scales_raw) == set(band_names), "RESPONSE_SCALES", f"scales must name exactly the bands {band_names}")
    scales = {name: parse_rational(scales_raw[name], "RESPONSE_SCALES") for name in band_names}

    sign = manifest.get("odd_response_sign")
    require(sign in (1, -1), "ODD_RESPONSE_SIGN", "odd_response_sign must be +1 or -1")

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
        "sign": int(sign),
        "axis_scales": axis_scales,
        "odd_axis_signs": [int(s) for s in odd_axis_signs_raw],
        "repair_ledger_rows": len(ledger),
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
        self.lam_kernel = F5(scales["kernel_band"]) * F5(params["sign"])
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


def certificate_payload(manifest: Mapping[str, Any], base_dir: Path | None = None) -> dict[str, Any]:
    base = base_dir or MODULE_DIR
    params = validate_manifest(manifest)
    carrier, group_row, plus, carrier_manifest = load_carrier(manifest, base)

    verts = standard_vertices()
    matched = orientation_matched_assignments(carrier, verts)
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
        f"physical currents need derived dimension eleven; records that commute give {derived_dimension}",
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
        "source_defined_operators_domain_pairing_refinement": True,
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
    require(all(gate.values()), "GATE", "physical-current gate did not pass")

    return {
        "schema": RECEIPT_SCHEMA,
        "issue": 566,
        "manifest_sha256": sha256_json(manifest),
        "carrier_manifest_sha256": sha256_json(carrier_manifest),
        "source_firewall": {
            "forbidden_dependency_hits": [],
            "uses_only": [
                "certified twelve-port carrier packet",
                "reversible response automorphism typing",
                "four exact response band scales",
                "one common odd-response sign",
                "irreversible strict-descent repair ledger (excluded from currents)",
            ],
        },
        "source_definedness": {
            "domain": "real port fields on the twelve primitive central atoms of the certified carrier",
            "operators": "K built from the derived oriented frame, the Galois-twisted kernel intertwiner, and the declared response scales",
            "inner_product": "standard Hermitian pairing on the charged double-triplet response space C^3 (+) C^3",
            "response_pairing": "Hilbert-Schmidt pullback -Re tr(K(f)K(f')) with exact band coefficients",
            "refinement_maps": "the declared carrier tower maps, each intertwined by the current lift",
            "all_source_defined": True,
        },
        "frame_realization": {
            "coordinate_model": "twelve unnormalized icosahedron vertices, cyclic permutations of (0, +/-1, +/-phi)",
            "arithmetic_field": "Q(sqrt5), no floating point in any proof decision",
            "orientation_matched_assignments": 60,
            "equivalence": "realizations form one proper orbit; exact band coefficients verified equal for an alternative orientation-matched realization",
            "canonical_assignment": list(psi),
            "axis_representatives": [carrier.ports[p] for p in frame.axis_reps],
        },
        "port_to_generator_map": {
            "model": params["model"],
            "response_band_scales": {k: str(v) for k, v in params["scales"].items()},
            "odd_response_sign": params["sign"],
            "injective": True,
            "image_real_dimension": image_rank,
            "skew_adjoint": True,
            "block_dimensions": {"even_block_u3": 9, "kernel_block_so3": 3},
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
        "refinement": {
            "carrier_tower": refinement_row,
            "naturality": naturality_maps,
            "natural": True,
        },
        "response_moduli": {
            "equivariant_lift_dimension": moduli_dimension,
            "burnside_rank_check": "sum of squared fixed-port counts over A5 equals 240 = 4 * 60",
            "open_source_data": [
                "unit band response scale",
                "quintet band response scale",
                "frame band response scale",
                "kernel band response scale",
                "one common odd-response sign",
            ],
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
            "physical_layer": "this receipt constructs current operators on a faithful charged response space with nonabelian closure and an inner A5 action",
            "separating_witness": "the abelian record model matches the coefficient module yet fails the physical-current gate with CENTER_NOT_ONE_DIMENSIONAL",
            "distinguished": True,
        },
        "physical_current_gate": {**gate, "passed": True},
        "claim_boundary": {
            "closes": "PORT-CURRENT-INNER on the declared echosahedral response branch",
            "declared_source_data": "the four A5-equivariant response band scales and the common odd-response sign are declared reversible-response source data, not derived from raw consensus dynamics",
            "does_not_close": [
                "derivation of the response scales or the response space from raw OPH consensus dynamics",
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


def negative_control_cases(manifest: Mapping[str, Any]) -> list[tuple[str, dict[str, Any], str]]:
    cases: list[tuple[str, dict[str, Any], str]] = []

    abelian = copy.deepcopy(manifest)
    abelian["response_model"] = "abelian_record"
    cases.append(("abelian_record_model", abelian, "CENTER_NOT_ONE_DIMENSIONAL"))

    rank_deficient = copy.deepcopy(manifest)
    rank_deficient["response_band_scales"]["kernel_band"] = "0"
    cases.append(("rank_deficient_kernel_band", rank_deficient, "IMAGE_RANK_DEFICIENT"))

    dead_center = copy.deepcopy(manifest)
    dead_center["response_band_scales"]["unit_band"] = "0"
    cases.append(("rank_deficient_dead_center", dead_center, "IMAGE_RANK_DEFICIENT"))

    non_equivariant = copy.deepcopy(manifest)
    non_equivariant["even_quintet_axis_scales"] = ["2", "1", "1", "1", "1", "1"]
    cases.append(("non_equivariant_axis_response", non_equivariant, "COVARIANCE_BROKEN"))

    uncommon_sign = copy.deepcopy(manifest)
    uncommon_sign["odd_axis_signs"] = [1, 1, 1, 1, 1, -1]
    cases.append(("odd_response_sign_not_common", uncommon_sign, "COVARIANCE_BROKEN"))

    symmetric = copy.deepcopy(manifest)
    symmetric["response_model"] = "symmetric_record_control"
    cases.append(("symmetric_record_pairing", symmetric, "SKEW_ADJOINTNESS_BROKEN"))

    conflated = copy.deepcopy(manifest)
    conflated["strict_descent_repairs"]["ledger"][0]["defines_currents"] = True
    cases.append(("repair_conflated_with_response", conflated, "REPAIR_RESPONSE_CONFLATION"))

    forbidden = copy.deepcopy(manifest)
    forbidden["downstream_hint"] = {"measured_coupling_target": "alpha_inverse"}
    cases.append(("inject_downstream_target", forbidden, "FORBIDDEN_DEPENDENCY"))

    return cases


def negative_control_payload(manifest: Mapping[str, Any], base_dir: Path | None = None) -> dict[str, Any]:
    results: list[dict[str, Any]] = []
    for name, mutant, expected_code in negative_control_cases(manifest):
        actual_code = "ACCEPTED"
        try:
            certificate_payload(mutant, base_dir)
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
                "conclusion": "coefficient records that commute never pass the physical-current gate",
            },
            "rank_deficient": {
                "kernel_band_scale_zero_image_dimension": 9,
                "unit_band_scale_zero_image_dimension": 11,
                "conclusion": "degenerate response scales cannot produce the twelve-dimensional physical current",
            },
            "equivariance": {
                "per_axis_rescaling": "breaks K(g.f) = Pi(g) K(f) Pi(g)* on any element moving the rescaled axis",
                "per_axis_sign_flip": "breaks the common odd-response sign and covariance",
            },
            "typing": {
                "repair_conflation": "an irreversible strict-descent repair declared as a current source fails closed",
                "firewall": "a measured-coupling target in the source manifest fails closed",
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
