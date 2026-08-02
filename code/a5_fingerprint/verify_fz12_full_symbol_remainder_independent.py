#!/usr/bin/env python3
"""Independent verifier for the target-free FZ-12 full-symbol packet.

The verifier deliberately does not import the producer.  It reconstructs the
thirty unit seam directions from the exact icosahedral carrier table over
``Q(sqrt(5))``, recomputes their even moment polynomials, and checks the
eighth-order identity and Taylor bounds encoded by the receipt.

This is a spatial-symbol certificate.  Nothing in this file identifies the
symbol with a physical frequency, field sector, clock, frame, or readout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable


HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
DEFAULT_BASE_GEOMETRY = RUNTIME / "a5_multipole_fixed_point_receipt.json"
DEFAULT_SOURCE = RUNTIME / "seam_current_edge_prediction_receipt.json"
DEFAULT_CUSTODY = RUNTIME / "fz12_custody_projection.json"
DEFAULT_REMAINDER_LEAN = HERE.parent.parent / "Lean" / "Screen" / "SeamCurrentEdge30Remainder.lean"
DEFAULT_RECEIPT = RUNTIME / "fz12_full_symbol_remainder_receipt.json"

BASE_GEOMETRY_BYTES = 8864
BASE_GEOMETRY_SHA256 = "d96b80c71a64f48bb7a2a7b2592bad5e122ce9c5f3ca3b4636af29693c426ecd"
SOURCE_BYTES = 9296
SOURCE_SHA256 = "0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915"
CUSTODY_BYTES = 3624
CUSTODY_SHA256 = "dfc6574a5bde9b576df6aeef2e03aba8602ec2723a0a4606e7942c405a57c643"
REMAINDER_LEAN_BYTES = 2968
REMAINDER_LEAN_SHA256 = "e5d90cf9d84a4b4394355f9704114937268b2bfd32c478c4f93e6c53cd493928"

BASE_GEOMETRY_SCHEMA = "oph.a5_multipole_fixed_point_receipt.v1"
BASE_GEOMETRY_STATUS = "EXACT_A5_FINGERPRINT_CERTIFICATE__PHYSICAL_MAP_OPEN"
SOURCE_SCHEMA = "oph.seam_current_edge_prediction_candidate.v1"
SOURCE_STATUS = (
    "EXACT_SOURCE_NATIVE_EDGE_RAY__PROSPECTIVE_PHYSICAL_BRANCH_UNARMED__"
    "PHYSICAL_PRODUCER_OPEN"
)
CUSTODY_SCHEMA = "oph.fz12.custody_projection.v1"
CUSTODY_STATUS = "FROZEN_FZ12_SOURCE_AND_CUSTODY_PINNED__NO_COMPARISON_DATA"

SCHEMA = "oph.fz12.full_symbol_remainder.v1"
STATUS = (
    "EXACT_TARGET_FREE_EDGE_SYMBOL_Q_LE_ONE_REMAINDER__"
    "PHYSICAL_FREQUENCY_AND_COMPARISON_OPEN"
)


class VerificationError(ValueError):
    """An input failed the independent full-symbol contract."""


def check(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def self_digest(value: dict[str, Any], field: str) -> str:
    body = {key: item for key, item in value.items() if key != field}
    return "sha256:" + sha256(canonical_json_bytes(body))


def load_fixed_parent(
    path: Path,
    *,
    byte_count: int,
    digest: str,
    schema: str,
    status: str,
    self_field: str,
) -> tuple[bytes, dict[str, Any]]:
    raw = path.read_bytes()
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as error:
        raise VerificationError(f"invalid parent JSON: {path.name}") from error
    check(isinstance(value, dict), f"parent root is not an object: {path.name}")
    check(raw == canonical_json_bytes(value), f"noncanonical parent: {path.name}")
    check(sha256(raw) == digest, f"parent raw hash drift: {path.name}")
    check(len(raw) == byte_count, f"parent byte-count drift: {path.name}")
    check(value.get("schema") == schema, f"parent schema drift: {path.name}")
    check(value.get("status") == status, f"parent status drift: {path.name}")
    check(
        value.get(self_field) == self_digest(value, self_field),
        f"parent self-digest drift: {path.name}",
    )
    return raw, value


def verify_base_range_evidence(value: dict[str, Any]) -> None:
    frame = value.get("cartesian_frame")
    check(isinstance(frame, dict), "base geometry cartesian frame missing")
    check(
        frame.get("vertex_normalization")
        == "I6 = 1 on all twelve unit vertices",
        "base geometry I6 normalization drift",
    )
    critical = value.get("critical_points")
    check(isinstance(critical, dict), "base geometry critical-point proof missing")
    check(
        critical.get("census")
        == {"maxima": 12, "minima": 20, "saddles": 30, "total": 62},
        "base geometry critical-point census drift",
    )
    check(
        "every tangent Hessian eigenvalue is nonzero"
        in str(critical.get("nondegeneracy", "")),
        "base geometry critical-point nondegeneracy drift",
    )
    orbits = critical.get("orbits")
    check(isinstance(orbits, list), "base geometry critical-point orbits missing")
    extrema = {
        row.get("orbit"): (row.get("type"), row.get("value"))
        for row in orbits
        if isinstance(row, dict)
    }
    check(
        extrema.get("vertex_pole") == ("maximum", "1+0*sqrt5")
        and extrema.get("vertex_ring") == ("maximum", "1+0*sqrt5"),
        "base geometry I6 maximum evidence drift",
    )
    check(
        extrema.get("face_high") == ("minimum", "-5/9+0*sqrt5")
        and extrema.get("face_low") == ("minimum", "-5/9+0*sqrt5"),
        "base geometry I6 minimum evidence drift",
    )


def load_fixed_lean_parent(path: Path) -> bytes:
    raw = path.read_bytes()
    check(len(raw) == REMAINDER_LEAN_BYTES, "remainder Lean byte-count drift")
    check(sha256(raw) == REMAINDER_LEAN_SHA256, "remainder Lean raw hash drift")
    try:
        source = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise VerificationError("remainder Lean source is not UTF-8") from error
    check(
        "theorem seamMoment8_eq (k : Vec3)" in source,
        "remainder Lean seamMoment8_eq theorem missing",
    )
    check(
        "seamMoment8 k = (10 / 3 : ℝ) * radiusSquared k ^ 4" in source
        and "- (8 / 15 : ℝ) * radiusSquared k * I6 k := by" in source,
        "remainder Lean seamMoment8_eq statement drift",
    )
    forbidden = ("sorry", "admit", "axiom seamMoment8_eq")
    check(
        not any(token in source for token in forbidden),
        "remainder Lean source contains an admitted theorem",
    )
    check(
        "#print axioms OPH.SeamCurrentEdge30Remainder.seamMoment8_eq" in source,
        "remainder Lean axiom audit missing",
    )
    return raw


# A scalar is represented as ``a + b sqrt(5)`` with exact rational a and b.
Q5 = tuple[Fraction, Fraction]
Vec3 = tuple[Q5, Q5, Q5]
Monomial = tuple[int, int, int]
Polynomial = dict[Monomial, Q5]

ZERO: Q5 = (Fraction(0), Fraction(0))
ONE: Q5 = (Fraction(1), Fraction(0))
PHI: Q5 = (Fraction(1, 2), Fraction(1, 2))


def q5(value: int | Fraction) -> Q5:
    return (Fraction(value), Fraction(0))


def q5_add(left: Q5, right: Q5) -> Q5:
    return (left[0] + right[0], left[1] + right[1])


def q5_neg(value: Q5) -> Q5:
    return (-value[0], -value[1])


def q5_sub(left: Q5, right: Q5) -> Q5:
    return q5_add(left, q5_neg(right))


def q5_mul(left: Q5, right: Q5) -> Q5:
    return (
        left[0] * right[0] + 5 * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def q5_scale(value: Q5, scalar: Fraction | int) -> Q5:
    factor = Fraction(scalar)
    return (factor * value[0], factor * value[1])


def q5_pow(value: Q5, exponent: int) -> Q5:
    check(exponent >= 0, "negative Q(sqrt(5)) exponent")
    result = ONE
    base = value
    power = exponent
    while power:
        if power & 1:
            result = q5_mul(result, base)
        base = q5_mul(base, base)
        power >>= 1
    return result


def vec_neg(value: Vec3) -> Vec3:
    return tuple(q5_neg(item) for item in value)  # type: ignore[return-value]


def vec_sub(left: Vec3, right: Vec3) -> Vec3:
    return tuple(q5_sub(a, b) for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def vec_scale(value: Vec3, scalar: Fraction | int) -> Vec3:
    return tuple(q5_scale(item, scalar) for item in value)  # type: ignore[return-value]


def dot(left: Vec3, right: Vec3) -> Q5:
    result = ZERO
    for a, b in zip(left, right, strict=True):
        result = q5_add(result, q5_mul(a, b))
    return result


def polynomial_add(left: Polynomial, right: Polynomial) -> Polynomial:
    result = dict(left)
    for monomial, coefficient in right.items():
        result[monomial] = q5_add(result.get(monomial, ZERO), coefficient)
        if result[monomial] == ZERO:
            del result[monomial]
    return result


def polynomial_scale(value: Polynomial, scalar: Fraction | int) -> Polynomial:
    result = {
        monomial: q5_scale(coefficient, scalar)
        for monomial, coefficient in value.items()
    }
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient != ZERO}


def polynomial_mul(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for (i, j, k), a in left.items():
        for (r, s, t), b in right.items():
            monomial = (i + r, j + s, k + t)
            result[monomial] = q5_add(
                result.get(monomial, ZERO), q5_mul(a, b)
            )
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient != ZERO}


def polynomial_pow(value: Polynomial, exponent: int) -> Polynomial:
    result: Polynomial = {(0, 0, 0): ONE}
    base = value
    power = exponent
    while power:
        if power & 1:
            result = polynomial_mul(result, base)
        base = polynomial_mul(base, base)
        power >>= 1
    return result


def linear_power(vector: Vec3, exponent: int) -> Polynomial:
    result: Polynomial = {}
    factorial = math.factorial(exponent)
    for i in range(exponent + 1):
        for j in range(exponent - i + 1):
            k = exponent - i - j
            multinomial = Fraction(
                factorial,
                math.factorial(i) * math.factorial(j) * math.factorial(k),
            )
            coefficient = q5_mul(
                q5_pow(vector[0], i),
                q5_mul(q5_pow(vector[1], j), q5_pow(vector[2], k)),
            )
            coefficient = q5_scale(coefficient, multinomial)
            if coefficient != ZERO:
                result[(i, j, k)] = coefficient
    return result


def evaluate_polynomial(value: Polynomial, point: Vec3) -> Q5:
    result = ZERO
    for (i, j, k), coefficient in value.items():
        term = q5_mul(
            coefficient,
            q5_mul(
                q5_pow(point[0], i),
                q5_mul(q5_pow(point[1], j), q5_pow(point[2], k)),
            ),
        )
        result = q5_add(result, term)
    return result


def port_vectors() -> tuple[Vec3, ...]:
    minus_one = q5(-1)
    plus_one = q5(1)
    minus_phi = q5_neg(PHI)
    return (
        (ZERO, minus_one, minus_phi),
        (minus_one, minus_phi, ZERO),
        (minus_phi, ZERO, minus_one),
        (plus_one, minus_phi, ZERO),
        (ZERO, plus_one, minus_phi),
        (minus_phi, ZERO, plus_one),
        (PHI, ZERO, minus_one),
        (ZERO, minus_one, PHI),
        (minus_one, PHI, ZERO),
        (PHI, ZERO, plus_one),
        (plus_one, PHI, ZERO),
        (ZERO, plus_one, PHI),
    )


SEAM_LEFT = (
    0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3,
    4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10,
)
SEAM_RIGHT = (
    1, 2, 3, 4, 6, 2, 3, 5, 7, 4, 5, 8, 6, 7, 9,
    6, 8, 10, 7, 8, 11, 9, 10, 9, 11, 10, 11, 10, 11, 11,
)


def unit_seam_directions() -> tuple[Vec3, ...]:
    ports = port_vectors()
    directions = tuple(
        vec_scale(vec_sub(ports[right], ports[left]), Fraction(1, 2))
        for left, right in zip(SEAM_LEFT, SEAM_RIGHT, strict=True)
    )
    check(len(directions) == 30, "independent seam table is not thirty rows")
    check(all(dot(direction, direction) == ONE for direction in directions), "independent seam direction is not unit norm")
    return directions


def moment_polynomial(directions: Iterable[Vec3], degree: int) -> Polynomial:
    result: Polynomial = {}
    for direction in directions:
        result = polynomial_add(result, linear_power(direction, degree))
    return result


def independent_geometry_checks() -> dict[str, Any]:
    directions = unit_seam_directions()
    r2: Polynomial = {
        (2, 0, 0): ONE,
        (0, 2, 0): ONE,
        (0, 0, 2): ONE,
    }
    m2 = moment_polynomial(directions, 2)
    m4 = moment_polynomial(directions, 4)
    m6 = moment_polynomial(directions, 6)
    m8 = moment_polynomial(directions, 8)

    check(m2 == polynomial_scale(r2, 10), "independent M2 identity failed")
    check(
        m4 == polynomial_scale(polynomial_pow(r2, 2), 6),
        "independent M4 identity failed",
    )

    # Rearranging M6 = (30/7) r^6 - (2/7) I6 defines the same normalized
    # homogeneous I6 polynomial without importing the producer's expression.
    i6 = polynomial_add(
        polynomial_scale(polynomial_pow(r2, 3), 15),
        polynomial_scale(m6, Fraction(-7, 2)),
    )
    expected_m8 = polynomial_add(
        polynomial_scale(polynomial_pow(r2, 4), Fraction(10, 3)),
        polynomial_scale(polynomial_mul(r2, i6), Fraction(-8, 15)),
    )
    check(m8 == expected_m8, "independent M8 identity failed")

    for port in port_vectors():
        radius_squared = dot(port, port)
        check(
            evaluate_polynomial(i6, port) == q5_pow(radius_squared, 3),
            "independent I6 vertex normalization failed",
        )

    projective_representatives: list[Vec3] = []
    projective_counts: list[int] = []
    for direction in directions:
        matched = False
        for index, representative in enumerate(projective_representatives):
            if direction == representative or direction == vec_neg(representative):
                projective_counts[index] += 1
                matched = True
                break
        if not matched:
            projective_representatives.append(direction)
            projective_counts.append(1)
    check(len(projective_representatives) == 15, "independent projective axis count failed")
    check(
        projective_counts == [2] * 15,
        "independent projective axis multiplicity failed",
    )

    # Recomputed failures behind the declared mutation controls.
    m2_sixty = moment_polynomial(directions + directions, 2)
    m2_fifteen = moment_polynomial(projective_representatives, 2)
    raw_differences = tuple(vec_scale(direction, 2) for direction in directions)
    m2_raw = moment_polynomial(raw_differences, 2)
    check(m2_sixty == polynomial_scale(r2, 20), "60-label control did not double M2")
    check(m2_fifteen == polynomial_scale(r2, 5), "15-axis control did not halve M2")
    check(m2_raw == polynomial_scale(r2, 40), "raw norm-two control did not quadruple M2")

    return {
        "direction_count": len(directions),
        "projective_axis_count": len(projective_representatives),
        "raw_difference_norm_squared": "4",
        "M2": "10",
        "M4": "6",
        "M6": "30/7-(2/7)I6(n)",
        "M8": "10/3-(8/15)I6(n)",
        "controls": {
            "double_directed_labels_60": True,
            "projective_axes_15": True,
            "raw_norm2_differences": True,
        },
    }


def expected_parent_pins(
    parents: list[tuple[bytes, dict[str, Any]]], lean_raw: bytes,
) -> list[dict[str, Any]]:
    metadata = [
        (
            "code/a5_fingerprint/runtime/a5_multipole_fixed_point_receipt.json",
            "exact target-free carrier geometry, normalized I6, and global I6 extrema",
            "receipt_sha256",
        ),
        (
            "code/a5_fingerprint/runtime/seam_current_edge_prediction_receipt.json",
            "frozen FZ-12 p-basis coefficient ray",
            "receipt_sha256",
        ),
        (
            "code/a5_fingerprint/runtime/fz12_custody_projection.json",
            "data-free FZ-12 freeze and custody projection",
            "projection_sha256",
        ),
    ]
    result = []
    for (raw, value), (path, role, self_field) in zip(parents, metadata, strict=True):
        result.append(
            {
                "path": path,
                "role": role,
                "bytes": len(raw),
                "sha256": "sha256:" + sha256(raw),
                "schema": value["schema"],
                "status": value["status"],
                "self_digest": value[self_field],
            }
        )
    result.append(
        {
            "path": "Lean/Screen/SeamCurrentEdge30Remainder.lean",
            "role": "kernel-checked exact eighth seam moment",
            "bytes": len(lean_raw),
            "sha256": "sha256:" + sha256(lean_raw),
            "theorem": "OPH.SeamCurrentEdge30Remainder.seamMoment8_eq",
            "sorry_free": True,
        }
    )
    return result


def expected_geometry_contract() -> dict[str, Any]:
    return {
        "source_direction_count": 30,
        "unoriented_axis_count": 15,
        "raw_carrier_difference_norm_squared": "4",
        "unit_direction_definition": "w_e = carrierSeamDifference(e) / 2",
        "unit_direction_norm_squared": "1",
        "direction_multiplicity": (
            "one unit signed direction for each of the thirty canonical seam rows"
        ),
        "equal_weight": "1/5",
        "symbol_normalization": "lambda_hat = a^2 Lambda_a",
        "q_definition": "q = a k",
        "direction_definition": "n is a unit direction on S^2",
        "q_domain": {"minimum": "0", "maximum": "1", "inclusive": True},
        "pointwise_bound": "abs(w_e . n) <= 1",
        "cosine_even_orientation_boundary": True,
    }


def expected_exact_moments() -> dict[str, Any]:
    return {
        "definition": "M_r(n) = sum_e (w_e . n)^r",
        "M2": "10",
        "M4": "6",
        "M6": "30/7 - (2/7) I6(n)",
        "M8": "10/3 - (8/15) I6(n)",
        "M8_general_carrier": (
            "M8(k) = (10/3) |k|^8 - (8/15) |k|^2 I6(k)"
        ),
        "I6_normalization": "I6 = 1 on every unit primitive-port direction",
        "I6_unit_sphere_range": {"minimum": "-5/9", "maximum": "1"},
    }


def expected_exact_symbol() -> dict[str, Any]:
    return {
        "full_symbol": (
            "lambda_hat(q,n) = (1/5) sum_e [1 - cos(q (w_e . n))]"
        ),
        "P6": (
            "q^2 - q^4/20 + (1/840 - I6(n)/12600) q^6"
        ),
        "P8": (
            "P6 + (-1/60480 + I6(n)/378000) q^8"
        ),
        "P6_coefficients": {
            "q2": "1",
            "q4": "-1/20",
            "q6_isotropic": "1/840",
            "q6_I6": "-1/12600",
        },
        "P8_increment_coefficients": {
            "q8_isotropic": "-1/60480",
            "q8_I6": "1/378000",
        },
    }


def expected_taylor_remainders() -> dict[str, Any]:
    return {
        "validity_domain": "0 <= q <= 1",
        "R6": {
            "definition": "R6 = lambda_hat - P6",
            "signed_lower": (
                "-q^8 (10/3 - (8/15) I6(n)) / (5*8!)"
            ),
            "signed_upper": "0",
            "global_absolute_upper": "(7/388800) q^8",
            "global_absolute_coefficient": "7/388800",
            "self_contained_absolute_upper": "q^8/6720",
        },
        "R8": {
            "definition": "R8 = lambda_hat - P8",
            "signed_lower": "0",
            "pointwise_upper": "q^10 M10(n) / (5*10!)",
            "moment_domination": "M10(n) <= M8(n)",
            "global_upper": "(7/34992000) q^10",
            "global_coefficient": "7/34992000",
        },
    }


def expected_positivity() -> dict[str, Any]:
    return {
        "lower_bound": "(19/20) q^2",
        "upper_bound": "q^2",
        "strict_positive_for": "0 < q <= 1",
        "physical_monotonicity_claimed": False,
        "group_velocity_claimed": False,
    }


def expected_exposure_boundary() -> dict[str, Any]:
    return {
        "comparison_inputs": [],
        "comparison_data_read": False,
        "public_measurement_read": False,
        "target_values_read": False,
        "comparison_permitted": False,
        "score_emitted": False,
        "evidence_claimed": False,
        "verdict_emitted": False,
    }


def expected_physical_boundary() -> dict[str, Any]:
    return {
        "physical_position_identified": False,
        "spatial_symbol_to_frequency_squared_proved": False,
        "physical_sector_identified": False,
        "photon_attachment_proved": False,
        "physical_clock_proved": False,
        "frame_and_boost_map_proved": False,
        "cofinal_gluing_proved": False,
        "finite_physical_scale_identified": False,
        "wave_packet_and_readout_proved": False,
        "interaction_kinematics_proved": False,
        "physical_group_velocity_claimed": False,
    }


def expected_negative_controls() -> dict[str, Any]:
    rows = [
        ("double_directed_labels_60", "thirty-direction multiplicity"),
        ("projective_axes_15", "signed-direction multiplicity"),
        ("raw_norm2_differences", "unit-direction normalization"),
        ("fz11_rank_six_sign", "negative edge-orbit rank-six sign"),
        ("mutated_M8", "exact eighth moment"),
        ("q_above_one", "certified q domain"),
        ("unequal_weights", "equal weight 1/5"),
        ("extra_hop_shell", "complete thirty-direction support"),
    ]
    return {
        "all_detectors_fired": True,
        "controls": [
            {
                "control": control,
                "detector_fired": True,
                "expected_failure": expected_failure,
            }
            for control, expected_failure in rows
        ],
    }


def verify_receipt(
    path: Path,
    parents: list[tuple[bytes, dict[str, Any]]],
    lean_raw: bytes,
) -> dict[str, Any]:
    raw = path.read_bytes()
    try:
        receipt = json.loads(raw)
    except json.JSONDecodeError as error:
        raise VerificationError("invalid full-symbol receipt JSON") from error
    check(isinstance(receipt, dict), "full-symbol receipt root is not an object")
    check(raw == canonical_json_bytes(receipt), "noncanonical full-symbol receipt")
    check(
        set(receipt)
        == {
            "schema",
            "status",
            "issue",
            "parent_pins",
            "geometry_contract",
            "exact_moments",
            "exact_symbol",
            "taylor_remainders",
            "positivity_and_monotonicity",
            "exposure_boundary",
            "physical_boundary",
            "negative_controls",
            "receipt_sha256",
        },
        "full-symbol top-level keys drift",
    )
    check(receipt.get("schema") == SCHEMA, "full-symbol schema drift")
    check(receipt.get("status") == STATUS, "full-symbol status drift")
    check(receipt.get("issue") == 667, "full-symbol issue drift")
    check(
        receipt.get("receipt_sha256") == self_digest(receipt, "receipt_sha256"),
        "full-symbol self-digest drift",
    )
    check(
        receipt.get("parent_pins") == expected_parent_pins(parents, lean_raw),
        "full-symbol parent pins drift",
    )

    independent_geometry_checks()
    check(
        receipt.get("geometry_contract") == expected_geometry_contract(),
        "full-symbol geometry contract drift",
    )
    check(
        receipt.get("exact_moments") == expected_exact_moments(),
        "full-symbol exact moments drift",
    )
    check(
        receipt.get("exact_symbol") == expected_exact_symbol(),
        "full-symbol exact symbol drift",
    )
    check(
        receipt.get("taylor_remainders") == expected_taylor_remainders(),
        "full-symbol Taylor remainder drift",
    )
    check(
        receipt.get("positivity_and_monotonicity") == expected_positivity(),
        "full-symbol positivity boundary drift",
    )
    check(
        receipt.get("exposure_boundary") == expected_exposure_boundary(),
        "full-symbol exposure boundary drift",
    )
    check(
        receipt.get("physical_boundary") == expected_physical_boundary(),
        "full-symbol physical boundary drift",
    )
    check(
        receipt.get("negative_controls") == expected_negative_controls(),
        "full-symbol negative controls drift",
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-geometry", type=Path, default=DEFAULT_BASE_GEOMETRY)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--custody", type=Path, default=DEFAULT_CUSTODY)
    parser.add_argument("--remainder-lean", type=Path, default=DEFAULT_REMAINDER_LEAN)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        base_geometry = load_fixed_parent(
            args.base_geometry,
            byte_count=BASE_GEOMETRY_BYTES,
            digest=BASE_GEOMETRY_SHA256,
            schema=BASE_GEOMETRY_SCHEMA,
            status=BASE_GEOMETRY_STATUS,
            self_field="receipt_sha256",
        )
        verify_base_range_evidence(base_geometry[1])
        source = load_fixed_parent(
            args.source,
            byte_count=SOURCE_BYTES,
            digest=SOURCE_SHA256,
            schema=SOURCE_SCHEMA,
            status=SOURCE_STATUS,
            self_field="receipt_sha256",
        )
        custody = load_fixed_parent(
            args.custody,
            byte_count=CUSTODY_BYTES,
            digest=CUSTODY_SHA256,
            schema=CUSTODY_SCHEMA,
            status=CUSTODY_STATUS,
            self_field="projection_sha256",
        )
        lean_raw = load_fixed_lean_parent(args.remainder_lean)
        verify_receipt(args.receipt, [base_geometry, source, custody], lean_raw)
    except (OSError, VerificationError) as error:
        print(f"FZ12_FULL_SYMBOL_REMAINDER_VERIFY_FAIL: {error}", file=sys.stderr)
        return 1
    print("FZ12_FULL_SYMBOL_REMAINDER_VERIFY_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
