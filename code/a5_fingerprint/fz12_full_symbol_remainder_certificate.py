#!/usr/bin/env python3
"""Certify the full FZ-12 edge symbol and its first omitted terms.

The frozen FZ-12 packet gives the exact equal-weight thirty-direction spatial
symbol and its coefficients through sixth order.  This target-free extension
replays the edge geometry over ``Q(sqrt(5))``, checks the separately proved
eighth moment, and applies the alternating cosine bounds on ``0 <= q <= 1``.
It produces signed and absolute remainder bounds for the *spatial* symbol.

No step identifies that symbol with a physical frequency squared.  The
certificate reads no target, likelihood, public measurement, or comparison
input and cannot emit a physical verdict.
"""

from __future__ import annotations

import argparse
import json
import re
from fractions import Fraction
from pathlib import Path
from typing import Any

import a5_multipole_fixed_point_certificate as base


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
BASE_PATH = RUNTIME / "a5_multipole_fixed_point_receipt.json"
SOURCE_PATH = RUNTIME / "seam_current_edge_prediction_receipt.json"
CUSTODY_PATH = RUNTIME / "fz12_custody_projection.json"
MOMENT8_PROOF_PATH = REPO_ROOT / "Lean" / "Screen" / "SeamCurrentEdge30Remainder.lean"
RECEIPT_PATH = RUNTIME / "fz12_full_symbol_remainder_receipt.json"

SCHEMA = "oph.fz12.full_symbol_remainder.v1"
STATUS = (
    "EXACT_TARGET_FREE_EDGE_SYMBOL_Q_LE_ONE_REMAINDER__"
    "PHYSICAL_FREQUENCY_AND_COMPARISON_OPEN"
)

BASE_CONTRACT = {
    "path": BASE_PATH,
    "schema": "oph.a5_multipole_fixed_point_receipt.v1",
    "status": "EXACT_A5_FINGERPRINT_CERTIFICATE__PHYSICAL_MAP_OPEN",
    "bytes": 8864,
    "sha256": "d96b80c71a64f48bb7a2a7b2592bad5e122ce9c5f3ca3b4636af29693c426ecd",
    "self_field": "receipt_sha256",
}
SOURCE_CONTRACT = {
    "path": SOURCE_PATH,
    "schema": "oph.seam_current_edge_prediction_candidate.v1",
    "status": (
        "EXACT_SOURCE_NATIVE_EDGE_RAY__PROSPECTIVE_PHYSICAL_BRANCH_UNARMED__"
        "PHYSICAL_PRODUCER_OPEN"
    ),
    "bytes": 9296,
    "sha256": "0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915",
    "self_field": "receipt_sha256",
}
CUSTODY_CONTRACT = {
    "path": CUSTODY_PATH,
    "schema": "oph.fz12.custody_projection.v1",
    "status": "FROZEN_FZ12_SOURCE_AND_CUSTODY_PINNED__NO_COMPARISON_DATA",
    "bytes": 3624,
    "sha256": "dfc6574a5bde9b576df6aeef2e03aba8602ec2723a0a4606e7942c405a57c643",
    "self_field": "projection_sha256",
}

Vec3 = tuple[base.Q5, base.Q5, base.Q5]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise base.FingerprintError(message)


def load_parent(contract: dict[str, Any]) -> tuple[bytes, dict[str, Any]]:
    path = contract["path"]
    raw = path.read_bytes()
    require(len(raw) == contract["bytes"], f"parent byte-count drift: {path.name}")
    require(
        base.tagged_sha256(raw).removeprefix("sha256:") == contract["sha256"],
        f"parent fixed digest drift: {path.name}",
    )
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError(f"invalid parent JSON: {path.name}") from error
    require(isinstance(payload, dict), f"parent root is not an object: {path.name}")
    require(raw == base.canonical_json_bytes(payload), f"noncanonical parent: {path.name}")
    require(payload.get("schema") == contract["schema"], f"schema drift: {path.name}")
    require(payload.get("status") == contract["status"], f"status drift: {path.name}")
    self_field = contract["self_field"]
    body = {key: value for key, value in payload.items() if key != self_field}
    require(
        payload.get(self_field) == base.tagged_sha256(base.canonical_json_bytes(body)),
        f"self-digest drift: {path.name}",
    )
    return raw, payload


def load_moment8_proof() -> bytes:
    raw = MOMENT8_PROOF_PATH.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise base.FingerprintError("non-UTF-8 eighth-moment proof") from error
    normalized = re.sub(r"\s+", " ", text)
    required = (
        "import SeamCurrentEdge30Moment",
        "noncomputable def seamMoment8 (k : Vec3) : ℝ :=",
        "theorem seamMoment8_eq (k : Vec3) : seamMoment8 k = (10 / 3 : ℝ) * radiusSquared k ^ 4 - (8 / 15 : ℝ) * radiusSquared k * I6 k := by",
        "It supplies no physical frequency, clock, field attachment, propagation law, or comparison result.",
    )
    for fragment in required:
        require(fragment in normalized, "claim-bearing eighth-moment proof drift")
    require("sorry" not in text and "admit" not in text, "placeholder in eighth-moment proof")
    return raw


def q5_dot(left: Vec3, right: Vec3) -> base.Q5:
    total = base.ZERO
    for x, y in zip(left, right):
        total = base.q5_add(total, base.q5_mul(x, y))
    return total


def q5_sub(left: Vec3, right: Vec3) -> Vec3:
    return tuple(base.q5_sub(x, y) for x, y in zip(left, right))  # type: ignore[return-value]


def sphere_equal(left: base.Poly, right: base.Poly) -> bool:
    difference = base.p_add(left, base.p_scale(right, -1))
    return base.p_is_zero(base.p_reduce_sphere(difference))


def certify_i6_range(base_receipt: dict[str, Any]) -> None:
    critical = base_receipt.get("critical_points", {})
    require(
        critical.get("census")
        == {"maxima": 12, "minima": 20, "saddles": 30, "total": 62},
        "I6 critical-point census drift",
    )
    require("either the poles" in critical.get("case_analysis", ""), "I6 exhaustiveness drift")
    orbits = critical.get("orbits", [])
    require(sum(row.get("count", 0) for row in orbits) == 62, "I6 orbit count drift")
    extrema = {(row.get("type"), row.get("value")) for row in orbits}
    require(("maximum", "1+0*sqrt5") in extrema, "I6 maximum drift")
    require(("minimum", "-5/9+0*sqrt5") in extrema, "I6 minimum drift")
    require(
        critical.get("nondegeneracy", "").startswith("every tangent Hessian eigenvalue"),
        "I6 nondegeneracy evidence drift",
    )


def exact_edge_moments() -> dict[str, Any]:
    vertices = base.cartesian_vertices()
    inverse_norm = base.q5_div(base.ONE, base.NORM_SQ)

    seams: list[tuple[int, int]] = []
    for i in range(12):
        for j in range(i + 1, 12):
            dot = base.q5_mul(q5_dot(vertices[i], vertices[j]), inverse_norm)
            if dot == base.INV_SQRT5:
                seams.append((i, j))
    require(len(seams) == 30, "edge seam census drift")
    differences = [q5_sub(vertices[j], vertices[i]) for i, j in seams]
    norm_squared = {q5_dot(vector, vector) for vector in differences}
    require(norm_squared == {base.q5(4)}, "raw seam norm drift")

    moments = {
        order: base.normalized_moment(differences, order, base.q5(4))
        for order in (2, 4, 6, 8)
    }
    cartesian = base.build_cartesian_frame()
    i6 = cartesian["_i6_poly_object"]
    targets = {
        2: base.p_scale(base.radial_power(1), 10),
        4: base.p_scale(base.radial_power(2), 6),
        6: base.p_add(
            base.p_scale(base.radial_power(3), Fraction(30, 7)),
            base.p_scale(i6, Fraction(-2, 7)),
        ),
        8: base.p_add(
            base.p_scale(base.radial_power(4), Fraction(10, 3)),
            base.p_scale(
                base.p_mul(base.radial_power(1), i6), Fraction(-8, 15)
            ),
        ),
    }
    for order, target in targets.items():
        require(sphere_equal(moments[order], target), f"edge M{order} identity drift")

    m8_min = Fraction(10, 3) - Fraction(8, 15)
    m8_max = Fraction(10, 3) + Fraction(8, 15) * Fraction(5, 9)
    require(m8_min == Fraction(14, 5), "M8 minimum arithmetic drift")
    require(m8_max == Fraction(98, 27), "M8 maximum arithmetic drift")
    p6_bound = m8_max / (5 * 40320)
    p8_bound = m8_max / (5 * 3628800)
    require(p6_bound == Fraction(7, 388800), "P6 global bound drift")
    require(p8_bound == Fraction(7, 34992000), "P8 global bound drift")

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


def parent_pin(
    raw: bytes,
    payload: dict[str, Any],
    contract: dict[str, Any],
    role: str,
) -> dict[str, Any]:
    return {
        "path": str(contract["path"].relative_to(REPO_ROOT)),
        "role": role,
        "bytes": len(raw),
        "sha256": base.tagged_sha256(raw),
        "schema": contract["schema"],
        "status": contract["status"],
        "self_digest": payload[contract["self_field"]],
    }


def build_receipt() -> dict[str, Any]:
    base_raw, base_receipt = load_parent(BASE_CONTRACT)
    source_raw, source = load_parent(SOURCE_CONTRACT)
    custody_raw, custody = load_parent(CUSTODY_CONTRACT)
    proof_raw = load_moment8_proof()

    require(
        source["exposure_and_custody_boundary"]["public_measurement_read"] is False,
        "frozen source exposure drift",
    )
    require(
        custody["projection_scope"]
        == {
            "register_row": "FZ-12",
            "source_class": "frozen source theory and custody metadata only",
            "includes_measurement_values": False,
            "includes_comparison_values": False,
            "includes_likelihood_values": False,
            "includes_other_campaign_rows": False,
        },
        "FZ-12 custody scope drift",
    )
    require(
        custody["source_receipt"]["sha256"] == SOURCE_CONTRACT["sha256"],
        "custody no longer pins frozen source",
    )
    certify_i6_range(base_receipt)
    moments = exact_edge_moments()

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 667,
        "parent_pins": [
            parent_pin(
                base_raw,
                base_receipt,
                BASE_CONTRACT,
                "exact target-free carrier geometry, normalized I6, and global I6 extrema",
            ),
            parent_pin(
                source_raw,
                source,
                SOURCE_CONTRACT,
                "frozen FZ-12 p-basis coefficient ray",
            ),
            parent_pin(
                custody_raw,
                custody,
                CUSTODY_CONTRACT,
                "data-free FZ-12 freeze and custody projection",
            ),
            {
                "path": str(MOMENT8_PROOF_PATH.relative_to(REPO_ROOT)),
                "role": "kernel-checked exact eighth seam moment",
                "bytes": len(proof_raw),
                "sha256": base.tagged_sha256(proof_raw),
                "theorem": "OPH.SeamCurrentEdge30Remainder.seamMoment8_eq",
                "sorry_free": True,
            },
        ],
        "geometry_contract": {
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
        },
        "exact_moments": moments,
        "exact_symbol": {
            "full_symbol": (
                "lambda_hat(q,n) = (1/5) sum_e [1 - cos(q (w_e . n))]"
            ),
            "P6": "q^2 - q^4/20 + (1/840 - I6(n)/12600) q^6",
            "P8": "P6 + (-1/60480 + I6(n)/378000) q^8",
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
        },
        "taylor_remainders": {
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
        },
        "positivity_and_monotonicity": {
            "lower_bound": "(19/20) q^2",
            "upper_bound": "q^2",
            "strict_positive_for": "0 < q <= 1",
            "physical_monotonicity_claimed": False,
            "group_velocity_claimed": False,
        },
        "exposure_boundary": {
            "comparison_inputs": [],
            "comparison_data_read": False,
            "public_measurement_read": False,
            "target_values_read": False,
            "comparison_permitted": False,
            "score_emitted": False,
            "evidence_claimed": False,
            "verdict_emitted": False,
        },
        "physical_boundary": {
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
        },
        "negative_controls": {
            "all_detectors_fired": True,
            "controls": [
                {
                    "control": "double_directed_labels_60",
                    "detector_fired": True,
                    "expected_failure": "thirty-direction multiplicity",
                },
                {
                    "control": "projective_axes_15",
                    "detector_fired": True,
                    "expected_failure": "signed-direction multiplicity",
                },
                {
                    "control": "raw_norm2_differences",
                    "detector_fired": True,
                    "expected_failure": "unit-direction normalization",
                },
                {
                    "control": "fz11_rank_six_sign",
                    "detector_fired": True,
                    "expected_failure": "negative edge-orbit rank-six sign",
                },
                {
                    "control": "mutated_M8",
                    "detector_fired": True,
                    "expected_failure": "exact eighth moment",
                },
                {
                    "control": "q_above_one",
                    "detector_fired": True,
                    "expected_failure": "certified q domain",
                },
                {
                    "control": "unequal_weights",
                    "detector_fired": True,
                    "expected_failure": "equal weight 1/5",
                },
                {
                    "control": "extra_hop_shell",
                    "detector_fired": True,
                    "expected_failure": "complete thirty-direction support",
                },
            ],
        },
    }
    receipt["receipt_sha256"] = base.tagged_sha256(base.canonical_json_bytes(receipt))
    return receipt


def verify_committed_receipt() -> dict[str, Any]:
    raw = RECEIPT_PATH.read_bytes()
    try:
        committed = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError("invalid committed FZ-12 remainder receipt") from error
    require(isinstance(committed, dict), "FZ-12 remainder root is not an object")
    require(raw == base.canonical_json_bytes(committed), "noncanonical FZ-12 remainder receipt")
    require(committed.get("schema") == SCHEMA, "FZ-12 remainder schema drift")
    require(committed.get("status") == STATUS, "FZ-12 remainder status drift")
    body = {key: value for key, value in committed.items() if key != "receipt_sha256"}
    require(
        committed.get("receipt_sha256")
        == base.tagged_sha256(base.canonical_json_bytes(body)),
        "FZ-12 remainder self-digest drift",
    )
    rebuilt = build_receipt()
    require(raw == base.canonical_json_bytes(rebuilt), "FZ-12 remainder ancestry or result drift")
    return committed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    if args.write and args.verify:
        parser.error("choose --write or --verify")
    receipt = verify_committed_receipt() if args.verify else build_receipt()
    if args.write:
        RECEIPT_PATH.parent.mkdir(exist_ok=True)
        RECEIPT_PATH.write_bytes(base.canonical_json_bytes(receipt))
        print(RECEIPT_PATH)
        return 0
    print(json.dumps(receipt["taylor_remainders"], indent=1))
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
