#!/usr/bin/env python3
"""Independent verifier for the exact carrier frequency-speed packet.

This verifier does not import the producer.  It reconstructs the twelve
primitive vertices and thirty seam directions over ``Q(sqrt(5))``, checks the
exact tight-frame identities, pins both frozen source receipts, and
audits the kernel-checked Lean theorem surface.

The verified frequency is the norm of an auxiliary sine feature map.  The
packet does not identify it with a physical frequency, particle sector,
clock, wave packet, signal front, or laboratory readout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
DEFAULT_FZ11 = RUNTIME / "spin_six_primitive_port_prediction_receipt.json"
DEFAULT_FZ12 = RUNTIME / "seam_current_edge_prediction_receipt.json"
DEFAULT_LEAN = ROOT / "Lean" / "Screen" / "CarrierFrequencySpeed.lean"
DEFAULT_RECEIPT = RUNTIME / "carrier_frequency_speed_receipt.json"

SCHEMA = "oph.carrier_frequency_speed.v1"
STATUS = (
    "EXACT_POSITIVE_TIGHT_FRAME_FREQUENCY_CONTRACTION__"
    "FZ11_FZ12_INSTANTIATED__PHYSICAL_BRIDGES_OPEN"
)

FZ11_BYTES = 4809
FZ11_SHA256 = "8ac97d7c46199717ed031610efdda65c40f6a251e78715d6bc05888d598e66d8"
FZ11_SCHEMA = "oph.spin_six_primitive_port_prediction.v1"
FZ11_STATUS = (
    "FROZEN_PROSPECTIVE_PRIMITIVE_TWELVE_PORT_BRANCH_PREDICTION__"
    "PHYSICAL_COMPARISON_UNARMED"
)

FZ12_BYTES = 9296
FZ12_SHA256 = "0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915"
FZ12_SCHEMA = "oph.seam_current_edge_prediction_candidate.v1"
FZ12_STATUS = (
    "EXACT_SOURCE_NATIVE_EDGE_RAY__PROSPECTIVE_PHYSICAL_BRANCH_UNARMED__"
    "PHYSICAL_PRODUCER_OPEN"
)

LEAN_BYTES = 11315
LEAN_SHA256 = "7403545b04045b68091492e5524943cd8dd56f5ff2176e4236878f395e0c0af9"


class VerificationError(ValueError):
    """A file or semantic field failed the independent contract."""


def check(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def reject_constant(value: str) -> None:
    raise VerificationError(f"non-finite JSON constant rejected: {value}")


def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        check(key not in value, f"duplicate JSON key rejected: {key}")
        value[key] = item
    return value


def strict_json_loads(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            raw.decode("ascii"),
            object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise VerificationError(f"invalid canonical JSON: {label}") from error
    check(isinstance(value, dict), f"JSON root is not an object: {label}")
    return value


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
    label: str,
    byte_count: int,
    digest: str,
    schema: str,
    status: str,
) -> tuple[bytes, dict[str, Any]]:
    raw = path.read_bytes()
    value = strict_json_loads(raw, label)
    check(raw == canonical_json_bytes(value), f"noncanonical parent: {label}")
    check(len(raw) == byte_count, f"parent byte-count drift: {label}")
    check(sha256(raw) == digest, f"parent raw hash drift: {label}")
    check(value.get("schema") == schema, f"parent schema drift: {label}")
    check(value.get("status") == status, f"parent status drift: {label}")
    check(
        value.get("receipt_sha256") == self_digest(value, "receipt_sha256"),
        f"parent self-digest drift: {label}",
    )
    return raw, value


# A scalar is represented as a + b sqrt(5), with rational a and b.
Q5 = tuple[Fraction, Fraction]
Vec3 = tuple[Q5, Q5, Q5]

ZERO: Q5 = (Fraction(0), Fraction(0))
ONE: Q5 = (Fraction(1), Fraction(0))
PHI: Q5 = (Fraction(1, 2), Fraction(1, 2))
INV_SQRT5: Q5 = (Fraction(0), Fraction(1, 5))


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


def q5_scale(value: Q5, scalar: int | Fraction) -> Q5:
    factor = Fraction(scalar)
    return (factor * value[0], factor * value[1])


def q5_inverse(value: Q5) -> Q5:
    denominator = value[0] * value[0] - 5 * value[1] * value[1]
    check(denominator != 0, "division by zero in Q(sqrt(5))")
    return (value[0] / denominator, -value[1] / denominator)


def sum_q5(values: Iterable[Q5]) -> Q5:
    total = ZERO
    for value in values:
        total = q5_add(total, value)
    return total


def dot(left: Vec3, right: Vec3) -> Q5:
    return sum_q5(q5_mul(a, b) for a, b in zip(left, right, strict=True))


def vec_sub(left: Vec3, right: Vec3) -> Vec3:
    return tuple(q5_sub(a, b) for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def vec_scale(value: Vec3, scalar: int | Fraction) -> Vec3:
    return tuple(q5_scale(item, scalar) for item in value)  # type: ignore[return-value]


def cartesian_vertices() -> list[Vec3]:
    vertices: list[Vec3] = []
    for first in (-1, 1):
        for second in (-1, 1):
            scalar_first = q5(first)
            phi_second = q5_scale(PHI, second)
            vertices.extend(
                [
                    (ZERO, scalar_first, phi_second),
                    (scalar_first, phi_second, ZERO),
                    (phi_second, ZERO, scalar_first),
                ]
            )
    check(len(set(vertices)) == 12, "primitive vertex uniqueness drift")
    return vertices


def outer_sum(vectors: list[Vec3]) -> list[list[Q5]]:
    return [
        [sum_q5(q5_mul(vector[i], vector[j]) for vector in vectors) for j in range(3)]
        for i in range(3)
    ]


def independent_support_checks() -> dict[str, Any]:
    vertices = cartesian_vertices()
    radius_squared = dot(vertices[0], vertices[0])
    check(radius_squared == (Fraction(5, 2), Fraction(1, 2)), "vertex radius drift")
    check(all(dot(v, v) == radius_squared for v in vertices), "unequal vertex radius")

    vertex_outer = outer_sum(vertices)
    for i in range(3):
        for j in range(3):
            expected = q5_scale(radius_squared, 4) if i == j else ZERO
            check(vertex_outer[i][j] == expected, "vertex second moment drift")

    inverse_radius_squared = q5_inverse(radius_squared)
    seams: list[Vec3] = []
    for i in range(12):
        for j in range(i + 1, 12):
            normalized_dot = q5_mul(
                dot(vertices[i], vertices[j]), inverse_radius_squared
            )
            if normalized_dot == INV_SQRT5:
                seams.append(
                    vec_scale(vec_sub(vertices[j], vertices[i]), Fraction(1, 2))
                )
    check(len(seams) == 30, "edge support count drift")
    check(all(dot(seam, seam) == ONE for seam in seams), "edge norm drift")

    edge_outer = outer_sum(seams)
    for i in range(3):
        for j in range(3):
            expected = q5(10) if i == j else ZERO
            check(edge_outer[i][j] == expected, "edge second moment drift")
    return {
        "vertex12": {
            "support_count": 12,
            "raw_radius_squared": "5/2+1/2*sqrt5",
            "unit_second_moment": "sum_i (u_i.dot.x)^2 = 4 |x|^2",
            "tight_constant": "4",
            "normalized_symbol_prefactor": "1/(2 a^2)",
        },
        "edge30": {
            "support_count": 30,
            "unit_norm_squared": "1",
            "unit_second_moment": "sum_e (w_e.dot.x)^2 = 10 |x|^2",
            "tight_constant": "10",
            "normalized_symbol_prefactor": "1/(5 a^2)",
        },
    }


def expected_generic_theorem() -> dict[str, Any]:
    return {
        "premises": [
            "a finite support with strictly positive weights",
            "a positive tight constant t",
            "sum_i rho_i (v_i dot x)^2 = t |x|^2 for every x",
            "the full normalized symbol is (2/(t a^2)) sum_i rho_i [1-cos(a v_i dot k)] with a != 0",
        ],
        "feature_map": "Phi_i(k)=sqrt(4 rho_i/(a^2 t)) sin(a v_i dot k/2)",
        "feature_identity": "Lambda_a(k)=||Phi_a(k)||^2",
        "contraction": "||Phi_a(k)-Phi_a(p)|| <= |k-p|",
        "frequency_definition": "Omega_a(k)=||Phi_a(k)||=sqrt(Lambda_a(k))",
        "frequency_bound": "|Omega_a(k)-Omega_a(p)| <= |k-p|",
        "certified_upper_constant": "1",
        "scope": "all momenta; no Taylor truncation and no angular orientation choice",
    }


def expected_physical_boundary() -> dict[str, bool]:
    return {
        "physical_position_proved": False,
        "physical_frequency_proved": False,
        "physical_clock_proved": False,
        "photon_or_other_field_sector_proved": False,
        "wave_packet_or_signal_front_proved": False,
        "carrier_frame_and_boost_proved": False,
        "finite_physical_scale_proved": False,
        "positive_scale_lower_bound_proved": False,
        "source_lag_and_detector_readout_proved": False,
        "branch_exclusivity_proved": False,
        "comparison_permitted": False,
    }


def expected_exposure_boundary() -> dict[str, Any]:
    return {
        "public_measurement_read": False,
        "comparison_data_read": False,
        "comparison_inputs": [],
        "score_emitted": False,
        "verdict_emitted": False,
    }


def expected_branch_bindings() -> dict[str, Any]:
    return {
        "FZ-11": "the vertex12 tight frame normalizes to the immutable primitive-port symbol",
        "FZ-12": "the edge30 tight frame normalizes to the immutable seam-current character symbol",
        "frozen_bytes_modified": False,
        "new_prediction_payload": False,
    }


def expected_falsification_value() -> dict[str, str]:
    return {
        "conditional_kill": "after physical position, frequency, clock, field-sector, wave-packet, carrier-frame, finite-scale, and detector/readout attachments, a resolved intrinsic superluminal group-speed component excludes the complete declared positive-weight cosine branch without a Taylor-remainder loophole",
        "null_rule": "a luminal null remains inconclusive while a can approach zero; this packet supplies no positive physical lower bound",
        "support_rule": "a subluminal signal is not unique to OPH and requires the frozen linked coefficient and angular tests for support",
    }


def expected_contract_mutation_controls() -> dict[str, Any]:
    return {
        "controls": [
            "replace the positive-weight premise by signed weights",
            "change the edge tight constant",
            "change the vertex tight constant",
            "change the certified unit upper constant",
            "promote the auxiliary norm to physical frequency",
            "expose comparison data",
            "mutate a frozen parent pin",
        ],
        "all_rejected_by_test_suite": True,
    }


def verify_lean(path: Path) -> bytes:
    raw = path.read_bytes()
    check(len(raw) == LEAN_BYTES, "Lean byte-count drift")
    check(sha256(raw) == LEAN_SHA256, "Lean raw hash drift")
    try:
        source = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise VerificationError("Lean source is not UTF-8") from error
    required = (
        "structure PositiveTightFrame",
        "theorem feature_norm_sq_eq_cosineSymbol",
        "theorem feature_dist_le",
        "theorem frequency_global_one_lipschitz",
        "theorem unit_port_second_moment_eq",
        "theorem edge30_cosineSymbol_eq_fz12",
        "theorem vertex12_cosineSymbol_eq_fz11",
        "theorem fz12_frequency_global_one_lipschitz",
        "theorem fz11_frequency_global_one_lipschitz",
        "requires separate position, field, clock, wave-packet, frame, scale, and",
        "#print axioms OPH.CarrierFrequencySpeed.frequency_global_one_lipschitz",
    )
    for fragment in required:
        check(fragment in source, f"Lean theorem surface missing: {fragment}")
    check(
        re.search(r"\b(?:sorry|admit)\b", source) is None,
        "Lean source contains a placeholder",
    )
    return raw


def expected_parent_pins() -> list[dict[str, Any]]:
    return [
        {
            "row": "FZ-11",
            "path": "code/a5_fingerprint/runtime/spin_six_primitive_port_prediction_receipt.json",
            "bytes": FZ11_BYTES,
            "sha256": "sha256:" + FZ11_SHA256,
            "schema": FZ11_SCHEMA,
        },
        {
            "row": "FZ-12",
            "path": "code/a5_fingerprint/runtime/seam_current_edge_prediction_receipt.json",
            "bytes": FZ12_BYTES,
            "sha256": "sha256:" + FZ12_SHA256,
            "schema": FZ12_SCHEMA,
        },
    ]


def expected_lean_pin() -> dict[str, Any]:
    return {
        "path": "Lean/Screen/CarrierFrequencySpeed.lean",
        "bytes": LEAN_BYTES,
        "sha256": "sha256:" + LEAN_SHA256,
        "theorems": [
            "feature_norm_sq_eq_cosineSymbol",
            "feature_dist_le",
            "frequency_global_one_lipschitz",
            "edge30_cosineSymbol_eq_fz12",
            "fz12_frequency_global_one_lipschitz",
            "unit_port_second_moment_eq",
            "vertex12_cosineSymbol_eq_fz11",
            "fz11_frequency_global_one_lipschitz",
        ],
        "sorry_free": True,
    }


def verify_payload(
    payload: dict[str, Any],
    *,
    check_files: bool = True,
    fz11_path: Path = DEFAULT_FZ11,
    fz12_path: Path = DEFAULT_FZ12,
    lean_path: Path = DEFAULT_LEAN,
) -> None:
    expected_keys = {
        "schema",
        "status",
        "issue",
        "generic_exact_theorem",
        "exact_support_instantiations",
        "branch_bindings",
        "falsification_value",
        "physical_boundary",
        "exposure_boundary",
        "lean",
        "parent_pins",
        "contract_mutation_controls",
        "receipt_sha256",
    }
    check(set(payload) == expected_keys, "top-level keys drift")
    check(
        payload.get("receipt_sha256") == self_digest(payload, "receipt_sha256"),
        "receipt self-digest drift",
    )
    check(payload.get("schema") == SCHEMA, "receipt schema drift")
    check(payload.get("status") == STATUS, "receipt status drift")
    check(payload.get("issue") == 704, "receipt issue drift")

    check(
        payload.get("generic_exact_theorem") == expected_generic_theorem(),
        "generic theorem contract drift",
    )
    check(
        payload.get("exact_support_instantiations") == independent_support_checks(),
        "support instantiation drift",
    )
    check(
        payload.get("branch_bindings") == expected_branch_bindings(),
        "branch binding drift",
    )
    check(
        payload.get("falsification_value") == expected_falsification_value(),
        "falsification boundary drift",
    )
    check(
        payload.get("physical_boundary") == expected_physical_boundary(),
        "physical boundary drift",
    )
    check(
        payload.get("exposure_boundary") == expected_exposure_boundary(),
        "exposure boundary drift",
    )
    check(payload.get("parent_pins") == expected_parent_pins(), "parent pins drift")
    check(payload.get("lean") == expected_lean_pin(), "Lean pin drift")
    check(
        payload.get("contract_mutation_controls")
        == expected_contract_mutation_controls(),
        "contract-mutation registry drift",
    )

    if check_files:
        load_fixed_parent(
            fz11_path,
            label="FZ-11",
            byte_count=FZ11_BYTES,
            digest=FZ11_SHA256,
            schema=FZ11_SCHEMA,
            status=FZ11_STATUS,
        )
        load_fixed_parent(
            fz12_path,
            label="FZ-12",
            byte_count=FZ12_BYTES,
            digest=FZ12_SHA256,
            schema=FZ12_SCHEMA,
            status=FZ12_STATUS,
        )
        verify_lean(lean_path)


def verify_receipt(
    receipt_path: Path,
    *,
    fz11_path: Path = DEFAULT_FZ11,
    fz12_path: Path = DEFAULT_FZ12,
    lean_path: Path = DEFAULT_LEAN,
) -> dict[str, Any]:
    raw = receipt_path.read_bytes()
    payload = strict_json_loads(raw, "carrier frequency-speed receipt")
    check(raw == canonical_json_bytes(payload), "receipt is not canonical JSON")
    verify_payload(
        payload,
        check_files=True,
        fz11_path=fz11_path,
        fz12_path=fz12_path,
        lean_path=lean_path,
    )
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fz11", type=Path, default=DEFAULT_FZ11)
    parser.add_argument("--fz12", type=Path, default=DEFAULT_FZ12)
    parser.add_argument("--lean", type=Path, default=DEFAULT_LEAN)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args()
    try:
        verify_receipt(
            args.receipt,
            fz11_path=args.fz11,
            fz12_path=args.fz12,
            lean_path=args.lean,
        )
    except (OSError, VerificationError) as error:
        print(f"CARRIER_FREQUENCY_SPEED_INDEPENDENT_FAIL: {error}", file=sys.stderr)
        return 1
    print("CARRIER_FREQUENCY_SPEED_INDEPENDENT_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
