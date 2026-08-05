#!/usr/bin/env python3
"""Exact target-free frequency-speed certificate for positive cosine carriers.

The certificate binds the generic tight-frame contraction proved in
``Lean/Screen/CarrierFrequencySpeed.lean`` to the immutable FZ-11 vertex and
FZ-12 edge supports.  It reads no comparison data and makes no physical
frequency, clock, field, frame, scale, trajectory, or signal-front claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import a5_multipole_fixed_point_certificate as base


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
LEAN_PATH = ROOT / "Lean" / "Screen" / "CarrierFrequencySpeed.lean"
FZ11_PATH = RUNTIME / "spin_six_primitive_port_prediction_receipt.json"
FZ12_PATH = RUNTIME / "seam_current_edge_prediction_receipt.json"
RECEIPT_PATH = RUNTIME / "carrier_frequency_speed_receipt.json"

SCHEMA = "oph.carrier_frequency_speed.v1"
STATUS = (
    "EXACT_POSITIVE_TIGHT_FRAME_FREQUENCY_CONTRACTION__"
    "FZ11_FZ12_INSTANTIATED__PHYSICAL_BRIDGES_OPEN"
)

PARENTS = (
    {
        "path": FZ11_PATH,
        "row": "FZ-11",
        "schema": "oph.spin_six_primitive_port_prediction.v1",
        "status": (
            "FROZEN_PROSPECTIVE_PRIMITIVE_TWELVE_PORT_BRANCH_PREDICTION__"
            "PHYSICAL_COMPARISON_UNARMED"
        ),
        "bytes": 4809,
        "sha256": "8ac97d7c46199717ed031610efdda65c40f6a251e78715d6bc05888d598e66d8",
    },
    {
        "path": FZ12_PATH,
        "row": "FZ-12",
        "schema": "oph.seam_current_edge_prediction_candidate.v1",
        "status": (
            "EXACT_SOURCE_NATIVE_EDGE_RAY__PROSPECTIVE_PHYSICAL_BRANCH_UNARMED__"
            "PHYSICAL_PRODUCER_OPEN"
        ),
        "bytes": 9296,
        "sha256": "0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915",
    },
)


def canonical_json_bytes(payload: dict[str, Any]) -> bytes:
    return (
        json.dumps(
            payload,
            allow_nan=False,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        )
        + "\n"
    ).encode("ascii")


def tagged_sha256(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def reject_constant(value: str) -> None:
    raise RuntimeError(f"non-finite JSON constant rejected: {value}")


def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        require(key not in value, f"duplicate JSON key rejected: {key}")
        value[key] = item
    return value


def load_strict_json(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            raw.decode("ascii"),
            object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RuntimeError(f"invalid canonical JSON: {label}") from error
    require(isinstance(value, dict), f"JSON root is not an object: {label}")
    return value


def self_digest(payload: dict[str, Any], field: str) -> str:
    body = {key: value for key, value in payload.items() if key != field}
    return tagged_sha256(canonical_json_bytes(body))


def load_parent(contract: dict[str, Any]) -> tuple[bytes, dict[str, Any]]:
    raw = contract["path"].read_bytes()
    require(len(raw) == contract["bytes"], f"{contract['row']} byte drift")
    require(
        hashlib.sha256(raw).hexdigest() == contract["sha256"],
        f"{contract['row']} hash drift",
    )
    payload = load_strict_json(raw, contract["row"])
    require(raw == canonical_json_bytes(payload), f"{contract['row']} noncanonical")
    require(
        payload.get("schema") == contract["schema"], f"{contract['row']} schema drift"
    )
    require(
        payload.get("status") == contract["status"], f"{contract['row']} status drift"
    )
    require(
        payload.get("receipt_sha256") == self_digest(payload, "receipt_sha256"),
        f"{contract['row']} self-digest drift",
    )
    return raw, payload


def q5_outer_sum(
    vectors: list[tuple[base.Q5, base.Q5, base.Q5]],
) -> list[list[base.Q5]]:
    return [
        [sum_q5(base.q5_mul(v[i], v[j]) for v in vectors) for j in range(3)]
        for i in range(3)
    ]


def sum_q5(values) -> base.Q5:
    total = base.ZERO
    for value in values:
        total = base.q5_add(total, value)
    return total


def dot_q5(left, right) -> base.Q5:
    return sum_q5(base.q5_mul(x, y) for x, y in zip(left, right))


def sub_q5(left, right):
    return tuple(base.q5_sub(x, y) for x, y in zip(left, right))


def scale_q5(scale: Fraction, vector):
    return tuple(base.q5_scale(x, scale) for x in vector)


def support_certificate() -> dict[str, Any]:
    vertices = base.cartesian_vertices()
    require(len(vertices) == 12, "vertex count drift")
    radius_squared = dot_q5(vertices[0], vertices[0])
    require(
        all(dot_q5(v, v) == radius_squared for v in vertices), "vertex radius drift"
    )
    vertex_outer = q5_outer_sum(vertices)
    for i in range(3):
        for j in range(3):
            expected = base.q5_scale(radius_squared, 4) if i == j else base.ZERO
            require(vertex_outer[i][j] == expected, "vertex tight-frame drift")
    inverse_radius_squared = base.q5_div(base.ONE, radius_squared)
    seams = []
    for i in range(12):
        for j in range(i + 1, 12):
            normalized_dot = base.q5_mul(
                dot_q5(vertices[i], vertices[j]), inverse_radius_squared
            )
            if normalized_dot == base.INV_SQRT5:
                seams.append(scale_q5(Fraction(1, 2), sub_q5(vertices[j], vertices[i])))
    require(len(seams) == 30, "edge support count drift")
    require(all(dot_q5(v, v) == base.ONE for v in seams), "edge unit norm drift")
    edge_outer = q5_outer_sum(seams)
    for i in range(3):
        for j in range(3):
            expected = base.q5(10) if i == j else base.ZERO
            require(edge_outer[i][j] == expected, "edge tight-frame drift")
    return {
        "vertex12": {
            "support_count": 12,
            "raw_radius_squared": base.q5_str(radius_squared),
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


def load_lean() -> bytes:
    raw = LEAN_PATH.read_bytes()
    text = raw.decode("utf-8")
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
    )
    for fragment in required:
        require(fragment in text, f"Lean claim surface missing: {fragment}")
    require("sorry" not in text and "admit" not in text, "Lean placeholder detected")
    return raw


def build_receipt() -> dict[str, Any]:
    loaded = [load_parent(contract) for contract in PARENTS]
    lean_raw = load_lean()
    supports = support_certificate()

    parent_pins = []
    for contract, (raw, payload) in zip(PARENTS, loaded):
        parent_pins.append(
            {
                "row": contract["row"],
                "path": contract["path"].relative_to(ROOT).as_posix(),
                "bytes": len(raw),
                "sha256": tagged_sha256(raw),
                "schema": payload["schema"],
            }
        )

    body: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 704,
        "generic_exact_theorem": {
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
        },
        "exact_support_instantiations": supports,
        "branch_bindings": {
            "FZ-11": "the vertex12 tight frame normalizes to the immutable primitive-port symbol",
            "FZ-12": "the edge30 tight frame normalizes to the immutable seam-current character symbol",
            "frozen_bytes_modified": False,
            "new_prediction_payload": False,
        },
        "falsification_value": {
            "conditional_kill": "after physical position, frequency, clock, field-sector, wave-packet, carrier-frame, finite-scale, and detector/readout attachments, a resolved intrinsic superluminal group-speed component excludes the complete declared positive-weight cosine branch without a Taylor-remainder loophole",
            "null_rule": "a luminal null remains inconclusive while a can approach zero; this packet supplies no positive physical lower bound",
            "support_rule": "a subluminal signal is not unique to OPH and requires the frozen linked coefficient and angular tests for support",
        },
        "physical_boundary": {
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
        },
        "exposure_boundary": {
            "public_measurement_read": False,
            "comparison_data_read": False,
            "comparison_inputs": [],
            "score_emitted": False,
            "verdict_emitted": False,
        },
        "lean": {
            "path": LEAN_PATH.relative_to(ROOT).as_posix(),
            "bytes": len(lean_raw),
            "sha256": tagged_sha256(lean_raw),
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
        },
        "parent_pins": parent_pins,
        "contract_mutation_controls": {
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
        },
    }
    body["receipt_sha256"] = tagged_sha256(canonical_json_bytes(body))
    return body


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = build_receipt()
    raw = canonical_json_bytes(receipt)
    if args.write:
        RECEIPT_PATH.write_bytes(raw)
    else:
        print(raw.decode("ascii"), end="")


if __name__ == "__main__":
    main()
