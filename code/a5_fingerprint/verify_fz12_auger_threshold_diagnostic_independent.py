#!/usr/bin/env python3
"""Independent verifier for the exposed FZ-12 Auger diagnostic."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
DEFAULT_SOURCE = RUNTIME / "seam_current_edge_prediction_receipt.json"
DEFAULT_CUSTODY = RUNTIME / "fz12_custody_projection.json"
DEFAULT_OBSERVATION = RUNTIME / "fz12_observation_map_receipt.json"
DEFAULT_RECEIPT = RUNTIME / "fz12_auger_threshold_diagnostic_receipt.json"

SOURCE_BYTES = 9296
SOURCE_SHA256 = "0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915"
CUSTODY_BYTES = 3624
CUSTODY_SHA256 = "dfc6574a5bde9b576df6aeef2e03aba8602ec2723a0a4606e7942c405a57c643"
OBSERVATION_BYTES = 5276
OBSERVATION_SHA256 = "0b34218b075a3d51fe0badb00e3bb889743ae19b6dcf161f7877683e25121d17"
RECEIPT_BYTES = 5659
RECEIPT_SHA256 = "1cfcd876c122edf7406c54b3caeda7d457500d8137f90b0142110913fad262dd"

SCHEMA = "oph.fz12.auger_threshold_diagnostic.v1"
STATUS = "EXPOSED_RETROSPECTIVE_FZ12_AUGER_THRESHOLD_BOUND__CONDITIONAL_DIAGNOSTIC_ONLY"


class VerificationError(ValueError):
    """An input failed the independent diagnostic contract."""


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


def independent_coefficients() -> dict[str, str]:
    c = Fraction(-1, 20)
    d_iso = Fraction(1, 840)
    d_i6 = Fraction(-1, 12600)
    delta2 = c
    delta4_iso = d_iso - 2 * c * c
    delta4_i6 = d_i6
    check(delta4_iso == Fraction(-2, 525), "independent isotropic basis map failed")
    return {
        "c_over_a2": str(c),
        "d_iso_over_a4": str(d_iso),
        "d_I6_over_a4": str(d_i6),
        "delta_gamma2_over_a2": str(delta2),
        "delta_gamma4_iso_over_a4": str(delta4_iso),
        "delta_gamma4_I6_over_a4": str(delta4_i6),
    }


def independent_decimals() -> dict[str, str]:
    with localcontext() as context:
        context.prec = 80
        a_eV_inverse = Decimal(20).sqrt() * Decimal("1e-29")
        a_m = a_eV_inverse * Decimal("1.973269804e-7")
        a_planck = a_m / Decimal("1.616255e-35")
    return {
        "a_upper_approx_eV_inverse": format(a_eV_inverse, ".9E").replace("E", "e"),
        "a_upper_approx_m": format(a_m, ".9E").replace("E", "e"),
        "a_upper_approx_planck_lengths": format(a_planck, ".10f"),
        "headline_m": format(a_m, ".4E").replace("E", "e"),
        "headline_planck_lengths": format(a_planck, ".4f"),
    }


def expected_parent_pins(
    parents: list[tuple[bytes, dict[str, Any]]],
) -> list[dict[str, Any]]:
    metadata = [
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
        (
            "code/a5_fingerprint/runtime/fz12_observation_map_receipt.json",
            "exact formal FZ-12 frequency and velocity map",
            "receipt_sha256",
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
    return result


def verify_receipt(
    path: Path, parents: list[tuple[bytes, dict[str, Any]]]
) -> dict[str, Any]:
    raw = path.read_bytes()
    try:
        receipt = json.loads(raw)
    except json.JSONDecodeError as error:
        raise VerificationError("invalid diagnostic receipt JSON") from error
    check(isinstance(receipt, dict), "diagnostic receipt root is not an object")
    check(raw == canonical_json_bytes(receipt), "noncanonical diagnostic receipt")
    check(sha256(raw) == RECEIPT_SHA256, "diagnostic receipt raw hash drift")
    check(len(raw) == RECEIPT_BYTES, "diagnostic receipt byte-count drift")
    check(
        set(receipt)
        == {
            "basis_translation",
            "conditional_bound",
            "diagnostic_id",
            "exposed_public_input",
            "issue",
            "open_physical_attachments",
            "parent_pins",
            "receipt_sha256",
            "schema",
            "scope_boundary",
            "status",
            "theory_input",
        },
        "diagnostic top-level keys drift",
    )
    check(receipt.get("schema") == SCHEMA, "diagnostic schema drift")
    check(receipt.get("status") == STATUS, "diagnostic status drift")
    check(receipt.get("issue") == 667, "diagnostic issue drift")
    check(
        receipt.get("diagnostic_id") == "FZ12-AUGER-THRESHOLD-EXPOSED-RETROSPECTIVE",
        "diagnostic identifier drift",
    )
    check(
        receipt.get("receipt_sha256") == self_digest(receipt, "receipt_sha256"),
        "diagnostic self-digest drift",
    )
    check(
        receipt.get("parent_pins") == expected_parent_pins(parents),
        "diagnostic parent pins drift",
    )

    check(
        receipt.get("theory_input")
        == {
            "basis": "momentum p basis",
            "dispersion": ("E^2 = p^2 + c p^4 + (d_iso + d_I6 I6) p^6 + O(p^8)"),
            "coefficients": {
                "c": "-a^2/20",
                "d_iso": "a^4/840",
                "d_I6": "-a^4/12600",
            },
            "physical_attachment_closed": False,
        },
        "theory input drift",
    )
    check(
        receipt.get("basis_translation")
        == {
            "auger_energy_basis": (
                "E^2 = p^2 + delta_gamma2 E^4 + "
                "(delta_gamma4_iso + delta_gamma4_I6 I6) E^6 + O(E^8)"
            ),
            "series_identity": (
                "E^4 = p^4 + 2 c p^6 + O(p^8), so delta_gamma2 = c and "
                "delta_gamma4 = d - 2 c^2; the isotropic c leaves the I6 "
                "coefficient unchanged"
            ),
            "exact_coefficients": independent_coefficients(),
            "paper_consumes_delta_gamma4": False,
            "mapping_scope": "leading EFT series through sixth order only",
            "higher_order_remainder_control": False,
        },
        "basis translation drift",
    )
    check(
        receipt.get("exposed_public_input")
        == {
            "exposure_class": "EXPOSED_RETROSPECTIVE",
            "source": (
                "Pierre Auger Collaboration, Testing effects of Lorentz invariance "
                "violation in the propagation of astroparticles with the Pierre "
                "Auger Observatory, JCAP 01 (2022) 023"
            ),
            "doi": "10.1088/1475-7516/2022/01/023",
            "arxiv": "2112.06773",
            "paper_parameterization": ("E_i^2 = m_i^2 + p_i^2 + delta_i,n E_i^(2+n)"),
            "selected_scenario": (
                "alternative source scenario with a subdominant proton component "
                "extending to 1e20 eV"
            ),
            "selected_scenario_maximum_proton_energy_eV": "1e20",
            "delta_gamma2_strict_lower_bound_eV_minus2": "-1e-58",
            "confidence_level_for_direct_electromagnetic_bound": None,
            "confidence_statement": (
                "the paper states no confidence level for this direct "
                "electromagnetic scenario bound"
            ),
            "reference_scenarios": (
                "no electromagnetic bound when the reference source scenarios "
                "contain no protons beyond 1e19 eV"
            ),
            "scenario_dependence": True,
            "public_value_used_for_diagnostic_only": True,
        },
        "exposed public input drift",
    )

    delta_lower = Fraction(-1, 10**58)
    a_squared_upper = -20 * delta_lower
    check(a_squared_upper == Fraction(20, 10**58), "independent bound algebra failed")
    check(
        receipt.get("conditional_bound")
        == {
            "inequality_derivation": (
                "delta_gamma2 = -a^2/20 > -1e-58 eV^-2 implies a^2 < 20e-58 eV^-2"
            ),
            "a_squared_strict_upper_bound_eV_minus2": "20e-58",
            "a_strict_upper_bound_exact_eV_inverse": "sqrt(20)*1e-29",
            **independent_decimals(),
            "unit_conversion": {
                "hbar_c_eV_m": "1.973269804e-7",
                "planck_length_m": "1.616255e-35",
                "constant_source": "CODATA 2022 nominal decimal values",
                "constant_uncertainties_propagated": False,
                "rounding_note": (
                    "the published coefficient limit has one-significant-digit "
                    "precision; additional conversion digits are reproducibility aids"
                ),
            },
            "conditional_on_all_open_attachments": True,
        },
        "conditional bound drift",
    )
    check(
        receipt.get("open_physical_attachments")
        == {
            "spatial_symbol_to_frequency_squared_proved": False,
            "photon_sector_attachment_proved": False,
            "electron_sector_attachment_proved": False,
            "positron_sector_attachment_proved": False,
            "shared_electromagnetic_interaction_kinematics_proved": False,
            "source_composition_and_propagation_nuisance_model_validated": False,
            "energy_frame_and_boost_convention_fixed": False,
            "carrier_orientation_map_proved": False,
            "finite_physical_carrier_scale_identified": False,
            "statement": (
                "the Auger threshold calculation couples photon, electron, positron, "
                "source-composition, frame, and interaction kinematics; the finite "
                "FZ-12 source theorem supplies none of those physical identifications. "
                "The source and propagation model may be validated or profiled under "
                "a frozen comparison contract rather than derived from OPH"
            ),
        },
        "open physical attachments drift",
    )
    check(
        receipt.get("scope_boundary")
        == {
            "public_measurement_read": True,
            "comparison_inputs": [
                "published scenario bound delta_gamma2 > -1e-58 eV^-2"
            ],
            "leading_EFT_mapping_only": True,
            "comparison_permitted": False,
            "comparison_budget_consumed": False,
            "score_emitted": False,
            "evidence_claimed": False,
            "verdict_emitted": False,
            "OPH_exclusion": False,
            "frozen_FZ12_modified": False,
            "diagnostic_type": "conditional upper bound on one open scale",
            "statement": (
                "this exposed scenario-dependent translation supplies no OPH test "
                "result until the physical sector and threshold-kinematics bridges "
                "are proved; it cannot confirm, support, falsify, or exclude OPH"
            ),
        },
        "scope boundary drift",
    )
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--custody", type=Path, default=DEFAULT_CUSTODY)
    parser.add_argument("--observation-map", type=Path, default=DEFAULT_OBSERVATION)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args(argv)
    try:
        parents = [
            load_fixed_parent(
                args.source,
                byte_count=SOURCE_BYTES,
                digest=SOURCE_SHA256,
                schema="oph.seam_current_edge_prediction_candidate.v1",
                status=(
                    "EXACT_SOURCE_NATIVE_EDGE_RAY__PROSPECTIVE_PHYSICAL_BRANCH_UNARMED__"
                    "PHYSICAL_PRODUCER_OPEN"
                ),
                self_field="receipt_sha256",
            ),
            load_fixed_parent(
                args.custody,
                byte_count=CUSTODY_BYTES,
                digest=CUSTODY_SHA256,
                schema="oph.fz12.custody_projection.v1",
                status=("FROZEN_FZ12_SOURCE_AND_CUSTODY_PINNED__NO_COMPARISON_DATA"),
                self_field="projection_sha256",
            ),
            load_fixed_parent(
                args.observation_map,
                byte_count=OBSERVATION_BYTES,
                digest=OBSERVATION_SHA256,
                schema="oph.fz12.formal_observation_map.v1",
                status=(
                    "EXACT_FORMAL_SERIES_MAP_FROM_FROZEN_FZ12_RAY__"
                    "PHYSICAL_SECTOR_CLOCK_FRAME_REMAINDER_AND_COMPARISON_OPEN"
                ),
                self_field="receipt_sha256",
            ),
        ]
        verify_receipt(args.receipt, parents)
    except (OSError, UnicodeError, VerificationError) as error:
        print(f"FZ12_AUGER_DIAGNOSTIC_VERIFY_FAIL: {error}", file=sys.stderr)
        return 1
    print("FZ12_AUGER_DIAGNOSTIC_VERIFY_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
