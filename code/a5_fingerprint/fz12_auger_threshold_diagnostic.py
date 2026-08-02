#!/usr/bin/env python3
"""Build the exposed-retrospective FZ-12 Pierre Auger diagnostic.

The calculation translates the frozen FZ-12 momentum-basis ray into the
energy-basis convention used by the Pierre Auger Collaboration and then
applies one published scenario-dependent lower bound on ``delta_gamma2``.
It is a conditional diagnostic of the open carrier scale.  It is not a
prediction comparison, score, evidential update, verdict, or OPH exclusion.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from typing import Any

import a5_multipole_fixed_point_certificate as base


HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
SOURCE_PATH = RUNTIME / "seam_current_edge_prediction_receipt.json"
CUSTODY_PATH = RUNTIME / "fz12_custody_projection.json"
OBSERVATION_MAP_PATH = RUNTIME / "fz12_observation_map_receipt.json"
RECEIPT_PATH = RUNTIME / "fz12_auger_threshold_diagnostic_receipt.json"

SCHEMA = "oph.fz12.auger_threshold_diagnostic.v1"
STATUS = "EXPOSED_RETROSPECTIVE_FZ12_AUGER_THRESHOLD_BOUND__CONDITIONAL_DIAGNOSTIC_ONLY"

SOURCE_SCHEMA = "oph.seam_current_edge_prediction_candidate.v1"
SOURCE_STATUS = (
    "EXACT_SOURCE_NATIVE_EDGE_RAY__PROSPECTIVE_PHYSICAL_BRANCH_UNARMED__"
    "PHYSICAL_PRODUCER_OPEN"
)
CUSTODY_SCHEMA = "oph.fz12.custody_projection.v1"
CUSTODY_STATUS = "FROZEN_FZ12_SOURCE_AND_CUSTODY_PINNED__NO_COMPARISON_DATA"
OBSERVATION_MAP_SCHEMA = "oph.fz12.formal_observation_map.v1"
OBSERVATION_MAP_STATUS = (
    "EXACT_FORMAL_SERIES_MAP_FROM_FROZEN_FZ12_RAY__"
    "PHYSICAL_SECTOR_CLOCK_FRAME_REMAINDER_AND_COMPARISON_OPEN"
)

PARENT_CONTRACTS = (
    {
        "path": SOURCE_PATH,
        "relative_path": (
            "code/a5_fingerprint/runtime/seam_current_edge_prediction_receipt.json"
        ),
        "bytes": 9296,
        "sha256": "0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915",
        "schema": SOURCE_SCHEMA,
        "status": SOURCE_STATUS,
        "self_field": "receipt_sha256",
        "role": "frozen FZ-12 p-basis coefficient ray",
    },
    {
        "path": CUSTODY_PATH,
        "relative_path": "code/a5_fingerprint/runtime/fz12_custody_projection.json",
        "bytes": 3624,
        "sha256": "dfc6574a5bde9b576df6aeef2e03aba8602ec2723a0a4606e7942c405a57c643",
        "schema": CUSTODY_SCHEMA,
        "status": CUSTODY_STATUS,
        "self_field": "projection_sha256",
        "role": "data-free FZ-12 freeze and custody projection",
    },
    {
        "path": OBSERVATION_MAP_PATH,
        "relative_path": "code/a5_fingerprint/runtime/fz12_observation_map_receipt.json",
        "bytes": 5276,
        "sha256": "0b34218b075a3d51fe0badb00e3bb889743ae19b6dcf161f7877683e25121d17",
        "schema": OBSERVATION_MAP_SCHEMA,
        "status": OBSERVATION_MAP_STATUS,
        "self_field": "receipt_sha256",
        "role": "exact formal FZ-12 frequency and velocity map",
    },
)

# Exposed published input.  The direct electromagnetic scenario bound in the
# paper is stated without a confidence level.
DELTA_GAMMA2_LOWER_BOUND = Fraction(-1, 10**58)  # eV^-2, strict lower bound

# Declared CODATA 2022 nominal decimal values used only for unit conversion.
HBAR_C_EV_M = Decimal("1.973269804e-7")
PLANCK_LENGTH_M = Decimal("1.616255e-35")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise base.FingerprintError(message)


def sha256_hex(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def load_parent(contract: dict[str, Any]) -> tuple[bytes, dict[str, Any]]:
    path = contract["path"]
    raw = path.read_bytes()
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError(f"invalid parent JSON: {path.name}") from error
    require(isinstance(value, dict), f"parent is not an object: {path.name}")
    require(
        raw == base.canonical_json_bytes(value), f"noncanonical parent: {path.name}"
    )
    require(
        sha256_hex(raw) == contract["sha256"], f"parent raw hash drift: {path.name}"
    )
    require(len(raw) == contract["bytes"], f"parent byte-count drift: {path.name}")
    require(
        value.get("schema") == contract["schema"], f"parent schema drift: {path.name}"
    )
    require(
        value.get("status") == contract["status"], f"parent status drift: {path.name}"
    )
    self_field = contract["self_field"]
    body = {key: item for key, item in value.items() if key != self_field}
    require(
        value.get(self_field) == base.tagged_sha256(base.canonical_json_bytes(body)),
        f"parent self-digest drift: {path.name}",
    )
    return raw, value


def decimal_outputs() -> dict[str, str]:
    with localcontext() as context:
        context.prec = 80
        sqrt_twenty = Decimal(20).sqrt()
        a_eV_inverse = sqrt_twenty * Decimal("1e-29")
        a_m = a_eV_inverse * HBAR_C_EV_M
        a_planck = a_m / PLANCK_LENGTH_M
    return {
        "a_upper_approx_eV_inverse": format(a_eV_inverse, ".9E").replace("E", "e"),
        "a_upper_approx_m": format(a_m, ".9E").replace("E", "e"),
        "a_upper_approx_planck_lengths": format(a_planck, ".10f"),
        "headline_m": format(a_m, ".4E").replace("E", "e"),
        "headline_planck_lengths": format(a_planck, ".4f"),
    }


def exact_basis_translation() -> dict[str, str]:
    c = Fraction(-1, 20)
    d_iso = Fraction(1, 840)
    d_i6 = Fraction(-1, 12600)
    delta2 = c
    delta4_iso = d_iso - 2 * c * c
    delta4_i6 = d_i6
    require(delta4_iso == Fraction(-2, 525), "Auger isotropic E-basis map drift")
    require(delta4_i6 == Fraction(-1, 12600), "Auger I6 E-basis map drift")
    return {
        "c_over_a2": str(c),
        "d_iso_over_a4": str(d_iso),
        "d_I6_over_a4": str(d_i6),
        "delta_gamma2_over_a2": str(delta2),
        "delta_gamma4_iso_over_a4": str(delta4_iso),
        "delta_gamma4_I6_over_a4": str(delta4_i6),
    }


def build_receipt() -> dict[str, Any]:
    loaded = [load_parent(contract) for contract in PARENT_CONTRACTS]
    source = loaded[0][1]
    custody = loaded[1][1]
    observation = loaded[2][1]

    require(
        source.get("conditional_physical_candidate", {}).get("coefficients")
        == {
            "C4_over_a2": "-1/20",
            "B0_over_a4": "1/840",
            "B6_over_a4": "-1/12600",
        },
        "frozen FZ-12 coefficient ray drift",
    )
    require(
        custody.get("source_receipt", {}).get("sha256")
        == PARENT_CONTRACTS[0]["sha256"],
        "custody no longer pins the FZ-12 source",
    )
    observation_ancestry = observation.get("ancestry", {})
    require(
        observation_ancestry.get("source_receipt", {}).get("sha256")
        == f"sha256:{PARENT_CONTRACTS[0]['sha256']}",
        "observation map source ancestry drift",
    )
    require(
        observation_ancestry.get("custody_projection", {}).get("sha256")
        == f"sha256:{PARENT_CONTRACTS[1]['sha256']}",
        "observation map custody ancestry drift",
    )
    require(
        observation.get("certificate_scope", {}).get("physical_prediction") is False,
        "observation map was promoted to a physical prediction",
    )
    require(
        observation.get("exposure_boundary", {}).get("comparison_permitted") is False,
        "observation map unexpectedly permits comparison",
    )

    translation = exact_basis_translation()
    a_squared_upper = -20 * DELTA_GAMMA2_LOWER_BOUND
    require(a_squared_upper == Fraction(20, 10**58), "scale-bound algebra drift")
    decimal = decimal_outputs()

    parent_pins = []
    for contract, (raw, value) in zip(PARENT_CONTRACTS, loaded, strict=True):
        parent_pins.append(
            {
                "path": contract["relative_path"],
                "role": contract["role"],
                "bytes": len(raw),
                "sha256": f"sha256:{sha256_hex(raw)}",
                "schema": value["schema"],
                "status": value["status"],
                "self_digest": value[contract["self_field"]],
            }
        )

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 667,
        "diagnostic_id": "FZ12-AUGER-THRESHOLD-EXPOSED-RETROSPECTIVE",
        "parent_pins": parent_pins,
        "theory_input": {
            "basis": "momentum p basis",
            "dispersion": "E^2 = p^2 + c p^4 + (d_iso + d_I6 I6) p^6 + O(p^8)",
            "coefficients": {
                "c": "-a^2/20",
                "d_iso": "a^4/840",
                "d_I6": "-a^4/12600",
            },
            "physical_attachment_closed": False,
        },
        "basis_translation": {
            "auger_energy_basis": (
                "E^2 = p^2 + delta_gamma2 E^4 + "
                "(delta_gamma4_iso + delta_gamma4_I6 I6) E^6 + O(E^8)"
            ),
            "series_identity": (
                "E^4 = p^4 + 2 c p^6 + O(p^8), so delta_gamma2 = c and "
                "delta_gamma4 = d - 2 c^2; the isotropic c leaves the I6 "
                "coefficient unchanged"
            ),
            "exact_coefficients": translation,
            "paper_consumes_delta_gamma4": False,
            "mapping_scope": "leading EFT series through sixth order only",
            "higher_order_remainder_control": False,
        },
        "exposed_public_input": {
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
        "conditional_bound": {
            "inequality_derivation": (
                "delta_gamma2 = -a^2/20 > -1e-58 eV^-2 implies a^2 < 20e-58 eV^-2"
            ),
            "a_squared_strict_upper_bound_eV_minus2": "20e-58",
            "a_strict_upper_bound_exact_eV_inverse": "sqrt(20)*1e-29",
            **decimal,
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
        "open_physical_attachments": {
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
        "scope_boundary": {
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
    }
    receipt["receipt_sha256"] = base.tagged_sha256(base.canonical_json_bytes(receipt))
    return receipt


def verify_committed_receipt() -> dict[str, Any]:
    raw = RECEIPT_PATH.read_bytes()
    try:
        committed = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError(
            "invalid committed Auger diagnostic JSON"
        ) from error
    require(isinstance(committed, dict), "Auger diagnostic root is not an object")
    require(
        raw == base.canonical_json_bytes(committed), "noncanonical Auger diagnostic"
    )
    body = {key: value for key, value in committed.items() if key != "receipt_sha256"}
    require(
        committed.get("receipt_sha256")
        == base.tagged_sha256(base.canonical_json_bytes(body)),
        "Auger diagnostic self-digest drift",
    )
    require(
        raw == base.canonical_json_bytes(build_receipt()),
        "Auger diagnostic ancestry, arithmetic, schema, or boundary drift",
    )
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
        RUNTIME.mkdir(exist_ok=True)
        RECEIPT_PATH.write_bytes(base.canonical_json_bytes(receipt))
        print(RECEIPT_PATH)
        return 0
    print(json.dumps(receipt["conditional_bound"], indent=1))
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
