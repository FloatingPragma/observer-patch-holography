#!/usr/bin/env python3
"""Build fail-closed BD receipts for OPH issues #368 and #369.

The BD literature fixes a visible massless-cohomology branch.  It does not
provide the nonzero-mode spectrum, stabilized moduli, or the route-specific
data needed for a threshold calculation and low-energy decoupling map.  This
builder therefore reproduces target-side proxies, emits the issue-368 open
threshold certificate, and emits the issue-369 transverse-rank obstruction.
It never fills missing BD data with benchmark MSSM defaults, and it cannot
promote even a fully populated packet; that requires a separate evaluator.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Iterable

from jsonschema import Draft202012Validator


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_PATH = Path(__file__).resolve()
DEFAULT_PACKET = REPO_ROOT / "code/particles/data/oph_bd_threshold_spectrum_inputs.json"
DEFAULT_OUT_DIR = REPO_ROOT / "code/particles/runs/uv/oph_bd_threshold_spectrum"
SCHEMA_PATH = REPO_ROOT / "code/particles/uv/oph_bd_threshold_spectrum_inputs.schema.json"
BASE_REQUIRED_UV_INPUTS = {
    "bundle_moduli_point_and_mass_matrix",
    "bulk_five_brane_sector_and_moduli_or_absence",
    "complex_structure_moduli_point_and_mass_matrix",
    "dilaton_point_and_mass",
    "heavy_nonzero_mode_spectrum",
    "hidden_sector_completion_and_spectrum",
    "kahler_moduli_point_metric_and_mass_matrix",
    "low_energy_spectrum_solver_receipt",
    "normalized_yukawa_boundary_data",
    "physical_moduli_jacobian",
    "string_compactification_and_matching_scales",
    "supersymmetric_sector_decoupling_map",
    "threshold_matching_receipt",
}
ROUTE_CONDITIONAL_UV_INPUTS = {
    "conventional_soft_terms_boundary_conditions",
    "conventional_susy_breaking_and_mediation_data",
    "bd_equivalent_non_susy_uv_deformation_receipt",
}
ALL_UV_INPUTS = BASE_REQUIRED_UV_INPUTS | ROUTE_CONDITIONAL_UV_INPUTS
UV_RECEIPT_ENVELOPE_ARTIFACT = "oph_bd_uv_input_receipt_envelope"


class PacketError(ValueError):
    """Raised when the frozen source packet is malformed or drifts upstream."""


def _canonical_bytes(payload: Any) -> bytes:
    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _repo_relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def _resolve_repo_path(raw_path: str) -> Path:
    candidate = Path(raw_path)
    if candidate.is_absolute() or ".." in candidate.parts:
        raise PacketError(f"source path must be repository-relative: {raw_path}")
    resolved = (REPO_ROOT / candidate).resolve()
    try:
        resolved.relative_to(REPO_ROOT.resolve())
    except ValueError as exc:
        raise PacketError(f"source path escapes repository: {raw_path}") from exc
    if not resolved.is_file():
        raise PacketError(f"source path does not exist: {raw_path}")
    if resolved.suffix.lower() == ".pdf":
        raise PacketError(f"binary PDF is not a calculation dependency: {raw_path}")
    return resolved


def _json_pointer(document: Any, pointer: str) -> Any:
    if pointer == "":
        return document
    if not pointer.startswith("/"):
        raise PacketError(f"invalid JSON pointer: {pointer}")
    current = document
    for raw_token in pointer[1:].split("/"):
        token = raw_token.replace("~1", "/").replace("~0", "~")
        if isinstance(current, list):
            current = current[int(token)]
        elif isinstance(current, dict):
            if token not in current:
                raise PacketError(f"JSON pointer does not exist: {pointer}")
            current = current[token]
        else:
            raise PacketError(f"JSON pointer crosses a scalar: {pointer}")
    return current


def _decimal(value: Any) -> Decimal:
    return value if isinstance(value, Decimal) else Decimal(str(value))


def _decimal_string(value: Decimal) -> str:
    rendered = format(value, "f")
    if "." in rendered:
        rendered = rendered.rstrip("0").rstrip(".")
    return rendered or "0"


def _fraction(raw: str) -> Decimal:
    numerator, separator, denominator = raw.partition("/")
    if not separator:
        return Decimal(numerator)
    return Decimal(numerator) / Decimal(denominator)


def _derive_one_higgs_geometry(packet: dict[str, Any]) -> dict[str, int]:
    """Derive and cross-check the literature-backed determinantal counts."""

    inventory = packet["literature_inventory"]["available"]
    data = inventory["one_higgs_determinantal_data"]
    domain = int(data["map_domain_dimension_complex"])
    codomain = int(data["map_codomain_dimension_complex"])
    rank = int(data["map_rank_on_one_higgs_locus"])
    if rank > min(domain, codomain):
        raise PacketError("one-Higgs matrix rank exceeds its dimensions")
    kernel = domain - rank
    cokernel = codomain - rank
    normal = kernel * cokernel
    codimension = int(inventory["one_higgs_locus_codimension_complex"])
    if normal != codimension:
        raise PacketError(
            "one-Higgs codimension does not equal the determinantal normal dimension"
        )
    bundle_tangent = int(inventory["bundle_moduli_count"]) - codimension
    extension_tangent = int(data["extension_tangent_dimension_complex"])
    constituent_tangent = sum(
        int(value) for value in data["constituent_bundle_deformations_complex"]
    )
    if extension_tangent + constituent_tangent != bundle_tangent:
        raise PacketError(
            "one-Higgs tangent decomposition does not match the bundle-moduli count"
        )
    external_source_ids = {source["id"] for source in packet["external_sources"]}
    if not set(data["source_ids"]) <= external_source_ids:
        raise PacketError("one-Higgs determinantal data name an unknown source")
    return {
        "bundle_tangent_complex": bundle_tangent,
        "cokernel_complex": cokernel,
        "constituent_tangent_complex": constituent_tangent,
        "domain_complex": domain,
        "codomain_complex": codomain,
        "extension_tangent_complex": extension_tangent,
        "kernel_complex": kernel,
        "normal_complex": normal,
        "rank": rank,
    }


def _moduli_receipt_presence(packet: dict[str, Any]) -> dict[str, bool]:
    """Report packet-slot population without treating it as scientific verification."""

    uv_inputs = packet["bd_uv_inputs"]
    return {
        "all_completion_receipt_slots_populated": all(
            receipt is not None for receipt in uv_inputs.values()
        ),
        "physical_jacobian_receipt_slot_populated": (
            uv_inputs["physical_moduli_jacobian"] is not None
        ),
        "selected_moduli_point_declared": (
            packet["bd_branch"]["selected_moduli_point"] is not None
        ),
    }


def _load_packet(packet_path: Path) -> dict[str, Any]:
    packet = json.loads(packet_path.read_text(encoding="utf-8"))
    schema_reference = packet.get("$schema")
    if not isinstance(schema_reference, str):
        raise PacketError("source packet must name its JSON schema")
    schema_path = (packet_path.parent / schema_reference).resolve()
    try:
        schema_path.relative_to(REPO_ROOT.resolve())
    except ValueError as exc:
        raise PacketError("source-packet schema escapes repository") from exc
    if not schema_path.is_file():
        raise PacketError("source-packet schema does not exist")
    if schema_path != SCHEMA_PATH.resolve():
        raise PacketError(
            "source packet must use the canonical issue-368 schema: "
            f"{_repo_relative(SCHEMA_PATH)}"
        )
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    validation_errors = sorted(
        Draft202012Validator(schema).iter_errors(packet),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )
    if validation_errors:
        error = validation_errors[0]
        location = "/" + "/".join(str(part) for part in error.absolute_path)
        raise PacketError(
            f"source-packet schema validation failed at {location}: {error.message}"
        )
    required = {
        "artifact",
        "continuation_issues",
        "schema_version",
        "issue",
        "bd_branch",
        "comparison_registry",
        "external_sources",
        "literature_inventory",
        "target_coordinates",
        "reference_inputs",
        "proxy_assumptions",
        "numerical_policy",
        "bd_uv_inputs",
    }
    missing = sorted(required - set(packet))
    if missing:
        raise PacketError(f"source packet is missing top-level keys: {missing}")
    if packet["artifact"] != "oph_bd_threshold_spectrum_source_packet":
        raise PacketError("unexpected source-packet artifact identifier")
    if packet["schema_version"] != 2:
        raise PacketError("unsupported source-packet schema version")
    if packet["issue"] != 368:
        raise PacketError("source packet is not for issue 368")
    if packet["continuation_issues"] != [369]:
        raise PacketError("source packet must bind the issue-369 continuation")
    inventory = packet["literature_inventory"]["available"]
    if inventory["moduli_count_field"] != "complex":
        raise PacketError("BD moduli counts must be typed over the complex numbers")
    _derive_one_higgs_geometry(packet)
    external_source_ids = [source["id"] for source in packet["external_sources"]]
    if len(external_source_ids) != len(set(external_source_ids)):
        raise PacketError("external source identifiers must be unique")
    if not isinstance(packet["bd_uv_inputs"], dict):
        raise PacketError("bd_uv_inputs must be an object")
    declared_uv_inputs = set(packet["bd_uv_inputs"])
    if declared_uv_inputs != ALL_UV_INPUTS:
        missing_uv = sorted(ALL_UV_INPUTS - declared_uv_inputs)
        extra_uv = sorted(declared_uv_inputs - ALL_UV_INPUTS)
        raise PacketError(
            "bd_uv_inputs does not match the frozen contract: "
            f"missing={missing_uv}, extra={extra_uv}"
        )

    precision = int(packet["numerical_policy"]["decimal_precision_digits"])
    if precision < 50:
        raise PacketError("decimal precision must be at least 50 digits")
    getcontext().prec = precision
    getcontext().rounding = ROUND_HALF_EVEN
    return packet


def _validate_uv_receipts(packet: dict[str, Any]) -> list[dict[str, Any]]:
    """Verify each supplied receipt's bytes and issue/branch/slot envelope."""

    ledger: list[dict[str, Any]] = []
    for name, receipt in sorted(packet["bd_uv_inputs"].items()):
        if receipt is None:
            continue
        if not isinstance(receipt, dict):
            raise PacketError(f"UV receipt {name} must be null or an object")
        required = {"path", "sha256", "status"}
        if name == "supersymmetric_sector_decoupling_map":
            required.add("route")
        missing = sorted(required - set(receipt))
        allowed = set(required)
        extra = sorted(set(receipt) - allowed)
        if missing or extra:
            raise PacketError(
                f"UV receipt {name} has invalid fields: "
                f"missing={missing}, extra={extra}"
            )
        if receipt["status"] != "hash_and_envelope_verified":
            raise PacketError(
                f"UV receipt {name} must have "
                "status=hash_and_envelope_verified before it can satisfy an "
                "input slot"
            )
        expected_hash = receipt["sha256"]
        if (
            not isinstance(expected_hash, str)
            or len(expected_hash) != 64
            or any(character not in "0123456789abcdef" for character in expected_hash)
        ):
            raise PacketError(f"UV receipt {name} has an invalid SHA-256 value")
        source_path = _resolve_repo_path(receipt["path"])
        if source_path.suffix.lower() != ".json":
            raise PacketError(f"UV receipt {name} must be a JSON envelope")
        actual_hash = _sha256_file(source_path)
        if actual_hash != expected_hash:
            raise PacketError(
                f"UV receipt hash mismatch for {name}: "
                f"packet={expected_hash}, source={actual_hash}"
            )
        try:
            envelope = json.loads(source_path.read_text(encoding="utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise PacketError(f"UV receipt {name} is not valid UTF-8 JSON") from exc
        if not isinstance(envelope, dict):
            raise PacketError(f"UV receipt {name} must contain a JSON object")
        expected_envelope = {
            "artifact": UV_RECEIPT_ENVELOPE_ARTIFACT,
            "bd_branch_identifier": packet["bd_branch"]["identifier"],
            "issue": 368,
            "schema_version": 1,
            "uv_input": name,
        }
        for field, expected_value in expected_envelope.items():
            if envelope.get(field) != expected_value:
                raise PacketError(
                    f"UV receipt {name} has wrong envelope field {field}: "
                    f"expected={expected_value!r}, actual={envelope.get(field)!r}"
                )
        row = {
            "byte_count": source_path.stat().st_size,
            "path": _repo_relative(source_path),
            "sha256": actual_hash,
            "status": receipt["status"],
            "uv_input": name,
        }
        if "route" in receipt:
            row["route"] = receipt["route"]
        ledger.append(row)
    return ledger


def _validate_local_sources(
    packet: dict[str, Any], packet_path: Path
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    ledger_by_path: dict[str, dict[str, Any]] = {}
    selectors: list[dict[str, Any]] = []
    external_source_ids = {source["id"] for source in packet["external_sources"]}

    coordinate_groups = (
        ("target", packet["target_coordinates"]),
        ("reference", packet["reference_inputs"]),
    )
    for group_name, coordinates in coordinate_groups:
        for coordinate_name, record in coordinates.items():
            source = record.get("local_source")
            if source is None:
                external_source_id = record.get("external_source_id")
                if external_source_id is not None:
                    if external_source_id not in external_source_ids:
                        raise PacketError(
                            f"unknown external source for {coordinate_name}: "
                            f"{external_source_id}"
                        )
                    selectors.append(
                        {
                            "claim_status": record["claim_status"],
                            "coordinate": coordinate_name,
                            "external_source_id": external_source_id,
                            "group": group_name,
                            "value": str(record["value"]),
                            "uncertainty": record.get("uncertainty"),
                        }
                    )
                continue
            source_path = _resolve_repo_path(source["path"])
            document = json.loads(source_path.read_text(encoding="utf-8"))
            actual = _json_pointer(document, source["json_pointer"])
            expected = record["value"]
            if _decimal(actual) != _decimal(expected):
                raise PacketError(
                    f"upstream value drift for {group_name}.{coordinate_name}: "
                    f"packet={expected}, source={actual}"
                )
            promotion_pointer = source.get("promotion_json_pointer")
            promotion_value = None
            if promotion_pointer is not None:
                promotion_value = _json_pointer(document, promotion_pointer)
                required_value = source["promotion_required_value"]
                if promotion_value is not required_value:
                    raise PacketError(
                        f"upstream promotion-status drift for {coordinate_name}: "
                        f"packet={required_value}, source={promotion_value}"
                    )

            relative = _repo_relative(source_path)
            ledger_by_path.setdefault(
                relative,
                {
                    "byte_count": source_path.stat().st_size,
                    "path": relative,
                    "sha256": _sha256_file(source_path),
                },
            )
            selectors.append(
                {
                    "claim_status": record["claim_status"],
                    "coordinate": coordinate_name,
                    "group": group_name,
                    "json_pointer": source["json_pointer"],
                    "path": relative,
                    "promotion_json_pointer": promotion_pointer,
                    "promotion_value": promotion_value,
                    "value": str(expected),
                }
            )

    packet_relative = _repo_relative(packet_path)
    ledger_by_path[packet_relative] = {
        "byte_count": packet_path.stat().st_size,
        "path": packet_relative,
        "sha256": _sha256_file(packet_path),
    }
    script_relative = _repo_relative(SCRIPT_PATH)
    ledger_by_path[script_relative] = {
        "byte_count": SCRIPT_PATH.stat().st_size,
        "path": script_relative,
        "sha256": _sha256_file(SCRIPT_PATH),
    }
    schema_reference = packet.get("$schema")
    if not isinstance(schema_reference, str):
        raise PacketError("source packet must name its JSON schema")
    schema_path = (packet_path.parent / schema_reference).resolve()
    try:
        schema_path.relative_to(REPO_ROOT.resolve())
    except ValueError as exc:
        raise PacketError("source-packet schema escapes repository") from exc
    if not schema_path.is_file():
        raise PacketError("source-packet schema does not exist")
    schema_relative = _repo_relative(schema_path)
    ledger_by_path[schema_relative] = {
        "byte_count": schema_path.stat().st_size,
        "path": schema_relative,
        "sha256": _sha256_file(schema_path),
    }
    return sorted(ledger_by_path.values(), key=lambda row: row["path"]), selectors


def _target_values(packet: dict[str, Any]) -> dict[str, Decimal]:
    return {
        name: _decimal(record["value"])
        for name, record in packet["target_coordinates"].items()
    }


def _reference_values(packet: dict[str, Any]) -> dict[str, Decimal]:
    values: dict[str, Decimal] = {}
    for name, record in packet["reference_inputs"].items():
        values[name] = _decimal(record["value"])
        if "uncertainty" in record:
            values[f"{name}_uncertainty"] = _decimal(record["uncertainty"])
    return values


def _low_energy_proxy(
    targets: dict[str, Decimal], pi: Decimal
) -> dict[str, str]:
    two = Decimal(2)
    y_t = two.sqrt() * targets["mt_gev"] / targets["v_gev"]
    lambda_h = targets["mH_gev"] ** 2 / (two * targets["v_gev"] ** 2)
    lambda_mssm_max = (pi / two) * (
        targets["alpha2_mz"] + targets["alphaY_mz"]
    )
    m_h_tree_max = (
        two * lambda_mssm_max * targets["v_gev"] ** 2
    ).sqrt()
    return {
        "Delta_lambda_large_tan_beta_coordinate": _decimal_string(
            lambda_h - lambda_mssm_max
        ),
        "lambda_H_tree_coordinate": _decimal_string(lambda_h),
        "lambda_MSSM_tree_ceiling_coordinate": _decimal_string(lambda_mssm_max),
        "m_h_tree_ceiling_GeV": _decimal_string(m_h_tree_max),
        "y_t_tree_coordinate": _decimal_string(y_t),
    }


def _stop_proxy_rows(
    packet: dict[str, Any],
    targets: dict[str, Decimal],
    pi: Decimal,
) -> list[dict[str, Any]]:
    two = Decimal(2)
    twelve = Decimal(12)
    lambda_tree_max = (pi / two) * (
        targets["alpha2_mz"] + targets["alphaY_mz"]
    )
    m_tree_max = (
        two * lambda_tree_max * targets["v_gev"] ** 2
    ).sqrt()
    prefactor = (
        Decimal(3)
        * targets["mt_gev"] ** 4
        / (two * pi**2 * targets["v_gev"] ** 2)
    )
    rows: list[dict[str, Any]] = []
    for tan_beta_raw in packet["proxy_assumptions"]["tan_beta_scan"]:
        tan_beta = _decimal(tan_beta_raw)
        cos_2beta = (Decimal(1) - tan_beta**2) / (
            Decimal(1) + tan_beta**2
        )
        tree_m2 = (m_tree_max * cos_2beta) ** 2
        delta_m2 = targets["mH_gev"] ** 2 - tree_m2
        bracket = delta_m2 / prefactor
        for mixing_raw in packet["proxy_assumptions"]["stop_mixing_scan"]:
            mixing = (
                Decimal(0)
                if mixing_raw == "0"
                else Decimal(6).sqrt()
            )
            mixing_term = mixing**2 * (
                Decimal(1) - mixing**2 / twelve
            )
            m_s = targets["mt_gev"] * (
                (bracket - mixing_term) / two
            ).exp()
            hierarchy_ratio = m_s / targets["mt_gev"]
            rows.append(
                {
                    "M_S_over_m_t": _decimal_string(hierarchy_ratio),
                    "X_t_over_M_S": _decimal_string(mixing),
                    "large_log": _decimal_string(
                        (m_s**2 / targets["mt_gev"] ** 2).ln()
                    ),
                    "proxy_validity_warning": (
                        "no_BD_soft_terms_and_no_controlled_running_scheme"
                    ),
                    "required_M_S_GeV": _decimal_string(m_s),
                    "tan_beta": _decimal_string(tan_beta),
                }
            )
    return rows


def _gauge_proxy_row(
    m_susy: Decimal,
    targets: dict[str, Decimal],
    references: dict[str, Decimal],
    beta_sm: list[Decimal],
    beta_mssm: list[Decimal],
    pi: Decimal,
) -> dict[str, Any]:
    alpha1_mz = Decimal(5) * targets["alphaY_mz"] / Decimal(3)
    inverse_mz = (
        Decimal(1) / alpha1_mz,
        Decimal(1) / targets["alpha2_mz"],
        Decimal(1) / references["alpha3_reference_mz"],
    )
    log_susy_over_mz = (
        m_susy / references["mZ_reference_gev"]
    ).ln()
    inverse_susy = [
        inverse_mz[index]
        - beta_sm[index] / (Decimal(2) * pi) * log_susy_over_mz
        for index in range(3)
    ]
    log_unification_over_susy = (
        (inverse_susy[0] - inverse_susy[1])
        * Decimal(2)
        * pi
        / (beta_mssm[0] - beta_mssm[1])
    )
    m_unification = m_susy * log_unification_over_susy.exp()
    alpha_u_inverse = (
        inverse_susy[0]
        - beta_mssm[0]
        / (Decimal(2) * pi)
        * log_unification_over_susy
    )
    inverse_alpha3_susy_pred = (
        alpha_u_inverse
        + beta_mssm[2]
        / (Decimal(2) * pi)
        * log_unification_over_susy
    )
    inverse_alpha3_mz_pred = (
        inverse_alpha3_susy_pred
        + beta_sm[2]
        / (Decimal(2) * pi)
        * log_susy_over_mz
    )
    alpha3_pred = Decimal(1) / inverse_alpha3_mz_pred
    alpha3 = references["alpha3_reference_mz"]
    alpha3_uncertainty = references[
        "alpha3_reference_mz_uncertainty"
    ]
    alpha3_low = alpha3 - alpha3_uncertainty
    alpha3_high = alpha3 + alpha3_uncertainty
    residual_central = Decimal(1) / alpha3 - inverse_alpha3_mz_pred
    residual_low = Decimal(1) / alpha3_high - inverse_alpha3_mz_pred
    residual_high = Decimal(1) / alpha3_low - inverse_alpha3_mz_pred

    return {
        "M_SUSY_GeV": _decimal_string(m_susy),
        "M_U_GeV": _decimal_string(m_unification),
        "alpha3_pred_mZ": _decimal_string(alpha3_pred),
        "alpha_U_inverse": _decimal_string(alpha_u_inverse),
        "log10_MU_GeV": _decimal_string(
            m_unification.log10()
        ),
        "required_combined_inverse_coupling_correction": {
            "central": _decimal_string(residual_central),
            "interval_from_alpha3_reference_only": [
                _decimal_string(residual_low),
                _decimal_string(residual_high),
            ],
            "matching_convention": (
                "alpha3_reference_mz^(-1)-alpha3_pred_mz^(-1), evaluated "
                "at mZ after the declared one-loop SM/MSSM step running"
            ),
        },
    }


def _gauge_proxy(
    packet: dict[str, Any],
    targets: dict[str, Decimal],
    references: dict[str, Decimal],
    pi: Decimal,
) -> dict[str, Any]:
    beta_sm = [
        _fraction(value)
        for value in packet["proxy_assumptions"]["sm_beta_coefficients"]
    ]
    beta_mssm = [
        _fraction(value)
        for value in packet["proxy_assumptions"][
            "mssm_beta_coefficients"
        ]
    ]
    rows = [
        _gauge_proxy_row(
            _decimal(m_susy),
            targets,
            references,
            beta_sm,
            beta_mssm,
            pi,
        )
        for m_susy in packet["proxy_assumptions"]["m_susy_scan_gev"]
    ]
    return {
        "alpha3_reference_mz": _decimal_string(
            references["alpha3_reference_mz"]
        ),
        "alpha3_reference_uncertainty": _decimal_string(
            references["alpha3_reference_mz_uncertainty"]
        ),
        "rows": rows,
    }


def _external_source_ledger(packet: dict[str, Any]) -> list[dict[str, Any]]:
    ledger: list[dict[str, Any]] = []
    for source in packet["external_sources"]:
        hash_key = (
            "sha256_source_archive"
            if "sha256_source_archive" in source
            else "sha256_pdf"
        )
        ledger.append(
            {
                "id": source["id"],
                "role": source["role"],
                "sha256": source[hash_key],
                "source_kind": (
                    "arxiv_source_archive"
                    if hash_key == "sha256_source_archive"
                    else "published_pdf"
                ),
                "url": source["url"],
                "verification_scope": (
                    "recorded_hash; this builder does not fetch external bytes"
                ),
            }
        )
    return ledger


def _missing_uv_inputs(packet: dict[str, Any]) -> list[str]:
    uv_inputs = packet["bd_uv_inputs"]
    missing = {
        name for name in BASE_REQUIRED_UV_INPUTS if uv_inputs.get(name) is None
    }
    decoupling = uv_inputs.get("supersymmetric_sector_decoupling_map")
    if decoupling is not None:
        route = decoupling.get("route")
        if route == "conventional_susy_breaking":
            for name in (
                "conventional_susy_breaking_and_mediation_data",
                "conventional_soft_terms_boundary_conditions",
            ):
                if uv_inputs.get(name) is None:
                    missing.add(name)
        elif route == "bd_equivalent_non_susy_uv_deformation":
            if (
                uv_inputs.get("bd_equivalent_non_susy_uv_deformation_receipt")
                is None
            ):
                missing.add("bd_equivalent_non_susy_uv_deformation_receipt")
        else:
            raise PacketError(
                "decoupling route must be conventional_susy_breaking or "
                "bd_equivalent_non_susy_uv_deformation"
            )
    return sorted(missing)


def build_receipts(
    packet_path: Path = DEFAULT_PACKET,
) -> dict[str, dict[str, Any]]:
    packet_path = packet_path.resolve()
    packet = _load_packet(packet_path)
    local_ledger, selectors = _validate_local_sources(packet, packet_path)
    uv_receipt_ledger = _validate_uv_receipts(packet)
    pi = _decimal(packet["numerical_policy"]["pi_decimal"])
    targets = _target_values(packet)
    references = _reference_values(packet)
    missing_uv_inputs = _missing_uv_inputs(packet)
    target_status_blockers = sorted(
        name
        for name, record in packet["target_coordinates"].items()
        if "not_promoted" in record["claim_status"]
        or "candidate_only" in record["claim_status"]
    )

    source_dag = {
        "artifact": "oph_bd_threshold_spectrum_source_dag",
        "external_sources": _external_source_ledger(packet),
        "flow": [
            "published_BD_sources -> massless_branch_inventory",
            "local_OPH_and_compare_only_sources -> proxy_target_coordinates",
            "proxy_target_coordinates -> proxy_targets_only",
            "BD_UV_inputs -> threshold_and_low_energy_forward_map",
            "threshold_and_low_energy_forward_map -> OPH_target_comparison",
            "physical_moduli_jacobian -> moduli_locking_rank_test",
        ],
        "forbidden_edge": (
            "proxy_targets -> BD_threshold_compatibility_or_full_witness_promotion"
        ),
        "continuation_issues": packet["continuation_issues"],
        "issue": 368,
        "local_dependencies": local_ledger + uv_receipt_ledger,
        "numeric_selectors": selectors,
        "status": "frozen_source_dag_with_open_BD_UV_inputs",
        "decoupling_routes": {
            "conventional_susy_breaking": [
                "conventional_susy_breaking_and_mediation_data",
                "conventional_soft_terms_boundary_conditions",
            ],
            "bd_equivalent_non_susy_uv_deformation": [
                "bd_equivalent_non_susy_uv_deformation_receipt"
            ],
            "scope": (
                "The route belongs to the BD-to-low-energy interface. The "
                "OPH target used here has no supersymmetric partner sector."
            ),
        },
    }

    scheme_lock = {
        "artifact": "oph_bd_threshold_spectrum_scheme_lock",
        "gauge_proxy": {
            "alpha1_normalization": packet["proxy_assumptions"][
                "alpha1_normalization"
            ],
            "loop_order": packet["proxy_assumptions"]["gauge_loop_order"],
            "matching": "single_common_SM_to_MSSM_step",
            "renormalization_scheme": "not_fixed_by_BD_sources",
        },
        "higgs_stop_proxy": {
            "formula_class": "leading_one_loop_top_stop_proxy",
            "mass_inputs": "declared_pole_or_surface_coordinates",
            "renormalization_scheme": "not_a_common_running_scheme",
            "route_scope": (
                "conventional SUSY-breaking route only; a non-SUSY route that "
                "certifies BD must be an independently consistent, "
                "OPH-equivalent deformation and derive its own masses and "
                "decoupling thresholds"
            ),
        },
        "numeric_policy": packet["numerical_policy"],
        "precision_scope": (
            "80 digits control arithmetic and deterministic serialization. "
            "Physical significance remains limited by the frozen input "
            "precision, comparison uncertainty, loop order, and open scheme."
        ),
        "promotion_allowed": False,
        "status": "OPEN_SCHEME_AND_BD_BOUNDARY_DATA_NOT_SUPPLIED",
    }

    beta_provenance = {
        "artifact": "oph_bd_threshold_spectrum_beta_provenance",
        "rows": [
            {
                "coefficients": packet["proxy_assumptions"][
                    "sm_beta_coefficients"
                ],
                "effective_theory": "Standard Model",
                "loop_order": 1,
                "source_status": "borrowed_standard_QFT_coefficients",
            },
            {
                "coefficients": packet["proxy_assumptions"][
                    "mssm_beta_coefficients"
                ],
                "effective_theory": "MSSM",
                "loop_order": 1,
                "source_status": "borrowed_standard_QFT_coefficients",
            },
        ],
        "bd_heavy_field_beta_contributions": None,
        "promotion_allowed": False,
        "status": "proxy_coefficients_frozen_BD_heavy_coefficients_missing",
    }

    proxy_targets = {
        "artifact": "oph_bd_threshold_spectrum_proxy_targets",
        "bd_branch_values_used": False,
        "gauge_unification_proxy": _gauge_proxy(
            packet, targets, references, pi
        ),
        "issue": 368,
        "low_energy_tree_coordinates": _low_energy_proxy(targets, pi),
        "proxy_scope": (
            "target-side arithmetic only; no BD heavy spectrum, soft terms, "
            "thresholds, or spectrum-generator output enters"
        ),
        "status": "PROXY_TARGETS_REPRODUCED_NOT_BD_COMPATIBILITY",
        "stop_threshold_proxy": _stop_proxy_rows(
            packet, targets, pi
        ),
    }

    threshold_certificate = {
        "acceptance_criteria": {
            "proxy_numbers_are_not_a_heavy_spectrum_certificate": True,
            "selected_string_witness_promoted_to_full_OPH_target": False,
            "threshold_compatibility_backed_by_reproducible_calculation": False,
            "threshold_gate_remains_explicitly_open": True,
        },
        "artifact": "oph_bd_threshold_spectrum_certificate",
        "bd_threshold_spectrum_certificate_receipt": False,
        "compatibility_evaluated": False,
        "first_blocked_gate": (
            "heavy_nonzero_mode_spectrum"
            if "heavy_nonzero_mode_spectrum" in missing_uv_inputs
            else missing_uv_inputs[0]
            if missing_uv_inputs
            else "unpromoted_target_surface"
        ),
        "issue": 368,
        "missing_source_objects": missing_uv_inputs,
        "promotion_allowed": False,
        "reason": (
            "The published BD sources provide a massless cohomology spectrum, "
            "not the stabilized UV data or a map that decouples the "
            "supersymmetric compactification into the non-supersymmetric OPH "
            "low-energy branch."
        ),
        "decoupling_route_policy": {
            "conventional_susy_breaking": (
                "requires mediation and soft boundary conditions"
            ),
            "bd_equivalent_non_susy_uv_deformation": (
                "requires an independently consistent UV deformation or "
                "continuation proven OPH-equivalent to the BD branch, "
                "including preservation of its cited visible data, worldsheet "
                "or modular consistency, anomaly cancellation, a stable "
                "vacuum, spectrum, couplings, and thresholds"
            ),
            "oph_incorporates_supersymmetry": False,
        },
        "status": "OPEN_SOURCE_DATA_INSUFFICIENT",
        "target_surface_blockers": target_status_blockers,
    }

    inventory = packet["literature_inventory"]["available"]
    one_higgs_geometry = _derive_one_higgs_geometry(packet)
    receipt_presence = _moduli_receipt_presence(packet)
    bundle_tangent_complex = one_higgs_geometry["bundle_tangent_complex"]
    published_slice_complex = (
        inventory["kahler_moduli_count"]
        + inventory["complex_structure_moduli_count"]
        + bundle_tangent_complex
    )
    published_slice_real = 2 * published_slice_complex
    target_coordinate_order = [
        "alpha2_mz",
        "alphaY_mz",
        "v_gev",
        "mH_gev",
        "mt_gev",
    ]
    if set(target_coordinate_order) != set(packet["target_coordinates"]):
        raise PacketError("issue-369 target-coordinate registry drifted")
    target_dimension_real = len(target_coordinate_order)
    promoted_coordinate_count = sum(
        "not_promoted" not in packet["target_coordinates"][name]["claim_status"]
        and "candidate_only"
        not in packet["target_coordinates"][name]["claim_status"]
        for name in target_coordinate_order
    )

    moduli_certificate = {
        "artifact": "oph_bd_moduli_locking_certificate",
        "bd_moduli_locking_certificate_receipt": False,
        "certificate_issue": 369,
        "compatibility_evaluated": False,
        "flat_direction_classification": {
            "families": [
                {
                    "classification": "OPH-visible published modulus",
                    "complex_dimension": inventory["kahler_moduli_count"],
                    "family": "complexified_Kahler",
                    "physical_meaning": (
                        "volumes and B-field axions; they enter G4, gauge "
                        "kinetics, KK/winding scales, and thresholds"
                    ),
                    "status": "published_proxy_flat_completion_survival_unresolved",
                },
                {
                    "classification": "OPH-visible published modulus",
                    "complex_dimension": inventory[
                        "complex_structure_moduli_count"
                    ],
                    "family": "complex_structure",
                    "physical_meaning": (
                        "periods, spectra, normalized Yukawas, and thresholds"
                    ),
                    "status": "published_proxy_flat_completion_survival_unresolved",
                },
                {
                    "classification": "OPH-visible published modulus",
                    "complex_dimension": bundle_tangent_complex,
                    "decomposition_complex": {
                        "constituent_bundle_deformations": one_higgs_geometry[
                            "constituent_tangent_complex"
                        ],
                        "one_Higgs_extension_tangents": one_higgs_geometry[
                            "extension_tangent_complex"
                        ],
                    },
                    "family": "bundle_tangent_to_one_Higgs_locus",
                    "physical_meaning": (
                        "bundle cohomology, Yukawas, singlet couplings, "
                        "operator rules, and heavy/vectorlike masses"
                    ),
                    "status": "published_proxy_flat_completion_survival_unresolved",
                },
                {
                    "classification": "OPH-visible branch-exit",
                    "complex_dimension": inventory[
                        "one_higgs_locus_codimension_complex"
                    ],
                    "family": "bundle_normal_to_one_Higgs_locus",
                    "physical_meaning": (
                        "changes the Higgs/lepton mass pairing and leaves the "
                        "one-massless-Higgs-pair stratum"
                    ),
                    "status": "not_a_flat_direction_inside_the_n_equals_1_slice",
                },
            ],
            "missing_uncounted_visible_families": [
                "dilaton_axion",
                "hidden_bundle_or_bulk_M5_moduli",
                "physical_SUSY_breaking_or_non_SUSY_decoupling_data",
                "heavy_spectrum_and_matching_scales",
            ],
            "completion_dimension_policy": (
                "Missing completion variables can add ambient coordinates, "
                "while stabilization and vacuum constraints can remove "
                "directions; their net effect is not certified."
            ),
            "oph_invisible_and_quotiented": [
                "diffeomorphism_and_bundle_or_B_field_gauge_orbits",
                "bundle_isomorphisms_and_extension_rescaling_already_removed_by_the_cohomological_moduli_count",
                "exact_field_basis_and_local_representative_changes",
                "exact_scheme_scale_or_duality_presentation_changes_with_identical_complete_OPH_readout",
            ],
            "policy": (
                "A direction is not declared OPH-invisible merely because the "
                "five-coordinate comparison map does not read it."
            ),
        },
        "full_transverse_rank_verified": False,
        "isolation_verified": False,
        "missing_source_objects": [
            name
            for name in missing_uv_inputs
            if "moduli" in name
            or name
            in {
                "dilaton_point_and_mass",
                "physical_moduli_jacobian",
                "string_compactification_and_matching_scales",
            }
        ],
        "one_higgs_stratum_certificate": {
            "cokernel_dimension_complex": one_higgs_geometry[
                "cokernel_complex"
            ],
            "kernel_dimension_complex": one_higgs_geometry["kernel_complex"],
            "map": (
                "M_minus:C^"
                f"{one_higgs_geometry['domain_complex']}"
                "->C^"
                f"{one_higgs_geometry['codomain_complex']}"
            ),
            "normal_dimension_complex": one_higgs_geometry["normal_complex"],
            "normal_space": "Hom(ker(M_minus),coker(M_minus))",
            "rank_on_n_equals_1": one_higgs_geometry["rank"],
            "source_ids": inventory["one_higgs_determinantal_data"][
                "source_ids"
            ],
            "scope": (
                "At a smooth point of the cited codimension-two pullback "
                "locus, the induced normal map has exact complex rank two. "
                "This certifies only the Higgs-multiplicity stratum, not "
                "full OPH moduli locking."
            ),
            "transverse_rank_complex": one_higgs_geometry["normal_complex"],
        },
        "physical_source_slice": {
            "all_completion_receipt_slots_populated": receipt_presence[
                "all_completion_receipt_slots_populated"
            ],
            "completed_constraint_locus_verified": False,
            "completed_dimension_certified": False,
            "completed_dimension_real": None,
            "definition": (
                "The zero locus of the hidden/anomaly, safety-realization, "
                "spectrum, threshold, vacuum, stabilization, and decoupling "
                "equations inside the published one-Higgs moduli space times "
                "completion data, quotiented only by independently verified "
                "presentation redundancies. Threshold data are derived on "
                "that locus, not benchmark tuning knobs."
            ),
            "published_precompletion_slice": {
                "bundle_tangent_dimension_complex": bundle_tangent_complex,
                "dimension_complex": published_slice_complex,
                "dimension_real": published_slice_real,
                "scope": (
                    "smooth published visible-sector one-Higgs solution "
                    "deformations after gauge and presentation quotients, "
                    "before unprovided completion constraints"
                ),
            },
        },
        "promotion_allowed": False,
        "published_parameter_counts": packet["literature_inventory"][
            "available"
        ],
        "rank_obstruction_certificate_receipt": True,
        "rank_test": {
            "physical_jacobian_receipt_slot_populated": receipt_presence[
                "physical_jacobian_receipt_slot_populated"
            ],
            "physical_jacobian_verified": False,
            "completed_slice_dimension_real": None,
            "emitted_proxy_jacobian_rank": 0,
            "emitted_proxy_jacobian_shape": [
                target_dimension_real,
                published_slice_real,
            ],
            "full_column_rank_required_dimension_real": None,
            "published_precompletion_proxy_flat_dimension_real": (
                published_slice_real
            ),
            "published_slice_five_coordinate_infinitesimal_nullity_lower_bound": (
                published_slice_real - target_dimension_real
            ),
            "published_slice_five_coordinate_rank_upper_bound": (
                target_dimension_real
            ),
            "rank_deficiency_caveat": (
                "The 137-dimensional bound applies only to maps on the full "
                "published pre-completion slice. Unknown stabilization and "
                "vacuum constraints can reduce the completed slice. For a "
                "nonconstant map, kernel vectors are only infinitesimal "
                "target-null directions unless constant rank or explicit "
                "target-preserving curves are certified. The emitted proxy "
                "is constant, so all 142 published directions are actual "
                "proxy-flat directions."
            ),
            "selected_moduli_point_declared": receipt_presence[
                "selected_moduli_point_declared"
            ],
            "selected_moduli_point_verified": False,
        },
        "source_packet_issue": 368,
        "status": "FAILED_NO_COMPLETED_RANK_ISOLATION_CERTIFICATE",
        "target_map": {
            "common_threshold_scheme_fixed": packet["comparison_registry"][
                "common_threshold_scheme_fixed"
            ],
            "complete_oph_target": packet["comparison_registry"][
                "complete_oph_target"
            ],
            "comparison_values": {
                name: packet["target_coordinates"][name]["value"]
                for name in target_coordinate_order
            },
            "continuous_coordinate_order": target_coordinate_order,
            "coordinate_roles": {
                name: packet["target_coordinates"][name]["claim_status"]
                for name in target_coordinate_order
            },
            "definition": (
                "F_BD,n=1 maps a constrained physical BD point and its "
                "derived threshold receipt to the ordered five-coordinate "
                "comparison registry. The first three coordinates have "
                "declared OPH-surface status; mH and mt are candidate-only. "
                "Discrete group, generation, Higgs-count, exotic, and safety "
                "gates define the stratum and do not add Jacobian rows."
            ),
            "emitted_proxy_depends_on_bd_branch_values": False,
            "emitted_proxy_is_physical_map": False,
            "physical_forward_map_evaluable": False,
            "promoted_coordinate_count": promoted_coordinate_count,
            "target_dimension_real": target_dimension_real,
            "threshold_data_policy": (
                "physical heavy thresholds and decoupling data are derived "
                "from a fixed UV completion; convention-only scheme changes "
                "are quotiented, while genuine boundary parameters are "
                "visible source coordinates"
            ),
        },
        "verdict": {
            "operator_safe_selected_candidate": (
                "RETRACTED_BY_MODULI_LOCKING_GATE"
            ),
            "passing_oph_string_witness": False,
            "recovered_oph_core": "UNAFFECTED",
            "structural_bd_zero_mode_audit": "RETAINED_AS_BENCHMARK",
        },
    }

    return {
        "source_dag.json": source_dag,
        "scheme_lock.json": scheme_lock,
        "beta_provenance.json": beta_provenance,
        "proxy_targets.json": proxy_targets,
        "threshold_spectrum_certificate.json": threshold_certificate,
        "moduli_locking_certificate.json": moduli_certificate,
    }


def build_manifest(
    receipts: dict[str, dict[str, Any]],
    packet_path: Path = DEFAULT_PACKET,
) -> dict[str, Any]:
    artifacts = {
        filename: {
            "byte_count": len(_canonical_bytes(payload)),
            "sha256": _sha256_bytes(_canonical_bytes(payload)),
        }
        for filename, payload in sorted(receipts.items())
    }
    packet = _load_packet(packet_path.resolve())
    return {
        "artifact": "oph_bd_threshold_spectrum_manifest",
        "artifacts": artifacts,
        "generator": {
            "path": _repo_relative(SCRIPT_PATH),
            "sha256": _sha256_file(SCRIPT_PATH),
        },
        "input_packet": {
            "path": _repo_relative(packet_path.resolve()),
            "sha256": _sha256_file(packet_path.resolve()),
        },
        "issue": 368,
        "issues": [368, 369],
        "numeric_backend": "python_decimal_base10",
        "precision_digits": int(
            packet["numerical_policy"]["decimal_precision_digits"]
        ),
        "reproduction_command": (
            "python3 code/particles/uv/"
            "build_oph_bd_threshold_spectrum_receipts.py --check"
        ),
        "status": "frozen_fail_closed_receipt_bundle",
    }


def _expected_files(
    packet_path: Path,
) -> dict[str, bytes]:
    receipts = build_receipts(packet_path)
    manifest = build_manifest(receipts, packet_path)
    return {
        **{
            filename: _canonical_bytes(payload)
            for filename, payload in receipts.items()
        },
        "manifest.json": _canonical_bytes(manifest),
    }


def _write_bundle(packet_path: Path, out_dir: Path) -> None:
    expected = _expected_files(packet_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    for filename, content in expected.items():
        (out_dir / filename).write_bytes(content)


def _check_bundle(packet_path: Path, out_dir: Path) -> list[str]:
    mismatches: list[str] = []
    expected = _expected_files(packet_path)
    for filename, expected_bytes in expected.items():
        path = out_dir / filename
        if not path.is_file():
            mismatches.append(f"missing:{filename}")
        elif path.read_bytes() != expected_bytes:
            mismatches.append(f"drift:{filename}")
    extras = sorted(
        path.name
        for path in out_dir.glob("*.json")
        if path.name not in expected
    )
    mismatches.extend(f"unexpected:{name}" for name in extras)
    return mismatches


def _summary(
    packet_path: Path, out_dir: Path, mismatches: Iterable[str]
) -> dict[str, Any]:
    mismatches = list(mismatches)
    threshold = build_receipts(packet_path)[
        "threshold_spectrum_certificate.json"
    ]
    try:
        display_path = _repo_relative(out_dir)
    except ValueError:
        display_path = str(out_dir)
    return {
        "bundle": display_path,
        "compatibility_evaluated": threshold["compatibility_evaluated"],
        "mismatches": mismatches,
        "promotion_allowed": threshold["promotion_allowed"],
        "status": threshold["status"],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet", type=Path, default=DEFAULT_PACKET)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare the committed bundle with an in-memory rebuild",
    )
    args = parser.parse_args(argv)
    packet_path = args.packet.resolve()
    out_dir = args.out_dir.resolve()

    try:
        if args.check:
            mismatches = _check_bundle(packet_path, out_dir)
        else:
            _write_bundle(packet_path, out_dir)
            mismatches = []
        print(
            json.dumps(
                _summary(packet_path, out_dir, mismatches),
                indent=2,
                sort_keys=True,
            )
        )
        return 1 if mismatches else 0
    except (PacketError, KeyError, ValueError, json.JSONDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
