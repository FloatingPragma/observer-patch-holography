#!/usr/bin/env python3
"""Fail-closed receipt-class and literal-witness audit for GitHub issue #518.

The registry is outside every subject artifact. A subject cannot select its own
receipt class, and a boolean such as ``theorem_grade`` never establishes a
literal property. The audit recomputes the promoted Petz and hierarchy rows,
checks the Thomson and Planck comparison arithmetic, quarantines identity and
schema rows, and keeps the four A5 layers separately typed.
"""

from __future__ import annotations

import argparse
import copy
from decimal import Decimal, getcontext
from fractions import Fraction
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REGISTRY = ROOT / "code" / "audit" / "issue_518_receipt_registry.json"
RECEIPT_CLASSES = {
    "identity",
    "schema",
    "producer",
    "recomputation",
    "prediction",
    "comparison",
}
PROMOTION_ELIGIBLE = {"producer", "recomputation", "prediction"}
NEVER_SUFFICIENT = {"identity", "schema", "comparison"}
A5_ROLE_CLASSES = {
    "coefficient_algebra": "recomputation",
    "physical_currents": "producer",
    "global_descent": "identity",
    "matter_realization": "schema",
}

getcontext().prec = 100


def _json_pointer(payload: Any, pointer: str) -> Any:
    if pointer == "":
        return payload
    if not pointer.startswith("/"):
        raise ValueError(f"invalid JSON pointer: {pointer!r}")
    value = payload
    for raw_token in pointer[1:].split("/"):
        token = raw_token.replace("~1", "/").replace("~0", "~")
        if isinstance(value, list):
            value = value[int(token)]
        else:
            value = value[token]
    return value


def _load_json(
    relative_path: str,
    *,
    root: Path,
    overrides: Mapping[str, Any] | None,
) -> Any:
    if overrides and relative_path in overrides:
        return copy.deepcopy(overrides[relative_path])
    path = root / relative_path
    return json.loads(path.read_text(encoding="utf-8"))


def _interval_endpoints(text: str) -> tuple[Decimal, Decimal]:
    match = re.fullmatch(r"\[([^,]+), ([^\]]+)\]", text)
    if not match:
        raise ValueError(f"not a closed interval string: {text!r}")
    return Decimal(match.group(1)), Decimal(match.group(2))


def _check_hierarchy_interval(payload: Mapping[str, Any], _row: Mapping[str, Any]) -> list[str]:
    failures: list[str] = []
    phi_lo = _interval_endpoints(str(payload["phi_at_lower"]))
    phi_hi = _interval_endpoints(str(payload["phi_at_upper"]))
    derivative = _interval_endpoints(str(payload["derivative_enclosure_union"]))
    controls = payload.get("perturbation_controls", {})
    if not (phi_lo[0] > 0 and phi_lo[1] > 0):
        failures.append("hierarchy_lower_endpoint_not_strictly_positive")
    if not (phi_hi[0] < 0 and phi_hi[1] < 0):
        failures.append("hierarchy_upper_endpoint_not_strictly_negative")
    if not derivative[1] < 0:
        failures.append("hierarchy_derivative_not_strictly_negative")
    if not all(value is True for value in controls.values()) or len(controls) < 3:
        failures.append("hierarchy_mutation_controls_not_all_fail_closed")
    if payload.get("unique_root_certified") is not True:
        failures.append("hierarchy_unique_root_flag_disagrees_with_literal_checks")
    return failures


def _check_petz(payload: Mapping[str, Any], _row: Mapping[str, Any]) -> list[str]:
    failures: list[str] = []
    matrix = [
        [Fraction(value) for value in row]
        for row in payload.get("recovery_matrix", [])
    ]
    if len(matrix) != 9 or any(len(row) != 3 for row in matrix):
        return ["petz_recovery_matrix_shape_not_9_by_3"]
    nonnegative = all(value >= 0 for row in matrix for value in row)
    column_sums = [
        sum((row[column] for row in matrix), Fraction(0))
        for column in range(3)
    ]
    trace_preserving = column_sums == [Fraction(1)] * 3
    marginal = [
        [
            sum((matrix[3 * b_out + d][b_in] for d in range(3)), Fraction(0))
            for b_in in range(3)
        ]
        for b_out in range(3)
    ]
    identity = [
        [Fraction(int(row == column)) for column in range(3)]
        for row in range(3)
    ]
    reference = [
        Fraction(payload.get("reference_sigma_B", {}).get(str(index), "0"))
        for index in range(3)
    ]
    full_support = sum(reference, Fraction(0)) == 1 and min(reference) > 0
    induced_l1 = max(
        sum((abs(row[column]) for row in matrix), Fraction(0))
        for column in range(3)
    )
    literal = {
        "classical_choi_diagonal_nonnegative": nonnegative,
        "trace_preserving": trace_preserving,
        "basis_marginal_recovery_identity": marginal == identity,
        "cptp": nonnegative and trace_preserving,
        "trace_norm_contractive": nonnegative and trace_preserving and induced_l1 <= 1,
        "reference_state_normalized": sum(reference, Fraction(0)) == 1,
        "full_support": full_support,
    }
    for key, computed in literal.items():
        if payload.get(key) is not computed:
            failures.append(f"petz_literal_mismatch:{key}")
    if not all(literal.values()) or payload.get("pass") is not True:
        failures.append("petz_literal_witness_does_not_pass")
    return failures


def _check_thomson(
    payload: Mapping[str, Any],
    row: Mapping[str, Any],
    *,
    root: Path,
    overrides: Mapping[str, Any] | None,
) -> list[str]:
    failures: list[str] = []
    rounding_tolerance = Decimal("1e-80")
    compare = payload["codata_mapped_endpoint_packet"]
    exact = compare["exact_one_loop_package"]
    anchor = Decimal(exact["source_anchor_alpha_inv_mz"])
    lepton = Decimal(exact["lepton_delta_alpha_inv"])
    quark = Decimal(exact["quark_delta_alpha_inv_screened"])
    transport = Decimal(exact["implemented_transport_delta_alpha_inv"])
    endpoint = Decimal(exact["implemented_endpoint_alpha_inv"])
    target = Decimal(compare["compare_alpha_inv"])
    required = Decimal(exact["required_transport_delta_alpha_inv"])
    missing = Decimal(exact["missing_source_transport_delta_alpha_inv"])
    if abs((lepton + quark) - transport) > rounding_tolerance:
        failures.append("thomson_transport_components_do_not_sum")
    if abs((anchor + transport) - endpoint) > rounding_tolerance:
        failures.append("thomson_endpoint_not_recomputed_from_components")
    if abs((target - anchor) - required) > rounding_tolerance:
        failures.append("thomson_required_transport_not_recomputed")
    if abs((target - endpoint) - missing) > rounding_tolerance:
        failures.append("thomson_missing_transport_not_recomputed")
    if payload.get("promotion_allowed") is not False:
        failures.append("thomson_comparison_improperly_promoted")

    support_path = row.get("supporting_artifacts", [None])[0]
    support = _load_json(support_path, root=root, overrides=overrides)
    if support.get("promotion_allowed") is not False:
        failures.append("thomson_interval_scaffold_improperly_promoted")
    if support.get("R_Q_certificate", {}).get("status") != "missing_source_artifact":
        failures.append("thomson_missing_source_map_not_recorded")
    if support.get("conclusion", {}).get("unique_fixed_point_in_I_P") is not False:
        failures.append("thomson_unique_fixed_point_claim_not_blocked")
    rounding = support.get("interval_backend", {}).get("rounding", "")
    if "not_theorem_grade" not in rounding:
        failures.append("thomson_interval_backend_boundary_missing")
    # Deliberately do not read issue_223_acceptance.theorem_grade_object_defined.
    return failures


def _check_quarantined_identity(
    _payload: Mapping[str, Any], row: Mapping[str, Any]
) -> list[str]:
    if row.get("promoted") is not False:
        return ["identity_receipt_improperly_promoted"]
    return []


def _canonical_hash(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    import hashlib

    return hashlib.sha256(encoded).hexdigest()


def _check_nogo_source_packet(source_packet: Mapping[str, Any]) -> list[str]:
    expected_keys = {
        "branch",
        "P",
        "alpha_U",
        "beta_EW",
        "m_rep",
        "capacity_structural_axioms",
        "naturality_type_signatures",
        "excluded_target_conditions",
    }
    failures: list[str] = []
    if set(source_packet) != expected_keys:
        failures.append("nogo_source_packet_schema_mismatch")
    if source_packet.get("branch") != "source_audit_pixel_branch":
        failures.append("nogo_source_packet_not_on_source_audit_branch")
    exclusions = set(source_packet.get("excluded_target_conditions", []))
    required_exclusions = {
        "B_EW(P,N)=0",
        "Pi_EW(P,N)=4P",
        "N_CRC^EW",
        "epsilon_n",
        "epsilon_h",
        "epsilon_H",
        "measured weak scale",
        "measured Lambda or Planck capacity",
    }
    if exclusions != required_exclusions:
        failures.append("nogo_target_exclusion_set_mismatch")
    return failures


def _check_capacity_nonidentifiability(
    payload: Mapping[str, Any],
    _row: Mapping[str, Any],
) -> list[str]:
    failures: list[str] = []
    source_packet = payload.get("source_packet", {})
    failures.extend(_check_nogo_source_packet(source_packet))
    source_hash = _canonical_hash(source_packet)
    if payload.get("source_packet_sha256") != source_hash:
        failures.append("capacity_nogo_source_packet_hash_mismatch")
    if payload.get("target_artifacts_consumed") != []:
        failures.append("capacity_nogo_consumes_target_artifact")
    block = payload.get("capacity_nonidentifiability", {})
    models = block.get("countermodels", [])
    if len(models) < 2:
        return failures + ["capacity_nogo_requires_two_countermodels"]

    fixed_points: list[Fraction] = []
    for model in models:
        try:
            lo, hi = [Fraction(value) for value in model["domain_log_capacity"]]
            image_lo, image_hi = [
                Fraction(value) for value in model["image_log_capacity"]
            ]
            fixed = Fraction(model["fixed_log_capacity"])
            weight = Fraction(model["averaging_weight_lambda"])
            derivative = 1 - weight
            expected_image = (
                derivative * lo + weight * fixed,
                derivative * hi + weight * fixed,
            )
            computed = {
                "derivative": str(derivative),
                "monotone": derivative >= 0,
                "strict_contraction": abs(derivative) < 1,
                "self_map": lo <= expected_image[0] <= expected_image[1] <= hi,
                "unique_fixed_point": 0 < weight <= 1 and lo < fixed < hi,
            }
            if model.get("antecedent_fingerprint") != source_hash:
                failures.append(
                    f"capacity_nogo_antecedent_mismatch:{model.get('completion_id')}"
                )
            if (image_lo, image_hi) != expected_image:
                failures.append(
                    f"capacity_nogo_image_recomputation_failed:{model.get('completion_id')}"
                )
            for key, value in computed.items():
                displayed = model.get(key)
                if key == "derivative":
                    displayed = str(Fraction(str(displayed)))
                if displayed != value:
                    failures.append(
                        f"capacity_nogo_literal_mismatch:{model.get('completion_id')}:{key}"
                    )
            if Fraction(model.get("fixed_point_residual", "1")) != 0:
                failures.append(
                    f"capacity_nogo_fixed_point_residual_nonzero:{model.get('completion_id')}"
                )
            fixed_points.append(fixed)
        except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
            failures.append(
                "capacity_nogo_countermodel_exception:"
                f"{model.get('completion_id')}:{type(exc).__name__}:{exc}"
            )

    if len(set(fixed_points)) < 2:
        failures.append("capacity_nogo_outputs_not_distinct")
    if block.get("same_antecedents") is not True:
        failures.append("capacity_nogo_same_antecedents_not_recorded")
    if block.get("distinct_fixed_log_capacities") is not True:
        failures.append("capacity_nogo_distinct_outputs_not_recorded")
    if block.get("result") != (
        "NO_UNIQUE_CAPACITY_PRODUCER_FROM_DECLARED_ANTECEDENTS"
    ):
        failures.append("capacity_nogo_result_mismatch")
    return failures


def _finite_square_defects(completion: Mapping[str, Any]) -> dict[str, int]:
    sets = completion["sets"]
    maps = completion["maps"]
    q_source = list(sets["Q_s"])
    q_higgs = set(sets["Q_H"])
    o_source = set(sets["O_s"])
    o_higgs = set(sets["O_H"])

    expected_domains = {
        "rho_sH": set(q_source),
        "n_s": set(q_source),
        "n_H": q_higgs,
        "h_s": set(q_source),
        "chi_sH": o_source,
        "h_H": q_higgs,
    }
    expected_codomains = {
        "rho_sH": q_higgs,
        "n_s": set(q_source),
        "n_H": q_higgs,
        "h_s": o_source,
        "chi_sH": o_higgs,
        "h_H": o_higgs,
    }
    for name, domain in expected_domains.items():
        if set(maps[name]) != domain:
            raise ValueError(f"{name} is not total on its declared domain")
        if not set(maps[name].values()).issubset(expected_codomains[name]):
            raise ValueError(f"{name} leaves its declared codomain")

    rho = maps["rho_sH"]
    n_s = maps["n_s"]
    n_h = maps["n_H"]
    h_s = maps["h_s"]
    chi = maps["chi_sH"]
    h_h = maps["h_H"]
    epsilon_n = max(
        (
            int(rho[n_s[x]] != n_h[rho[x]])
            for x in q_source
        ),
        default=0,
    )
    epsilon_h = max(
        (
            int(chi[h_s[x]] != h_h[rho[x]])
            for x in q_source
        ),
        default=0,
    )
    return {
        "epsilon_n": epsilon_n,
        "epsilon_h": epsilon_h,
        "epsilon_H": max(epsilon_n, epsilon_h),
    }


def _check_naturality_nonidentifiability(
    payload: Mapping[str, Any],
    _row: Mapping[str, Any],
) -> list[str]:
    failures: list[str] = []
    source_packet = payload.get("source_packet", {})
    failures.extend(_check_nogo_source_packet(source_packet))
    source_hash = _canonical_hash(source_packet)
    if payload.get("source_packet_sha256") != source_hash:
        failures.append("naturality_nogo_source_packet_hash_mismatch")
    if payload.get("target_artifacts_consumed") != []:
        failures.append("naturality_nogo_consumes_target_artifact")
    block = payload.get("naturality_nonidentifiability", {})
    models = block.get("countermodels", [])
    if len(models) < 3:
        return failures + ["naturality_nogo_requires_three_countermodels"]

    defects: list[dict[str, int]] = []
    for model in models:
        try:
            computed = _finite_square_defects(model)
            if model.get("antecedent_fingerprint") != source_hash:
                failures.append(
                    f"naturality_nogo_antecedent_mismatch:{model.get('completion_id')}"
                )
            if model.get("evaluated_defects") != computed:
                failures.append(
                    f"naturality_nogo_defect_recomputation_failed:{model.get('completion_id')}"
                )
            defects.append(computed)
        except (KeyError, TypeError, ValueError) as exc:
            failures.append(
                "naturality_nogo_countermodel_exception:"
                f"{model.get('completion_id')}:{type(exc).__name__}:{exc}"
            )

    epsilon_h_values = [item["epsilon_H"] for item in defects]
    if 0 not in epsilon_h_values or not any(value > 0 for value in epsilon_h_values):
        failures.append("naturality_nogo_does_not_exhibit_distinct_defects")
    if not any(item["epsilon_n"] > 0 for item in defects):
        failures.append("naturality_nogo_normal_square_countermodel_missing")
    if not any(item["epsilon_h"] > 0 for item in defects):
        failures.append("naturality_nogo_obstruction_square_countermodel_missing")
    if block.get("evaluated_epsilon_H_values") != epsilon_h_values:
        failures.append("naturality_nogo_epsilon_summary_mismatch")
    if block.get("same_antecedents") is not True:
        failures.append("naturality_nogo_same_antecedents_not_recorded")
    if block.get("result") != (
        "NO_UNIQUE_NATURALITY_DEFECT_FROM_DECLARED_ANTECEDENTS"
    ):
        failures.append("naturality_nogo_result_mismatch")
    return failures


def _check_dark_negative(
    payload: Mapping[str, Any],
    _row: Mapping[str, Any],
    *,
    root: Path,
) -> list[str]:
    failures: list[str] = []
    if payload.get("negative_control") is not True:
        failures.append("dark_fixture_not_marked_negative_control")
    if payload.get("active_theorem_bundle") is not False:
        failures.append("dark_fixture_marked_active")
    marker = "dark_sector_calibration_negative"
    for manifest in (root / "code").glob("**/manifest*.json"):
        if marker in manifest.read_text(encoding="utf-8"):
            failures.append(f"dark_fixture_present_in_active_manifest:{manifest.relative_to(root)}")
    return failures


def _check_planck(payload: Mapping[str, Any], row: Mapping[str, Any]) -> list[str]:
    failures: list[str] = []
    input_pointer = str(row.get("planck_input_pointer", ""))
    output_pointer = str(row.get("planck_output_pointer", ""))
    if input_pointer == output_pointer:
        failures.append("planck_input_output_pointers_overlap")
        return failures
    inputs = _json_pointer(payload, input_pointer)
    outputs = _json_pointer(payload, output_pointer)
    n_crc = float(inputs["N_CRC"])
    c = 299_792_458.0
    hbar = 1.054_571_817e-34
    grav = 6.674_30e-11
    mpc = 3.085_677_581_491_367e22
    for name, combo in inputs["combos"].items():
        omega_l = float(combo["OmegaL"][0])
        h0 = float(combo["H0"][0])
        lam = 3.0 * omega_l * (h0 * 1e3 / mpc) ** 2 / c**2
        lp2 = hbar * grav / c**3
        computed_n = 3.0 * math.pi / (lam * lp2)
        output = outputs[name]
        if not math.isclose(float(output["N_Lambda"]), computed_n, rel_tol=2e-15):
            failures.append(f"planck_N_Lambda_recomputation_failed:{name}")
        if not math.isclose(float(output["gap_ratio"]), n_crc / computed_n, rel_tol=2e-15):
            failures.append(f"planck_gap_ratio_recomputation_failed:{name}")
        if not math.isclose(float(output["gap_ln"]), math.log(n_crc / computed_n), rel_tol=2e-15):
            failures.append(f"planck_gap_log_recomputation_failed:{name}")
    return failures


def _check_a5_coefficient(payload: Mapping[str, Any], _row: Mapping[str, Any]) -> list[str]:
    types = {
        entry.get("lie_type")
        for entry in payload.get("compact_lie_trichotomy", [])
    }
    required = {"su(3)+su(2)+u(1)", "su(2)+su(2)+u(1)^6", "u(1)^12"}
    return [] if types == required else ["a5_coefficient_trichotomy_incomplete"]


def _check_a5_current(payload: Mapping[str, Any], row: Mapping[str, Any]) -> list[str]:
    failures: list[str] = []
    if payload.get("conditional_algebraic_gate", {}).get("passed") is not True:
        failures.append("a5_conditional_current_algebra_gate_failed")
    gate = payload.get("physical_source_gate", {})
    binding = payload.get("semantic_response_binding")
    if gate.get("passed") is not True:
        failures.append("a5_physical_current_source_gate_not_passed")
    if gate.get("passed") is True and not isinstance(binding, Mapping):
        failures.append("a5_physical_current_gate_passed_without_semantic_binding")
    if isinstance(binding, Mapping) and binding.get("sector_structure_recomputed") is not True:
        failures.append("a5_physical_current_binding_not_recomputed")
    if row.get("promoted") is not False:
        failures.append("a5_physical_current_improperly_promoted")
    return failures


def _check_a5_global(payload: Mapping[str, Any], row: Mapping[str, Any]) -> list[str]:
    failures: list[str] = []
    if not payload:
        failures.append("a5_cover_kernel_identity_missing")
    if row.get("promoted") is not False:
        failures.append("a5_cover_kernel_improperly_promoted_to_global_descent")
    return failures


def _check_a5_matter(payload: Mapping[str, Any], row: Mapping[str, Any]) -> list[str]:
    failures: list[str] = []
    if payload.get("conditional_algebraic_gate", {}).get("mar_class_nonempty_witnessed") is not True:
        failures.append("a5_matter_schema_nonempty_witness_missing")
    gate = payload.get("physical_source_gate", {})
    if gate.get("passed") not in {True, False}:
        failures.append("a5_matter_physical_source_gate_status_missing")
    if gate.get("passed") is True:
        balance = payload.get("block_determinant_balance", {})
        if not (
            balance.get("declared_equals_derived") is True
            or balance.get("declared_matches_derived_pair_up_to_conjugation") is True
        ):
            failures.append("a5_matter_gate_passed_without_derived_balance")
    else:
        closure = payload.get("issue_closure_condition", {})
        if (
            closure.get("met_locally") is not False
            or closure.get("physical_source_realization_gate_passed") is not False
        ):
            failures.append("a5_matter_open_gate_not_reflected_in_closure_status")
    if row.get("promoted") is not False:
        failures.append("a5_matter_schema_improperly_promoted")
    return failures


def audit_registry(
    registry: Mapping[str, Any],
    *,
    root: Path = ROOT,
    overrides: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit a registry payload and return a machine-readable fail-closed report."""
    global_failures: list[str] = []
    rows = list(registry.get("rows", []))
    if set(registry.get("receipt_classes", [])) != RECEIPT_CLASSES:
        global_failures.append("receipt_class_vocabulary_mismatch")
    if registry.get("a5_role_classes") != A5_ROLE_CLASSES:
        global_failures.append("a5_role_class_map_mismatch")
    ids = [row.get("receipt_id") for row in rows]
    if len(ids) != len(set(ids)):
        global_failures.append("duplicate_receipt_id")
    try:
        hierarchy_manifest = _load_json(
            "code/particles/hierarchy/manifest.json",
            root=root,
            overrides=overrides,
        )
        authority = hierarchy_manifest.get("promotion_authority", {})
        if authority.get("registry") != "code/audit/issue_518_receipt_registry.json":
            global_failures.append("hierarchy_manifest_missing_independent_promotion_registry")
        closed_text = "\n".join(
            hierarchy_manifest.get("claim_boundary", {}).get("closed_by_bundle", [])
        )
        if "bridge-defined fixed-point identity N_CRC^EW" in closed_text:
            global_failures.append("backsolved_capacity_identity_present_in_closed_bundle")
        if "RG/Higgs naturality defect epsilon_H=0" in closed_text:
            global_failures.append("selected_naturality_identity_present_in_closed_bundle")
        if (
            "target-free scalar packet and generic Banach/readback conditions "
            "admit distinct capacity fixed points"
        ) not in closed_text:
            global_failures.append("capacity_nonidentifiability_missing_from_closed_bundle")
        if (
            "map type signatures admit both commuting and noncommuting "
            "comparison-square completions"
        ) not in closed_text:
            global_failures.append("naturality_nonidentifiability_missing_from_closed_bundle")
        open_text = "\n".join(
            hierarchy_manifest.get("claim_boundary", {}).get(
                "not_closed_by_bundle",
                [],
            )
        )
        if "source-derived public-record/readback map F" not in open_text:
            global_failures.append("physical_capacity_producer_blocker_missing")
        if "source-derived definitions of rho_sH" not in open_text:
            global_failures.append("physical_naturality_map_blocker_missing")
    except (KeyError, ValueError, TypeError, OSError, json.JSONDecodeError) as exc:
        global_failures.append(
            f"hierarchy_manifest_audit_exception:{type(exc).__name__}:{exc}"
        )
    try:
        nogo = _load_json(
            (
                "code/particles/hierarchy/certificates/"
                "antecedent_only_nonidentifiability_receipt.json"
            ),
            root=root,
            overrides=overrides,
        )
        source = nogo["source_packet"]
        pixel = _load_json(
            (
                "code/particles/hierarchy/certificates/"
                "R_P_source_audit_pixel_certificate.json"
            ),
            root=root,
            overrides=overrides,
        )
        rounds = _load_json(
            (
                "code/particles/hierarchy/certificates/"
                "R_m_rep_24_certificate.json"
            ),
            root=root,
            overrides=overrides,
        )
        if source.get("P") != pixel.get("P_cand"):
            global_failures.append("nogo_source_P_not_bound_to_source_audit_receipt")
        if source.get("alpha_U") != pixel.get("alpha_U_P_cand"):
            global_failures.append(
                "nogo_source_alpha_U_not_bound_to_source_audit_receipt"
            )
        if source.get("m_rep") != rounds.get("result", {}).get("m_rep"):
            global_failures.append("nogo_m_rep_not_bound_to_round_count_receipt")
    except (KeyError, ValueError, TypeError, OSError, json.JSONDecodeError) as exc:
        global_failures.append(
            f"nogo_source_binding_exception:{type(exc).__name__}:{exc}"
        )

    row_reports: list[dict[str, Any]] = []
    for row in rows:
        failures: list[str] = []
        receipt_id = str(row.get("receipt_id", "<missing>"))
        receipt_class = row.get("receipt_class")
        if receipt_class not in RECEIPT_CLASSES:
            failures.append("unknown_receipt_class")
        promoted = row.get("promoted") is True
        if promoted and receipt_class not in PROMOTION_ELIGIBLE:
            failures.append("receipt_class_not_promotion_eligible")
        if promoted and row.get("open_gates"):
            failures.append("promoted_row_has_open_gates")
        antecedents = set(row.get("antecedent_artifacts", []))
        targets = set(row.get("target_artifacts", []))
        if promoted and antecedents.intersection(targets):
            failures.append("promoted_row_consumes_target_artifact")
        if promoted and not row.get("mutation_controls"):
            failures.append("promoted_row_has_no_mutation_controls")

        role = row.get("role")
        if role in A5_ROLE_CLASSES and receipt_class != A5_ROLE_CLASSES[role]:
            failures.append("a5_role_receipt_class_mismatch")

        try:
            raw_payload = _load_json(str(row["artifact"]), root=root, overrides=overrides)
            payload = _json_pointer(raw_payload, str(row.get("json_pointer", "")))
            kind = row.get("literal_check")
            if kind == "hierarchy_interval":
                failures.extend(_check_hierarchy_interval(payload, row))
            elif kind == "petz":
                failures.extend(_check_petz(payload, row))
            elif kind == "thomson":
                failures.extend(
                    _check_thomson(
                        payload,
                        row,
                        root=root,
                        overrides=overrides,
                    )
                )
            elif kind == "quarantined_identity":
                failures.extend(_check_quarantined_identity(payload, row))
            elif kind == "capacity_nonidentifiability":
                failures.extend(_check_capacity_nonidentifiability(payload, row))
            elif kind == "naturality_nonidentifiability":
                failures.extend(_check_naturality_nonidentifiability(payload, row))
            elif kind == "dark_negative_control":
                failures.extend(_check_dark_negative(payload, row, root=root))
            elif kind == "planck_separation":
                failures.extend(_check_planck(payload, row))
            elif kind == "a5_coefficient":
                failures.extend(_check_a5_coefficient(payload, row))
            elif kind == "a5_current":
                failures.extend(_check_a5_current(payload, row))
            elif kind == "a5_global_descent":
                failures.extend(_check_a5_global(payload, row))
            elif kind == "a5_matter":
                failures.extend(_check_a5_matter(payload, row))
            else:
                failures.append("unknown_literal_check")
        except (KeyError, ValueError, TypeError, OSError, json.JSONDecodeError) as exc:
            failures.append(f"literal_check_exception:{type(exc).__name__}:{exc}")

        row_reports.append(
            {
                "receipt_id": receipt_id,
                "receipt_class": receipt_class,
                "promoted": promoted,
                "pass": not failures,
                "failures": failures,
            }
        )

    all_failures = global_failures + [
        f"{row['receipt_id']}:{failure}"
        for row in row_reports
        for failure in row["failures"]
    ]
    return {
        "artifact": "issue_518_receipt_promotion_audit",
        "issue": 518,
        "class_source": "independent_registry_only",
        "subject_declared_classes_ignored": True,
        "rows": row_reports,
        "promoted_receipts": [
            row["receipt_id"]
            for row in row_reports
            if row["promoted"] and row["pass"]
        ],
        "global_failures": global_failures,
        "failures": all_failures,
        "pass": not all_failures,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", default=str(DEFAULT_REGISTRY))
    args = parser.parse_args()
    registry = json.loads(Path(args.registry).read_text(encoding="utf-8"))
    report = audit_registry(registry)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
