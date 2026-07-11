#!/usr/bin/env python3
"""Deterministically derive the hardened blind-analysis lock from a sealed catalog.

This module is pure preregistration plumbing: it verifies the public/private
catalog, classifies every dynamic and diagnostic circuit, derives exact latent
per-circuit hypothesis tables, and builds the self-hashed analysis lock.  It
never reads credentials, connects to IBM, or submits a job.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from itertools import permutations, product
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

import blind_preregister
import cayley_blind_likelihood_analysis as analysis
import generative_repair_kernel as kernels
from record_gated_cayley_circuits import build_recipe


DYNAMIC_SHOTS = 192
DIAGNOSTIC_SHOTS = 512
CONTAMINATION_PROBABILITY = 0.08
CONTAMINATION_SENSITIVITY = (0.04, 0.15)
LABEL_COMPONENT_COUNT = 256


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _uniform_table(valid_codes: Sequence[int]) -> dict[str, float]:
    vector = np.full(analysis.JOINT_CARDINALITY, 1.0 / analysis.JOINT_CARDINALITY)
    return analysis.vector_to_probability_mapping(vector, valid_codes, include_zeros=True)


def _joint_table(
    heated: int,
    decision: int,
    final: int,
    valid_codes: Sequence[int],
    *,
    acceptance_probability: float | None = None,
    rejected_final: int | None = None,
) -> dict[str, float]:
    if acceptance_probability is None:
        return {analysis.outcome_key(heated, decision, final, valid_codes): 1.0}
    probability = float(acceptance_probability)
    if not 0.0 <= probability <= 1.0 or rejected_final is None:
        raise ValueError("invalid MH acceptance table")
    table: dict[str, float] = {}
    if probability > 0.0:
        table[analysis.outcome_key(heated, 1, final, valid_codes)] = probability
    if probability < 1.0:
        table[analysis.outcome_key(heated, 0, rejected_final, valid_codes)] = 1.0 - probability
    return table


def _cayley_candidates(parameters: Mapping[str, Any]) -> tuple[dict[str, dict[str, float]], list[int]]:
    model_key = str(parameters["model"])
    model = kernels.builtin_cayley_models()[model_key]
    encoding = tuple(int(value) for value in parameters["state_encoding"])
    valid_codes = sorted(encoding)
    protocol_by_model = {
        analysis.REPAIR_MODEL: "record_gated",
        "lazy_heat": "open_loop_heat",
        "delayed_record": "delayed_record",
        "shuffled_record": "shuffled_record",
        "inverted_record": "inverted_record",
    }
    candidates: dict[str, dict[str, float]] = {}
    for hypothesis, protocol in protocol_by_model.items():
        recipe = build_recipe(
            model,
            protocol,
            int(parameters["initial_state"]),
            int(parameters["disturbance_slot"]),
            int(parameters["second_slot"]),
            encoding,
        )
        candidates[hypothesis] = _joint_table(
            encoding[recipe.heated_state],
            recipe.decision_record,
            encoding[recipe.final_state],
            valid_codes,
        )
    # Filled after every fixed-design stratum has been assembled.
    candidates[analysis.STATE_PREPARATION_MODEL] = dict(candidates[analysis.REPAIR_MODEL])
    candidates["calibrated_noise"] = _uniform_table(valid_codes)
    return candidates, valid_codes


def _mh_candidates(parameters: Mapping[str, Any]) -> tuple[dict[str, dict[str, float]], list[int]]:
    spectrum = kernels.builtin_spectra()[str(parameters["spectrum"])]
    permutation = tuple(int(value) for value in parameters["label_permutation"])
    source = int(parameters["semantic_source"])
    target = int(parameters["semantic_target"])
    physical_source = permutation[source]
    physical_target = permutation[target]
    valid_codes = sorted(permutation)
    candidates: dict[str, dict[str, float]] = {}
    for model_name, kappa in (("kappa_1", 1.0), ("kappa_0", 0.0), ("kappa_2", 2.0)):
        acceptance = kernels.metropolis_acceptance(
            spectrum,
            float(parameters["beta"]),
            kappa,
        )[source, target]
        candidates[model_name] = _joint_table(
            physical_source,
            1,
            physical_target,
            valid_codes,
            acceptance_probability=float(acceptance),
            rejected_final=physical_source,
        )
    candidates["mh_calibrated_noise"] = _uniform_table(valid_codes)
    return candidates, valid_codes


def _semantic_counterpart_key(descriptor: Mapping[str, Any]) -> str:
    parameters = dict(descriptor["parameters"])
    parameters.pop("state_encoding", None)
    parameters.pop("label_permutation", None)
    return analysis.sha256_json(
        {"family": descriptor["family"], "parameters": parameters}
    )


def _transform_table(
    probabilities: Mapping[str, Any],
    source_valid_codes: Sequence[int],
    target_valid_codes: Sequence[int],
    code_mapping: Mapping[int, int],
) -> dict[str, float]:
    source = analysis.probability_mapping_to_vector(probabilities, source_valid_codes)
    target = np.zeros_like(source)
    for index, probability in enumerate(source):
        if probability <= 0.0:
            continue
        heated, decision, final = analysis.joint_tuple(index)
        target[
            analysis.joint_index(
                code_mapping[heated],
                decision,
                code_mapping[final],
            )
        ] += probability
    return analysis.vector_to_probability_mapping(target, target_valid_codes)


def _endpoint_for_cayley(parameters: Mapping[str, Any]) -> tuple[str, str]:
    model_key = str(parameters["model"])
    protocol = str(parameters["protocol"])
    correct_model = {
        "record_gated": analysis.REPAIR_MODEL,
        "open_loop_heat": "lazy_heat",
        "delayed_record": "delayed_record",
        "shuffled_record": "shuffled_record",
        "inverted_record": "inverted_record",
    }[protocol]
    if model_key == "s3" and protocol == "record_gated":
        return analysis.PRIMARY_ENDPOINT, correct_model
    if model_key == "z5" and protocol == "record_gated":
        return "secondary_z5", correct_model
    return f"cayley_control_{model_key}_{protocol}", correct_model


def _candidate_provenance(
    opaque_id: str,
    descriptor: Mapping[str, Any],
    candidates: Mapping[str, Mapping[str, Any]],
    logical_circuit_sha256: str,
) -> dict[str, dict[str, str]]:
    derivations = {
        model: analysis.sha256_json(
            {
                "opaque_id": opaque_id,
                "family": descriptor["family"],
                "parameters": descriptor["parameters"],
                "hypothesis": model,
                "logical_circuit_sha256": logical_circuit_sha256,
                "analysis_source_sha256": analysis.analysis_code_sha256(),
            }
        )
        for model in candidates
    }
    return analysis.build_candidate_provenance(
        candidates,
        derivations,
        logical_circuit_sha256,
    )


def ideal_primary_label_power(lock: Mapping[str, Any]) -> dict[str, Any]:
    """Expected finite-shot LR power against the frozen global mapping mixture."""

    analysis.validate_analysis_lock(lock)
    rows = [row for row in lock["expected_rows"] if row["endpoint"] == analysis.PRIMARY_ENDPOINT]
    components = lock["label_layout_model"]["components"]
    component_ids = [component["component_id"] for component in components]
    weights = {component["component_id"]: float(component["prior_weight"]) for component in components}
    contributions: dict[str, list[dict[str, Any]]] = {}
    for calibration_id in sorted({row["calibration_id"] for row in rows}):
        selected_rows = [row for row in rows if row["calibration_id"] == calibration_id]
        variants = analysis._calibration_variants(lock["calibrations"][calibration_id])
        variant_contributions: list[dict[str, Any]] = []
        for channel in variants:
            repair_score = 0.0
            component_scores = {component_id: 0.0 for component_id in component_ids}
            for row in selected_rows:
                valid_codes = row["valid_codes"]
                latent_repair = analysis.probability_mapping_to_vector(
                    row["candidate_probabilities"][analysis.REPAIR_MODEL], valid_codes
                )
                true_probabilities = analysis.convolve_calibration(latent_repair, channel)
                shots = int(row["shots"])
                repair_score += shots * float(
                    np.sum(true_probabilities * np.log(true_probabilities))
                )
                for component in components:
                    latent_component = analysis.probability_mapping_to_vector(
                        component["row_probabilities"][row["row_id"]], valid_codes
                    )
                    predicted = analysis.convolve_calibration(latent_component, channel)
                    component_scores[component["component_id"]] += shots * float(
                        np.sum(true_probabilities * np.log(predicted))
                    )
            label_marginal = analysis._logsumexp(
                [
                    np.log(weights[component_id]) + component_scores[component_id]
                    for component_id in component_ids
                ]
            )
            variant_contributions.append(
                {
                    "repair_log_score": repair_score,
                    "component_log_scores": component_scores,
                    "label_log_marginal": label_marginal,
                    "log_likelihood_ratio": repair_score - label_marginal,
                }
            )
        contributions[calibration_id] = variant_contributions

    per_backend = {
        calibration_id: {
            "nominal_expected_log_likelihood_ratio": values[0]["log_likelihood_ratio"],
            "sensitivity_expected_log_likelihood_ratio_lower_bound": min(
                value["log_likelihood_ratio"] for value in values
            ),
        }
        for calibration_id, values in contributions.items()
    }
    pooled_values: list[float] = []
    calibration_ids = sorted(contributions)
    nominal_component_scores = {component_id: 0.0 for component_id in component_ids}
    nominal_repair = 0.0
    for calibration_id in calibration_ids:
        nominal = contributions[calibration_id][0]
        nominal_repair += nominal["repair_log_score"]
        for component_id in component_ids:
            nominal_component_scores[component_id] += nominal["component_log_scores"][component_id]
    for indices in product(
        *(range(len(contributions[calibration_id])) for calibration_id in calibration_ids)
    ):
        repair_score = 0.0
        component_scores = {component_id: 0.0 for component_id in component_ids}
        for calibration_id, index in zip(calibration_ids, indices):
            value = contributions[calibration_id][index]
            repair_score += value["repair_log_score"]
            for component_id in component_ids:
                component_scores[component_id] += value["component_log_scores"][component_id]
        marginal = analysis._logsumexp(
            [
                np.log(weights[component_id]) + component_scores[component_id]
                for component_id in component_ids
            ]
        )
        pooled_values.append(repair_score - marginal)
    ranked_nominal = sorted(
        nominal_component_scores.items(), key=lambda item: (-item[1], item[0])
    )
    reference_id = lock["label_layout_model"]["reference_component_id"]
    receipt = {
        "method": "finite-shot expected complete-joint conditional likelihood ratio",
        "primary_row_count": len(rows),
        "component_count": len(components),
        "reference_component_id": reference_id,
        "reference_component_is_unique_top": bool(
            ranked_nominal[0][0] == reference_id
            and ranked_nominal[0][1] > ranked_nominal[1][1]
        ),
        "nominal_top_to_second_log_score_gap": ranked_nominal[0][1] - ranked_nominal[1][1],
        "per_backend": per_backend,
        "pooled_nominal_expected_log_likelihood_ratio": pooled_values[0],
        "pooled_sensitivity_expected_log_likelihood_ratio_lower_bound": min(pooled_values),
    }
    receipt["receipt_sha256"] = analysis.sha256_json(receipt)
    return receipt


def build_locked_analysis(
    public_manifest: Mapping[str, Any],
    reveal: Mapping[str, Any],
) -> dict[str, Any]:
    """Return a complete hardened analysis lock for exactly one sealed catalog."""

    blind_preregister.verify_bundle(public_manifest, reveal, rebuild_circuits=False)
    payload = reveal["sealed_payload"]
    public_by_id = {row["opaque_id"]: row for row in public_manifest["circuits"]}
    private_by_id = payload["circuits"]
    if set(public_by_id) != set(private_by_id):
        raise ValueError("public/private circuit sets differ")

    private_slots_by_role = {
        slot["role_opaque_id"]: slot for slot in payload["backend_slots"]
    }
    rows: list[dict[str, Any]] = []
    diagnostics_by_role: dict[str, dict[str, int]] = defaultdict(dict)
    correct_model_by_endpoint: dict[str, str] = {}
    row_descriptor_by_id: dict[str, Mapping[str, Any]] = {}
    counterpart_index: dict[tuple[str, str], str] = {}

    for opaque_id in sorted(private_by_id):
        descriptor = private_by_id[opaque_id]
        public_row = public_by_id[opaque_id]
        role_id = str(descriptor["backend_role_opaque_id"])
        slot = private_slots_by_role[role_id]
        family = str(descriptor["family"])
        if family == "readout_calibration":
            if int(public_row["shots"]) != DIAGNOSTIC_SHOTS:
                raise ValueError("diagnostic calibration shot count changed")
            if descriptor["parameters"].get("evidentiary_role") != "diagnostic_only":
                raise ValueError("basis calibration is not marked diagnostic-only")
            diagnostics_by_role[role_id][opaque_id] = int(
                descriptor["parameters"]["basis_code"]
            )
            continue
        if int(public_row["shots"]) != DYNAMIC_SHOTS:
            raise ValueError("dynamic circuit shot count changed")

        parameters = descriptor["parameters"]
        if family == "cayley":
            candidates, valid_codes = _cayley_candidates(parameters)
            endpoint, correct_model = _endpoint_for_cayley(parameters)
            stratum_id = "sp_" + analysis.sha256_json(
                {
                    "backend_role": role_id,
                    "model": parameters["model"],
                    "protocol": parameters["protocol"],
                    "initial_state": parameters["initial_state"],
                }
            )[:32]
            family_fields = {"state_preparation_stratum_id": stratum_id}
        elif family == "mh":
            candidates, valid_codes = _mh_candidates(parameters)
            spectrum_key = str(parameters["spectrum"])
            endpoint = f"mh_{spectrum_key}"
            correct_model = analysis.MH_REFERENCE_MODEL
            family_fields = {
                "identifiability_role": (
                    "abelian_negative_control"
                    if spectrum_key in ("z3_cyclic_control", "z5_cyclic_control")
                    else "identifiable"
                )
            }
        else:
            raise ValueError(f"unsupported sealed circuit family {family!r}")
        previous_correct = correct_model_by_endpoint.setdefault(endpoint, correct_model)
        if previous_correct != correct_model:
            raise ValueError("one endpoint maps to multiple generating models")

        calibration_id = "cal_" + role_id
        row = {
            "row_id": opaque_id,
            "opaque_id": opaque_id,
            "logical_circuit_sha256": public_row["logical_circuit_sha256"],
            "family": family,
            "endpoint": endpoint,
            "backend_role": str(slot["role"]),
            "backend_name": str(slot["backend"]),
            "layout_id": str(slot["layout_opaque_id"]),
            "physical_layout": [int(value) for value in slot["layout"]],
            "calibration_id": calibration_id,
            "shots": DYNAMIC_SHOTS,
            "max_leakage_fraction": (
                analysis.INDIVIDUAL_CIRCUIT_CATASTROPHIC_LEAKAGE_FRACTION
            ),
            "valid_codes": valid_codes,
            "candidate_probabilities": candidates,
            **family_fields,
        }
        rows.append(row)
        row_descriptor_by_id[opaque_id] = descriptor
        counterpart_index[(role_id, _semantic_counterpart_key(descriptor))] = opaque_id

    if len(rows) != 3072 or sum(len(values) for values in diagnostics_by_role.values()) != 32:
        raise ValueError("production catalog must contain 3,072 dynamic and 32 diagnostic circuits")

    # Derive the fixed-design state-preparation-only null once per opaque stratum.
    strata: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row["family"] == "cayley":
            strata[row["state_preparation_stratum_id"]].append(row)
    for stratum_rows in strata.values():
        valid_codes = stratum_rows[0]["valid_codes"]
        aggregate = sum(
            analysis.probability_mapping_to_vector(
                row["candidate_probabilities"][analysis.REPAIR_MODEL], valid_codes
            )
            for row in stratum_rows
        ) / len(stratum_rows)
        state_null = analysis.state_preparation_only_joint_null(
            analysis.vector_to_probability_mapping(aggregate, valid_codes),
            valid_codes,
        )
        for row in stratum_rows:
            row["candidate_probabilities"][analysis.STATE_PREPARATION_MODEL] = state_null

    for row in rows:
        descriptor = row_descriptor_by_id[row["opaque_id"]]
        row["candidate_provenance"] = _candidate_provenance(
            row["opaque_id"],
            descriptor,
            row["candidate_probabilities"],
            row["logical_circuit_sha256"],
        )

    calibrations: dict[str, dict[str, Any]] = {}
    for role_id, expected_codes in diagnostics_by_role.items():
        descriptor_receipt = {
            opaque_id: {
                "basis_code": basis_code,
                "logical_circuit_sha256": public_by_id[opaque_id][
                    "logical_circuit_sha256"
                ],
                "shots": public_by_id[opaque_id]["shots"],
            }
            for opaque_id, basis_code in sorted(expected_codes.items())
        }
        protocol_hash = analysis.sha256_json(
            {"role": role_id, "diagnostics": descriptor_receipt}
        )
        derivation_hash = analysis.sha256_json(
            {
                "mode": "uniform_joint_contamination_prior",
                "nominal": CONTAMINATION_PROBABILITY,
                "sensitivity": list(CONTAMINATION_SENSITIVITY),
                "diagnostic_protocol_sha256": protocol_hash,
            }
        )
        calibrations["cal_" + role_id] = analysis.contamination_calibration_packet(
            contamination_probability=CONTAMINATION_PROBABILITY,
            sensitivity_probabilities=CONTAMINATION_SENSITIVITY,
            receipt_sha256=protocol_hash,
            derivation_sha256=derivation_hash,
            diagnostic_opaque_ids=sorted(expected_codes),
            diagnostic_expected_codes=expected_codes,
        )

    # Frozen global primary-S3 label alternatives.  The first 128 exact S3
    # state permutations are crossed with the two preregistered layout-role
    # assignments, giving 256 equal components.  Every map is a bijection
    # between the row's six valid physical codes; no label control manufactures
    # leakage.  Identity/no-swap is the exact true component.
    role_ids = sorted(private_slots_by_role)
    if len(role_ids) != 2:
        raise ValueError("label/layout mixture requires exactly two roles")
    other_role = {role_ids[0]: role_ids[1], role_ids[1]: role_ids[0]}
    components: list[dict[str, Any]] = []
    candidate_component_specs = [
        (swap_layout, permutation)
        for permutation in permutations(range(6))
        for swap_layout in (False, True)
    ]
    weight = 1.0 / LABEL_COMPONENT_COUNT
    rows_by_id = {row["row_id"]: row for row in rows}
    primary_rows_by_id = {
        row_id: row
        for row_id, row in rows_by_id.items()
        if row["endpoint"] == analysis.PRIMARY_ENDPOINT
    }
    if len(primary_rows_by_id) != 432:
        raise ValueError("primary S3 catalog must contain 432 circuits")
    selected_component_specs: list[tuple[bool, tuple[int, ...]]] = []
    component_table_hashes: set[str] = set()
    for swap_layout, permutation in candidate_component_specs:
        row_probabilities: dict[str, dict[str, float]] = {}
        row_derivations: dict[str, str] = {}
        for row_id, row in primary_rows_by_id.items():
            descriptor = row_descriptor_by_id[row_id]
            source_row = row
            source_descriptor = descriptor
            if swap_layout:
                role_id = str(descriptor["backend_role_opaque_id"])
                counterpart_id = counterpart_index[
                    (other_role[role_id], _semantic_counterpart_key(descriptor))
                ]
                source_row = rows_by_id[counterpart_id]
                source_descriptor = row_descriptor_by_id[counterpart_id]
            source_encoding = tuple(
                int(value) for value in source_descriptor["parameters"]["state_encoding"]
            )
            target_encoding = tuple(
                int(value) for value in descriptor["parameters"]["state_encoding"]
            )
            if swap_layout:
                source_order = sorted(source_encoding)
                target_order = sorted(target_encoding)
            else:
                source_order = list(source_encoding)
                target_order = list(target_encoding)
            code_mapping = {
                int(source_order[index]): int(target_order[permutation[index]])
                for index in range(6)
            }
            probabilities = _transform_table(
                source_row["candidate_probabilities"][analysis.REPAIR_MODEL],
                source_row["valid_codes"],
                row["valid_codes"],
                code_mapping,
            )
            row_probabilities[row_id] = probabilities
            row_derivations[row_id] = analysis.sha256_json(
                {
                    "row_id": row_id,
                    "source_row_id": source_row["row_id"],
                    "swap_layout": swap_layout,
                    "semantic_permutation": list(permutation),
                    "code_mapping": code_mapping,
                }
            )
        table_set_hash = analysis.sha256_json(row_probabilities)
        if table_set_hash in component_table_hashes:
            continue
        component_table_hashes.add(table_set_hash)
        selected_component_specs.append((swap_layout, permutation))
        component_id = "lbl_" + analysis.sha256_json(
            {
                "swap_layout": swap_layout,
                "semantic_permutation": list(permutation),
            }
        )[:32]
        components.append(
            analysis.build_label_layout_component(
                component_id=component_id,
                prior_weight=weight,
                row_probabilities=row_probabilities,
                component_derivation_sha256=analysis.sha256_json(
                    {
                        "component_id": component_id,
                        "swap_layout": swap_layout,
                        "semantic_permutation": list(permutation),
                        "row_derivations": row_derivations,
                    }
                ),
                row_derivation_sha256=row_derivations,
            )
        )
        if len(components) == LABEL_COMPONENT_COUNT:
            break
    if len(components) != LABEL_COMPONENT_COUNT:
        raise ValueError("could not construct 256 distinct global label/layout components")
    label_layout_model = {
        "mapping_scope": "global_shared_across_primary_rows",
        "component_set_derivation_sha256": analysis.sha256_json(
            {
                "component_specs": selected_component_specs,
                "catalog_precommitment_sha256": public_manifest[
                    "catalog_precommitment_sha256"
                ],
            }
        ),
        "components": components,
        "reference_component_id": components[0]["component_id"],
    }

    rows_by_endpoint: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        if row["endpoint"] != analysis.PRIMARY_ENDPOINT:
            rows_by_endpoint[row["endpoint"]].append(row["row_id"])
    secondary_tests = [
        {
            "test_id": "gof_" + analysis.sha256_json(endpoint)[:24],
            "row_ids": sorted(row_ids),
            "model": correct_model_by_endpoint[endpoint],
        }
        for endpoint, row_ids in sorted(rows_by_endpoint.items())
    ]

    lock = analysis.build_analysis_lock(
        expected_rows=sorted(rows, key=lambda row: row["row_id"]),
        calibrations=calibrations,
        secondary_tests=secondary_tests,
        catalog_precommitment_sha256=str(public_manifest["catalog_precommitment_sha256"]),
        label_layout_model=label_layout_model,
    )
    covered = set(lock["catalog_coverage"]["dynamic_analysis_opaque_ids"]) | set(
        lock["catalog_coverage"]["diagnostic_calibration_opaque_ids"]
    )
    if covered != set(public_by_id):
        raise ValueError("analysis lock does not classify the complete manifest catalog")
    if len(covered) != 3104:
        raise ValueError("production manifest coverage must contain exactly 3,104 circuits")
    power = ideal_primary_label_power(lock)
    if not power["reference_component_is_unique_top"]:
        raise ValueError("ideal power preflight does not uniquely recover the reference mapping")
    if any(
        details["sensitivity_expected_log_likelihood_ratio_lower_bound"]
        <= np.log(analysis.PER_BACKEND_LR_THRESHOLD)
        for details in power["per_backend"].values()
    ):
        raise ValueError("ideal per-backend label-mixture power does not exceed the frozen gate")
    if power["pooled_sensitivity_expected_log_likelihood_ratio_lower_bound"] <= np.log(
        analysis.POOLED_LR_THRESHOLD
    ):
        raise ValueError("ideal pooled label-mixture power does not exceed the frozen gate")
    body = dict(lock)
    body.pop("analysis_lock_sha256")
    body["ideal_primary_label_power"] = power
    lock = dict(body)
    lock["analysis_lock_sha256"] = analysis.sha256_json(body)
    analysis.validate_analysis_lock(lock)
    return lock


def main() -> int:
    raise SystemExit(
        "This module is invoked by blind_preregister_orchestrate.py; it does not access IBM."
    )


if __name__ == "__main__":
    main()
