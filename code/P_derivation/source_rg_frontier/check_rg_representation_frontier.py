#!/usr/bin/env python3
"""Independent resolver for the issue #32 RG representation frontier.

This checker imports no producer code.  It resolves every committed source
path, re-hashes every input, reconstructs the representation indices from the
finite matter data, recomputes each conditional coefficient row, verifies the
non-selection witnesses, and refuses any claim that matching objects exist.
"""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator


PACKAGE = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE.parents[2]
OUTPUT = PACKAGE / "outputs" / "rg_representation_frontier.json"
SCHEMA_PATH = PACKAGE / "schemas" / "rg_representation_frontier_v1.schema.json"
POLICY_PATH = PACKAGE / "data" / "source_rg_policy_v1.json"


def fail(code: str, message: str) -> None:
    raise SystemExit(f"{code}: {message}")


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        fail(code, message)


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def frac(value: Fraction | int) -> str:
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def resolve_pin(pin: dict[str, Any], *, json_file: bool = True) -> Any:
    rel = pin["path"]
    path = (REPO_ROOT / rel).resolve()
    require(path.is_relative_to(REPO_ROOT.resolve()), "PATH_TRAVERSAL", f"path escapes repository: {rel}")
    require(path.is_file(), "PIN_MISSING", f"pinned file missing: {rel}")
    raw = path.read_bytes()
    require(len(raw) == pin["bytes"], "PIN_SIZE", f"byte size changed: {rel}")
    require(hashlib.sha256(raw).hexdigest() == pin["byte_sha256"], "PIN_HASH", f"byte hash changed: {rel}")
    if not json_file:
        return raw
    parsed = json.loads(raw)
    require(
        canonical_sha256(parsed) == pin["canonical_json_sha256"],
        "PIN_CANONICAL_HASH",
        f"canonical JSON hash changed: {rel}",
    )
    return parsed


def coefficients(n_g: int, n_h: int) -> dict[str, str]:
    return {
        "b_Y": frac(Fraction(20, 9) * n_g + Fraction(1, 6) * n_h),
        "b_2": frac(-Fraction(22, 3) + Fraction(4, 3) * n_g + Fraction(1, 6) * n_h),
        "b_3": frac(-11 + Fraction(4, 3) * n_g),
    }


def rederive_indices(matter: dict[str, Any], family: dict[str, Any]) -> dict[str, Fraction]:
    fields = matter["realized_package"]["fields"]
    states = family["generation"]["states"]
    require({row["label"] for row in states} == set(fields), "FIELD_LABELS", "field labels drifted")
    su3 = Fraction(0)
    su2 = Fraction(0)
    u1 = Fraction(0)
    for row in states:
        name = row["label"]
        color = int(row["color"])
        weak = int(row["weak"])
        charge = Fraction(row["hypercharge"])
        multiplicity = int(row["weyl_states"])
        require(fields[name]["dimension"] == multiplicity, "FIELD_DIMENSION", f"{name} dimension drifted")
        require(Fraction(fields[name]["charge"]) == charge, "FIELD_CHARGE", f"{name} charge drifted")
        require(color * weak == multiplicity, "REP_DIMENSION", f"{name} representation dimension drifted")
        if color == 3:
            su3 += Fraction(1, 2) * weak
        else:
            require(color == 1, "COLOR_REP", f"unknown color dimension for {name}")
        if weak == 2:
            su2 += Fraction(1, 2) * color
        else:
            require(weak == 1, "WEAK_REP", f"unknown weak dimension for {name}")
        u1 += multiplicity * charge * charge
    return {"su3": su3, "su2": su2, "u1": u1}


def check(payload: dict[str, Any]) -> None:
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    errors = sorted(Draft202012Validator(schema).iter_errors(payload), key=lambda err: list(err.path))
    require(not errors, "SCHEMA", errors[0].message if errors else "unknown schema error")

    require(payload["promotion_allowed"] is False, "PROMOTION", "frontier must never promote")
    require(payload["issue"] == 32, "ISSUE", "wrong issue")
    require(
        payload["status"] == "PARTIAL_EXACT_REPRESENTATION_INDICES__SOURCE_MATCHING_OPEN",
        "STATUS",
        "frontier must remain partial and open",
    )

    policy = resolve_pin(payload["source_inputs"]["policy"])
    require(policy == json.loads(POLICY_PATH.read_text(encoding="utf-8")), "POLICY", "resolved policy differs")
    policy_by_role = {entry["role"]: entry for entry in policy["source_inputs"]}
    rows = payload["source_inputs"]["artifacts"]
    require({row["role"] for row in rows} == set(policy_by_role), "SOURCE_ROLES", "source role set drifted")
    resolved: dict[str, dict[str, Any]] = {}
    for row in rows:
        spec = policy_by_role[row["role"]]
        require(row["path"] == spec["path"], "SOURCE_PATH", f"path mismatch for {row['role']}")
        require(row["status"] == spec["status"], "SOURCE_STATUS", f"status mismatch for {row['role']}")
        resolved[row["role"]] = resolve_pin(row)
    expected_input_digest = canonical_sha256(
        {
            "policy": payload["source_inputs"]["policy"],
            "artifacts": rows,
        }
    )
    require(
        payload["source_inputs"]["input_bundle_digest"] == expected_input_digest,
        "INPUT_DIGEST",
        "input bundle digest mismatch",
    )

    forbidden_paths = set(policy["target_firewall"]["forbidden_source_paths"])
    require(
        all(row["path"] not in forbidden_paths for row in rows),
        "TARGET_PATH",
        "target or external-validation path entered the source bundle",
    )
    forbidden_tokens = [
        token.casefold() for token in policy["target_firewall"]["forbidden_structured_tokens"]
    ]
    source_text = json.dumps(resolved, sort_keys=True).casefold()
    hits = [token for token in forbidden_tokens if token in source_text]
    require(not hits, "TARGET_TOKEN", f"target token entered source inputs: {hits}")
    require(payload["source_inputs"]["target_paths_present"] is False, "TARGET_FLAG", "target flag must be false")

    global_form = resolved["global_charge_lattice"]
    require(global_form["hypercharge_convention"] == "q = 6Y", "HYPERCHARGE", "q=6Y convention drifted")
    matter = resolved["finite_matter_representation"]
    family = resolved["conditional_rank_three_screen_context"]
    indices = rederive_indices(matter, family)
    require(indices == {"su3": Fraction(2), "su2": Fraction(2), "u1": Fraction(10, 3)}, "INDICES", "indices drifted")
    stated = payload["representation_indices"]["matter_copy"]
    require(Fraction(stated["sum_weyl_T_SU3"]) == indices["su3"], "SU3_INDEX", "SU(3) index mismatch")
    require(Fraction(stated["sum_weyl_T_SU2"]) == indices["su2"], "SU2_INDEX", "SU(2) index mismatch")
    require(Fraction(stated["sum_weyl_dim3_dim2_Y2"]) == indices["u1"], "U1_INDEX", "U(1) index mismatch")
    require(stated["rank"] == 15, "MATTER_RANK", "matter rank must be fifteen")
    require(stated["charge_conjugation_invariant"] is True, "CONJUGATION", "quadratic indices must be conjugation invariant")

    scalar = payload["representation_indices"]["conditional_scalar_doublet"]
    require(
        scalar["status"] == "representation_exact_if_present__existence_and_count_not_selected",
        "SCALAR_STATUS",
        "scalar must stay conditional",
    )
    require(Fraction(scalar["sum_complex_scalar_T_SU2"]) == Fraction(1, 2), "SCALAR_SU2", "scalar SU(2) index mismatch")
    require(Fraction(scalar["sum_complex_scalar_dim3_dim2_Y2"]) == Fraction(1, 2), "SCALAR_U1", "scalar U(1) index mismatch")

    qft = payload["qft_import_boundary"]
    require(qft["status"] == "imported_not_derived_from_A1_A3", "QFT_IMPORT", "QFT functional must remain imported")
    require(qft["oph_native_one_loop_beta_theorem"] is False, "QFT_PROMOTION", "native QFT theorem is absent")
    require(qft["external_593_packet_consumed"] is False, "EXTERNAL_593", "#593 output must not enter source ancestry")

    law = payload["parametric_one_loop_law"]
    require(
        law["convention"] == "d g_i/d ln(mu)=b_i g_i^3/(16 pi^2)",
        "BETA_CONVENTION",
        "one-loop sign convention mismatch",
    )
    require(
        law["hypercharge_normalization"] == "SM gprime convention with q=6Y",
        "HYPERCHARGE_NORMALIZATION",
        "hypercharge normalization mismatch",
    )
    require(law["coefficients"] == {
        "b_Y": "(20/9) N_g + (1/6) N_H",
        "b_2": "-22/3 + (4/3) N_g + (1/6) N_H",
        "b_3": "-11 + (4/3) N_g",
    }, "BETA_LAW", "parametric beta law mismatch")

    invisible = payload["invisible_sector_boundary"]
    require(
        invisible["exact_result"]
        == {
            "delta_b_Y": "0",
            "delta_b_2": "0",
            "delta_b_3": "0",
            "status": "proved_after_imported_one_loop_gauge_functional",
        },
        "STERILE_GAUGE_SHIFT",
        "zero-index direct-sum gauge shift drifted",
    )
    require(
        invisible["full_WZ_decoupling_proved"] is False,
        "STERILE_OVERPROMOTION",
        "gauge-index invariance is not a full W/Z decoupling theorem",
    )
    require(
        invisible["status"]
        == "GAUGE_INDEX_INVARIANCE_PROVED__FULL_WZ_DECOUPLING_OPEN",
        "STERILE_STATUS",
        "invisible-sector boundary status drifted",
    )
    require(
        "zero Yukawa and scalar vertices"
        in invisible["missing_for_full_WZ_decoupling"][0],
        "STERILE_VERTEX_BOUNDARY",
        "full decoupling must retain the open vertex condition",
    )
    conditional = payload["conditional_evaluations"]
    require(len(conditional) == 1, "CONDITIONAL_COUNT", "one declared evaluation expected")
    row = conditional[0]
    require(row["N_g"] == 3 and row["N_H"] == 1, "CONDITIONAL_INPUT", "declared evaluation changed")
    require(row["coefficients"] == coefficients(3, 1), "CONDITIONAL_BETA", "declared evaluation arithmetic mismatch")
    require(row["promotion_allowed"] is False, "CONDITIONAL_PROMOTION", "declared evaluation must not promote")
    require(
        row["status"] == "conditional_declared_completion_not_OPH_selected",
        "CONDITIONAL_STATUS",
        "declared completion must stay unselected",
    )

    multiplicity = resolved["mandatory_multiplicity_countermodels"]
    require(
        multiplicity["family_multiplicity_window"]["verdict"]["count_inside_window"]
        == "not_source_selected",
        "FAMILY_SELECTION",
        "family count is not selected",
    )
    require(
        multiplicity["scalar_response_multiplicity"]["verdict"]["scalar_existence"]
        == "not_source_determined",
        "SCALAR_SELECTION",
        "scalar existence is not selected",
    )
    census = payload["nonidentifiability_witnesses"]["census_nonuniqueness"]
    expected_pairs = [(3, 0), (3, 1), (3, 2), (4, 1), (5, 1)]
    actual_pairs = [(item["N_g"], item["N_H"]) for item in census["witnesses"]]
    require(actual_pairs == expected_pairs, "CENSUS_WITNESSES", "countermodel grid drifted")
    for item in census["witnesses"]:
        require(
            item["coefficients"] == coefficients(item["N_g"], item["N_H"]),
            "CENSUS_BETA",
            f"countermodel arithmetic mismatch for {(item['N_g'], item['N_H'])}",
        )
        require(
            item["physical_selection_status"] == "not_source_selected",
            "CENSUS_PROMOTION",
            "countermodel row was promoted",
        )
    require(
        len({canonical_sha256(item["coefficients"]) for item in census["witnesses"]}) == 5,
        "CENSUS_FIBER",
        "countermodel beta vectors must be distinct",
    )

    witness_blocks = payload["nonidentifiability_witnesses"]
    parent_digests = {
        block["parent_projection_digest"]
        for block in witness_blocks.values()
    }
    require(len(parent_digests) == 1, "PARENT_PROJECTION", "witnesses must share one finite parent projection")
    schemes = witness_blocks["coordinate_nonuniqueness"]["finite_scheme_redefinitions"]
    require([row["c"] for row in schemes] == ["0", "1/8"], "SCHEME_WITNESSES", "scheme witnesses drifted")
    require(
        witness_blocks["coordinate_nonuniqueness"]["one_loop_coefficient_invariant"] is True,
        "SCHEME_ONE_LOOP",
        "finite scheme witness must preserve the one-loop coefficient",
    )
    thresholds = witness_blocks["threshold_nonuniqueness"]
    require(thresholds["physical_threshold_claim"] is False, "THRESHOLD_PROMOTION", "abstract witness is not a threshold claim")
    require(
        thresholds["abstract_extension_witnesses"][0]["dimensionless_mass_parameters"]
        != thresholds["abstract_extension_witnesses"][1]["dimensionless_mass_parameters"],
        "THRESHOLD_WITNESSES",
        "threshold extension witnesses must differ",
    )

    required_matching = {
        "ordered_eft_intervals",
        "threshold_locations",
        "decoupling_maps",
        "scheme_maps",
        "jacobians",
        "finite_order_term_masks",
        "certified_vector_remainders",
    }
    require(set(payload["matching_objects"]) == required_matching, "MATCHING_OBJECTS", "matching object set drifted")
    for name, item in payload["matching_objects"].items():
        require(item["status"] == "not_emitted", "MATCHING_PROMOTION", f"{name} must remain not emitted")
        require("reason" in item and item["reason"], "MATCHING_REASON", f"{name} lacks an explicit blocker")

    formal = payload["formal_certificate"]
    lean_path = (REPO_ROOT / formal["path"]).resolve()
    lean_raw = lean_path.read_bytes()
    require(len(lean_raw) == formal["bytes"], "LEAN_SIZE", "Lean certificate size drifted")
    require(hashlib.sha256(lean_raw).hexdigest() == formal["byte_sha256"], "LEAN_HASH", "Lean certificate hash drifted")
    lean_text = lean_raw.decode("utf-8")
    for theorem in formal["theorems"]:
        require(f"theorem {theorem}" in lean_text, "LEAN_THEOREM", f"Lean theorem missing: {theorem}")

    acceptance = {row["criterion"]: row for row in payload["acceptance_map"]}
    require(
        acceptance["target_clean_source_emits_complete_matching_packet"]["status"] == "open",
        "ACCEPTANCE_OVERCLAIM",
        "complete matching acceptance must stay open",
    )
    require(
        all(row["status"] != "complete" for row in payload["acceptance_map"]),
        "ACCEPTANCE_COMPLETE",
        "no #32 acceptance criterion is complete",
    )

    digest_payload = deepcopy(payload)
    stated_digest = digest_payload.pop("subject_digest")
    require(canonical_sha256(digest_payload) == stated_digest, "SUBJECT_DIGEST", "subject digest mismatch")


def main() -> int:
    require(OUTPUT.is_file(), "OUTPUT_MISSING", f"frontier output missing: {OUTPUT}")
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    check(payload)
    print(
        "RG representation frontier OK: finite source pins, exact indices, "
        "parametric one-loop law, non-selection witnesses, target firewall, "
        "Lean binding, and open matching gates all verified"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
