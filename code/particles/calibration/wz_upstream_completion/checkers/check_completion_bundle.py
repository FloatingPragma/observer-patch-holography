#!/usr/bin/env python3
"""Fail-closed *specification linter* for the W/Z completion templates.

This module validates the fixed templates and evaluates their candidate
conjunction. It does not resolve external artifacts or derive evidence from
them. Consequently it is not a production aggregate verifier and can never
authorize scientific promotion. It deliberately imports no diagram generator,
self-energy engine, or target data.
"""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import json
from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]

SCHEMAS = {
    "action": "sm_eft_action_packet_v1.schema.json",
    "matching": "eft_matching_packet_v1.schema.json",
    "yukawa": "full_yukawa_packet_v1.schema.json",
    "fj": "fj_equivalence_receipt_v1.schema.json",
    "renorm": "renormalization_packet_v1.schema.json",
    "brst": "gauge_brst_raw_receipt_v1.schema.json",
    "pole_w": "physical_pole_interpretation_receipt_v1.schema.json",
    "pole_z": "physical_pole_interpretation_receipt_v1.schema.json",
    "law": "source_law_covariance_packet_v1.schema.json",
    "clock": "operational_clock_packet_v1.schema.json",
}
TEMPLATES = {
    "action": "sm_eft_action_packet_TEMPLATE.json",
    "matching": "eft_matching_packet_TEMPLATE.json",
    "yukawa": "full_yukawa_packet_TEMPLATE.json",
    "fj": "fj_equivalence_receipt_TEMPLATE.json",
    "renorm": "renormalization_packet_TEMPLATE.json",
    "brst": "gauge_brst_raw_receipt_TEMPLATE.json",
    "pole_w": "physical_pole_interpretation_TEMPLATE.json",
    "pole_z": "physical_pole_interpretation_Z_TEMPLATE.json",
    "law": "source_law_covariance_packet_TEMPLATE.json",
    "clock": "operational_clock_packet_TEMPLATE.json",
}

ZERO_HASH = "0" * 64
STATUS = "DRAFT_SUFFICIENCY_STACK_DEFINED__SIMULATION_RECEIPTS_OPEN__NO_OPH_NATIVE_POLE_PROMOTION"
VERIFIER_POLICY = "SPECIFICATION_ONLY__EXTERNAL_ARTIFACT_RESOLUTION_NOT_IMPLEMENTED"


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def validate_all():
    data = {}
    for key in SCHEMAS:
        schema = load(ROOT / "schemas" / SCHEMAS[key])
        obj = load(ROOT / "templates" / TEMPLATES[key])
        errors = sorted(
            Draft202012Validator(schema).iter_errors(obj),
            key=lambda e: list(e.path),
        )
        if errors:
            raise AssertionError(
                f"{key} schema: " + "; ".join(e.message for e in errors)
            )
        data[key] = obj
    return data


def _all_true(mapping: dict) -> bool:
    return all(value is True for value in mapping.values())


def _matching_reasons(m: dict) -> list[str]:
    reasons: list[str] = []
    if m["claim_lane"] != "OPH_NATIVE_PHYSICAL":
        reasons.append("matching lane not OPH-native physical")
    if not m["source_root_unique"]:
        reasons.append("matching source root not unique")
    for interval in m["intervals"]:
        if not interval["ordering_checked"] or not (
            interval["q_high_over_Estar"] > interval["q_low_over_Estar"]
        ):
            reasons.append(f"matching interval ordering invalid: {interval['id']}")
        if not (
            interval["beta"]["derived_from_census"]
            and interval["beta"]["independently_checked"]
        ):
            reasons.append(f"beta functions not independently census-derived: {interval['id']}")
        # The template names a pure one-Higgs, three-generation SM interval.
        # For that exact census, independently recompute the gauge coefficients.
        if interval["eft_id"] == "SM_ONE_HIGGS_THREE_GENERATION_TEMPLATE":
            got = interval["beta"]["gauge_coefficients"]
            expected = {"gprime": Fraction(41, 6), "g": Fraction(-19, 6), "gs": Fraction(-7, 1)}
            if any(Fraction(got[k]) != v for k, v in expected.items()):
                reasons.append("pure-SM interval has non-SM one-loop gauge coefficients")
        if not all(t["checked"] for t in interval["thresholds"]):
            reasons.append(f"unchecked threshold map: {interval['id']}")
    if not _all_true(m["evidence"]):
        reasons.append("EFT matching evidence incomplete")
    return reasons


def _law_reasons(law: dict) -> list[str]:
    reasons: list[str] = []
    if not law["active_mode_consistent"]:
        reasons.append("source-law active mode inconsistent")
    if law["mode"] == "deterministic_delta":
        required = law["deterministic"]
        if not _all_true({
            "global_root_unique": required["global_root_unique"],
            "selector_unique": required["selector_unique"],
            "primitives_exact": required["primitives_exact"],
            "matching_deterministic": required["matching_deterministic"],
            "covariance_zero_exact": required["covariance_zero_exact"],
        }):
            reasons.append("deterministic delta law lacks uniqueness/exactness proof")
    else:
        if not _all_true({k: v for k, v in law["stochastic"].items() if k != "ensemble_hash"}):
            reasons.append("stochastic source law/weights/covariance incomplete")
    if not _all_true(law["evidence"]):
        reasons.append("source-law/covariance evidence incomplete")
    return reasons


def _pole_reasons(p: dict, expected_boson: str) -> list[str]:
    reasons: list[str] = []
    if p["boson"] != expected_boson:
        reasons.append(f"{expected_boson} pole receipt has wrong boson tag")
    if p["current_amplitude"]["positive_residue_required"] is not False:
        reasons.append(f"{expected_boson} unstable pole incorrectly requires positive residue")
    simple = p["simple_zero"]
    if not (
        simple["rank_at_pole"] == simple["matrix_dimension"] - 1
        and simple["left_kernel_dimension"] == 1
        and simple["right_kernel_dimension"] == 1
    ):
        reasons.append(f"{expected_boson} pole lacks rank-n-minus-one Laurent hypothesis")
    if not (
        p["evidence"]["promotion_ready"]
        and p["evidence"]["independent_checker"]
        and p["current_amplitude"]["same_pole"]
        and p["current_amplitude"]["nonzero_vertex_contractions"]
        and p["simple_zero"]["derivative_excludes_zero"]
        and p["contour"]["rouche_passed"]
        and p["contour"]["boundary_nonzero"]
        and p["null_vectors"]["laurent_denominator_excludes_zero"]
    ):
        reasons.append(f"{expected_boson} physical-current pole interpretation incomplete")
    return reasons


def _common_hash_reasons(d: dict) -> list[str]:
    reasons: list[str] = []
    action_hashes = [
        d["action"]["hashes"]["action_ast"],
        d["matching"]["intervals"][0]["hashes"]["action_ast"],
        d["yukawa"]["hashes"]["action_ast"],
        d["fj"]["hashes"]["action_ast"],
        d["renorm"]["hashes"]["action_ast"],
        d["brst"]["hashes"]["action_ast"],
    ]
    nonzero = [h for h in action_hashes if h != ZERO_HASH]
    if nonzero and len(set(nonzero)) != 1:
        reasons.append("action hash mismatch across receipts")
    term_masks = [
        d["matching"]["intervals"][0]["beta"]["monomial_mask_hash"],
        d["fj"]["hashes"]["term_mask"],
        d["renorm"]["hashes"]["term_mask"],
        d["brst"]["hashes"]["term_mask"],
        d["pole_w"]["hashes"]["term_mask"],
        d["pole_z"]["hashes"]["term_mask"],
    ]
    nonzero = [h for h in term_masks if h != ZERO_HASH]
    if nonzero and len(set(nonzero)) != 1:
        reasons.append("term-mask hash mismatch across receipts")
    return reasons


def _contains_placeholder_hash(obj) -> bool:
    if isinstance(obj, dict):
        return any(_contains_placeholder_hash(v) for v in obj.values())
    if isinstance(obj, list):
        return any(_contains_placeholder_hash(v) for v in obj)
    return obj == ZERO_HASH


def promotion_reasons(d):
    """Return reasons the *declared candidate conjunction* is not satisfied.

    The values consumed here are unverified producer declarations. An empty
    list is useful only for testing predicate wiring; it is never promotion.
    """
    r: list[str] = []
    a, m, y, f, n, b = (
        d["action"], d["matching"], d["yukawa"], d["fj"], d["renorm"], d["brst"]
    )
    law, clock = d["law"], d["clock"]

    if a["claim_lane"] != "OPH_NATIVE_PHYSICAL":
        r.append("action lane not OPH-native physical")
    if a["geometry"]["origin"] != "oph_receipt":
        r.append("spacetime/action geometry is imported or unproved")
    if not _all_true(a["evidence"]):
        r.append("action parents incomplete")
    if not a["source_ancestry"]["target_blacklist_passed"]:
        r.append("target ancestry not cleared")

    r.extend(_matching_reasons(m))

    if not _all_true(y["evidence"]):
        r.append("full Yukawa/open-channel evidence incomplete")
    if y["approximation"]["mode"] != "full" and not y["approximation"]["remainder_bound_supplied"]:
        r.append("light-fermion approximation lacks remainder")

    if not (_all_true(f["map"]) and _all_true(f["checks"]) and _all_true(f["evidence"])):
        r.append("direct/converted FJ equivalence incomplete")
    if f["engines"]["shared_generated_expressions"] or f["engines"]["shared_integral_backend"]:
        r.append("FJ engines not independent")

    if not (
        n["evidence"]["promotion_ready"]
        and n["evidence"]["independent_checker"]
        and n["uv_poles"]["exact_cancellation"]
        and n["st_restoration"]["linearized_identity_passed"]
        and n["gamma5"]["finite_restoration_declared"]
    ):
        r.append("renormalization/ST incomplete")
    if n["counterterm_generation"]["handwritten_vertex_list"]:
        r.append("handwritten CT vertex list forbidden")
    if not (
        n["counterterm_generation"]["from_bare_substitution"]
        and n["counterterm_generation"]["external_brst_sources"]
    ):
        r.append("counterterms not generated from complete bare action")

    if not b["diagram_universe"]["complete"]:
        r.append("diagram universe incomplete")
    if not _all_true(b["brst"]):
        r.append("BRST construction incomplete")
    if not _all_true(b["identities"]):
        r.append("BRST/ST/Ward/Nielsen/FJ identities incomplete")
    if b["independence"]["shared_generated_expressions"] or b["independence"]["shared_integral_backend"]:
        r.append("pole engines not independent")
    if not _all_true({
        "raw": b["evidence"]["raw_diagrams_present"],
        "checker": b["evidence"]["small_checker_passed"],
        "mutations": b["evidence"]["mutation_suite_passed"],
    }):
        r.append("BRST raw/checker/mutation evidence incomplete")
    if not (b["precision"]["nested_balls"] and b["precision"]["zero_containment"]):
        r.append("complex-ball precision certification incomplete")

    r.extend(_pole_reasons(d["pole_w"], "W"))
    r.extend(_pole_reasons(d["pole_z"], "Z"))
    r.extend(_law_reasons(law))

    gap = clock["dimensionless_gap_over_Estar"]
    if not (gap["positive_gap_certified"] and gap["hi"] >= gap["lo"] > 0):
        r.append("source clock gap not positively certified")
    if not _all_true(clock["naturality"]):
        r.append("operational clock naturality/calibration incomplete")
    if not (clock["source_ancestry"]["target_blacklist_passed"] and _all_true(clock["evidence"])):
        r.append("source clock ancestry/evidence incomplete")

    r.extend(_common_hash_reasons(d))
    if any(_contains_placeholder_hash(section) for section in d.values()):
        r.append("placeholder hash present")
    return list(dict.fromkeys(r))


def aggregate_status(d):
    """Return the non-promoting status for this specification-only checker."""
    reasons = promotion_reasons(d)
    return {
        "status": STATUS,
        "verifier_policy": VERIFIER_POLICY,
        "schemas_valid": True,
        "distinct_schema_documents": len(set(SCHEMAS.values())),
        "receipt_instances": len(SCHEMAS),
        "candidate_conjunction_satisfied": not reasons,
        "production_artifacts_resolved": False,
        "production_digests_recomputed": False,
        "evidence_booleans_derived_by_verifier": False,
        "promotion_allowed": False,
        "failed_reasons": reasons,
        "promotion_blockers": [
            "specification-only checker",
            "external artifact resolver absent",
            "proof-bearing production schemas absent",
        ],
    }


def main():
    d = validate_all()
    out = aggregate_status(d)
    print(json.dumps(out, indent=2))
    (ROOT / "outputs" / "aggregate_status.json").write_text(
        json.dumps(out, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
