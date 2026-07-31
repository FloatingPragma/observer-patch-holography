#!/usr/bin/env python3
"""Build the exact, non-promoting representation frontier for issue #32.

The producer reads only finite OPH source receipts.  It derives exact
representation indices and then applies one explicitly imported one-loop QFT
functional.  Family and scalar multiplicities remain variables.  No EFT
interval, mass threshold, decoupling map, scheme conversion, Jacobian, term
mask, or remainder enclosure is fabricated.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import re
from pathlib import Path
from typing import Any


PACKAGE = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE.parents[2]
POLICY_PATH = PACKAGE / "data" / "source_rg_policy_v1.json"
DEFAULT_OUTPUT = PACKAGE / "outputs" / "rg_representation_frontier.json"
SCHEMA = "oph.rg_representation_frontier.v1"
LEAN_PATH = REPO_ROOT / "Lean" / "Screen" / "RGRepresentationFrontier.lean"


class FrontierError(ValueError):
    """A fail-closed producer error with a stable machine code."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise FrontierError(code, message)


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


def file_pin(path: Path, *, role: str | None = None, status: str | None = None) -> dict[str, Any]:
    raw = path.read_bytes()
    parsed = json.loads(raw)
    pin: dict[str, Any] = {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "bytes": len(raw),
        "byte_sha256": hashlib.sha256(raw).hexdigest(),
        "canonical_json_sha256": canonical_sha256(parsed),
    }
    if role is not None:
        pin["role"] = role
    if status is not None:
        pin["status"] = status
    return pin


def resolve_source_inputs(policy: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    allowed_roots = tuple(policy["target_firewall"]["allowed_source_roots"])
    forbidden_paths = set(policy["target_firewall"]["forbidden_source_paths"])
    rows: list[dict[str, Any]] = []
    payloads: dict[str, dict[str, Any]] = {}
    for entry in policy["source_inputs"]:
        rel = entry["path"]
        require(rel not in forbidden_paths, "TARGET_PATH", f"forbidden source path {rel}")
        require(
            any(rel == root or rel.startswith(root + "/") for root in allowed_roots),
            "PATH_ALLOWLIST",
            f"source path outside the closed allowlist: {rel}",
        )
        path = (REPO_ROOT / rel).resolve()
        require(
            path.is_relative_to(REPO_ROOT.resolve()),
            "PATH_TRAVERSAL",
            f"source path escapes the repository: {rel}",
        )
        require(path.is_file(), "SOURCE_MISSING", f"source artifact is missing: {rel}")
        payload = json.loads(path.read_text(encoding="utf-8"))
        payloads[entry["role"]] = payload
        rows.append(
            file_pin(
                path,
                role=entry["role"],
                status=entry["status"],
            )
        )
    require(len(rows) == 4, "SOURCE_COUNT", "exactly four finite source artifacts are required")
    return rows, payloads


def target_firewall(source_rows: list[dict[str, Any]], source_payloads: dict[str, dict[str, Any]], policy: dict[str, Any]) -> None:
    forbidden_paths = set(policy["target_firewall"]["forbidden_source_paths"])
    require(
        all(row["path"] not in forbidden_paths for row in source_rows),
        "TARGET_PATH",
        "target or external-validation artifact entered the source bundle",
    )
    forbidden_tokens = tuple(
        token.casefold() for token in policy["target_firewall"]["forbidden_structured_tokens"]
    )
    structured = json.dumps(source_payloads, sort_keys=True).casefold()
    hits = [
        token
        for token in forbidden_tokens
        if re.search(
            rf"(?<![a-z0-9_]){re.escape(token)}(?![a-z0-9_])",
            structured,
        )
    ]
    require(not hits, "TARGET_TOKEN", f"forbidden structured target tokens found: {hits}")


def representation_indices(
    matter: dict[str, Any],
    family_context: dict[str, Any],
    global_form: dict[str, Any],
) -> dict[str, Any]:
    require(matter.get("issue") == 314, "MATTER_ISSUE", "matter receipt must be issue #314")
    require(
        matter.get("schema") == "oph.super_tannakian_matter_receipt.v5",
        "MATTER_SCHEMA",
        "unexpected matter receipt schema",
    )
    require(
        matter.get("conditional_algebraic_gate", {}).get("passed") is True
        and matter.get("physical_source_gate", {}).get("passed") is False
        and matter.get("physical_source_gate", {}).get(
            "upstream_current_representation_source_bound"
        )
        is False
        and matter.get("physical_source_gate", {}).get("matter_lift_source_bound")
        is False,
        "MATTER_SCOPE",
        "the matter representation must remain a declared conditional fixture",
    )
    require(
        global_form.get("hypercharge_convention") == "q = 6Y",
        "HYPERCHARGE_CONVENTION",
        "global-form packet must freeze q=6Y",
    )
    require(
        "do not select a physical global form"
        in str(global_form.get("description", "")).lower(),
        "GLOBAL_FORM_SCOPE",
        "the charge lattice must retain its open physical-global-form boundary",
    )
    family_scope = family_context.get("conditional_structural_scope", {})
    require(
        family_scope.get("physical_current_source_bound") is False
        and family_scope.get("physical_matter_lift_source_bound") is False
        and family_scope.get("physical_global_form_selected") is False,
        "FAMILY_SCOPE",
        "the rank-three context must retain its conditional structural scope",
    )

    source_fields = matter["realized_package"]["fields"]
    states = family_context["generation"]["states"]
    require(len(states) == 5, "FIELD_COUNT", "five irreducible matter blocks are required")
    require(
        {row["label"] for row in states} == set(source_fields),
        "FIELD_LABELS",
        "conditional context must mirror the five #314 matter blocks",
    )

    su3_sum = Fraction(0)
    su2_sum = Fraction(0)
    u1_sum = Fraction(0)
    normalized_states: list[dict[str, Any]] = []
    for state in sorted(states, key=lambda row: row["label"]):
        label = state["label"]
        color = int(state["color"])
        weak = int(state["weak"])
        charge = Fraction(state["hypercharge"])
        multiplicity = int(state["weyl_states"])
        require(
            source_fields[label]["dimension"] == multiplicity,
            "FIELD_DIMENSION",
            f"{label} dimension disagrees with the #314 matter receipt",
        )
        require(
            Fraction(source_fields[label]["charge"]) == charge,
            "FIELD_CHARGE",
            f"{label} hypercharge disagrees with the #314 matter receipt",
        )
        require(
            multiplicity == color * weak,
            "REPRESENTATION_DIMENSION",
            f"{label} multiplicity is not color times weak dimension",
        )
        if color == 3:
            su3_sum += Fraction(1, 2) * weak
        else:
            require(color == 1, "COLOR_REP", f"unsupported color dimension for {label}")
        if weak == 2:
            su2_sum += Fraction(1, 2) * color
        else:
            require(weak == 1, "WEAK_REP", f"unsupported weak dimension for {label}")
        u1_sum += multiplicity * charge * charge
        normalized_states.append(
            {
                "label": label,
                "color_dimension": color,
                "weak_dimension": weak,
                "hypercharge": frac(charge),
                "weyl_multiplicity": multiplicity,
            }
        )

    require(su3_sum == 2, "SU3_INDEX", f"expected per-copy SU(3) index 2, got {su3_sum}")
    require(su2_sum == 2, "SU2_INDEX", f"expected per-copy SU(2) index 2, got {su2_sum}")
    require(u1_sum == Fraction(10, 3), "U1_INDEX", f"expected U(1) index 10/3, got {u1_sum}")

    scalar = matter["scalar_and_channel_selection"]
    require(
        scalar["admissible_scalar_charges"] == [3, -3],
        "SCALAR_CHARGE",
        "conditional scalar charges must be the conjugate q=+/-3 pair",
    )
    require(
        "not derived" in scalar["scalar_content_status"],
        "SCALAR_STATUS",
        "scalar multiplicity must remain explicitly non-derived",
    )

    return {
        "matter_copy": {
            "states": normalized_states,
            "rank": sum(row["weyl_multiplicity"] for row in normalized_states),
            "sum_weyl_T_SU3": frac(su3_sum),
            "sum_weyl_T_SU2": frac(su2_sum),
            "sum_weyl_dim3_dim2_Y2": frac(u1_sum),
            "charge_conjugation_invariant": True,
            "status": "exact_finite_representation_index",
        },
        "conditional_scalar_doublet": {
            "color_dimension": 1,
            "weak_dimension": 2,
            "absolute_hypercharge": "1/2",
            "sum_complex_scalar_T_SU3": "0",
            "sum_complex_scalar_T_SU2": "1/2",
            "sum_complex_scalar_dim3_dim2_Y2": "1/2",
            "status": "representation_exact_if_present__existence_and_count_not_selected",
        },
    }


def beta_coefficients(n_g: int, n_h: int) -> dict[str, str]:
    require(n_g >= 0 and n_h >= 0, "NEGATIVE_MULTIPLICITY", "multiplicities must be nonnegative")
    b_y = Fraction(20, 9) * n_g + Fraction(1, 6) * n_h
    b_2 = -Fraction(22, 3) + Fraction(4, 3) * n_g + Fraction(1, 6) * n_h
    b_3 = -11 + Fraction(4, 3) * n_g
    return {"b_Y": frac(b_y), "b_2": frac(b_2), "b_3": frac(b_3)}


def parametric_law() -> dict[str, Any]:
    return {
        "convention": "d g_i/d ln(mu)=b_i g_i^3/(16 pi^2)",
        "hypercharge_normalization": "SM gprime convention with q=6Y",
        "imported_functional": "b=-(11/3)C2(G)+(2/3)sum_Weyl T(R)+(1/3)sum_complex_scalar T(R)",
        "variables": {
            "N_g": "number of rank-fifteen matter copies; source selection open",
            "N_H": "number of complex color-singlet weak doublets with |Y|=1/2; source selection open",
        },
        "coefficients": {
            "b_Y": "(20/9) N_g + (1/6) N_H",
            "b_2": "-22/3 + (4/3) N_g + (1/6) N_H",
            "b_3": "-11 + (4/3) N_g",
        },
        "opposite_sign_warning": "Some OPH receipts print B_i=-b_i when testing asymptotic freedom.",
        "gut_normalization": {
            "relation": "b_1_GUT=(3/5)b_Y",
            "status": "derived_conversion_not_the_active_convention",
        },
        "status": "exact_after_imported_one_loop_QFT_functional",
    }


def multiplicity_witnesses(multiplicity: dict[str, Any]) -> dict[str, Any]:
    family = multiplicity["family_multiplicity_window"]
    scalar = multiplicity["scalar_response_multiplicity"]
    family_counts = [row["n_g"] for row in family["in_window_non_selection"]["members"]]
    scalar_counts = [
        row["copy_count"] for row in scalar["countermodel_battery"]["configurations"]
    ]
    require(family_counts == [3, 4, 5], "FAMILY_COUNTERMODELS", "expected family witnesses 3,4,5")
    require(scalar_counts == [0, 2, 2], "SCALAR_COUNTERMODELS", "expected zero and two-scalar witnesses")
    require(
        family["verdict"]["count_inside_window"] == "not_source_selected",
        "FAMILY_SELECTION",
        "family count must remain unselected",
    )
    require(
        scalar["verdict"]["scalar_existence"] == "not_source_determined",
        "SCALAR_SELECTION",
        "scalar existence must remain unselected",
    )

    census_rows = []
    for n_g, n_h, witness in (
        (3, 0, "zero_scalar_completion"),
        (3, 1, "declared_one_doublet_completion"),
        (3, 2, "duplicate_or_inert_doublet_completion"),
        (4, 1, "four_family_reduct"),
        (5, 1, "five_family_reduct"),
    ):
        census_rows.append(
            {
                "N_g": n_g,
                "N_H": n_h,
                "witness": witness,
                "coefficients": beta_coefficients(n_g, n_h),
                "physical_selection_status": "not_source_selected",
            }
        )
    require(
        len({canonical_sha256(row["coefficients"]) for row in census_rows}) == len(census_rows),
        "COEFFICIENT_FIBER",
        "countermodel coefficient rows must be pairwise distinct",
    )

    parent_projection = {
        "finite_matter_module": "fixed",
        "global_charge_lattice": "fixed",
        "mass_coordinates_present": False,
        "renormalization_scheme_present": False,
    }
    parent_digest = canonical_sha256(parent_projection)
    return {
        "census_nonuniqueness": {
            "parent_projection_digest": parent_digest,
            "witnesses": census_rows,
            "verdict": "current finite reduct does not select a complete beta vector",
            "scope": "registered family and scalar multiplicity interfaces only",
        },
        "coordinate_nonuniqueness": {
            "parent_projection_digest": parent_digest,
            "finite_scheme_redefinitions": [
                {
                    "id": "identity",
                    "map": "g'=g+0*g^3+O(g^5)",
                    "c": "0",
                    "jacobian_at_origin": "1",
                },
                {
                    "id": "cubic_shift",
                    "map": "g'=g+(1/8)*g^3+O(g^5)",
                    "c": "1/8",
                    "jacobian_at_origin": "1",
                },
            ],
            "one_loop_coefficient_invariant": True,
            "higher_order_coefficients_and_finite_matching_not_selected": True,
            "verdict": "finite parents do not select scheme maps or higher-loop coordinates",
        },
        "threshold_nonuniqueness": {
            "parent_projection_digest": parent_digest,
            "abstract_extension_witnesses": [
                {
                    "id": "mass_assignment_A",
                    "dimensionless_mass_parameters": ["m", "2*m"],
                },
                {
                    "id": "mass_assignment_B",
                    "dimensionless_mass_parameters": ["2*m", "m"],
                },
            ],
            "physical_threshold_claim": False,
            "verdict": "the finite representation reduct contains no mass coordinate and therefore emits no threshold location or ordering",
        },
    }


def matching_objects() -> dict[str, dict[str, str]]:
    reasons = {
        "ordered_eft_intervals": "active light/heavy census and source mass ordering are absent",
        "threshold_locations": "finite parents emit no mass coordinate or physical threshold",
        "decoupling_maps": "no source threshold crossing or heavy-field action is emitted",
        "scheme_maps": "no operational renormalization-scheme selector is emitted",
        "jacobians": "the finite scheme and threshold maps required for differentiation are absent",
        "finite_order_term_masks": "the complete action and loop-order source law are absent",
        "certified_vector_remainders": "no source domain, norm, higher-order kernel, or enclosure is emitted",
    }
    return {
        key: {"status": "not_emitted", "reason": reason}
        for key, reason in reasons.items()
    }


def invisible_sector_boundary() -> dict[str, Any]:
    return {
        "equivalence_scope": (
            "strict one-loop gauge coefficients modulo direct summands with "
            "zero SU(3), SU(2), and U(1) representation indices"
        ),
        "exact_result": {
            "delta_b_Y": "0",
            "delta_b_2": "0",
            "delta_b_3": "0",
            "status": "proved_after_imported_one_loop_gauge_functional",
        },
        "full_WZ_decoupling_proved": False,
        "missing_for_full_WZ_decoupling": [
            "source proof of zero Yukawa and scalar vertices",
            "source proof of zero mass mixing with the retained sector",
            "all-order or declared-finite-order decoupling map",
        ],
        "census_rule": (
            "enumerate every W/Z-coupled field and operator; quotient only by "
            "summands whose complete interaction and mixing vertices are "
            "proved zero"
        ),
        "status": "GAUGE_INDEX_INVARIANCE_PROVED__FULL_WZ_DECOUPLING_OPEN",
    }


def build_frontier() -> dict[str, Any]:
    policy = json.loads(POLICY_PATH.read_text(encoding="utf-8"))
    require(policy.get("schema") == "oph.source_rg_policy.v1", "POLICY_SCHEMA", "wrong policy schema")
    source_rows, source_payloads = resolve_source_inputs(policy)
    target_firewall(source_rows, source_payloads, policy)

    indices = representation_indices(
        source_payloads["finite_matter_representation"],
        source_payloads["conditional_rank_three_screen_context"],
        source_payloads["global_charge_lattice"],
    )
    witnesses = multiplicity_witnesses(
        source_payloads["mandatory_multiplicity_countermodels"]
    )
    policy_pin = file_pin(POLICY_PATH)
    input_digest = canonical_sha256(
        {
            "policy": policy_pin,
            "artifacts": source_rows,
        }
    )
    conditional = {
        "N_g": 3,
        "N_H": 1,
        "coefficients": beta_coefficients(3, 1),
        "gut_normalized_b_1": "41/10",
        "status": "conditional_declared_completion_not_OPH_selected",
        "promotion_allowed": False,
        "resolved_bounded_context": [
            "#634 finite local action domain; continuum promotion and physical coefficients absent"
        ],
        "required_open_attachments": [
            "#569 physical family attachment",
            "#636 physical scalar action and kinetic normalization",
            "#637 complete source Yukawa matrices",
            "#631 local physical screen/electroweak carrier",
            "#632 complete W/Z-coupled census modulo proved zero-vertex decoupling",
        ],
        "explicit_non_dependencies": [
            "#630 scalar/Yukawa/FJ integration",
            "#638 source-to-FJ coordinate map",
        ],
    }
    lean_raw = LEAN_PATH.read_bytes()
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": 32,
        "status": "PARTIAL_EXACT_REPRESENTATION_INDICES__SOURCE_MATCHING_OPEN",
        "promotion_allowed": False,
        "source_inputs": {
            "policy": policy_pin,
            "artifacts": source_rows,
            "input_bundle_digest": input_digest,
            "target_paths_present": False,
        },
        "qft_import_boundary": {
            **policy["imported_qft_premises"][0],
            "oph_native_one_loop_beta_theorem": False,
            "external_593_packet_consumed": False,
        },
        "formal_certificate": {
            "path": LEAN_PATH.relative_to(REPO_ROOT).as_posix(),
            "bytes": len(lean_raw),
            "byte_sha256": hashlib.sha256(lean_raw).hexdigest(),
            "theorems": [
                "representation_indices",
                "declared_three_one_evaluation",
                "family_copy_shift",
                "scalar_copy_shift",
                "zero_gauge_index_direct_sum_invariance",
            ],
            "scope": "axiom-free rational arithmetic after the imported one-loop QFT functional; no physical multiplicity or matching selection",
        },
        "representation_indices": indices,
        "parametric_one_loop_law": parametric_law(),
        "invisible_sector_boundary": invisible_sector_boundary(),
        "conditional_evaluations": [conditional],
        "nonidentifiability_witnesses": witnesses,
        "matching_objects": matching_objects(),
        "acceptance_map": [
            {
                "criterion": "target_clean_source_emits_complete_matching_packet",
                "status": "open",
                "landed": [
                    "target-clean finite representation provenance",
                    "exact parametric one-loop gauge coefficient law under an imported QFT functional",
                    "independent arithmetic replay and non-selection witnesses",
                ],
                "missing": list(matching_objects()),
            },
            {
                "criterion": "dependencies_imports_and_blockers_explicit",
                "status": "partial",
                "landed": [
                    "resolved finite-parent hashes",
                    "QFT import boundary",
                    "family, scalar, threshold, and scheme non-selection boundaries",
                ],
                "missing": [
                    "physical family attachment",
                    "physical scalar action and kinetic normalization",
                    "complete source Yukawa matrices",
                    "physical local carrier attachment",
                    "complete W/Z-coupled census modulo proved zero-vertex decoupling",
                ],
            },
            {
                "criterion": "companion_surfaces_aligned",
                "status": "partial",
                "landed": [
                    "machine contract, P trunk, particle ledgers, completion plan, and issue comments",
                ],
                "missing": [
                    "paper prose changes wait for a positive matching theorem",
                ],
            },
        ],
    }
    payload["subject_digest"] = canonical_sha256(payload)
    return payload


def serialize(payload: dict[str, Any]) -> bytes:
    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check-byte-exact", action="store_true")
    args = parser.parse_args()

    payload = build_frontier()
    built = serialize(payload)
    if args.check_byte_exact:
        require(args.output.is_file(), "OUTPUT_MISSING", f"committed output missing: {args.output}")
        require(
            args.output.read_bytes() == built,
            "OUTPUT_DRIFT",
            "committed frontier is not byte-identical to a clean rebuild",
        )
        print(f"byte-exact: {args.output}")
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(built)
    print(f"saved: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
