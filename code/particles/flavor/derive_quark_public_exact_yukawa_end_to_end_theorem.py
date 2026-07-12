#!/usr/bin/env python3
"""Emit the public exact end-to-end Yukawa theorem above the descended sigma datum."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PUBLIC_SIGMA_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_public_physical_sigma_datum_descent.json"
EXACT_PDG_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_exact_pdg_end_to_end_theorem.json"
EXACT_YUKAWA_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_exact_yukawa_end_to_end_theorem.json"
SCHEME_OBSTRUCTION_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_running_mass_scheme_convention_obstruction.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_end_to_end_theorem.json"

MISSING_FOR_PROMOTION = [
    "QUARK_SIGMA_SOURCE_QUOTIENT",
    "QUARK_SIGMA_SOURCE_SELECTOR",
    "QUARK_EDGE_STATISTICS_CORRECTION_THEOREM",
    "QUARK_SIGMA_REFINEMENT_COMPATIBILITY",
    "NO_TARGET_LEAK_DAG_QUARK_SIGMA_SOURCE",
]

MASS_YUKAWA_CONSISTENCY_BLOCKER = "QUARK_EXACT_MASS_YUKAWA_SURFACE_CONSISTENCY"
PHYSICAL_YUKAWA_NORMALIZATION_BLOCKER = "QUARK_COMMON_SCALE_DIMENSIONLESS_YUKAWA_CERTIFICATE"
MASS_YUKAWA_REL_TOL = 1.0e-6
MASS_YUKAWA_ABS_TOL_GEV = 1.0e-8
MASS_SINGULAR_VALUE_INDEX = {
    "u": ("singular_values_u", 0),
    "c": ("singular_values_u", 1),
    "t": ("singular_values_u", 2),
    "d": ("singular_values_d", 0),
    "s": ("singular_values_d", 1),
    "b": ("singular_values_d", 2),
}


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _mass_yukawa_consistency(
    exact_masses: dict[str, Any],
    forward_yukawas: dict[str, Any],
) -> dict[str, Any]:
    """Check that the declared sextet is the singular-value spectrum of Y_u/Y_d."""

    comparisons: dict[str, dict[str, Any]] = {}
    malformed_inputs: list[str] = []
    all_consistent = True
    for particle, (singular_key, index) in MASS_SINGULAR_VALUE_INDEX.items():
        try:
            declared_mass = float(exact_masses[particle])
        except (KeyError, TypeError, ValueError):
            malformed_inputs.append(f"exact_running_values_gev.{particle}")
            all_consistent = False
            continue
        singular_values = forward_yukawas.get(singular_key)
        if not isinstance(singular_values, list) or len(singular_values) != 3:
            malformed_inputs.append(singular_key)
            all_consistent = False
            continue
        try:
            singular_value = float(singular_values[index])
        except (TypeError, ValueError):
            malformed_inputs.append(f"{singular_key}[{index}]")
            all_consistent = False
            continue

        absolute_residual = abs(declared_mass - singular_value)
        scale = max(abs(declared_mass), abs(singular_value), MASS_YUKAWA_ABS_TOL_GEV)
        relative_residual = absolute_residual / scale
        consistent = math.isclose(
            declared_mass,
            singular_value,
            rel_tol=MASS_YUKAWA_REL_TOL,
            abs_tol=MASS_YUKAWA_ABS_TOL_GEV,
        )
        all_consistent = all_consistent and consistent
        comparisons[particle] = {
            "declared_mass_gev": declared_mass,
            "forward_yukawa_singular_value": singular_value,
            "absolute_residual_gev": absolute_residual,
            "relative_residual": relative_residual,
            "consistent_within_tolerance": consistent,
        }

    if malformed_inputs:
        status = "invalid_mass_or_yukawa_spectrum_surface"
    elif all_consistent:
        status = "consistent_single_exact_sextet_matrix_pair"
    else:
        status = "mismatched_mass_sextet_and_forward_yukawa_singular_values"
    return {
        "status": status,
        "consistent": all_consistent and not malformed_inputs,
        "single_exact_sextet_matrix_pair_claim_allowed": all_consistent and not malformed_inputs,
        "relative_tolerance": MASS_YUKAWA_REL_TOL,
        "absolute_tolerance_gev": MASS_YUKAWA_ABS_TOL_GEV,
        "comparisons": comparisons,
        "malformed_inputs": malformed_inputs,
    }


def build_artifact(
    public_sigma_theorem: dict[str, Any],
    exact_pdg_theorem: dict[str, Any],
    exact_yukawa_theorem: dict[str, Any],
    scheme_obstruction: dict[str, Any] | None = None,
) -> dict[str, Any]:
    forward = dict(exact_yukawa_theorem["forward_yukawa_artifact"])
    exact_masses = dict(exact_pdg_theorem["exact_running_values_gev"])
    mass_yukawa_consistency = _mass_yukawa_consistency(exact_masses, forward)
    exact_pair_consistent = mass_yukawa_consistency["consistent"] is True
    scheme_obstruction = scheme_obstruction or {}
    dimensional_audit = dict(scheme_obstruction.get("stored_matrix_dimensional_audit") or {})
    physical_yukawa_normalization_closed = (
        dimensional_audit.get("certified_physical_yukawa_matrices") is True
        and dimensional_audit.get("exact_forward_yukawa_claim_allowed") is True
    )
    sigma_non_circularity = dict(public_sigma_theorem.get("non_circularity_status") or {})
    promotion_allowed = (
        exact_pair_consistent
        and physical_yukawa_normalization_closed
        and public_sigma_theorem.get("public_promotion_allowed") is True
        and sigma_non_circularity.get("promotion_allowed", public_sigma_theorem.get("public_promotion_allowed")) is True
        and public_sigma_theorem.get("proof_status") == "closed_source_only_public_physical_sigma_datum_descent"
    )
    if not exact_pair_consistent:
        proof_status = "blocked_inconsistent_exact_mass_yukawa_surfaces"
    elif not physical_yukawa_normalization_closed:
        proof_status = "blocked_mixed_scheme_dimensionful_mass_textures_not_physical_yukawas"
    elif promotion_allowed:
        proof_status = "closed_source_only_public_exact_yukawa_end_to_end_theorem"
    else:
        proof_status = "blocked_target_derived_sigma_source_missing"
    return {
        "artifact": "oph_quark_public_exact_yukawa_end_to_end_theorem",
        "generated_utc": _timestamp(),
        "proof_status": proof_status,
        "target_name": (
            "source_only_public_exact_forward_quark_yukawas"
            if promotion_allowed
            else (
                "inconsistent_selected_class_mass_yukawa_support_surfaces"
                if not exact_pair_consistent
                else (
                    "selected_class_mixed_scheme_mass_texture_audit"
                    if not physical_yukawa_normalization_closed
                    else "selected_class_conditional_forward_quark_yukawas"
                )
            )
        ),
        "theorem_scope": public_sigma_theorem["theorem_scope"],
        "claim_tier": (
            "source_only_public_yukawa_theorem"
            if promotion_allowed
            else (
                "inconsistent_exact_support_surfaces"
                if not exact_pair_consistent
                else (
                    "mixed_scheme_dimensionful_mass_texture_audit"
                    if not physical_yukawa_normalization_closed
                    else "selected_class_conditional_on_source_sigma"
                )
            )
        ),
        "public_promotion_allowed": promotion_allowed,
        "source_only_sigma_emitted": bool(public_sigma_theorem.get("source_only_sigma_emitted")),
        "downstream_algebra_closed": exact_pair_consistent and physical_yukawa_normalization_closed,
        "mass_texture_algebra_closed": exact_pair_consistent,
        "physical_yukawa_construction_closed": physical_yukawa_normalization_closed,
        "display_allowed_as_selected_class_exact_witness": (
            exact_pair_consistent and physical_yukawa_normalization_closed
        ),
        "display_allowed_as_selected_class_mass_texture_audit": exact_pair_consistent,
        "mass_yukawa_consistency": mass_yukawa_consistency,
        "scheme_and_dimensional_normalization": {
            "artifact": scheme_obstruction.get("artifact"),
            "proof_status": scheme_obstruction.get("proof_status"),
            "physical_yukawa_normalization_closed": physical_yukawa_normalization_closed,
            "matrix_classification": dimensional_audit.get("current_classification"),
            "common_renormalization_scale": dimensional_audit.get("common_renormalization_scale"),
            "dimensionless_normalization_supplied": dimensional_audit.get(
                "dimensionless_normalization_supplied"
            ),
            "promotion_requirements": dimensional_audit.get("promotion_requirements", []),
        },
        "non_circularity_status": {
            "promotion_allowed": promotion_allowed,
            "public_sigma_promotion_allowed": public_sigma_theorem.get("public_promotion_allowed"),
            "public_sigma_proof_status": public_sigma_theorem.get("proof_status"),
            "target_derived_sigma_datum_used": sigma_non_circularity.get("target_derived_sigma_datum_used"),
            "source_sigma_selector_closed": sigma_non_circularity.get("source_sigma_selector_closed"),
            "exact_mass_yukawa_pair_consistent": exact_pair_consistent,
            "physical_yukawa_normalization_closed": physical_yukawa_normalization_closed,
            "missing_source_object": None
            if promotion_allowed
            else "quark_sigma_source_datum_no_target_leak_required",
            "strict_audit_label": "source_only_public_yukawa_theorem"
            if promotion_allowed
            else (
                "inconsistent_exact_support_surfaces"
                if not exact_pair_consistent
                else (
                    "mixed_scheme_dimensionful_mass_texture_audit"
                    if not physical_yukawa_normalization_closed
                    else "selected_class_conditional_on_source_sigma"
                )
            ),
        },
        "supporting_theorem_artifacts": {
            "public_sigma_datum_descent": public_sigma_theorem["artifact"],
            "local_exact_pdg_wrapper": exact_pdg_theorem["artifact"],
            "local_exact_yukawa_wrapper": exact_yukawa_theorem["artifact"],
        },
        "theorem_statement": (
            "The stored declared mass sextet and forward Yukawa matrices do not form one exact output pair: the "
            "singular values of Y_u and Y_d disagree with one or more declared running masses. They remain separate "
            "audit surfaces until rebuilt from one common input surface and certified by the mass-Yukawa consistency "
            "check. Independently, public source-only promotion still requires a no-target source sigma theorem."
            if not exact_pair_consistent
            else (
                (
                    "The stored matrices reproduce the declared mixed-convention mass coordinates as dimensionful "
                    "mass textures. They are not physical Yukawa matrices: the artifact supplies neither one common "
                    "renormalization scale nor threshold transport, top conversion, a running Higgs expectation value, "
                    "or the dimensionless normalization y_q=sqrt(2)m_q/v. Public physical-Yukawa promotion is blocked "
                    "independently of the source-spread obstruction."
                )
                if not physical_yukawa_normalization_closed
                else (
                    "Conditional on a no-target source theorem emitting the physical sigma datum on the public quark "
                    "frame class selected by P, the selected bridge-fiber descent makes that datum representative-independent, "
                    "the affine mean law emits the absolute sector scales algebraically, the ordered three-point readout "
                    "yields the scheme-labelled quark coordinates, and a common-scale dimensionless conversion emits "
                    "physical Yukawa matrices Y_u and Y_d on that selected class."
                )
            )
        ),
        "selected_public_physical_frame_class": public_sigma_theorem["selected_public_physical_frame_class"],
        "descended_physical_sigma_datum": public_sigma_theorem["descended_physical_sigma_datum"],
        "public_exact_outputs": {
            "pair_status": mass_yukawa_consistency["status"],
            "single_exact_mass_texture_pair_claim_allowed": exact_pair_consistent,
            "single_exact_physical_yukawa_pair_claim_allowed": (
                exact_pair_consistent and physical_yukawa_normalization_closed
            ),
            "exact_running_values_gev": exact_masses,
            "forward_yukawa_artifact": forward,
        },
        "lemma_chain": [
            public_sigma_theorem["artifact"],
            "oph_quark_absolute_readout_algebraic_collapse",
            exact_pdg_theorem["artifact"],
            exact_yukawa_theorem["artifact"],
        ],
        "minimal_exact_blocker_set": (
            []
            if promotion_allowed
            else (
                [
                    *MISSING_FOR_PROMOTION,
                    MASS_YUKAWA_CONSISTENCY_BLOCKER,
                    PHYSICAL_YUKAWA_NORMALIZATION_BLOCKER,
                ]
                if not exact_pair_consistent
                else [*MISSING_FOR_PROMOTION, PHYSICAL_YUKAWA_NORMALIZATION_BLOCKER]
                if not physical_yukawa_normalization_closed
                else MISSING_FOR_PROMOTION
            )
        ),
        "notes": [
            (
                "This is the theorem wrapper that upgrades the closed declared-carrier exact Yukawa chain to the selected public class."
                if promotion_allowed
                else (
                    "The mass and matrix artifacts are retained for separate audit only; they may not be described as one exact sextet/Yukawa witness while their spectra disagree."
                    if not exact_pair_consistent
                    else (
                        "The numerical matrices are retained as mixed-scheme mass-texture audit witnesses; they are not certified physical Yukawa matrices."
                        if not physical_yukawa_normalization_closed
                        else "This artifact displays the selected-class exact witness but is not promotable under the strict non-circularity audit because the public sigma datum descends from an exact target surface."
                    )
                )
            ),
            "The matrices Y_u and Y_d are the numerical matrices emitted on the local chain; the consistency block records whether their singular values match the separately declared sextet.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the public exact Yukawa end-to-end theorem artifact.")
    parser.add_argument("--public-sigma-theorem", default=str(PUBLIC_SIGMA_JSON))
    parser.add_argument("--exact-pdg-theorem", default=str(EXACT_PDG_JSON))
    parser.add_argument("--exact-yukawa-theorem", default=str(EXACT_YUKAWA_JSON))
    parser.add_argument("--scheme-obstruction", default=str(SCHEME_OBSTRUCTION_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(
        _load_json(Path(args.public_sigma_theorem)),
        _load_json(Path(args.exact_pdg_theorem)),
        _load_json(Path(args.exact_yukawa_theorem)),
        _load_json(Path(args.scheme_obstruction)),
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
