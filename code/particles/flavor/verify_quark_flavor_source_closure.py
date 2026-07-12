#!/usr/bin/env python3
"""Verify exact S3/D12 lemmas and emit the conditional flavor-source contract.

The exact output covers centered-plane algebra and normalized-trace arithmetic.
The source-carrier, channel assignment, physical scale, and RG packet remain
explicit hypotheses.  This verifier contains no quark comparison values.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp


CODE_ROOT = Path(__file__).resolve().parents[2]
FLAVOR_RUNS = CODE_ROOT / "particles" / "runs" / "flavor"
DEFAULT_OUTPUT = FLAVOR_RUNS / "quark_flavor_source_closure_contract.json"
EXPECTED_PROOF_BUNDLE_SHA256 = (
    "86e476d4b9c7910aefa0d259395e63e40b0b71c1534cf192b4f578cb21320759"
)

PIXEL_CERTIFICATE = (
    CODE_ROOT
    / "particles"
    / "hierarchy"
    / "certificates"
    / "R_P_source_audit_pixel_certificate.json"
)
FAMILY_KERNEL = FLAVOR_RUNS / "family_transport_kernel.json"
ODD_RESPONSE = FLAVOR_RUNS / "quark_odd_response_law.json"
MEAN_LAW = FLAVOR_RUNS / "quark_sector_mean_split.json"
SCHEME_OBSTRUCTION = FLAVOR_RUNS / "quark_running_mass_scheme_convention_obstruction.json"
AXIOM_NO_GO = FLAVOR_RUNS / "quark_axiom_level_yukawa_moduli_nonidentifiability.json"
RSCC_CANDIDATE = FLAVOR_RUNS / "quark_rscc_completion_candidate.json"
RSCC_AUDIT = FLAVOR_RUNS / "quark_rscc_completion_candidate_audit.json"
RSCC_ARITHMETIC = FLAVOR_RUNS / "quark_rscc_module_arithmetic.json"


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _exact_algebra() -> dict[str, Any]:
    x, rho, amplitude = sp.symbols("x rho amplitude", real=True)

    def center(vector: list[sp.Expr | int]) -> sp.Matrix:
        column = sp.Matrix(vector)
        return column - sp.ones(3, 1) * sum(column) / 3

    linear = sp.simplify(center([-1, x, 1]))
    quadratic = sp.simplify(center([1, x**2, 1]))
    gram = sp.Matrix(
        [
            [linear.dot(linear), linear.dot(quadratic)],
            [quadratic.dot(linear), quadratic.dot(quadratic)],
        ]
    )
    determinant = sp.factor(gram.det())
    expected_determinant = sp.Rational(4, 3) * (1 - x**2) ** 2

    up_ray = sp.Matrix(
        [
            -(2 * rho + 1) / (3 * (1 + rho)),
            (rho - 1) / (3 * (1 + rho)),
            (rho + 2) / (3 * (1 + rho)),
        ]
    )
    down_ray = sp.Matrix(
        [
            -(rho + 2) / (3 * (1 + rho)),
            (1 - rho) / (3 * (1 + rho)),
            (2 * rho + 1) / (3 * (1 + rho)),
        ]
    )
    alpha, beta = sp.symbols("alpha beta")
    up_solution = sp.solve(
        list(2 * amplitude * up_ray - alpha * linear - beta * quadratic),
        [alpha, beta],
        dict=True,
    )[0]
    down_solution = sp.solve(
        list(2 * amplitude * down_ray - alpha * linear - beta * quadratic),
        [alpha, beta],
        dict=True,
    )[0]
    expected_up_beta = sp.factor(
        amplitude
        * (-rho * x + rho - x - 1)
        / ((1 + rho) * (x**2 - 1))
    )
    expected_down_beta = sp.factor(
        amplitude
        * (-rho * x - rho - x + 1)
        / ((1 + rho) * (x**2 - 1))
    )

    r = sp.symbols("r", positive=True, real=True)
    rho_from_r = 3 / (2 + r)
    x_from_r = (r - 1) / (r + 1)
    rho_x_residual = sp.factor(rho_from_r - 3 * (1 - x_from_r) / (3 - x_from_r))

    checks = {
        "LQ_gram_determinant": sp.simplify(determinant - expected_determinant) == 0,
        "up_linear_coordinate": sp.simplify(up_solution[alpha] - amplitude) == 0,
        "down_linear_coordinate": sp.simplify(down_solution[alpha] - amplitude) == 0,
        "up_ray_quadratic_coordinate": sp.simplify(up_solution[beta] - expected_up_beta) == 0,
        "down_ray_quadratic_coordinate": sp.simplify(
            down_solution[beta] - expected_down_beta
        )
        == 0,
        "rho_x_relation": rho_x_residual == 0,
    }
    if not all(checks.values()):
        raise AssertionError(f"exact flavor algebra failed: {checks}")

    dimensions = {
        "family_module": 5,
        "family_matrix_slots": 25,
        "up_odd_module": 2,
        "down_odd_module": 5,
        "up_even_oriented_module": 10,
        "down_even_electroweak_module": 4,
    }
    normalized_trace_coefficients = {
        "isotropic_heat_slot": sp.Rational(1, dimensions["family_matrix_slots"]),
        "up_odd_rank_one": sp.Rational(1, dimensions["up_odd_module"]),
        "down_odd_rank_one": sp.Rational(1, dimensions["down_odd_module"]),
        "up_even_rank_one": sp.Rational(1, dimensions["up_even_oriented_module"]),
        "down_even_rank_one": sp.Rational(1, dimensions["down_even_electroweak_module"]),
    }
    expected_coefficients = {
        "isotropic_heat_slot": sp.Rational(1, 25),
        "up_odd_rank_one": sp.Rational(1, 2),
        "down_odd_rank_one": sp.Rational(1, 5),
        "up_even_rank_one": sp.Rational(1, 10),
        "down_even_rank_one": sp.Rational(1, 4),
    }
    if normalized_trace_coefficients != expected_coefficients:
        raise AssertionError("normalized-trace denominator check failed")

    return {
        "proof_status": "exact_symbolic_verification_passed",
        "checks": checks,
        "centered_basis": {
            "L": [str(sp.factor(entry)) for entry in linear],
            "Q": [str(sp.factor(entry)) for entry in quadratic],
            "gram_determinant": str(determinant),
            "basis_domain": "x != +/-1",
        },
        "ray_decomposition": {
            "up_linear_coordinate": str(sp.factor(up_solution[alpha])),
            "up_quadratic_coordinate": str(sp.factor(up_solution[beta])),
            "down_linear_coordinate": str(sp.factor(down_solution[alpha])),
            "down_quadratic_coordinate": str(sp.factor(down_solution[beta])),
        },
        "rho_x_relation": "rho=3*(1-x)/(3-x)",
        "normalized_trace_theorem": {
            "statement": (
                "A complex-linear U(d)-conjugation-invariant functional with ell(I)=1 "
                "is Tr/d. A U(d)xU(d)-isotropic Hilbert-Schmidt width with total "
                "width t is (t/d^2)I."
            ),
            "dimensions": dimensions,
            "rank_one_or_slot_coefficients": {
                key: str(value) for key, value in normalized_trace_coefficients.items()
            },
            "conditional_denominator_tuple": [5, 2, 5, 10, 4],
            "scope_boundary": (
                "The coefficients follow only after the physical heat, odd, and even "
                "channel modules, orientations, and signs are supplied. Representation "
                "theory does not select those modules."
            ),
        },
    }


def build_artifact(proof_bundle: Path | None = None) -> dict[str, Any]:
    pixel = _read_json(PIXEL_CERTIFICATE)
    kernel = _read_json(FAMILY_KERNEL)
    odd = _read_json(ODD_RESPONSE)
    mean = _read_json(MEAN_LAW)
    scheme = _read_json(SCHEME_OBSTRUCTION)
    no_go = _read_json(AXIOM_NO_GO)
    rscc_candidate = _read_json(RSCC_CANDIDATE)
    rscc_audit = _read_json(RSCC_AUDIT)
    rscc_arithmetic = _read_json(RSCC_ARITHMETIC)

    receipts = {
        "F1_pixel_gauge_source_root": {
            "closed": False,
            "current_status": pixel.get("status"),
            "rscc_candidate_status": rscc_candidate["provenance"][
                "pixel_branch_status"
            ],
            "rscc_stage5_quark_ancestor": rscc_candidate["provenance"][
                "pixel_source_uses_internal_stage5_quark_model"
            ],
            "required": "unique no-target interval or Krawczyk source-root certificate",
        },
        "F2_refinement_natural_flavor_carrier": {
            "closed": False,
            "current_kernel_status": kernel.get("status"),
            "current_kernel_proof_status": kernel.get("proof_status"),
            "odd_response_status": odd.get("proof_status"),
            "odd_response_missing_object": odd.get("upstream_missing_object"),
            "required": (
                "source-emitted labeled simple-spectrum family generator, charged seed, "
                "and common-refinement intertwiners"
            ),
            "rscc_contribution": (
                "proposes F=1+2*V_std but supplies no physical regular-heat-to-F "
                "attachment, labeled generator, or refinement transport"
            ),
        },
        "F3_physical_channel_identification_functor": {
            "closed": False,
            "required": (
                "identify the physical S3 heat carrier, tau_f, rho/x dictionary, and edge "
                "response equations, then derive the 5,2,5,10,4 response modules, signs, "
                "orientation, and exclusion of competing source-natural assignments"
            ),
            "rscc_arithmetic_status": rscc_arithmetic["proof_status"],
            "rscc_physical_functor_closed": False,
        },
        "F4_affine_sector_mean_law_on_source_carrier": {
            "closed": False,
            "current_algebra_status": mean.get("proof_status"),
            "domain_requirement": (
                "1-x^2-x^2/(1+rho) != 0, in addition to x != +/-1 and rho > 0"
            ),
            "reason_open": (
                "the current-family algebra is closed only after candidate/target-derived "
                "spread inputs; it is not derived on a source-emitted carrier"
            ),
            "rscc_status": (
                "explicit formulas proposed, but inherited rays, even responses, and "
                "A/B affine law remain postulated"
            ),
        },
        "F5_dimensionful_or_dimensionless_physical_readout": {
            "closed": False,
            "required": (
                "a declared source-scale mass normalization or common-scale dimensionless "
                "Yukawa normalization including m_q(mu)=y_q(mu)*v(mu)/sqrt(2) and the "
                "running Higgs expectation value in the same scheme"
            ),
            "rscc_status": (
                "g_ch/v and m/E_star are dimensionally improved candidate ratios; the "
                "optional D10 GeV display mixes nonmatching P/alpha and scale branches"
            ),
        },
        "F6_source_only_RG_threshold_scheme_packet": {
            "closed": False,
            "single_running_sextet_claim_allowed": scheme["row_partition"][
                "single_running_quark_sextet_claim_allowed"
            ],
            "required": (
                "locally-Lipschitz beta system with existence/non-blowup over every required "
                "interval, frozen thresholds and matching maps, top conversion, common "
                "scheme, and no-target provenance"
            ),
            "rscc_packet_supplied": False,
        },
    }
    all_receipts_closed = all(receipt["closed"] for receipt in receipts.values())

    dependency_paths = [
        PIXEL_CERTIFICATE,
        FAMILY_KERNEL,
        ODD_RESPONSE,
        MEAN_LAW,
        SCHEME_OBSTRUCTION,
        AXIOM_NO_GO,
        RSCC_CANDIDATE,
        RSCC_AUDIT,
        RSCC_ARITHMETIC,
    ]
    bundle_receipt: dict[str, Any] = {
        "expected_sha256": EXPECTED_PROOF_BUNDLE_SHA256,
        "path_checked": None,
        "matches": None,
    }
    if proof_bundle is not None:
        observed_sha256 = _sha256(proof_bundle)
        bundle_receipt.update(
            {
                "path_checked": str(proof_bundle),
                "observed_sha256": observed_sha256,
                "matches": observed_sha256 == EXPECTED_PROOF_BUNDLE_SHA256,
            }
        )
    return {
        "artifact": "oph_quark_flavor_source_closure_contract_v1",
        "claim_class": "exact_algebra_plus_conditional_physical_completion_contract",
        "promotion_allowed": False,
        "submitted_proof_bundle_receipt": bundle_receipt,
        "rscc_formula_level_candidate": {
            "artifact": rscc_candidate["artifact"],
            "audit_artifact": rscc_audit["artifact"],
            "claim_class": rscc_candidate["claim_class"],
            "promotion_allowed": False,
            "old_candidate_flavor_decimals_consumed_directly": rscc_candidate[
                "provenance"
            ]["old_candidate_flavor_decimals_consumed"],
            "exact_arithmetic_status": rscc_arithmetic["proof_status"],
            "maximum_mixed_chart_relative_residual_percent": rscc_audit[
                "descriptive_mixed_chart_comparison"
            ]["max_abs_relative_error_percent"],
            "negative_control_beats_full_rscc": rscc_audit["negative_control"][
                "beats_full_rscc_maximum_error"
            ],
            "all_F1_to_F6_receipts_remain_open": True,
            "status_statement": (
                "RSCC turns parts of the unspecified source packet into an explicit, "
                "falsifiable formula-level proposal. Its module/effect/sign ledger, "
                "Gaussian truncation, inherited affine readout, source root, scale, "
                "and RG packet are not derived, so it discharges no physical receipt."
            ),
        },
        "exact_algebra": _exact_algebra(),
        "necessity_theorem": {
            "artifact": no_go.get("artifact"),
            "proof_status": no_go.get("proof_status"),
            "counterfamily_parameter_space": no_go["counterfamily"]["parameter_space"],
            "conclusion": (
                "Any successful flavor source functional must be nonconstant on the "
                "independent positive up/down centered-spread rescaling orbit."
            ),
        },
        "conditional_refinement_descent": {
            "theorem_status": "exact_given_hypotheses",
            "hypotheses": [
                "a source-labeled simple-spectrum family generator K_r",
                "a source charged seed Y_ch,r in the same labeled basis",
                "exact refinement intertwiners K_s U_sr = U_sr K_r",
                "monomial source-label transport with U_sr e_i,r equal to a phase times e_P(i),s",
                "a uniform positive spectral-gap lower bound",
                "positive g_ch, common dimensional typing of g_ch and Y_ch entries, and a uniform nonzero lower bound on every selected charged-seed edge",
            ],
            "conclusion": (
                "ordered spectral gaps, trace/gap functionals, normalized edge magnitudes, "
                "and phase-invariant projector witnesses descend and are perturbatively stable"
            ),
            "current_hypotheses_discharged": False,
        },
        "conditional_RG_uniqueness": {
            "theorem_status": "standard_ODE_theorem_given_hypotheses",
            "statement": (
                "A frozen locally-Lipschitz beta system that exists without blowup across the "
                "declared intervals, with single-valued threshold matching, has a unique "
                "transported image. This proves uniqueness after an RG packet is supplied, "
                "not derivation of that packet."
            ),
            "current_packet_supplied": False,
        },
        "flavor_source_closure_receipts": receipts,
        "all_six_receipts_closed": all_receipts_closed,
        "conditional_mass_evaluator_unique_if_all_receipts_close": True,
        "current_repository_emits_physical_quark_sextet": False,
        "current_numeric_formula_role": (
            "post-hoc target-informed repository-template and RSCC "
            "formula-discovery diagnostics"
        ),
        "sufficient_new_content_not_proved_minimal": (
            "a unique, target-free Flavor Source Closure functor from observer-like "
            "self-reading OPH overlap/record receipts to the carrier, response modules, "
            "physical scale, and RG packet, nonconstant on the rescaling orbit"
        ),
        "dependency_files": {
            str(path.relative_to(CODE_ROOT)): {"sha256": _sha256(path)}
            for path in dependency_paths
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--proof-bundle", type=Path)
    parser.add_argument("--print-json", action="store_true")
    args = parser.parse_args()

    artifact = build_artifact(args.proof_bundle)
    text = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    if args.print_json:
        print(text, end="")
    else:
        print(f"saved: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
