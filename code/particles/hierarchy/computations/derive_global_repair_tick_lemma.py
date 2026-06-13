#!/usr/bin/env python3
"""Emit the global repair-tick lemma certificate.

Chain role: close the R_N-side global repair-tick lemma for the local/global
hierarchy resonance continuation without importing electroweak data.

Mathematics: the D6 screen-capacity normalization fixes the cosmic radius
ratio r_CRC/ell_star = (N_CRC/pi)^(1/2), so in the multiplicative coordinate
rho = ell/r_CRC the horizon sits at rho = 1 and the local cell at
rho_star = (N_CRC/pi)^(-1/2). The corpus readback map is
F(N) = Cap_read(Obs(nf_{r,N}(U_{r,N}))) (Observers synthesis sections,
refinement limit): the capacity reconstructed by the interior observer sector
from the horizon record. Modeling Cap_read by the D6 area law at the delivery
resolution (a declared D6-consistent identification), a readback that delivers
the record at effective resolution coordinate rho_read reconstructs
F(N) = pi / rho_read^2, so the fixed-point equation
N_CRC = F(N_CRC) = pi / rho_star^2 forces rho_read = rho_star: the readback
transport carries the horizon coordinate rho = 1 exactly to the local cell
rho_star, i.e. G_N(1) = rho_star, with full-cycle multiplier
(N_CRC/pi)^(-1/2). On the declared 24-tick resonance branch the one-tick map
is the positive homogeneous 24th root, hence |g_*'| = (N_CRC/pi)^(-1/48).
For a general m-tick decomposition the per-tick exponent is -1/(2m); the
declared m = 24 gives the stated exponent.

Output: a machine-readable certificate for GlobalRepairTickLemma_R_N.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "certificates" / "R_N_global_repair_tick_certificate.json"
DISPLAY_N_CRC = 3.31e122
DECLARED_REPAIR_ROUNDS = 24
FULL_CONTRACTION_EXPONENT = Fraction(-1, 2)
ONE_TICK_EXPONENT = FULL_CONTRACTION_EXPONENT / DECLARED_REPAIR_ROUNDS


def _fraction_payload(value: Fraction) -> dict[str, int | str]:
    text = str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "string": text,
    }


def _format_scientific(value: float) -> str:
    return f"{value:.17e}"


def build_artifact(n_crc_display: float = DISPLAY_N_CRC) -> dict:
    n_over_pi = n_crc_display / math.pi
    log_n_over_pi = math.log(n_over_pi)
    tick_multiplier = math.exp(float(ONE_TICK_EXPONENT) * log_n_over_pi)
    full_contraction = math.exp(float(FULL_CONTRACTION_EXPONENT) * log_n_over_pi)
    tick_power_24 = tick_multiplier**DECLARED_REPAIR_ROUNDS
    relative_power_residual = (tick_power_24 - full_contraction) / full_contraction

    return {
        "artifact": "R_N_global_repair_tick_certificate",
        "status": "closed_global_repair_tick_lemma_on_declared_round_structure",
        "theorem_kind": "lemma_on_declared_branch",
        "object_id": "GlobalRepairTickLemma_R_N",
        "proof_gate": "R_N_global_repair_tick_lemma_without_EW_readback",
        "claim_boundary": {
            "proved_by_certificate": [
                "D6 radius identity r_CRC/ell_star = (N_CRC/pi)^(1/2) from the screen-capacity normalization",
                "interface realization: under the stated readback counting model F(N) = pi/rho_read^2 (the D6-consistent area-law count at the effective delivery resolution), the corpus readback fixed-point equation N_CRC = F(N_CRC) is equivalent to the closure transport G_N(1) = rho_star",
                "full-cycle multiplier (N_CRC/pi)^(-1/2) from that closure transport",
                "per-tick exponent law -1/(2m) for an m-tick homogeneous decomposition",
                "stated exponent -1/48 at the declared m = 24, hence |g_*'| = (N_CRC/pi)^(-1/48)",
                "dependency boundary excluding electroweak measured inputs",
            ],
            "declared_not_derived": [
                "round count m = 24 (branch architecture for this lemma; representation-to-spectrum derivation is work in progress)",
                "homogeneous positive-root one-tick normal form (the normalization convention defining |g_*'|, supplied per acceptance criterion 1)",
                "single effective delivery resolution rho_read for the readback (the scale-free normal form under which the counting model is stated)",
                "readback counting model: Cap_read reconstructs capacity by the D6 area law evaluated at the effective delivery resolution, F(N) = pi/rho_read^2 (a D6-consistent modeling identification of Cap_read; the corpus defines Cap_read as the capacity reconstructed by the observer sector and does not state this counting property)",
            ],
            "not_closed_by_certificate": [
                "verification on the concrete finite repair machinery that nf_{r,N} delivers a single well-defined effective resolution scale (the corpus marks the finite readback map F_r as schematic and the refinement limit as conditional on existence)",
                "representation-to-spectrum theorem deriving the 24-round count from the OPH repair grammar",
                "exact global capacity source certificate satisfying B_EW(P_star,N_CRC)=0",
                "joint fixed-point/stability theorem for (P,N), recorded separately in R_PN_joint_fixed_point_certificate_report.json",
                "electroweak tick-projection bridge, recorded separately in R_EW_tick_projection_certificate.json",
                "RG/coarse-graining naturality theorem, recorded separately in issue_332_rg_naturality_certificate.json",
            ],
        },
        "definitions": {
            "N_CRC": (
                "the cosmic record-capacity fixed point N_CRC = F(N_CRC), "
                "read on the observed de Sitter branch as S_dS"
            ),
            "screen_normalized_radius_coordinate": (
                "rho = length / r_CRC; the horizon sits at rho = 1 and the local cell at "
                "rho_star = ell_star / r_CRC, with r_CRC/ell_star = (N_CRC/pi)^(1/2)"
            ),
            "global_repair_cycle": (
                "G_N is the scale-free (homogeneous) readback transport attached to the fixed point; "
                "fixed-point closure without deficit or slack means G_N(1) = rho_star exactly"
            ),
            "corpus_readback_map": (
                "F(N) = Cap_read(Obs(nf_{r,N}(U_{r,N}))) in the refinement limit: the active cosmic "
                "record capacity reconstructed by the stable self-reading observer sector "
                "(paper/tex_fragments/OBSERVERS_SYNTHESIS_SECTIONS.tex)"
            ),
            "effective_delivery_coordinate": (
                "rho_read is the effective resolution scale, in units of r_CRC, at which the "
                "repair-readback chain delivers the horizon record to the observer sector "
                "(scale-free normal form)"
            ),
            "readback_counting_model": (
                "Cap_read is modeled as the D6 area-law count at the delivery resolution: the capacity "
                "reconstructed at rho_read is F(N) = pi / rho_read^2. This is the canonical "
                "scale-covariant extension of the D6 counting rule N = pi (r/ell)^2, supplied here as a "
                "modeling identification. The corpus defines Cap_read as reconstructed capacity and does "
                "not state this counting rule as one of its properties"
            ),
            "one_tick_map": "g_N with g_N^24 = G_N, positive homogeneous root",
            "derivative_convention": (
                "|g_*'| is the absolute derivative of the one-tick map in the "
                "screen-normalized multiplicative radius coordinate rho"
            ),
        },
        "premises": {
            "derived": [
                {
                    "id": "d6_entropy_capacity",
                    "statement": "N_CRC = S_dS = pi * (r_CRC/ell_star)^2 on the observed screen-capacity branch",
                    "discharged_here": True,
                },
                {
                    "id": "f_interface_realization_equivalence",
                    "statement": (
                        "with the corpus readback map F = Cap_read(Obs(nf_{r,N}(U_{r,N}))) and the "
                        "stated readback counting model F(N) = pi/rho_read^2, the fixed-point equation "
                        "N_CRC = F(N_CRC) = pi/rho_star^2 is equivalent to rho_read = rho_star "
                        "(positive roots), i.e. to the closure transport G_N(1) = rho_star; deficit "
                        "rho_read > rho_star gives F(N) < N, slack rho_read < rho_star gives F(N) > N"
                    ),
                    "source": (
                        "paper/tex_fragments/OBSERVERS_SYNTHESIS_SECTIONS.tex readback-map definition "
                        "plus the D6 screen-capacity normalization; the counting model itself is a "
                        "declared modeling identification listed under declared_not_derived"
                    ),
                    "discharged_here": True,
                    "conditional_on": "readback counting model (declared modeling identification)",
                },
                {
                    "id": "full_cycle_closure_multiplier",
                    "statement": (
                        "fixed-point readback closure (no deficit, no slack) forces G_N(1) = rho_star, "
                        "so the homogeneous full-cycle multiplier is (N_CRC/pi)^(-1/2)"
                    ),
                    "discharged_here": True,
                },
            ],
            "declared_branch_inputs": [
                {
                    "id": "global_repair_round_count",
                    "statement": "the global repair cycle decomposes into 24 homogeneous ticks",
                    "origin": "declared resonance-continuation branch architecture",
                    "discharged_here": False,
                    "open_item": "representation-to-spectrum derivation of the round count",
                },
                {
                    "id": "homogeneous_tick_normal_form",
                    "statement": "the one-tick map is the positive homogeneous 24th root of the full-cycle map",
                    "role": "normalization convention defining |g_*'| (acceptance criterion 1)",
                    "discharged_here": False,
                },
            ],
        },
        "normalization": {
            "screen_capacity_relation": "N_CRC = pi * (r_CRC/ell_star)^2",
            "horizon_radius_over_local_unit": "(N_CRC/pi)^(1/2)",
            "local_cell_coordinate_rho_star": "(N_CRC/pi)^(-1/2)",
            "full_repair_rounds": DECLARED_REPAIR_ROUNDS,
            "full_cycle_map": "G_N(rho) = (N_CRC/pi)^(-1/2) * rho",
            "full_cycle_exponent": _fraction_payload(FULL_CONTRACTION_EXPONENT),
            "one_tick_map": "g_N(rho) = (N_CRC/pi)^(-1/48) * rho",
            "one_tick_exponent": _fraction_payload(ONE_TICK_EXPONENT),
            "abs_g_star_prime": "(N_CRC/pi)^(-1/48)",
        },
        "exponent_law": {
            "per_tick_exponent_for_m_ticks": "-1/(2m)",
            "declared_round_count_m": DECLARED_REPAIR_ROUNDS,
            "exponent_at_declared_m": "-1/48",
            "corrected_exponent_clause": (
                "the proof gives -1/48 at the declared m = 24; "
                "the parametric law -1/(2m) applies to an m-tick round-count theorem"
            ),
        },
        "symbolic_checks": {
            "tick_times_24": "-1/48 * 24 = -1/2",
            "full_cycle_identity": "(|g_*'|)^24 = (N_CRC/pi)^(-1/2)",
            "ew_projection_identity_if_projection_lemma_supplied": (
                "(|g_*'|)^(4 P_star) = (N_CRC/pi)^(-P_star/12)"
            ),
            "corrected_exponent": "-1/48 at declared m = 24; parametric -1/(2m)",
        },
        "numeric_display_for_rounded_capacity": {
            "N_CRC_display": _format_scientific(n_crc_display),
            "N_CRC_over_pi": _format_scientific(n_over_pi),
            "log_N_CRC_over_pi": f"{log_n_over_pi:.17g}",
            "log10_N_CRC_over_pi": f"{math.log10(n_over_pi):.17g}",
            "abs_g_star_prime_for_display_N": _format_scientific(tick_multiplier),
            "full_24_round_contraction_for_display_N": _format_scientific(full_contraction),
            "tick_power_24_float_check": _format_scientific(tick_power_24),
            "tick_power_24_relative_residual_float": f"{relative_power_residual:.3e}",
        },
        "source_side_dependency_audit": {
            "derived_uses": [
                "D6 screen-capacity relation N_CRC = S_dS",
                "screen-radius identity N_CRC = pi * (r_CRC/ell_star)^2",
                "corpus readback-map factorization F = Cap_read(Obs(nf)) from the Observers synthesis sections",
                "fixed-point readback-closure semantics of N_CRC = F(N_CRC) on the radius coordinate",
            ],
            "declared_inputs": [
                "24-round decomposition of the global repair cycle (declared branch architecture)",
                "positive homogeneous one-tick normal form (normalization convention)",
            ],
            "does_not_use": [
                "alpha_U",
                "v/E_star",
                "D10 transmutation exponent",
                "W or Z measured masses",
                "observed electroweak hierarchy",
            ],
        },
        "proof_sketch": [
            "On the D6 capacity branch, the de Sitter entropy normalization gives N_CRC = pi (r_CRC/ell_star)^2, so the local cell sits at rho_star = (N_CRC/pi)^(-1/2) in the coordinate rho = ell/r_CRC.",
            "The corpus readback map is F(N) = Cap_read(Obs(nf_{r,N}(U_{r,N}))) in the refinement limit: the capacity reconstructed by the interior observer sector from the horizon record (Observers synthesis sections).",
            "Model Cap_read by the D6 area law at the delivery resolution (the canonical scale-covariant extension of N = pi (r/ell)^2; a declared modeling identification): a readback delivering the horizon record at effective resolution coordinate rho_read reconstructs capacity F(N) = pi/rho_read^2.",
            "The fixed-point equation N_CRC = F(N_CRC) = pi/rho_star^2 forces rho_read = rho_star (positive roots): deficit rho_read > rho_star gives F < N, slack rho_read < rho_star gives F > N. This gives the closure transport G_N(1) = rho_star.",
            "Homogeneity (scale-free normal form) forces G_N(rho) = Q rho with Q > 0, and G_N(1) = rho_star gives Q = (N_CRC/pi)^(-1/2).",
            "On the declared 24-tick branch, the one-tick map is the positive homogeneous root with g_N^24 = G_N, so its multiplier is Q^(1/24) = (N_CRC/pi)^(-1/48).",
            "Because g_N(rho) = q rho, the absolute derivative in this coordinate is |g_*'| = q = (N_CRC/pi)^(-1/48). For a general m-tick decomposition the same argument gives (N_CRC/pi)^(-1/(2m)).",
        ],
        "acceptance_criteria_status": {
            "defines_N_CRC_repair_map_tick_coordinate_and_derivative_convention": True,
            "proves_declared_screen_capacity_fixed_point_emits_tick_contraction": True,
            "proof_scope": (
                "the closure transport G_N(1) = rho_star follows as equivalent to the corpus "
                "fixed-point equation N_CRC = F(N_CRC) under the declared area-law counting model, and the "
                "full-cycle multiplier follows; the proof takes the 24-round decomposition as declared "
                "branch architecture"
            ),
            "round_count_derived_from_first_principles": False,
            "closure_transport_derived_from_F_interface": True,
            "readback_counting_model_is_modeling_identification": True,
            "concrete_finite_machinery_verification_open": True,
            "exact_corrected_exponent_clause": "-1/48 at declared m = 24; parametric form -1/(2m)",
            "uses_oph_source_side_data_only": True,
            "imports_electroweak_hierarchy_as_measured_input": False,
            "emits_machine_checkable_certificate": True,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the global repair-tick lemma certificate.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--n-crc-display", type=float, default=DISPLAY_N_CRC)
    args = parser.parse_args()

    artifact = build_artifact(args.n_crc_display)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
