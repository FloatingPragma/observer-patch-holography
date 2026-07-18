#!/usr/bin/env python3
"""Certify non-identifiability of the scoped current OPH charged signature.

The exact theorem is deliberately scoped to D9 charged-channel admissibility,
fixed D10 P and v(P), and the invariant 12/24 structural receipts.  It excludes
mass-dependent electromagnetic endpoint transport.  On this scope, the three
spectral coordinates of a charged Yukawa matrix are free. Separate audits
show why the MaxEnt, phase, family-attachment, and Z6 proposals do not
select those coordinates.

No charged reference packet is read. Target-informed constants
appear only as examples being audited, never as theorem inputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_current_oph_moduli_independence.json"
)

RECEIPTS = (
    ROOT / "particles" / "runs" / "calibration" / "d10_ew_forward_transmutation_certificate.json",
    ROOT / "particles" / "hierarchy" / "certificates" / "R_screen_sieve_icosahedral_certificate.json",
    ROOT / "particles" / "hierarchy" / "certificates" / "R_m_rep_24_certificate.json",
    ROOT / "particles" / "runs" / "leptons" / "charged_z3_maxent_balance_audit.json",
    ROOT / "particles" / "runs" / "leptons" / "charged_z3_phase_holonomy_no_go.json",
    ROOT / "particles" / "runs" / "leptons" / "charged_absolute_scale_underdetermination_theorem.json",
    ROOT / "particles" / "runs" / "leptons" / "charged_trace_lift_theorem.json",
)


def koide_from_action_gap(action_gap: float) -> float:
    """Natural trace base: E_c/E_0 = 2 exp(-action_gap)."""

    return (1.0 + 2.0 * math.exp(-action_gap)) / 3.0


def balanced_roots(delta: float) -> tuple[float, float, float]:
    """Sorted roots of the balanced Hermitian C3 circulant."""

    return tuple(sorted(
        1.0 + math.sqrt(2.0) * math.cos(delta + 2.0 * math.pi * k / 3.0)
        for k in range(3)
    ))


def square_root_mass_invariant(values: Iterable[float]) -> float:
    masses = tuple(float(value) for value in values)
    return sum(masses) / sum(math.sqrt(value) for value in masses) ** 2


def determinant_one(values: Iterable[float]) -> tuple[float, ...]:
    data = tuple(float(value) for value in values)
    geometric_mean = math.prod(data) ** (1.0 / len(data))
    return tuple(value / geometric_mean for value in data)


def centered_logs(values: Iterable[float]) -> tuple[float, ...]:
    logs = tuple(math.log(float(value)) for value in values)
    mean = sum(logs) / len(logs)
    return tuple(value - mean for value in logs)


def matrix_rank(rows: Iterable[Iterable[float]], tolerance: float = 1.0e-12) -> int:
    """Small deterministic Gaussian-elimination rank used by the receipt."""

    matrix = [list(map(float, row)) for row in rows]
    if not matrix:
        return 0
    row_count = len(matrix)
    column_count = len(matrix[0])
    rank = 0
    for column in range(column_count):
        pivot = max(range(rank, row_count), key=lambda row: abs(matrix[row][column]))
        if abs(matrix[pivot][column]) <= tolerance:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for row in range(row_count):
            if row == rank:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                value - factor * pivot_entry
                for value, pivot_entry in zip(matrix[row], matrix[rank], strict=True)
            ]
        rank += 1
        if rank == row_count:
            break
    return rank


def spectral_family(mu: float, x: float, y: float) -> dict[str, Any]:
    """Three independent log-singular-value coordinates.

    The centered plane is parametrized by (x,y,-x-y); mu is the common affine
    coordinate.  D9 charged-channel admissibility permits every such positive
    diagonal matrix, before optional left/right unitary changes of frame.
    """

    logs = (mu + x, mu + y, mu - x - y)
    singular_values = tuple(math.exp(value) for value in logs)
    centered = tuple(value - mu for value in logs)
    return {
        "mu": mu,
        "x": x,
        "y": y,
        "log_singular_values": logs,
        "centered_log_singular_values": centered,
        "singular_values": singular_values,
        "product_singular_values": math.prod(singular_values),
        "log_abs_determinant": sum(logs),
        "positive_simple_spectrum": len({round(value, 14) for value in singular_values}) == 3,
    }


def _receipt_bindings() -> list[dict[str, Any]]:
    bindings: list[dict[str, Any]] = []
    for path in RECEIPTS:
        raw = path.read_bytes()
        payload = json.loads(raw)
        bindings.append(
            {
                "path": path.relative_to(ROOT).as_posix(),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact": payload.get("artifact", payload.get("certificate_id")),
                "status": payload.get("status", payload.get("proof_status")),
            }
        )
    return bindings


def _source_scope_checks() -> dict[str, bool]:
    d10, screen, round_count, maxent, phase, absolute, trace_lift = (
        json.loads(path.read_text(encoding="utf-8")) for path in RECEIPTS
    )
    return {
        "d10_forward_map_closed": d10.get("status") == "closed_forward_p_to_t_map",
        "conditional_screen_branch_has_twelve_uniform_unit_defects": (
            screen.get("orbit_stabilizer", {}).get("orbit_size") == 12
            and screen.get("strict_unit_defect_minimum", {}).get("charges") == [1] * 12
        ),
        "oriented_round_count_is_24": round_count.get("result", {}).get("m_rep") == 24,
        "maxent_does_not_promote_koide": maxent.get("public_koide_promotion_allowed") is False,
        "phase_receipt_does_not_promote_2_over_9": phase.get("public_phase_promotion_allowed") is False,
        "centered_scale_receipt_does_not_promote": absolute.get("public_promotion_allowed") is False,
        "trace_lift_receipt_does_not_promote": (
            trace_lift.get("promotion", {}).get("public_promotion_allowed") is False
        ),
    }


def _phase_countermodel(delta: float) -> dict[str, Any]:
    roots = balanced_roots(delta)
    masses = determinant_one(root * root for root in roots)
    return {
        "delta_radians": delta,
        "roots": roots,
        "positive_simple_spectrum": all(root > 0.0 for root in roots)
        and len({round(root, 14) for root in roots}) == 3,
        "determinant_normalized_masses": masses,
        "mass_ratios_to_lightest": tuple(value / masses[0] for value in masses),
        "Q": square_root_mass_invariant(masses),
        "loop_phase_3delta_mod_2pi": (3.0 * delta) % (2.0 * math.pi),
    }


def build_artifact() -> dict[str, Any]:
    phase_a = _phase_countermodel(0.1)
    phase_b = _phase_countermodel(2.0 / 9.0)
    family_base = spectral_family(0.0, -2.0, 0.25)
    family_x = spectral_family(0.0, -1.5, 0.25)
    family_y = spectral_family(0.0, -2.0, 0.75)
    family_mu = spectral_family(0.7, -2.0, 0.25)
    scale_base = determinant_one(value * value for value in balanced_roots(2.0 / 9.0))
    scale_shifted = tuple(2.0 * value for value in scale_base)
    exponents_a = (7, 4, 3)
    exponents_b = tuple(value + 6 for value in exponents_a)
    shape_columns = [[1.0, 0.0], [0.0, 1.0], [-1.0, -1.0]]
    full_columns = [[1.0, 1.0, 0.0], [1.0, 0.0, 1.0], [1.0, -1.0, -1.0]]
    scope_checks = _source_scope_checks()
    return {
        "artifact": "oph_charged_current_oph_moduli_independence",
        "status": "CLOSED_SCOPED_CHARGED_SIGNATURE_NONIDENTIFIABILITY_CERTIFICATE",
        "runtime_charged_reference_packet_consumed": False,
        "historical_target_informed_examples_present": True,
        "public_charged_mass_promotion_allowed": False,
        "source_receipt_bindings": _receipt_bindings(),
        "source_scope_checks": scope_checks,
        "source_scope_checks_pass": all(scope_checks.values()),
        "scope": {
            "included": [
                "G_SM=(SU(3)xSU(2)xU(1))/Z6 and D9 charged-channel admissibility",
                "N_c=N_g=3, one Higgs doublet, and fixed D10 P and v(P)",
                "the invariant twelve-port screen and count-only 24-slot register",
                "the named MaxEnt, phase, centered-scale, and trace-lift audit receipts",
            ],
            "excluded": [
                "mass-dependent electromagnetic vacuum-polarization transport",
                "the Thomson endpoint and its moving P(kappa) fixed point",
                "empirical hadronic payloads, target-anchored charged readouts, and thresholds",
                "any physical charged Yukawa matrix supplied as an input",
            ],
            "why_excluded": (
                "Electromagnetic endpoint transport changes under a common charged rescaling. "
                "It can therefore break the scale freedom, but its closure depends on empirical "
                "inputs or source inputs whose derivation is work in progress; it is not an "
                "antecedent of the charged common-shift theorem."
            ),
        },
        "theorem": {
            "statement": (
                "At the charged-readout layer, D9 identifies the allowed charged Yukawa channel "
                "but leaves its generation factor Mat_3(C) unselected. Fixed P and v(P), plus "
                "the invariant 12/24 receipts, add no map choosing a point or singular-value orbit "
                "in that space. The explicit family diag(exp(mu+x),exp(mu+y),exp(mu-x-y)) "
                "therefore gives two independent centered shape directions and one independent "
                "common affine direction compatible with this scoped signature. The premises do "
                "not uniquely entail, or supply a natural selector for, a charged mass triple."
            ),
            "coordinate_count": (
                "Three real spectral coordinates are required: two centered log-shape coordinates "
                "and one common affine/determinant coordinate, followed by a physical family-line "
                "attachment and a mass-scheme map."
            ),
            "categorical_form": (
                "D9 isolates the charged isotypic channel C tensor Mat_3(C), not a point or "
                "singular-value orbit in Mat_3(C). The nonzero elements of its determinant line "
                "form a C* torsor. A bare one-dimensional line has no natural nonzero normalized "
                "section: naturality under every scalar automorphism would force a putative "
                "section to vanish. The physical scale uses log|det Y_e|, equivalently the sum "
                "of log singular values."
            ),
        },
        "independent_spectral_coordinate_family": {
            "formula": "Y(mu,x,y)=diag(exp(mu+x),exp(mu+y),exp(mu-x-y))",
            "shape_jacobian_rows": shape_columns,
            "shape_jacobian_rank": matrix_rank(shape_columns),
            "full_log_jacobian_rows": full_columns,
            "full_log_jacobian_rank": matrix_rank(full_columns),
            "base": family_base,
            "x_direction": family_x,
            "y_direction": family_y,
            "mu_direction": family_mu,
            "preservation_statement": (
                "All four matrices are positive diagonal elements of the D9-allowed Mat_3(C) "
                "generation factor. Varying x or y preserves the product of singular values; "
                "varying mu preserves centered log gaps and multiplies that product by exp(3 mu)."
            ),
        },
        "distinct_failed_selector_audits": {
            "maxent_balance": {
                "law": "E_c/E_0=2*exp(-(S_c-S_0)); Q=(1+E_c/E_0)/3",
                "equal_action": {
                    "S_c_minus_S_0": 0.0,
                    "Q": koide_from_action_gap(0.0),
                },
                "one_bit_compensated": {
                    "S_c_minus_S_0": "ln(2)",
                    "Q": koide_from_action_gap(math.log(2.0)),
                },
                "conclusion": (
                    "This is an algebraic selector audit, not a construction of two full OPH "
                    "models. The current relative-MaxEnt receipt leaves the charged-block action "
                    "gap and the ensemble-to-deterministic-Yukawa bridge unspecified."
                ),
            },
            "balanced_phase": {
                "countermodel_A": phase_a,
                "countermodel_B": phase_b,
                "conclusion": (
                    "Both determinant-normalized spectra are positive, simple, and exactly "
                    "Koide-balanced, but their ratios and loop phases differ."
                ),
            },
            "common_affine_scale": {
                "Y0_singular_values": scale_base,
                "Y1_equals_2Y0_singular_values": scale_shifted,
                "centered_logs_Y0": centered_logs(scale_base),
                "centered_logs_Y1": centered_logs(scale_shifted),
                "determinant_Y0": math.prod(scale_base),
                "determinant_Y1": math.prod(scale_shifted),
                "conclusion": (
                    "Y_e -> exp(kappa)Y_e preserves centered log operators, their projectors, "
                    "centered log gaps, ratios, and the scoped centered refinement data while "
                    "multiplying |det(Y_e)| by exp(3*kappa). It does not preserve absolute "
                    "spectral gaps or the excluded electromagnetic endpoint transport."
                ),
            },
            "family_attachment": {
                "freedom": "Y_e -> U_L Y_e U_R^dagger for U_L,U_R in U(3)",
                "finite_subfamily": "the six permutation attachments",
                "preserved": [
                    "singular values",
                    "Koide invariant",
                    "determinant magnitude",
                    "D9 gauge representations and anomaly data",
                ],
                "changed": (
                    "attachment relative to a separately fixed physical family observable; "
                    "without such an observable these unitary changes may be basis redundancy"
                ),
                "boundary": (
                    "Attachment is unnecessary for an unordered mass triple, but necessary to "
                    "identify its lines consistently with neutrino/mixing or other family data."
                ),
            },
            "z6_defect_lengths": {
                "exponents_A": exponents_a,
                "exponents_B": exponents_b,
                "residues_mod_6_A": tuple(value % 6 for value in exponents_a),
                "residues_mod_6_B": tuple(value % 6 for value in exponents_b),
                "sum_A": sum(exponents_a),
                "sum_B": sum(exponents_b),
                "determinant_power_change": sum(exponents_b) - sum(exponents_a),
                "conclusion": (
                    "Conditional audit: even granting the proposed identification "
                    "of suppression exponents with Z6 residue lifts, quotient data would not "
                    "choose the integer lifts. Current OPH does not first derive that identification."
                ),
            },
        },
        "extended_centered_candidate_scope": {
            "statement": (
                "If one additionally grants a fixed theorem-grade centered charged readout, the "
                "two shape coordinates cease to be free, but the exact common-shift family "
                "Y_e -> exp(kappa)Y_e is free. The on-disk centered operator is a "
                "template-dependent candidate and does not satisfy its promotion gate."
            ),
            "remaining_spectral_modulus": "mu=(1/3)log|det Y_e|",
        },
        "independent_inputs_not_supplied": [
            "charged-block action gap ln(2) plus ensemble-to-Yukawa isometry",
            "physical C3 family connection and loop law fixing delta=2/9",
            "physical root-to-family attachment",
            "absolute excitation/determinant law det(Y_e)=sqrt(2)/6^14",
            "renormalization coordinate and QED/electroweak pole-mass conversion",
        ],
        "smallest_honest_positive_routes": [
            "derive a physical uncentered charged response Y_e(P) directly",
            "derive a quotient-visible refinement-stable W5 port/repair record for the two shape coordinates and close the determinant line separately",
        ],
        "claim_boundary": (
            "This certificate proves non-identifiability only for the explicitly scoped charged "
            "signature and makes no claim outside that scope. A closing extension must "
            "add a source object that is nonconstant along the displayed spectral family, bind "
            "its receipt scope, and freeze it before charged-mass comparison."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": artifact["status"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
