#!/usr/bin/env python3
"""Verify the conditional charged face-incidence fixed-point theorem.

This module consumes no charged reference masses.  It proves consequences of
the explicitly declared three-coordinate affine repair kernel used by the
retrospective face-incidence candidate.  The kernel and numerical branch are
historically target-informed and are not derived from the OPH axioms.

The exact conclusion is therefore conditional: the declared map is a strict
contraction with one fixed point, and the displayed determinant-normalized
spectrum follows after all physical attachment and scheme hypotheses are
granted.  Explicit source-multiplier witnesses show that the listed symmetry,
analytic, diagonal, and contraction properties do not select the unit source
normalization.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import mpmath as mp


mp.mp.dps = 100

HERE = Path(__file__).resolve()
CODE_ROOT = HERE.parents[2]
FACE_RECEIPT = (
    CODE_ROOT / "particles" / "runs" / "leptons"
    / "charged_icosahedral_face_carrier_frontier.json"
)
DEFAULT_OUT = (
    CODE_ROOT / "particles" / "runs" / "leptons"
    / "charged_face_incidence_conditional_theorem.json"
)

V1_ARCHIVE_SHA256 = "4dff09dbc66368a99dbddac4e238e3dcd42906dada7d2f084ff059dda32b04ef"
THEOREM_ARCHIVE_SHA256 = "cc550b311c5f2bc416d78f65b4ec9e1bf32f10f9b73fc4135fd05de626b48d4f"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def text(value: mp.mpf | mp.mpc, digits: int = 80) -> str:
    return mp.nstr(value, digits)


def build_artifact(face_receipt: dict[str, Any], face_receipt_sha256: str) -> dict[str, Any]:
    if face_receipt.get("checks_pass") is not True:
        raise ValueError("conditional theorem requires the closed geometric face receipt")
    if face_receipt.get("public_charged_mass_promotion_allowed") is not False:
        raise ValueError("geometric receipt must retain the physical-value boundary")

    # Frozen historically target-informed tuple declared by the submitted map.
    p_value = mp.mpf("1.6309681897")
    v_gev = mp.mpf("246.6174823334856")
    alpha_u = mp.mpf("0.041124336195630495")

    n_c = 3
    n_g = 3
    beta_ew = 4
    vertices, edges, faces = 12, 30, 20
    euler = vertices - edges + faces
    vertex_degree = 5
    face_degree = 3
    phi = (1 + mp.sqrt(5)) / 2
    alpha_p = (p_value - phi) / mp.sqrt(mp.pi)

    s_kappa = alpha_u / (edges + faces)
    q_kappa = (
        alpha_u / (edges + 1)
        + alpha_u**2 / ((edges + 1) * (vertices - n_g + 1))
    )
    s_chi = alpha_u**2 / (2 ** (vertices - n_g))
    q_chi = alpha_u / ((vertex_degree + euler) * (vertices - 1))
    s_zeta = alpha_p**2 / (faces + 1)
    q_zeta = (
        alpha_u / (edges - n_g)
        + alpha_u**2 / ((edges - n_g) * vertex_degree)
    )

    kappa = s_kappa / (1 + q_kappa)
    chi = s_chi / (1 - q_chi)
    zeta = s_zeta / (1 - q_zeta)
    contraction = max(q_kappa, q_chi, q_zeta)

    def repair(theta: tuple[mp.mpf, mp.mpf, mp.mpf]) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
        k, c, z = theta
        return (
            s_kappa - q_kappa * k,
            s_chi + q_chi * c,
            s_zeta + q_zeta * z,
        )

    fixed_point = (kappa, chi, zeta)
    fixed_image = repair(fixed_point)
    residual = [fixed_image[index] - fixed_point[index] for index in range(3)]

    delta0 = mp.mpf(beta_ew) / (2 * n_c * n_g)
    delta = delta0 + zeta
    rho = mp.sqrt(2) * mp.e ** (-chi)
    roots = [
        1 + rho * mp.cos(delta + 2 * mp.pi * index / 3)
        for index in range(3)
    ]
    reversed_roots = [
        1 + rho * mp.cos(-delta + 2 * mp.pi * index / 3)
        for index in range(3)
    ]
    root_squares = sorted(root * root for root in roots)
    shape_geometric_mean = mp.power(mp.fprod(root_squares), mp.mpf(1) / 3)
    determinant_exponent = vertices + euler
    bare_determinant = v_gev**3 / (2 * mp.power(6, determinant_exponent))
    endpoint_geometric_mean = mp.power(bare_determinant, mp.mpf(1) / 3) * mp.e**kappa
    masses_gev = [
        endpoint_geometric_mean * square / shape_geometric_mean
        for square in root_squares
    ]
    determinant_residual = (
        mp.fprod(masses_gev)
        / (bare_determinant * mp.e ** (3 * kappa))
        - 1
    )

    def multiplier_witness(
        lambda_kappa: int,
        lambda_chi: int,
        lambda_zeta: int,
    ) -> dict[str, Any]:
        witness_kappa = lambda_kappa * s_kappa / (1 + q_kappa)
        witness_chi = lambda_chi * s_chi / (1 - q_chi)
        witness_zeta = lambda_zeta * s_zeta / (1 - q_zeta)
        witness_rho = mp.sqrt(2) * mp.e ** (-witness_chi)
        witness_delta = delta0 + witness_zeta
        witness_roots = [
            1 + witness_rho * mp.cos(witness_delta + 2 * mp.pi * index / 3)
            for index in range(3)
        ]
        witness_squares = sorted(root * root for root in witness_roots)
        witness_shape_mean = mp.power(mp.fprod(witness_squares), mp.mpf(1) / 3)
        witness_mean = (
            mp.power(bare_determinant, mp.mpf(1) / 3)
            * mp.e**witness_kappa
        )
        witness_masses_mev = [
            1000 * witness_mean * square / witness_shape_mean
            for square in witness_squares
        ]
        return {
            "source_multipliers": [lambda_kappa, lambda_chi, lambda_zeta],
            "fixed_point": [
                text(witness_kappa),
                text(witness_chi),
                text(witness_zeta),
            ],
            "masses_mev": [text(value) for value in witness_masses_mev],
            "same_jacobian_and_contraction": True,
        }

    countermodels = [
        multiplier_witness(2, 1, 1),
        multiplier_witness(1, 2, 1),
        multiplier_witness(1, 1, 2),
    ]
    base_mass_text = [text(1000 * value) for value in masses_gev]
    countermodels_change_masses = all(
        row["masses_mev"] != base_mass_text for row in countermodels
    )
    orientation_residual = [
        left - right
        for left, right in zip(sorted(roots), sorted(reversed_roots), strict=True)
    ]
    bare_ratio = 1 / mp.sqrt(2)
    endpoint_ratio = mp.e ** (-chi) / mp.sqrt(2)

    checks = {
        "structural_incidence_is_12_30_20": (
            vertices, edges, faces, euler, vertex_degree, face_degree
        ) == (12, 30, 20, 2, 5, 3),
        "strict_contraction": contraction < 1,
        "fixed_point_residual_below_1e_90": max(abs(value) for value in residual) < mp.mpf("1e-90"),
        "positive_simple_roots": (
            min(roots) > 0
            and min(
                abs(roots[left] - roots[right])
                for left in range(3)
                for right in range(left)
            ) > mp.mpf("1e-20")
        ),
        "orientation_reversal_preserves_spectrum": max(
            abs(value) for value in orientation_residual
        ) < mp.mpf("1e-90"),
        "determinant_identity": abs(determinant_residual) < mp.mpf("1e-90"),
        "explicit_source_multiplier_witnesses_change_masses": countermodels_change_masses,
        "bare_and_endpoint_amplitude_ratios_are_distinct": endpoint_ratio != bare_ratio,
    }

    return {
        "artifact": "oph_charged_face_incidence_conditional_theorem",
        "status": "CONDITIONAL_FIXED_POINT_AND_READOUT_CLOSED_SOURCE_LAW_RIGIDITY_OPEN",
        "runtime_charged_reference_consumed": False,
        "historical_charged_target_informed": True,
        "global_source_only": False,
        "branch_tuple_coherent": False,
        "mass_scheme_certified": False,
        "public_prediction_promotion_allowed": False,
        "runtime_dependency": {"mpmath": mp.__version__},
        "provenance": {
            "geometric_face_receipt": "charged_icosahedral_face_carrier_frontier.json",
            "geometric_face_receipt_sha256": face_receipt_sha256,
            "candidate_v1_archive_sha256": V1_ARCHIVE_SHA256,
            "conditional_theorem_archive_sha256": THEOREM_ARCHIVE_SHA256,
        },
        "declared_hybrid_inputs": {
            "P": text(p_value),
            "v_gev": text(v_gev),
            "alpha_U": text(alpha_u),
            "warning": (
                "P and v belong to the truncated D10 comparison probe, while alpha_U "
                "belongs to the distinct canonical public pixel."
            ),
        },
        "conditional_hypotheses": [
            "an equivariant physical generation attachment to the A5/C3 face-corner bundle",
            "a source-derived bare block balance and an explicit bare-to-endpoint amplitude-repair bridge",
            "a charged connection deriving the 2/9 base phase and the endpoint phase correction",
            "a normalized primitive charged Z6 determinant character attached to the physical determinant line",
            "the displayed affine finite-incidence repair kernel is selected by an OPH source action and is refinement-natural",
            "one coherent receipt-bound P-v-alpha_U source branch is frozen before charged comparison",
            "the terminal charged response is attached to a declared running coordinate and controlled pole-scheme map",
        ],
        "repair_map": {
            "formula": "T(k,c,z)=(s_k-q_k*k, s_c+q_c*c, s_z+q_z*z)",
            "sources": {
                "s_kappa": text(s_kappa),
                "s_chi": text(s_chi),
                "s_zeta": text(s_zeta),
            },
            "feedback": {
                "q_kappa": text(q_kappa),
                "q_chi": text(q_chi),
                "q_zeta": text(q_zeta),
            },
            "jacobian_diagonal": [
                text(-q_kappa),
                text(q_chi),
                text(q_zeta),
            ],
            "contraction_constant": text(contraction),
            "unique_iteratively_stable_fixed_point": [
                text(kappa),
                text(chi),
                text(zeta),
            ],
            "fixed_point_residual": [text(value) for value in residual],
            "stability_scope": (
                "Banach stability applies to iteration of this declared three-variable "
                "map, not to regulator refinement, radiative stability, or source-law selection."
            ),
        },
        "bare_to_endpoint_stage_boundary": {
            "bare_balanced_amplitude_ratio": text(bare_ratio),
            "endpoint_amplitude_ratio": text(endpoint_ratio),
            "endpoint_over_bare": text(mp.e ** (-chi)),
            "conclusion": (
                "Exact bare balance and the nonzero endpoint amplitude defect cannot be "
                "simultaneous statements about one un-staged operator. A physical repair "
                "bridge between the two stages remains required."
            ),
        },
        "conditional_spectrum": {
            "delta0": text(delta0),
            "delta": text(delta),
            "rho": text(rho),
            "roots": [text(value) for value in roots],
            "masses_mev": [text(1000 * value) for value in masses_gev],
            "determinant_relative_residual": text(determinant_residual),
        },
        "scoped_source_normalization_non_rigidity": {
            "statement": (
                "The three explicit integer source-multiplier witnesses retain the "
                "displayed diagonal Jacobian and contraction while changing the masses. "
                "Thus those enumerated admissibility properties do not select unit source "
                "normalization. This is not a non-entailment theorem for every OPH axiom."
            ),
            "witnesses": countermodels,
        },
        "proved": [
            "strict contraction, unique fixed point, and global iteration convergence for the declared affine map",
            "determinant-normalized circulant spectrum and mass readout conditional on all listed physical hypotheses",
            "non-rigidity of unit source normalization under the enumerated candidate properties",
        ],
        "open": [
            "OPH derivation and rigidity of the repair sources and feedback kernel",
            "bare-to-endpoint block-balance repair",
            "physical generation attachment and charged phase connection",
            "normalized determinant-line character",
            "refinement naturality and coherent no-target-leak source branch",
            "running, threshold, and pole-scheme attachment",
        ],
        "checks": checks,
        "checks_pass": all(checks.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--face-receipt", type=Path, default=FACE_RECEIPT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    face_raw = args.face_receipt.read_bytes()
    artifact = build_artifact(json.loads(face_raw), sha256(face_raw))
    output = (json.dumps(artifact, indent=2, sort_keys=True) + "\n").encode()
    if args.check:
        actual = args.out.read_bytes() if args.out.exists() else None
        ok = actual == output
        print(json.dumps({"status": "OK" if ok else "DRIFT"}, indent=2))
        return 0 if ok else 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_bytes(output)
    print(json.dumps({"status": artifact["status"], "checks_pass": artifact["checks_pass"]}, indent=2))
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
