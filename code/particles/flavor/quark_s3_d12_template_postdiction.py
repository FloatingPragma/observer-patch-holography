#!/usr/bin/env python3
"""Evaluate the quarantined S3/D12 quark formula-discovery witness.

This module deliberately contains no quark comparison values.  It reconstructs
the frozen algebraic ansatz from repository artifacts, while preserving their
actual ancestry: most numerical flavor inputs descend from the handwritten
``family_transport_kernel`` template and the pixel candidate has an internal
Stage-5 quark-model ancestor.  The output is therefore a dimensionless,
target-informed diagnostic.  It is not an OPH source-only mass prediction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import mpmath as mp


mp.mp.dps = 80

CODE_ROOT = Path(__file__).resolve().parents[2]
RUNS = CODE_ROOT / "particles" / "runs"
DEFAULT_OUTPUT = RUNS / "flavor" / "quark_s3_d12_template_postdiction.json"

PIXEL_CERTIFICATE = (
    CODE_ROOT
    / "particles"
    / "hierarchy"
    / "certificates"
    / "R_P_source_audit_pixel_certificate.json"
)
FAMILY_KERNEL = RUNS / "flavor" / "family_transport_kernel.json"
EDGE_COCYCLE = RUNS / "flavor" / "overlap_edge_transport_cocycle.json"
ODD_RESPONSE = RUNS / "flavor" / "quark_odd_response_law.json"
CHARGED_BUDGET = RUNS / "flavor" / "charged_budget_transport.json"
P_DERIVATION = CODE_ROOT / "P_derivation" / "paper_math.py"

CLAIM_CLASS = "post_hoc_target_informed_repository_template_ansatz_not_physical_postdiction"


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _m(value: object) -> mp.mpf:
    return mp.mpf(str(value))


def _fmt(value: mp.mpf, digits: int = 40) -> str:
    return mp.nstr(value, n=digits, strip_zeros=False)


def _center(values: Iterable[mp.mpf]) -> list[mp.mpf]:
    vector = list(values)
    mean = mp.fsum(vector) / len(vector)
    return [entry - mean for entry in vector]


def load_repository_inputs() -> tuple[dict[str, object], dict[str, Any]]:
    pixel = _read_json(PIXEL_CERTIFICATE)
    kernel = _read_json(FAMILY_KERNEL)
    edge = _read_json(EDGE_COCYCLE)
    odd = _read_json(ODD_RESPONSE)
    charged = _read_json(CHARGED_BUDGET)

    suppression = edge["derived_pairwise_suppression"]
    charged_certificate = charged["charged_dirac_scalarization_certificate"]
    inputs: dict[str, object] = {
        "P": pixel["P_cand"],
        "alpha_U": pixel["alpha_U_P_cand"],
        "S_13": suppression[0][2],
        "S_23": suppression[1][2],
        "delta_21": edge["hermitian_descendant_riesz_margin"]["gamma"],
        "DeltaS_13": odd["Delta_S_q"][0][2],
        "g_ch_template": charged_certificate["current_family_seed_data"]["seed_value"],
        "N_c": 3,
        "d_std_S3": 2,
    }

    source_files = [
        PIXEL_CERTIFICATE,
        FAMILY_KERNEL,
        EDGE_COCYCLE,
        ODD_RESPONSE,
        CHARGED_BUDGET,
        P_DERIVATION,
    ]
    p_source = P_DERIVATION.read_text(encoding="utf-8")
    provenance = {
        "dependency_files": {
            path.relative_to(CODE_ROOT).as_posix(): {"sha256": _sha256(path)} for path in source_files
        },
        "family_transport_kernel_status": kernel.get("status"),
        "family_transport_kernel_proof_status": kernel.get("proof_status"),
        "family_transport_kernel_note": kernel.get("metadata", {}).get("note"),
        "charged_scale_source_status": charged_certificate.get("stored_shared_absolute_scale_status"),
        "pixel_source_status": pixel.get("status"),
        "evaluation_imports_quark_reference_values": False,
        "formula_discovery_target_informed": True,
        "numerical_flavor_template_consumed": kernel.get("status") == "template",
        "p_source_uses_stage5_quark_model": (
            "quarks = self.diagonal_quark_masses(d10.v)" in p_source
            and 'mass_source="internal_stage5_continuation"' in p_source
        ),
        "dimensionful_scale_emitted": False,
        "common_scale_scheme_declared": False,
        "full_rg_threshold_map_present": False,
        "public_prediction_allowed": False,
    }
    return inputs, provenance


def evaluate(
    inputs: dict[str, object] | None = None,
    *,
    scale_mode: str = "repository_template_gch",
) -> dict[str, Any]:
    if inputs is None:
        inputs, provenance = load_repository_inputs()
    else:
        provenance = {
            "evaluation_imports_quark_reference_values": False,
            "formula_discovery_target_informed": True,
            "numerical_flavor_template_consumed": True,
            "p_source_uses_stage5_quark_model": True,
            "dimensionful_scale_emitted": False,
            "common_scale_scheme_declared": False,
            "full_rg_threshold_map_present": False,
            "public_prediction_allowed": False,
            "note": "Explicit inputs supplied; caller is responsible for retaining their provenance.",
        }

    P = _m(inputs["P"])
    alpha = _m(inputs["alpha_U"])
    S13 = _m(inputs["S_13"])
    S23 = _m(inputs["S_23"])
    delta21 = _m(inputs["delta_21"])
    delta_s13 = _m(inputs["DeltaS_13"])
    n_c = _m(inputs["N_c"])
    d_std = _m(inputs["d_std_S3"])
    n_f = n_c + d_std
    beta_ew = n_c + 1
    phi = (1 + mp.sqrt(5)) / 2
    width = mp.pi * alpha

    # The S3 heat-spectrum identity is exact once tau is supplied.  This
    # identification of tau is a target-informed, unproved ansatz.
    tau_f = P / 4 - width / n_f
    r = mp.exp(-3 * tau_f)
    rho = 3 / (2 + r)
    x = (r - 1) / (r + 1)

    # These edge and mean laws are candidate formulas, not traced OPH theorems.
    edge_u = S13 + rho * delta21 / (1 + rho)
    edge_d = S23 + delta21 / (2 * (1 + rho - x * x))
    mean_u = edge_u - width * delta_s13
    mean_d = edge_d - width * (1 - delta_s13)
    a_u = edge_u + width / d_std
    a_d = edge_d - width / n_f

    linear_basis = _center([-mp.mpf(1), x, mp.mpf(1)])
    quadratic_basis = _center([mp.mpf(1), x * x, mp.mpf(1)])
    # Exact closed ray coordinates in the L/Q basis.  Their symbolic proof is
    # emitted by verify_quark_flavor_source_closure.py.
    b_u_ray = (
        a_u * (-rho * x + rho - x - 1) / ((1 + rho) * (x * x - 1))
    )
    b_d_ray = (
        a_d * (-rho * x - rho - x + 1) / ((1 + rho) * (x * x - 1))
    )
    q_u = -width * rho / (2 * n_f)
    q_d = -width / beta_ew
    b_u = b_u_ray + q_u
    b_d = b_d_ray + q_d

    coefficient_a = 1 / (2 * (1 + rho - x * x))
    coefficient_b = 1 / (2 * (1 - x * x - x * x / (1 + rho)))
    sigma_seed = (mean_u + mean_d) / 2
    eta = (mean_u - mean_d) / 2

    if scale_mode == "repository_template_gch":
        g_ch = _m(inputs["g_ch_template"])
        scale_formula = "dimensionless family-template mean plus minimum gap"
    elif scale_mode == "strict_closure_conjecture":
        g_ch = mp.exp(-width / phi)
        scale_formula = "dimensionless exp(-pi*alpha_U/phi) conjecture"
    else:
        raise ValueError(
            "scale_mode must be 'repository_template_gch' or 'strict_closure_conjecture'"
        )

    g_u = g_ch * mp.exp(-(coefficient_a * sigma_seed - coefficient_b * eta))
    g_d = g_ch * mp.exp(-(coefficient_a * sigma_seed + coefficient_b * eta))
    centered_u = [a_u * linear_basis[i] + b_u * quadratic_basis[i] for i in range(3)]
    centered_d = [a_d * linear_basis[i] + b_d * quadratic_basis[i] for i in range(3)]
    up = [g_u * mp.exp(entry) for entry in centered_u]
    down = [g_d * mp.exp(entry) for entry in centered_d]

    tolerance = mp.mpf("1e-60")
    assert abs(mp.fsum(linear_basis)) < tolerance
    assert abs(mp.fsum(quadratic_basis)) < tolerance
    assert abs(mp.fsum(centered_u)) < tolerance
    assert abs(mp.fsum(centered_d)) < tolerance

    return {
        "artifact": "oph_quark_s3_d12_template_postdiction_v2",
        "claim_class": CLAIM_CLASS,
        "promotion_allowed": False,
        "scale_mode": scale_mode,
        "unit_status": {
            "emitted_coordinates": "dimensionless",
            "bundle_comparison_interpretation": "silently treats one coordinate unit as 1 GeV",
            "dimensionful_bridge": "missing",
            "physical_mass_or_yukawa_claim_allowed": False,
        },
        "repository_inputs": inputs,
        "provenance": provenance,
        "exact_mathematics": {
            "S3_transposition_Cayley_Laplacian_spectrum": {"0": 1, "3": 4, "6": 1},
            "heat_gap_ratio_given_tau": "r=exp(-3*tau)",
            "LQ_Gram_determinant": "4*(1-x^2)^2/3",
            "LQ_scope": "universal basis of the centered three-vector plane for x != +/-1",
        },
        "unproved_formula_choices": {
            "heat_time": "tau_f=P/4-pi*alpha_U/5",
            "edge_laws": [
                "edge_u=S_13+rho*delta_21/(1+rho)",
                "edge_d=S_23+delta_21/(2*(1+rho-x^2))",
            ],
            "mean_allocation": "DeltaS_13 and its complement",
            "odd_even_denominator_tuple": [5, 2, 5, 10, 4],
            "affine_mean_law": "repository candidate A_ud/B_ud law",
            "scale_formula": scale_formula,
        },
        "derived_scalars": {
            "tau_f": _fmt(tau_f),
            "r": _fmt(r),
            "rho": _fmt(rho),
            "x": _fmt(x),
            "edge_u": _fmt(edge_u),
            "edge_d": _fmt(edge_d),
            "mean_u": _fmt(mean_u),
            "mean_d": _fmt(mean_d),
            "a_u": _fmt(a_u),
            "a_d": _fmt(a_d),
            "b_u_ray": _fmt(b_u_ray),
            "b_d_ray": _fmt(b_d_ray),
            "q_u": _fmt(q_u),
            "q_d": _fmt(q_d),
            "b_u": _fmt(b_u),
            "b_d": _fmt(b_d),
            "A_ud": _fmt(coefficient_a),
            "B_ud": _fmt(coefficient_b),
            "g_ch": _fmt(g_ch),
            "g_u": _fmt(g_u),
            "g_d": _fmt(g_d),
        },
        "basis": {
            "L=ctr(-1,x,1)": [_fmt(entry) for entry in linear_basis],
            "Q=ctr(1,x^2,1)": [_fmt(entry) for entry in quadratic_basis],
        },
        "centered_log_coordinates": {
            "up": [_fmt(entry) for entry in centered_u],
            "down": [_fmt(entry) for entry in centered_d],
        },
        "dimensionless_output_coordinates": {
            "u": _fmt(up[0]),
            "c": _fmt(up[1]),
            "t": _fmt(up[2]),
            "d": _fmt(down[0]),
            "s": _fmt(down[1]),
            "b": _fmt(down[2]),
        },
        "blocking_conditions": [
            "replace the handwritten family kernel with an OPH-emitted observer-like self-reading carrier receipt",
            "derive invariant, normalization-covariant edge and odd/even response scalars",
            "remove the Stage-5 quark-model ancestor from the pixel source branch",
            "emit a dimensionful normalization or same-scale dimensionless Yukawa normalization",
            "supply one common renormalization scheme and the full RG/threshold map",
            "freeze every formula and selector prospectively before comparison",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scale-mode",
        choices=["repository_template_gch", "strict_closure_conjecture"],
        default="repository_template_gch",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--allow-template-ancestry",
        action="store_true",
        help="Acknowledge that this command evaluates a handwritten-template diagnostic.",
    )
    parser.add_argument("--print-json", action="store_true")
    args = parser.parse_args()

    if not args.allow_template_ancestry:
        parser.error(
            "refusing to evaluate a template-ancestry candidate without "
            "--allow-template-ancestry"
        )
    artifact = evaluate(scale_mode=args.scale_mode)
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
