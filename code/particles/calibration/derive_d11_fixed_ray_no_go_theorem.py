#!/usr/bin/env python3
"""Emit the fixed-ray point statement for the one-scalar D11 branch.

Chain role: record where the reference Higgs/top central pair sits relative to
the one-scalar D11 fixed ray on the declared D10/D11 surface, with the
reference uncertainties propagated through the same Jacobian.

Mathematics: the fixed ray forces `pi_y = pi_lambda`, hence `w_HT = 0`. The
compare-only central pair (Q007TP4 top, S126M Higgs) inverted through the
declared linear Jacobian has `w_HT = -0.00249`, so that central point is off
the ray. Propagating the declared reference errors through the same Jacobian
gives `sigma(w_HT) = 0.00366`, a 0.68 sigma pull, so the one-scalar ray stays
data-compatible on this surface. The two-coordinate split is a chosen
extension; it is recorded as such.

Inputs: the declared D10/D11 calibration surface, the non-promoting one-scalar
fixed-ray seed, and the compare-only exact inverse slice used only as a
witness (it carries the reference uncertainties and codomain ids).

Output: a machine-readable point-statement artifact with the uncertainty
pull, the codomain dependence, and the chosen extension stated explicitly.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.artifact_paths import canonicalize_artifact_paths
D11_SURFACE_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_declared_calibration_surface.json"
D11_FORWARD_SEED_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed.json"
D11_EXACT_ADAPTER_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_reference_exact_adapter.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d11_fixed_ray_no_go_theorem.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(d11_surface: dict, forward_seed: dict, exact_adapter: dict) -> dict:
    core = dict(d11_surface["core"])
    jacobian = dict(d11_surface["jacobian"])
    sigma = float(forward_seed["sigma_D11_HT"])

    y_core = float(core["y_t_core_mt"])
    lambda_core = float(core["lambda_core_mt"])
    mt_core = float(core["mt_pole_core_gev"])
    mh_core = float(core["mH_core_gev"])
    d_mt = float(jacobian["d_mt_pole_d_y_t"])
    d_mh = float(jacobian["d_mH_d_lambda"])

    exact_targets = dict(exact_adapter["exact_reference_targets"])
    exact_slice = dict(exact_adapter["inverse_slice_coordinates"])
    exact_readback = dict(exact_adapter["normalized_readback"])
    reference_sigma = dict(exact_adapter["reference_uncertainties_gev"])
    reference_codomain = dict(exact_adapter["reference_codomain"])
    aux_codomain = dict(exact_adapter["auxiliary_direct_top_codomain"])

    pi_y_exact = float(exact_readback["pi_y"])
    pi_lambda_exact = float(exact_readback["pi_lambda"])
    sigma_exact = 0.5 * (pi_y_exact + pi_lambda_exact)
    eta_exact = 0.5 * (pi_y_exact - pi_lambda_exact)
    w_exact = pi_y_exact - pi_lambda_exact

    compatibility_functional = (
        (float(exact_targets["mt_pole_gev"]) - mt_core) / (d_mt * y_core)
        + (9.0 / 16.0) * (float(exact_targets["mH_gev"]) - mh_core) / (d_mh * lambda_core)
    )

    # Linear error propagation of the declared reference uncertainties through
    # the same Jacobian that defines the compare-only inverse slice.
    sigma_pi_y = float(reference_sigma["mt_pole_gev"]) / (d_mt * y_core)
    sigma_pi_lambda = (9.0 / 16.0) * float(reference_sigma["mH_gev"]) / (d_mh * lambda_core)
    sigma_w = math.sqrt(sigma_pi_y**2 + sigma_pi_lambda**2)
    pull_sigma = w_exact / sigma_w
    within_one_sigma = abs(pull_sigma) < 1.0

    # The fixed ray evaluated at the reference Higgs central value predicts one
    # top pole mass; its distance from the reference top value in reference
    # sigmas is the same 0.68 sigma statement read on the top axis.
    mt_on_ray_at_mh_central = mt_core + d_mt * y_core * pi_lambda_exact
    mt_ray_offset_gev = mt_on_ray_at_mh_central - float(exact_targets["mt_pole_gev"])
    mt_ray_offset_sigma = mt_ray_offset_gev / float(reference_sigma["mt_pole_gev"])

    # Codomain dependence: the same test against the auxiliary direct-top
    # extraction codomain (Q007TP) instead of the cross-section codomain (Q007TP4).
    pi_y_aux = (float(aux_codomain["mt_gev"]) - mt_core) / (d_mt * y_core)
    w_aux = pi_y_aux - pi_lambda_exact
    sigma_pi_y_aux = float(aux_codomain["sigma_gev"]) / (d_mt * y_core)
    sigma_w_aux = math.sqrt(sigma_pi_y_aux**2 + sigma_pi_lambda**2)
    pull_aux = w_aux / sigma_w_aux

    primary_id = str(reference_codomain["mt_pole_summary_id"])
    higgs_id = str(reference_codomain["mH_summary_id"])
    aux_id = str(aux_codomain["summary_id"])

    return {
        "artifact": "oph_d11_fixed_ray_no_go_theorem",
        "generated_utc": _timestamp(),
        "theorem_id": "D11FixedRayNoGoTheorem",
        "proof_status": "central_pair_off_fixed_ray_on_declared_linear_surface_within_one_sigma",
        "status": "closed_point_statement",
        "scope": "declared_d10_d11_surface_against_compare_only_central_pair",
        "theorem_statement": (
            f"On the declared linear D10/D11 surface the one-scalar fixed ray pi_y = pi_lambda has w_HT = 0 identically. "
            f"The {primary_id}/{higgs_id} central pair (mt = {float(exact_targets['mt_pole_gev'])}, mH = {float(exact_targets['mH_gev'])}) "
            f"inverted through the declared Jacobian has w_HT = {w_exact:.5f}, so that central point is off the ray. "
            f"The declared reference uncertainties propagated through the same Jacobian give sigma(w_HT) = {sigma_w:.5f}, "
            f"a {abs(pull_sigma):.2f} sigma pull; the fixed ray stays data-compatible on this surface."
        ),
        "excluded_object": (
            f"the central reference point ({primary_id} top, {higgs_id} Higgs) as an exact point of the fixed ray "
            "pi_y = pi_lambda on the declared linear Jacobian surface"
        ),
        "quantifiers": (
            "one point, one declared linear surface, one reference codomain; the statement covers no other codomain, "
            "no nonlinear surface, and makes no claim that the fixed ray is excluded by the data"
        ),
        "instantiation": "repository objects: the declared D11 surface, the forward seed, and the compare-only exact adapter",
        "live_alternatives": [
            f"the one-scalar fixed ray itself, at {abs(pull_sigma):.2f} sigma from the {primary_id}/{higgs_id} central pair on this surface",
            f"the same ray against the {aux_id} direct-top codomain, at {abs(pull_aux):.2f} sigma",
            "the two-coordinate split Theta_D11_HT(mu_t) = (delta_y_t, delta_lambda), a chosen extension",
        ],
        "source_artifacts": {
            "d11_declared_surface": str(D11_SURFACE_JSON),
            "d11_forward_seed": str(D11_FORWARD_SEED_JSON),
            "d11_reference_exact_adapter": str(D11_EXACT_ADAPTER_JSON),
        },
        "current_fixed_ray_branch": {
            "sigma_D11_HT": sigma,
            "pi_y": sigma,
            "pi_lambda": sigma,
            "eta_HT": 0.0,
            "w_HT": 0.0,
            "law": "delta_y_t / y_t_core = -(9/16) * delta_lambda / lambda_core = sigma_D11_HT",
        },
        "exact_compare_witness": {
            "mt_pole_exact_gev": float(exact_targets["mt_pole_gev"]),
            "mH_exact_gev": float(exact_targets["mH_gev"]),
            "delta_y_t_exact": float(exact_slice["delta_y_t_mt"]),
            "delta_lambda_exact": float(exact_slice["delta_lambda_mt"]),
            "pi_y_exact": pi_y_exact,
            "pi_lambda_exact": pi_lambda_exact,
            "Sigma_HT_exact": sigma_exact,
            "eta_HT_exact": eta_exact,
            "w_HT_exact": w_exact,
        },
        "fixed_ray_point_test": {
            "compatibility_functional_formula": (
                "((mt_target - mt_pole_core_gev) / (d_mt_pole_d_y_t * y_t_core_mt)) "
                "+ (9/16) * ((mH_target - mH_core_gev) / (d_mH_d_lambda * lambda_core_mt))"
            ),
            "compatibility_functional_value": compatibility_functional,
            "fixed_ray_condition": "the central pair is an exact point of the one-scalar ray iff compatibility_functional = 0",
            "central_pair_on_fixed_ray": False,
            "central_pair_within_one_sigma_of_fixed_ray": within_one_sigma,
        },
        "uncertainty_pull": {
            "reference_sigma_mt_pole_gev": float(reference_sigma["mt_pole_gev"]),
            "reference_sigma_mH_gev": float(reference_sigma["mH_gev"]),
            "sigma_pi_y": sigma_pi_y,
            "sigma_pi_lambda": sigma_pi_lambda,
            "sigma_w_HT": sigma_w,
            "pull_sigma": pull_sigma,
            "abs_pull_sigma": abs(pull_sigma),
            "mt_on_fixed_ray_at_mH_central_gev": mt_on_ray_at_mh_central,
            "mt_fixed_ray_offset_gev": mt_ray_offset_gev,
            "mt_fixed_ray_offset_sigma": mt_ray_offset_sigma,
            "propagation": "linear, through the declared D11 Jacobian, uncorrelated reference errors",
        },
        "codomain_dependence": {
            "primary_codomain": {
                "mt_pole_summary_id": primary_id,
                "mH_summary_id": higgs_id,
                "mt_pole_description": str(reference_codomain["mt_pole_description"]),
                "w_HT": w_exact,
                "pull_sigma": pull_sigma,
            },
            "auxiliary_codomain": {
                "mt_summary_id": aux_id,
                "mt_description": str(aux_codomain["description"]),
                "mt_gev": float(aux_codomain["mt_gev"]),
                "sigma_gev": float(aux_codomain["sigma_gev"]),
                "pi_y": pi_y_aux,
                "w_HT": w_aux,
                "sigma_w_HT": sigma_w_aux,
                "pull_sigma": pull_aux,
            },
            "note": "The sign and size of w_HT depend on which top extraction codomain is paired with the Higgs value; the two codomains are distinct PDG entries.",
        },
        "smallest_supported_extension": {
            "object": "Theta_D11_HT(mu_t) = (delta_y_t, delta_lambda)",
            "equivalent_coordinates": "(Sigma_HT, eta_HT)",
            "extension_role": "chosen_extension",
            "one_extra_scalar_beyond_fixed_ray": True,
            "forced_by_data_on_declared_surface": False,
            "exact_extension_values": {
                "Sigma_HT_exact": sigma_exact,
                "eta_HT_exact": eta_exact,
                "w_HT_exact": w_exact,
            },
            "notes": [
                "The fixed ray freezes eta_HT = 0; the chosen extension adds one scalar and lands the central pair exactly.",
                "The exact center Sigma_HT_exact differs from sigma_D11_HT, so a wedge-only patch around the old center misses the central pair.",
            ],
        },
        "proof": [
            "On the one-scalar branch, pi_y = pi_lambda = sigma_D11_HT, so w_HT vanishes identically.",
            "The compare-only central pair on the same declared D11 surface has pi_y_exact != pi_lambda_exact, hence w_HT_exact != 0; the central point is off the ray.",
            "Linear propagation of the declared reference uncertainties through the same Jacobian gives sigma(w_HT); the pull |w_HT_exact| / sigma(w_HT) is below one, so the central-value offset is inside the reference band.",
            "The two-coordinate forward readout object Theta_D11_HT(mu_t) = (delta_y_t, delta_lambda), equivalently (Sigma_HT, eta_HT), is the chosen extension that lands the central pair exactly; the data on this surface do not force it.",
        ],
        "notes": [
            "This is a point statement about the one-scalar fixed ray on the declared linear surface. It does not promote any Higgs or top coordinate on the declared D10/D11 surface.",
            "The compare-only exact inverse slice stays compare-only and is used here only as a witness for the position of the central pair relative to the fixed ray.",
            "The target-anchored Higgs and top exactifiers stay compare-only validation artifacts.",
            "The artifact file name and theorem id are kept for path stability; the content is the point statement above.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D11 fixed-ray point-statement artifact.")
    parser.add_argument("--d11-surface", default=str(D11_SURFACE_JSON))
    parser.add_argument("--forward-seed", default=str(D11_FORWARD_SEED_JSON))
    parser.add_argument("--exact-adapter", default=str(D11_EXACT_ADAPTER_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    d11_surface = json.loads(Path(args.d11_surface).read_text(encoding="utf-8"))
    forward_seed = json.loads(Path(args.forward_seed).read_text(encoding="utf-8"))
    exact_adapter = json.loads(Path(args.exact_adapter).read_text(encoding="utf-8"))
    artifact = canonicalize_artifact_paths(
        build_artifact(d11_surface, forward_seed, exact_adapter)
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
