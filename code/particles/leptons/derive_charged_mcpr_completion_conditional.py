#!/usr/bin/env python3
"""Evaluate the MCPR charged-response completion on the source-audit branch.

MCPR (minimal complete public response) is the submitted completion package
that closes the charged-lepton response by stipulating a typed, incidence
complete response architecture on an oriented icosahedral flag: eight finite
registers, eight primitive path classes, normalized rank-one trace weights, a
public-block amplitude normalization, an oriented process phase, and a Z6
determinant character.  Conditional on that architecture, the response map has
a unique globally attracting fixed point and emits one dimensionless charged
triple with zero runtime charged reference input.

This module recomputes the entire chain from the source-audit pixel
certificate and freezes the result as a conditional completed-model candidate.
The 2026-07-13 symmetry-breaking audit classifies the architecture itself as
stipulated: the register dimensions, path table, signs, amplitude choice,
phase offset, and determinant exponent are declared inputs of the model, so
the output is a declared-model coordinate on the OPH+MCPR theory.  The
matching frontier objects that a source-only promotion requires (W5-ORB,
A5-FAM, DET-CAN, the interacting kernel packet, the operational scale, and
prospective provenance) are recorded as open gates.

No electron, muon, tau, square-root-mass, pole-mass, or charged-Yukawa
numerical reference is read by this module.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import json
from pathlib import Path
from typing import Any

import mpmath as mp

WORKING_DPS = 120


def _scoped_dps(func):
    """Run func at WORKING_DPS and restore the global precision on exit.

    The module never mutates the global mpmath precision at import; every
    entry point that computes with mpmath carries its own scope so parallel
    modules keep their own default precision.
    """

    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        with mp.workdps(WORKING_DPS):
            return func(*args, **kwargs)

    return wrapper

HERE = Path(__file__).resolve()
CODE_ROOT = HERE.parents[2]
SOURCE_CERTIFICATE = (
    CODE_ROOT
    / "particles"
    / "hierarchy"
    / "certificates"
    / "R_P_source_audit_pixel_certificate.json"
)
DEFAULT_OUT = (
    CODE_ROOT
    / "particles"
    / "runs"
    / "leptons"
    / "charged_mcpr_completion_conditional.json"
)

# Display constant for the optional MeV column.  The clock lane classifies
# this decimal as a calibration checksum; it never enters the response map.
DISPLAY_E_STAR_GEV = "1.220890e19"

EXPECTED_REGISTER_DIMENSIONS = [50, 31, 10, 512, 77, 21, 27, 5]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _text(value: mp.mpf | mp.mpc, digits: int = 60) -> str:
    return mp.nstr(value, digits)


def icosahedral_incidence_solve() -> dict[str, Any]:
    """Solve 5V=2E, 3F=2E, V-E+F=2 over positive integers."""

    solutions: list[tuple[int, int, int]] = []
    for vertices in range(1, 100):
        for edges in range(1, 300):
            if 5 * vertices != 2 * edges:
                continue
            if (2 * edges) % 3:
                continue
            faces = 2 * edges // 3
            if vertices - edges + faces == 2:
                solutions.append((vertices, edges, faces))
    return {
        "solutions": solutions,
        "unique": solutions == [(12, 30, 20)],
        "V": 12,
        "E": 30,
        "F": 20,
        "euler": 2,
        "valence": 5,
        "face_size": 3,
    }


def register_dimensions(geom: dict[str, Any]) -> dict[str, int]:
    """Stipulated incidence registers of the MCPR flag architecture.

    Every dimension is a finite incidence count relative to one oriented face
    with one marked corner.  The audit classifies the selection of exactly
    these eight registers as a declared model input.
    """

    v, e, f = geom["V"], geom["E"], geom["F"]
    euler, valence = geom["euler"], geom["valence"]
    n_g = 3
    return {
        "kappa_0_edges_plus_faces": e + f,
        "kappa_1_edges_plus_identity": e + 1,
        "kappa_2_transverse_vertices_plus_identity": v - n_g + 1,
        "chi_0_subsets_of_transverse_vertices": 2 ** (v - n_g),
        "chi_1_collar_labels_times_nonanchor_vertices": (valence + euler) * (v - 1),
        "zeta_0_faces_plus_identity": f + 1,
        "zeta_1_edges_outside_marked_face": e - n_g,
        "zeta_2_anchor_star": valence,
    }


@_scoped_dps
def response_fixed_point(
    p_value: mp.mpf, alpha_u: mp.mpf, dims: dict[str, int]
) -> dict[str, Any]:
    """Solve the three-channel MCPR response map exactly."""

    phi = (1 + mp.sqrt(5)) / 2
    alpha_p = (p_value - phi) / mp.sqrt(mp.pi)

    s_kappa = alpha_u / dims["kappa_0_edges_plus_faces"]
    q_kappa = alpha_u / dims["kappa_1_edges_plus_identity"] + alpha_u**2 / (
        dims["kappa_1_edges_plus_identity"]
        * dims["kappa_2_transverse_vertices_plus_identity"]
    )
    s_chi = alpha_u**2 / dims["chi_0_subsets_of_transverse_vertices"]
    q_chi = alpha_u / dims["chi_1_collar_labels_times_nonanchor_vertices"]
    s_zeta = alpha_p**2 / dims["zeta_0_faces_plus_identity"]
    q_zeta = alpha_u / dims["zeta_1_edges_outside_marked_face"] + alpha_u**2 / (
        dims["zeta_1_edges_outside_marked_face"] * dims["zeta_2_anchor_star"]
    )

    kappa = s_kappa / (1 + q_kappa)
    chi = s_chi / (1 - q_chi)
    zeta = s_zeta / (1 - q_zeta)

    image = (
        s_kappa - q_kappa * kappa,
        s_chi + q_chi * chi,
        s_zeta + q_zeta * zeta,
    )
    residual = [image[i] - (kappa, chi, zeta)[i] for i in range(3)]
    contraction = max(abs(q_kappa), abs(q_chi), abs(q_zeta))

    return {
        "alpha_p": alpha_p,
        "sources": (s_kappa, s_chi, s_zeta),
        "feedback": (q_kappa, q_chi, q_zeta),
        "kappa": kappa,
        "chi": chi,
        "zeta": zeta,
        "residual": residual,
        "contraction": contraction,
    }


@_scoped_dps
def charged_triple(
    p_value: mp.mpf,
    alpha_u: mp.mpf,
    fixed_point: dict[str, Any],
    geom: dict[str, Any],
) -> dict[str, Any]:
    """Map the fixed point to the dimensionless charged triple."""

    n_c = 3
    n_g = 3
    beta_ew = n_c + 1
    kappa = fixed_point["kappa"]
    chi = fixed_point["chi"]
    zeta = fixed_point["zeta"]

    delta_0 = mp.mpf(beta_ew) / (2 * n_c * n_g)
    delta = delta_0 + zeta
    rho = mp.sqrt(2) * mp.e ** (-chi)
    roots = [1 + rho * mp.cos(delta + 2 * mp.pi * j / 3) for j in range(3)]
    reversed_roots = [1 + rho * mp.cos(-delta + 2 * mp.pi * j / 3) for j in range(3)]
    squares = sorted(root * root for root in roots)
    reversed_squares = sorted(root * root for root in reversed_roots)
    shape_geomean = mp.fprod(squares) ** (mp.mpf(1) / 3)

    determinant_exponent = geom["V"] + geom["euler"]
    det_m_over_v_cubed = mp.e ** (3 * kappa) / (2 * mp.power(6, determinant_exponent))
    geometric_mean_over_v = det_m_over_v_cubed ** (mp.mpf(1) / 3)
    masses_over_v = [geometric_mean_over_v * square / shape_geomean for square in squares]

    v_over_e_star = p_value ** (-mp.mpf(1) / 2) * mp.e ** (
        -2 * mp.pi / (beta_ew * alpha_u)
    )
    masses_over_e_star = [value * v_over_e_star for value in masses_over_v]

    sqrt_mass_invariant_direct = mp.fsum(masses_over_v) / (
        mp.fsum(mp.sqrt(value) for value in masses_over_v) ** 2
    )
    sqrt_mass_invariant_formula = (1 + mp.e ** (-2 * chi)) / 3

    return {
        "delta_0": delta_0,
        "delta": delta,
        "rho": rho,
        "roots": roots,
        "squares": squares,
        "reversed_squares": reversed_squares,
        "shape_geomean": shape_geomean,
        "determinant_exponent": determinant_exponent,
        "det_m_over_v_cubed": det_m_over_v_cubed,
        "masses_over_v": masses_over_v,
        "v_over_e_star": v_over_e_star,
        "masses_over_e_star": masses_over_e_star,
        "ratios": (
            masses_over_v[1] / masses_over_v[0],
            masses_over_v[2] / masses_over_v[0],
            masses_over_v[2] / masses_over_v[1],
        ),
        "sqrt_mass_invariant_direct": sqrt_mass_invariant_direct,
        "sqrt_mass_invariant_formula": sqrt_mass_invariant_formula,
    }


def audit_boundary() -> dict[str, Any]:
    """Claim boundary from the 2026-07-13 symmetry-breaking audit."""

    return {
        "verdict": "SOURCE_ONLY_PHYSICAL_PREDICTION_OPEN",
        "architecture_provenance": (
            "The eight register dimensions, path table, signs, source "
            "monomials, amplitude prescription, phase offset, determinant "
            "exponent, and physical readout are declared model inputs. "
            "Arithmetic consistency checks establish contraction, fixed-point "
            "identities, positivity, simplicity, and determinant arithmetic "
            "after those choices are supplied; they do not establish that the "
            "source selects the architecture."
        ),
        "amplitude_boundary": (
            "The public-block amplitude 1/sqrt(2) corresponds to the tracial "
            "weight t=1/2 on C (+) M_2(C). Traciality and C3 invariance admit "
            "the full one-parameter family; the Hilbert-space normalized "
            "trace has t=1/3 and gives amplitude 1/2. A physical "
            "state-selection and response-isometry theorem is open."
        ),
        "phase_boundary": (
            "The identity (N_c+1)/(2 N_c N_g) = 2/9 does not select the "
            "oriented phase. A holonomy reading requires link variables, a "
            "loop, and a nonlocal readout map; a common one-link phase on a "
            "triangular C3 cycle exponentiates to three times itself."
        ),
        "determinant_boundary": (
            "A product of fourteen normalized rank-one traces is neither a "
            "Z6 group character nor a Yukawa determinant. The canonical "
            "determinant norm carries the kinetic factor "
            "det M_L = (v/sqrt2)^3 |det Y_e| / sqrt(det Z_L det Z_E). The "
            "6^-14 weight is a candidate positive path weight."
        ),
        "maxent_shape_silence": (
            "A unique G-invariant MaxEnt state has zero expectation on every "
            "nontrivial irreducible module, and the current homogeneous "
            "twelve-port branch has the unique uniform minimizer. Charged "
            "family shape therefore requires a symmetry-breaking selected "
            "orbit, an extremal conditioned branch, or a source non-singlet."
        ),
        "completion_route": {
            "W5_ORB": {
                "status": "open",
                "statement": (
                    "an OPH-derived A5-invariant effective action on W5 whose "
                    "unique minimizing orbit is simple-spectrum with "
                    "Delta_Q > 0; first decisive quantity is the W5 Hessian "
                    "scalar h_5 at the homogeneous state"
                ),
            },
            "A5_FAM": {
                "status": "open",
                "statement": (
                    "physical chiral family fibers carrying 3 or 3' of A5, "
                    "with the multiplicity-one attachment "
                    "Sym^2_0(3) ~= W5 and a source-fixed sign"
                ),
            },
            "DET_CAN": {
                "status": "open",
                "statement": (
                    "a normed, sector-isolated determinant-line descent with "
                    "the kinetic determinant factor included"
                ),
            },
            "QFT_POLE": {
                "status": "open",
                "statement": (
                    "a BRST-complete interacting charged kernel with "
                    "renormalization fixed without charged targets and "
                    "singularity invariants equal to the source invariants"
                ),
            },
            "SCALE": {
                "status": "open",
                "statement": (
                    "a source-derived operational clock ratio for MeV output"
                ),
            },
            "PROV": {
                "status": "open",
                "statement": (
                    "canonical taint classes, transitive hashes, and a "
                    "prospective freeze before held-out evaluation"
                ),
            },
        },
    }


@_scoped_dps
def build_artifact(source: dict[str, Any], source_meta: dict[str, Any]) -> dict[str, Any]:
    p_value = mp.mpf(str(source["P_cand"]))
    alpha_u = mp.mpf(str(source["alpha_U_P_cand"]))

    geom = icosahedral_incidence_solve()
    dims = register_dimensions(geom)
    fixed_point = response_fixed_point(p_value, alpha_u, dims)
    triple = charged_triple(p_value, alpha_u, fixed_point, geom)

    min_root_gap = min(
        abs(triple["roots"][i] - triple["roots"][j])
        for i in range(3)
        for j in range(i)
    )
    product_v = mp.fprod(triple["masses_over_v"])

    checks = {
        "icosahedral_incidence_unique": geom["unique"],
        "register_dimensions_are_50_31_10_512_77_21_27_5": list(dims.values())
        == EXPECTED_REGISTER_DIMENSIONS,
        "strict_contraction": fixed_point["contraction"] < 1,
        "fixed_point_residual_below_1e_100": max(
            abs(x) for x in fixed_point["residual"]
        )
        < mp.mpf("1e-100"),
        "positive_simple_regular_C3_roots": min(triple["roots"]) > 0
        and min_root_gap > mp.mpf("1e-30"),
        "orientation_reversal_isospectral": max(
            abs(a - b)
            for a, b in zip(triple["squares"], triple["reversed_squares"], strict=True)
        )
        < mp.mpf("1e-100"),
        "determinant_identity_in_v_units": abs(
            product_v / triple["det_m_over_v_cubed"] - 1
        )
        < mp.mpf("1e-100"),
        "sqrt_mass_invariant_identity": abs(
            triple["sqrt_mass_invariant_direct"]
            - triple["sqrt_mass_invariant_formula"]
        )
        < mp.mpf("1e-100"),
        "source_certificate_v_over_Estar_consistency": abs(
            triple["v_over_e_star"] / mp.mpf(str(source["v_over_E_star_P_cand"])) - 1
        )
        < mp.mpf("2e-15"),
    }
    checks_pass = all(bool(v) for v in checks.values())

    e_star = mp.mpf(DISPLAY_E_STAR_GEV)
    display = {
        "display_scale_status": (
            "downstream unit display only; the decimal is classified by the "
            "clock lane as a calibration checksum and enters no solve path"
        ),
        "E_star_GeV": DISPLAY_E_STAR_GEV,
        "masses_MeV": [
            _text(1000 * value * e_star) for value in triple["masses_over_e_star"]
        ],
    }

    return {
        "artifact": "oph_charged_mcpr_completion_conditional",
        "schema_version": 1,
        "status": "CONDITIONAL_COMPLETED_MODEL_CANDIDATE_NOT_SOURCE_ONLY",
        "row_class": "conditional_on_declared_mcpr_response_architecture",
        "source_only": False,
        "charged_reference_data_consumed": False,
        "public_charged_mass_promotion_allowed": False,
        "claim_class": {
            "current_OPH5_alone_unique": False,
            "conditional_on_declared_architecture_unique": True,
            "new_continuous_or_fitted_parameters": 0,
            "historically_blind_prospective_prediction": False,
            "runtime_charged_reference_consumed": False,
        },
        "source_input": {
            "path": str(source_meta.get("path", "")),
            "sha256": source_meta.get("sha256", ""),
            "artifact": source.get("artifact"),
            "certificate_status": source.get("status"),
            "P": _text(p_value, 45),
            "alpha_U": _text(alpha_u, 20),
        },
        "incidence": {k: geom[k] for k in ("V", "E", "F", "euler", "valence", "face_size", "unique")},
        "response_register_dimensions": dims,
        "response_map": {
            "alpha_P": _text(fixed_point["alpha_p"]),
            "sources": [_text(x) for x in fixed_point["sources"]],
            "feedback": [_text(x) for x in fixed_point["feedback"]],
            "contraction_constant": _text(fixed_point["contraction"]),
            "fixed_point": {
                "kappa": _text(fixed_point["kappa"]),
                "chi": _text(fixed_point["chi"]),
                "zeta": _text(fixed_point["zeta"]),
            },
        },
        "regular_C3_shape": {
            "delta_0": _text(triple["delta_0"]),
            "delta": _text(triple["delta"]),
            "rho": _text(triple["rho"]),
            "roots_unsorted": [_text(x) for x in triple["roots"]],
            "root_squares_sorted": [_text(x) for x in triple["squares"]],
            "sqrt_mass_invariant": _text(triple["sqrt_mass_invariant_direct"]),
        },
        "determinant_character": {
            "primitive_Z6_event_count": triple["determinant_exponent"],
            "formula_det_M_over_v_cubed": "exp(3*kappa)/(2*6^14)",
            "det_M_over_v_cubed": _text(triple["det_m_over_v_cubed"]),
            "audit_classification": "candidate positive path weight",
        },
        "conditional_prediction": {
            "v_over_E_star": _text(triple["v_over_e_star"]),
            "masses_over_v": [_text(x) for x in triple["masses_over_v"]],
            "masses_over_E_star": [_text(x) for x in triple["masses_over_e_star"]],
            "ratios": {
                "m_mu_over_m_e": _text(triple["ratios"][0]),
                "m_tau_over_m_e": _text(triple["ratios"][1]),
                "m_tau_over_m_mu": _text(triple["ratios"][2]),
            },
        },
        "optional_scale_display": display,
        "audit_boundary": audit_boundary(),
        "checks": {k: bool(v) for k, v in checks.items()},
        "checks_pass": checks_pass,
        "no_target_dependency_statement": (
            "No electron, muon, tau, square-root-mass-invariant, pole-mass, "
            "or charged-Yukawa numerical reference is read by this module."
        ),
    }


def build() -> dict[str, Any]:
    source = json.loads(SOURCE_CERTIFICATE.read_text(encoding="utf-8"))
    meta = {
        "path": "hierarchy/certificates/R_P_source_audit_pixel_certificate.json",
        "sha256": sha256_file(SOURCE_CERTIFICATE),
    }
    return build_artifact(source, meta)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": artifact["status"],
                "checks_pass": artifact["checks_pass"],
                "ratios": artifact["conditional_prediction"]["ratios"],
                "masses_MeV_display": artifact["optional_scale_display"]["masses_MeV"],
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
