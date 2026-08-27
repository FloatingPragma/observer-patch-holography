#!/usr/bin/env python3
"""Certify the local mathematical capacity--alpha branch with intervals.

The declared finite coordinate relations are

    P(alpha) = phi + sqrt(pi) alpha,
    N(alpha) = pi exp(6 pi / (P(alpha) alpha_U(P(alpha)))),

where ``alpha_U(P)`` is the root selected by the committed finite-screen pixel
closure.  This module certifies, on a whole interval around the CODATA-fed
comparison coordinate, that the selected root is a C1 branch and that
``d log N / d log alpha`` is sign-definite.  It then encloses that derivative
and its reciprocal.  All elementary operations use outward-rounded
``mpmath.iv`` arithmetic.  The two implicit denominators are inherited from
and numerically exposed by ``interval_contraction_certificate``.

This is a theorem about a declared mathematical branch.  It supplies no
physical identification between the capacity variables, no epochwise law,
and no reason for cosmological evolution or a laboratory readout to select the
branch.  Those are the undischarged premises B1, B2, and the physical part of
B3 recorded in the receipt.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
RER_ROOT = HERE.parents[2]
P_DERIVATION = RER_ROOT / "code" / "P_derivation"
INTERVAL_ENGINE = P_DERIVATION / "interval_contraction_certificate.py"
PAPER_MATH = P_DERIVATION / "paper_math.py"
VERIFIER = HERE / "verify_capacity_alpha_interval_certificate.py"
DEFAULT_OUTPUT = HERE / "runtime" / "capacity_alpha_interval_certificate.json"

sys.path.insert(0, str(P_DERIVATION))
try:
    import interval_contraction_certificate as icc
finally:
    sys.path.pop(0)


SCHEMA = "oph.capacity_alpha_interval_certificate.v1"
ARTIFACT_DATE = "2026-08-27"
ALPHA_INVERSE_COMPARISON = "137.035999177"
DEFAULT_LOG_ALPHA_HALF_WIDTH = "0.00001"

PREMISES = {
    "B1_same_capacity": {
        "statement": (
            "The electroweak bridge coordinate N_EW is the same physical "
            "capacity N that occurs in the conditional fixed-capacity w-law."
        ),
        "status": "undischarged_physical_identification",
    },
    "B2_epochwise_bridge": {
        "statement": (
            "The finite relation N = pi exp(6 pi/(P alpha_U(P))) is a physical "
            "epoch-by-epoch law on the domain used by the cosmological model."
        ),
        "status": "undischarged_physical_epochwise_law",
    },
    "B3_physical_solver_tangent": {
        "statement": (
            "Physical evolution selects this committed root branch and "
            "co-variation convention, while a local optical-clock alpha is the "
            "same homogeneous cosmological alpha with the required time map."
        ),
        "status": "undischarged_physical_branch_selection_and_readout",
        "mathematical_subresult": (
            "C1 regularity and a sign-definite tangent of the selected declared "
            "root branch are certified by this artifact; this does not select "
            "the branch physically."
        ),
    },
}


class CapacityCertificateError(RuntimeError):
    """Raised when an interval proof obligation cannot be certified."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sign(ivb: Any, value: Any, *, label: str) -> str:
    zero = ivb.zero
    if value.b < zero.a:
        return "strictly_negative"
    if value.a > zero.b:
        return "strictly_positive"
    raise CapacityCertificateError(f"{label} interval contains zero")


def _configure_guidance_scales(
    chain: Any,
    p_center: Any,
    *,
    mp_dps: int,
    iv_dps: int,
) -> None:
    """Set bracket-inflation guidance; later interval checks carry the proof.

    These centered differences choose an efficient initial bracket only.  A
    bad estimate can make certification slower or fail closed, but it cannot
    make the endpoint signs or implicit-denominator checks pass incorrectly.
    """

    point = chain.point
    mpb = point.b
    hp = mpb.num(10) ** (-14)
    u_plus, _ = point.solve_alpha_u(p_center + hp)
    u_minus, _ = point.solve_alpha_u(p_center - hp)
    dudp = abs(u_plus - u_minus) / (mpb.two * hp)

    u_center, mu_u_center = point.solve_alpha_u(p_center)
    mz_center = point.solve_mz(u_center, p_center, mu_u_center)
    mz_plus_p = point.solve_mz(
        u_center,
        p_center + hp,
        icc.mu_u_of_p(mpb, p_center + hp),
    )
    mz_minus_p = point.solve_mz(
        u_center,
        p_center - hp,
        icc.mu_u_of_p(mpb, p_center - hp),
    )
    dmz_rel_dp = abs(mz_plus_p - mz_minus_p) / (mpb.two * hp) / mz_center

    hu = mpb.num(10) ** (-14)
    mz_plus_u = point.solve_mz(u_center + hu, p_center, mu_u_center)
    mz_minus_u = point.solve_mz(u_center - hu, p_center, mu_u_center)
    dmz_rel_du = abs(mz_plus_u - mz_minus_u) / (mpb.two * hu) / mz_center

    floor = mpb.num(10) ** (-(min(mp_dps, iv_dps) - 14))
    chain.scales = {
        "dudp_abs": dudp + mpb.num(10) ** (-6),
        "dmz_rel_dp": dmz_rel_dp + mpb.num(10) ** (-6),
        "dmz_rel_du": dmz_rel_du + mpb.num(10) ** (-6),
        "u_floor": floor,
        "mz_floor": floor,
    }


def certify_branch(
    *,
    mp_dps: int = 60,
    iv_dps: int = 60,
    su2_cutoff: int = 120,
    su3_cutoff: int = 90,
    log_alpha_half_width: str = DEFAULT_LOG_ALPHA_HALF_WIDTH,
) -> dict[str, Any]:
    """Return the outward-rounded local branch and tangent certificate."""

    if mp_dps < 28 or iv_dps < 28:
        raise ValueError("at least 28 point and interval digits are required")
    if su2_cutoff < 4 or su3_cutoff < 3:
        raise ValueError("representation cutoffs are too small for certification")

    mpb = icc.MpBackend(mp_dps)
    ivb = icc.IvBackend(iv_dps)
    point = icc.PointSolver(mpb, su2_cutoff, su3_cutoff)
    chain = icc.IntervalChain(ivb, point, su2_cutoff, su3_cutoff)
    d = chain.dual

    log_radius_iv = ivb.num(log_alpha_half_width)
    if not log_radius_iv.a > ivb.zero.b:
        raise ValueError("log-alpha half-width must be positive")

    alpha_center_iv = ivb.one / ivb.num(ALPHA_INVERSE_COMPARISON)
    alpha_lo = alpha_center_iv * ivb.exp(-log_radius_iv)
    alpha_hi = alpha_center_iv * ivb.exp(log_radius_iv)
    alpha_box = ivb.hull(alpha_lo.a, alpha_hi.b)
    if not alpha_box.a > ivb.zero.b:
        raise CapacityCertificateError("alpha domain is not positive")

    alpha_center_mp = mpb.one / mpb.num(ALPHA_INVERSE_COMPARISON)
    p_center = mpb.phi + mpb.sqrt_pi * alpha_center_mp
    _configure_guidance_scales(
        chain, p_center, mp_dps=mp_dps, iv_dps=iv_dps
    )

    alpha_dual = icc.Dual(alpha_box, ivb.one)
    p_dual = d.phi + alpha_dual * d.sqrt_pi
    mu_u_dual = icc.mu_u_of_p(d, p_dual)
    if not p_dual.x.a > ivb.zero.b:
        raise CapacityCertificateError("P domain is not positive")

    p_mid = icc._mid_mp(mpb, p_dual.x)
    u_candidate, _ = point.solve_alpha_u(p_mid)
    p_width = icc._width_mp(mpb, p_dual.x)
    pad = (
        chain.scales["dudp_abs"] * p_width * mpb.num(8)
        + chain.scales["u_floor"]
    )
    u_dual, u_checks = chain.verified_alpha_u(
        p_dual, mu_u_dual, u_candidate, pad
    )
    if not u_dual.x.a > ivb.zero.b:
        raise CapacityCertificateError("alpha_U branch is not positive")

    p_times_u = p_dual * u_dual
    if ivb.contains_zero(p_times_u.x):
        raise CapacityCertificateError("P alpha_U denominator contains zero")

    log_n = d.log(d.pi) + d.num(6) * d.pi / p_times_u
    tangent_direct = alpha_box * log_n.d

    dp_dalpha = p_dual.d
    du_dp = u_dual.d / dp_dalpha
    dln_u_dln_p = (p_dual.x / u_dual.x) * du_dp
    dln_pu_dln_p = ivb.one + dln_u_dln_p
    dln_p_dln_alpha = alpha_box * dp_dalpha / p_dual.x
    log_n_over_pi = d.num(6).x * ivb.pi / p_times_u.x
    dln_n_dln_p = -log_n_over_pi * dln_pu_dln_p
    tangent_chain_rule = dln_n_dln_p * dln_p_dln_alpha
    tangent = ivb.hull(
        min(tangent_direct.a, tangent_chain_rule.a),
        max(tangent_direct.b, tangent_chain_rule.b),
    )
    tangent_sign = _sign(ivb, tangent, label="d log N / d log alpha")
    reciprocal = ivb.one / tangent
    reciprocal_sign = _sign(
        ivb, reciprocal, label="d log alpha / d log N"
    )

    log_alpha_box = ivb.log(alpha_box)
    log_alpha_center = ivb.log(alpha_center_iv)
    centered_log_displacement = log_alpha_box - log_alpha_center
    certified_log_radius = ivb.sup_abs(centered_log_displacement)

    if tangent_sign == "strictly_negative":
        slope_abs_lower = -tangent.b
        slope_abs_upper = -tangent.a
    else:
        slope_abs_lower = tangent.a
        slope_abs_upper = tangent.b
    if not slope_abs_lower > ivb.zero.b:
        raise CapacityCertificateError("absolute tangent lower bound is not positive")

    return {
        "configuration": {
            "mp_dps": mp_dps,
            "iv_dps": iv_dps,
            "su2_cutoff": su2_cutoff,
            "su3_cutoff": su3_cutoff,
            "edge_sum_tail_bounds_included": True,
            "requested_log_alpha_half_width": log_alpha_half_width,
        },
        "comparison_coordinate": {
            "alpha_inverse": ALPHA_INVERSE_COMPARISON,
            "alpha": icc._iv_pair(alpha_center_iv),
            "P_point_display_only": icc._mp_str(mpb, p_center),
            "ancestry": (
                "CODATA-fed comparison coordinate; this is not the independent "
                "source-forward P coordinate"
            ),
        },
        "certified_domain": {
            "alpha": icc._iv_pair(alpha_box),
            "log_alpha": icc._iv_pair(log_alpha_box),
            "guaranteed_symmetric_log_alpha_inner_radius": log_alpha_half_width,
            "outer_enclosure_max_abs_log_alpha_displacement": (
                icc._iv_pair(certified_log_radius)["hi"]
            ),
            "P": icc._iv_pair(p_dual.x),
            "alpha_U": icc._iv_pair(u_dual.x),
            "P_times_alpha_U": icc._iv_pair(p_times_u.x),
            "positive_denominators": True,
        },
        "implicit_function_certificate": {
            "alpha_U_pixel_closure": u_checks,
            "selected_branch_C1_on_domain": True,
            "selection_scope": (
                "the root selected by the committed scan-and-bisect solver; "
                "uniqueness is proved inside the verified bracket, while global "
                "uniqueness over the full scan window is not claimed"
            ),
        },
        "local_derivatives": {
            "d_log_alpha_U_d_log_P": icc._iv_pair(dln_u_dln_p),
            "d_log_P_alpha_U_d_log_P": icc._iv_pair(dln_pu_dln_p),
            "log_N_over_pi": icc._iv_pair(log_n_over_pi),
            "d_log_N_d_log_P": icc._iv_pair(dln_n_dln_p),
            "d_log_P_d_log_alpha": icc._iv_pair(dln_p_dln_alpha),
            "d_log_N_d_log_alpha": {
                **icc._iv_pair(tangent),
                "sign": tangent_sign,
                "zero_excluded": True,
                "encloses_direct_dual_evaluation": icc._iv_pair(tangent_direct),
                "encloses_factorized_chain_rule_evaluation": icc._iv_pair(
                    tangent_chain_rule
                ),
            },
            "d_log_alpha_d_log_N": {
                **icc._iv_pair(reciprocal),
                "sign": reciprocal_sign,
                "computed_only_after_zero_exclusion": True,
            },
        },
        "mean_value_certificate": {
            "valid_for_every_pair_of_alpha_values_in_certified_domain": True,
            "guaranteed_symmetric_log_alpha_inner_radius": log_alpha_half_width,
            "outer_enclosure_max_abs_log_alpha_displacement": (
                icc._iv_pair(certified_log_radius)["hi"]
            ),
            "abs_slope_lower": icc._iv_pair(slope_abs_lower)["lo"],
            "abs_slope_upper": icc._iv_pair(slope_abs_upper)["hi"],
            "forward_bound": (
                "m |Delta log alpha| <= |Delta log N| <= "
                "L |Delta log alpha|, with m=abs_slope_lower and "
                "L=abs_slope_upper"
            ),
            "inverse_bound": (
                "|Delta log alpha| <= (1/m) |Delta log N| for points on "
                "the certified branch image"
            ),
            "argument": (
                "The implicit-function denominators are sign-definite on the "
                "whole box, so the selected branch is C1.  The displayed "
                "logarithmic derivative is sign-definite there; the ordinary "
                "mean-value theorem in log alpha gives the two-sided bound and "
                "local invertibility."
            ),
        },
    }


def build_certificate(**kwargs: Any) -> dict[str, Any]:
    branch = certify_branch(**kwargs)
    return {
        "artifact": "oph_capacity_alpha_interval_certificate",
        "schema": SCHEMA,
        "date": ARTIFACT_DATE,
        "classification": {
            "epistemic_status": "rigorous_interval_certificate_for_declared_mathematical_branch",
            "mathematical_branch_differentiability": "attained_on_certified_domain",
            "physical_epoch_evolution": "undischarged",
            "physical_attachment": "requires_B1_B2_and_physical_part_of_B3",
            "OPH_evidence_status": "none_from_this_conditional_certificate",
        },
        "premises": PREMISES,
        "bindings": {
            "producer_path": str(Path(__file__).resolve().relative_to(RER_ROOT)),
            "producer_sha256": sha256(Path(__file__).resolve()),
            "interval_engine_path": str(INTERVAL_ENGINE.relative_to(RER_ROOT)),
            "interval_engine_sha256": sha256(INTERVAL_ENGINE),
            "paper_math_path": str(PAPER_MATH.relative_to(RER_ROOT)),
            "paper_math_sha256": sha256(PAPER_MATH),
            "verifier_path": str(VERIFIER.relative_to(RER_ROOT)),
            "verifier_sha256": sha256(VERIFIER) if VERIFIER.exists() else None,
        },
        "branch_certificate": branch,
        "nonclaims": [
            "B1 is not derived or discharged.",
            "B2 is not derived or discharged.",
            "The physical branch-selection and readout part of B3 is not derived or discharged.",
            "Mathematical differentiability of the declared branch does not imply physical epoch evolution.",
            "The CODATA-fed comparison coordinate is not an independent prediction.",
            "No large cosmological evolution is extrapolated beyond the certified local domain.",
        ],
    }


def write_certificate(path: Path, certificate: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mp-dps", type=int, default=60)
    parser.add_argument("--iv-dps", type=int, default=60)
    parser.add_argument("--su2-cutoff", type=int, default=120)
    parser.add_argument("--su3-cutoff", type=int, default=90)
    parser.add_argument(
        "--log-alpha-half-width", default=DEFAULT_LOG_ALPHA_HALF_WIDTH
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    certificate = build_certificate(
        mp_dps=args.mp_dps,
        iv_dps=args.iv_dps,
        su2_cutoff=args.su2_cutoff,
        su3_cutoff=args.su3_cutoff,
        log_alpha_half_width=args.log_alpha_half_width,
    )
    write_certificate(args.output, certificate)
    tangent = certificate["branch_certificate"]["local_derivatives"][
        "d_log_N_d_log_alpha"
    ]
    reciprocal = certificate["branch_certificate"]["local_derivatives"][
        "d_log_alpha_d_log_N"
    ]
    print(f"d log N / d log alpha in [{tangent['lo']}, {tangent['hi']}]")
    print(f"d log alpha / d log N in [{reciprocal['lo']}, {reciprocal['hi']}]")
    print(f"saved: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
