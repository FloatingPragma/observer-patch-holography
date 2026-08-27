#!/usr/bin/env python3
"""Independently replay and verify the capacity--alpha interval receipt.

The replay shares only the generic outward-rounded interval/implicit-function
engine with the producer.  It reconstructs the branch directly and evaluates
the logarithmic tangent through the factorized analytic chain rule, rather
than calling the producer or its certificate-building function.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from decimal import Decimal, InvalidOperation, localcontext
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
RER_ROOT = HERE.parents[2]
P_DERIVATION = RER_ROOT / "code" / "P_derivation"
DEFAULT_RECEIPT = HERE / "runtime" / "capacity_alpha_interval_certificate.json"
SCHEMA = "oph.capacity_alpha_interval_certificate.v1"

EXPECTED_CLASSIFICATION = {
    "epistemic_status": "rigorous_interval_certificate_for_declared_mathematical_branch",
    "mathematical_branch_differentiability": "attained_on_certified_domain",
    "physical_epoch_evolution": "undischarged",
    "physical_attachment": "requires_B1_B2_and_physical_part_of_B3",
    "OPH_evidence_status": "none_from_this_conditional_certificate",
}

EXPECTED_NONCLAIMS = {
    "B1 is not derived or discharged.",
    "B2 is not derived or discharged.",
    "The physical branch-selection and readout part of B3 is not derived or discharged.",
    "Mathematical differentiability of the declared branch does not imply physical epoch evolution.",
    "The CODATA-fed comparison coordinate is not an independent prediction.",
    "No large cosmological evolution is extrapolated beyond the certified local domain.",
}

EXPECTED_PREMISES = {
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

sys.path.insert(0, str(P_DERIVATION))
try:
    import interval_contraction_certificate as icc
finally:
    sys.path.pop(0)


class VerificationError(RuntimeError):
    """Raised on any malformed, stale, or mathematically failed receipt."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _configure_scales_independently(
    chain: Any,
    p_center: Any,
    *,
    mp_dps: int,
    iv_dps: int,
) -> None:
    point = chain.point
    mpb = point.b
    hp = mpb.num(10) ** (-14)
    u_hi, _ = point.solve_alpha_u(p_center + hp)
    u_lo, _ = point.solve_alpha_u(p_center - hp)
    dudp = abs(u_hi - u_lo) / (mpb.two * hp)
    u0, mu0 = point.solve_alpha_u(p_center)
    mz0 = point.solve_mz(u0, p_center, mu0)
    mz_hi_p = point.solve_mz(
        u0, p_center + hp, icc.mu_u_of_p(mpb, p_center + hp)
    )
    mz_lo_p = point.solve_mz(
        u0, p_center - hp, icc.mu_u_of_p(mpb, p_center - hp)
    )
    dmz_rel_dp = abs(mz_hi_p - mz_lo_p) / (mpb.two * hp) / mz0
    hu = mpb.num(10) ** (-14)
    mz_hi_u = point.solve_mz(u0 + hu, p_center, mu0)
    mz_lo_u = point.solve_mz(u0 - hu, p_center, mu0)
    dmz_rel_du = abs(mz_hi_u - mz_lo_u) / (mpb.two * hu) / mz0
    floor = mpb.num(10) ** (-(min(mp_dps, iv_dps) - 14))
    chain.scales = {
        "dudp_abs": dudp + mpb.num(10) ** (-6),
        "dmz_rel_dp": dmz_rel_dp + mpb.num(10) ** (-6),
        "dmz_rel_du": dmz_rel_du + mpb.num(10) ** (-6),
        "u_floor": floor,
        "mz_floor": floor,
    }


def independent_replay(configuration: dict[str, Any]) -> dict[str, Any]:
    """Reconstruct the certificate without importing the producer module."""

    mp_dps = int(configuration["mp_dps"])
    iv_dps = int(configuration["iv_dps"])
    su2_cutoff = int(configuration["su2_cutoff"])
    su3_cutoff = int(configuration["su3_cutoff"])
    log_half_width = str(configuration["requested_log_alpha_half_width"])

    mpb = icc.MpBackend(mp_dps)
    ivb = icc.IvBackend(iv_dps)
    point = icc.PointSolver(mpb, su2_cutoff, su3_cutoff)
    chain = icc.IntervalChain(ivb, point, su2_cutoff, su3_cutoff)
    d = chain.dual

    alpha0_iv = ivb.one / ivb.num("137.035999177")
    radius = ivb.num(log_half_width)
    alpha_box = ivb.hull(
        (alpha0_iv * ivb.exp(-radius)).a,
        (alpha0_iv * ivb.exp(radius)).b,
    )
    alpha0_mp = mpb.one / mpb.num("137.035999177")
    p0 = mpb.phi + mpb.sqrt_pi * alpha0_mp
    _configure_scales_independently(
        chain, p0, mp_dps=mp_dps, iv_dps=iv_dps
    )

    alpha = icc.Dual(alpha_box, ivb.one)
    p = d.phi + alpha * d.sqrt_pi
    mu_u = icc.mu_u_of_p(d, p)
    p_mid = icc._mid_mp(mpb, p.x)
    u_candidate, _ = point.solve_alpha_u(p_mid)
    pad = (
        chain.scales["dudp_abs"] * icc._width_mp(mpb, p.x) * mpb.num(8)
        + chain.scales["u_floor"]
    )
    u, checks = chain.verified_alpha_u(p, mu_u, u_candidate, pad)

    # Independent factorized chain rule:
    # dlnN/dlnalpha = [-6pi/(Pu)] [1+dlnu/dlnP]
    #                  [alpha/P dP/dalpha].
    du_dp = u.d / p.d
    dln_u_dln_p = p.x * du_dp / u.x
    log_n_over_pi = ivb.num(6) * ivb.pi / (p.x * u.x)
    dln_n_dln_p = -log_n_over_pi * (ivb.one + dln_u_dln_p)
    dln_p_dln_alpha = alpha_box * p.d / p.x
    tangent = dln_n_dln_p * dln_p_dln_alpha
    if ivb.contains_zero(tangent):
        raise VerificationError("independent tangent replay contains zero")
    reciprocal = ivb.one / tangent

    return {
        "alpha": icc._iv_pair(alpha_box),
        "alpha_U": icc._iv_pair(u.x),
        "tangent": icc._iv_pair(tangent),
        "reciprocal": icc._iv_pair(reciprocal),
        "R_u": checks["R_u_enclosure"],
        "h_m": checks["m_z_ift_for_input_seed"]["h_m_enclosure"],
    }


def _decimal_interval(block: dict[str, Any]) -> tuple[Decimal, Decimal]:
    try:
        lo = Decimal(str(block["lo"]))
        hi = Decimal(str(block["hi"]))
    except (KeyError, ValueError, TypeError, InvalidOperation) as exc:
        raise VerificationError("malformed decimal interval") from exc
    if not lo.is_finite() or not hi.is_finite():
        raise VerificationError("non-finite decimal interval")
    if not lo <= hi:
        raise VerificationError("interval endpoints are reversed")
    return lo, hi


def _decimal_value(value: Any, label: str) -> Decimal:
    try:
        parsed = Decimal(str(value))
    except (ValueError, TypeError, InvalidOperation) as exc:
        raise VerificationError(f"malformed decimal value: {label}") from exc
    if not parsed.is_finite():
        raise VerificationError(f"non-finite decimal value: {label}")
    return parsed


def _require_strict_sign(block: dict[str, Any], label: str) -> str:
    lo, hi = _decimal_interval(block)
    if hi < 0:
        return "strictly_negative"
    if lo > 0:
        return "strictly_positive"
    raise VerificationError(f"{label} contains zero")


def _require_contains(
    outer: dict[str, Any], inner: dict[str, Any], label: str
) -> None:
    outer_lo, outer_hi = _decimal_interval(outer)
    inner_lo, inner_hi = _decimal_interval(inner)
    if not (outer_lo <= inner_lo and inner_hi <= outer_hi):
        raise VerificationError(f"recorded {label} does not contain replay")


def verify_certificate(
    artifact: dict[str, Any], *, replay: bool = True
) -> dict[str, Any]:
    if artifact.get("schema") != SCHEMA:
        raise VerificationError("unknown certificate schema")
    if artifact.get("artifact") != "oph_capacity_alpha_interval_certificate":
        raise VerificationError("wrong artifact type")

    bindings = artifact.get("bindings", {})
    for path_key, hash_key in (
        ("producer_path", "producer_sha256"),
        ("interval_engine_path", "interval_engine_sha256"),
        ("paper_math_path", "paper_math_sha256"),
        ("verifier_path", "verifier_sha256"),
    ):
        rel = bindings.get(path_key)
        expected = bindings.get(hash_key)
        if not isinstance(rel, str) or not isinstance(expected, str):
            raise VerificationError(f"missing binding {path_key}/{hash_key}")
        path = (RER_ROOT / rel).resolve()
        if not path.is_relative_to(RER_ROOT) or not path.is_file():
            raise VerificationError(f"invalid bound path: {rel}")
        if _sha256(path) != expected:
            raise VerificationError(f"stale or mutated bound file: {rel}")

    premises = artifact.get("premises", {})
    if premises != EXPECTED_PREMISES:
        raise VerificationError("physical premise statements or statuses changed")

    classification = artifact.get("classification", {})
    if classification != EXPECTED_CLASSIFICATION:
        raise VerificationError("classification or physical attachment changed")
    nonclaims = artifact.get("nonclaims")
    if (
        not isinstance(nonclaims, list)
        or len(nonclaims) != len(EXPECTED_NONCLAIMS)
        or set(nonclaims) != EXPECTED_NONCLAIMS
    ):
        raise VerificationError("required physical nonclaims changed")

    branch = artifact.get("branch_certificate", {})
    configuration = branch.get("configuration", {})
    domain = branch.get("certified_domain", {})
    implicit = branch.get("implicit_function_certificate", {})
    if implicit.get("selected_branch_C1_on_domain") is not True:
        raise VerificationError("selected branch C1 verdict missing")
    if domain.get("positive_denominators") is not True:
        raise VerificationError("positive explicit-denominator verdict missing")

    requested_radius = _decimal_value(
        configuration.get("requested_log_alpha_half_width"),
        "requested log-alpha half-width",
    )
    inner_radius = _decimal_value(
        domain.get("guaranteed_symmetric_log_alpha_inner_radius"),
        "guaranteed symmetric inner radius",
    )
    outer_radius = _decimal_value(
        domain.get("outer_enclosure_max_abs_log_alpha_displacement"),
        "outer log-alpha enclosure radius",
    )
    if not Decimal(0) < requested_radius == inner_radius <= outer_radius:
        raise VerificationError("configuration/domain inner-radius semantics fail")
    log_lo, log_hi = _decimal_interval(domain.get("log_alpha", {}))
    try:
        iv_dps = int(configuration.get("iv_dps", 0))
        with localcontext() as ctx:
            ctx.prec = max(80, iv_dps + 16)
            center_log = (Decimal(1) / Decimal("137.035999177")).ln()
            required_lo = center_log - inner_radius
            required_hi = center_log + inner_radius
    except (ValueError, TypeError, InvalidOperation) as exc:
        raise VerificationError("invalid interval precision or log domain") from exc
    if not (log_lo <= required_lo and required_hi <= log_hi):
        raise VerificationError("certified domain does not contain inner log ball")

    derivatives = branch.get("local_derivatives", {})
    tangent = derivatives.get("d_log_N_d_log_alpha", {})
    reciprocal = derivatives.get("d_log_alpha_d_log_N", {})
    tangent_sign = _require_strict_sign(tangent, "recorded tangent")
    reciprocal_sign = _require_strict_sign(reciprocal, "recorded reciprocal")
    if tangent_sign != reciprocal_sign:
        raise VerificationError("tangent and reciprocal signs disagree")
    if tangent.get("zero_excluded") is not True:
        raise VerificationError("tangent zero-exclusion verdict missing")
    if reciprocal.get("computed_only_after_zero_exclusion") is not True:
        raise VerificationError("reciprocal ordering verdict missing")

    roots = implicit.get("alpha_U_pixel_closure", {})
    if roots.get("endpoint_signs_verified") is not True:
        raise VerificationError("alpha_U endpoint-sign verdict missing")
    if roots.get("R_u_sign_definite") is not True:
        raise VerificationError("R_u sign-definiteness missing")
    r_u_sign = _require_strict_sign(roots.get("R_u_enclosure", {}), "R_u")
    expected_r_orientation = (
        "decreasing" if r_u_sign == "strictly_negative" else "increasing"
    )
    if roots.get("orientation") != expected_r_orientation:
        raise VerificationError("alpha_U orientation contradicts R_u")
    for key in ("m_z_ift_for_R_u", "m_z_ift_for_input_seed"):
        mz = roots.get(key, {})
        if mz.get("endpoint_signs_verified") is not True:
            raise VerificationError(f"{key} endpoint-sign verdict missing")
        if mz.get("h_m_sign_definite") is not True:
            raise VerificationError(f"{key} h_m sign-definiteness missing")
        h_m_sign = _require_strict_sign(
            mz.get("h_m_enclosure", {}), f"{key} h_m"
        )
        expected_m_orientation = (
            "decreasing" if h_m_sign == "strictly_negative" else "increasing"
        )
        if mz.get("orientation") != expected_m_orientation:
            raise VerificationError(f"{key} orientation contradicts h_m")

    mean_value = branch.get("mean_value_certificate", {})
    if mean_value.get(
        "valid_for_every_pair_of_alpha_values_in_certified_domain"
    ) is not True:
        raise VerificationError("mean-value domain verdict missing")
    lower = _decimal_value(mean_value.get("abs_slope_lower"), "abs slope lower")
    upper = _decimal_value(mean_value.get("abs_slope_upper"), "abs slope upper")
    if not Decimal(0) < lower <= upper:
        raise VerificationError("invalid mean-value slope bounds")
    mean_inner = _decimal_value(
        mean_value.get("guaranteed_symmetric_log_alpha_inner_radius"),
        "mean-value inner radius",
    )
    mean_outer = _decimal_value(
        mean_value.get("outer_enclosure_max_abs_log_alpha_displacement"),
        "mean-value outer radius",
    )
    if mean_inner != inner_radius or mean_outer != outer_radius:
        raise VerificationError("mean-value/domain radius mismatch")

    tangent_lo, tangent_hi = _decimal_interval(tangent)
    if tangent_sign == "strictly_negative":
        expected_lower = tangent_hi.copy_negate()
        expected_upper = tangent_lo.copy_negate()
    else:
        expected_lower = tangent_lo
        expected_upper = tangent_hi
    if lower != expected_lower or upper != expected_upper:
        raise VerificationError("mean-value slopes do not match tangent enclosure")
    for key, label in (
        ("encloses_direct_dual_evaluation", "direct-dual tangent subblock"),
        (
            "encloses_factorized_chain_rule_evaluation",
            "factorized tangent subblock",
        ),
    ):
        _require_contains(tangent, tangent.get(key, {}), label)

    replay_block = None
    if replay:
        replay_block = independent_replay(branch["configuration"])
        _require_contains(tangent, replay_block["tangent"], "tangent")
        # Reciprocal enclosure is order-reversing; the recorded interval was
        # computed from a hull that contains the replay tangent and must
        # therefore contain its reciprocal too.
        _require_contains(reciprocal, replay_block["reciprocal"], "reciprocal")
        _require_contains(
            branch["certified_domain"]["alpha"], replay_block["alpha"], "alpha domain"
        )
        _require_strict_sign(replay_block["R_u"], "replayed R_u")
        _require_strict_sign(replay_block["h_m"], "replayed h_m")

    return {
        "verified": True,
        "schema": SCHEMA,
        "bindings_verified": True,
        "physical_premises_remain_undischarged": True,
        "independent_interval_replay": replay_block is not None,
        "replay": replay_block,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", nargs="?", type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument("--skip-replay", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        artifact = json.loads(args.receipt.read_text(encoding="utf-8"))
        result = verify_certificate(artifact, replay=not args.skip_replay)
    except (OSError, json.JSONDecodeError, VerificationError, KeyError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
