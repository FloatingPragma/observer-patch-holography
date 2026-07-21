#!/usr/bin/env python3
"""Global at-most-one fixed-point certificate for the OPH P/alpha closure map.

The stage-2 contraction certificate (``interval_contraction_certificate.py``)
proves existence and uniqueness of the closure-map fixed point on a small
alpha interval around each fixed point. This module upgrades the uniqueness
half from interval-local to global on the declared physical domain:

- the declared domain is the solver scan window of
  ``paper_math.PaperMathContext.solve_closure``: alpha in [0.005, 0.01]
  (equivalently alpha_inv in [100, 200]); both certified fixed points
  (source mode near alpha_inv = 136.9948 and gauge-width mode near
  alpha_inv = 137.0356, each on its own map) lie inside it;
- the domain is covered by an adaptive interval subdivision; on each piece
  the map derivative g'(alpha) is enclosed by the same forward-mode interval
  AD chain the stage-2 certificate uses (verified implicit brackets,
  sign-definite residual derivatives, SU(2)/SU(3) edge-sum tail majorants
  folded in);
- a piece certifies by, in order of strength: contraction (sup |g'| < 1),
  sign-definite residual derivative (g' - 1 < 0 or g' - 1 > 0, so the
  residual g(alpha) - alpha is strictly monotone on the piece), or a
  residual enclosure bounded away from zero (no fixed point in the piece);
- pieces that certify none of these are bisected up to a depth cap and a
  deterministic work budget; anything left is reported as the exceptional
  set and the conclusion degrades to partial;
- the per-piece verdicts are synthesized into a global at-most-one
  statement: sup |g'| < 1 on every piece gives |g'| < 1 on the convex
  domain and hence at most one fixed point by the mean value theorem;
  otherwise strictly monotone residual runs plus root-free gaps bound the
  root count, with certified endpoint residual signs closing the count.

Claim boundary: the certificate covers the declared numerical map at the
declared representation cutoffs (the edge-sum tail majorants extend every
enclosure to the infinite-cutoff sums); it certifies at-most-one on the
declared domain, not any relation to the measured fine-structure constant,
and the stage-3 landing verdict (closure row CL-1) is unchanged.

mpmath contexts are private instances; the global ``mpmath.mp`` and
``mpmath.iv`` precision settings are never touched. The construction is
deterministic: the subdivision order, the work budget, and every verdict are
functions of the declared parameters only, and no wall-clock quantity is
recorded in the artifact.
"""

from __future__ import annotations

import argparse
import json
import time
from decimal import Decimal
from pathlib import Path
from typing import Any, Callable

from mpmath.libmp.libmpf import ComplexResult

import interval_contraction_certificate as icc


ARTIFACT_DATE = "2026-07-14"
DEFAULT_OUT = (
    Path(__file__).resolve().parent
    / "runtime"
    / f"p_global_uniqueness_certificate_{ARTIFACT_DATE}.json"
)

DOMAIN_LO = "0.005"
DOMAIN_HI = "0.01"
DOMAIN_DECLARATION = (
    "declared solver scan window of paper_math.PaperMathContext.solve_closure "
    "(alpha_min = 0.005, alpha_max = 0.01, i.e. alpha_inv in [100, 200]); "
    "adopted as the declared physical domain of the closure map"
)

# Interval evaluation on a wide box can fail structurally (bracket inflation
# exhausted) or arithmetically (log of an interval touching zero, division by
# an interval containing zero). Both are handled as an undecided piece and
# resolved by bisection. ComplexResult subclasses ValueError; it is named
# explicitly for clarity.
EVAL_ERRORS = (icc.CertificateError, ComplexResult, ValueError, ArithmeticError)

VERDICT_CONTRACTION = "contraction"
VERDICT_DECREASING = "residual_decreasing"
VERDICT_INCREASING = "residual_increasing"
VERDICT_ROOT_FREE_POS = "root_free_residual_positive"
VERDICT_ROOT_FREE_NEG = "root_free_residual_negative"
VERDICT_UNDECIDED = "undecided"
VERDICT_EVAL_FAILED = "evaluation_failed"

DECISIVE_VERDICTS = {
    VERDICT_CONTRACTION,
    VERDICT_DECREASING,
    VERDICT_INCREASING,
    VERDICT_ROOT_FREE_POS,
    VERDICT_ROOT_FREE_NEG,
}
DECREASING_FAMILY = {VERDICT_CONTRACTION, VERDICT_DECREASING}
ROOT_FREE = {VERDICT_ROOT_FREE_POS, VERDICT_ROOT_FREE_NEG}


def _endpoint_mp(mpb: icc.MpBackend, endpoint: Any, side: int) -> Any:
    """Convert an mpmath.iv endpoint (a thin interval) to an mp number."""
    return mpb.ctx.mpf(endpoint._mpi_[side])


def prepare_chain(
    mode: str,
    mp_dps: int,
    iv_dps: int,
    su2_cutoff: int,
    su3_cutoff: int,
) -> tuple[icc.MpBackend, icc.PointSolver, icc.IvBackend, icc.IntervalChain, Any]:
    """Build point solver and interval chain with the stage-2 scale seeding."""
    mpb = icc.MpBackend(mp_dps)
    point = icc.PointSolver(mpb, su2_cutoff, su3_cutoff)
    ivb = icc.IvBackend(iv_dps)
    chain = icc.IntervalChain(ivb, point, su2_cutoff, su3_cutoff)

    alpha_star = point.fixed_point(mode)
    inner_star = point.inner(alpha_star, mode)

    # Finite-difference scale estimates (bracket-inflation guidance only; the
    # verified sign checks carry the rigor, exactly as in the stage-2 module).
    hp = mpb.num(10) ** (-14)
    p_star = inner_star["p"]
    u_plus, _ = point.solve_alpha_u(p_star + hp)
    u_minus, _ = point.solve_alpha_u(p_star - hp)
    dudp = abs(u_plus - u_minus) / (2 * hp)
    mu_u_star = inner_star["mu_u"]
    u_star = inner_star["alpha_u"]
    mz_star = inner_star["mz"]
    mzp = point.solve_mz(u_star, p_star + hp, icc.mu_u_of_p(mpb, p_star + hp))
    mzm = point.solve_mz(u_star, p_star - hp, icc.mu_u_of_p(mpb, p_star - hp))
    dmz_rel_dp = abs(mzp - mzm) / (2 * hp) / mz_star
    hu = mpb.num(10) ** (-14)
    mzpu = point.solve_mz(u_star + hu, p_star, mu_u_star)
    mzmu = point.solve_mz(u_star - hu, p_star, mu_u_star)
    dmz_rel_du = abs(mzpu - mzmu) / (2 * hu) / mz_star
    floor = mpb.num(10) ** (-(min(mp_dps, iv_dps) - 14))
    chain.scales = {
        "dudp_abs": dudp + mpb.num(10) ** (-6),
        "dmz_rel_dp": dmz_rel_dp + mpb.num(10) ** (-6),
        "dmz_rel_du": dmz_rel_du + mpb.num(10) ** (-6),
        "u_floor": floor,
        "mz_floor": floor,
    }
    return mpb, point, ivb, chain, alpha_star


def classify_piece(
    chain: icc.IntervalChain,
    ivb: icc.IvBackend,
    mpb: icc.MpBackend,
    lo_mp: Any,
    hi_mp: Any,
    mode: str,
) -> dict[str, Any]:
    """Evaluate one subdivision piece and return its verdict record."""
    box = ivb.hull(ivb.thin(lo_mp).a, ivb.thin(hi_mp).b)
    record: dict[str, Any] = {"lo": lo_mp, "hi": hi_mp}
    try:
        g_box, _ = chain.g_dual(box, mode)
    except EVAL_ERRORS as exc:
        record["verdict"] = VERDICT_EVAL_FAILED
        record["reason"] = f"{type(exc).__name__}: {exc}"
        return record

    gp = g_box.d
    record["gprime"] = icc._iv_pair(gp)
    lipschitz = ivb.sup_abs(gp)
    if lipschitz < ivb.one.a:
        record["verdict"] = VERDICT_CONTRACTION
        record["sup_abs_gprime"] = _endpoint_mp(mpb, lipschitz, 1)
        return record
    if gp.b < ivb.one.a:
        record["verdict"] = VERDICT_DECREASING
        return record
    if gp.a > ivb.one.b:
        record["verdict"] = VERDICT_INCREASING
        return record

    residual = g_box.x - box
    record["residual"] = icc._iv_pair(residual)
    if residual.a > ivb.zero.b:
        record["verdict"] = VERDICT_ROOT_FREE_POS
        return record
    if residual.b < ivb.zero.a:
        record["verdict"] = VERDICT_ROOT_FREE_NEG
        return record

    record["verdict"] = VERDICT_UNDECIDED
    record["reason"] = "derivative enclosure straddles the critical value and the residual enclosure contains zero"
    return record


def synthesize(
    pieces: list[dict[str, Any]],
    endpoint_sign: Callable[[Any], str | None],
) -> dict[str, Any]:
    """Combine sorted per-piece verdicts into a global at-most-one statement.

    ``pieces`` must tile the domain in ascending order. ``endpoint_sign`` maps
    an alpha value (mp number) to a certified residual sign ("+", "-") or
    None when the sign cannot be certified; it is only consulted when the
    verdict mix requires run-endpoint bookkeeping.
    """
    verdict_counts: dict[str, int] = {}
    for piece in pieces:
        verdict_counts[piece["verdict"]] = verdict_counts.get(piece["verdict"], 0) + 1

    unresolved = [p for p in pieces if p["verdict"] not in DECISIVE_VERDICTS]
    if unresolved:
        return {
            "at_most_one_on_domain": False,
            "status": "partial",
            "argument": (
                "the subdivision left pieces with no certified verdict; roots "
                "cannot be excluded or counted on the exceptional set"
            ),
            "verdict_counts": verdict_counts,
        }

    verdicts = {p["verdict"] for p in pieces}
    if verdicts <= {VERDICT_CONTRACTION}:
        return {
            "at_most_one_on_domain": True,
            "status": "certified",
            "argument": (
                "sup |g'| < 1 certified on every piece of a finite tiling of the "
                "domain, so |g'| < 1 holds on the convex domain; two distinct "
                "fixed points x < y would give |x - y| = |g(x) - g(y)| = "
                "|g'(xi)| |x - y| < |x - y| by the mean value theorem (g is C^1 "
                "on each piece: the implicit-function denominators are "
                "interval-verified sign-definite); at most one fixed point on "
                "the domain follows"
            ),
            "verdict_counts": verdict_counts,
        }
    if verdicts <= DECREASING_FAMILY:
        return {
            "at_most_one_on_domain": True,
            "status": "certified",
            "argument": (
                "g' - 1 < 0 certified on every piece, so the residual "
                "g(alpha) - alpha is strictly decreasing on the whole domain "
                "and has at most one zero"
            ),
            "verdict_counts": verdict_counts,
        }
    if verdicts <= {VERDICT_INCREASING}:
        return {
            "at_most_one_on_domain": True,
            "status": "certified",
            "argument": (
                "g' - 1 > 0 certified on every piece, so the residual is "
                "strictly increasing on the whole domain and has at most one zero"
            ),
            "verdict_counts": verdict_counts,
        }

    # Mixed verdicts: roots live only in monotone runs (root-free pieces
    # exclude them, endpoints included); each maximal monotone run is strictly
    # monotone by continuity, so it carries at most one root, and certified
    # endpoint residual signs exclude runs whose signs forbid a crossing.
    runs: list[dict[str, Any]] = []
    for piece in pieces:
        if piece["verdict"] in ROOT_FREE:
            continue
        direction = "decreasing" if piece["verdict"] in DECREASING_FAMILY else "increasing"
        if runs and runs[-1]["direction"] == direction and runs[-1]["hi"] == piece["lo"]:
            runs[-1]["hi"] = piece["hi"]
            runs[-1]["pieces"] += 1
        else:
            runs.append(
                {"lo": piece["lo"], "hi": piece["hi"], "direction": direction, "pieces": 1}
            )

    candidate_runs = []
    for run in runs:
        s_lo = endpoint_sign(run["lo"])
        s_hi = endpoint_sign(run["hi"])
        run["sign_lo"] = s_lo
        run["sign_hi"] = s_hi
        if run["direction"] == "decreasing":
            excluded = s_lo == "-" or s_hi == "+"
        else:
            excluded = s_lo == "+" or s_hi == "-"
        run["root_excluded_by_endpoint_signs"] = excluded
        if not excluded:
            candidate_runs.append(run)

    at_most_one = len(candidate_runs) <= 1
    return {
        "at_most_one_on_domain": at_most_one,
        "status": "certified" if at_most_one else "not_established",
        "argument": (
            "root-free pieces carry no root (residual enclosure excludes zero "
            "on the closed piece); each maximal strictly monotone run carries "
            "at most one root; certified endpoint residual signs exclude all "
            f"but {len(candidate_runs)} candidate run(s)"
        ),
        "verdict_counts": verdict_counts,
        "monotone_runs": len(runs),
        "candidate_runs": len(candidate_runs),
    }


def sweep_mode(
    mode: str,
    mp_dps: int,
    iv_dps: int,
    su2_cutoff: int,
    su3_cutoff: int,
    domain_lo: str,
    domain_hi: str,
    initial_pieces: int,
    depth_cap: int,
    eval_budget: int,
) -> dict[str, Any]:
    mpb, point, ivb, chain, alpha_star = prepare_chain(
        mode, mp_dps, iv_dps, su2_cutoff, su3_cutoff
    )

    lo_dom = mpb.num(domain_lo)
    hi_dom = mpb.num(domain_hi)
    if not (lo_dom < hi_dom):
        raise icc.CertificateError("domain is empty")

    # Initial uniform tiling; boundaries are shared mp numbers so the pieces
    # tile the domain exactly.
    bounds = [
        lo_dom + (hi_dom - lo_dom) * mpb.num(k) / mpb.num(initial_pieces)
        for k in range(initial_pieces + 1)
    ]
    stack: list[tuple[Any, Any, int]] = [
        (bounds[k], bounds[k + 1], 0) for k in range(initial_pieces)
    ]
    stack.reverse()  # process left to right (LIFO)

    pieces: list[dict[str, Any]] = []
    evaluations = 0
    max_depth_used = 0
    while stack:
        lo, hi, depth = stack.pop()
        max_depth_used = max(max_depth_used, depth)
        if evaluations >= eval_budget:
            pieces.append(
                {
                    "lo": lo,
                    "hi": hi,
                    "depth": depth,
                    "verdict": VERDICT_UNDECIDED,
                    "reason": "work budget exhausted before evaluation",
                }
            )
            continue
        record = classify_piece(chain, ivb, mpb, lo, hi, mode)
        evaluations += 1
        record["depth"] = depth
        if record["verdict"] in DECISIVE_VERDICTS:
            pieces.append(record)
            continue
        if depth < depth_cap and evaluations < eval_budget:
            mid = (lo + hi) / mpb.two
            stack.append((mid, hi, depth + 1))
            stack.append((lo, mid, depth + 1))
            continue
        pieces.append(record)

    def endpoint_sign(x_mp: Any) -> str | None:
        try:
            g_thin, _ = chain.g_dual(ivb.thin(x_mp), mode)
        except EVAL_ERRORS:
            return None
        residual = g_thin.x - ivb.thin(x_mp)
        if residual.a > ivb.zero.b:
            return "+"
        if residual.b < ivb.zero.a:
            return "-"
        return None

    conclusion = synthesize(pieces, endpoint_sign)

    # Worst-piece statistics.
    worst_l = None
    gp_lo_min = None
    gp_hi_max = None
    min_width = None
    for piece in pieces:
        if "sup_abs_gprime" in piece:
            if worst_l is None or piece["sup_abs_gprime"] > worst_l:
                worst_l = piece["sup_abs_gprime"]
        if "gprime" in piece:
            g_lo = mpb.ctx.mpf(piece["gprime"]["lo"])
            g_hi = mpb.ctx.mpf(piece["gprime"]["hi"])
            gp_lo_min = g_lo if gp_lo_min is None or g_lo < gp_lo_min else gp_lo_min
            gp_hi_max = g_hi if gp_hi_max is None or g_hi > gp_hi_max else gp_hi_max
        width = piece["hi"] - piece["lo"]
        min_width = width if min_width is None or width < min_width else min_width

    def fmt(x: Any, dps: int = 20) -> str:
        return icc._mp_str(mpb, x, dps)

    # Run-length compressed verdict table (verdict blocks in domain order).
    blocks: list[dict[str, Any]] = []
    for piece in pieces:
        if blocks and blocks[-1]["verdict"] == piece["verdict"]:
            blocks[-1]["alpha_hi"] = fmt(piece["hi"])
            blocks[-1]["pieces"] += 1
        else:
            blocks.append(
                {
                    "verdict": piece["verdict"],
                    "alpha_lo": fmt(piece["lo"]),
                    "alpha_hi": fmt(piece["hi"]),
                    "pieces": 1,
                }
            )

    exceptional = [
        {
            "alpha_lo": fmt(p["lo"]),
            "alpha_hi": fmt(p["hi"]),
            "width": fmt(p["hi"] - p["lo"], 8),
            "depth": p["depth"],
            "verdict": p["verdict"],
            "reason": p.get("reason", ""),
        }
        for p in pieces
        if p["verdict"] not in DECISIVE_VERDICTS
    ]
    exceptional_width = mpb.zero
    for p in pieces:
        if p["verdict"] not in DECISIVE_VERDICTS:
            exceptional_width = exceptional_width + (p["hi"] - p["lo"])

    block: dict[str, Any] = {
        "mode": mode,
        "map_definition": icc.MAP_DEFINITIONS[mode],
        "outer_equation": "P = phi + alpha*sqrt(pi)",
        "backend": (
            "mpmath.iv binary interval arithmetic with outward rounding on every "
            "elementary operation; private MPIntervalContext/MPContext instances; "
            "forward-mode interval AD with implicit nodes by the implicit function "
            "theorem, identical to the stage-2 contraction certificate chain"
        ),
        "iv_dps": iv_dps,
        "point_dps": mp_dps,
        "su2_cutoff": su2_cutoff,
        "su3_cutoff": su3_cutoff,
        "edge_sum_tail_bounds": {
            "included": True,
            "method": (
                "geometric majorants on the SU(2)/SU(3) edge-sum truncation tails, "
                "added one-sidedly inside every interval evaluation, so each piece "
                "verdict covers both the declared finite-cutoff sums and their "
                "infinite-cutoff limits"
            ),
        },
        "subdivision": {
            "initial_pieces": initial_pieces,
            "total_pieces": len(pieces),
            "map_evaluations": evaluations,
            "eval_budget": eval_budget,
            "depth_cap": depth_cap,
            "max_depth_used": max_depth_used,
            "min_piece_width": fmt(min_width, 8) if min_width is not None else None,
        },
        "verdict_counts": conclusion.pop("verdict_counts"),
        "worst_piece": {
            "max_sup_abs_gprime": fmt(worst_l) if worst_l is not None else None,
            "gprime_global_hull": {
                "lo": fmt(gp_lo_min) if gp_lo_min is not None else None,
                "hi": fmt(gp_hi_max) if gp_hi_max is not None else None,
            },
        },
        "verdict_blocks": blocks,
        "exceptional_set": exceptional,
        "exceptional_set_total_width": fmt(exceptional_width, 8),
        "conclusion": conclusion,
        "fixed_point_point_estimate_display_only": {
            "alpha": fmt(alpha_star, 30),
            "alpha_inv": fmt(mpb.one / alpha_star, 30),
        },
    }
    return block


def build_certificate(
    mp_dps: int = 40,
    iv_dps: int = 40,
    su2_cutoff: int = 120,
    su3_cutoff: int = 90,
    domain_lo: str = DOMAIN_LO,
    domain_hi: str = DOMAIN_HI,
    initial_pieces: int = 256,
    depth_cap: int = 8,
    eval_budget: int = 2000,
    modes: tuple[str, ...] = icc.MODES,
) -> dict[str, Any]:
    mode_blocks = {}
    for mode in modes:
        mode_blocks[mode] = sweep_mode(
            mode,
            mp_dps,
            iv_dps,
            su2_cutoff,
            su3_cutoff,
            domain_lo,
            domain_hi,
            initial_pieces,
            depth_cap,
            eval_budget,
        )

    all_at_most_one = all(
        block["conclusion"]["at_most_one_on_domain"] for block in mode_blocks.values()
    )
    artifact: dict[str, Any] = {
        "artifact": "oph_p_global_uniqueness_certificate",
        "date": ARTIFACT_DATE,
        "claim_status": "global_at_most_one_certificate_for_declared_closure_map",
        "claim_boundary": (
            "At-most-one fixed point on the declared physical domain is certified "
            "by adaptive interval subdivision for the declared numerical map at "
            "the stated representation cutoffs (the edge-sum tail majorants extend "
            "every enclosure to the infinite-cutoff sums). Combined with the "
            "stage-2 contraction certificate (existence on a subinterval, same "
            "declared map and cutoffs), each readout map has exactly one fixed "
            "point on the domain. The declared one-loop RG/matching conventions, "
            "the tree-level m_Z closure, the Stage-5 continuation masses, and the "
            "exact one-loop kernel are certified as declared numerical structure, "
            "not as physical endpoint theorems. This is not an exact fine-structure "
            "derivation; the stage-3 landing verdict of closure row CL-1 is "
            "unchanged."
        ),
        "protocol_stage": (
            "global uniqueness supplement to stage 2 of the basin-then-contract "
            "protocol for the P coordinate; see the P-closure issue (#545)"
        ),
        "domain": {
            "alpha": {"lo": domain_lo, "hi": domain_hi},
            "alpha_inv": {
                "lo": format(Decimal(1) / Decimal(domain_hi), "f"),
                "hi": format(Decimal(1) / Decimal(domain_lo), "f"),
            },
            "declaration": DOMAIN_DECLARATION,
        },
        "stage2_certificate": (
            "runtime/p_interval_contraction_certificate_2026-07-14.json "
            "(existence and local uniqueness on the contraction interval)"
        ),
        "inner_root_scope": (
            "g is the declared map: the inner alpha_U and m_Z roots are the ones "
            "the declared scan-and-bisect solver selects; each piece evaluation "
            "encloses those roots in verified sign-change brackets with "
            "sign-definite residual derivatives; global uniqueness of the inner "
            "roots over the solver scan windows is not claimed"
        ),
        "promotion_allowed": False,
        "exact_alpha_promoted": False,
        "consumer_policy": {
            "may_feed_live_particle_predictions": False,
            "may_feed_compare_or_audit_surfaces": True,
            "hidden_external_alpha_allowed": False,
            "default_thomson_endpoint_allowed": False,
        },
        "conclusion": {
            "at_most_one_on_domain_all_modes": all_at_most_one,
            "statement": (
                "each declared readout map has at most one fixed point on alpha "
                "in [" + domain_lo + ", " + domain_hi + "]; with stage-2 existence "
                "this is exactly one fixed point per map on the declared domain"
                if all_at_most_one
                else "at-most-one is not established on the full declared domain; "
                "see the per-mode exceptional sets"
            ),
        },
        "modes": mode_blocks,
    }
    return artifact


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Emit the global at-most-one fixed-point certificate for the OPH "
            "P/alpha closure map on the declared physical domain."
        )
    )
    parser.add_argument("--mp-dps", type=int, default=40, help="Point-solver working digits.")
    parser.add_argument("--iv-dps", type=int, default=40, help="Interval arithmetic working digits.")
    parser.add_argument("--su2-cutoff", type=int, default=120, help="SU(2) representation cutoff.")
    parser.add_argument("--su3-cutoff", type=int, default=90, help="SU(3) representation cutoff.")
    parser.add_argument("--domain-lo", default=DOMAIN_LO, help="Domain lower alpha endpoint.")
    parser.add_argument("--domain-hi", default=DOMAIN_HI, help="Domain upper alpha endpoint.")
    parser.add_argument(
        "--initial-pieces",
        type=int,
        default=256,
        help="Uniform initial subdivision piece count.",
    )
    parser.add_argument(
        "--depth-cap",
        type=int,
        default=8,
        help="Maximum bisection depth below the initial subdivision.",
    )
    parser.add_argument(
        "--eval-budget",
        type=int,
        default=2000,
        help="Deterministic work budget: maximum interval map evaluations per mode.",
    )
    parser.add_argument(
        "--mode",
        choices=icc.MODES + ("both",),
        default="both",
        help="Which readout map(s) to certify.",
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output JSON path.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    modes = icc.MODES if args.mode == "both" else (args.mode,)
    t0 = time.time()
    artifact = build_certificate(
        mp_dps=args.mp_dps,
        iv_dps=args.iv_dps,
        su2_cutoff=args.su2_cutoff,
        su3_cutoff=args.su3_cutoff,
        domain_lo=args.domain_lo,
        domain_hi=args.domain_hi,
        initial_pieces=args.initial_pieces,
        depth_cap=args.depth_cap,
        eval_budget=args.eval_budget,
        modes=modes,
    )
    wall_seconds = round(time.time() - t0, 2)
    text = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    for mode in modes:
        block = artifact["modes"][mode]
        print(f"{mode}:")
        print(f"  pieces                = {block['subdivision']['total_pieces']}")
        print(f"  verdict counts        = {block['verdict_counts']}")
        print(f"  worst sup|g'|         = {block['worst_piece']['max_sup_abs_gprime']}")
        print(f"  exceptional set width = {block['exceptional_set_total_width']}")
        print(f"  at most one on domain = {block['conclusion']['at_most_one_on_domain']}")
    print(f"conclusion (all modes): {artifact['conclusion']['at_most_one_on_domain_all_modes']}")
    print(f"wall seconds: {wall_seconds} (not recorded in the artifact)")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
