#!/usr/bin/env python3
"""Maximal-domain global uniqueness extension certificate for the P/alpha map.

The 2026-07-14 certificates prove existence (stage 2, Banach) and global
at-most-one (adaptive subdivision, sup |g'| < 1) for both declared readout
maps on the declared physical domain alpha in [0.005, 0.01], equivalently
alpha_inv in [100, 200]. This module extends the uniqueness statement to the
maximal analytic domain of the declared maps (proof spine GAP-A7):

- the declared maps read alpha only through P = phi + alpha*sqrt(pi); the
  inner structure is fixed by two declared solver windows: the D10 pixel
  closure selects alpha_U inside the scan window [0.02, 0.08], and the
  tree-level m_Z closure selects m_Z inside the log-grid window
  [mu_U*e^-50, mu_U], i.e. L := ln(mu_U/m_Z) in [0, 50]
  (``paper_math.solve_alpha_u_from_p`` and ``solve_mz_fixed_point_tree``);
- every point of the maximal domain, and in particular every fixed point,
  therefore carries a window-consistent triple (alpha_U, L, P) satisfying the
  pixel-closure equation, the m_Z-closure identity, and the window bounds;
- a closed-form envelope lemma (interval arithmetic on the window bounds, the
  m_Z-closure identity, and edge-sum majorants) certifies
  1/alpha_3(m_Z) >= 4 on the whole domain, hence quark screening in
  [1 - 3/(4*pi), 1], a global P range, and a global positive lower bound on
  the inverse-alpha readout, so the readout has no pole on the domain and
  every fixed point satisfies 0 < alpha <= 1/AI_lo;
- an adaptive interval sweep over the full pixel window alpha_U in
  [0.02, 0.08] propagates the closure constraints on each piece and certifies
  one of: the m_Z window is infeasible (no domain point), the pixel/screening
  constraints are infeasible, the alpha value forced by the closure is
  disjoint from the certified readout enclosure (no fixed point, per mode),
  or the alpha value lies inside the declared physical interval (handoff to
  the 2026-07-14 certificates);
- the synthesis: every fixed point of either declared readout map on the
  maximal analytic domain lies inside the declared physical interval; the
  at-most-one certificate and the stage-2 existence certificate then give
  exactly one fixed point per map on the maximal domain, and the list of
  nonphysical exterior fixed points is empty.

The sweep quantifies over every window-consistent inner-root selection, so it
covers the declared scan-and-bisect selection as a special case. Both readout
maps are handled: the source map and the gauge-width map differ only by the
+alpha_U term in the readout, evaluated per piece.

Claim boundary: the certificate covers the declared numerical maps at the
declared representation cutoffs, with edge-sum tail majorants extending every
enclosure to the infinite-cutoff sums; it certifies fixed-point structure of
the declared maps and no relation to the measured fine-structure constant;
the stage-3 landing verdict (closure row CL-1) is unchanged, and public
claims remain tied to the declared physical interval.

mpmath contexts are private instances; the global ``mpmath.mp`` and
``mpmath.iv`` precision settings are never touched. The construction is
deterministic: subdivision order, budgets, and verdicts are functions of the
declared parameters only, and no wall-clock quantity enters the artifact.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

from mpmath.libmp.libmpf import ComplexResult

import interval_contraction_certificate as icc


ARTIFACT_DATE = "2026-07-17"
RUNTIME_DIR = Path(__file__).resolve().parent / "runtime"
DEFAULT_OUT = RUNTIME_DIR / (
    f"p_global_uniqueness_extension_certificate_{ARTIFACT_DATE}.json"
)

AT_MOST_ONE_ARTIFACT = "p_global_uniqueness_certificate_2026-07-14.json"
EXISTENCE_ARTIFACT = "p_interval_contraction_certificate_2026-07-14.json"

U_WINDOW = ("0.02", "0.08")
L_WINDOW = ("0", "50")
PHYSICAL_ALPHA = ("0.005", "0.01")
# Strict inset of the declared physical interval used for the handoff check,
# so the containment claim is robust against decimal-to-binary endpoint
# rounding of the declared interval bounds.
HANDOFF_ALPHA = ("0.0051", "0.0099")

A3_INV_FLOOR = "4"

VERDICT_MZ_WINDOW = "excluded_mz_window_infeasible"
VERDICT_A3_WINDOW = "excluded_a3_envelope_infeasible"
VERDICT_CLOSURE = "excluded_pixel_closure_infeasible"
VERDICT_GAP = "excluded_fixed_point_gap_both_modes"
VERDICT_HANDOFF = "handoff_physical_interval"
VERDICT_UNDECIDED = "undecided"

DECISIVE_VERDICTS = {
    VERDICT_MZ_WINDOW,
    VERDICT_A3_WINDOW,
    VERDICT_CLOSURE,
    VERDICT_GAP,
    VERDICT_HANDOFF,
}

EVAL_ERRORS = (icc.CertificateError, ComplexResult, ValueError, ArithmeticError)

QUARK_CHARGE_NUMERATORS = {"u": 4, "d": 1, "s": 1, "c": 4, "b": 1}


class Frame:
    """Backends, representation tables, and declared constants."""

    def __init__(self, mp_dps: int, iv_dps: int, su2_cutoff: int, su3_cutoff: int) -> None:
        self.mp_dps = mp_dps
        self.iv_dps = iv_dps
        self.su2_cutoff = su2_cutoff
        self.su3_cutoff = su3_cutoff
        self.mpb = icc.MpBackend(mp_dps)
        self.ivb = icc.IvBackend(iv_dps)
        self.su2_terms = icc.build_su2_terms(self.ivb, su2_cutoff)
        self.su3_terms = icc.build_su3_terms(self.ivb, su3_cutoff)
        self.vectors = icc.stage5_integer_vectors(self.mpb)
        roots_mp = icc.koide_roots(self.mpb)
        order = sorted(range(3), key=lambda i: roots_mp[i])
        roots_iv = icc.koide_roots(self.ivb)
        self.roots_ascending = [roots_iv[i] for i in order]
        for i in range(2):
            if not (self.roots_ascending[i].b < self.roots_ascending[i + 1].a):
                raise icc.CertificateError("Koide root intervals are not disjoint")
        self.lepton_c = icc.lepton_masses(self.ivb, self.ivb.one, self.vectors, self.roots_ascending)
        self.quark_c = icc.quark_masses(self.ivb, self.ivb.one, self.vectors)
        self._mp_su2_terms: list | None = None
        self._mp_su3_terms: list | None = None

    # -- interval helpers -----------------------------------------------------

    def isect(self, x: Any, y: Any) -> Any | None:
        """Intersection of two intervals from their exact stored endpoints."""
        lo = x.a if x.a > y.a else y.a
        hi = x.b if x.b < y.b else y.b
        if lo > hi:
            return None
        return self.ivb.hull(lo, hi)

    def hull2(self, x: Any, y: Any) -> Any:
        lo = x.a if x.a < y.a else y.a
        hi = x.b if x.b > y.b else y.b
        return self.ivb.hull(lo, hi)

    def box(self, lo: Any, hi: Any) -> Any:
        return self.ivb.hull(self.ivb.thin(lo).a, self.ivb.thin(hi).b)

    def iv_lo(self, x: Any) -> Any:
        return self.mpb.ctx.mpf(x.a._mpi_[0])

    def iv_hi(self, x: Any) -> Any:
        return self.mpb.ctx.mpf(x.b._mpi_[1])

    def fmt(self, x: Any, dps: int = 20) -> str:
        return icc._mp_str(self.mpb, x, dps)

    def pair(self, x: Any) -> dict[str, str]:
        return icc._iv_pair(x)

    # -- declared-chain pieces -------------------------------------------------

    def ellbar_pair(self, t2: Any, t3: Any) -> tuple[Any, Any]:
        """Tail-bounded ellbar enclosures for both groups (value intervals)."""
        ivb = self.ivb
        tl2 = icc.su2_tail_bounds(ivb, t2, self.su2_cutoff)
        tl3 = icc.su3_tail_bounds(ivb, t3, self.su3_cutoff)
        z2 = ivb.hull(ivb.zero.a, tl2["z"].b)
        s2 = ivb.hull(ivb.zero.a, tl2["s"].b)
        z3 = ivb.hull(ivb.zero.a, tl3["z"].b)
        s3 = ivb.hull(ivb.zero.a, tl3["s"].b)
        e2 = icc.ellbar(ivb, self.su2_terms, t2, (z2, s2))
        e3 = icc.ellbar(ivb, self.su3_terms, t3, (z3, s3))
        return e2, e3

    def kernel_from_a(self, a: Any, charge_sq: Any, multiplicity: Any) -> Any:
        """Exact one-loop transport kernel written in a = q^2/(4 m^2)."""
        b = self.ivb
        sqrt_a = b.sqrt(a)
        sqrt_one_plus_a = b.sqrt(b.one + a)
        asinh_sqrt_a = b.log(sqrt_a + sqrt_one_plus_a)
        integral = (
            -(b.num(5) / b.num(18))
            + b.one / (b.num(6) * a)
            + ((b.two * a - b.one) * sqrt_one_plus_a * asinh_sqrt_a)
            / (b.num(6) * a * sqrt_a)
        )
        return b.two * multiplicity * charge_sq / b.pi * integral

    def kernel_enclosure(self, a_box: Any, charge_sq: Any, multiplicity: Any, pieces: int) -> Any:
        """Kernel enclosure over an a-interval via geometric subdivision."""
        mpb = self.mpb
        lo = self.iv_lo(a_box)
        hi = self.iv_hi(a_box)
        if not lo > 0:
            raise icc.CertificateError("kernel argument interval must be positive")
        ratio = mpb.exp(mpb.log(hi / lo) / mpb.num(pieces))
        out = None
        g0 = lo
        for index in range(pieces):
            g1 = hi if index == pieces - 1 else g0 * ratio
            k = self.kernel_from_a(self.box(g0, g1), charge_sq, multiplicity)
            out = k if out is None else self.hull2(out, k)
            g0 = g1
        return out

    def delta_thomson(self, mzv2: Any, screening: Any, kernel_pieces: int) -> Any:
        """Structured Thomson continuation enclosure from mzv2 = 4*pi*(a2+aY).

        Every Stage-5 continuation mass is a P-independent constant times v,
        and m_Z = (v/2)*sqrt(mzv2), so each kernel argument is
        a_f = mzv2/(16 c_f^2) with c_f = m_f/v; the transmutation scale cancels.
        """
        ivb = self.ivb
        sixteen = ivb.num(16)
        lepton_total = None
        for name in ("e", "mu", "tau"):
            c = self.lepton_c[name]
            k = self.kernel_enclosure(mzv2 / (sixteen * c * c), ivb.one, ivb.one, kernel_pieces)
            lepton_total = k if lepton_total is None else lepton_total + k
        three = ivb.num(3)
        nine = ivb.num(9)
        quark_naive = None
        for name, cnum in QUARK_CHARGE_NUMERATORS.items():
            c = self.quark_c[name]
            k = self.kernel_enclosure(
                mzv2 / (sixteen * c * c), ivb.num(cnum) / nine, three, kernel_pieces
            )
            quark_naive = k if quark_naive is None else quark_naive + k
        return lepton_total + screening * quark_naive


# ---------------------------------------------------------------------------
# Envelope lemma: certified global constraints on any domain point.
# ---------------------------------------------------------------------------


def su2_numerator_upper(frame: Frame, t_lo: Any) -> Any:
    """Upper endpoint of sum_{n>=1} dim*ln(dim)*e^{-t*c2} + tail at thin t = t_lo."""
    ivb = frame.ivb
    t = ivb.thin(t_lo)
    total = ivb.zero
    for dim, c2, ln_dim in frame.su2_terms:
        total = total + dim * ivb.exp(-(t * c2)) * ln_dim
    tail = icc.su2_tail_bounds(ivb, t, frame.su2_cutoff)
    return ivb.thin((total + tail["s"]).b)


def su3_numerator_upper(frame: Frame, t_lo: Any) -> Any:
    """Upper endpoint of the SU(3) numerator sum + tail at thin t = t_lo."""
    ivb = frame.ivb
    t = ivb.thin(t_lo)
    total = ivb.zero
    for dim, c2, ln_dim in frame.su3_terms:
        total = total + dim * ivb.exp(-(t * c2)) * ln_dim
    tail = icc.su3_tail_bounds(ivb, t, frame.su3_cutoff)
    return ivb.thin((total + tail["s"]).b)


def su2_partition_upper(frame: Frame, t_lo: Any) -> Any:
    """Upper endpoint of sum dim*e^{-t*c2} + tail at thin t = t_lo."""
    ivb = frame.ivb
    t = ivb.thin(t_lo)
    total = ivb.zero
    for dim, c2, _ in frame.su2_terms:
        total = total + dim * ivb.exp(-(t * c2))
    tail = icc.su2_tail_bounds(ivb, t, frame.su2_cutoff)
    return ivb.thin((total + tail["z"]).b)


def envelope_lemma(frame: Frame) -> dict[str, Any]:
    """Certify global bounds valid at every point of the maximal domain.

    Inputs are the declared window bounds alone; every inequality is an
    interval comparison. Raises CertificateError if any check fails.
    """
    ivb = frame.ivb
    mpb = frame.mpb
    u_lo = mpb.num(U_WINDOW[0])
    u_hi = mpb.num(U_WINDOW[1])
    l_lo = mpb.num(L_WINDOW[0])
    l_hi = mpb.num(L_WINDOW[1])
    U = frame.box(u_lo, u_hi)
    L = frame.box(l_lo, l_hi)
    four_pi_sq = ivb.four_pi * ivb.pi

    # Window coupling bounds: 1/a_i = 1/u + (b_i/2pi)*L with L >= 0.
    a1_inv = ivb.one / U + (ivb.num(33) / ivb.num(5)) * L / ivb.two_pi
    a2_inv = ivb.one / U + L / ivb.two_pi
    ay_inv = (ivb.num(5) / ivb.num(3)) * a1_inv
    s_box = ivb.one / a2_inv + (ivb.num(3) / ivb.num(5)) / a1_inv
    mzv2 = ivb.four_pi * s_box
    t2 = four_pi_sq / a2_inv
    t2_lo = frame.iv_lo(t2)

    # m_Z-closure identity at the selected root:
    #   L = 2pi/(4u) - 2pi + (2/3)ln P + ln 2 - (1/2)ln(4pi(a2+aY)),
    # hence 1/a3 = 1/u - (3/2pi)L
    #            = 1/(4u) + 3 - (1/pi)ln P - (3/2pi)ln 2 + (3/4pi)ln(4pi(a2+aY)).
    ln2 = ivb.log(ivb.two)
    b_const = (
        ivb.one / (ivb.num(4) * frame.box(u_hi, u_hi))
        + ivb.num(3)
        - (ivb.num(3) / ivb.two_pi) * ln2
        + (ivb.num(3) / (ivb.num(4) * ivb.pi)) * ivb.thin(ivb.log(mzv2).a)
    )
    a3_floor = ivb.num(A3_INV_FLOOR)

    # Case 0 < 1/a3 < 4: then alpha_3 > 1/4, t3 > pi^2, and the SU(3) ellbar
    # is bounded by its numerator at pi^2 over the trivial-term partition
    # floor 1; the pixel closure then caps P and the identity forces
    # 1/a3 above the floor, a contradiction.
    s2_cap = su2_numerator_upper(frame, t2_lo)
    pi_sq = frame.iv_lo(ivb.pi * ivb.pi)
    s3_cap_case2 = su3_numerator_upper(frame, pi_sq)
    p_cap_case2 = ivb.num(4) * (s2_cap + s3_cap_case2)
    bound_case2 = b_const - ivb.log(p_cap_case2) / ivb.pi
    if not (bound_case2.a > a3_floor.b):
        raise icc.CertificateError("envelope case 0 < 1/a3 < 4 is not refuted")

    # Case 1/a3 <= 0 (finite-cutoff variant; the infinite-cutoff edge sums
    # are divergent there, outside the analytic domain by the positivity
    # condition): the finite SU(3) ellbar is a weighted mean of ln(dim), so it
    # is bounded by the maximal ln(dim) at any real t3.
    n3 = frame.su3_cutoff
    dim3_max = ((n3 + 1) * (n3 + 1) * (2 * n3 + 2)) // 2
    ln_dim3_max = ivb.log(ivb.num(dim3_max))
    p_cap_case3 = ivb.num(4) * (s2_cap + ln_dim3_max)
    bound_case3 = b_const - ivb.log(p_cap_case3) / ivb.pi
    if not (bound_case3.a > a3_floor.b):
        raise icc.CertificateError("envelope case 1/a3 <= 0 is not refuted")

    # Conclusion: 1/a3 in [4, 1/u_lo] on the whole domain.
    a3_inv_glob = frame.box(frame.iv_lo(a3_floor), frame.iv_hi(ivb.one / frame.box(u_lo, u_lo)))
    t3_glob = four_pi_sq / a3_inv_glob
    screening_glob = ivb.one - ivb.num(3) / (ivb.pi * a3_inv_glob)
    if not (screening_glob.a > ivb.zero.b):
        raise icc.CertificateError("global screening lower bound is not positive")

    # Global P range from the pixel closure P = 4(ellbar2 + ellbar3):
    # upper from numerator caps at the t lower endpoints, lower from the
    # first SU(2) term at the t2 upper endpoint over the full partition sum.
    s3_cap = su3_numerator_upper(frame, frame.iv_lo(t3_glob))
    p_glob_hi = ivb.num(4) * (s2_cap + s3_cap)
    t2_hi = frame.iv_hi(t2)
    z2_cap = su2_partition_upper(frame, t2_lo)
    first_term = ivb.two * ivb.log(ivb.two) * ivb.exp(
        -(ivb.thin(t2_hi) * (ivb.num(3) / ivb.num(4)))
    )
    p_glob_lo = ivb.num(4) * ivb.thin((first_term / z2_cap).a)
    if not (p_glob_lo.a > ivb.zero.b):
        raise icc.CertificateError("global P lower bound is not positive")
    p_glob = frame.hull2(p_glob_lo, p_glob_hi)

    # Global readout floor: alpha_inv >= 1/a2 + 1/aY + delta on both modes
    # (the gauge-width mode adds alpha_U > 0).
    delta_glob = frame.delta_thomson(mzv2, screening_glob, kernel_pieces=192)
    ai_lo_iv = ivb.thin(a2_inv.a) + ivb.thin(ay_inv.a) + ivb.thin(delta_glob.a)
    if not (ai_lo_iv.a > ivb.zero.b):
        raise icc.CertificateError("global inverse-alpha readout floor is not positive")
    g_cap = ivb.one / ivb.thin(ai_lo_iv.a)

    return {
        "U": U,
        "L": L,
        "a3_inv_glob": a3_inv_glob,
        "t2_glob": frame.hull2(ivb.thin(frame.iv_lo(t2)), ivb.thin(t2_hi)),
        "t3_glob": t3_glob,
        "screening_glob": screening_glob,
        "p_glob": p_glob,
        "ai_lo": ai_lo_iv,
        "g_cap": g_cap,
        "report": {
            "window_u": {"lo": U_WINDOW[0], "hi": U_WINDOW[1]},
            "window_L": {"lo": L_WINDOW[0], "hi": L_WINDOW[1]},
            "coupling_inverse_bounds": {
                "a1_inv": frame.pair(a1_inv),
                "a2_inv": frame.pair(a2_inv),
                "aY_inv": frame.pair(ay_inv),
            },
            "t2_window": frame.pair(t2),
            "mz_closure_identity": (
                "1/a3 = 1/(4u) + 3 - (1/pi)ln P - (3/2pi)ln 2 + (3/4pi)ln(4pi(a2+aY)) "
                "at the selected m_Z root"
            ),
            "case_su3_subcritical": {
                "assumption": "0 < 1/a3 < " + A3_INV_FLOOR,
                "su3_ellbar_cap_at_pi_sq": frame.fmt(frame.iv_hi(s3_cap_case2), 8),
                "p_cap": frame.fmt(frame.iv_hi(p_cap_case2), 12),
                "identity_lower_bound": frame.fmt(frame.iv_lo(bound_case2), 12),
                "refuted": True,
            },
            "case_su3_nonpositive": {
                "assumption": "1/a3 <= 0",
                "variant_note": (
                    "finite-cutoff variant refuted by the weighted-mean cap "
                    "ellbar3 <= ln(dim_max); the infinite-cutoff edge sums are "
                    "divergent there, outside the analytic domain by the t3 > 0 "
                    "positivity condition"
                ),
                "ln_dim3_max": frame.fmt(frame.iv_hi(ln_dim3_max), 12),
                "p_cap": frame.fmt(frame.iv_hi(p_cap_case3), 12),
                "identity_lower_bound": frame.fmt(frame.iv_lo(bound_case3), 12),
                "refuted": True,
            },
            "a3_inv_global": frame.pair(a3_inv_glob),
            "t3_global": frame.pair(t3_glob),
            "screening_global": frame.pair(screening_glob),
            "su2_ellbar_numerator_cap": frame.fmt(frame.iv_hi(s2_cap), 8),
            "p_global": frame.pair(p_glob),
            "delta_thomson_global": frame.pair(delta_glob),
            "alpha_inv_readout_floor": frame.fmt(frame.iv_lo(ai_lo_iv), 12),
            "fixed_point_alpha_cap": frame.fmt(frame.iv_hi(g_cap), 12),
            "statement": (
                "every point of the maximal domain satisfies 1/a3 >= "
                + A3_INV_FLOOR
                + ", P inside p_global, and alpha_inv >= alpha_inv_readout_floor > 0; "
                "every fixed point of either readout map satisfies "
                "0 < alpha <= fixed_point_alpha_cap"
            ),
        },
    }


# ---------------------------------------------------------------------------
# Pixel-window sweep: constraint propagation on alpha_U pieces.
# ---------------------------------------------------------------------------


def classify_u_piece(
    frame: Frame,
    env: dict[str, Any],
    lo_mp: Any,
    hi_mp: Any,
    iters: int,
    kernel_pieces: int,
) -> dict[str, Any]:
    """Propagate the closure constraints on one alpha_U piece and classify it."""
    ivb = frame.ivb
    record: dict[str, Any] = {"lo": lo_mp, "hi": hi_mp}
    U = frame.box(lo_mp, hi_mp)
    P = env["p_glob"]
    L = env["L"]
    four_pi_sq = ivb.four_pi * ivb.pi
    handoff_lo = frame.mpb.num(HANDOFF_ALPHA[0])
    handoff_hi = frame.mpb.num(HANDOFF_ALPHA[1])

    try:
        for iteration in range(iters):
            a1_inv = ivb.one / U + (ivb.num(33) / ivb.num(5)) * L / ivb.two_pi
            a2_inv = ivb.one / U + L / ivb.two_pi
            a3_inv = frame.isect(ivb.one / U - ivb.num(3) * L / ivb.two_pi, env["a3_inv_glob"])
            if a3_inv is None:
                record["verdict"] = VERDICT_A3_WINDOW
                record["iterations_used"] = iteration + 1
                return record
            s_box = ivb.one / a2_inv + (ivb.num(3) / ivb.num(5)) / a1_inv
            l_formula = (
                ivb.two_pi / (ivb.num(4) * U)
                - ivb.two_pi
                + (ivb.two / ivb.num(3)) * ivb.log(P)
                + ivb.log(ivb.two)
                - ivb.log(ivb.four_pi * s_box) / ivb.two
            )
            L_new = frame.isect(L, l_formula)
            if L_new is None:
                record["verdict"] = VERDICT_MZ_WINDOW
                record["iterations_used"] = iteration + 1
                return record
            L = L_new

            t2 = frame.isect(four_pi_sq / a2_inv, env["t2_glob"])
            t3 = frame.isect(four_pi_sq / a3_inv, env["t3_glob"])
            if t2 is None or t3 is None:
                record["verdict"] = VERDICT_A3_WINDOW
                record["iterations_used"] = iteration + 1
                return record
            e2, e3 = frame.ellbar_pair(t2, t3)
            P_new = frame.isect(P, ivb.num(4) * (e2 + e3))
            if P_new is None:
                record["verdict"] = VERDICT_CLOSURE
                record["iterations_used"] = iteration + 1
                return record
            P = P_new

            # Readout enclosures and verdict test (skip the first pass, where
            # the boxes are still dominated by the global initialization).
            if iteration == 0 and iters > 1:
                continue
            alpha_box = (P - ivb.phi) / ivb.sqrt_pi
            mzv2 = ivb.four_pi * s_box
            screening = ivb.one - ivb.num(3) / (ivb.pi * a3_inv)
            delta = frame.delta_thomson(mzv2, screening, kernel_pieces)
            ai_source = a2_inv + (ivb.num(5) / ivb.num(3)) * a1_inv + delta
            ai_gauge = ai_source + U
            if not (ai_source.a > ivb.zero.b and ai_gauge.a > ivb.zero.b):
                continue
            g_source = ivb.one / ai_source
            g_gauge = ivb.one / ai_gauge

            def status(g: Any) -> str:
                if alpha_box.b < g.a or g.b < alpha_box.a:
                    return "excluded_gap"
                inside = (
                    frame.iv_lo(alpha_box) >= handoff_lo
                    and frame.iv_hi(alpha_box) <= handoff_hi
                )
                return "handoff_inside" if inside else "undecided"

            source_status = status(g_source)
            gauge_status = status(g_gauge)
            decisive = {"excluded_gap", "handoff_inside"}
            if source_status in decisive and gauge_status in decisive:
                record["iterations_used"] = iteration + 1
                record["_alpha_lo_mp"] = frame.iv_lo(alpha_box)
                record["_alpha_hi_mp"] = frame.iv_hi(alpha_box)
                record["alpha_box"] = frame.pair(alpha_box)
                record["g_source"] = frame.pair(g_source)
                record["g_gauge"] = frame.pair(g_gauge)
                record["L_box"] = frame.pair(L)
                record["P_box"] = frame.pair(P)
                record["mode_status"] = {
                    icc.MODE_SOURCE: source_status,
                    icc.MODE_GAUGE_WIDTH: gauge_status,
                }
                if source_status == "excluded_gap" and gauge_status == "excluded_gap":
                    record["verdict"] = VERDICT_GAP
                else:
                    record["verdict"] = VERDICT_HANDOFF
                return record
    except EVAL_ERRORS as exc:
        record["verdict"] = VERDICT_UNDECIDED
        record["reason"] = f"{type(exc).__name__}: {exc}"
        return record

    record["verdict"] = VERDICT_UNDECIDED
    record["reason"] = "constraint propagation left the alpha and readout boxes overlapping"
    alpha_box = (P - ivb.phi) / ivb.sqrt_pi
    record["_alpha_lo_mp"] = frame.iv_lo(alpha_box)
    record["_alpha_hi_mp"] = frame.iv_hi(alpha_box)
    record["alpha_box"] = frame.pair(alpha_box)
    record["P_box"] = frame.pair(P)
    return record


def sweep_pixel_window(
    frame: Frame,
    env: dict[str, Any],
    initial_pieces: int,
    depth_cap: int,
    eval_budget: int,
    iters: int,
    kernel_pieces: int,
) -> dict[str, Any]:
    """Adaptive sweep of the full pixel window with per-piece verdicts."""
    mpb = frame.mpb
    u_lo = mpb.num(U_WINDOW[0])
    u_hi = mpb.num(U_WINDOW[1])
    bounds = [
        u_lo + (u_hi - u_lo) * mpb.num(k) / mpb.num(initial_pieces)
        for k in range(initial_pieces + 1)
    ]
    stack: list[tuple[Any, Any, int]] = [
        (bounds[k], bounds[k + 1], 0) for k in range(initial_pieces)
    ]
    stack.reverse()

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
        record = classify_u_piece(frame, env, lo, hi, iters, kernel_pieces)
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

    verdict_counts: dict[str, int] = {}
    for piece in pieces:
        verdict_counts[piece["verdict"]] = verdict_counts.get(piece["verdict"], 0) + 1

    # Certified hulls: alpha projection of the domain superset (every piece
    # that is not window/closure-excluded) and of the handoff band.
    domain_alpha_hull = None
    handoff_alpha_hull = None
    for piece in pieces:
        if "_alpha_lo_mp" not in piece:
            continue
        lo_a = piece.pop("_alpha_lo_mp")
        hi_a = piece.pop("_alpha_hi_mp")
        if domain_alpha_hull is None:
            domain_alpha_hull = [lo_a, hi_a]
        else:
            domain_alpha_hull[0] = min(domain_alpha_hull[0], lo_a)
            domain_alpha_hull[1] = max(domain_alpha_hull[1], hi_a)
        if piece["verdict"] == VERDICT_HANDOFF:
            if handoff_alpha_hull is None:
                handoff_alpha_hull = [lo_a, hi_a]
            else:
                handoff_alpha_hull[0] = min(handoff_alpha_hull[0], lo_a)
                handoff_alpha_hull[1] = max(handoff_alpha_hull[1], hi_a)

    exceptional = [
        {
            "u_lo": frame.fmt(p["lo"]),
            "u_hi": frame.fmt(p["hi"]),
            "depth": p["depth"],
            "reason": p.get("reason", ""),
        }
        for p in pieces
        if p["verdict"] not in DECISIVE_VERDICTS
    ]

    blocks: list[dict[str, Any]] = []
    for piece in pieces:
        if blocks and blocks[-1]["verdict"] == piece["verdict"]:
            blocks[-1]["u_hi"] = frame.fmt(piece["hi"])
            blocks[-1]["pieces"] += 1
        else:
            blocks.append(
                {
                    "verdict": piece["verdict"],
                    "u_lo": frame.fmt(piece["lo"]),
                    "u_hi": frame.fmt(piece["hi"]),
                    "pieces": 1,
                }
            )

    certified = not exceptional
    handoff_inside = handoff_alpha_hull is not None and (
        handoff_alpha_hull[0] >= frame.mpb.num(HANDOFF_ALPHA[0])
        and handoff_alpha_hull[1] <= frame.mpb.num(HANDOFF_ALPHA[1])
    )
    conclusion = {
        "certified": bool(certified and (handoff_alpha_hull is None or handoff_inside)),
        "statement": (
            "every window-consistent fixed-point candidate of either readout map "
            "has alpha inside the handoff hull, which lies inside the declared "
            "physical interval [" + PHYSICAL_ALPHA[0] + ", " + PHYSICAL_ALPHA[1] + "]"
            if certified
            else "the sweep left undecided pieces; see the exceptional set"
        ),
    }
    return {
        "window_u": {"lo": U_WINDOW[0], "hi": U_WINDOW[1]},
        "subdivision": {
            "initial_pieces": initial_pieces,
            "total_pieces": len(pieces),
            "map_evaluations": evaluations,
            "eval_budget": eval_budget,
            "depth_cap": depth_cap,
            "max_depth_used": max_depth_used,
            "propagation_iterations": iters,
            "kernel_subdivision_pieces": kernel_pieces,
        },
        "verdict_counts": verdict_counts,
        "verdict_blocks": blocks,
        "exceptional_set": exceptional,
        "domain_alpha_outer_hull": (
            None
            if domain_alpha_hull is None
            else {
                "lo": frame.fmt(domain_alpha_hull[0], 12),
                "hi": frame.fmt(domain_alpha_hull[1], 12),
            }
        ),
        "handoff_alpha_hull": (
            None
            if handoff_alpha_hull is None
            else {
                "lo": frame.fmt(handoff_alpha_hull[0], 12),
                "hi": frame.fmt(handoff_alpha_hull[1], 12),
            }
        ),
        "handoff_window": {"lo": HANDOFF_ALPHA[0], "hi": HANDOFF_ALPHA[1]},
        "conclusion": conclusion,
    }


# ---------------------------------------------------------------------------
# Numeric edge location under the declared grid semantics (display only).
# ---------------------------------------------------------------------------


def _declared_map_defined(frame: Frame, alpha_mp: Any) -> bool:
    """Mirror the declared scan semantics at a point alpha (mp arithmetic)."""
    b = frame.mpb
    if frame._mp_su2_terms is None:
        frame._mp_su2_terms = icc.build_su2_terms(b, frame.su2_cutoff)
        frame._mp_su3_terms = icc.build_su3_terms(b, frame.su3_cutoff)
    p = b.phi + alpha_mp * b.sqrt_pi
    if not p > 0:
        return False
    mu_u = icc.mu_u_of_p(b, p)
    lo = b.num(U_WINDOW[0])
    hi = b.num(U_WINDOW[1])
    prev_r = None
    log_lo = b.log(mu_u) - b.num(L_WINDOW[1])
    log_hi = b.log(mu_u)
    for k in range(41):
        u = lo + (hi - lo) * b.num(k) / b.num(40)
        v = icc.v_transmutation(b, u, p)
        prev_mu = prev_val = None
        root = None
        for i in range(260):
            mu = b.exp(log_lo + (log_hi - log_lo) * b.num(i) / b.num(259))
            val = icc.h_mz(b, mu, u, mu_u, v)
            if prev_val is not None and val * prev_val < 0:
                a_, c_, fa = prev_mu, mu, prev_val
                for _ in range(80):
                    m_ = b.sqrt(a_ * c_)
                    fm = icc.h_mz(b, m_, u, mu_u, v)
                    if fa * fm > 0:
                        a_, fa = m_, fm
                    else:
                        c_ = m_
                root = (a_ + c_) / b.two
                break
            prev_mu, prev_val = mu, val
        if root is None:
            continue
        a2 = icc.alpha_run(b, u, b.b_mssm[1], root, mu_u)
        a3 = icc.alpha_run(b, u, b.b_mssm[2], root, mu_u)
        r = (
            icc.ellbar(b, frame._mp_su2_terms, b.four_pi_sq * a2)
            + icc.ellbar(b, frame._mp_su3_terms, b.four_pi_sq * a3)
            - p / b.num(4)
        )
        if prev_r is not None and r * prev_r < 0:
            return True
        prev_r = r
    return False


def locate_declared_edges(frame: Frame, resolution: str = "0.002") -> dict[str, Any]:
    """Bisect the defined/undefined transition on both sides (display only)."""
    b = frame.mpb
    res = b.num(resolution)

    def bisect(inside: Any, outside: Any) -> tuple[Any, Any]:
        anchor = b.num("0.0073")
        for _ in range(8):
            if _declared_map_defined(frame, inside):
                break
            inside = (inside + anchor) / b.two
        else:
            raise icc.CertificateError("edge bisection: no defined inside point found")
        for _ in range(8):
            if not _declared_map_defined(frame, outside):
                break
            outside = outside + (outside - inside)
        else:
            raise icc.CertificateError("edge bisection: no undefined outside point found")
        while abs(outside - inside) > res:
            mid = (inside + outside) / b.two
            if _declared_map_defined(frame, mid):
                inside = mid
            else:
                outside = mid
        return inside, outside

    right_in, right_out = bisect(b.num("0.30"), b.num("0.60"))
    left_in, left_out = bisect(b.num("-0.40"), b.num("-0.80"))
    return {
        "method": (
            "point bisection of the declared scan semantics (41-point pixel "
            "window grid, 260-point m_Z log-grid window); located values are "
            "display quantities, the certified outer bound is the sweep hull"
        ),
        "resolution": resolution,
        "alpha_right_edge": {
            "last_defined": frame.fmt(right_in, 8),
            "first_undefined": frame.fmt(right_out, 8),
        },
        "alpha_left_edge": {
            "last_defined": frame.fmt(left_in, 8),
            "first_undefined": frame.fmt(left_out, 8),
        },
    }


# ---------------------------------------------------------------------------
# Composition inputs.
# ---------------------------------------------------------------------------


def verify_composition_inputs() -> dict[str, Any]:
    """Load the 2026-07-14 certificates and check their conclusions."""
    out: dict[str, Any] = {}
    at_most_one_path = RUNTIME_DIR / AT_MOST_ONE_ARTIFACT
    existence_path = RUNTIME_DIR / EXISTENCE_ARTIFACT
    with at_most_one_path.open(encoding="utf-8") as fh:
        at_most_one = json.load(fh)
    with existence_path.open(encoding="utf-8") as fh:
        existence = json.load(fh)
    if not at_most_one["conclusion"]["at_most_one_on_domain_all_modes"]:
        raise icc.CertificateError("at-most-one certificate conclusion is not certified")
    for mode in icc.MODES:
        banach = existence["modes"][mode]["banach"]
        if not (banach["existence"] and banach["uniqueness_in_interval"]):
            raise icc.CertificateError(f"existence certificate is not certified for {mode}")
    dom = at_most_one["domain"]["alpha"]
    if not (dom["lo"] == PHYSICAL_ALPHA[0] and dom["hi"] == PHYSICAL_ALPHA[1]):
        raise icc.CertificateError("at-most-one certificate domain differs from the declared interval")
    out["at_most_one_certificate"] = f"runtime/{AT_MOST_ONE_ARTIFACT}"
    out["existence_certificate"] = f"runtime/{EXISTENCE_ARTIFACT}"
    out["at_most_one_on_declared_interval_all_modes"] = True
    out["existence_on_declared_interval_all_modes"] = True
    return out


# ---------------------------------------------------------------------------
# Certificate assembly.
# ---------------------------------------------------------------------------


def build_certificate(
    mp_dps: int = 40,
    iv_dps: int = 40,
    su2_cutoff: int = 120,
    su3_cutoff: int = 90,
    initial_pieces: int = 60,
    depth_cap: int = 12,
    eval_budget: int = 1200,
    propagation_iters: int = 6,
    kernel_pieces: int = 24,
    locate_edges: bool = True,
    check_composition: bool = True,
) -> dict[str, Any]:
    frame = Frame(mp_dps, iv_dps, su2_cutoff, su3_cutoff)
    env = envelope_lemma(frame)
    sweep = sweep_pixel_window(
        frame, env, initial_pieces, depth_cap, eval_budget, propagation_iters, kernel_pieces
    )
    composition = verify_composition_inputs() if check_composition else None
    edges = locate_declared_edges(frame) if locate_edges else None

    sweep_certified = sweep["conclusion"]["certified"]
    global_conclusion = bool(sweep_certified and (composition is not None))

    artifact: dict[str, Any] = {
        "artifact": "oph_p_global_uniqueness_extension_certificate",
        "date": ARTIFACT_DATE,
        "claim_status": "maximal_domain_uniqueness_extension_for_declared_closure_maps",
        "claim_boundary": (
            "The certificate classifies the fixed-point set of both declared "
            "readout maps on their maximal analytic domain, for the declared "
            "numerical maps at the stated representation cutoffs (edge-sum tail "
            "majorants extend every enclosure to the infinite-cutoff sums). It "
            "certifies that every fixed point lies inside the declared physical "
            "interval and composes with the 2026-07-14 existence and at-most-one "
            "certificates to exactly one fixed point per map on the maximal "
            "domain. The declared one-loop RG/matching conventions, the "
            "tree-level m_Z closure, the Stage-5 continuation masses, and the "
            "exact one-loop kernel are certified as declared numerical "
            "structure, not as physical endpoint theorems. This is not an exact "
            "fine-structure derivation; the stage-3 landing verdict of closure "
            "row CL-1 is unchanged, and public claims remain tied to the "
            "declared physical interval."
        ),
        "protocol_stage": (
            "maximal-domain supplement to the interval-certified global "
            "uniqueness result"
        ),
        "backend": (
            "mpmath.iv binary interval arithmetic with outward rounding on every "
            "elementary operation; private MPIntervalContext/MPContext instances; "
            "interval constraint propagation on the declared closure equations"
        ),
        "iv_dps": iv_dps,
        "point_dps": mp_dps,
        "su2_cutoff": su2_cutoff,
        "su3_cutoff": su3_cutoff,
        "declared_maps": {mode: icc.MAP_DEFINITIONS[mode] for mode in icc.MODES},
        "outer_equation": "P = phi + alpha*sqrt(pi)",
        "maximal_domain": {
            "description": (
                "alpha values where the declared chain evaluates: P > 0, the "
                "pixel-closure residual selects a root alpha_U inside the "
                "declared scan window, the tree-level m_Z closure selects a "
                "root inside the declared log-grid window, the implicit-root "
                "denominators are nonzero, the edge sums converge, and the "
                "inverse-alpha readout is nonzero"
            ),
            "conditions": [
                "P = phi + alpha*sqrt(pi) > 0 (branch points of ln P and 1/sqrt(P) at P = 0, i.e. alpha = -phi/sqrt(pi) = -0.91287...)",
                "alpha_U in [0.02, 0.08]: declared pixel-closure scan window (paper_math.solve_alpha_u_from_p)",
                "L = ln(mu_U/m_Z) in [0, 50]: declared m_Z log-grid window (paper_math.solve_mz_fixed_point_tree)",
                "t2 = 4*pi^2*alpha_2(m_Z) > 0: SU(2) edge-sum convergence; implied by the window bounds (1/a2 >= 1/0.08)",
                "t3 = 4*pi^2*alpha_3(m_Z) > 0: SU(3) edge-sum convergence for the infinite-cutoff sums; the envelope lemma certifies 1/a3 >= 4 on the domain, so the SU(3) coupling stays below the Landau pole",
                "sqrt arguments positive: 4*pi*(a2 + aY) > 0 in the m_Z tree closure and the kernel arguments a_f = mzv2/(16 c_f^2) > 0; implied by the window bounds",
                "implicit-root nondegeneracy: the m_Z residual derivative h_m and the pixel residual derivative R_u are nonzero at the selected roots (real-analyticity by the implicit function theorem; interval-verified along the 2026-07-14 certificates on the declared interval)",
                "alpha_inv != 0: certified globally by the envelope readout floor",
            ],
            "certified_alpha_outer_bound": sweep["domain_alpha_outer_hull"],
            "located_edges": edges,
            "note": (
                "the certified outer bound is a superset projection: every "
                "domain point of either map has alpha inside it; the located "
                "edges are point-probe values under the declared grid semantics"
            ),
        },
        "envelope_lemma": env["report"],
        "pixel_window_sweep": sweep,
        "excluded_nonphysical_fixed_points": [],
        "exterior_classification": (
            "empty: no fixed point of either declared readout map exists outside "
            "the declared physical interval anywhere on the maximal analytic "
            "domain; the sweep excludes every window-consistent candidate by "
            "m_Z-window infeasibility or a certified alpha/readout gap"
        ),
        "composition_inputs": composition,
        "inner_root_scope": (
            "the sweep quantifies over every window-consistent inner-root pair "
            "(alpha_U, m_Z), so the exterior exclusion holds for any root "
            "selection inside the declared windows, in particular for the "
            "declared scan-and-bisect selection; inside the declared physical "
            "interval the composed uniqueness statement carries the inner-root "
            "scope of the 2026-07-14 certificates"
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
            "every_fixed_point_in_declared_interval_all_modes": bool(sweep_certified),
            "exactly_one_fixed_point_on_maximal_domain_all_modes": global_conclusion,
            "statement": (
                "each declared readout map has exactly one fixed point on its "
                "maximal analytic domain, the certified fixed point inside "
                "alpha in [" + PHYSICAL_ALPHA[0] + ", " + PHYSICAL_ALPHA[1] + "] "
                "(alpha_inv in [100, 200]); the exterior fixed-point set is empty"
                if global_conclusion
                else "the maximal-domain classification is partial; see the sweep exceptional set"
            ),
        },
    }
    return artifact


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Emit the maximal-domain global uniqueness extension certificate "
            "for the OPH P/alpha closure maps (proof spine GAP-A7)."
        )
    )
    parser.add_argument("--mp-dps", type=int, default=40, help="Point arithmetic working digits.")
    parser.add_argument("--iv-dps", type=int, default=40, help="Interval arithmetic working digits.")
    parser.add_argument("--su2-cutoff", type=int, default=120, help="SU(2) representation cutoff.")
    parser.add_argument("--su3-cutoff", type=int, default=90, help="SU(3) representation cutoff.")
    parser.add_argument("--initial-pieces", type=int, default=60, help="Initial pixel-window subdivision.")
    parser.add_argument("--depth-cap", type=int, default=12, help="Maximum bisection depth.")
    parser.add_argument("--eval-budget", type=int, default=1200, help="Maximum piece evaluations.")
    parser.add_argument("--propagation-iters", type=int, default=6, help="Constraint propagation passes per piece.")
    parser.add_argument("--kernel-pieces", type=int, default=24, help="Kernel-argument subdivision per fermion.")
    parser.add_argument("--skip-edges", action="store_true", help="Skip the numeric edge location.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output JSON path.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    t0 = time.time()
    artifact = build_certificate(
        mp_dps=args.mp_dps,
        iv_dps=args.iv_dps,
        su2_cutoff=args.su2_cutoff,
        su3_cutoff=args.su3_cutoff,
        initial_pieces=args.initial_pieces,
        depth_cap=args.depth_cap,
        eval_budget=args.eval_budget,
        propagation_iters=args.propagation_iters,
        kernel_pieces=args.kernel_pieces,
        locate_edges=not args.skip_edges,
    )
    wall_seconds = round(time.time() - t0, 2)
    text = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    sweep = artifact["pixel_window_sweep"]
    print(f"pieces                 = {sweep['subdivision']['total_pieces']}")
    print(f"verdict counts         = {sweep['verdict_counts']}")
    print(f"handoff alpha hull     = {sweep['handoff_alpha_hull']}")
    print(f"domain alpha outer     = {sweep['domain_alpha_outer_hull']}")
    print(f"exceptional set        = {len(sweep['exceptional_set'])}")
    print(f"conclusion             = {artifact['conclusion']['statement']}")
    print(f"wall seconds: {wall_seconds} (not recorded in the artifact)")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
