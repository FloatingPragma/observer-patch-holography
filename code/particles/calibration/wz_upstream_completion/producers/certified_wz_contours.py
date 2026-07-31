#!/usr/bin/env python3
"""Certified principal-sheet W and Z zero exclusion on the synthetic fixture.

The verifier upgrades the sampled boundary diagnostic to a directed
interval certificate on declared off-axis boxes.  Every quantity is a
rectangular complex interval with outward rounding, every logarithm
and square root carries a sheet gate refusing the principal cut, and
the argument principle is applied with endpoint-ratio argument
chaining, so a reported winding is a proved root count for the
holomorphic principal-sheet function on the declared box.

Interior holomorphy is discharged by an exact rational cut-exclusion
certificate rather than a cell sweep.  On a box strictly inside the
open upper half plane, every branch-carrying subexpression of the
compiled evaluation formula is classified once and for all:

* ``log(s/mu2)``: the argument stays in the open upper half plane;
* massless-partner charts (one vanishing mass): both logarithm
  arguments have strictly negative imaginary part on the box, so they
  exclude the closed negative real axis;
* two-massive charts: the discriminant ``(s - s_plus)(s - s_minus)``
  has nonvanishing imaginary part on the box whenever the
  threshold-sum line ``Re s = m1 + m2`` misses the box real range,
  an exact Fraction comparison; the Feynman roots are never real
  because a real root of the quadratic with ``Im s > 0`` is forced
  into ``{0, 1}`` while the quadratic takes the nonzero values ``m1``
  and ``m2`` there, so every root-chart logarithm argument excludes
  its cut and every divisor excludes zero.

The boundary winding chain replaces segment hull evaluation by a
centered form: a segment image is enclosed by the point value at the
exact rational midpoint plus the interval derivative over the segment
hull times the segment offset rectangle.  Point enclosures carry
rounding widths only, so the centered enclosure width scales
quadratically with segment length and the adaptive subdivision
terminates at moderate depth.  The derivative of the compiled block
is evaluated from the closed forms ``d/dx root_term(x) =
log((x-1)/x)`` and ``2 s x + b = +/- sqrt(disc)`` for the two Feynman
roots, sheet gated through the same interval layer.

What is certified, exactly:

* directed complex-interval evaluation of the compiled one-loop
  transverse blocks, fixture-exact coefficients, with the exact
  finite ``d = 4 - 2 eps`` dimensional-prefactor correction included
  in the compiled term list, on boxes that lie strictly in the upper
  half plane, so no evaluation meets the physical cut on the real
  axis;
* interior holomorphy of the evaluation formula on the full box by
  the exact rational cut-exclusion certificate above;
* boundary exclusion: every boundary segment centered-form enclosure
  excludes zero with the declared argument-width gate, subdividing
  adaptively to the declared depth cap;
* the winding number by rigorous interval argument chaining on point
  enclosures, with the total enclosure inside the declared tolerance
  of an integer multiple of two pi;
* enclosure nesting across the declared precision ladder, recorded
  for the summed winding total and per quantity at a declared probe
  point: the inverse propagator, its derivative, and every distinct
  loop-function value; the overall status fails closed when any
  nesting comparison fails.

What is not certified and stays false in the receipt: the second-sheet
continuation through the physical cut, any pole enclosure, Laurent and
residue data, BMHV restoration, and any physical-current or unit
claim.  The certified theorem is negative and structural: the masked
one-loop inverse propagator has no principal-sheet zero in the
declared upper-half boxes, consistent with first-sheet analyticity,
and the resonance zero is a second-sheet object.  The winding-one
value of the earlier axis-crossing sampled diagnostic is recorded as a
branch-gluing artifact superseded by this certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import deque
from fractions import Fraction
from pathlib import Path
from typing import Any

from mpmath import iv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

import wz_pole_receipts as wzp  # noqa: E402
from complex_interval import CInterval, SheetError  # noqa: E402

VECTOR_PATH = ROOT / "outputs" / "fj_direct_vector_blocks.json"
OUT_PATH = ROOT / "outputs" / "certified_wz_contours.json"
CHECKPOINT_PATH = ROOT / "outputs" / "certified_wz_contours_checkpoint.json"

SCHEMA = "oph.certified_wz_contours.v2"
STATUS_CERTIFIED = "PRINCIPAL_SHEET_ZERO_EXCLUSION_CERTIFIED"
STATUS_FAILED = "CERTIFICATION_INCOMPLETE"

PRECISIONS = (128, 192, 256)
INITIAL_SEGMENTS_PER_EDGE = 8
MAX_SUBDIVISION_DEPTH = 12
ARG_WIDTH_GATE_NUM = Fraction(157, 100)
WINDING_TOLERANCE_NUM = Fraction(157, 100)
HOLOMORPHY_METHOD = "exact_rational_cut_exclusion"
BOUNDARY_METHOD = "centered_form_segment_enclosure"
PROGRESS_STRIDE = 200

BOXES = {
    "W": {
        "re": (Fraction(109, 1000), Fraction(114, 1000)),
        "im": (Fraction(1, 5000), Fraction(1, 500)),
    },
    "Z": {
        "re": (Fraction(172, 1000), Fraction(176, 1000)),
        "im": (Fraction(1, 5000), Fraction(3, 2000)),
    },
}


def sha256_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def progress(message: str) -> None:
    print(f"[certified_wz_contours] {message}", file=sys.stderr, flush=True)


def _iv_fraction(value: Fraction) -> Any:
    return iv.mpf(value.numerator) / iv.mpf(value.denominator)


def a0_interval(m2: Fraction, mu2: Fraction) -> CInterval:
    if m2 == 0:
        return CInterval.from_fraction(0)
    ratio = _iv_fraction(Fraction(m2, mu2))
    value = _iv_fraction(m2) * (iv.mpf(1) - iv.log(ratio))
    return CInterval(value, iv.mpf(0))


def a0p_interval(m2: Fraction, mu2: Fraction) -> CInterval:
    if m2 == 0:
        return CInterval.from_fraction(0)
    ratio = _iv_fraction(Fraction(m2, mu2))
    return CInterval(-iv.log(ratio), iv.mpf(0))


def _root_term(x: CInterval, gate: Any, xm1: CInterval | None = None) -> CInterval:
    """``x log((x-1)/x) - log(x-1)`` with an optional stable ``x-1``.

    Callers whose chart forms ``x-1`` by cancellation pass the
    rationalized enclosure through ``xm1``; the direct subtraction is
    the default."""

    one = CInterval.from_fraction(1)
    if xm1 is None:
        xm1 = x - one
    return x * (xm1 / x).log(gate) - xm1.log(gate)


def _root_term_slope(
    x: CInterval, gate: Any, xm1: CInterval | None = None
) -> CInterval:
    """Derivative of the root term in its own variable.

    ``d/dx [x log((x-1)/x) - log(x-1)] = log((x-1)/x)``: the algebraic
    parts ``x/(x-1) - 1 - 1/(x-1)`` cancel exactly."""

    one = CInterval.from_fraction(1)
    if xm1 is None:
        xm1 = x - one
    return (xm1 / x).log(gate)


def _feynman_roots(
    s: CInterval, m1: Fraction, m2: Fraction, gate: Any
) -> tuple[CInterval, CInterval, CInterval, CInterval]:
    """Both Feynman roots, the gated square root, and a stable ``x1-1``.

    The quadratic-formula subtraction cancels catastrophically when a
    root approaches one or zero, so both roots use the rationalized
    forms ``x1 - 1 = -2 m2/(sqrt(disc) + s - m1 + m2)`` and
    ``x2 = 2 m1/(s + m1 - m2 + sqrt(disc))``, exact identities of the
    quadratic ``s x^2 - (s + m1 - m2) x + m1``."""

    one = CInterval.from_fraction(1)
    b = -(s + CInterval.from_fraction(m1 - m2))
    disc = b * b - CInterval.from_fraction(4) * s * CInterval.from_fraction(m1)
    sq = disc.sqrt(gate)
    x1m1 = CInterval.from_fraction(-2 * m2) / (
        sq + s + CInterval.from_fraction(m2 - m1)
    )
    x1 = one + x1m1
    x2 = CInterval.from_fraction(2 * m1) / (
        s + CInterval.from_fraction(m1 - m2) + sq
    )
    return x1, x2, sq, x1m1


def b0_interval(
    s: CInterval, m1: Fraction, m2: Fraction, mu2: Fraction, gate: Any
) -> CInterval:
    """Principal-branch B0 finite part on an off-axis rectangle.

    Every logarithm and square root is sheet gated; the caller treats
    SheetError as a subdivision signal.  The box lies strictly in the
    upper half plane, so no infinitesimal prescription is needed and
    the principal branch is the declared continuation.  Massless
    arguments are the exact algebraic limits of the root form on that
    half plane: a Feynman root at zero contributes i pi, because the
    vanishing root approaches zero with negative imaginary part there,
    and a root at one contributes zero, because its paired logarithms
    cancel in the limit."""

    one = CInterval.from_fraction(1)
    two = CInterval.from_fraction(2)
    mu = CInterval.from_fraction(mu2)
    total = two - (s / mu).log(gate)
    i_pi = CInterval(iv.mpf(0), iv.pi)
    if m1 == 0 and m2 == 0:
        return total + i_pi
    if m1 == 0:
        x = one - CInterval.from_fraction(m2) / s
        return total + i_pi + _root_term(x, gate)
    if m2 == 0:
        x = CInterval.from_fraction(m1) / s
        return total + _root_term(x, gate)
    x1, x2, _, x1m1 = _feynman_roots(s, m1, m2, gate)
    return total + _root_term(x1, gate, xm1=x1m1) + _root_term(x2, gate)


def b0p_interval(
    s: CInterval, m1: Fraction, m2: Fraction, mu2: Fraction, gate: Any
) -> CInterval:
    """Derivative of the principal-branch B0 finite part in ``s``.

    The chain rule composes ``d/dx root_term = log((x-1)/x)`` with the
    implicit root slopes.  For the two-massive chart the quadratic
    ``s x^2 - (s + m1 - m2) x + m1`` gives ``x' = x(1-x)/(2 s x + b)``
    and ``2 s x + b`` equals ``+sqrt(disc)`` on the first root and
    ``-sqrt(disc)`` on the second, so no additional branch decision
    enters.  Massless charts differentiate their explicit roots."""

    one = CInterval.from_fraction(1)
    minus_inv_s = -(one / s)
    if m1 == 0 and m2 == 0:
        return minus_inv_s
    if m1 == 0:
        x = one - CInterval.from_fraction(m2) / s
        x_slope = CInterval.from_fraction(m2) / (s * s)
        return minus_inv_s + _root_term_slope(x, gate) * x_slope
    if m2 == 0:
        x = CInterval.from_fraction(m1) / s
        x_slope = -(CInterval.from_fraction(m1) / (s * s))
        return minus_inv_s + _root_term_slope(x, gate) * x_slope
    x1, x2, sq, x1m1 = _feynman_roots(s, m1, m2, gate)
    x1_slope = x1 * (-x1m1) / sq
    x2_slope = -(x2 * (one - x2) / sq)
    return (
        minus_inv_s
        + _root_term_slope(x1, gate, xm1=x1m1) * x1_slope
        + _root_term_slope(x2, gate) * x2_slope
    )


def _poly_derivative(
    coefficients: list[tuple[int, Fraction]],
) -> list[tuple[int, Fraction]]:
    return [
        (power - 1, fraction * power)
        for power, fraction in coefficients
        if power > 0
    ]


class IntervalEvaluator:
    """Interval evaluation of one compiled transverse block."""

    def __init__(self, compiled: list[dict[str, Any]], gate: Any) -> None:
        self.compiled = compiled
        self.gate = gate
        mu2 = wzp.FIXTURE["mu_ren2"]
        self.mu2 = mu2
        self.loop_cache: dict[tuple, CInterval] = {}
        self.derived = [
            {
                "num": _poly_derivative(term["num"]),
                "den": _poly_derivative(term["den"]),
            }
            for term in compiled
        ]

    def _loop(self, head: str, args: tuple, s_key: tuple, s: CInterval) -> CInterval:
        key = (head, args) + (s_key if head == "B0" else ())
        cached = self.loop_cache.get(key)
        if cached is not None:
            return cached
        if head == "A0":
            value = a0_interval(args[0], self.mu2)
        elif head == "A0p":
            value = a0p_interval(args[0], self.mu2)
        elif head == "B0":
            value = b0_interval(s, args[0], args[1], self.mu2, self.gate)
        else:
            raise SheetError(f"unexpected loop head {head}")
        self.loop_cache[key] = value
        return value

    def _loop_slope(
        self, head: str, args: tuple, s_key: tuple, s: CInterval
    ) -> CInterval | None:
        if head in ("A0", "A0p"):
            return None
        if head == "B0":
            key = ("B0p", args) + s_key
            cached = self.loop_cache.get(key)
            if cached is not None:
                return cached
            value = b0p_interval(s, args[0], args[1], self.mu2, self.gate)
            self.loop_cache[key] = value
            return value
        raise SheetError(f"unexpected loop head {head}")

    def _poly(self, coefficients: list, s: CInterval) -> CInterval:
        total = CInterval.from_fraction(0)
        for power, fraction in coefficients:
            term = CInterval.from_fraction(fraction) * s.pow_int(int(power))
            total = total + term
        return total

    def transverse(self, s: CInterval, s_key: tuple) -> CInterval:
        total = CInterval.from_fraction(0)
        for term in self.compiled:
            denominator = self._poly(term["den"], s)
            coefficient = self._poly(term["num"], s) / denominator
            if term["head"] is None:
                total = total + coefficient
            else:
                total = total + coefficient * self._loop(
                    term["head"], term["args"], s_key, s
                )
        pi_sq = iv.pi * iv.pi
        loop_factor = CInterval(iv.mpf(1) / (iv.mpf(16) * pi_sq), iv.mpf(0))
        return total * loop_factor

    def transverse_derivative(self, s: CInterval, s_key: tuple) -> CInterval:
        total = CInterval.from_fraction(0)
        for term, derived in zip(self.compiled, self.derived):
            den = self._poly(term["den"], s)
            num = self._poly(term["num"], s)
            num_d = self._poly(derived["num"], s)
            den_d = self._poly(derived["den"], s)
            coefficient = num / den
            coefficient_slope = (num_d * den - num * den_d) / (den * den)
            if term["head"] is None:
                total = total + coefficient_slope
                continue
            loop_value = self._loop(term["head"], term["args"], s_key, s)
            total = total + coefficient_slope * loop_value
            loop_slope = self._loop_slope(term["head"], term["args"], s_key, s)
            if loop_slope is not None:
                total = total + coefficient * loop_slope
        pi_sq = iv.pi * iv.pi
        loop_factor = CInterval(iv.mpf(1) / (iv.mpf(16) * pi_sq), iv.mpf(0))
        return total * loop_factor

    def inverse_propagator(
        self, s: CInterval, s_key: tuple, tree_mass: Fraction
    ) -> CInterval:
        return s - CInterval.from_fraction(tree_mass) - self.transverse(s, s_key)

    def inverse_propagator_derivative(
        self, s: CInterval, s_key: tuple, tree_mass: Fraction
    ) -> CInterval:
        one = CInterval.from_fraction(1)
        return one - self.transverse_derivative(s, s_key)


def _segment_hull(
    start: tuple[Fraction, Fraction], end: tuple[Fraction, Fraction]
) -> CInterval:
    return CInterval.box(
        min(start[0], end[0]),
        max(start[0], end[0]),
        min(start[1], end[1]),
        max(start[1], end[1]),
    )


def _boundary_segments(
    box: dict[str, tuple[Fraction, Fraction]], per_edge: int
) -> list[tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]]:
    re_lo, re_hi = box["re"]
    im_lo, im_hi = box["im"]
    corners = [
        (re_lo, im_lo),
        (re_hi, im_lo),
        (re_hi, im_hi),
        (re_lo, im_hi),
    ]
    segments = []
    for index in range(4):
        start = corners[index]
        end = corners[(index + 1) % 4]
        for step in range(per_edge):
            t0 = Fraction(step, per_edge)
            t1 = Fraction(step + 1, per_edge)
            p0 = (
                start[0] + (end[0] - start[0]) * t0,
                start[1] + (end[1] - start[1]) * t0,
            )
            p1 = (
                start[0] + (end[0] - start[0]) * t1,
                start[1] + (end[1] - start[1]) * t1,
            )
            segments.append((p0, p1))
    return segments


def certify_interior_holomorphy(
    compiled: list[dict[str, Any]],
    box: dict[str, tuple[Fraction, Fraction]],
) -> dict[str, Any]:
    """Exact rational cut exclusion for the evaluation formula.

    The certificate discharges holomorphy of the compiled formula on
    the whole box from a finite list of exact comparisons; no interval
    subdivision enters.  Each branch-carrying subexpression is
    classified by chart, and the recorded facts are the ones the
    docstring derivation consumes: the box sits in the open upper half
    plane and to the right of zero, massless-partner logarithm
    arguments have strictly negative imaginary part there, and each
    two-massive discriminant misses the closed negative real axis
    because the threshold-sum line misses the box real range while the
    Feynman roots stay off the real axis.  Denominator polynomials are
    checked to exclude zero on the box hull."""

    re_lo, re_hi = box["re"]
    im_lo, im_hi = box["im"]
    base_facts = {
        "box_in_open_upper_half_plane": im_lo > 0,
        "box_real_part_positive": re_lo > 0,
    }
    rows: list[dict[str, Any]] = []
    denominators_ok = True
    hull = CInterval.box(re_lo, re_hi, im_lo, im_hi)
    seen_dens: set = set()
    for term in compiled:
        den_key = tuple(term["den"])
        if den_key not in seen_dens:
            seen_dens.add(den_key)
            probe = IntervalEvaluator([], None)
            den_value = probe._poly(term["den"], hull)
            if den_value.contains_zero():
                denominators_ok = False
                rows.append(
                    {
                        "kind": "denominator",
                        "den": [[p, str(f)] for p, f in term["den"]],
                        "excludes_zero_on_box": False,
                    }
                )
    seen_pairs: set = set()
    pairs_ok = True
    for term in compiled:
        if term["head"] != "B0" or term["args"] in seen_pairs:
            continue
        seen_pairs.add(term["args"])
        m1, m2 = term["args"]
        if m1 == 0 and m2 == 0:
            row = {
                "m1": str(m1),
                "m2": str(m2),
                "chart": "both_masses_zero",
                "certificate": "log(s/mu2) only; box in open upper half plane",
                "holomorphic": bool(im_lo > 0),
            }
        elif m1 == 0 or m2 == 0:
            row = {
                "m1": str(m1),
                "m2": str(m2),
                "chart": "one_mass_zero",
                "certificate": (
                    "explicit root chart; both logarithm arguments have "
                    "strictly negative imaginary part on the box"
                ),
                "holomorphic": bool(im_lo > 0 and re_lo > 0),
            }
        else:
            threshold_sum = m1 + m2
            outside = threshold_sum < re_lo or threshold_sum > re_hi
            row = {
                "m1": str(m1),
                "m2": str(m2),
                "chart": "two_massive",
                "threshold_sum": str(threshold_sum),
                "threshold_sum_outside_box_re_range": bool(outside),
                "certificate": (
                    "discriminant imaginary part nonvanishing off the "
                    "threshold-sum line; Feynman roots never real because "
                    "the quadratic equals m1 at x=0 and m2 at x=1"
                ),
                "holomorphic": bool(outside and im_lo > 0 and re_lo > 0),
            }
        pairs_ok = pairs_ok and row["holomorphic"]
        rows.append(row)
    holomorphic = bool(
        base_facts["box_in_open_upper_half_plane"]
        and base_facts["box_real_part_positive"]
        and denominators_ok
        and pairs_ok
    )
    return {
        "method": HOLOMORPHY_METHOD,
        "base_facts": base_facts,
        "denominators_exclude_zero_on_box": denominators_ok,
        "loop_charts": rows,
        "holomorphic_on_box": holomorphic,
    }


def certify_winding(
    evaluator: Any,
    box: dict[str, tuple[Fraction, Fraction]],
    tree_mass: Fraction,
    gate: Any,
    fixed_segments: list | None = None,
) -> dict[str, Any]:
    """Rigorous argument-principle winding with endpoint-ratio chaining.

    Each boundary segment is enclosed by the centered form: the point
    value at the exact rational midpoint plus the interval derivative
    over the segment hull times the offset rectangle.  The enclosure
    must exclude zero and, rotated by the conjugate midpoint value so
    the test chart faces the positive real axis, show argument width
    below the gate; that confines the image of the segment to a cone
    narrower than pi, so the argument cannot wrap within the segment.
    The rotation leaves angular width invariant while keeping the
    two-argument arctangent away from its own chart cut, so a segment
    whose image crosses the negative real axis is certified in place
    rather than subdivided against the chart.  The argument increment
    of the segment is the principal argument of the endpoint ratio,
    evaluated as ``f(end) * conj(f(start))`` on point rectangles whose
    enclosure widths carry rounding only; cone confinement bounds the
    within-segment increment by the gate, strictly below pi, so the
    principal value is the unique admissible lift and no branch
    decision enters.  The winding is the summed increment enclosure
    divided by two pi, certified when the residual enclosure lies
    inside the declared tolerance.

    With ``fixed_segments`` the adaptive search is replaced by a
    replay of the given partition: every segment must pass the gates
    at the working precision without subdivision, which ties the
    precision ladder to one declared partition so per-segment
    enclosure nesting is a well-defined comparison."""

    point_values: dict[tuple, CInterval] = {}

    def point_value(point: tuple) -> CInterval:
        cached = point_values.get(point)
        if cached is not None:
            return cached
        rect = CInterval.box(point[0], point[0], point[1], point[1])
        evaluator.loop_cache.clear()
        value = evaluator.inverse_propagator(rect, ("pt", str(point)), tree_mass)
        point_values[point] = value
        return value

    def conjugate(value: CInterval) -> CInterval:
        return CInterval(value.re, -value.im)

    def segment_enclosure(start: tuple, end: tuple) -> CInterval | None:
        mid = ((start[0] + end[0]) / 2, (start[1] + end[1]) / 2)
        hull = _segment_hull(start, end)
        try:
            center = point_value(mid)
            evaluator.loop_cache.clear()
            slope = evaluator.inverse_propagator_derivative(
                hull, ("seg", str(start), str(end)), tree_mass
            )
            offset = CInterval.box(
                min(start[0], end[0]) - mid[0],
                max(start[0], end[0]) - mid[0],
                min(start[1], end[1]) - mid[1],
                max(start[1], end[1]) - mid[1],
            )
            image = center + slope * offset
            rotated = image * conjugate(center)
            ok = (
                (not image.contains_zero())
                and (not rotated.contains_zero())
                and rotated.arg().delta <= gate
            )
        except SheetError:
            return None
        return image if ok else None

    certified: list[tuple] = []
    enclosures: list[CInterval] = []
    max_depth_used = 0
    processed = 0
    if fixed_segments is not None:
        for index, (start, end) in enumerate(fixed_segments):
            processed += 1
            if processed % PROGRESS_STRIDE == 0:
                progress(
                    f"  winding replay: {processed}/{len(fixed_segments)}"
                )
            image = segment_enclosure(start, end)
            if image is None:
                return {
                    "certified": False,
                    "method": BOUNDARY_METHOD,
                    "partition": "replayed_base_partition",
                    "reason": f"fixed segment {index} fails the gates",
                    "segments": len(certified),
                }
            certified.append((start, end))
            enclosures.append(image)
    else:
        worklist = deque(
            (segment, 0)
            for segment in _boundary_segments(box, INITIAL_SEGMENTS_PER_EDGE)
        )
        while worklist:
            (start, end), depth = worklist.popleft()
            processed += 1
            if processed % PROGRESS_STRIDE == 0:
                progress(
                    f"  winding: {processed} segment evaluations, "
                    f"{len(certified)} certified, depth<= {max_depth_used}"
                )
            image = segment_enclosure(start, end)
            if image is not None:
                certified.append((start, end))
                enclosures.append(image)
                continue
            if depth >= MAX_SUBDIVISION_DEPTH:
                return {
                    "certified": False,
                    "method": BOUNDARY_METHOD,
                    "partition": "adaptive",
                    "reason": "subdivision depth cap reached",
                    "segments": len(certified),
                }
            mid = ((start[0] + end[0]) / 2, (start[1] + end[1]) / 2)
            worklist.appendleft(((mid, end), depth + 1))
            worklist.appendleft(((start, mid), depth + 1))
            max_depth_used = max(max_depth_used, depth + 1)

    two_pi = iv.pi * iv.mpf(2)
    total = iv.mpf(0)
    for start, end in certified:
        try:
            ratio = point_value(end) * conjugate(point_value(start))
            increment = ratio.arg()
        except SheetError:
            return {
                "certified": False,
                "method": BOUNDARY_METHOD,
                "reason": "endpoint ratio enclosure meets the chart cut",
                "segments": len(certified),
            }
        if abs(increment.a) > gate or abs(increment.b) > gate:
            return {
                "certified": False,
                "method": BOUNDARY_METHOD,
                "reason": "endpoint ratio increment outside the cone gate",
                "segments": len(certified),
            }
        total = total + increment
    winding = int(
        round((float(total.a) + float(total.b)) / 2 / float(two_pi.b))
    )
    tolerance = _iv_fraction(WINDING_TOLERANCE_NUM)
    residual = total - two_pi * iv.mpf(winding)
    within = bool(residual.a > -tolerance.a and residual.b < tolerance.a)
    return {
        "certified": within,
        "method": BOUNDARY_METHOD,
        "partition": (
            "adaptive" if fixed_segments is None else "replayed_base_partition"
        ),
        "winding": winding,
        "segments": len(certified),
        "max_depth_used": max_depth_used,
        "total_variation_interval": [str(total.a), str(total.b)],
        "reason": None if within else "total variation outside tolerance",
        "_raw_total": total,
        "_segments": certified,
        "_segment_enclosures": enclosures,
    }


def probe_quantities(
    evaluator: IntervalEvaluator,
    box: dict[str, tuple[Fraction, Fraction]],
    tree_mass: Fraction,
) -> dict[str, CInterval]:
    """Evaluate the ladder probe set at the exact box center.

    The probe set carries one enclosure per certified quantity kind:
    the inverse propagator, its derivative, and every distinct loop
    value the evaluation touches, so precision-ladder nesting is
    recorded per quantity rather than only for the summed total."""

    center = (
        (box["re"][0] + box["re"][1]) / 2,
        (box["im"][0] + box["im"][1]) / 2,
    )
    rect = CInterval.box(center[0], center[0], center[1], center[1])
    evaluator.loop_cache.clear()
    quantities = {
        "inverse_propagator": evaluator.inverse_propagator(
            rect, ("probe",), tree_mass
        ),
        "inverse_propagator_derivative": evaluator.inverse_propagator_derivative(
            rect, ("probe",), tree_mass
        ),
    }
    for key, value in evaluator.loop_cache.items():
        head, args = key[0], key[1]
        label = f"{head}({','.join(str(a) for a in args)})"
        quantities[label] = value
    return quantities


def certify_particle(
    name: str,
    compiled: list[dict[str, Any]],
    tree_mass: Fraction,
    precision: int,
    fixed_segments: list | None = None,
) -> dict[str, Any]:
    iv.prec = precision
    tight_gate = _iv_fraction(ARG_WIDTH_GATE_NUM).a
    box = BOXES[name]
    progress(f"{name} @ {precision} bits: interior certificate")
    holomorphy = certify_interior_holomorphy(compiled, box)
    progress(
        f"{name} @ {precision} bits: interior "
        f"holomorphic_on_box={holomorphy['holomorphic_on_box']}; "
        "boundary winding"
    )
    evaluator = IntervalEvaluator(compiled, tight_gate)
    winding = certify_winding(
        evaluator, box, tree_mass, tight_gate, fixed_segments=fixed_segments
    )
    progress(
        f"{name} @ {precision} bits: winding certified="
        f"{winding['certified']} winding={winding.get('winding')} "
        f"segments={winding.get('segments')}"
    )
    probes = probe_quantities(evaluator, box, tree_mass)
    certified = bool(
        holomorphy["holomorphic_on_box"]
        and winding["certified"]
        and winding.get("winding") == 0
    )
    return {
        "particle": name,
        "precision_bits": precision,
        "box": {
            "re": [str(box["re"][0]), str(box["re"][1])],
            "im": [str(box["im"][0]), str(box["im"][1])],
            "half_plane": "upper",
        },
        "tree_mass_sq": str(tree_mass),
        "interior_holomorphy": holomorphy,
        "boundary_winding": winding,
        "_probes": probes,
        "zero_exclusion_certified": certified,
        "reading": (
            "the masked one-loop inverse propagator has no zero in the "
            "declared upper-half box on the principal sheet, consistent "
            "with first-sheet analyticity; the resonance zero belongs to "
            "the second sheet"
            if certified
            else "not certified"
        ),
    }


def _write_checkpoint(rows: dict[str, Any]) -> None:
    serializable = {
        key: {k: v for k, v in row.items() if k != "_raw_total"}
        for key, row in rows.items()
    }
    CHECKPOINT_PATH.write_text(
        json.dumps(serializable, sort_keys=True, indent=1, default=str) + "\n",
        encoding="utf-8",
    )


def build_receipt() -> dict[str, Any]:
    vector_raw = VECTOR_PATH.read_bytes()
    vector = json.loads(vector_raw.decode("utf-8"))
    g1 = wzp.FIXTURE["g1"]
    g2 = wzp.FIXTURE["g2"]
    v = wzp.FIXTURE["v"]
    masses = {
        "W": Fraction(g2**2 * v**2, 4),
        "Z": Fraction((g1**2 + g2**2) * v**2, 4),
    }
    progress("compiling transverse blocks with prefactor corrections")
    expressions = {
        "W": wzp.block_transverse(vector, "WpWm"),
        "Z": wzp.block_transverse(vector, "ZZ"),
    }
    corrections = {
        name: wzp.dimensional_prefactor_finite_correction(expression)
        for name, expression in expressions.items()
    }
    compiled = {
        name: wzp.compile_block(expressions[name])
        + wzp.compile_block(corrections[name])
        for name in expressions
    }

    def interval_nested(outer: Any, inner: Any) -> bool:
        return bool(inner.a >= outer.a and inner.b <= outer.b)

    def cinterval_nested(outer: CInterval, inner: CInterval) -> bool:
        return interval_nested(outer.re, inner.re) and interval_nested(
            outer.im, inner.im
        )

    results: dict[str, Any] = {}
    nesting: dict[str, Any] = {}
    checkpoint_rows: dict[str, Any] = {}
    all_certified = True
    all_nested = True
    for name in ("W", "Z"):
        rows = {}
        raw_totals = []
        probe_ladder: list[dict[str, CInterval]] = []
        enclosure_ladder: list[list[CInterval]] = []
        base_segments: list | None = None
        for precision in PRECISIONS:
            row = certify_particle(
                name,
                compiled[name],
                masses[name],
                precision,
                fixed_segments=base_segments,
            )
            winding_row = row["boundary_winding"]
            segments = winding_row.pop("_segments", None)
            enclosures = winding_row.pop("_segment_enclosures", None)
            if base_segments is None and segments is not None:
                base_segments = segments
            if enclosures is not None:
                enclosure_ladder.append(enclosures)
            checkpoint_rows[f"{name}:{precision}"] = {
                k: v
                for k, v in winding_row.items()
                if not k.startswith("_")
            } | {"zero_exclusion_certified": row["zero_exclusion_certified"]}
            _write_checkpoint(checkpoint_rows)
            raw = winding_row.pop("_raw_total", None)
            probe_ladder.append(row.pop("_probes"))
            rows[str(precision)] = row
            all_certified = all_certified and row["zero_exclusion_certified"]
            if raw is not None:
                raw_totals.append(raw)
        segment_nesting = {
            "segments": len(enclosure_ladder[0]) if enclosure_ladder else 0,
            "ladders_compared": max(len(enclosure_ladder) - 1, 0),
            "all_nested": bool(
                len(enclosure_ladder) == len(PRECISIONS)
                and all(
                    len(enclosure_ladder[k]) == len(enclosure_ladder[0])
                    for k in range(len(enclosure_ladder))
                )
                and all(
                    cinterval_nested(
                        enclosure_ladder[k][j], enclosure_ladder[k + 1][j]
                    )
                    for k in range(len(enclosure_ladder) - 1)
                    for j in range(len(enclosure_ladder[0]))
                )
            ),
        }
        total_nested = (
            all(
                interval_nested(raw_totals[k], raw_totals[k + 1])
                for k in range(len(raw_totals) - 1)
            )
            if len(raw_totals) == len(PRECISIONS)
            else False
        )
        probe_names = sorted(probe_ladder[0].keys())
        per_quantity = {
            label: all(
                label in probe_ladder[k + 1]
                and cinterval_nested(
                    probe_ladder[k][label], probe_ladder[k + 1][label]
                )
                for k in range(len(probe_ladder) - 1)
            )
            for label in probe_names
        }
        box = BOXES[name]
        nesting[name] = {
            "enclosures_nested_with_precision": bool(total_nested),
            "comparison": "exact endpoint comparison of the raw enclosures",
            "probe_point": [
                str((box["re"][0] + box["re"][1]) / 2),
                str((box["im"][0] + box["im"][1]) / 2),
            ],
            "per_quantity_probe_nesting": per_quantity,
            "per_segment_enclosure_nesting": segment_nesting,
        }
        all_nested = (
            all_nested
            and total_nested
            and all(per_quantity.values())
            and segment_nesting["all_nested"]
        )
        results[name] = rows
    all_certified = all_certified and all_nested

    payload = {
        "schema": SCHEMA,
        "type": "EXTERNAL_SM_EFT_VALIDATION",
        "pins": {
            "vector_blocks_sha256": sha256_bytes(vector_raw),
            "diagnostic_module_sha256": sha256_bytes(
                (ROOT / "producers" / "wz_pole_receipts.py").read_bytes()
            ),
            "interval_module_sha256": sha256_bytes(
                (ROOT / "producers" / "complex_interval.py").read_bytes()
            ),
            "verifier_module_sha256": sha256_bytes(Path(__file__).read_bytes()),
        },
        "fixture": {k: str(f) for k, f in wzp.FIXTURE.items()},
        "dimensional_prefactor_finite_correction": {
            name: str(corrections[name]) for name in corrections
        },
        "serialized_gates": {
            "precisions_bits": list(PRECISIONS),
            "initial_segments_per_edge": INITIAL_SEGMENTS_PER_EDGE,
            "max_subdivision_depth": MAX_SUBDIVISION_DEPTH,
            "holomorphy_method": HOLOMORPHY_METHOD,
            "boundary_method": BOUNDARY_METHOD,
            "arg_width_gate": str(ARG_WIDTH_GATE_NUM),
            "winding_tolerance": str(WINDING_TOLERANCE_NUM),
        },
        "sheet_statement": (
            "the boxes lie strictly in the upper half plane; the principal "
            "branch is the declared continuation and every logarithm and "
            "square root is cut-gated per evaluation; the second-sheet "
            "continuation through the physical cut is not constructed"
        ),
        "results": results,
        "precision_nesting": nesting,
        "promotion": {
            "complex_ball_certified": all_certified,
            "sheet_certified_on_declared_boxes": all_certified,
            "principal_sheet_zero_exclusion_certified": all_certified,
            "root_count_certified_on_declared_boxes": all_certified,
            "pole_enclosure_certified": False,
            "second_sheet_certified": False,
            "laurent_residue_certified": False,
            "bmhv_restoration_certified": False,
            "physical_current_claim": False,
            "oph_native": False,
            "unit_claim": False,
        },
        "artifact_finding": (
            "the earlier sampled diagnostic reported winding one on boxes "
            "crossing the real axis; that value is a branch-gluing "
            "artifact of the discontinuous epsilon evaluation across the "
            "physical cut and is superseded by this certificate: the "
            "principal-sheet function has no zero in the declared "
            "off-axis boxes, and the resonance zero is a second-sheet "
            "object whose chart is not constructed here"
        ),
        "status": STATUS_CERTIFIED if all_certified else STATUS_FAILED,
    }
    payload["receipt_sha256"] = sha256_bytes(
        canonical_json(
            {k: v for k, v in payload.items() if k != "receipt_sha256"}
        ).encode("utf-8")
    )
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    payload = build_receipt()
    if args.verify:
        stored = json.loads(OUT_PATH.read_text(encoding="utf-8"))
        if stored != payload:
            print("CERTIFIED_CONTOUR_DRIFT", file=sys.stderr)
            return 1
        print("CERTIFIED_CONTOUR_VERIFIED")
        return 0
    OUT_PATH.write_text(
        json.dumps(payload, sort_keys=True, indent=1) + "\n", encoding="utf-8"
    )
    if CHECKPOINT_PATH.exists():
        CHECKPOINT_PATH.unlink()
    print(json.dumps({"status": payload["status"]}))
    return 0 if payload["status"] == STATUS_CERTIFIED else 1


if __name__ == "__main__":
    raise SystemExit(main())
