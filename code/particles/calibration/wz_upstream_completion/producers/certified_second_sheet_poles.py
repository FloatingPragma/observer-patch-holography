#!/usr/bin/env python3
"""Certified second-sheet W and Z pole enclosures on the synthetic fixture.

The producer completes the pole side of the contour row: a declared
continuation path through the physical cut, one simple root in a
certified contour, Laurent residue data, and a derivative ball
excluding zero at the pole.

Declared continuation path.  The principal-sheet formula is continued
through a declared cut window between exact rational thresholds.  Per
open channel the continuation adds the realized-jump analytic
function of the principal-branch formula:

* two-massive channels add ``2 pi i sqrt(lambda(s))/s`` with
  ``lambda(s) = s^2 - 2(m1+m2)s + (m1-m2)^2`` under the same sheet
  gate as the Feynman-root square root;
* massless-partner channels add the rational function
  ``-2 pi i m/s``, the exact realized jump of the explicit root
  chart, with the both-massless channel adding zero;

and closed channels stay on the principal branch.  The identification
of the glued formula with the continuation through the window is a
declared convention backed by recorded boundary-value consistency
probes on a delta ladder at two exact rational window points per open
channel; the analyticity of the glued formula on the lower-half
corridor connecting the window to the pole box is machine-checked by
the same exact rational cut-exclusion certificate as the
zero-exclusion receipt, mirrored to the open lower half plane.

Certified on each declared pole box, per precision:

* interior holomorphy of the continued formula on the pole box and on
  the declared corridor box by exact rational cut exclusion;
* boundary winding one by the centered-form segment enclosures and
  endpoint-ratio argument chaining imported from the zero-exclusion
  verifier: winding one counts one root with multiplicity one, so the
  root is simple;
* an interval-Newton contraction certifying a unique root inside a
  declared refinement box and returning a tight pole ball;
* a derivative ball over the refinement box excluding zero: the
  Laurent denominator of the propagator residue;
* the scalar Laurent residue ball ``1/G'`` over the refinement box;
* per-quantity probe nesting across the declared precision ladder,
  failing closed on any comparison failure.

What is not certified and stays false in the receipt: BMHV
restoration, physical-current normalization, any unit claim, and any
OPH-native statement.  The receipt consumes the transverse blocks
with the exact finite dimensional-prefactor correction and records
the zero-exclusion receipt as the principal-sheet counterpart.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from fractions import Fraction
from math import isqrt
from pathlib import Path
from typing import Any

from mpmath import iv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "producers"))

import wz_pole_receipts as wzp  # noqa: E402
import certified_wz_contours as cwc  # noqa: E402
from complex_interval import CInterval, SheetError  # noqa: E402

VECTOR_PATH = ROOT / "outputs" / "fj_direct_vector_blocks.json"
ZERO_EXCLUSION_PATH = ROOT / "outputs" / "certified_wz_contours.json"
OUT_PATH = ROOT / "outputs" / "certified_second_sheet_poles.json"
CHECKPOINT_PATH = ROOT / "outputs" / "certified_second_sheet_poles_checkpoint.json"

SCHEMA = "oph.certified_second_sheet_poles.v1"
STATUS_CERTIFIED = "SECOND_SHEET_POLE_CERTIFIED_ON_DECLARED_CONTINUATION"
STATUS_FAILED = "CERTIFICATION_INCOMPLETE"

PRECISIONS = cwc.PRECISIONS
NEWTON_STEPS = 4
PROBE_DELTAS = (Fraction(1, 10**6), Fraction(1, 10**8))
PROBE_RESIDUAL_GATE = Fraction(1, 1000)

# Declared boxes: the pole box carries the winding-one certificate, the
# seed box starts the interval-Newton refinement, and the corridor box
# connects the cut window to the pole box inside the open lower half
# plane.  All endpoints are exact rationals.
GEOMETRY = {
    "W": {
        "pole_box": {
            "re": (Fraction(1119, 10000), Fraction(1121, 10000)),
            "im": (Fraction(-12, 100000), Fraction(-5, 100000)),
        },
        "newton_seed": {
            "re": (Fraction(1119625, 10**7), Fraction(1119632, 10**7)),
            "im": (Fraction(-882, 10**7), Fraction(-878, 10**7)),
        },
        "corridor_box": {
            "re": (Fraction(9, 100), Fraction(1130, 10000)),
            "im": (Fraction(-13, 100000), Fraction(-1, 10**9)),
        },
    },
    "Z": {
        "pole_box": {
            "re": (Fraction(1742, 10000), Fraction(1746, 10000)),
            "im": (Fraction(-15, 10000), Fraction(-11, 10000)),
        },
        "newton_seed": {
            "re": (Fraction(1744045, 10**7), Fraction(1744055, 10**7)),
            "im": (Fraction(-12735, 10**7), Fraction(-12725, 10**7)),
        },
        "corridor_box": {
            "re": (Fraction(16, 100), Fraction(1746, 10000)),
            "im": (Fraction(-16, 10000), Fraction(-1, 10**9)),
        },
    },
}


def sha256_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def rational_sqrt(value: Fraction) -> Fraction:
    """Exact rational square root; fails closed on non-square input."""

    pn, qn = value.numerator, value.denominator
    a, b = isqrt(pn), isqrt(qn)
    if a * a != pn or b * b != qn:
        raise SystemExit(f"fixture mass {value} has no rational square root")
    return Fraction(a, b)


class SecondSheetEvaluator(cwc.IntervalEvaluator):
    """Continued transverse evaluation on the declared second sheet.

    Open channels receive the realized-jump addition on top of the
    principal-branch loop value; closed channels are untouched.  The
    additions and their derivatives run through the same sheet-gated
    interval layer, so a refusing enclosure stays a subdivision
    signal."""

    def __init__(
        self,
        compiled: list[dict[str, Any]],
        gate: Any,
        open_pairs: set,
    ) -> None:
        super().__init__(compiled, gate)
        self.open_pairs = open_pairs

    def _two_pi_i(self) -> CInterval:
        return CInterval(iv.mpf(0), 2 * iv.pi)

    def addition(self, args: tuple, s: CInterval) -> CInterval:
        m1, m2 = args
        if m1 == 0 or m2 == 0:
            m = m1 if m2 == 0 else m2
            return self._two_pi_i() * CInterval.from_fraction(-m) / s
        sq = cwc._feynman_roots(s, m1, m2, self.gate)[2]
        return self._two_pi_i() * sq / s

    def addition_slope(self, args: tuple, s: CInterval) -> CInterval:
        m1, m2 = args
        if m1 == 0 or m2 == 0:
            m = m1 if m2 == 0 else m2
            return self._two_pi_i() * CInterval.from_fraction(m) / (s * s)
        sq = cwc._feynman_roots(s, m1, m2, self.gate)[2]
        sq_slope = (s - CInterval.from_fraction(m1 + m2)) / sq
        return self._two_pi_i() * (sq_slope * s - sq) / (s * s)

    def _loop(self, head: str, args: tuple, s_key: tuple, s: CInterval) -> CInterval:
        value = super()._loop(head, args, s_key, s)
        if head == "B0" and args in self.open_pairs:
            key = ("B0II_add", args) + s_key
            cached = self.loop_cache.get(key)
            if cached is None:
                cached = self.addition(args, s)
                self.loop_cache[key] = cached
            value = value + cached
        return value

    def _loop_slope(
        self, head: str, args: tuple, s_key: tuple, s: CInterval
    ) -> CInterval | None:
        slope = super()._loop_slope(head, args, s_key, s)
        if head == "B0" and args in self.open_pairs:
            slope = slope + self.addition_slope(args, s)
        return slope


def classify_channels(
    compiled: list[dict[str, Any]], tree_mass: Fraction
) -> dict[str, Any]:
    """Exact threshold classification and the declared cut window."""

    pairs = sorted({term["args"] for term in compiled if term["head"] == "B0"})
    rows = []
    open_pairs = set()
    open_ths = []
    closed_ths = []
    for m1, m2 in pairs:
        threshold = (rational_sqrt(m1) + rational_sqrt(m2)) ** 2
        is_open = threshold < tree_mass
        if is_open:
            open_pairs.add((m1, m2))
            open_ths.append(threshold)
        else:
            closed_ths.append(threshold)
        if m1 == 0 and m2 == 0:
            addition = "0"
        elif m1 == 0 or m2 == 0:
            m = m1 if m2 == 0 else m2
            addition = f"-2*pi*i*({m})/s"
        else:
            addition = "2*pi*i*sqrt(lambda(s))/s"
        rows.append(
            {
                "m1": str(m1),
                "m2": str(m2),
                "threshold": str(threshold),
                "open": bool(is_open),
                "addition": addition if is_open else "principal branch",
            }
        )
    window = (max(open_ths), min(closed_ths + [tree_mass]))
    return {
        "channels": rows,
        "open_pairs": open_pairs,
        "window": window,
    }


def consistency_probes(
    open_pairs: set, window: tuple[Fraction, Fraction], gate: Any
) -> dict[str, Any]:
    """Boundary-value probes backing the declared continuation.

    Per open channel and window point, the principal-branch jump
    across the window is compared against the declared addition just
    below the axis on a delta ladder.  Acceptance reads the ladder:
    the finest-delta relative residual must sit under the declared
    gate and must not exceed the coarse-delta residual, so finite
    delta truncation is separated from a genuine mismatch.  Channels
    whose declared addition is zero are checked for a vanishing jump
    on the same ladder.  Probe rows are evidence for the declaration,
    typed as consistency checks."""

    x_points = [
        window[0] + (window[1] - window[0]) * Fraction(1, 3),
        window[0] + (window[1] - window[0]) * Fraction(2, 3),
    ]
    dummy = SecondSheetEvaluator([], gate, set())
    rows = []
    all_pass = True
    deltas = sorted(PROBE_DELTAS, reverse=True)
    for m1, m2 in sorted(open_pairs):
        for x in x_points:
            jump_mids = []
            declared_mids = []
            for delta in deltas:
                up = cwc.b0_interval(
                    CInterval.box(x, x, delta, delta),
                    m1,
                    m2,
                    wzp.FIXTURE["mu_ren2"],
                    gate,
                )
                down_point = CInterval.box(x, x, -delta, -delta)
                down = cwc.b0_interval(
                    down_point, m1, m2, wzp.FIXTURE["mu_ren2"], gate
                )
                jump_mids.append(complex((up - down).midpoint()))
                declared_mids.append(
                    complex(dummy.addition((m1, m2), down_point).midpoint())
                )
            if m1 == 0 and m2 == 0:
                magnitudes = [abs(j) for j in jump_mids]
                passed = bool(
                    magnitudes[-1] < float(PROBE_RESIDUAL_GATE)
                    and magnitudes[-1] <= magnitudes[0]
                )
                row_kind = "vanishing_jump"
                residuals = [f"{m:.3e}" for m in magnitudes]
            else:
                residuals_num = [
                    abs(j - d) / max(abs(d), 1e-30)
                    for j, d in zip(jump_mids, declared_mids)
                ]
                passed = bool(
                    residuals_num[-1] < float(PROBE_RESIDUAL_GATE)
                    and residuals_num[-1] <= residuals_num[0]
                )
                row_kind = "declared_addition_match"
                residuals = [f"{r:.3e}" for r in residuals_num]
            all_pass = all_pass and passed
            rows.append(
                {
                    "m1": str(m1),
                    "m2": str(m2),
                    "x": str(x),
                    "kind": row_kind,
                    "deltas": [str(d) for d in deltas],
                    "residuals_coarse_to_fine": residuals,
                    "passed": passed,
                }
            )
    return {"rows": rows, "all_passed": bool(all_pass)}


def interior_certificate(
    compiled: list[dict[str, Any]],
    box: dict[str, tuple[Fraction, Fraction]],
    open_pairs: set,
) -> dict[str, Any]:
    """Exact rational cut exclusion mirrored to the lower half plane.

    The chart facts mirror the zero-exclusion certificate: on a box
    with ``im_hi < 0`` the massless-partner logarithm arguments sit
    strictly inside one open half plane, the two-massive discriminant
    misses the closed negative real axis whenever the threshold-sum
    line misses the box real range, and the Feynman roots stay off the
    real axis.  Open-channel additions reuse the same discriminant
    fact; massless additions are rational with the box off zero."""

    re_lo, re_hi = box["re"]
    im_lo, im_hi = box["im"]
    base_facts = {
        "box_in_open_lower_half_plane": im_hi < 0,
        "box_real_part_positive": re_lo > 0,
    }
    rows = []
    ok = base_facts["box_in_open_lower_half_plane"] and base_facts[
        "box_real_part_positive"
    ]
    hull = CInterval.box(re_lo, re_hi, im_lo, im_hi)
    probe = cwc.IntervalEvaluator([], None)
    seen_dens: set = set()
    for term in compiled:
        den_key = tuple(term["den"])
        if den_key in seen_dens:
            continue
        seen_dens.add(den_key)
        if probe._poly(term["den"], hull).contains_zero():
            ok = False
            rows.append(
                {
                    "kind": "denominator",
                    "den": [[p, str(f)] for p, f in term["den"]],
                    "excludes_zero_on_box": False,
                }
            )
    seen = set()
    for term in compiled:
        if term["head"] != "B0" or term["args"] in seen:
            continue
        seen.add(term["args"])
        m1, m2 = term["args"]
        if m1 == 0 and m2 == 0:
            holo = im_hi < 0
            certificate = "log(s/mu2) only; box in open lower half plane"
        elif m1 == 0 or m2 == 0:
            holo = im_hi < 0 and re_lo > 0
            certificate = (
                "explicit root chart; both logarithm arguments confined "
                "to one open half plane on the box"
            )
        else:
            threshold_sum = m1 + m2
            outside = threshold_sum < re_lo or threshold_sum > re_hi
            holo = bool(outside and im_hi < 0 and re_lo > 0)
            certificate = (
                "discriminant imaginary part nonvanishing off the "
                "threshold-sum line; Feynman roots never real"
            )
        row = {
            "m1": str(m1),
            "m2": str(m2),
            "holomorphic": bool(holo),
            "certificate": certificate,
        }
        if (m1, m2) in open_pairs:
            row["addition_holomorphic"] = bool(holo)
        rows.append(row)
        ok = ok and holo
    return {
        "method": cwc.HOLOMORPHY_METHOD + "_lower_half_plane",
        "base_facts": base_facts,
        "loop_charts": rows,
        "holomorphic_on_box": bool(ok),
    }


def _intersect(a: CInterval, b: CInterval) -> CInterval | None:
    re_lo = max(a.re.a, b.re.a)
    re_hi = min(a.re.b, b.re.b)
    im_lo = max(a.im.a, b.im.a)
    im_hi = min(a.im.b, b.im.b)
    if re_lo > re_hi or im_lo > im_hi:
        return None
    return CInterval(iv.mpf([re_lo, re_hi]), iv.mpf([im_lo, im_hi]))


def _contains(outer: CInterval, inner: CInterval) -> bool:
    return bool(
        inner.re.a >= outer.re.a
        and inner.re.b <= outer.re.b
        and inner.im.a >= outer.im.a
        and inner.im.b <= outer.im.b
    )


def newton_refine(
    evaluator: SecondSheetEvaluator,
    seed: dict[str, tuple[Fraction, Fraction]],
    tree_mass: Fraction,
) -> dict[str, Any]:
    """Interval-Newton contraction to a certified unique root ball.

    ``N(X) = mid(X) - G(mid)/G'(X)``: strict containment of ``N(X)``
    in ``X`` certifies existence and uniqueness of a root in ``X`` and
    the iteration contracts to a tight pole ball.  The derivative
    enclosure over the final ball is the Laurent denominator ball."""

    region = CInterval.box(
        seed["re"][0], seed["re"][1], seed["im"][0], seed["im"][1]
    )
    contracted = False
    width_floor = 1e-12
    try:
        for step_index in range(NEWTON_STEPS):
            widths = (
                float(region.re.b) - float(region.re.a),
                float(region.im.b) - float(region.im.a),
            )
            if contracted and max(widths) < width_floor:
                break
            re_mid = (
                Fraction(float(region.re.a)) + Fraction(float(region.re.b))
            ) / 2
            im_mid = (
                Fraction(float(region.im.a)) + Fraction(float(region.im.b))
            ) / 2
            mid_rect = CInterval.box(re_mid, re_mid, im_mid, im_mid)
            evaluator.loop_cache.clear()
            g_mid = evaluator.inverse_propagator(
                mid_rect, ("nmid", step_index), tree_mass
            )
            evaluator.loop_cache.clear()
            gp_region = evaluator.inverse_propagator_derivative(
                region, ("nreg", step_index), tree_mass
            )
            step = g_mid / gp_region
            candidate = mid_rect - step
            if _contains(region, candidate):
                contracted = True
            intersected = _intersect(region, candidate)
            if intersected is None:
                return {"contracted": False, "reason": "empty intersection"}
            region = intersected
        evaluator.loop_cache.clear()
        gp_final = evaluator.inverse_propagator_derivative(
            region, ("nfin",), tree_mass
        )
        evaluator.loop_cache.clear()
        g_final = evaluator.inverse_propagator(region, ("gfin",), tree_mass)
    except SheetError as error:
        return {"contracted": False, "reason": f"sheet gate: {error}"}
    denominator_excludes_zero = not gp_final.contains_zero()
    residual_contains_zero = g_final.contains_zero()
    residue = None
    if denominator_excludes_zero:
        one = CInterval.from_fraction(1)
        residue = one / gp_final
    return {
        "contracted": bool(contracted),
        "pole_ball": {
            "re": [str(region.re.a), str(region.re.b)],
            "im": [str(region.im.a), str(region.im.b)],
        },
        "derivative_ball": {
            "re": [str(gp_final.re.a), str(gp_final.re.b)],
            "im": [str(gp_final.im.a), str(gp_final.im.b)],
        },
        "denominator_ball_excludes_zero": bool(denominator_excludes_zero),
        "null_vectors": {
            "projection": "transverse_scalar_block_1x1",
            "left_residual_contains_zero": bool(residual_contains_zero),
            "right_residual_contains_zero": bool(residual_contains_zero),
            "laurent_denominator_excludes_zero": bool(
                denominator_excludes_zero
            ),
        },
        "residue_ball": None
        if residue is None
        else {
            "re": [str(residue.re.a), str(residue.re.b)],
            "im": [str(residue.im.a), str(residue.im.b)],
        },
        "_region": region,
        "_derivative": gp_final,
        "_residue": residue,
    }


def certify_pole(
    name: str,
    compiled: list[dict[str, Any]],
    tree_mass: Fraction,
    open_pairs: set,
    precision: int,
    fixed_segments: list | None = None,
) -> dict[str, Any]:
    iv.prec = precision
    gate = cwc._iv_fraction(cwc.ARG_WIDTH_GATE_NUM).a
    geometry = GEOMETRY[name]
    cwc.progress(f"{name} @ {precision} bits: second-sheet interior certificates")
    interior = interior_certificate(compiled, geometry["pole_box"], open_pairs)
    corridor = interior_certificate(compiled, geometry["corridor_box"], open_pairs)
    evaluator = SecondSheetEvaluator(compiled, gate, open_pairs)
    cwc.progress(f"{name} @ {precision} bits: pole-box winding")
    winding = cwc.certify_winding(
        evaluator,
        geometry["pole_box"],
        tree_mass,
        gate,
        fixed_segments=fixed_segments,
    )
    cwc.progress(
        f"{name} @ {precision} bits: winding certified="
        f"{winding['certified']} winding={winding.get('winding')} "
        f"segments={winding.get('segments')}"
    )
    newton = newton_refine(
        SecondSheetEvaluator(compiled, gate, open_pairs),
        geometry["newton_seed"],
        tree_mass,
    )
    probes = cwc.probe_quantities(evaluator, geometry["pole_box"], tree_mass)
    null_vectors = newton.get("null_vectors", {})
    certified = bool(
        interior["holomorphic_on_box"]
        and corridor["holomorphic_on_box"]
        and winding["certified"]
        and winding.get("winding") == 1
        and newton.get("contracted")
        and null_vectors.get("left_residual_contains_zero")
        and null_vectors.get("right_residual_contains_zero")
        and null_vectors.get("laurent_denominator_excludes_zero")
    )
    return {
        "particle": name,
        "precision_bits": precision,
        "pole_box": {
            "re": [str(f) for f in geometry["pole_box"]["re"]],
            "im": [str(f) for f in geometry["pole_box"]["im"]],
            "half_plane": "lower",
        },
        "tree_mass_sq": str(tree_mass),
        "interior_holomorphy": interior,
        "corridor_holomorphy": corridor,
        "boundary_winding": winding,
        "interval_newton": newton,
        "_probes": probes,
        "simple_root_certified": bool(
            winding["certified"] and winding.get("winding") == 1
        ),
        "pole_certified": certified,
        "reading": (
            "the continued inverse propagator has exactly one simple zero "
            "in the declared lower-half pole box on the declared "
            "continuation; the interval-Newton ball encloses it and the "
            "derivative ball excludes zero"
            if certified
            else "not certified"
        ),
    }


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
    cwc.progress("compiling transverse blocks with prefactor corrections")
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

    iv.prec = max(PRECISIONS)
    gate = cwc._iv_fraction(cwc.ARG_WIDTH_GATE_NUM).a
    continuation: dict[str, Any] = {}
    channel_data: dict[str, Any] = {}
    for name in ("W", "Z"):
        data = classify_channels(compiled[name], masses[name])
        channel_data[name] = data
        cwc.progress(f"{name}: continuation consistency probes")
        probes = consistency_probes(data["open_pairs"], data["window"], gate)
        continuation[name] = {
            "window": [str(data["window"][0]), str(data["window"][1])],
            "channels": data["channels"],
            "consistency_probes": probes,
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
    all_certified = all(
        continuation[name]["consistency_probes"]["all_passed"]
        for name in ("W", "Z")
    )
    for name in ("W", "Z"):
        rows = {}
        raw_totals = []
        probe_ladder: list[dict[str, CInterval]] = []
        enclosure_ladder: list[list[CInterval]] = []
        ball_ladder: list[dict[str, CInterval | None]] = []
        base_segments: list | None = None
        for precision in PRECISIONS:
            row = certify_pole(
                name,
                compiled[name],
                masses[name],
                channel_data[name]["open_pairs"],
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
            newton_row = row["interval_newton"]
            ball_ladder.append(
                {
                    "pole_ball": newton_row.pop("_region", None),
                    "derivative_ball": newton_row.pop("_derivative", None),
                    "residue_ball": newton_row.pop("_residue", None),
                }
            )
            checkpoint_rows[f"{name}:{precision}"] = {
                "winding": winding_row.get("winding"),
                "certified": row["pole_certified"],
            }
            CHECKPOINT_PATH.write_text(
                json.dumps(checkpoint_rows, sort_keys=True, indent=1) + "\n",
                encoding="utf-8",
            )
            raw = winding_row.pop("_raw_total", None)
            probe_ladder.append(row.pop("_probes"))
            rows[str(precision)] = row
            all_certified = all_certified and row["pole_certified"]
            if raw is not None:
                raw_totals.append(raw)
        total_nested = (
            all(
                interval_nested(raw_totals[k], raw_totals[k + 1])
                for k in range(len(raw_totals) - 1)
            )
            if len(raw_totals) == len(PRECISIONS)
            else False
        )
        per_quantity = {
            label: all(
                label in probe_ladder[k + 1]
                and cinterval_nested(
                    probe_ladder[k][label], probe_ladder[k + 1][label]
                )
                for k in range(len(probe_ladder) - 1)
            )
            for label in sorted(probe_ladder[0].keys())
        }
        segment_nesting = {
            "segments": len(enclosure_ladder[0]) if enclosure_ladder else 0,
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
        ball_nesting = {
            label: all(
                ball_ladder[k][label] is not None
                and ball_ladder[k + 1][label] is not None
                and cinterval_nested(
                    ball_ladder[k][label], ball_ladder[k + 1][label]
                )
                for k in range(len(ball_ladder) - 1)
            )
            for label in ("pole_ball", "derivative_ball", "residue_ball")
        }
        nesting[name] = {
            "enclosures_nested_with_precision": bool(total_nested),
            "per_quantity_probe_nesting": per_quantity,
            "per_segment_enclosure_nesting": segment_nesting,
            "newton_ball_nesting": ball_nesting,
        }
        all_certified = (
            all_certified
            and total_nested
            and all(per_quantity.values())
            and segment_nesting["all_nested"]
            and all(ball_nesting.values())
        )
        results[name] = rows

    zero_exclusion_sha = (
        sha256_bytes(ZERO_EXCLUSION_PATH.read_bytes())
        if ZERO_EXCLUSION_PATH.exists()
        else None
    )
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
            "zero_exclusion_module_sha256": sha256_bytes(
                (ROOT / "producers" / "certified_wz_contours.py").read_bytes()
            ),
            "verifier_module_sha256": sha256_bytes(Path(__file__).read_bytes()),
            "zero_exclusion_receipt_sha256": zero_exclusion_sha,
        },
        "fixture": {k: str(f) for k, f in wzp.FIXTURE.items()},
        "dimensional_prefactor_finite_correction": {
            name: str(corrections[name]) for name in corrections
        },
        "declared_continuation": continuation,
        "serialized_gates": {
            "precisions_bits": list(PRECISIONS),
            "initial_segments_per_edge": cwc.INITIAL_SEGMENTS_PER_EDGE,
            "max_subdivision_depth": cwc.MAX_SUBDIVISION_DEPTH,
            "holomorphy_method": cwc.HOLOMORPHY_METHOD + "_lower_half_plane",
            "boundary_method": cwc.BOUNDARY_METHOD,
            "arg_width_gate": str(cwc.ARG_WIDTH_GATE_NUM),
            "winding_tolerance": str(cwc.WINDING_TOLERANCE_NUM),
            "newton_steps": NEWTON_STEPS,
            "probe_residual_gate": str(PROBE_RESIDUAL_GATE),
        },
        "sheet_statement": (
            "the pole boxes lie strictly in the open lower half plane on "
            "the declared continuation through the recorded cut window; "
            "the identification with the continuation is a declared "
            "convention backed by the recorded boundary-value probes; "
            "the principal-sheet zero exclusion is the separate receipt "
            "pinned above"
        ),
        "results": results,
        "precision_nesting": nesting,
        "promotion": {
            "second_sheet_pole_certified_on_declared_continuation": all_certified,
            "simple_root_certified": all_certified,
            "laurent_denominator_ball_certified": all_certified,
            "scalar_residue_ball_certified": all_certified,
            "matrix_rank_laurent_certified": False,
            "bmhv_restoration_certified": False,
            "physical_current_claim": False,
            "oph_native": False,
            "unit_claim": False,
        },
        "status": STATUS_CERTIFIED if all_certified else STATUS_FAILED,
    }
    payload["receipt_sha256"] = sha256_bytes(
        cwc.canonical_json(
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
            print("SECOND_SHEET_POLE_DRIFT", file=sys.stderr)
            return 1
        print("SECOND_SHEET_POLE_VERIFIED")
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
