#!/usr/bin/env python3
"""Scalar W/Z pole enclosures on a declared algebraic continuation chart.

This producer certifies a deliberately narrow statement: the scalar
transverse inverse-propagator function built from the frozen synthetic
fixture has one simple zero on an explicitly declared algebraic chart.
It does not certify that chart as *the* physical or unique resonance
sheet, does not certify the full neutral/charged matrix rank-``n-1``
Laurent hypothesis, and therefore does not discharge issue #593's
combined pole-plus-Laurent acceptance row.

Declared chart.  Each channel carries an explicit base-chart label and
an explicit crossed/principal sheet action.  The crossing window is
declared input, not inferred from the tree mass.  In particular, the W
window is below the W-gamma threshold, so W-gamma receives the declared
``principal`` action—no added chart-relative correction.  This label is
not a certified physical-sheet assignment.

The algebraic corrections are relative to the implemented base charts:

* two-massive crossed channels add ``2 pi i sqrt(lambda(s))/s``;
* either ordering of a one-mass channel uses the same mass-exchange-
  symmetric lower-principal base chart, and its crossed-chart correction
  is ``2 pi i (1 - m/s)``;
* the both-massless base formula is already the upper-half continuation,
  so its relative correction is zero.

Finite-delta boundary-value probes are serialized as non-certifying
consistency diagnostics only.  They do not prove the continuation
identity and do not gate the scalar pole certificate.  The algebraic
chart identity is declared; independent certification of that identity
remains false.

Certified on each declared scalar pole box, per precision:

* interior holomorphy of the continued formula on the pole box and on
  the declared corridor box by exact rational cut exclusion;
* boundary winding one by the centered-form segment enclosures and
  endpoint-ratio argument chaining imported from the zero-exclusion
  verifier: winding one counts one root with multiplicity one, so the
  root is simple;
* an interval-Newton image strictly inside its source box, with exact
  dyadic endpoint evidence and a strictly positive recorded margin,
  certifying a unique root inside the declared refinement box;
* a derivative ball over the refinement box excluding zero: the
  Laurent denominator of the propagator residue;
* the scalar Laurent residue ball ``1/G'`` over the refinement box;
* replayable raw partition, interval, Newton, gate, and precision-
  nesting evidence.

What is not certified and stays false in the receipt: independently
certified continuation identity, standard/physical second-sheet
identification, a sign bridge from the pinned engine convention
``G = s - m_tree^2 - Pi_engine`` to the separately written theorem
convention, full matrix rank-``n-1`` Laurent data, issue #593 row closure,
BMHV restoration, physical-current normalization, unitarity, any unit
claim, and any OPH-native statement.
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

SCHEMA = "oph.certified_second_sheet_poles.v2"
STATUS_CERTIFIED = "SCALAR_POLE_CERTIFIED_ON_DECLARED_ALGEBRAIC_CHART"
STATUS_FAILED = "CERTIFICATION_INCOMPLETE"

PRECISIONS = cwc.PRECISIONS
NEWTON_STEPS = 4
PROBE_DELTAS = (Fraction(1, 10**6), Fraction(1, 10**8))
PROBE_RESIDUAL_GATE = Fraction(1, 1000)

# These windows are declared chart input.  A channel is crossed exactly
# when its threshold is at or below the lower window endpoint.  The
# open interval then contains no threshold and lies strictly left of
# the pole box.  The W-gamma threshold is the upper W endpoint and gets
# the declared "principal" action (no chart-relative correction), not a
# certified physical-sheet identification.
DECLARED_WINDOWS = {
    "W": (Fraction(9, 100), Fraction(1, 9)),
    "Z": (Fraction(4, 25), Fraction(25, 144)),
}

# Declared boxes: the pole box carries the winding-one certificate, the
# seed box starts the interval-Newton refinement, and the corridor box
# connects the cut window to the pole box inside the open lower half
# plane.  All endpoints are exact rationals.
GEOMETRY = {
    "W": {
        "pole_box": {
            "re": (Fraction(1119, 10000), Fraction(14, 125)),
            "im": (Fraction(-9, 10000), Fraction(-7, 10000)),
        },
        "newton_seed": {
            "re": (Fraction(111953, 10**6), Fraction(111955, 10**6)),
            "im": (Fraction(-803, 10**6), Fraction(-802, 10**6)),
        },
        "corridor_box": {
            "re": (Fraction(9, 100), Fraction(113, 1000)),
            "im": (Fraction(-1, 1000), Fraction(-1, 10**9)),
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


def _mpf_tuple_fraction(raw: tuple[int, int, int, int]) -> Fraction:
    """Convert an mpmath binary endpoint tuple to its exact dyadic value."""

    sign, mantissa, exponent, _bit_count = raw
    value = Fraction(mantissa)
    if exponent >= 0:
        value *= 2**exponent
    else:
        value /= 2 ** (-exponent)
    return -value if sign else value


def _endpoint_fraction(interval: Any, upper: bool = False) -> Fraction:
    """Return one directed interval endpoint without a float round trip."""

    return _mpf_tuple_fraction(interval._mpi_[1 if upper else 0])


def _complex_interval_evidence(
    interval: CInterval,
) -> dict[str, dict[str, str]]:
    return cwc.serialize_cinterval(interval)


class SecondSheetEvaluator(cwc.IntervalEvaluator):
    """Continued transverse evaluation on the declared algebraic chart.

    Crossed channels receive the chart-relative correction on top of
    the implemented base value; principal channels are untouched.  The
    one-mass base chart is explicitly mass-exchange symmetric, so both
    mass orderings receive the same correction.  All additions and
    derivatives run through the same sheet-gated interval layer."""

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
        if m1 == 0 and m2 == 0:
            return CInterval.from_fraction(0)
        if m1 == 0 or m2 == 0:
            mass = m2 if m1 == 0 else m1
            one = CInterval.from_fraction(1)
            return self._two_pi_i() * (
                one - CInterval.from_fraction(mass) / s
            )
        sq = cwc._feynman_roots(s, m1, m2, self.gate)[2]
        return self._two_pi_i() * sq / s

    def addition_slope(self, args: tuple, s: CInterval) -> CInterval:
        m1, m2 = args
        if m1 == 0 and m2 == 0:
            return CInterval.from_fraction(0)
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


def _base_chart_and_correction(
    m1: Fraction, m2: Fraction
) -> tuple[str, str]:
    """Return the implemented base chart and its relative crossing correction."""

    if m1 == 0 and m2 == 0:
        return "both_massless_upper_half_continued", "0"
    if m1 == 0 or m2 == 0:
        mass = m2 if m1 == 0 else m1
        return (
            "one_mass_symmetric_lower_principal",
            f"2*pi*i*(1-({mass})/s)",
        )
    return "two_mass_principal_root_chart", "2*pi*i*sqrt(lambda(s))/s"


def classify_channels(
    compiled: list[dict[str, Any]], name: str, tree_mass: Fraction
) -> dict[str, Any]:
    """Replay the exact sheet vector from an explicit declared window.

    The lower endpoint is the greatest crossed threshold.  Thresholds
    at or below it are crossed; all others remain principal.  The open
    window itself must contain no threshold and must end strictly left
    of the pole box.  ``tree_mass`` is recorded as fixture context only
    and never decides a channel sheet.
    """

    window = DECLARED_WINDOWS[name]
    pole_re_lo = GEOMETRY[name]["pole_box"]["re"][0]
    if not (window[0] < window[1] < pole_re_lo):
        raise SystemExit(f"{name} declared crossing window has invalid geometry")

    pairs = sorted({term["args"] for term in compiled if term["head"] == "B0"})
    rows = []
    crossed_pairs = set()
    crossed_thresholds = []
    principal_thresholds = []
    for m1, m2 in pairs:
        threshold = (rational_sqrt(m1) + rational_sqrt(m2)) ** 2
        if window[0] < threshold < window[1]:
            raise SystemExit(f"{name} threshold {threshold} lies inside the window")
        crossed = threshold <= window[0]
        if crossed:
            crossed_pairs.add((m1, m2))
            crossed_thresholds.append(threshold)
        else:
            principal_thresholds.append(threshold)
        base_chart, crossing_correction = _base_chart_and_correction(m1, m2)
        rows.append(
            {
                "m1": str(m1),
                "m2": str(m2),
                "threshold": str(threshold),
                "sheet_action": "crossed" if crossed else "principal",
                "base_chart": base_chart,
                "crossing_correction": crossing_correction,
                "applied_correction": crossing_correction if crossed else "0",
                "correction_reference": "relative_to_declared_base_chart",
            }
        )
    if not crossed_thresholds or max(crossed_thresholds) != window[0]:
        raise SystemExit(f"{name} window lower endpoint is not the crossed frontier")
    if principal_thresholds and min(principal_thresholds) < window[1]:
        raise SystemExit(f"{name} principal threshold lies left of window end")
    return {
        "channels": rows,
        "crossed_pairs": crossed_pairs,
        "window": window,
        "tree_mass_sq_context": tree_mass,
        "sheet_vector": [
            {
                "m1": row["m1"],
                "m2": row["m2"],
                "sheet_action": row["sheet_action"],
            }
            for row in rows
        ],
    }


def consistency_probes(
    crossed_pairs: set, window: tuple[Fraction, Fraction], gate: Any
) -> dict[str, Any]:
    """Non-certifying finite-delta diagnostics for the declared chart.

    These point probes compare upper/lower finite-offset values with
    the declared chart-relative correction.  They can expose a sign or
    branch mismatch, but finite deltas and midpoint arithmetic do
    not certify the continuation identity.  Their aggregate never
    gates the scalar root theorem.
    """

    x_points = [
        window[0] + (window[1] - window[0]) * Fraction(1, 3),
        window[0] + (window[1] - window[0]) * Fraction(2, 3),
    ]
    dummy = SecondSheetEvaluator([], gate, set())
    rows = []
    all_pass = True
    deltas = sorted(PROBE_DELTAS, reverse=True)
    for m1, m2 in sorted(crossed_pairs):
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
                    "within_diagnostic_gate": passed,
                }
            )
    return {
        "role": "non_certifying_finite_delta_diagnostic",
        "gates_scalar_certificate": False,
        "residual_gate": str(PROBE_RESIDUAL_GATE),
        "rows": rows,
        "all_within_diagnostic_gate": bool(all_pass),
    }


def interior_certificate(
    compiled: list[dict[str, Any]],
    box: dict[str, tuple[Fraction, Fraction]],
    crossed_pairs: set,
) -> dict[str, Any]:
    """Exact rational cut exclusion mirrored to the lower half plane.

    The chart facts mirror the zero-exclusion certificate: on a box
    with ``im_hi < 0`` the massless-partner logarithm arguments sit
    strictly inside one open half plane, the two-massive discriminant
    misses the closed negative real axis whenever the threshold-sum
    line misses the box real range, and the Feynman roots stay off the
    real axis.  Crossed-channel corrections reuse the same
    discriminant fact; one-mass corrections are rational with the box
    off zero."""

    re_lo, re_hi = box["re"]
    im_lo, im_hi = box["im"]
    base_facts = {
        "box_in_open_lower_half_plane": im_hi < 0,
        "box_real_part_positive": re_lo > 0,
    }
    rows = []
    hull = CInterval.box(re_lo, re_hi, im_lo, im_hi)
    denominator_specs: dict[str, list] = {}
    for term in compiled:
        coefficients = [
            [int(power), str(fraction)] for power, fraction in term["den"]
        ]
        denominator_id = sha256_bytes(
            cwc.canonical_json(coefficients).encode("utf-8")
        )
        denominator_specs.setdefault(denominator_id, term["den"])
    denominator_probe = cwc.IntervalEvaluator([], None)
    denominator_witnesses = []
    denominator_enclosures: dict[str, CInterval] = {}
    for denominator_id in sorted(denominator_specs):
        coefficients_raw = denominator_specs[denominator_id]
        coefficients = [
            [int(power), str(fraction)]
            for power, fraction in coefficients_raw
        ]
        enclosure = denominator_probe._poly(coefficients_raw, hull)
        denominator_enclosures[denominator_id] = enclosure
        excludes_zero = not enclosure.contains_zero()
        denominator_witnesses.append(
            {
                "denominator_id": denominator_id,
                "coefficients": coefficients,
                "enclosure": _complex_interval_evidence(enclosure),
                "excludes_zero": bool(excludes_zero),
                "zero_exclusion_abs2_lower": str(
                    cwc._zero_exclusion_abs2_lower(enclosure)
                ),
            }
        )
    denominators_ok = all(
        witness["excludes_zero"] for witness in denominator_witnesses
    )
    pairs_ok = True
    seen = set()
    for term in compiled:
        if term["head"] != "B0" or term["args"] in seen:
            continue
        seen.add(term["args"])
        m1, m2 = term["args"]
        if m1 == 0 and m2 == 0:
            holo = im_hi < 0
            row = {
                "m1": str(m1),
                "m2": str(m2),
                "chart": "both_masses_zero",
                "certificate": "log(s/mu2) only; box in open lower half plane",
                "holomorphic": bool(holo),
            }
        elif m1 == 0 or m2 == 0:
            holo = im_hi < 0 and re_lo > 0
            sign = "positive"
            base_chart, _correction = _base_chart_and_correction(m1, m2)
            row = {
                "m1": str(m1),
                "m2": str(m2),
                "chart": base_chart,
                "root_chart_log_argument_imaginary_sign": sign,
                "certificate": (
                    "declared mass-exchange-symmetric one-mass algebraic "
                    f"base chart; root-chart arguments have strictly {sign} "
                    "imaginary part and chart-relative rational corrections "
                    "are holomorphic because the box excludes s=0"
                ),
                "holomorphic": bool(holo),
            }
        else:
            symmetry_line = m1 + m2
            outside = symmetry_line < re_lo or symmetry_line > re_hi
            holo = bool(outside and im_hi < 0 and re_lo > 0)
            row = {
                "m1": str(m1),
                "m2": str(m2),
                "chart": "two_massive",
                "discriminant_symmetry_line": str(symmetry_line),
                "discriminant_symmetry_line_outside_box_re_range": bool(
                    outside
                ),
                "certificate": (
                    "discriminant imaginary part nonvanishing off the "
                    "discriminant-symmetry line Re(s)=m1+m2 (not the "
                    "physical threshold); Feynman roots never real"
                ),
                "holomorphic": bool(holo),
            }
        if (m1, m2) in crossed_pairs:
            row["chart_correction_holomorphic"] = bool(holo)
        rows.append(row)
        pairs_ok = pairs_ok and holo
    holomorphic = bool(
        base_facts["box_in_open_lower_half_plane"]
        and base_facts["box_real_part_positive"]
        and denominators_ok
        and pairs_ok
    )
    return {
        "method": cwc.HOLOMORPHY_METHOD + "_lower_half_plane",
        "base_facts": base_facts,
        "coefficient_denominators_exclude_zero_on_box": denominators_ok,
        "coefficient_denominator_witnesses": denominator_witnesses,
        "loop_charts": rows,
        "holomorphic_on_box": holomorphic,
        "_coefficient_denominator_enclosures": denominator_enclosures,
    }


def _intersect(a: CInterval, b: CInterval) -> CInterval | None:
    re_lo = max(a.re.a, b.re.a)
    re_hi = min(a.re.b, b.re.b)
    im_lo = max(a.im.a, b.im.a)
    im_hi = min(a.im.b, b.im.b)
    if re_lo > re_hi or im_lo > im_hi:
        return None
    return CInterval(iv.mpf([re_lo, re_hi]), iv.mpf([im_lo, im_hi]))


def _strict_inclusion_evidence(
    outer: CInterval, inner: CInterval
) -> dict[str, Any]:
    """Exact endpoint margins for ``inner`` lying in ``interior(outer)``."""

    margins = {
        "re_lower": _endpoint_fraction(inner.re) - _endpoint_fraction(outer.re),
        "re_upper": _endpoint_fraction(outer.re, upper=True)
        - _endpoint_fraction(inner.re, upper=True),
        "im_lower": _endpoint_fraction(inner.im) - _endpoint_fraction(outer.im),
        "im_upper": _endpoint_fraction(outer.im, upper=True)
        - _endpoint_fraction(inner.im, upper=True),
    }
    minimum = min(margins.values())
    return {
        "margins": {key: str(value) for key, value in margins.items()},
        "minimum_margin": str(minimum),
        "strict": bool(minimum > 0),
    }


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
    strict_inclusion_certified = False
    proof_step_index: int | None = None
    proof_margin: dict[str, Any] | None = None
    width_floor = Fraction(1, 10**12)
    iterations: list[dict[str, Any]] = []
    try:
        for step_index in range(NEWTON_STEPS):
            widths = (
                _endpoint_fraction(region.re, upper=True)
                - _endpoint_fraction(region.re),
                _endpoint_fraction(region.im, upper=True)
                - _endpoint_fraction(region.im),
            )
            if strict_inclusion_certified and max(widths) < width_floor:
                break
            re_mid = (
                _endpoint_fraction(region.re)
                + _endpoint_fraction(region.re, upper=True)
            ) / 2
            im_mid = (
                _endpoint_fraction(region.im)
                + _endpoint_fraction(region.im, upper=True)
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
            inclusion = _strict_inclusion_evidence(region, candidate)
            if inclusion["strict"] and not strict_inclusion_certified:
                strict_inclusion_certified = True
                proof_step_index = step_index
                proof_margin = inclusion
            intersected = _intersect(region, candidate)
            if intersected is None:
                return {
                    "contracted": False,
                    "strict_interior_inclusion_certified": False,
                    "reason": "empty intersection",
                    "iterations": iterations,
                }
            iterations.append(
                {
                    "step_index": step_index,
                    "source_box": _complex_interval_evidence(region),
                    "midpoint": [str(re_mid), str(im_mid)],
                    "midpoint_value": _complex_interval_evidence(g_mid),
                    "derivative_on_source": _complex_interval_evidence(gp_region),
                    "newton_image": _complex_interval_evidence(candidate),
                    "strict_inclusion": inclusion,
                    "intersection_box": _complex_interval_evidence(intersected),
                }
            )
            region = intersected
        evaluator.loop_cache.clear()
        gp_final = evaluator.inverse_propagator_derivative(
            region, ("nfin",), tree_mass
        )
        evaluator.loop_cache.clear()
        g_final = evaluator.inverse_propagator(region, ("gfin",), tree_mass)
    except SheetError as error:
        return {
            "contracted": False,
            "strict_interior_inclusion_certified": False,
            "reason": f"sheet gate: {error}",
            "iterations": iterations,
        }
    denominator_excludes_zero = not gp_final.contains_zero()
    residual_contains_zero = g_final.contains_zero()
    residue = None
    if denominator_excludes_zero:
        one = CInterval.from_fraction(1)
        residue = one / gp_final
    return {
        "contracted": bool(strict_inclusion_certified),
        "strict_interior_inclusion_certified": bool(
            strict_inclusion_certified
        ),
        "proof_step_index": proof_step_index,
        "strict_inclusion_margin": proof_margin,
        "seed_box": {
            "re": [str(seed["re"][0]), str(seed["re"][1])],
            "im": [str(seed["im"][0]), str(seed["im"][1])],
        },
        "iterations": iterations,
        "pole_ball": _complex_interval_evidence(region),
        "residual_ball": _complex_interval_evidence(g_final),
        "derivative_ball": _complex_interval_evidence(gp_final),
        "denominator_ball_excludes_zero": bool(denominator_excludes_zero),
        "null_vectors": {
            "projection": "transverse_scalar_block_1x1",
            "full_matrix_rank_n_minus_1_certified": False,
            "left_residual_contains_zero": bool(residual_contains_zero),
            "right_residual_contains_zero": bool(residual_contains_zero),
            "laurent_denominator_excludes_zero": bool(
                denominator_excludes_zero
            ),
        },
        "residue_ball": None
        if residue is None
        else _complex_interval_evidence(residue),
        "_region": region,
        "_derivative": gp_final,
        "_residual": g_final,
        "_residue": residue,
    }


def declared_chart_segment_quantities(
    evaluator: SecondSheetEvaluator,
    segments: list,
    crossed_pairs: set[tuple[Fraction, Fraction]],
) -> tuple[list[dict[str, Any]], dict[str, CInterval]]:
    """Emit unambiguous base/correction/combined B0 and B0' evidence."""

    rows: list[dict[str, Any]] = []
    quantities: dict[str, CInterval] = {}
    for start, end in segments:
        segment_id = cwc._segment_id(start, end)
        hull = cwc._segment_hull(start, end)
        for m1, m2 in sorted(crossed_pairs):
            args = (m1, m2)
            base = cwc.b0_interval(
                hull, m1, m2, wzp.FIXTURE["mu_ren2"], evaluator.gate
            )
            correction = evaluator.addition(args, hull)
            combined = base + correction
            base_slope = cwc.b0p_interval(
                hull, m1, m2, wzp.FIXTURE["mu_ren2"], evaluator.gate
            )
            correction_slope = evaluator.addition_slope(args, hull)
            combined_slope = base_slope + correction_slope
            rendered = f"{m1},{m2}"
            values = {
                f"segment:{segment_id}:base_chart:B0({rendered})": base,
                f"segment:{segment_id}:chart_correction:B0({rendered})": correction,
                f"segment:{segment_id}:combined:B0({rendered})": combined,
                f"segment:{segment_id}:base_chart:B0p({rendered})": base_slope,
                f"segment:{segment_id}:chart_correction:B0p({rendered})": correction_slope,
                f"segment:{segment_id}:combined:B0p({rendered})": combined_slope,
            }
            quantities.update(values)
            base_chart, crossing_formula = _base_chart_and_correction(m1, m2)
            rows.append(
                {
                    "segment_id": segment_id,
                    "m1": str(m1),
                    "m2": str(m2),
                    "base_chart": base_chart,
                    "crossing_correction_formula": crossing_formula,
                    "base_B0": _complex_interval_evidence(base),
                    "chart_correction_B0": _complex_interval_evidence(correction),
                    "combined_B0": _complex_interval_evidence(combined),
                    "base_B0p": _complex_interval_evidence(base_slope),
                    "chart_correction_B0p": _complex_interval_evidence(
                        correction_slope
                    ),
                    "combined_B0p": _complex_interval_evidence(
                        combined_slope
                    ),
                }
            )
    return rows, quantities


def certify_pole(
    name: str,
    compiled: list[dict[str, Any]],
    tree_mass: Fraction,
    crossed_pairs: set,
    precision: int,
    fixed_segments: list | None = None,
) -> dict[str, Any]:
    iv.prec = precision
    gate = cwc._iv_fraction(cwc.ARG_WIDTH_GATE_NUM).a
    geometry = GEOMETRY[name]
    cwc.progress(f"{name} @ {precision} bits: declared-chart interior certificates")
    interior = interior_certificate(compiled, geometry["pole_box"], crossed_pairs)
    corridor = interior_certificate(
        compiled, geometry["corridor_box"], crossed_pairs
    )
    denominator_enclosures = {
        f"pole_box:{identifier}": enclosure
        for identifier, enclosure in interior.pop(
            "_coefficient_denominator_enclosures"
        ).items()
    } | {
        f"corridor_box:{identifier}": enclosure
        for identifier, enclosure in corridor.pop(
            "_coefficient_denominator_enclosures"
        ).items()
    }
    evaluator = SecondSheetEvaluator(compiled, gate, crossed_pairs)
    cwc.progress(f"{name} @ {precision} bits: pole-box winding")
    winding = cwc.certify_winding(
        evaluator,
        geometry["pole_box"],
        tree_mass,
        gate,
        fixed_segments=fixed_segments,
        collect_evidence=True,
    )
    cwc.progress(
        f"{name} @ {precision} bits: winding certified="
        f"{winding['certified']} winding={winding.get('winding')} "
        f"segments={winding.get('segments')}"
    )
    raw_segments = winding.get("_segments") or []
    if not winding.get("certified") or not raw_segments:
        raise SystemExit(f"{name} @ {precision}: no certified partition evidence")
    chart_rows, chart_quantities = declared_chart_segment_quantities(
        evaluator, raw_segments, crossed_pairs
    )
    newton = newton_refine(
        SecondSheetEvaluator(compiled, gate, crossed_pairs),
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
        and Fraction(winding["winding_tolerance_slack"]) > 0
        and newton.get("strict_interior_inclusion_certified")
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
        "declared_chart_segment_quantities": chart_rows,
        "interval_newton": newton,
        "_probes": probes,
        "_declared_chart_quantity_enclosures": chart_quantities,
        "_coefficient_denominator_enclosures": denominator_enclosures,
        "simple_scalar_root_certified": bool(
            interior["holomorphic_on_box"]
            and winding["certified"]
            and winding.get("winding") == 1
        ),
        "scalar_pole_certified": certified,
        "reading": (
            "the scalar transverse inverse propagator has exactly one "
            "simple zero in the declared lower-half pole box on the "
            "declared algebraic chart; this does not certify a unique or "
            "physical resonance sheet or full matrix Laurent data"
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
        data = classify_channels(compiled[name], name, masses[name])
        channel_data[name] = data
        cwc.progress(f"{name}: continuation consistency probes")
        diagnostics = consistency_probes(
            data["crossed_pairs"], data["window"], gate
        )
        continuation[name] = {
            "window": [str(data["window"][0]), str(data["window"][1])],
            "tree_mass_sq_context": str(data["tree_mass_sq_context"]),
            "chart_definition": {
                "identity_status": "declared_not_independently_certified",
                "correction_reference": "relative_to_declared_base_chart",
                "one_mass_mass_exchange_symmetric": True,
                "window_is_declared_input": True,
            },
            "channels": data["channels"],
            "sheet_vector": data["sheet_vector"],
            "consistency_diagnostics": diagnostics,
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
    for name in ("W", "Z"):
        rows = {}
        raw_totals = []
        probe_ladder: list[dict[str, CInterval]] = []
        denominator_ladder: list[dict[str, CInterval]] = []
        enclosure_ladder: list[list[CInterval]] = []
        segment_record_ladder: list[list[dict[str, Any]]] = []
        partition_id_ladder: list[list[str]] = []
        ball_ladder: list[dict[str, CInterval | None]] = []
        base_segments: list | None = None
        for precision in PRECISIONS:
            row = certify_pole(
                name,
                compiled[name],
                masses[name],
                channel_data[name]["crossed_pairs"],
                precision,
                fixed_segments=base_segments,
            )
            winding_row = row["boundary_winding"]
            segments = winding_row.pop("_segments", None)
            enclosures = winding_row.pop("_segment_enclosures", None)
            segment_records = winding_row.pop("_segment_records", None)
            winding_quantities = winding_row.pop("_quantity_enclosures", {})
            raw = winding_row.pop("_raw_total", None)
            if base_segments is None and segments is not None:
                base_segments = segments
            if enclosures is not None:
                enclosure_ladder.append(enclosures)
            if segment_records is not None:
                segment_record_ladder.append(segment_records)
            partition_id_ladder.append(
                [
                    evidence["segment_id"]
                    for evidence in winding_row.get("segment_evidence", [])
                ]
            )
            newton_row = row["interval_newton"]
            ball_ladder.append(
                {
                    "pole_ball": newton_row.pop("_region", None),
                    "derivative_ball": newton_row.pop("_derivative", None),
                    "residual_ball": newton_row.pop("_residual", None),
                    "residue_ball": newton_row.pop("_residue", None),
                }
            )
            checkpoint_rows[f"{name}:{precision}"] = {
                "winding": winding_row.get("winding"),
                "certified": row["scalar_pole_certified"],
            }
            CHECKPOINT_PATH.write_text(
                json.dumps(checkpoint_rows, sort_keys=True, indent=1) + "\n",
                encoding="utf-8",
            )
            raw_probes = row.pop("_probes")
            probes: dict[str, CInterval] = {}
            for label, value in raw_probes.items():
                if label.startswith("probe:center:integral:B0II_add("):
                    normalized = label.replace(
                        "probe:center:integral:B0II_add(",
                        "probe:center:chart_correction:B0(",
                        1,
                    )
                elif label.startswith("probe:center:integral:B0("):
                    normalized = label.replace(
                        "probe:center:integral:B0(",
                        "probe:center:base_chart:B0(",
                        1,
                    )
                elif label.startswith("probe:center:derivative:B0p("):
                    normalized = label.replace(
                        "probe:center:derivative:B0p(",
                        "probe:center:base_chart:B0p(",
                        1,
                    )
                elif label.startswith(
                    (
                        "probe:center:inverse_propagator",
                        "probe:center:integral:A0(",
                        "probe:center:integral:A0p(",
                    )
                ):
                    normalized = label
                else:
                    raise SystemExit(f"unclassified center-probe quantity {label}")
                probes[normalized] = value
            center_re = sum(GEOMETRY[name]["pole_box"]["re"], Fraction()) / 2
            center_im = sum(GEOMETRY[name]["pole_box"]["im"], Fraction()) / 2
            center_rect = CInterval.box(
                center_re, center_re, center_im, center_im
            )
            center_evaluator = SecondSheetEvaluator(
                compiled[name],
                cwc._iv_fraction(cwc.ARG_WIDTH_GATE_NUM).a,
                channel_data[name]["crossed_pairs"],
            )
            for m1, m2 in sorted(channel_data[name]["crossed_pairs"]):
                rendered = f"{m1},{m2}"
                base = cwc.b0_interval(
                    center_rect,
                    m1,
                    m2,
                    wzp.FIXTURE["mu_ren2"],
                    center_evaluator.gate,
                )
                correction = center_evaluator.addition((m1, m2), center_rect)
                base_slope = cwc.b0p_interval(
                    center_rect,
                    m1,
                    m2,
                    wzp.FIXTURE["mu_ren2"],
                    center_evaluator.gate,
                )
                correction_slope = center_evaluator.addition_slope(
                    (m1, m2), center_rect
                )
                probes[
                    f"probe:center:base_chart:B0({rendered})"
                ] = base
                probes[
                    f"probe:center:chart_correction:B0({rendered})"
                ] = correction
                probes[
                    f"probe:center:combined:B0({rendered})"
                ] = base + correction
                probes[
                    f"probe:center:base_chart:B0p({rendered})"
                ] = base_slope
                probes[
                    f"probe:center:chart_correction:B0p({rendered})"
                ] = correction_slope
                probes[
                    f"probe:center:combined:B0p({rendered})"
                ] = base_slope + correction_slope
            normalized_winding_quantities: dict[str, CInterval] = {}
            for label, value in winding_quantities.items():
                if ":integral:B0II_add(" in label:
                    normalized = label.replace(
                        ":integral:B0II_add(",
                        ":chart_correction:B0(",
                        1,
                    )
                elif ":integral:B0(" in label:
                    normalized = label.replace(
                        ":integral:B0(", ":base_chart:B0(", 1
                    )
                elif ":derivative:B0p(" in label:
                    normalized = label.replace(
                        ":derivative:B0p(", ":base_chart:B0p(", 1
                    )
                elif label.startswith("probe:center:integral:A0("):
                    normalized = label
                else:
                    raise SystemExit(
                        f"unclassified winding quantity {label}"
                    )
                normalized_winding_quantities[normalized] = value
            winding_quantities = normalized_winding_quantities
            chart_quantities = row.pop(
                "_declared_chart_quantity_enclosures"
            )
            quantities = winding_quantities | probes | chart_quantities
            denominator_ladder.append(
                row.pop("_coefficient_denominator_enclosures")
            )
            probe_ladder.append(quantities)
            row["quantity_enclosures"] = {
                label: _complex_interval_evidence(value)
                for label, value in sorted(quantities.items())
            }
            rows[str(precision)] = row
            all_certified = all_certified and row["scalar_pole_certified"]
            if raw is not None:
                raw_totals.append(raw)
        partition_ids_match = bool(
            len(partition_id_ladder) == len(PRECISIONS)
            and bool(partition_id_ladder[0])
            and all(
                partition_id_ladder[index] == partition_id_ladder[0]
                for index in range(1, len(partition_id_ladder))
            )
        )
        segment_counts_match = bool(
            len(segment_record_ladder) == len(PRECISIONS)
            and bool(segment_record_ladder[0])
            and all(
                len(segment_record_ladder[index])
                == len(segment_record_ladder[0])
                for index in range(len(segment_record_ladder))
            )
        )

        def segment_cinterval_field_nested(field: str) -> bool:
            return bool(
                segment_counts_match
                and all(
                    cinterval_nested(
                        segment_record_ladder[index][segment][field],
                        segment_record_ladder[index + 1][segment][field],
                    )
                    for index in range(len(segment_record_ladder) - 1)
                    for segment in range(len(segment_record_ladder[0]))
                )
            )

        def segment_real_interval_field_nested(field: str) -> bool:
            return bool(
                segment_counts_match
                and all(
                    interval_nested(
                        segment_record_ladder[index][segment][field],
                        segment_record_ladder[index + 1][segment][field],
                    )
                    for index in range(len(segment_record_ladder) - 1)
                    for segment in range(len(segment_record_ladder[0]))
                )
            )

        center_values_nested = segment_cinterval_field_nested("_center_value")
        derivative_hulls_nested = segment_cinterval_field_nested(
            "_derivative_hull"
        )
        offsets_nested = segment_cinterval_field_nested("_offset")
        images_nested = segment_cinterval_field_nested("_image")
        rotated_images_nested = segment_cinterval_field_nested(
            "_rotated_image"
        )
        rotated_arguments_nested = segment_real_interval_field_nested(
            "_rotated_argument"
        )
        start_values_nested = segment_cinterval_field_nested("_start_value")
        end_values_nested = segment_cinterval_field_nested("_end_value")
        endpoint_values_nested = bool(
            start_values_nested and end_values_nested
        )
        endpoint_ratios_nested = segment_cinterval_field_nested(
            "_endpoint_ratio"
        )
        endpoint_increments_nested = segment_real_interval_field_nested(
            "_endpoint_increment"
        )
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
            label: bool(
                all(
                    set(probes) == set(probe_ladder[0])
                    for probes in probe_ladder
                )
                and all(
                    cinterval_nested(
                        probe_ladder[k][label],
                        probe_ladder[k + 1][label],
                    )
                    for k in range(len(probe_ladder) - 1)
                )
            )
            for label in probe_names
        }
        segment_nesting = {
            "segments": len(enclosure_ladder[0]) if enclosure_ladder else 0,
            "ladders_compared": max(len(enclosure_ladder) - 1, 0),
            "partition_ids_match": partition_ids_match,
            "center_values_nested": center_values_nested,
            "derivative_hulls_nested": derivative_hulls_nested,
            "offsets_nested": offsets_nested,
            "images_nested": images_nested,
            "rotated_images_nested": rotated_images_nested,
            "rotated_arguments_nested": rotated_arguments_nested,
            "endpoint_values_nested": endpoint_values_nested,
            "endpoint_ratios_nested": endpoint_ratios_nested,
            "endpoint_increments_nested": endpoint_increments_nested,
            "all_nested": bool(
                partition_ids_match
                and len(enclosure_ladder) == len(PRECISIONS)
                and segment_counts_match
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
                and center_values_nested
                and derivative_hulls_nested
                and offsets_nested
                and images_nested
                and rotated_images_nested
                and rotated_arguments_nested
                and endpoint_values_nested
                and endpoint_ratios_nested
                and endpoint_increments_nested
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
            for label in (
                "pole_ball",
                "derivative_ball",
                "residual_ball",
                "residue_ball",
            )
        }
        denominator_key_sets_match = bool(
            len(denominator_ladder) == len(PRECISIONS)
            and bool(denominator_ladder[0])
            and all(
                set(denominator_ladder[index])
                == set(denominator_ladder[0])
                for index in range(1, len(denominator_ladder))
            )
        )
        denominator_nested = bool(
            denominator_key_sets_match
            and all(
                cinterval_nested(
                    denominator_ladder[k][identifier],
                    denominator_ladder[k + 1][identifier],
                )
                for identifier in denominator_ladder[0]
                for k in range(len(denominator_ladder) - 1)
            )
        )
        comparison_pairs = []
        for ladder_index in range(len(PRECISIONS) - 1):
            quantity_pair = {
                label: cinterval_nested(
                    probe_ladder[ladder_index][label],
                    probe_ladder[ladder_index + 1][label],
                )
                for label in probe_names
            }
            segment_pair = bool(
                partition_id_ladder[ladder_index]
                == partition_id_ladder[ladder_index + 1]
                and len(segment_record_ladder[ladder_index])
                == len(segment_record_ladder[ladder_index + 1])
                and all(
                    cinterval_nested(
                        segment_record_ladder[ladder_index][segment_index][field],
                        segment_record_ladder[ladder_index + 1][segment_index][
                            field
                        ],
                    )
                    for segment_index in range(
                        len(segment_record_ladder[ladder_index])
                    )
                    for field in (
                        "_center_value",
                        "_derivative_hull",
                        "_offset",
                        "_image",
                        "_rotated_image",
                        "_start_value",
                        "_end_value",
                        "_endpoint_ratio",
                    )
                )
                and all(
                    interval_nested(
                        segment_record_ladder[ladder_index][segment_index][field],
                        segment_record_ladder[ladder_index + 1][segment_index][
                            field
                        ],
                    )
                    for segment_index in range(
                        len(segment_record_ladder[ladder_index])
                    )
                    for field in (
                        "_rotated_argument",
                        "_endpoint_increment",
                    )
                )
            )
            newton_pair = {
                label: bool(
                    ball_ladder[ladder_index][label] is not None
                    and ball_ladder[ladder_index + 1][label] is not None
                    and cinterval_nested(
                        ball_ladder[ladder_index][label],
                        ball_ladder[ladder_index + 1][label],
                    )
                )
                for label in ball_nesting
            }
            denominator_pair = bool(
                set(denominator_ladder[ladder_index])
                == set(denominator_ladder[ladder_index + 1])
                and all(
                    cinterval_nested(
                        denominator_ladder[ladder_index][identifier],
                        denominator_ladder[ladder_index + 1][identifier],
                    )
                    for identifier in denominator_ladder[ladder_index]
                )
            )
            comparison_pairs.append(
                {
                    "outer_precision": PRECISIONS[ladder_index],
                    "inner_precision": PRECISIONS[ladder_index + 1],
                    "total_variation_nested": interval_nested(
                        raw_totals[ladder_index],
                        raw_totals[ladder_index + 1],
                    ),
                    "partition_identical": (
                        rows[str(PRECISIONS[ladder_index])][
                            "boundary_winding"
                        ]["partition_sha256"]
                        == rows[str(PRECISIONS[ladder_index + 1])][
                            "boundary_winding"
                        ]["partition_sha256"]
                    ),
                    "segment_enclosures_nested": segment_pair,
                    "quantity_enclosure_nesting": quantity_pair,
                    "newton_ball_nesting": newton_pair,
                    "coefficient_denominators_nested": denominator_pair,
                }
            )
        nesting[name] = {
            "enclosures_nested_with_precision": bool(total_nested),
            "per_quantity_enclosure_nesting": per_quantity,
            "quantity_key_sets_match": bool(
                all(
                    set(quantities) == set(probe_ladder[0])
                    for quantities in probe_ladder
                )
            ),
            "quantity_enclosures_all_nested": bool(
                all(
                    set(quantities) == set(probe_ladder[0])
                    for quantities in probe_ladder
                )
                and all(per_quantity.values())
            ),
            "per_segment_enclosure_nesting": segment_nesting,
            "newton_ball_nesting": ball_nesting,
            "coefficient_denominator_nesting": {
                "records": (
                    len(denominator_ladder[0])
                    if denominator_ladder
                    else 0
                ),
                "key_sets_match": denominator_key_sets_match,
                "all_nested": denominator_nested,
            },
            "comparison_pairs": comparison_pairs,
        }
        all_certified = (
            all_certified
            and total_nested
            and all(per_quantity.values())
            and segment_nesting["all_nested"]
            and all(ball_nesting.values())
            and denominator_nested
            and all(
                pair["partition_identical"]
                and pair["total_variation_nested"]
                and pair["segment_enclosures_nested"]
                and all(pair["quantity_enclosure_nesting"].values())
                and all(pair["newton_ball_nesting"].values())
                and pair["coefficient_denominators_nested"]
                for pair in comparison_pairs
            )
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
        "claim_scope": {
            "certified_object": (
                "scalar_transverse_inverse_propagator_on_declared_"
                "algebraic_chart"
            ),
            "declared_chart_only": True,
            "continuation_identity_independently_certified": False,
            "standard_second_sheet_identification_certified": False,
            "full_matrix_rank_n_minus_1_laurent_certified": False,
            "issue_593_pole_laurent_row_discharged": False,
            "finite_delta_diagnostics_are_noncertifying": True,
            "independent_numerical_replay_certified": False,
            "engine_inverse_propagator_convention": (
                "G(s)=s-m_tree^2-Pi_engine(s)"
            ),
            "theorem_self_energy_sign_bridge_certified": False,
            "issue_593_precision_ladder_row_discharged": False,
            "issue_593_independent_third_verifier_row_discharged": False,
            "issue_593_full_acceptance_satisfied": False,
        },
        "source_subject": {
            "relative_path": "outputs/fj_direct_vector_blocks.json",
            "artifact_role": "direct_FJ_vector_blocks_input",
            "schema": vector.get("schema"),
            "target": vector.get("target"),
            "units": vector.get("units"),
            "bytes": len(vector_raw),
            "sha256": sha256_bytes(vector_raw),
        },
        "pins": {
            "vector_blocks_sha256": sha256_bytes(vector_raw),
            "diagnostic_producer_sha256": sha256_bytes(
                (ROOT / "producers" / "wz_pole_receipts.py").read_bytes()
            ),
            "interval_arithmetic_module_sha256": sha256_bytes(
                (ROOT / "producers" / "complex_interval.py").read_bytes()
            ),
            "principal_zero_exclusion_producer_sha256": sha256_bytes(
                (ROOT / "producers" / "certified_wz_contours.py").read_bytes()
            ),
            "producer_module_sha256": sha256_bytes(Path(__file__).read_bytes()),
            "checker_module_sha256": sha256_bytes(
                (ROOT / "checkers" / "check_certified_second_sheet_poles.py").read_bytes()
            ),
            "principal_zero_exclusion_receipt_sha256": zero_exclusion_sha,
        },
        "fixture": {k: str(f) for k, f in wzp.FIXTURE.items()},
        "ckm_fixture": {
            f"V{i}{j}": str(wzp.CKM_FIXTURE[(i, j)])
            for i in (1, 2, 3)
            for j in (1, 2, 3)
        },
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
            "newton_strict_interior_required": True,
            "interval_encoding": "exact_dyadic_rational_endpoints",
        },
        "sheet_statement": (
            "the pole boxes lie strictly in the open lower half plane on a "
            "declared algebraic chart with an explicit per-channel sheet "
            "vector; in the W chart W-gamma has the declared action principal, "
            "meaning no added chart-relative correction rather than a "
            "certified physical-sheet assignment; "
            "the one-mass base chart is mass-exchange symmetric and every "
            "correction is relative to its recorded base chart; "
            "finite-delta probes are non-certifying diagnostics, and no "
            "independent continuation identity or physical/unique resonance "
            "sheet is certified"
        ),
        "results": results,
        "precision_nesting": nesting,
        "promotion": {
            "scalar_pole_on_declared_algebraic_chart_certified": all_certified,
            "simple_scalar_root_certified": all_certified,
            "scalar_laurent_denominator_ball_certified": all_certified,
            "scalar_residue_ball_certified": all_certified,
            "continuation_identity_independently_certified": False,
            "standard_second_sheet_identification_certified": False,
            "theorem_self_energy_sign_bridge_certified": False,
            "full_matrix_rank_n_minus_1_laurent_certified": False,
            "issue_593_precision_ladder_row_discharged": False,
            "issue_593_independent_third_verifier_row_discharged": False,
            "issue_593_pole_laurent_row_discharged": False,
            "issue_593_full_acceptance_satisfied": False,
            "bmhv_restoration_certified": False,
            "physical_current_claim": False,
            "oph_native": False,
            "unit_claim": False,
            "unitarity_claim": False,
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
    try:
        payload = build_receipt()
        if args.verify:
            stored = json.loads(OUT_PATH.read_text(encoding="utf-8"))
            if stored != payload:
                print("SECOND_SHEET_POLE_DRIFT", file=sys.stderr)
                return 1
            print("SECOND_SHEET_POLE_VERIFIED")
            return 0
        OUT_PATH.write_text(
            json.dumps(payload, sort_keys=True, indent=1) + "\n",
            encoding="utf-8",
        )
        print(json.dumps({"status": payload["status"]}))
        return 0 if payload["status"] == STATUS_CERTIFIED else 1
    finally:
        if CHECKPOINT_PATH.exists():
            CHECKPOINT_PATH.unlink()


if __name__ == "__main__":
    raise SystemExit(main())
