#!/usr/bin/env python3
"""Certified principal-sheet W and Z zero exclusion on the synthetic fixture.

The verifier upgrades the sampled boundary diagnostic to a directed
interval certificate on declared off-axis boxes.  Every quantity is a
rectangular complex interval with outward rounding, every logarithm
and square root carries a sheet gate refusing the principal cut, and
the argument principle is applied with interval argument chaining, so
a reported winding is a proved root count for the holomorphic
principal-sheet function on the declared box.

What is certified, exactly:

* directed complex-interval evaluation of the compiled one-loop
  transverse blocks, fixture-exact coefficients, on boxes that lie
  strictly in the upper half plane, so no evaluation meets the
  physical cut on the real axis;
* interior holomorphy: on a declared grid refinement of each box,
  every logarithm and square-root argument of every loop function
  excludes the principal cut and every denominator excludes zero;
* boundary exclusion: every boundary segment enclosure excludes zero
  with the declared argument-width gate, subdividing adaptively to
  the declared depth cap;
* the winding number by rigorous interval argument chaining, with the
  total variation enclosure inside the declared tolerance of an
  integer multiple of two pi;
* enclosure nesting across the declared precision ladder.

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

SCHEMA = "oph.certified_wz_contours.v1"
STATUS_CERTIFIED = "PRINCIPAL_SHEET_ZERO_EXCLUSION_CERTIFIED"
STATUS_FAILED = "CERTIFICATION_INCOMPLETE"

PRECISIONS = (128, 192, 256)
INITIAL_SEGMENTS_PER_EDGE = 8
MAX_SUBDIVISION_DEPTH = 12
HOLOMORPHY_GRID = 8
ARG_WIDTH_GATE_NUM = Fraction(157, 100)
HOLOMORPHY_ARG_GATE_NUM = Fraction(31, 5)
WINDING_TOLERANCE_NUM = Fraction(157, 100)

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


def _root_term(x: CInterval, gate: Any) -> CInterval:
    one = CInterval.from_fraction(1)
    return x * ((x - one) / x).log(gate) - (x - one).log(gate)


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
    four = CInterval.from_fraction(4)
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
    a = s
    b = -(s + CInterval.from_fraction(m1 - m2))
    c = CInterval.from_fraction(m1)
    disc = b * b - four * a * c
    sq = disc.sqrt(gate)
    x1 = (-b + sq) / (two * a)
    x2 = (-b - sq) / (two * a)
    return total + _root_term(x1, gate) + _root_term(x2, gate)


class IntervalEvaluator:
    """Interval evaluation of one compiled transverse block."""

    def __init__(self, compiled: list[dict[str, Any]], gate: Any) -> None:
        self.compiled = compiled
        self.gate = gate
        mu2 = wzp.FIXTURE["mu_ren2"]
        self.mu2 = mu2
        self.loop_cache: dict[tuple, CInterval] = {}

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

    def inverse_propagator(
        self, s: CInterval, s_key: tuple, tree_mass: Fraction
    ) -> CInterval:
        return s - CInterval.from_fraction(tree_mass) - self.transverse(s, s_key)


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
    evaluator: IntervalEvaluator,
    box: dict[str, tuple[Fraction, Fraction]],
    tree_mass: Fraction,
    grid: int,
) -> dict[str, Any]:
    """Cut exclusion and denominator exclusion on a cover of the box.

    Holomorphy needs only cut and zero exclusion, so the evaluator is
    run under the declared loose argument gate; the tight gate belongs
    to the boundary winding chain.  A refusing cell is subdivided into
    quarters up to the declared depth cap, so enclosure overestimation
    on coarse cells is refined away rather than reported as failure."""

    re_lo, re_hi = box["re"]
    im_lo, im_hi = box["im"]
    worklist = []
    for i in range(grid):
        for j in range(grid):
            worklist.append(
                (
                    re_lo + (re_hi - re_lo) * Fraction(i, grid),
                    re_lo + (re_hi - re_lo) * Fraction(i + 1, grid),
                    im_lo + (im_hi - im_lo) * Fraction(j, grid),
                    im_lo + (im_hi - im_lo) * Fraction(j + 1, grid),
                    0,
                )
            )
    certified_cells = 0
    max_depth_used = 0
    while worklist:
        a, b, c, d, depth = worklist.pop()
        cell = CInterval.box(a, b, c, d)
        try:
            evaluator.loop_cache.clear()
            evaluator.inverse_propagator(
                cell, ("cell", str(a), str(c), depth), tree_mass
            )
            certified_cells += 1
            continue
        except SheetError:
            pass
        if depth >= MAX_SUBDIVISION_DEPTH:
            return {
                "grid": grid,
                "certified_cells": certified_cells,
                "holomorphic_on_box": False,
                "reason": "interior subdivision depth cap reached",
            }
        mid_re = (a + b) / 2
        mid_im = (c + d) / 2
        max_depth_used = max(max_depth_used, depth + 1)
        for cell_bounds in (
            (a, mid_re, c, mid_im),
            (mid_re, b, c, mid_im),
            (a, mid_re, mid_im, d),
            (mid_re, b, mid_im, d),
        ):
            worklist.append(cell_bounds + (depth + 1,))
    return {
        "grid": grid,
        "certified_cells": certified_cells,
        "max_depth_used": max_depth_used,
        "holomorphic_on_box": True,
    }


def certify_winding(
    evaluator: IntervalEvaluator,
    box: dict[str, tuple[Fraction, Fraction]],
    tree_mass: Fraction,
    gate: Any,
) -> dict[str, Any]:
    """Rigorous argument-principle winding with endpoint chaining.

    Each boundary segment enclosure must exclude zero with argument
    width below the gate, which confines the image of the segment to a
    cone narrower than pi, so the argument cannot wrap within the
    segment.  The argument increment of the segment is then the unique
    representative of the endpoint argument difference with modulus
    below the gate, evaluated on point rectangles whose enclosure
    widths carry rounding only.  The winding is the summed increment
    enclosure divided by two pi, certified when the residual enclosure
    lies inside the declared tolerance."""

    worklist = [
        (segment, 0)
        for segment in _boundary_segments(box, INITIAL_SEGMENTS_PER_EDGE)
    ]
    certified: list[tuple] = []
    max_depth_used = 0
    while worklist:
        (start, end), depth = worklist.pop(0)
        hull = _segment_hull(start, end)
        try:
            evaluator.loop_cache.clear()
            image = evaluator.inverse_propagator(
                hull, ("seg", str(start), str(end)), tree_mass
            )
            ok = (not image.contains_zero()) and image.arg().delta <= gate
        except SheetError:
            ok = False
        if ok:
            certified.append((start, end))
            continue
        if depth >= MAX_SUBDIVISION_DEPTH:
            return {
                "certified": False,
                "reason": "subdivision depth cap reached",
                "segments": len(certified),
            }
        mid = ((start[0] + end[0]) / 2, (start[1] + end[1]) / 2)
        worklist.insert(0, ((mid, end), depth + 1))
        worklist.insert(0, ((start, mid), depth + 1))
        max_depth_used = max(max_depth_used, depth + 1)

    point_args: dict[tuple, Any] = {}

    def endpoint_arg(point: tuple) -> Any:
        cached = point_args.get(point)
        if cached is not None:
            return cached
        rect = CInterval.box(point[0], point[0], point[1], point[1])
        evaluator.loop_cache.clear()
        value = evaluator.inverse_propagator(
            rect, ("pt", str(point)), tree_mass
        )
        angle = value.arg()
        point_args[point] = angle
        return angle

    two_pi = iv.pi * iv.mpf(2)
    total = iv.mpf(0)
    for start, end in certified:
        difference = endpoint_arg(end) - endpoint_arg(start)
        candidates = [
            difference,
            difference - two_pi,
            difference + two_pi,
        ]
        admissible = [
            candidate
            for candidate in candidates
            if abs(candidate.a) <= gate and abs(candidate.b) <= gate
        ]
        if len(admissible) != 1:
            return {
                "certified": False,
                "reason": "segment increment lift ambiguous",
                "segments": len(certified),
            }
        total = total + admissible[0]
    winding = int(round(float((total.a + total.b) / 2) / float(two_pi.b)))
    tolerance = _iv_fraction(WINDING_TOLERANCE_NUM)
    residual = total - two_pi * iv.mpf(winding)
    within = bool(residual.a > -tolerance.a and residual.b < tolerance.a)
    return {
        "certified": within,
        "winding": winding,
        "segments": len(certified),
        "max_depth_used": max_depth_used,
        "total_variation_interval": [str(total.a), str(total.b)],
        "reason": None if within else "total variation outside tolerance",
        "_raw_total": total,
    }


def certify_particle(
    name: str,
    compiled: list[dict[str, Any]],
    tree_mass: Fraction,
    precision: int,
) -> dict[str, Any]:
    iv.prec = precision
    tight_gate = _iv_fraction(ARG_WIDTH_GATE_NUM).a
    loose_gate = _iv_fraction(HOLOMORPHY_ARG_GATE_NUM).a
    box = BOXES[name]
    holomorphy = certify_interior_holomorphy(
        IntervalEvaluator(compiled, loose_gate), box, tree_mass,
        HOLOMORPHY_GRID,
    )
    winding = certify_winding(
        IntervalEvaluator(compiled, tight_gate), box, tree_mass, tight_gate
    )
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
    compiled = {
        "W": wzp.compile_block(wzp.block_transverse(vector, "WpWm")),
        "Z": wzp.compile_block(wzp.block_transverse(vector, "ZZ")),
    }

    results: dict[str, Any] = {}
    nesting: dict[str, Any] = {}
    all_certified = True
    for name in ("W", "Z"):
        rows = {}
        raw_totals = []
        for precision in PRECISIONS:
            row = certify_particle(name, compiled[name], masses[name], precision)
            raw = row["boundary_winding"].pop("_raw_total", None)
            rows[str(precision)] = row
            all_certified = all_certified and row["zero_exclusion_certified"]
            if raw is not None:
                raw_totals.append(raw)
        nested = (
            all(
                raw_totals[k + 1].a >= raw_totals[k].a
                and raw_totals[k + 1].b <= raw_totals[k].b
                for k in range(len(raw_totals) - 1)
            )
            if len(raw_totals) == len(PRECISIONS)
            else False
        )
        nesting[name] = {
            "enclosures_nested_with_precision": bool(nested),
            "comparison": "exact endpoint comparison of the raw enclosures",
        }
        results[name] = rows

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
        "serialized_gates": {
            "precisions_bits": list(PRECISIONS),
            "initial_segments_per_edge": INITIAL_SEGMENTS_PER_EDGE,
            "max_subdivision_depth": MAX_SUBDIVISION_DEPTH,
            "holomorphy_grid": HOLOMORPHY_GRID,
            "arg_width_gate": str(ARG_WIDTH_GATE_NUM),
            "holomorphy_arg_gate": str(HOLOMORPHY_ARG_GATE_NUM),
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
    print(json.dumps({"status": payload["status"]}))
    return 0 if payload["status"] == STATUS_CERTIFIED else 1


if __name__ == "__main__":
    raise SystemExit(main())
