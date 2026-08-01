#!/usr/bin/env python3
"""Independent verifier for the issue-654 v3 persistence packet.

This verifier deliberately does not import the v3 producer.  It rebuilds the
exact moment template, Taylor-tail bounds, and critical-axis separation from
the immutable Cartesian carrier implementation, then checks custody and every
open promotion boundary in the serialized receipt.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from math import factorial
from pathlib import Path
from typing import Any, Iterable

import a5_multipole_fixed_point_certificate as base


HERE = Path(__file__).resolve().parent
RECEIPT_PATH = HERE / "runtime" / "a5_multipole_persistence_receipt_v3.json"
PARENT_PATH = HERE / "runtime" / "a5_multipole_fixed_point_receipt_v2.json"

EXPECTED_SCHEMA = "oph.a5_multipole_persistence_receipt.v3"
EXPECTED_STATUS = (
    "EXACT_THROUGH_EIGHTH_ORDER_I6_TEMPLATE__"
    "FULL_COSINE_X10_PLUS_TAIL_BOUNDS__"
    "GLOBAL_PERSISTENCE_COVER_OPEN__PHYSICAL_MAP_OPEN"
)
PARENT_SCHEMA = "oph.a5_multipole_fixed_point_receipt.v2"
PARENT_STATUS = (
    "EXACT_A5_FINGERPRINT_CORE__QUANTITATIVE_PERSISTENCE_OPEN__"
    "PHYSICAL_MAP_OPEN"
)

Q5 = base.Q5
ZERO = base.ZERO
ONE = base.ONE


class PersistenceVerificationError(ValueError):
    """The independent verifier rejected the v3 receipt."""


def check(condition: bool, message: str) -> None:
    if not condition:
        raise PersistenceVerificationError(message)


def self_hash(value: dict[str, Any]) -> str:
    body = {key: item for key, item in value.items() if key != "receipt_sha256"}
    return "sha256:" + hashlib.sha256(base.canonical_json_bytes(body)).hexdigest()


def load_object(path: Path) -> tuple[bytes, dict[str, Any]]:
    raw = path.read_bytes()
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as error:
        raise PersistenceVerificationError(f"invalid JSON: {path}") from error
    check(isinstance(value, dict), f"JSON root is not an object: {path}")
    return raw, value


def q5_sum(values: Iterable[Q5]) -> Q5:
    out = ZERO
    for value in values:
        out = base.q5_add(out, value)
    return out


def dot(x: tuple[Q5, Q5, Q5], y: tuple[Q5, Q5, Q5]) -> Q5:
    return q5_sum(base.q5_mul(x[i], y[i]) for i in range(3))


def cross(
    x: tuple[Q5, Q5, Q5], y: tuple[Q5, Q5, Q5]
) -> tuple[Q5, Q5, Q5]:
    return (
        base.q5_sub(base.q5_mul(x[1], y[2]), base.q5_mul(x[2], y[1])),
        base.q5_sub(base.q5_mul(x[2], y[0]), base.q5_mul(x[0], y[2])),
        base.q5_sub(base.q5_mul(x[0], y[1]), base.q5_mul(x[1], y[0])),
    )


def add_rays(
    vertices: list[tuple[Q5, Q5, Q5]], indices: Iterable[int]
) -> tuple[Q5, Q5, Q5]:
    selected = tuple(indices)
    return tuple(
        q5_sum(vertices[index][axis] for index in selected) for axis in range(3)
    )  # type: ignore[return-value]


def decompose_on_i6(polynomial: base.Poly, i6: base.Poly) -> tuple[Q5, Q5]:
    reduced = base.p_reduce_sphere(polynomial)
    i6_reduced = base.p_reduce_sphere(i6)
    probe = next(mono for mono in i6_reduced if sum(mono) > 0)
    amplitude = base.q5_div(reduced.get(probe, ZERO), i6_reduced[probe])
    residue = base.p_add(
        reduced,
        base.p_scale_q5(i6_reduced, base.q5_neg(amplitude)),
    )
    constant = residue.pop((0, 0, 0), ZERO)
    check(base.p_is_zero(residue), "independent moment leaves constant-plus-I6 line")
    return constant, amplitude


def expected_template() -> dict[str, Any]:
    vertices = base.cartesian_vertices()
    check(len(vertices) == 12, "independent carrier count drift")
    check(
        all(dot(vertex, vertex) == base.NORM_SQ for vertex in vertices),
        "independent carrier norm drift",
    )
    i6 = base.build_cartesian_frame()["_i6_poly_object"]
    constant6, amplitude6 = decompose_on_i6(
        base.normalized_moment(vertices, 6, base.NORM_SQ), i6
    )
    constant8, amplitude8 = decompose_on_i6(
        base.normalized_moment(vertices, 8, base.NORM_SQ), i6
    )
    check(constant6 == base.q5(Fraction(12, 7)), "sixth-moment constant drift")
    check(amplitude6 == base.q5(Fraction(64, 175)), "sixth-moment I6 drift")
    check(constant8 == base.q5(Fraction(4, 3)), "eighth-moment constant drift")
    check(amplitude8 == base.q5(Fraction(256, 375)), "eighth-moment I6 drift")
    leading = amplitude6[0] / factorial(6)
    eighth = amplitude8[0] / factorial(8)
    check(leading == Fraction(4, 7875), "independent x6 coefficient drift")
    check(eighth == Fraction(2, 118125), "independent x8 coefficient drift")
    check(leading / eighth == 30, "independent coefficient ratio drift")
    return {
        "kernel": "Q_x(n) = sum_i [1 - cos(x u_i.n)], x = |a k|",
        "kernel_status": "DECLARED_EQUAL_WEIGHT_COSINE_BRANCH",
        "declared_branch_premise": True,
        "source_selected": False,
        "architecture_forced": False,
        "physical_source_selection_owner": 655,
        "declared_range": "0 < x <= 1",
        "domain_assumptions": {
            "angular_argument": "n in S^2 with ||n|| = 1",
            "carrier_directions": "u_i in S^2 with ||u_i|| = 1, i = 1,...,12",
            "carrier_unit_norm_verified": True,
            "normalization_method": (
                "the exact Cartesian vertices have common squared norm "
                "5/2+1/2*sqrt5; each direction is divided by the square root "
                "of that common squared norm"
            ),
        },
        "moment_decomposition": {
            "sum_i (u_i.n)^6": "12/7 + (64/175) I6(n)",
            "sum_i (u_i.n)^8": "4/3 + (256/375) I6(n)",
        },
        "anisotropic_part_through_eighth_order": (
            "A(x) I6(n), A(x) = 2 x^6 (30 - x^2) / 118125"
        ),
        "x6_coefficient": str(leading),
        "x8_subtracted_coefficient": str(eighth),
        "positivity_certificate": {
            "lower_factor_on_declared_range": "30 - x^2 >= 29",
            "strictly_positive": True,
        },
        "normalization_role": (
            "multiplication by the positive A(x) and subtraction of an "
            "isotropic constant do not change stationary directions"
        ),
    }


def expected_tail(row: dict[str, Any]) -> dict[str, Any]:
    try:
        x_max = Fraction(row["x_max"])
    except (KeyError, TypeError, ValueError, ZeroDivisionError) as error:
        raise PersistenceVerificationError("invalid serialized x_max") from error
    check(0 < x_max <= 1, "serialized tail range leaves 0 < x_max <= 1")
    port_count = 12
    gradient_geometric_factor = Fraction(110, 109)
    hessian_geometric_factor = Fraction(90, 89)
    lower = Fraction(2) * (30 - x_max * x_max) / 118125
    gradient_raw = (
        Fraction(port_count, factorial(9))
        * gradient_geometric_factor
        * x_max**10
    )
    hessian_raw = (
        Fraction(port_count, factorial(8))
        * hessian_geometric_factor
        * x_max**10
    )
    denominator = lower * x_max**6
    gradient = gradient_raw / denominator
    intrinsic_hessian = (hessian_raw + gradient_raw) / denominator
    if x_max == 1:
        check(gradient == Fraction(6875, 101152), "independent C1 value drift")
        check(
            intrinsic_hessian == Fraction(383125, 562658),
            "independent C2 value drift",
        )
    check(gradient < Fraction(7, 100), "independent C1 envelope fails")
    check(intrinsic_hessian < Fraction(7, 10), "independent C2 envelope fails")
    return {
        "x_max": str(x_max),
        "range": f"0 < x <= {x_max}",
        "tail_start": "x^10",
        "normalization_lower_bound": (
            f"A(x) >= {lower} x^6 for 0 < x <= {x_max}"
        ),
        "gradient_tail_derivation": (
            "12 x^10 [1/9! + 1/11! + ...] / A(x), with the "
            "bracket bounded by (110/109)/9!"
        ),
        "intrinsic_hessian_tail_derivation": (
            "12 x^10 [(90/89)/8! + (110/109)/9!] / A(x); "
            "the second summand bounds the spherical radial correction"
        ),
        "normalized_C1_gradient_bound": str(gradient),
        "normalized_C2_intrinsic_hessian_bound": str(intrinsic_hessian),
        "derivation_parameters": {
            "port_count": 12,
            "gradient_successive_ratio_upper": "1/110",
            "gradient_geometric_sum_factor": "110/109",
            "euclidean_hessian_successive_ratio_upper": "1/90",
            "euclidean_hessian_geometric_sum_factor": "90/89",
            "normalization_lower_coefficient": str(lower),
            "radial_hessian_correction_included": True,
            "radial_hessian_correction_bound": "the normalized gradient tail",
        },
        "clean_envelopes": {
            "gradient_strictly_below": "7/100",
            "intrinsic_hessian_strictly_below": "7/10",
        },
        "arithmetic": "exact Fraction arithmetic; no transcendental evaluation",
    }


def expected_axis_separation() -> dict[str, Any]:
    vertices = base.cartesian_vertices()
    inv_norm = base.q5_div(ONE, base.NORM_SQ)

    def unit_dot(i: int, j: int) -> Q5:
        return base.q5_mul(dot(vertices[i], vertices[j]), inv_norm)

    edges = [
        (i, j)
        for i in range(12)
        for j in range(i + 1, 12)
        if unit_dot(i, j) == base.INV_SQRT5
    ]
    edge_set = set(edges)
    faces = [
        (i, j, k)
        for i in range(12)
        for j in range(i + 1, 12)
        for k in range(j + 1, 12)
        if (i, j) in edge_set and (i, k) in edge_set and (j, k) in edge_set
    ]
    check(len(edges) == 30 and len(faces) == 20, "independent cell census drift")
    typed_rays = (
        [("vertex", vertex) for vertex in vertices]
        + [("edge", add_rays(vertices, edge)) for edge in edges]
        + [("face", add_rays(vertices, face)) for face in faces]
    )
    axes: list[tuple[str, tuple[Q5, Q5, Q5]]] = []
    for orbit, ray in typed_rays:
        parallel = [
            entry for entry in axes if cross(ray, entry[1]) == (ZERO, ZERO, ZERO)
        ]
        if parallel:
            check(len(parallel) == 1, "independent axis deduplication ambiguity")
            check(parallel[0][0] == orbit, "independent orbit-axis collision")
            continue
        axes.append((orbit, ray))
    counts = {
        orbit: sum(kind == orbit for kind, _ in axes)
        for orbit in ("vertex", "face", "edge")
    }
    check(counts == {"vertex": 6, "face": 10, "edge": 15}, "axis census drift")
    best: Q5 | None = None
    pair: tuple[int, int] | None = None
    for left, (_, x) in enumerate(axes):
        for right in range(left + 1, len(axes)):
            y = axes[right][1]
            cosine_squared = base.q5_div(
                base.q5_pow(dot(x, y), 2),
                base.q5_mul(dot(x, x), dot(y, y)),
            )
            check(
                base.q5_sign(cosine_squared) >= 0
                and base.q5_sign(base.q5_sub(ONE, cosine_squared)) > 0,
                "independent distinct axes coincide",
            )
            if best is None or base.q5_sign(base.q5_sub(cosine_squared, best)) > 0:
                best = cosine_squared
                pair = (left, right)
    check(best == base.q5(Fraction(1, 2), Fraction(1, 6)), "axis cosine drift")
    minimum = base.q5_sub(ONE, best)
    check(
        minimum == base.q5(Fraction(1, 2), Fraction(-1, 6)),
        "axis sine separation drift",
    )
    radius = base.q5(Fraction(1, 64))
    check(
        base.q5_sign(base.q5_sub(minimum, base.q5_scale(radius, 4))) > 0,
        "independent neighborhood separation fails",
    )
    check(pair is not None, "independent maximizing pair missing")
    return {
        "unoriented_axis_count": 31,
        "oriented_critical_point_count": 62,
        "axis_counts_by_orbit": counts,
        "maximum_squared_axis_cosine": "1/2+1/6*sqrt5",
        "minimum_squared_axis_sine": "1/2-sqrt5/6",
        "maximizing_pair_indices": list(pair),
        "maximizing_pair_orbit_types": [axes[pair[0]][0], axes[pair[1]][0]],
        "declared_local_neighborhood": "sin^2(angle to an axis) <= 1/64",
        "separation_check": (
            "four times the neighborhood sine-squared radius is strictly "
            "below the exact minimum axis sine-squared separation"
        ),
        "promotion_from_separation_permitted": False,
    }


def verify_receipt(
    receipt_path: Path = RECEIPT_PATH,
    parent_path: Path = PARENT_PATH,
) -> dict[str, Any]:
    raw, receipt = load_object(receipt_path)
    check(raw == base.canonical_json_bytes(receipt), "v3 receipt is not canonical JSON")
    check(receipt.get("schema") == EXPECTED_SCHEMA, "v3 schema drift")
    check(receipt.get("status") == EXPECTED_STATUS, "v3 status drift")
    check(receipt.get("issue") == 654, "v3 issue drift")
    check(receipt.get("receipt_sha256") == self_hash(receipt), "v3 self-digest drift")

    parent_raw, parent = load_object(parent_path)
    check(parent.get("schema") == PARENT_SCHEMA, "parent schema drift")
    check(parent.get("status") == PARENT_STATUS, "parent status drift")
    check(parent.get("receipt_sha256") == self_hash(parent), "parent self-digest drift")
    expected_relation = {
        "schema": PARENT_SCHEMA,
        "status": PARENT_STATUS,
        "path": "code/a5_fingerprint/runtime/a5_multipole_fixed_point_receipt_v2.json",
        "reason": "append-only quantitative extension; v1 and v2 remain immutable",
    }
    check("supersedes" not in receipt, "v3 must extend, not supersede, v2")
    check(receipt.get("extends") == expected_relation, "v3 extension relation drift")
    expected_pin = {
        "path": "code/a5_fingerprint/runtime/a5_multipole_fixed_point_receipt_v2.json",
        "bytes": len(parent_raw),
        "sha256": base.tagged_sha256(parent_raw),
        "receipt_sha256": parent["receipt_sha256"],
    }
    check(receipt.get("parent_pin") == expected_pin, "v3 parent pin drift")

    check(
        receipt.get("through_eighth_order_i6_template") == expected_template(),
        "through-eighth-order I6 template mismatch",
    )
    tail = receipt.get("normalized_tail_bounds")
    check(isinstance(tail, dict), "normalized tail section missing")
    check(tail == expected_tail(tail), "normalized tail-bound mismatch")
    check(
        receipt.get("critical_axis_separation") == expected_axis_separation(),
        "critical-axis separation mismatch",
    )

    expected_persistence = {
        "declared_full_cosine_kernel_for_tail_bounds": True,
        "through_eighth_order_i6_template": True,
        "exact_C1_C2_tail_bounds": True,
        "exact_critical_axis_separation": True,
        "global_interval_gradient_cover": False,
        "local_interval_newton_uniqueness_boxes": False,
        "finite_exactly_62_persistence_range": False,
        "status": "PARTIAL__GLOBAL_COVER_AND_LOCAL_UNIQUENESS_OPEN",
        "claim_boundary": (
            "the packet bounds the full-cosine tail on 0 < |a k| <= 1 "
            "but does not infer an exactly-62 critical-point theorem until "
            "the global off-neighborhood cover and local uniqueness boxes "
            "are both certified"
        ),
    }
    check(
        receipt.get("quantitative_persistence") == expected_persistence,
        "quantitative persistence boundary drift",
    )
    check(
        receipt.get("comparison_boundary")
        == {
            "public_measurement_read": False,
            "comparison_permitted": False,
            "physical_map_open": True,
        },
        "comparison or physical-map boundary drift",
    )
    controls = receipt.get("fail_closed_controls", {}).get("controls")
    check(isinstance(controls, list) and len(controls) == 11, "control census drift")
    check(all(row.get("detector_fired") is True for row in controls), "control drift")
    check(
        receipt.get("fail_closed_controls", {}).get("all_detectors_fired") is True,
        "control aggregate drift",
    )
    return {
        "status": "PASS",
        "receipt_sha256": receipt["receipt_sha256"],
        "parent_sha256": expected_pin["sha256"],
        "independent_template_recomputed": True,
        "independent_tail_bounds_recomputed": True,
        "independent_axis_separation_recomputed": True,
        "promotion_boundaries_fail_closed": True,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=RECEIPT_PATH)
    parser.add_argument("--parent", type=Path, default=PARENT_PATH)
    args = parser.parse_args(argv)
    print(json.dumps(verify_receipt(args.receipt, args.parent), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
