#!/usr/bin/env python3
"""Append-only quantitative persistence packet for issue 654.

The v1 fingerprint and v2 hardening receipts are immutable parents.  This
producer adds the exact part of the finite-kernel persistence calculation
without changing either parent:

* it identifies the equal-weight cosine kernel's anisotropic template through
  eighth order as one positive multiple of the normalized ``I6`` for
  ``0 < |a k| <= 1``;
* it bounds the normalized gradient and intrinsic Hessian of the remaining
  Taylor tail using exact rational arithmetic; and
* it computes the exact separation of the 31 antipodal critical axes.

The global interval cover and the local interval-Newton uniqueness boxes are
kept fail-closed.  Until both are present, this packet does not assert an
exactly-62 theorem for the non-polynomial cosine kernel.
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
import a5_multipole_fixed_point_hardening_certificate as hardening


HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
PARENT_PATH = RUNTIME / "a5_multipole_fixed_point_receipt_v2.json"
RECEIPT_PATH = RUNTIME / "a5_multipole_persistence_receipt_v3.json"

SCHEMA = "oph.a5_multipole_persistence_receipt.v3"
STATUS = (
    "EXACT_THROUGH_EIGHTH_ORDER_I6_TEMPLATE__"
    "FULL_COSINE_X10_PLUS_TAIL_BOUNDS__"
    "GLOBAL_PERSISTENCE_COVER_OPEN__PHYSICAL_MAP_OPEN"
)

Q5 = base.Q5
ZERO = base.ZERO
ONE = base.ONE
q5 = base.q5
q5_add = base.q5_add
q5_sub = base.q5_sub
q5_mul = base.q5_mul
q5_div = base.q5_div
q5_pow = base.q5_pow
q5_scale = base.q5_scale
q5_neg = base.q5_neg
q5_sign = base.q5_sign
q5_str = base.q5_str
require = base.require


def artifact_self_hash(value: dict[str, Any]) -> str:
    body = {key: item for key, item in value.items() if key != "receipt_sha256"}
    return "sha256:" + hashlib.sha256(base.canonical_json_bytes(body)).hexdigest()


def load_v2_parent(path: Path = PARENT_PATH) -> tuple[bytes, dict[str, Any]]:
    raw = path.read_bytes()
    try:
        parent = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError("invalid v2 fingerprint parent JSON") from error
    require(isinstance(parent, dict), "v2 fingerprint parent is not an object")
    require(parent.get("schema") == hardening.SCHEMA, "v2 parent schema drift")
    require(parent.get("status") == hardening.STATUS, "v2 parent status drift")
    require(
        parent.get("receipt_sha256") == artifact_self_hash(parent),
        "v2 parent self-digest drift",
    )
    require(
        raw == base.canonical_json_bytes(hardening.build_receipt()),
        "v2 parent fails byte-exact producer replay",
    )
    return raw, parent


def q5_sum(values: Iterable[Q5]) -> Q5:
    total = ZERO
    for value in values:
        total = q5_add(total, value)
    return total


def dot(x: tuple[Q5, Q5, Q5], y: tuple[Q5, Q5, Q5]) -> Q5:
    return q5_sum(q5_mul(x[axis], y[axis]) for axis in range(3))


def add_rays(
    vertices: list[tuple[Q5, Q5, Q5]], indices: Iterable[int]
) -> tuple[Q5, Q5, Q5]:
    selected = tuple(indices)
    return tuple(
        q5_sum(vertices[index][axis] for index in selected) for axis in range(3)
    )  # type: ignore[return-value]


def q5_lt(x: Q5, y: Q5) -> bool:
    return q5_sign(q5_sub(x, y)) < 0


def q5_gt(x: Q5, y: Q5) -> bool:
    return q5_sign(q5_sub(x, y)) > 0


def decompose_on_i6(
    polynomial: base.Poly,
    i6: base.Poly,
) -> tuple[Q5, Q5]:
    """Return ``constant, amplitude`` on the unit sphere.

    The caller supplies a carrier moment polynomial known to live in the
    constant plus ``I6`` line.  The residue check makes that premise an exact
    calculation rather than a label.
    """

    reduced = base.p_reduce_sphere(polynomial)
    i6_reduced = base.p_reduce_sphere(i6)
    probe = next(mono for mono in i6_reduced if sum(mono) > 0)
    amplitude = q5_div(reduced.get(probe, ZERO), i6_reduced[probe])
    residue = base.p_add(
        reduced,
        base.p_scale_q5(i6_reduced, q5_neg(amplitude)),
    )
    constant = residue.pop((0, 0, 0), ZERO)
    require(base.p_is_zero(residue), "moment leaves the constant-plus-I6 line")
    return constant, amplitude


def through_eighth_order_i6_template_certificate(
    *,
    t6_amplitude_override: Fraction | None = None,
    t8_amplitude_override: Fraction | None = None,
) -> dict[str, Any]:
    """Certify the exact anisotropic coefficient through eighth order."""

    vertices = base.cartesian_vertices()
    require(
        all(dot(vertex, vertex) == base.NORM_SQ for vertex in vertices),
        "carrier directions do not share the declared unit normalization",
    )
    cartesian = base.build_cartesian_frame()
    i6 = cartesian["_i6_poly_object"]
    moment6 = base.normalized_moment(vertices, 6, base.NORM_SQ)
    moment8 = base.normalized_moment(vertices, 8, base.NORM_SQ)
    constant6, amplitude6 = decompose_on_i6(moment6, i6)
    constant8, amplitude8 = decompose_on_i6(moment8, i6)

    require(constant6 == q5(Fraction(12, 7)), "sixth-moment constant drift")
    require(amplitude6 == q5(Fraction(64, 175)), "sixth-moment I6 drift")
    require(constant8 == q5(Fraction(4, 3)), "eighth-moment constant drift")
    require(amplitude8 == q5(Fraction(256, 375)), "eighth-moment I6 drift")
    require(amplitude6[1] == amplitude8[1] == 0, "moment amplitudes leave Q")

    if t6_amplitude_override is not None:
        require(
            t6_amplitude_override == amplitude6[0],
            "injected sixth-moment amplitude disagrees with exact geometry",
        )
    if t8_amplitude_override is not None:
        require(
            t8_amplitude_override == amplitude8[0],
            "injected eighth-moment amplitude disagrees with exact geometry",
        )

    leading = amplitude6[0] / factorial(6)
    eighth = amplitude8[0] / factorial(8)
    require(leading == Fraction(4, 7875), "cosine x^6 amplitude drift")
    require(eighth == Fraction(2, 118125), "cosine x^8 amplitude drift")
    require(leading / eighth == 30, "cosine amplitude ratio drift")

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


def cosine_tail_bounds(
    *,
    x_max: Fraction = Fraction(1),
    port_count: int = 12,
    gradient_geometric_factor_override: Fraction | None = None,
    hessian_geometric_factor_override: Fraction | None = None,
    normalization_lower_coefficient_override: Fraction | None = None,
    include_radial_hessian_correction: bool = True,
) -> dict[str, Any]:
    """Exact uniform C1/C2 bounds for the normalized Taylor tail.

    For the terms with powers ``2m >= 10``, differentiation costs one or
    two factorials.  On ``x <= 1`` the successive gradient terms have ratio
    at most ``1/110`` and the successive Hessian terms ratio at most ``1/90``.
    The intrinsic spherical Hessian is bounded by the Euclidean Hessian plus
    the radial-gradient correction.
    """

    require(port_count == 12, "tail proof requires the complete twelve-port orbit")
    require(x_max > 0, "tail range must be positive")
    require(x_max <= 1, "tail geometric-ratio proof is declared only for x <= 1")

    amplitude_lower_coefficient = Fraction(2) * (30 - x_max * x_max) / 118125
    require(amplitude_lower_coefficient > 0, "anisotropic normalization can vanish")
    if normalization_lower_coefficient_override is not None:
        require(
            normalization_lower_coefficient_override
            == amplitude_lower_coefficient,
            "injected normalization lower coefficient disagrees with A(x)",
        )

    gradient_geometric_factor = Fraction(110, 109)
    hessian_geometric_factor = Fraction(90, 89)
    if gradient_geometric_factor_override is not None:
        require(
            gradient_geometric_factor_override == gradient_geometric_factor,
            "injected gradient geometric factor disagrees with the ratio bound",
        )
    if hessian_geometric_factor_override is not None:
        require(
            hessian_geometric_factor_override == hessian_geometric_factor,
            "injected Hessian geometric factor disagrees with the ratio bound",
        )
    require(
        include_radial_hessian_correction,
        "intrinsic Hessian bound requires the radial-gradient correction",
    )

    gradient_raw_factor = (
        Fraction(port_count, factorial(9))
        * gradient_geometric_factor
        * x_max**10
    )
    euclidean_hessian_raw_factor = (
        Fraction(port_count, factorial(8))
        * hessian_geometric_factor
        * x_max**10
    )
    normalization_lower = amplitude_lower_coefficient * x_max**6
    gradient = gradient_raw_factor / normalization_lower
    intrinsic_hessian = (
        euclidean_hessian_raw_factor + gradient_raw_factor
    ) / normalization_lower

    if x_max == 1:
        require(gradient == Fraction(6875, 101152), "canonical C1 bound drift")
        require(
            intrinsic_hessian == Fraction(383125, 562658),
            "canonical C2 bound drift",
        )
    require(gradient < Fraction(7, 100), "C1 bound exceeds declared envelope")
    require(
        intrinsic_hessian < Fraction(7, 10),
        "C2 bound exceeds declared envelope",
    )

    return {
        "x_max": str(x_max),
        "range": f"0 < x <= {x_max}",
        "tail_start": "x^10",
        "normalization_lower_bound": (
            f"A(x) >= {amplitude_lower_coefficient} x^6 for "
            f"0 < x <= {x_max}"
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
            "port_count": port_count,
            "gradient_successive_ratio_upper": "1/110",
            "gradient_geometric_sum_factor": str(gradient_geometric_factor),
            "euclidean_hessian_successive_ratio_upper": "1/90",
            "euclidean_hessian_geometric_sum_factor": str(
                hessian_geometric_factor
            ),
            "normalization_lower_coefficient": str(
                amplitude_lower_coefficient
            ),
            "radial_hessian_correction_included": True,
            "radial_hessian_correction_bound": "the normalized gradient tail",
        },
        "clean_envelopes": {
            "gradient_strictly_below": "7/100",
            "intrinsic_hessian_strictly_below": "7/10",
        },
        "arithmetic": "exact Fraction arithmetic; no transcendental evaluation",
    }


def critical_axis_separation(
    vertices: list[tuple[Q5, Q5, Q5]] | None = None,
    *,
    neighborhood_sine_squared: Q5 | None = None,
    promote_separation_to_global_persistence: bool = False,
) -> dict[str, Any]:
    """Compute the 31 unoriented vertex, face, and edge axes exactly."""

    if vertices is None:
        vertices = base.cartesian_vertices()
    faces, edges = hardening.derive_adjacency_faces_edges(vertices)
    typed_rays = [
        ("vertex", vertex) for vertex in vertices
    ] + [
        ("edge", add_rays(vertices, edge)) for edge in edges
    ] + [
        ("face", add_rays(vertices, face)) for face in faces
    ]
    require(len(typed_rays) == 62, "critical ray census drift")

    axes: list[tuple[str, tuple[Q5, Q5, Q5]]] = []
    for orbit, ray in typed_rays:
        parallel = [
            entry
            for entry in axes
            if hardening.cross(ray, entry[1]) == (ZERO, ZERO, ZERO)
        ]
        if parallel:
            require(
                len(parallel) == 1 and dot(ray, parallel[0][1]) != ZERO,
                "axis deduplication is ambiguous",
            )
            require(
                parallel[0][0] == orbit,
                "different critical orbit types share an axis",
            )
            continue
        axes.append((orbit, ray))

    orbit_counts = {
        orbit: sum(1 for kind, _ in axes if kind == orbit)
        for orbit in ("vertex", "face", "edge")
    }
    require(
        orbit_counts == {"vertex": 6, "face": 10, "edge": 15},
        "unoriented critical-axis census drift",
    )
    require(len(axes) == 31, "critical axis count drift")

    max_cosine_squared: Q5 | None = None
    maximizing_pair: tuple[int, int] | None = None
    for left, (_, x) in enumerate(axes):
        for right in range(left + 1, len(axes)):
            y = axes[right][1]
            cosine_squared = q5_div(
                q5_pow(dot(x, y), 2),
                q5_mul(dot(x, x), dot(y, y)),
            )
            require(
                q5_sign(cosine_squared) >= 0
                and q5_sign(q5_sub(ONE, cosine_squared)) > 0,
                "two distinct critical axes coincide",
            )
            if max_cosine_squared is None or q5_gt(
                cosine_squared, max_cosine_squared
            ):
                max_cosine_squared = cosine_squared
                maximizing_pair = (left, right)

    require(max_cosine_squared is not None, "axis separation pair missing")
    expected_cosine_squared = q5(Fraction(1, 2), Fraction(1, 6))
    expected_sine_squared = q5(Fraction(1, 2), Fraction(-1, 6))
    require(
        max_cosine_squared == expected_cosine_squared,
        "minimum critical-axis separation drift",
    )
    minimum_sine_squared = q5_sub(ONE, max_cosine_squared)
    require(
        minimum_sine_squared == expected_sine_squared,
        "minimum sine-squared separation drift",
    )

    if neighborhood_sine_squared is None:
        neighborhood_sine_squared = q5(Fraction(1, 64))
    require(
        q5_gt(
            minimum_sine_squared,
            q5_scale(neighborhood_sine_squared, 4),
        ),
        "declared critical-axis neighborhoods are not separated",
    )
    require(
        not promote_separation_to_global_persistence,
        "axis separation alone cannot promote the global persistence theorem",
    )
    require(maximizing_pair is not None, "maximizing axis pair missing")
    maximizing_pair_orbits = [
        axes[maximizing_pair[0]][0],
        axes[maximizing_pair[1]][0],
    ]

    return {
        "unoriented_axis_count": 31,
        "oriented_critical_point_count": 62,
        "axis_counts_by_orbit": orbit_counts,
        "maximum_squared_axis_cosine": q5_str(max_cosine_squared),
        "minimum_squared_axis_sine": "1/2-sqrt5/6",
        "maximizing_pair_indices": list(maximizing_pair or ()),
        "maximizing_pair_orbit_types": maximizing_pair_orbits,
        "declared_local_neighborhood": "sin^2(angle to an axis) <= "
        + (
            str(neighborhood_sine_squared[0])
            if neighborhood_sine_squared[1] == 0
            else q5_str(neighborhood_sine_squared)
        ),
        "separation_check": (
            "four times the neighborhood sine-squared radius is strictly "
            "below the exact minimum axis sine-squared separation"
        ),
        "promotion_from_separation_permitted": False,
    }


def fail_closed_controls() -> dict[str, Any]:
    controls = []

    try:
        cosine_tail_bounds(x_max=Fraction(2))
    except base.FingerprintError:
        controls.append(
            {
                "control": "extend the Taylor-tail proof to x = 2",
                "expected_failure": "declared geometric-ratio range",
                "detector_fired": True,
            }
        )
    else:
        raise base.FingerprintError("x-range mutation escaped the tail detector")

    try:
        cosine_tail_bounds(port_count=13)
    except base.FingerprintError:
        controls.append(
            {
                "control": "inject a thirteenth port into the tail bound",
                "expected_failure": "complete twelve-port normalization",
                "detector_fired": True,
            }
        )
    else:
        raise base.FingerprintError("port-count mutation escaped the tail detector")

    try:
        through_eighth_order_i6_template_certificate(
            t6_amplitude_override=Fraction(65, 175)
        )
    except base.FingerprintError:
        controls.append(
            {
                "control": "mutate the sixth-moment I6 amplitude",
                "expected_failure": "exact carrier moment decomposition",
                "detector_fired": True,
            }
        )
    else:
        raise base.FingerprintError("sixth-moment mutation escaped the template detector")

    try:
        through_eighth_order_i6_template_certificate(
            t8_amplitude_override=Fraction(257, 375)
        )
    except base.FingerprintError:
        controls.append(
            {
                "control": "mutate the eighth-moment I6 amplitude",
                "expected_failure": "exact carrier moment decomposition",
                "detector_fired": True,
            }
        )
    else:
        raise base.FingerprintError("eighth-moment mutation escaped the template detector")

    tail_mutations = (
        (
            "mutate the gradient-tail geometric sum factor",
            {"gradient_geometric_factor_override": Fraction(111, 109)},
            "exact gradient-tail ratio bound",
        ),
        (
            "mutate the Hessian-tail geometric sum factor",
            {"hessian_geometric_factor_override": Fraction(91, 89)},
            "exact Hessian-tail ratio bound",
        ),
        (
            "mutate the anisotropic normalization lower coefficient",
            {"normalization_lower_coefficient_override": Fraction(59, 118125)},
            "exact positive normalization",
        ),
        (
            "drop the spherical radial-Hessian correction",
            {"include_radial_hessian_correction": False},
            "intrinsic Hessian conversion",
        ),
    )
    for label, kwargs, expected_failure in tail_mutations:
        try:
            cosine_tail_bounds(**kwargs)
        except base.FingerprintError:
            controls.append(
                {
                    "control": label,
                    "expected_failure": expected_failure,
                    "detector_fired": True,
                }
            )
        else:
            raise base.FingerprintError(f"{label} escaped the tail detector")

    vertices = base.cartesian_vertices()
    mutated = list(vertices)
    antipode = hardening.derive_antipode(vertices)
    first, opposite = 0, antipode[0]

    def rotate_x(
        vector: tuple[Q5, Q5, Q5],
    ) -> tuple[Q5, Q5, Q5]:
        x, y, z = vector
        return (
            x,
            q5_add(q5_scale(y, Fraction(3, 5)), q5_scale(z, Fraction(-4, 5))),
            q5_add(q5_scale(y, Fraction(4, 5)), q5_scale(z, Fraction(3, 5))),
        )

    mutated[first] = rotate_x(vertices[first])
    mutated[opposite] = rotate_x(vertices[opposite])
    try:
        critical_axis_separation(mutated)
    except base.FingerprintError:
        controls.append(
            {
                "control": "rotate one antipodal pair on the common sphere",
                "expected_failure": "exact icosahedral axis separation",
                "detector_fired": True,
            }
        )
    else:
        raise base.FingerprintError("geometry mutation escaped the axis detector")

    try:
        critical_axis_separation(
            neighborhood_sine_squared=q5(Fraction(1, 4))
        )
    except base.FingerprintError:
        controls.append(
            {
                "control": "enlarge the local axis neighborhood beyond separation",
                "expected_failure": "exact disjoint-neighborhood inequality",
                "detector_fired": True,
            }
        )
    else:
        raise base.FingerprintError("neighborhood-radius mutation escaped the axis detector")

    try:
        critical_axis_separation(
            promote_separation_to_global_persistence=True
        )
    except base.FingerprintError:
        controls.append(
            {
                "control": "promote axis separation to global persistence",
                "expected_failure": "global cover and local uniqueness remain open",
                "detector_fired": True,
            }
        )
    else:
        raise base.FingerprintError("separation-promotion mutation escaped the detector")

    require(len(controls) == 11, "persistence control count drift")
    require(all(row["detector_fired"] for row in controls), "a control failed")
    return {"controls": controls, "all_detectors_fired": True}


def build_receipt(*, parent_path: Path = PARENT_PATH) -> dict[str, Any]:
    parent_raw, parent = load_v2_parent(parent_path)
    receipt = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 654,
        "extends": {
            "schema": parent["schema"],
            "status": parent["status"],
            "path": "code/a5_fingerprint/runtime/a5_multipole_fixed_point_receipt_v2.json",
            "reason": (
                "append-only quantitative extension; v1 and v2 remain immutable"
            ),
        },
        "parent_pin": {
            "path": "code/a5_fingerprint/runtime/a5_multipole_fixed_point_receipt_v2.json",
            "bytes": len(parent_raw),
            "sha256": base.tagged_sha256(parent_raw),
            "receipt_sha256": parent["receipt_sha256"],
        },
        "through_eighth_order_i6_template": (
            through_eighth_order_i6_template_certificate()
        ),
        "normalized_tail_bounds": cosine_tail_bounds(),
        "critical_axis_separation": critical_axis_separation(),
        "quantitative_persistence": {
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
        },
        "fail_closed_controls": fail_closed_controls(),
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
            "physical_map_open": True,
        },
    }
    receipt["receipt_sha256"] = artifact_self_hash(receipt)
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    receipt = build_receipt()
    if args.write:
        RUNTIME.mkdir(exist_ok=True)
        RECEIPT_PATH.write_bytes(base.canonical_json_bytes(receipt))
        print(RECEIPT_PATH)
        return 0
    print(json.dumps(receipt["quantitative_persistence"], indent=2))
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
