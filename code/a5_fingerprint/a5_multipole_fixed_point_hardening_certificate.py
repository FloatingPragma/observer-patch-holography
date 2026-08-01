#!/usr/bin/env python3
"""Append-only hardening packet for the issue-654 A5 fingerprint.

The original v1 fingerprint receipt is an immutable parent of the custodied
FZ-11 prediction.  This producer therefore does not rewrite it.  It emits a
v2 receipt which:

* states explicitly that the quantitative Morse-persistence radius remains
  open;
* serializes the ordered exact carrier geometry, antipodes, positive face
  orientation, and equal-trace declaration;
* pins and validates the source-bound issue-599 response artifact instead of
  presenting the representative ``R = -J`` as a consequence of incidence
  alone;
* verifies every one of the 62 stationary directions directly against the original
  Cartesian tangent-gradient equation, after the meridian calculation has
  squared an equation;
* replaces two weak negative controls with an exact on-sphere angular
  perturbation and an injected response-sign mutation evaluated through the
  actual kernel matrices; and
* types M1--M3 as local nonphysical theorem rows rather than pretending that
  they are independently registered nature-facing predictions.

No interval gradient bound or finite persistence range is asserted here.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import a5_multipole_fixed_point_certificate as base


HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
LEGACY_RECEIPT_PATH = RUNTIME / "a5_multipole_fixed_point_receipt.json"
RECEIPT_PATH = RUNTIME / "a5_multipole_fixed_point_receipt_v2.json"
SEMANTIC_RESPONSE_PATH = (
    HERE.parent / "a5_closure" / "manifests" /
    "charged_response_semantic_artifact.json"
)

SCHEMA = "oph.a5_multipole_fixed_point_receipt.v2"
STATUS = (
    "EXACT_A5_FINGERPRINT_CORE__QUANTITATIVE_PERSISTENCE_OPEN__"
    "PHYSICAL_MAP_OPEN"
)
SEMANTIC_RESPONSE_SCHEMA = "oph.charged_response_semantic_artifact.v3"

# ``base.cartesian_vertices`` is ordered by the two signs first and the cyclic
# coordinate placement second.  The semantic response artifact uses p00..p11.
PORT_IDS = (
    "p00", "p04", "p08", "p01", "p05", "p10",
    "p02", "p06", "p09", "p03", "p07", "p11",
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
    body = {key: item for key, item in value.items() if key != "artifact_sha256"}
    encoded = json.dumps(
        body,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


def _load_json(path: Path) -> tuple[bytes, dict[str, Any]]:
    raw = path.read_bytes()
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError(f"invalid JSON parent: {path}") from error
    require(isinstance(value, dict), f"JSON parent is not an object: {path}")
    return raw, value


def load_legacy_parent() -> tuple[bytes, dict[str, Any]]:
    raw, value = _load_json(LEGACY_RECEIPT_PATH)
    require(value.get("schema") == base.SCHEMA, "legacy fingerprint schema drift")
    require(value.get("status") == base.STATUS, "legacy fingerprint status drift")
    require(
        raw == base.canonical_json_bytes(base.build_receipt()),
        "legacy fingerprint fails byte-exact replay",
    )
    claimed = value.get("receipt_sha256")
    body = {key: item for key, item in value.items() if key != "receipt_sha256"}
    require(
        claimed == base.tagged_sha256(base.canonical_json_bytes(body)),
        "legacy fingerprint self-digest drift",
    )
    return raw, value


def parse_artifact_q5(value: str) -> Q5:
    text = value.strip()
    if "sqrt(5)" not in text:
        return q5(Fraction(text))
    parts = text.split(" + ")
    require(len(parts) == 2, f"malformed artifact Q(sqrt5) scalar: {value!r}")
    radical = parts[1]
    require(
        radical.endswith("*sqrt(5)"),
        f"malformed artifact radical: {value!r}",
    )
    return q5(Fraction(parts[0]), Fraction(radical[: -len("*sqrt(5)")]))


def load_semantic_response(
    path: Path = SEMANTIC_RESPONSE_PATH,
) -> tuple[bytes, dict[str, Any]]:
    raw, artifact = _load_json(path)
    require(
        artifact.get("schema") == SEMANTIC_RESPONSE_SCHEMA,
        "semantic response schema drift",
    )
    require(artifact.get("issue") == 599, "semantic response issue drift")
    require(
        artifact.get("artifact_sha256") == artifact_self_hash(artifact),
        "semantic response self-digest drift",
    )
    response = artifact.get("source_response")
    require(isinstance(response, dict), "semantic source_response missing")
    require(
        response.get("operator") == "negative_graph_antipode_involution",
        "semantic response fails to select the representative R=-J",
    )
    require(
        response.get("impulse_readback_response_executed") is True
        and response.get("physical_perturb_readback_source_bound") is True,
        "semantic response is not source-bound and executed",
    )
    require(
        response.get("sector_eigenvalues")
        == {
            "frame_band": 1,
            "kernel_band": 1,
            "quintet_band": -1,
            "unit_band": -1,
        },
        "semantic response sector signs drift",
    )
    return raw, artifact


def q5_sum(values) -> Q5:
    total = ZERO
    for value in values:
        total = q5_add(total, value)
    return total


def dot(x: tuple[Q5, Q5, Q5], y: tuple[Q5, Q5, Q5]) -> Q5:
    return q5_sum(q5_mul(x[axis], y[axis]) for axis in range(3))


def add_vectors(*vectors: tuple[Q5, Q5, Q5]) -> tuple[Q5, Q5, Q5]:
    return tuple(q5_sum(vector[axis] for vector in vectors) for axis in range(3))  # type: ignore[return-value]


def cross(x: tuple[Q5, Q5, Q5], y: tuple[Q5, Q5, Q5]) -> tuple[Q5, Q5, Q5]:
    return (
        q5_sub(q5_mul(x[1], y[2]), q5_mul(x[2], y[1])),
        q5_sub(q5_mul(x[2], y[0]), q5_mul(x[0], y[2])),
        q5_sub(q5_mul(x[0], y[1]), q5_mul(x[1], y[0])),
    )


def determinant(
    x: tuple[Q5, Q5, Q5],
    y: tuple[Q5, Q5, Q5],
    z: tuple[Q5, Q5, Q5],
) -> Q5:
    return dot(x, cross(y, z))


def derive_adjacency_faces_edges(verts):
    inv_norm = q5_div(ONE, base.NORM_SQ)

    def unit_dot(i: int, j: int) -> Q5:
        return q5_mul(dot(verts[i], verts[j]), inv_norm)

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
    require(len(edges) == 30, "edge census drift")
    require(len(faces) == 20, "face census drift")
    return faces, edges


def derive_antipode(verts) -> list[int]:
    out = []
    for i, vertex in enumerate(verts):
        partners = [
            j for j, other in enumerate(verts)
            if dot(vertex, other) == q5_neg(base.NORM_SQ)
        ]
        require(len(partners) == 1, "antipode is not unique")
        out.append(partners[0])
    require(all(out[out[i]] == i for i in range(12)), "antipode is not involutive")
    return out


def serialized_geometry(
    verts,
    semantic_artifact: dict[str, Any],
) -> dict[str, Any]:
    require(len(PORT_IDS) == len(verts) == 12, "ordered port count drift")
    port_frame = semantic_artifact.get("port_vertex_frame")
    require(isinstance(port_frame, dict), "semantic port frame missing")
    for port_id, vertex in zip(PORT_IDS, verts, strict=True):
        supplied = port_frame.get(port_id)
        require(
            isinstance(supplied, list) and len(supplied) == 3,
            f"semantic coordinate missing for {port_id}",
        )
        require(
            tuple(parse_artifact_q5(str(item)) for item in supplied) == vertex,
            f"semantic coordinate drift at {port_id}",
        )

    antipode = derive_antipode(verts)
    carrier_binding = semantic_artifact.get("carrier_binding")
    require(isinstance(carrier_binding, dict), "semantic carrier binding missing")
    require(
        carrier_binding.get("port_order") == [f"p{index:02d}" for index in range(12)],
        "semantic canonical port order drift",
    )
    supplied_antipode = carrier_binding.get("antipode")
    require(isinstance(supplied_antipode, dict), "semantic antipode map missing")
    derived_by_port = {
        PORT_IDS[index]: PORT_IDS[partner]
        for index, partner in enumerate(antipode)
    }
    require(derived_by_port == supplied_antipode, "semantic antipode map drift")
    response = semantic_artifact["source_response"]
    canonical_ports = carrier_binding["port_order"]
    response_antipode = response.get("antipode_port_map")
    require(
        isinstance(response_antipode, list) and len(response_antipode) == 12,
        "semantic response antipode index map missing",
    )
    require(
        {
            port: canonical_ports[int(response_antipode[index])]
            for index, port in enumerate(canonical_ports)
        }
        == derived_by_port,
        "semantic response antipode index map disagrees with the exact geometry",
    )

    faces, _ = derive_adjacency_faces_edges(verts)
    oriented_faces = []
    for i, j, k in faces:
        if q5_sign(determinant(verts[i], verts[j], verts[k])) < 0:
            j, k = k, j
        det = determinant(verts[i], verts[j], verts[k])
        require(q5_sign(det) > 0, "face orientation is degenerate")
        oriented_faces.append(
            {
                "indices": [i, j, k],
                "ports": [PORT_IDS[i], PORT_IDS[j], PORT_IDS[k]],
                "determinant": q5_str(det),
            }
        )
    require(len(oriented_faces) == 20, "oriented face census drift")

    orientation = semantic_artifact.get("orientation_convention")
    require(isinstance(orientation, dict), "semantic orientation convention missing")
    require(
        orientation.get("faces")
        == "every oriented face has positive determinant in the assigned frame",
        "semantic face orientation convention drift",
    )

    rows = [
        {
            "index": index,
            "port": port_id,
            "coordinates_qsqrt5": [q5_str(value) for value in vertex],
            "norm_squared": q5_str(dot(vertex, vertex)),
            "antipode_index": antipode[index],
            "antipode_port": PORT_IDS[antipode[index]],
            "trace_weight": "1",
        }
        for index, (port_id, vertex) in enumerate(zip(PORT_IDS, verts, strict=True))
    ]
    require(
        all(dot(vertex, vertex) == base.NORM_SQ for vertex in verts),
        "serialized vertex norm drift",
    )
    geometry_body = {
        "ordered_ports": rows,
        "oriented_faces": oriented_faces,
        "equal_trace_declaration": {
            "declared": True,
            "weights": ["1"] * 12,
            "scope": "equal scalar activation of the twelve finite ports",
        },
        "orientation_source": orientation,
    }
    return {
        **geometry_body,
        "geometry_sha256": base.tagged_sha256(base.canonical_json_bytes(geometry_body)),
    }


def homogeneous_i6(verts) -> base.Poly:
    """The degree-six homogeneous harmonic whose unit restriction is I6."""

    r2 = base.radial_power(1)
    r4 = base.radial_power(2)
    r6 = base.radial_power(3)
    m2 = base.moment_sum(verts, 2)
    m4 = base.moment_sum(verts, 4)
    m6 = base.moment_sum(verts, 6)
    n = base.NORM_SQ
    total = base.p_zero()
    total = base.p_add(
        total,
        base.p_scale_q5(m6, q5_scale(q5_div(ONE, q5_pow(n, 3)), 231)),
    )
    total = base.p_add(
        total,
        base.p_scale_q5(
            base.p_mul(m4, r2),
            q5_scale(q5_div(ONE, q5_pow(n, 2)), -315),
        ),
    )
    total = base.p_add(
        total,
        base.p_scale_q5(
            base.p_mul(m2, r4),
            q5_scale(q5_div(ONE, n), 105),
        ),
    )
    total = base.p_add(total, base.p_scale(r6, -5 * len(verts)))
    return base.p_scale(total, Fraction(25, 132 * 16))


def poly_derivative(poly: base.Poly, axis: int) -> base.Poly:
    out: base.Poly = {}
    for mono, coeff in poly.items():
        power = mono[axis]
        if power == 0:
            continue
        derived = list(mono)
        derived[axis] -= 1
        out[tuple(derived)] = q5_scale(coeff, power)
    return out


def tangent_gradient_cross(
    poly: base.Poly,
    ray: tuple[Q5, Q5, Q5],
) -> tuple[Q5, Q5, Q5]:
    gradient = tuple(
        base.p_eval(poly_derivative(poly, axis), ray) for axis in range(3)
    )
    return cross(ray, gradient)  # type: ignore[arg-type]


def squared_meridian_sign_enumeration() -> dict[str, Any]:
    """Replay every sign branch introduced by the meridian squaring.

    For a nonzero, nonpolar latitude ``t = c^2``, ``sin(5 phi) = 0``
    leaves ``q = cos(5 phi)`` equal to ``+1`` or ``-1``.  Each choice of
    ``q`` represents five azimuths.  Squaring the remaining equation forgets
    the relative sign between ``c`` and ``q``.  Consequently each algebraic
    ``t`` root supplies twenty squared candidates: two signs of ``c``, two
    signs of ``q``, and five azimuths.  The original equation accepts exactly
    one q sign for each c sign.
    """

    latitudes = [
        ("vertex_ring", q5(Fraction(1, 5))),
        ("face_high", q5(Fraction(1, 3), Fraction(2, 15))),
        ("face_low", q5(Fraction(1, 3), Fraction(-2, 15))),
        ("edge_high", q5(Fraction(1, 2), Fraction(1, 10))),
        ("edge_low", q5(Fraction(1, 2), Fraction(-1, 10))),
    ]
    rows = []
    total_candidates = 0
    total_accepted = 0
    total_rejected = 0
    for name, t_value in latitudes:
        # P6'(c) = c * g(c^2).
        g_value = q5_scale(
            q5_add(
                q5_add(q5_scale(q5_pow(t_value, 2), 1386), q5_scale(t_value, -1260)),
                q5(210),
            ),
            Fraction(1, 16),
        )
        six_t_minus_one = q5_sub(q5_scale(t_value, 6), ONE)
        require(q5_sign(g_value) != 0, f"P6'/c vanishes at {name}")
        require(q5_sign(six_t_minus_one) != 0, f"6t-1 vanishes at {name}")

        # Exact magnitude replay of the squared meridian equation:
        # c^2 g(t)^2 = (21/8)^2 (1-t)^3 (6t-1)^2.
        lhs_squared = q5_mul(t_value, q5_pow(g_value, 2))
        rhs_squared = q5_scale(
            q5_mul(
                q5_pow(q5_sub(ONE, t_value), 3),
                q5_pow(six_t_minus_one, 2),
            ),
            Fraction(441, 64),
        )
        require(lhs_squared == rhs_squared, f"squared magnitude drift at {name}")

        branches = []
        accepted = 0
        rejected = 0
        for c_sign in (-1, 1):
            required_q_sign = (
                c_sign * q5_sign(g_value) * q5_sign(six_t_minus_one)
            )
            for q_sign in (-1, 1):
                passes_original = q_sign == required_q_sign
                directions = 5
                if passes_original:
                    accepted += directions
                else:
                    rejected += directions
                branches.append(
                    {
                        "c_sign": c_sign,
                        "cos_5phi_sign": q_sign,
                        "azimuth_count": directions,
                        "passes_original_unsquared_equation": passes_original,
                    }
                )
        require(accepted == 10 and rejected == 10, f"sign census drift at {name}")
        rows.append(
            {
                "orbit_latitude": name,
                "c_squared": q5_str(t_value),
                "P6_prime_over_c": q5_str(g_value),
                "six_c_squared_minus_one": q5_str(six_t_minus_one),
                "branches": branches,
                "squared_candidate_directions": accepted + rejected,
                "accepted_by_original_equation": accepted,
                "rejected_as_squaring_extraneous": rejected,
            }
        )
        total_candidates += accepted + rejected
        total_accepted += accepted
        total_rejected += rejected

    require(total_candidates == 100, "squared meridian candidate total drift")
    require(total_accepted == 50, "unsquared meridian acceptance total drift")
    require(total_rejected == 50, "squaring-extraneous total drift")
    return {
        "candidate_convention": (
            "five nonzero c^2 roots times two c signs times two cos(5phi) "
            "signs times five azimuths"
        ),
        "latitudes": rows,
        "squared_candidate_directions": total_candidates,
        "accepted_by_original_equation": total_accepted,
        "rejected_as_squaring_extraneous": total_rejected,
    }


def exact_stationary_ray_replay(
    verts,
    hom_i6: base.Poly | None = None,
) -> dict[str, Any]:
    if hom_i6 is None:
        hom_i6 = homogeneous_i6(verts)
    cartesian = base.build_cartesian_frame()
    sphere_i6 = cartesian["_i6_poly_object"]
    require(
        base.p_is_zero(
            base.p_reduce_sphere(
                base.p_add(hom_i6, base.p_scale(sphere_i6, -1))
            )
        ),
        "homogeneous I6 does not restrict to the registered sphere polynomial",
    )

    faces, edges = derive_adjacency_faces_edges(verts)
    rays: list[tuple[str, tuple[Q5, Q5, Q5]]] = []
    rays.extend(("vertex", vertex) for vertex in verts)
    rays.extend(
        ("face", add_vectors(verts[i], verts[j], verts[k]))
        for i, j, k in faces
    )
    rays.extend(
        ("edge", add_vectors(verts[i], verts[j]))
        for i, j in edges
    )
    require(len(rays) == 62, "stationary ray input count drift")

    for name, ray in rays:
        require(
            tangent_gradient_cross(hom_i6, ray) == (ZERO, ZERO, ZERO),
            f"{name} ray fails the original tangent-gradient equation",
        )

    # No same-direction duplicates are allowed.  Antipodal directions are
    # distinct points of S2 and are counted separately.
    antipodal_pairs = 0
    for i, (_, left) in enumerate(rays):
        for _, right in rays[i + 1 :]:
            if cross(left, right) != (ZERO, ZERO, ZERO):
                continue
            alignment = q5_sign(dot(left, right))
            require(alignment < 0, "stationary ray list contains a duplicate direction")
            antipodal_pairs += 1
    require(antipodal_pairs == 31, "stationary antipodal-pair count drift")

    pole = verts[0]
    pole_norm = dot(pole, pole)
    latitude_counts: dict[str, int] = {}
    for _, ray in rays:
        c_squared = q5_div(q5_pow(dot(ray, pole), 2), q5_mul(dot(ray, ray), pole_norm))
        key = q5_str(c_squared)
        latitude_counts[key] = latitude_counts.get(key, 0) + 1

    expected = {
        q5_str(q5(1)): 2,
        q5_str(q5(Fraction(1, 5))): 10,
        q5_str(q5(Fraction(1, 3), Fraction(2, 15))): 10,
        q5_str(q5(Fraction(1, 3), Fraction(-2, 15))): 10,
        q5_str(q5(Fraction(1, 2), Fraction(1, 10))): 10,
        q5_str(q5(Fraction(1, 2), Fraction(-1, 10))): 10,
        q5_str(q5(0)): 10,
    }
    require(latitude_counts == expected, "original-equation latitude census drift")
    require(
        all(q5_sign(value) > 0 and q5_sign(q5_sub(ONE, value)) > 0 for value in (
            q5(Fraction(1, 5)),
            q5(Fraction(1, 3), Fraction(2, 15)),
            q5(Fraction(1, 3), Fraction(-2, 15)),
            q5(Fraction(1, 2), Fraction(1, 10)),
            q5(Fraction(1, 2), Fraction(-1, 10)),
        )),
        "a squared-factor latitude left the physical interval",
    )

    sign_enumeration = squared_meridian_sign_enumeration()
    require(
        sign_enumeration["accepted_by_original_equation"] == 50,
        "meridian acceptance does not match the five ten-ray latitude orbits",
    )

    return {
        "equation": "x cross grad(H6)(x) = 0 for the homogeneous degree-six I6 polynomial",
        "coordinate_role": (
            "coordinate-invariant replay of the original tangent-gradient equation; "
            "it does not reuse the squared meridian equation"
        ),
        "accepted_directions": {
            "vertex": 12,
            "face": 20,
            "edge": 30,
            "total": 62,
        },
        "antipodal_pairs": antipodal_pairs,
        "accepted_latitude_counts": latitude_counts,
        "squared_meridian_sign_enumeration": sign_enumeration,
        "exclusion_argument": (
            "the squared meridian equation creates 100 sign-and-azimuth "
            "candidates over its five c^2 roots; exact replay of the original "
            "sign equation accepts 50 and rejects 50.  The accepted meridian "
            "directions, ten equatorial roots, and two poles are exactly the 62 "
            "distinct Cartesian tangent-gradient zeros"
        ),
    }


def response_signs_for_multiplier(verts, multiplier: int) -> list[int]:
    require(multiplier in (-1, 1), "response multiplier must be +-1")
    antipode = derive_antipode(verts)
    inv_norm = q5_div(ONE, base.NORM_SQ)
    dots = [
        [q5_mul(dot(verts[i], verts[j]), inv_norm) for j in range(12)]
        for i in range(12)
    ]
    kernels = [
        [
            [base.legendre_eval_q5(level, dots[i][j]) for j in range(12)]
            for i in range(12)
        ]
        for level in range(4)
    ]
    response = [
        [q5(multiplier) if antipode[i] == j else ZERO for j in range(12)]
        for i in range(12)
    ]

    def mat_mul(left, right):
        return [
            [
                q5_sum(q5_mul(left[i][k], right[k][j]) for k in range(12))
                for j in range(12)
            ]
            for i in range(12)
        ]

    signs = []
    for kernel in kernels:
        product = mat_mul(response, kernel)
        matched = []
        for sign in (-1, 1):
            if all(
                product[i][j] == q5_scale(kernel[i][j], sign)
                for i in range(12)
                for j in range(12)
            ):
                matched.append(sign)
        require(len(matched) == 1, "response mutation has no unique band sign")
        signs.append(matched[0])
    return signs


def hardened_controls(verts, legacy_controls: dict[str, Any]) -> dict[str, Any]:
    retained = [
        copy.deepcopy(row)
        for row in legacy_controls["controls"]
        if row["control"]
        not in {
            "perturbed vertex geometry (one pair stretched)",
            "R = +J response",
        }
    ]
    require(len(retained) == 6, "legacy control filtering drift")

    # Apply an exact rational orthogonal rotation to one antipodal pair.  All
    # twelve points remain on the same sphere and the pair remains antipodal.
    antipode = derive_antipode(verts)
    first, opposite = 0, antipode[0]

    def rotate_x(vector):
        x, y, z = vector
        return (
            x,
            q5_add(q5_scale(y, Fraction(3, 5)), q5_scale(z, Fraction(-4, 5))),
            q5_add(q5_scale(y, Fraction(4, 5)), q5_scale(z, Fraction(3, 5))),
        )

    perturbed = list(verts)
    perturbed[first] = rotate_x(verts[first])
    perturbed[opposite] = rotate_x(verts[opposite])
    require(
        perturbed[opposite]
        == tuple(q5_neg(value) for value in perturbed[first]),
        "on-sphere perturbation lost antipodality",
    )
    require(
        all(dot(vertex, vertex) == base.NORM_SQ for vertex in perturbed),
        "on-sphere perturbation changed a radius",
    )
    second = base.normalized_moment(perturbed, 2, base.NORM_SQ)
    isotropy_broken = second != base.p_scale(base.radial_power(1), 4)
    require(isotropy_broken, "on-sphere angular perturbation escaped detector")
    retained.append(
        {
            "control": "one antipodal pair rotated on the exact common sphere",
            "expected_failure": "second-moment isotropy",
            "common_radius_preserved": True,
            "antipode_preserved": True,
            "detector_fired": True,
        }
    )

    production_signs = response_signs_for_multiplier(verts, -1)
    mutated_signs = response_signs_for_multiplier(verts, 1)
    require(production_signs == [-1, 1, -1, 1], "production R=-J signs drift")
    require(mutated_signs == [1, -1, 1, -1], "mutated R=+J signs drift")
    retained.append(
        {
            "control": "response multiplier mutation from -J to +J",
            "expected_failure": "source-bound band sign vector",
            "production_signs_levels_0_to_3": production_signs,
            "mutated_signs_levels_0_to_3": mutated_signs,
            "implementation": "both signs recomputed by matrix multiplication against all four Legendre Gram kernels",
            "detector_fired": mutated_signs != production_signs,
        }
    )
    require(len(retained) == 8, "hardened control count drift")
    require(all(row["detector_fired"] for row in retained), "hardened control failed")
    return {"controls": retained, "all_detectors_fired": True}


def local_rows(legacy_rules: dict[str, Any]) -> dict[str, Any]:
    rows = copy.deepcopy(legacy_rules["frozen_rows"])
    for row in rows.values():
        row["type"] = "local_nonphysical_conditional_theorem_row"
        row["registration_scope"] = (
            "this issue-654 certificate only; not a separately custodied or "
            "nature-facing frozen prediction"
        )
    rows["OPH-A5-M2"]["premise_ancestry"] = (
        "pure leading-band regime; exact Morse nondegeneracy gives qualitative "
        "small-perturbation persistence, while the explicit norm, gradient lower "
        "bound, and finite parameter range remain open"
    )
    rows["OPH-A5-M3"]["premise_ancestry"] = (
        "the carrier incidence derives the antipode J and its parity action; "
        "the separately pinned source-bound issue-599 impulse/readback artifact "
        "selects the conventional representative R=-J and its four signs; the "
        "physical response channel remains open"
    )
    return {
        "local_certificate_rows": rows,
        "external_physical_rows": copy.deepcopy(legacy_rules["unfrozen_rows"]),
        "exact_local_decision_rules": copy.deepcopy(
            legacy_rules["blind_decision_rules"]
        ),
        "comparison_boundary": copy.deepcopy(legacy_rules["comparison_boundary"]),
        "typing_note": (
            "empirical tolerances are not frozen by this ineligible mathematical "
            "packet; any nature-facing campaign requires its own pre-exposure "
            "contract"
        ),
    }


def build_receipt(
    *,
    semantic_response_path: Path = SEMANTIC_RESPONSE_PATH,
) -> dict[str, Any]:
    legacy_raw, legacy = load_legacy_parent()
    semantic_raw, semantic = load_semantic_response(semantic_response_path)
    verts = base.cartesian_vertices()

    cartesian = copy.deepcopy(legacy["cartesian_frame"])
    cartesian["serialized_geometry"] = serialized_geometry(verts, semantic)

    critical = copy.deepcopy(legacy["critical_points"])
    critical["original_equation_replay"] = exact_stationary_ray_replay(verts)
    critical["quantitative_persistence"] = {
        "status": "OPEN",
        "missing": [
            "declared perturbation family or normalized C2 ball",
            "interval lower bound for the tangent-gradient norm off the 62 neighborhoods",
            "local uniqueness boxes with Hessian margins",
            "finite certified kernel or momentum range",
        ],
        "claim_boundary": (
            "exact Morse nondegeneracy proves existence of some sufficiently "
            "small neighborhood abstractly; this receipt supplies no explicit "
            "radius or finite non-pure-kernel range"
        ),
    }

    response = copy.deepcopy(legacy["band_response"])
    response["response_provenance"] = {
        "antipode_J": "derived exactly from the carrier incidence geometry",
        "representative_R_minus_J": (
            "selected by the pinned source-bound issue-599 impulse/readback artifact; "
            "the common overall sign is conventional"
        ),
        "semantic_artifact_self_hash": semantic["artifact_sha256"],
    }

    receipt = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 654,
        "supersedes": {
            "schema": legacy["schema"],
            "status": legacy["status"],
            "path": "code/a5_fingerprint/runtime/a5_multipole_fixed_point_receipt.json",
            "reason": (
                "append-only hardening; the v1 bytes remain immutable because FZ-11 "
                "pins them"
            ),
        },
        "parent_pins": [
            {
                "path": "code/a5_fingerprint/runtime/a5_multipole_fixed_point_receipt.json",
                "bytes": len(legacy_raw),
                "sha256": base.tagged_sha256(legacy_raw),
                "receipt_sha256": legacy["receipt_sha256"],
            },
            {
                "path": "code/a5_closure/manifests/charged_response_semantic_artifact.json",
                "bytes": len(semantic_raw),
                "sha256": base.tagged_sha256(semantic_raw),
                "artifact_sha256": semantic["artifact_sha256"],
            },
        ],
        "normalization": legacy["normalization"],
        "cartesian_frame": cartesian,
        "pole_frame": copy.deepcopy(legacy["pole_frame"]),
        "critical_points": critical,
        "band_response": response,
        "kernel_independence": copy.deepcopy(legacy["kernel_independence"]),
        "kinetic_stencil_conditional": copy.deepcopy(
            legacy["kinetic_stencil_conditional"]
        ),
        "fail_closed_controls": hardened_controls(
            verts, legacy["fail_closed_controls"]
        ),
        "decision_rules_and_local_rows": local_rows(
            legacy["decision_rules_and_ledger"]
        ),
    }
    receipt["receipt_sha256"] = base.tagged_sha256(
        base.canonical_json_bytes(
            {key: value for key, value in receipt.items() if key != "receipt_sha256"}
        )
    )
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
    print(json.dumps(receipt["critical_points"]["original_equation_replay"], indent=1))
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
