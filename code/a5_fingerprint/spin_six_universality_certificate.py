#!/usr/bin/env python3
"""Issue #655, universality half: the spin-six residue is carrier-universal.

The certificate proves, in exact arithmetic, that the leading directional
artifact of icosahedral carrier kinetics is universal in shape and order,
with no stencil selection premise:

* **Invariant-theory core.** For the sixty-element proper icosahedral
  rotation group, the dimension of the invariant subspace of the
  spin-``L`` harmonic space is
  ``m_L = (1/60)[(2L+1) + 15 U_{2L}(0) + 20 U_{2L}(1/2)
  + 12 U_{2L}(phi/2) + 12 U_{2L}((phi-1)/2)]``,
  with ``U`` the Chebyshev polynomials of the second kind evaluated over
  ``Q(sqrt5)``. The certificate computes the exact table and certifies
  ``m_1 = ... = m_5 = 0`` and ``m_6 = 1``.
* **Universality.** Any finite-range invariant kinetic operator
  ``lambda_a(k) = sum_d w(d) [1 - cos(a k . d)]`` with an invariant
  direction multiset and invariant weights has, at order
  ``a^{2m-2} k^{2m}``, an angular part that is an invariant polynomial of
  degree ``2m``. For ``2m = 2, 4`` the invariant table forces isotropy;
  for ``2m = 6`` the single invariant line forces every carrier's first
  directional artifact to be one multiple of the same normalized ``I6``.
  The first anisotropy therefore sits at ``a^4 k^6`` for every member of
  the class, and one binary refinement step suppresses it by exactly
  ``1/16`` at that order.
* **Constructive verification.** The three fundamental orbits (twelve
  vertex, twenty face, thirty edge directions) are certified explicitly:
  second and fourth moment sums isotropic, sixth moment sums equal to an
  isotropic constant plus an orbit-dependent multiple of the same ``I6``
  pinned by the issue #654 certificate.

Boundary. The theorem concerns the registered finite carrier class under
an invariance and finite-range premise; the physical carrier realization,
lattice scale, and residue amplitude are open, and no comparison is
opened here. The Standard Model with General Relativity carries exact
rotational invariance in this sector: it produces no rotational residue
of any shape at any order, so a certified detection of an
``l = 6``-shaped, ``l <= 5``-clean, ``a^4``-scaling residue is outside
that baseline, while a null result bounds the carrier scale.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import a5_multipole_fixed_point_certificate as base

HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "spin_six_universality_receipt.json"

SCHEMA = "oph.spin_six_universality_receipt.v1"
STATUS = "SPIN_SIX_RESIDUE_UNIVERSAL_ON_INVARIANT_CARRIER_CLASS__AMPLITUDE_OPEN"

q5 = base.q5
Q5 = base.Q5
ZERO = base.ZERO
ONE = base.ONE
PHI = base.PHI
q5_add = base.q5_add
q5_sub = base.q5_sub
q5_mul = base.q5_mul
q5_scale = base.q5_scale
q5_div = base.q5_div
q5_pow = base.q5_pow
q5_neg = base.q5_neg
q5_str = base.q5_str
require = base.require
p_zero = base.p_zero
p_add = base.p_add
p_scale = base.p_scale
p_scale_q5 = base.p_scale_q5
p_pow = base.p_pow
p_linear = base.p_linear
p_reduce_sphere = base.p_reduce_sphere
p_is_zero = base.p_is_zero
radial_power = base.radial_power
moment_sum = base.moment_sum


# ---------------------------------------------------------------------------
# Invariant dimensions from exact Chebyshev character sums
# ---------------------------------------------------------------------------


def chebyshev_u(n: int, x: Q5) -> Q5:
    """Chebyshev polynomial of the second kind, exact over Q(sqrt5)."""

    prev, cur = ONE, q5_scale(x, 2)
    if n == 0:
        return prev
    for _ in range(n - 1):
        prev, cur = cur, q5_sub(q5_mul(q5_scale(x, 2), cur), prev)
    return cur


def invariant_dimension(level: int) -> int:
    """dim of icosahedral-rotation invariants in the spin-level space.

    Conjugacy classes of the sixty proper rotations: identity (1),
    fifteen half-turns, twenty order-three rotations, twelve rotations by
    2 pi / 5 and twelve by 4 pi / 5. The spin-L character at rotation
    angle theta is U_{2L}(cos(theta/2)); the half-angle cosines are
    0, 1/2, phi/2, and (phi - 1)/2, all in Q(sqrt5).
    """

    half_cosines = [
        (15, ZERO),
        (20, q5(Fraction(1, 2))),
        (12, q5_scale(PHI, Fraction(1, 2))),
        (12, q5_scale(q5_sub(PHI, ONE), Fraction(1, 2))),
    ]
    total = q5(2 * level + 1)
    for count, cosine in half_cosines:
        total = q5_add(total, q5_scale(chebyshev_u(2 * level, cosine), count))
    require(total[1] == 0, f"character sum leaves Q at level {level}")
    value = total[0] / 60
    require(value.denominator == 1, f"non-integer multiplicity at {level}")
    require(value >= 0, f"negative multiplicity at {level}")
    return int(value)


def invariant_table_certificate() -> dict[str, Any]:
    table = {level: invariant_dimension(level) for level in range(17)}
    require(table[0] == 1, "trivial invariant drift")
    require(
        all(table[level] == 0 for level in range(1, 6)),
        "an invariant exists below level six",
    )
    require(table[6] == 1, "level-six invariant line is not one-dimensional")
    require(table[15] == 1, "first odd invariant is not at level fifteen")
    require(
        [level for level, dim in table.items() if dim > 0]
        == [0, 6, 10, 12, 15, 16],
        "allowed-level list drift on the window",
    )
    return {
        "formula": (
            "m_L = (1/60)[(2L+1) + 15 U_{2L}(0) + 20 U_{2L}(1/2) + "
            "12 U_{2L}(phi/2) + 12 U_{2L}((phi-1)/2)]"
        ),
        "table": {str(level): dim for level, dim in table.items()},
        "no_invariant_below_six": True,
        "level_six_dimension_one": True,
        "first_odd_invariant_level": 15,
    }


# ---------------------------------------------------------------------------
# Constructive orbit verification
# ---------------------------------------------------------------------------


def orbit_directions() -> dict[str, dict[str, Any]]:
    verts = base.cartesian_vertices()
    inv_norm = q5_div(ONE, base.NORM_SQ)

    def unit_dot(i: int, j: int) -> Q5:
        dot = ZERO
        for axis in range(3):
            dot = q5_add(dot, q5_mul(verts[i][axis], verts[j][axis]))
        return q5_mul(dot, inv_norm)

    inv_sqrt5 = base.INV_SQRT5
    faces = []
    for i in range(12):
        for j in range(i + 1, 12):
            if unit_dot(i, j) != inv_sqrt5:
                continue
            for k in range(j + 1, 12):
                if unit_dot(i, k) == inv_sqrt5 and unit_dot(j, k) == inv_sqrt5:
                    faces.append((i, j, k))
    edges = [
        (i, j)
        for i in range(12)
        for j in range(i + 1, 12)
        if unit_dot(i, j) == inv_sqrt5
    ]
    require(len(faces) == 20 and len(edges) == 30, "orbit census drift")

    def vec_sum(indices) -> tuple[Q5, Q5, Q5]:
        total = [ZERO, ZERO, ZERO]
        for idx in indices:
            for axis in range(3):
                total[axis] = q5_add(total[axis], verts[idx][axis])
        return tuple(total)

    face_dirs = [vec_sum(f) for f in faces]
    edge_dirs = [vec_sum(e) for e in edges]
    face_norm = q5_add(q5(6), q5_scale(PHI, 9))
    edge_norm = q5_add(q5(4), q5_scale(PHI, 4))
    for d in face_dirs:
        nsq = ZERO
        for axis in range(3):
            nsq = q5_add(nsq, q5_mul(d[axis], d[axis]))
        require(nsq == face_norm, "face direction norm drift")
    for d in edge_dirs:
        nsq = ZERO
        for axis in range(3):
            nsq = q5_add(nsq, q5_mul(d[axis], d[axis]))
        require(nsq == edge_norm, "edge direction norm drift")

    return {
        "vertex_12": {"dirs": verts, "norm_sq": base.NORM_SQ},
        "face_20": {"dirs": face_dirs, "norm_sq": face_norm},
        "edge_30": {"dirs": edge_dirs, "norm_sq": edge_norm},
    }


def orbit_moment_certificate(i6_reduced) -> dict[str, Any]:
    orbits = orbit_directions()
    probe = next(m for m in i6_reduced if sum(m) > 0)
    rows = []
    for name, data in orbits.items():
        dirs, norm_sq = data["dirs"], data["norm_sq"]
        count = len(dirs)
        for k in (1, 3, 5):
            require(
                p_is_zero(moment_sum(dirs, k)),
                f"odd moment {k} survives on {name}",
            )
        inv2 = q5_div(ONE, norm_sq)
        m2 = p_scale_q5(moment_sum(dirs, 2), inv2)
        m4 = p_scale_q5(moment_sum(dirs, 4), q5_mul(inv2, inv2))
        m6 = p_scale_q5(moment_sum(dirs, 6), q5_mul(inv2, q5_mul(inv2, inv2)))
        c2 = q5_scale(q5(count), Fraction(1, 3))
        require(m2 == p_scale_q5(radial_power(1), c2), f"{name} m2 anisotropic")
        c4 = q5_scale(q5(count), Fraction(1, 5))
        require(m4 == p_scale_q5(radial_power(2), c4), f"{name} m4 anisotropic")
        m6_red = p_reduce_sphere(m6)
        beta = q5_div(m6_red.get(probe, ZERO), i6_reduced[probe])
        residue = p_add(m6_red, p_scale_q5(i6_reduced, q5_neg(beta)))
        alpha = residue.pop((0, 0, 0), ZERO)
        require(
            p_is_zero(residue),
            f"{name} sixth moment leaves the invariant line",
        )
        rows.append(
            {
                "orbit": name,
                "count": count,
                "m2": f"({q5_str(c2)}) r^2",
                "m4": f"({q5_str(c4)}) r^4",
                "m6_isotropic": q5_str(alpha),
                "m6_i6_multiple": q5_str(beta),
                "i6_multiple_nonzero": beta != ZERO,
            }
        )
    require(
        all(row["i6_multiple_nonzero"] for row in rows),
        "an orbit has an accidentally vanishing level-six multiple",
    )
    return {
        "statement": (
            "for each fundamental orbit the second and fourth moment sums "
            "are exactly isotropic and the sixth moment sum is one "
            "isotropic constant plus one nonzero multiple of the same "
            "normalized I6"
        ),
        "orbits": rows,
    }


# ---------------------------------------------------------------------------
# The universality theorem and its boundary
# ---------------------------------------------------------------------------


def universality_statement() -> dict[str, Any]:
    return {
        "carrier_class": (
            "every kinetic operator lambda_a(k) = sum_d w(d)[1 - cos(a k.d)] "
            "with a finite direction multiset and weights invariant under "
            "the proper icosahedral rotation group"
        ),
        "argument": (
            "the order a^{2m-2} k^{2m} angular content is an invariant "
            "polynomial of degree 2m; its harmonic components sit at even "
            "levels l <= 2m; the exact invariant table forces zero at "
            "l = 2, 4 (orders k^4 and the l <= 4 part of k^6) and a "
            "one-dimensional space at l = 6, so every member's first "
            "directional artifact is one multiple of the same normalized "
            "I6 at order a^4 k^6"
        ),
        "consequences": {
            "shape_universality": (
                "the leading rotational residue of every carrier in the "
                "class has the exact I6 angular dependence with the "
                "62-point census of the issue #654 certificate"
            ),
            "clean_low_levels": (
                "zero anisotropic content at l = 1..5 at every order in a"
            ),
            "refinement_law": (
                "the leading residue scales as a^4, so one binary "
                "refinement step suppresses it by exactly 1/16 at that "
                "order, for every member of the class"
            ),
            "selection_premise_weakened": (
                "the issue #654 stencil receipt's declared equal-weight "
                "premise weakens to invariance plus finite range; stencil "
                "selection affects only the residue amplitude, never the "
                "shape, the forbidden levels, or the refinement exponent"
            ),
        },
        "baseline_contrast": (
            "the Standard Model with General Relativity is exactly "
            "rotationally invariant in this sector and produces no "
            "rotational residue of any shape at any order; a certified "
            "detection of an l = 6-shaped, l <= 5-clean, a^4-scaling "
            "residue therefore lies outside that baseline, and a null "
            "result bounds the carrier scale"
        ),
        "open_premises": {
            "physical_carrier_realization": (
                "a physical lattice-type carrier at finite scale a with "
                "invariant finite-range kinetics; owned by the source "
                "realization lanes"
            ),
            "amplitude_and_scale": (
                "no residue amplitude or carrier scale is derived; the "
                "amplitude question is the remaining half of issue #655 "
                "and the branch arithmetic of the kinetic-form dichotomy"
            ),
        },
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
            "eligibility": "INELIGIBLE_OPEN_PHYSICAL_MAP for issue #639",
        },
    }


def build_receipt() -> dict[str, Any]:
    cartesian = base.build_cartesian_frame()
    i6 = cartesian.pop("_i6_poly_object")
    cartesian.pop("_vertices_object")
    i6_reduced = p_reduce_sphere(i6)
    invariants = invariant_table_certificate()
    orbits = orbit_moment_certificate(i6_reduced)
    statement = universality_statement()
    parent = base.RECEIPT_PATH.read_bytes()
    receipt = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 655,
        "parent_pins": [
            {
                "path": (
                    "code/a5_fingerprint/runtime/"
                    "a5_multipole_fixed_point_receipt.json"
                ),
                "bytes": len(parent),
                "sha256": base.tagged_sha256(parent),
            }
        ],
        "invariant_table": invariants,
        "orbit_verification": orbits,
        "universality": statement,
    }
    receipt["receipt_sha256"] = base.tagged_sha256(
        base.canonical_json_bytes(
            {k: v for k, v in receipt.items() if k != "receipt_sha256"}
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
    print(json.dumps(receipt["invariant_table"]["table"], indent=1))
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
