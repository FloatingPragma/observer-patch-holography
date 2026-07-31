#!/usr/bin/env python3
"""Issue #655, universality half: spin six is the least populatable spin.

The certificate certifies the invariant table and the orbit moment sums
in exact arithmetic; the universality statement follows by the recorded
decomposition argument. The unconditional content, with no stencil
selection premise:

* **Invariant-theory core.** For the sixty-element proper icosahedral
  rotation group, the dimension of the invariant subspace of the
  spin-``L`` harmonic space is
  ``m_L = (1/60)[(2L+1) + 15 U_{2L}(0) + 20 U_{2L}(1/2)
  + 12 U_{2L}(phi/2) + 12 U_{2L}((phi-1)/2)]``,
  with ``U`` the Chebyshev polynomials of the second kind evaluated over
  ``Q(sqrt5)``. The certificate computes the exact table and certifies
  ``m_1 = ... = m_5 = 0`` and ``m_6 = 1``.
* **Universality.** Any finite-range invariant scalar cosine symbol
  ``lambda_a(k) = a^-2 sum_d c(d) [1 - cos(a k . d)]`` with an invariant
  direction multiset and scale-independent invariant coefficients has, at
  order
  ``a^{2m-2} k^{2m}``, an angular part that is an invariant polynomial of
  degree ``2m``. The invariant table forces exact isotropy through spin
  five at every order; every directional term below spin ten is one
  multiple of the same normalized ``I6``; and the possible artifact
  spins are exactly the even invariant levels ``{6, 10, 12, 16, ...}``,
  so spin six is the least symmetry-allowed nonzero anisotropic spin.
  Whenever the
  weighted sixth-moment ``I6`` coefficient is nonzero, a condition
  certified here for the equal-weight stencil and each fundamental
  orbit, the first artifact sits at ``a^4 k^6`` and one binary
  refinement step suppresses it by exactly ``1/16`` at that order. A
  tuned vanishing coefficient moves the first artifact to higher order
  in ``a``; for members whose directions share one radius, such as the
  exhibited control, the spin-six content cancels at every order and the
  first artifact moves to a higher even invariant spin.
* **Constructive verification.** The three fundamental orbits (twelve
  vertex, twenty face, thirty edge directions) are certified explicitly:
  second and fourth moment sums isotropic, sixth moment sums equal to an
  isotropic constant plus an orbit-dependent multiple of the same ``I6``
  pinned by the issue #654 certificate.

Boundary. The theorem concerns the registered finite carrier class under
an invariance, finite-range, and self-similar-coefficient premise; the
operator is a finite-difference generator on continuum fields, so no periodic lattice
is invoked (no three-dimensional periodic lattice carries this point
group), and the physical carrier realization, scale, sector assignment,
and residue amplitude are open, with no comparison opened here. Minimal
locally Lorentz-invariant Standard Model plus General Relativity contains
no intrinsic preferred spatial tensor in local vacuum. A physical
comparison still needs a sector map, a carrier frame, orientation
transport, an environmental model, and an amplitude floor before a null
can reject the branch.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import factorial
from pathlib import Path
from typing import Any

import a5_multipole_fixed_point_certificate as base

HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "spin_six_universality_receipt.json"

SCHEMA = "oph.spin_six_universality_receipt.v1"
STATUS = (
    "SPIN_SIX_LEAST_POPULATABLE_SPIN__I6_RIGID_BELOW_SPIN_TEN__"
    "GENERIC_LEADING_ORDER__AMPLITUDE_OPEN"
)


def load_parent_receipt() -> tuple[bytes, dict[str, Any]]:
    """Load the fingerprint receipt with schema, status, digest, and bytes checked."""

    parent_bytes = base.RECEIPT_PATH.read_bytes()
    try:
        parent = json.loads(parent_bytes)
    except json.JSONDecodeError as error:
        raise base.FingerprintError("invalid fingerprint parent JSON") from error
    require(
        parent_bytes == base.canonical_json_bytes(parent),
        "fingerprint parent is not canonical JSON",
    )
    require(parent.get("schema") == base.SCHEMA, "fingerprint parent schema drift")
    require(parent.get("status") == base.STATUS, "fingerprint parent status drift")
    claimed = parent.get("receipt_sha256")
    body = {key: value for key, value in parent.items() if key != "receipt_sha256"}
    require(
        claimed == base.tagged_sha256(base.canonical_json_bytes(body)),
        "fingerprint parent self-digest drift",
    )
    return parent_bytes, parent

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


def double_factorial(n: int) -> int:
    """Exact double factorial, with the standard 0!! = (-1)!! = 1."""

    require(n >= -1, "double factorial domain drift")
    result = 1
    while n > 0:
        result *= n
        n -= 2
    return result


def monomial_level_six_coefficient(even_power: int) -> Fraction:
    """Coefficient of P6 in t^even_power on [-1, 1]."""

    require(even_power >= 6 and even_power % 2 == 0, "invalid even power")
    return Fraction(
        13 * factorial(even_power),
        double_factorial(even_power - 6)
        * double_factorial(even_power + 7),
    )


def all_order_level_six_certificate() -> dict[str, Any]:
    """Pin the positive P6 coefficient used by the all-order argument."""

    values = {
        power: monomial_level_six_coefficient(power)
        for power in range(6, 18, 2)
    }
    require(values[6] == Fraction(16, 231), "t^6 to P6 coefficient drift")
    require(values[8] == Fraction(64, 495), "t^8 to P6 coefficient drift")
    require(all(value > 0 for value in values.values()), "P6 coefficient sign drift")
    return {
        "formula": (
            "[P6] t^(2m) = 13 (2m)! / ((2m-6)!! (2m+7)!!), "
            "strictly positive for every m >= 3"
        ),
        "checked_even_powers": {
            str(power): str(value) for power, value in values.items()
        },
        "strictly_positive": True,
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


def tuned_cancellation_control(i6_reduced) -> dict[str, Any]:
    """A positive-weight member with zero spin-six content at order k^6.

    Unit-direction weights 1 on the twelve vertex directions and 27/25 on
    the twenty face directions cancel the weighted sixth-moment I6
    coefficient exactly: 64/175 - (27/25)(64/189) = 0. The member
    certifies that the a^4 k^6 leading-order clause needs the
    nonvanishing-coefficient qualifier; its first artifact sits at a
    higher even invariant spin, while the forbidden spins one through
    five and the below-spin-ten I6 rigidity stay untouched.
    """

    vertex_multiple = q5(Fraction(64, 175))
    face_multiple = q5(Fraction(-64, 189))
    weight = q5(Fraction(27, 25))
    cancelled = q5_add(vertex_multiple, q5_mul(weight, face_multiple))
    require(cancelled == ZERO, "tuned cancellation drift")

    # single-radius proportionality lemma: for a unit orbit O the
    # spin-six content of sum_{d in O} (d.n)^{2m} is c_{2m,6} s_O with
    # one orbit scalar s_O = the I6 multiple of sum_{d in O} P6(d.n)
    # and c_{2m,6} > 0, so cancelling s_v + (27/25) s_f kills spin six
    # at every order for this member. Certify the orbit scalars.
    orbits_lemma = orbit_directions()
    i6_probe = next(m for m in i6_reduced if sum(m) > 0)
    scalars = {}
    for name in ("vertex_12", "face_20"):
        data = orbits_lemma[name]
        inv2 = q5_div(ONE, data["norm_sq"])
        total = p_zero()
        for power, coeff in enumerate(base.LEGENDRE[6]):
            if coeff == 0:
                continue
            raw = moment_sum(data["dirs"], power)
            if power % 2 == 1:
                require(p_is_zero(raw), "odd moment leak in lemma")
                continue
            total = p_add(
                total,
                p_scale_q5(raw, q5_scale(q5_pow(inv2, power // 2), coeff)),
            )
        reduced = p_reduce_sphere(total)
        scalars[name] = q5_div(reduced.get(i6_probe, ZERO), i6_reduced[i6_probe])
    require(scalars["vertex_12"] == q5(Fraction(132, 25)), "s_v drift")
    require(scalars["face_20"] == q5(Fraction(-44, 9)), "s_f drift")
    require(
        q5_add(scalars["vertex_12"], q5_mul(weight, scalars["face_20"]))
        == ZERO,
        "every-order cancellation scalar drift",
    )

    # independent recomputation from the orbit sums
    orbits = orbit_directions()
    probe = next(m for m in i6_reduced if sum(m) > 0)
    total = p_zero()
    for name, w in (("vertex_12", Fraction(1)), ("face_20", Fraction(27, 25))):
        data = orbits[name]
        inv2 = q5_div(ONE, data["norm_sq"])
        m6 = p_scale_q5(
            moment_sum(data["dirs"], 6), q5_mul(inv2, q5_mul(inv2, inv2))
        )
        total = p_add(total, p_scale(m6, w))
    reduced = p_reduce_sphere(total)
    beta = q5_div(reduced.get(probe, ZERO), i6_reduced[probe])
    require(beta == ZERO, "tuned member retains spin-six content at k^6")
    residue = p_add(reduced, p_scale_q5(i6_reduced, q5_neg(beta)))
    residue.pop((0, 0, 0), None)
    require(p_is_zero(residue), "tuned member leaks below spin ten at k^6")

    # The same-radius member has a nonzero spin-ten harmonic.  Each P10
    # orbit sum is already harmonic, so a nonzero reduced polynomial proves
    # that the first allowed anisotropic spin after the cancelled I6 line is
    # actually populated.
    p10_total = p_zero()
    for name, w in (("vertex_12", Fraction(1)), ("face_20", Fraction(27, 25))):
        data = orbits[name]
        inv2 = q5_div(ONE, data["norm_sq"])
        orbit_p10 = p_zero()
        for power, coeff in enumerate(base.LEGENDRE[10]):
            if coeff == 0:
                continue
            raw = moment_sum(data["dirs"], power)
            if power % 2 == 1:
                require(p_is_zero(raw), "odd moment leak in P10 control")
                continue
            orbit_p10 = p_add(
                orbit_p10,
                p_scale_q5(
                    raw,
                    q5_scale(q5_pow(inv2, power // 2), coeff),
                ),
            )
        p10_total = p_add(p10_total, p_scale(orbit_p10, w))
    p10_reduced = p_reduce_sphere(p10_total)
    require(not p_is_zero(p10_reduced), "tuned member loses its spin-ten line")
    p10_probe = sorted(p10_reduced)[0]

    # A multi-radius control shows why a sixth-moment cancellation need not
    # remove I6 at every higher order.  Vertex radius 1 with weight 1 and
    # face radius 2 with weight 27/1600 cancel the k^6 I6 coefficient, while
    # the k^8 coefficient remains -256/125.
    beta_by_orbit_power: dict[tuple[str, int], Q5] = {}
    for name in ("vertex_12", "face_20"):
        data = orbits[name]
        inv2 = q5_div(ONE, data["norm_sq"])
        for power in (6, 8):
            normalized = p_scale_q5(
                moment_sum(data["dirs"], power),
                q5_pow(inv2, power // 2),
            )
            reduced_power = p_reduce_sphere(normalized)
            beta_by_orbit_power[(name, power)] = q5_div(
                reduced_power.get(probe, ZERO), i6_reduced[probe]
            )
    require(
        beta_by_orbit_power[("vertex_12", 8)] == q5(Fraction(256, 375)),
        "vertex eighth-moment I6 coefficient drift",
    )
    require(
        beta_by_orbit_power[("face_20", 8)] == q5(Fraction(-256, 405)),
        "face eighth-moment I6 coefficient drift",
    )
    multi_weight = q5(Fraction(27, 1600))
    multi_beta6 = q5_add(
        beta_by_orbit_power[("vertex_12", 6)],
        q5_mul(q5_scale(multi_weight, 2**6), beta_by_orbit_power[("face_20", 6)]),
    )
    multi_beta8 = q5_add(
        beta_by_orbit_power[("vertex_12", 8)],
        q5_mul(q5_scale(multi_weight, 2**8), beta_by_orbit_power[("face_20", 8)]),
    )
    require(multi_beta6 == ZERO, "multi-radius k^6 cancellation drift")
    require(multi_beta8 == q5(Fraction(-256, 125)), "multi-radius k^8 control drift")

    return {
        "member": (
            "unit directions, weight 1 on the twelve vertex directions and "
            "27/25 on the twenty face directions"
        ),
        "weights_positive": True,
        "k6_i6_coefficient": "0 (exact cancellation 64/175 - (27/25)(64/189))",
        "single_radius_lemma": (
            "for a unit orbit the spin-six content of every moment is one "
            "positive multiple of the orbit scalar s_O, the I6 multiple "
            "of the orbit's P6 sum; certified s_v = 132/25, s_f = -44/9, "
            "and s_v + (27/25) s_f = 0, so this member's spin-six content "
            "vanishes at every order and its first artifact sits at a "
            "higher even invariant spin"
        ),
        "spin_ten_control": {
            "nonzero": True,
            "probe_monomial": "-".join(str(value) for value in p10_probe),
            "probe_coefficient": q5_str(p10_reduced[p10_probe]),
            "reading": (
                "the same-radius positive tuned member has a nonzero P10 "
                "orbit sum after its I6 line is cancelled"
            ),
        },
        "multi_radius_control": {
            "member": (
                "vertex radius 1 with weight 1 plus face radius 2 with "
                "weight 27/1600"
            ),
            "k6_i6_coefficient": q5_str(multi_beta6),
            "k8_i6_coefficient": q5_str(multi_beta8),
            "reading": (
                "a multi-radius member can cancel I6 at k^6 while "
                "retaining it at k^8"
            ),
        },
        "reading": (
            "the leading-order and 1/16 clauses hold under the certified "
            "nonvanishing-coefficient condition and fail on tuned members "
            "such as this one; the forbidden spins and the below-spin-ten "
            "I6 rigidity hold for every member"
        ),
    }


# ---------------------------------------------------------------------------
# The universality theorem and its boundary
# ---------------------------------------------------------------------------


def universality_statement() -> dict[str, Any]:
    return {
        "carrier_class": (
            "every scalar cosine symbol lambda_a(k) = a^-2 sum_d "
            "c(d)[1 - cos(a k.d)] with a finite dimensionless direction "
            "multiset and scale-independent coefficients invariant under "
            "the proper icosahedral rotation group; nonnegative coefficients "
            "and the second-moment normalization make it a positive normalized "
            "kinetic symbol"
        ),
        "argument": (
            "the order a^{2m-2} k^{2m} angular content is an invariant "
            "polynomial of degree 2m; its harmonic components sit at even "
            "levels l <= 2m; the exact invariant table forces zero at "
            "l = 2, 4, 8 and a one-dimensional space at l = 6, so the class "
            "is exactly isotropic through spin five at every order, every "
            "directional term below spin ten is one multiple of the same "
            "normalized I6, and the possible artifact spins are exactly "
            "the even invariant levels {6, 10, 12, 16, ...}"
        ),
        "consequences": {
            "least_populatable_spin": (
                "spin six is the least symmetry-allowed nonzero anisotropic "
                "spin, and every anisotropic harmonic component below spin ten "
                "carries the exact I6 angular dependence with the "
                "sign-symmetric 62-point census of the issue #654 "
                "certificate: 12 and 20 extrema of opposite index and 30 "
                "saddles"
            ),
            "clean_low_levels": (
                "zero anisotropic content at l = 1..5 at every order in a, "
                "for every member"
            ),
            "generic_leading_order": (
                "whenever the weighted sixth-moment I6 coefficient is "
                "nonzero, certified here for the equal-weight stencil and "
                "each fundamental orbit, the first artifact sits at "
                "a^4 k^6 and one binary refinement step suppresses it by "
                "exactly 1/16 at fixed physical k for the same dimensionless "
                "stencil under binary refinement; a tuned vanishing coefficient "
                "moves the first artifact to higher order in a, and for "
                "single-radius members to a higher even invariant spin"
            ),
            "selection_premise_weakened": (
                "the issue #654 stencil receipt's declared equal-weight "
                "premise weakens to invariance plus finite range; stencil "
                "selection affects the residue amplitude and, on tuned "
                "members, the leading order, never the forbidden levels or "
                "the below-spin-ten I6 rigidity"
            ),
        },
        "baseline_contrast": (
            "minimal locally Lorentz-invariant Standard Model plus General "
            "Relativity has no intrinsic preferred spatial tensor in local "
            "vacuum and therefore supplies no coefficient on this template; "
            "nonminimal effective operators, media, source effects, gravity, "
            "and instruments can imitate or contaminate an anisotropy and must "
            "be included in any physical comparison"
        ),
        "open_premises": {
            "physical_carrier_realization": (
                "a physical carrier with icosahedrally invariant "
                "finite-range kinetics at finite scale a, in a "
                "quasiperiodic or graph realization; no periodic "
                "three-dimensional lattice carries this point group, and "
                "the realization with its propagation-sector assignment is "
                "owned by the source lanes"
            ),
            "sector_and_frame": (
                "a scalar or polarization-independent physical sector, a "
                "carrier rest frame, and coherent orientation transport into "
                "the comparison frame are not established by this theorem"
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
    all_order = all_order_level_six_certificate()
    orbits = orbit_moment_certificate(i6_reduced)
    tuned = tuned_cancellation_control(i6_reduced)
    statement = universality_statement()
    parent, _ = load_parent_receipt()
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
        "all_order_level_six_coefficient": all_order,
        "orbit_verification": orbits,
        "tuned_cancellation_control": tuned,
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
