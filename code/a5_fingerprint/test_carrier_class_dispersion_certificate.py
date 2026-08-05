"""Adversarial tests for the carrier-class dispersion band certificate.

The high-precision checks rebuild the physical claims from the raw cosine
symbol with mpmath, independently of the exact polynomial pipeline: the
dispersion coefficients of concrete members are extracted by small-k series
fitting and compared against the certified class map.
"""

from __future__ import annotations

import copy
import json
import math
from fractions import Fraction
from pathlib import Path

import pytest

import a5_multipole_fixed_point_certificate as base
import carrier_class_dispersion_certificate as cert

HERE = Path(__file__).resolve().parent


def q5_float(x) -> float:
    a, b = x
    return float(a) + float(b) * math.sqrt(5.0)


@pytest.fixture(scope="module")
def receipt() -> dict:
    return json.loads(cert.RECEIPT_PATH.read_text())


@pytest.fixture(scope="module")
def rebuilt() -> dict:
    return cert.build_receipt()


def test_receipt_matches_fresh_rebuild(receipt, rebuilt):
    assert receipt == rebuilt


def test_receipt_self_digest(receipt):
    body = {k: v for k, v in receipt.items() if k != "receipt_sha256"}
    assert receipt["receipt_sha256"] == base.tagged_sha256(
        base.canonical_json_bytes(body)
    )


def test_mutated_receipt_fails_digest(receipt):
    mutated = copy.deepcopy(receipt)
    mutated["kernel_factorization"]["universal_constant"] = "64/34"
    body = {k: v for k, v in mutated.items() if k != "receipt_sha256"}
    assert mutated["receipt_sha256"] != base.tagged_sha256(
        base.canonical_json_bytes(body)
    )


def test_status_and_schema(receipt):
    assert receipt["schema"] == cert.SCHEMA
    assert receipt["status"] == cert.STATUS
    assert "ISOTROPIC_FLOOR_10_21" in receipt["status"]
    assert "full spatial symbol and hence complete through k^8" in receipt[
        "class_definition"
    ]
    assert "declared positive-weight scalar cosine class" in receipt[
        "kill_surface"
    ]


# ---------------------------------------------------------------------------
# Independent numerical route: series-fit the raw cosine symbol
# ---------------------------------------------------------------------------


def q5_mp(x, mp):
    a, b = x
    return mp.mpf(a.numerator) / a.denominator + (
        mp.mpf(b.numerator) / b.denominator
    ) * mp.sqrt(5)


def unit_orbit_mp(mp) -> dict[str, list[tuple]]:
    """Orbit directions at working mpmath precision, built from the exact
    Q(sqrt5) coordinates."""

    import spin_six_universality_certificate as universality

    orbits = universality.orbit_directions()
    out = {}
    for name, data in orbits.items():
        dirs, norm_sq = data["dirs"], data["norm_sq"]
        norm = mp.sqrt(q5_mp(norm_sq, mp))
        out[name] = [
            tuple(q5_mp(c, mp) / norm for c in v) for v in dirs
        ]
    return out


def i6_mp(direction, mp):
    """I6 from its definition (25/132) sum_v P6(v_hat . n), independently
    of the polynomial pipeline."""

    orbits = unit_orbit_mp(mp)
    n = [mp.mpf(c) for c in direction]
    norm = mp.sqrt(sum(c * c for c in n))
    n = [c / norm for c in n]

    def p6(t):
        return (231 * t**6 - 315 * t**4 + 105 * t**2 - 5) / 16

    total = sum(
        p6(v[0] * n[0] + v[1] * n[1] + v[2] * n[2])
        for v in orbits["vertex_12"]
    )
    return 25 * total / 132


def fit_dispersion(member, direction, mp):
    """Fit k^2, k^4, k^6 coefficients of the normalized symbol at small k.

    member: list of (weight, radius, orbit_name) shells.
    Returns (c2, c4, c6) with the k^2 coefficient normalized to one.
    """

    orbits = unit_orbit_mp(mp)
    n = [mp.mpf(c) for c in direction]
    norm = mp.sqrt(sum(c * c for c in n))
    n = [c / norm for c in n]

    def symbol(k):
        total = mp.mpf(0)
        for weight, radius, orbit in member:
            w = mp.mpf(weight)
            r = mp.mpf(radius)
            for u in orbits[orbit]:
                dot = r * (u[0] * n[0] + u[1] * n[1] + u[2] * n[2])
                total += w * (1 - mp.cos(k * dot))
        return total

    # Solve the 3x3 Vandermonde system for the even series coefficients
    # from three tiny sample points; with mp.dps = 80 and k ~ 1e-8 the
    # k^8 truncation leak on the fitted c4 is below 1e-32.
    ks = [mp.mpf("1e-8"), mp.mpf("2e-8"), mp.mpf("3e-8")]
    rows = [[k**2, k**4, k**6] for k in ks]
    values = [symbol(k) for k in ks]
    sol = mp.lu_solve(mp.matrix(rows), mp.matrix(values))
    c2, c4, c6 = sol[0], sol[1], sol[2]
    return c4 / c2, c6 / c2


def fit_normalized_even_coefficients(member, direction, mp, max_power=8):
    """Fit the normalized even series through ``max_power`` from raw cosines."""

    orbits = unit_orbit_mp(mp)
    n = [mp.mpf(c) for c in direction]
    norm = mp.sqrt(sum(c * c for c in n))
    n = [c / norm for c in n]

    def as_mp(value):
        if isinstance(value, Fraction):
            return mp.mpf(value.numerator) / value.denominator
        return mp.mpf(value)

    def symbol(k):
        total = mp.mpf(0)
        for weight, radius, orbit in member:
            w = as_mp(weight)
            r = as_mp(radius)
            for u in orbits[orbit]:
                dot = r * (u[0] * n[0] + u[1] * n[1] + u[2] * n[2])
                total += w * (1 - mp.cos(k * dot))
        return total

    terms = max_power // 2
    ks = [mp.mpf(index) * mp.mpf("1e-9") for index in range(1, terms + 1)]
    rows = [[k ** (2 * order) for order in range(1, terms + 1)] for k in ks]
    values = [symbol(k) for k in ks]
    solution = mp.lu_solve(mp.matrix(rows), mp.matrix(values))
    return tuple(solution[index] / solution[0] for index in range(1, terms))


def test_vertex_branch_series_fit_matches_fz11():
    mp = pytest.importorskip("mpmath").mp
    mp.dps = 80
    member = [(1.0, 1.0, "vertex_12")]
    direction = (0.3, -0.7, 0.53)
    c4, c6 = fit_dispersion(member, direction, mp)
    # Normalized symbol: divide by the k^2 coefficient, so C4 = c4/c2 with
    # a = 1: expect -1/20 and 1/840 + (2/7875) I6(n).
    assert abs(c4 - mp.mpf(-1) / 20) < mp.mpf("1e-30")
    expected_c6 = mp.mpf(1) / 840 + mp.mpf(2) / 7875 * mp.mpf(
        i6_mp(direction, mp)
    )
    assert abs(c6 - expected_c6) < mp.mpf("1e-12")


def test_edge_branch_series_fit_matches_fz12():
    mp = pytest.importorskip("mpmath").mp
    mp.dps = 80
    member = [(1.0, 1.0, "edge_30")]
    direction = (0.9, 0.1, -0.41)
    c4, c6 = fit_dispersion(member, direction, mp)
    assert abs(c4 - mp.mpf(-1) / 20) < mp.mpf("1e-30")
    expected_c6 = mp.mpf(1) / 840 - mp.mpf(1) / 12600 * mp.mpf(
        i6_mp(direction, mp)
    )
    assert abs(c6 - expected_c6) < mp.mpf("1e-12")


def test_face_branch_series_fit():
    mp = pytest.importorskip("mpmath").mp
    mp.dps = 80
    member = [(1.0, 1.0, "face_20")]
    direction = (-0.2, 0.8, 0.4)
    c4, c6 = fit_dispersion(member, direction, mp)
    assert abs(c4 - mp.mpf(-1) / 20) < mp.mpf("1e-30")
    expected_c6 = mp.mpf(1) / 840 - mp.mpf(2) / 14175 * mp.mpf(
        i6_mp(direction, mp)
    )
    assert abs(c6 - expected_c6) < mp.mpf("1e-12")


def test_two_shell_mixture_obeys_floor_strictly():
    """A concrete two-shell member: the fitted isotropic ratio equals the
    certified moment formula and sits strictly above 10/21."""

    mp = pytest.importorskip("mpmath").mp
    mp.dps = 80
    member = [(1.0, 1.0, "vertex_12"), (0.5, 2.0, "edge_30")]
    # Two directions separate the isotropic and I6 parts of c6.
    dir_a = (0.0, 0.0, 1.0)
    dir_b = (0.62, -0.33, 0.71)
    c4_a, c6_a = fit_dispersion(member, dir_a, mp)
    c4_b, c6_b = fit_dispersion(member, dir_b, mp)
    assert abs(c4_a - c4_b) < mp.mpf("1e-25")
    i6_a = i6_mp(dir_a, mp)
    i6_b = i6_mp(dir_b, mp)
    b0 = (c6_a * i6_b - c6_b * i6_a) / (i6_b - i6_a)
    ratio = b0 / c4_a**2
    # Moment prediction: weights per direction over shells
    # (w |O| r^m): mu_m = 12 * 1 + 0.5 * 30 * 2^m.
    mu2 = 12 + 15 * 4.0
    mu4 = 12 + 15 * 16.0
    mu6 = 12 + 15 * 64.0
    predicted = Fraction(10, 21) * Fraction(mu2 * mu6 / mu4**2)
    assert abs(ratio - mp.mpf(float(predicted))) < mp.mpf("1e-10")
    assert ratio > mp.mpf(10) / 21 + mp.mpf("0.01")
    assert c4_a < 0
    # Class-wide band membership through the same independent route.
    b6 = (c6_a - b0) / i6_a
    band_ratio = b6 / b0
    assert mp.mpf(-16) / 135 < band_ratio < mp.mpf(16) / 75


def test_eighth_order_by_four_term_series_fit():
    """Fit four even coefficients of the raw vertex symbol and confirm the
    eighth-order pair (D0, D6) and the 12/5 cross-order lock."""

    mp = pytest.importorskip("mpmath").mp
    mp.dps = 140
    orbits = unit_orbit_mp(mp)
    direction = (0.3, -0.7, 0.53)
    n = [mp.mpf(c) for c in direction]
    norm = mp.sqrt(sum(c * c for c in n))
    n = [c / norm for c in n]

    def symbol(k):
        total = mp.mpf(0)
        for u in orbits["vertex_12"]:
            dot = u[0] * n[0] + u[1] * n[1] + u[2] * n[2]
            total += 1 - mp.cos(k * dot)
        return total

    ks = [mp.mpf("1e-9"), mp.mpf("2e-9"), mp.mpf("3e-9"), mp.mpf("4e-9")]
    rows = [[k**2, k**4, k**6, k**8] for k in ks]
    values = [symbol(k) for k in ks]
    sol = mp.lu_solve(mp.matrix(rows), mp.matrix(values))
    c2, c4, c6, c8 = sol[0], sol[1], sol[2], sol[3]
    i6 = i6_mp(direction, mp)
    # Normalized: divide by c2.
    assert abs(c4 / c2 - mp.mpf(-1) / 20) < mp.mpf("1e-25")
    expected_c6 = mp.mpf(1) / 840 + mp.mpf(2) / 7875 * i6
    assert abs(c6 / c2 - expected_c6) < mp.mpf("1e-14")
    expected_c8 = -mp.mpf(1) / 60480 - mp.mpf(1) / 118125 * i6
    assert abs(c8 / c2 - expected_c8) < mp.mpf("1e-20")
    # Cross-order lock at the vertex point: (D6/D0)/(B6/B0) = 12/5.
    lock = (mp.mpf(-1) / 118125 / (mp.mpf(-1) / 60480)) / (
        mp.mpf(2) / 7875 / (mp.mpf(1) / 840)
    )
    assert abs(lock - mp.mpf(12) / 5) < mp.mpf("1e-30")


def test_eighth_order_receipt_content(receipt):
    eighth = receipt["eighth_order"]
    assert eighth["general_isotropic_coefficient"] == (
        "D0 = -(a^6/60480)(mu8/mu2)"
    )
    assert eighth["general_anisotropic_ratio"] == (
        "D6/D0 = (64/125) <I6(seed)> with positive weights w |O| r^8"
    )
    assert eighth["cross_order_lock"] == "12/5"
    assert eighth["cross_order_ratio_identity"] == "D6/D0 = (12/5)(B6/B0)"
    assert eighth["cross_order_identity"] == "5 D6 B0 = 12 B6 D0"
    assert eighth["zero_mixture_control"] == {
        "member": "vertex weight 25, face weight 27, one common radius",
        "B6_over_B0": "0",
        "D6_over_D0": "0",
        "division_free_identity": "5 D6 B0 = 12 B6 D0 = 0",
    }
    assert eighth["multi_radius_negative_control"] == {
        "member": "vertex radius 1 plus face radius 2, equal per-direction weights",
        "B6_over_B0": "-25168/218025",
        "D6_over_D0": "-407872/1443375",
        "lock_residual": "-57344/10360225",
        "verdict": "nonzero; the cross-order lock is single-radius only",
    }
    assert eighth["eighth_order_band"] == ["-64/225", "64/125"]
    assert eighth["kernel_universal_constant"] == "256/75"
    assert eighth["harmonic_invariant_multiplicities_l0_l2_l4_l6_l8"] == [
        1,
        0,
        0,
        1,
        0,
    ]
    rows = {row["orbit"]: row for row in eighth["branches"]}
    assert rows["vertex_12"]["D6_over_D0"] == "64/125"
    assert rows["edge_30"]["D6_over_D0"] == "-4/25"
    assert rows["face_20"]["D6_over_D0"] == "-64/225"
    assert all(
        row["radius_scope"] ==
        "unit radius; at radius r, D0/a^6 = -r^6/60480"
        for row in rows.values()
    )


def test_division_free_lock_at_zero_anisotropy_mixture():
    """The quotient of the two anisotropy ratios is 0/0 on the valid
    25:27 member, so the certificate must retain a polynomial identity."""

    numerator = 25 * 12 * Fraction(1) + 27 * 20 * Fraction(-5, 9)
    denominator = 25 * 12 + 27 * 20
    mean_i6 = numerator / denominator
    b6_over_b0 = Fraction(16, 75) * mean_i6
    d6_over_d0 = Fraction(64, 125) * mean_i6
    assert b6_over_b0 == 0
    assert d6_over_d0 == 0
    assert d6_over_d0 == Fraction(12, 5) * b6_over_b0


def test_zero_mixture_raw_symbol_is_isotropic_through_k8():
    """Rebuild the valid 25:27 per-direction mixture from raw cosines.

    Orbit cardinalities are part of the sum, so 25 vertex weight units and
    27 face weight units cancel the rank-six numerator exactly.
    """

    mp = pytest.importorskip("mpmath").mp
    mp.dps = 140
    member = [
        (25.0, 1.0, "vertex_12"),
        (27.0, 1.0, "face_20"),
    ]
    first = (0.31, -0.72, 0.61)
    second = (-0.42, 0.19, 0.88)
    c4_first, c6_first, c8_first = fit_normalized_even_coefficients(
        member, first, mp
    )
    c4_second, c6_second, c8_second = fit_normalized_even_coefficients(
        member, second, mp
    )
    assert abs(c4_first - mp.mpf(-1) / 20) < mp.mpf("1e-30")
    assert abs(c4_second - mp.mpf(-1) / 20) < mp.mpf("1e-30")
    assert abs(c6_first - mp.mpf(1) / 840) < mp.mpf("1e-12")
    assert abs(c6_second - mp.mpf(1) / 840) < mp.mpf("1e-12")
    assert abs(c6_first - c6_second) < mp.mpf("1e-12")
    assert abs(c8_first + mp.mpf(1) / 60480) < mp.mpf("1e-18")
    assert abs(c8_second + mp.mpf(1) / 60480) < mp.mpf("1e-18")
    assert abs(c8_first - c8_second) < mp.mpf("1e-18")


def test_cross_order_lock_is_not_widened_to_multi_radius_members(receipt):
    b6_over_b0 = Fraction(16, 75) * (
        Fraction(12) * 1**6 * 1
        + Fraction(20) * 2**6 * Fraction(-5, 9)
    ) / (Fraction(12) * 1**6 + Fraction(20) * 2**6)
    d6_over_d0 = Fraction(64, 125) * (
        Fraction(12) * 1**8 * 1
        + Fraction(20) * 2**8 * Fraction(-5, 9)
    ) / (Fraction(12) * 1**8 + Fraction(20) * 2**8)
    residual = d6_over_d0 - Fraction(12, 5) * b6_over_b0
    assert residual == Fraction(-57344, 10360225)
    assert residual != 0
    control = receipt["eighth_order"]["multi_radius_negative_control"]
    assert Fraction(control["B6_over_B0"]) == b6_over_b0
    assert Fraction(control["D6_over_D0"]) == d6_over_d0
    assert Fraction(control["lock_residual"]) == residual


def test_common_radius_scales_d0_by_radius_to_the_sixth():
    mp = pytest.importorskip("mpmath").mp
    mp.dps = 140
    radius = Fraction(7, 3)
    member = [(1, radius, "vertex_12")]
    first = (0.37, -0.49, 0.79)
    second = (-0.21, 0.91, 0.35)
    _, _, c8_first = fit_normalized_even_coefficients(member, first, mp)
    _, _, c8_second = fit_normalized_even_coefficients(member, second, mp)
    i6_first = i6_mp(first, mp)
    i6_second = i6_mp(second, mp)
    fitted_d0 = (c8_first * i6_second - c8_second * i6_first) / (
        i6_second - i6_first
    )
    expected_d0 = -(
        mp.mpf(radius.numerator) / radius.denominator
    ) ** 6 / 60480
    assert abs(fitted_d0 - expected_d0) < mp.mpf("1e-15")
    assert Fraction(-16807, 6298560) == -(radius**6) / 60480


def test_kernel_constant_by_float_rotation_sum():
    """Rebuild sum_g ((gu).n)^6 numerically from the certified rotation
    list and confirm the 60/7 + (64/35) I6 I6 factorization."""

    rotations, _ = cert.rotation_group_certificate()
    floats = [
        [[q5_float(entry) for entry in row] for row in g]
        for g in rotations
    ]
    u = (0.31, 0.42, 0.85)
    nrm = math.sqrt(sum(c * c for c in u))
    u = tuple(c / nrm for c in u)
    n = (-0.5, 0.62, 0.61)
    nrm = math.sqrt(sum(c * c for c in n))
    n = tuple(c / nrm for c in n)
    total = 0.0
    for g in floats:
        gu = (
            g[0][0] * u[0] + g[0][1] * u[1] + g[0][2] * u[2],
            g[1][0] * u[0] + g[1][1] * u[1] + g[1][2] * u[2],
            g[2][0] * u[0] + g[2][1] * u[1] + g[2][2] * u[2],
        )
        total += (gu[0] * n[0] + gu[1] * n[1] + gu[2] * n[2]) ** 6
    import mpmath

    mp = mpmath.mp
    mp.dps = 30
    expected = 60.0 / 7.0 + (64.0 / 35.0) * float(i6_mp(u, mp)) * float(
        i6_mp(n, mp)
    )
    assert abs(total - expected) < 1e-9


def test_eighth_order_kernel_constant_by_float_rotation_sum():
    """Independently evaluate the degree-eight kernel factorization."""

    rotations, _ = cert.rotation_group_certificate()
    floats = [
        [[q5_float(entry) for entry in row] for row in rotation]
        for rotation in rotations
    ]
    u = (0.24, -0.72, 0.65)
    u_norm = math.sqrt(sum(component * component for component in u))
    u = tuple(component / u_norm for component in u)
    n = (-0.61, 0.14, 0.78)
    n_norm = math.sqrt(sum(component * component for component in n))
    n = tuple(component / n_norm for component in n)
    total = 0.0
    for rotation in floats:
        rotated = (
            rotation[0][0] * u[0] + rotation[0][1] * u[1] + rotation[0][2] * u[2],
            rotation[1][0] * u[0] + rotation[1][1] * u[1] + rotation[1][2] * u[2],
            rotation[2][0] * u[0] + rotation[2][1] * u[1] + rotation[2][2] * u[2],
        )
        total += sum(rotated[index] * n[index] for index in range(3)) ** 8
    import mpmath

    mp = mpmath.mp
    mp.dps = 30
    expected = 20.0 / 3.0 + (256.0 / 75.0) * float(i6_mp(u, mp)) * float(
        i6_mp(n, mp)
    )
    assert abs(total - expected) < 1e-9


# ---------------------------------------------------------------------------
# Receipt content checks
# ---------------------------------------------------------------------------


def test_branch_table_is_the_frozen_pair_plus_face(receipt):
    rows = {
        row["orbit"]: row
        for row in receipt["single_orbit_branches"]["branches"]
    }
    assert rows["vertex_12"]["B6_over_B0"] == "16/75"
    assert rows["edge_30"]["B6_over_B0"] == "-1/15"
    assert rows["face_20"]["B6_over_B0"] == "-16/135"
    for row in rows.values():
        assert row["C4_over_a2"] == "-1/20"
        assert row["B0_over_a4"] == "1/840"
        assert row["B0_over_C4_squared"] == "10/21"


def test_band_endpoints_and_zero(receipt):
    band = receipt["rank_six_band"]
    assert band["band"] == ["-16/135", "16/75"]
    assert band["edge_point"] == "-1/15"
    assert band["census_recomputed"] == 62


def test_signed_control_violates_floor(receipt):
    control = receipt["radial_floor"]["signed_control"]
    ratio = Fraction(control["ratio"])
    assert ratio < Fraction(10, 21)


def test_floor_from_band_fractions():
    assert Fraction(16, 75) * Fraction(-5, 9) == Fraction(-16, 135)
    assert Fraction(16, 75) * Fraction(-5, 16) == Fraction(-1, 15)
    assert Fraction(2, 7875) == Fraction(16, 75) * Fraction(1, 840)
    assert Fraction(-1, 12600) == Fraction(-1, 15) * Fraction(1, 840)
    assert Fraction(-2, 14175) == Fraction(-16, 135) * Fraction(1, 840)


def test_frozen_receipts_untouched():
    """The class certificate reads the frozen receipts and changes nothing:
    both file hashes equal the custody hashes of their register rows."""

    import hashlib

    register = json.loads(
        (HERE.parent.parent / "claims" / "frozen_prediction_register.json")
        .read_text()
    )
    registered = {row["id"]: row for row in register["rows"]}
    for row_id, filename in (
        ("FZ-11", "spin_six_primitive_port_prediction_receipt.json"),
        ("FZ-12", "seam_current_edge_prediction_receipt.json"),
    ):
        frozen = (HERE / "runtime" / filename).read_bytes()
        assert registered[row_id]["content_sha256"] == hashlib.sha256(
            frozen
        ).hexdigest()


def test_frozen_parent_pins_cover_both_branches(receipt):
    parents = {
        row["row"]: row
        for row in receipt["frozen_receipt_crosscheck"]["parent_pins"]
    }
    for row_id, filename in (
        ("FZ-11", "spin_six_primitive_port_prediction_receipt.json"),
        ("FZ-12", "seam_current_edge_prediction_receipt.json"),
    ):
        raw = (HERE / "runtime" / filename).read_bytes()
        assert parents[row_id] == {
            "row": row_id,
            "path": f"code/a5_fingerprint/runtime/{filename}",
            "bytes": len(raw),
            "sha256": base.tagged_sha256(raw),
        }


@pytest.mark.parametrize(
    ("row_id", "section", "field", "replacement"),
    [
        ("FZ-11", "scale_free_relations", "B0_over_C4_squared", "11/21"),
        ("FZ-11", "scale_free_relations", "B6_over_B0", "17/75"),
        ("FZ-11", "scale_free_relations", "B6_over_C4_squared", "31/315"),
        ("FZ-11", "coefficients", "C4_over_a2", "-1/21"),
        ("FZ-11", "coefficients", "B0_over_a4", "1/841"),
        ("FZ-11", "coefficients", "B6_over_a4", "2/7876"),
        ("FZ-12", "scale_free_relations", "B0_over_C4_squared", "11/21"),
        ("FZ-12", "scale_free_relations", "B6_over_B0", "-2/15"),
        ("FZ-12", "scale_free_relations", "B6_over_C4_squared", "-1/63"),
        ("FZ-12", "coefficients", "C4_over_a2", "-1/21"),
        ("FZ-12", "coefficients", "B0_over_a4", "1/841"),
        ("FZ-12", "coefficients", "B6_over_a4", "-1/12601"),
    ],
)
def test_typed_frozen_branch_mutations_fail_closed(
    row_id, section, field, replacement
):
    fz11 = json.loads(
        (HERE / "runtime" / "spin_six_primitive_port_prediction_receipt.json")
        .read_text()
    )
    fz12 = json.loads(
        (HERE / "runtime" / "seam_current_edge_prediction_receipt.json")
        .read_text()
    )
    target = fz11["exact_prediction"] if row_id == "FZ-11" else fz12[
        "conditional_physical_candidate"
    ]
    target[section][field] = replacement
    with pytest.raises(base.FingerprintError):
        cert._validate_frozen_branch_values(fz11, fz12)


def test_tampered_kernel_constant_rejected(monkeypatch):
    """Force a wrong universal constant and confirm the producer fails
    closed."""

    original = cert.kernel_factorization_certificate

    def tampered(rotations, i6_raw, i6_reduced):
        d, payload = original(rotations, i6_raw, i6_reduced)
        wrong = base.q5_mul(d, base.q5(Fraction(2)))
        return wrong, payload

    monkeypatch.setattr(
        cert, "kernel_factorization_certificate", tampered
    )
    with pytest.raises(base.FingerprintError):
        cert.build_receipt()
