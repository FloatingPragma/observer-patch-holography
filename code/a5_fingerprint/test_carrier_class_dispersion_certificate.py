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
    the FZ-11 file hash equals the custody hash of the register row."""

    frozen = (
        HERE / "runtime" / "spin_six_primitive_port_prediction_receipt.json"
    ).read_bytes()
    import hashlib

    digest = hashlib.sha256(frozen).hexdigest()
    register = json.loads(
        (HERE.parent.parent / "claims" / "frozen_prediction_register.json")
        .read_text()
    )
    row = next(r for r in register["rows"] if r["id"] == "FZ-11")
    assert row["content_sha256"] == digest


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
