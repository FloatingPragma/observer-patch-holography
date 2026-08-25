"""Exact checks for the integer-k comb template and its receipt.

Run with pytest. The tests exercise the producer's exact limits
(Schwarzschild), the imported ratio ladder, tooth monotonicity, the declared
KMS net-response factor, verifier agreement, receipt determinism, and mutation guards
on the frozen constants ln(k)/(8*pi).

What is not proved here. These tests certify internal consistency of
the build-stage instrument only; they touch no event data and no
comparison dataset, and they do not register a prediction.
"""

from __future__ import annotations

import inspect
import json
import math
import os
import subprocess
import sys
from decimal import Decimal, getcontext

import pytest

import integer_k_comb_template as producer

HERE = os.path.dirname(os.path.abspath(__file__))
RECEIPT = os.path.join(HERE, "runtime", producer.RECEIPT_BASENAME)


def setup_module(_module) -> None:
    getcontext().prec = producer.WORKING_PRECISION


def _rel_close(a: Decimal, b: Decimal, tol: str = "1E-45") -> bool:
    if b == 0:
        return a == 0
    return abs(a - b) / abs(b) < Decimal(tol)


def _load_receipt() -> dict:
    with open(RECEIPT, "rb") as fh:
        return json.loads(fh.read().decode("ascii"))


def test_schwarzschild_limit_omega_h_zero() -> None:
    assert producer.omega_h_si(Decimal(62), Decimal(0)) == 0
    assert producer.omega_h_si(Decimal(1), Decimal(0)) == 0


def test_schwarzschild_limit_kappa() -> None:
    # kappa(M, 0) = c^3/(4*G*M) exactly: check kappa * 4*G*M == c^3.
    c = Decimal(producer.C_LIGHT_M_PER_S)
    for mass in (Decimal(1), Decimal(62), Decimal("10.5")):
        kappa = producer.kappa_si(mass, Decimal(0))
        assert _rel_close(kappa * 4 * producer.gm_si(mass), c ** 3)


def test_g_of_chi_surface_gravity_identity() -> None:
    # Statement-pinned g(chi) equals 4*G*M*kappa/c^3 at every spin.
    c = Decimal(producer.C_LIGHT_M_PER_S)
    mass = Decimal(62)
    for chi_s in ("0", "0.1", "0.3", "0.67", "0.9", "0.998"):
        chi = Decimal(chi_s)
        lhs = producer.g_of_chi(chi)
        rhs = 4 * producer.gm_si(mass) * producer.kappa_si(mass, chi) / c ** 3
        assert _rel_close(lhs, rhs)


def test_integer_division_entropy_sign_and_eligibility() -> None:
    assert producer.integer_division_after(60, 3) == 20
    signed_change, entropy_loss = producer.transition_entropy_nats(60, 3)
    log_three = Decimal(3).ln()
    assert _rel_close(signed_change, -log_three)
    assert _rel_close(entropy_loss, log_three)
    assert signed_change < 0 < entropy_loss
    with pytest.raises(ValueError, match="divide"):
        producer.integer_division_after(10, 3)
    with pytest.raises(ValueError, match="at least 2"):
        producer.integer_division_after(10, 1)


def test_detector_frame_mass_and_frequency_scaling() -> None:
    pi = producer.compute_pi()
    source_mass = Decimal("62")
    redshift = Decimal("0.25")
    chi = Decimal("0.67")
    mass_det = producer.detector_frame_mass_solar(source_mass, redshift)
    assert mass_det == Decimal("77.50")
    source_frequency = producer.tooth_frequency_hz(source_mass, chi, 2, 3, pi)
    detector_frequency = producer.detector_frame_tooth_frequency_hz(
        source_mass, redshift, chi, 2, 3, pi
    )
    assert _rel_close(detector_frequency, source_frequency / (1 + redshift))
    source_rotation = producer.rotation_line_hz(source_mass, chi, 2, pi)
    detector_rotation = producer.rotation_line_hz(mass_det, chi, 2, pi)
    assert _rel_close(detector_rotation, source_rotation / (1 + redshift))
    with pytest.raises(ValueError, match="nonnegative"):
        producer.detector_frame_mass_solar(source_mass, Decimal("-0.1"))


def test_ratio_ladder_values() -> None:
    pi = producer.compute_pi()
    assert float(producer.ladder_ratio(4)) == 2.0
    assert abs(float(producer.ladder_ratio(3)) - math.log(3) / math.log(2)) < 1e-12
    assert abs(float(producer.ladder_ratio(5)) - math.log(5) / math.log(2)) < 1e-12
    # Offset-subtracted ratio equals x_k / x_2 identically.
    for k in range(2, 13):
        xk = producer.universal_position(k, pi)
        x2 = producer.universal_position(2, pi)
        assert _rel_close(xk / x2, producer.ladder_ratio(k))


def test_universal_positions() -> None:
    pi = producer.compute_pi()
    # Statement display values (rounded there to four or five digits).
    statement = {2: 0.02758, 3: 0.04371, 4: 0.05516, 5: 0.06404}
    for k, val in statement.items():
        assert abs(float(producer.universal_position(k, pi)) - val) < 5e-6
    # Independent float path.
    for k in range(2, 13):
        expected = math.log(k) / (8 * math.pi)
        got = float(producer.universal_position(k, pi))
        assert abs(got - expected) < 1e-15 * expected


def test_tooth_monotonicity() -> None:
    pi = producer.compute_pi()
    mass = Decimal(producer.DECLARED_REFERENCE_MASS_SOLAR)
    chi = Decimal(producer.DECLARED_REFERENCE_CHI)
    xs = [producer.universal_position(k, pi) for k in range(2, 13)]
    dfs = [producer.tooth_offset_hz(mass, chi, k, pi) for k in range(2, 13)]
    fs = [
        producer.tooth_frequency_hz(mass, chi, 2, k, pi) for k in range(2, 13)
    ]
    for seq in (xs, dfs, fs):
        for lo, hi in zip(seq, seq[1:]):
            assert lo < hi


def test_kms_net_response_factor_algebra() -> None:
    weights = [producer.kms_weight(k) for k in range(2, 13)]
    for lo, hi in zip(weights, weights[1:]):
        assert lo < hi
    for k, w in zip(range(2, 13), weights):
        assert w < 1
        assert _rel_close(w, Decimal(k - 1) / Decimal(k))
    assert float(weights[0]) == 0.5


def test_linewidth_model() -> None:
    pi = producer.compute_pi()
    # Mass independence is structural, but the spin factor remains.
    params = inspect.signature(producer.linewidth_fraction).parameters
    assert "mass_solar" not in params
    assert list(params) == ["a", "chi", "k", "pi"]
    # The original 1.8%--18% band is the Schwarzschild g(0)=1 limit.
    wide_schw = producer.linewidth_fraction(Decimal(1), Decimal(0), 2, pi)
    narrow_schw = producer.linewidth_fraction(Decimal(10), Decimal(0), 2, pi)
    assert 0.18 < float(wide_schw) < 0.19
    assert 0.018 < float(narrow_schw) < 0.019
    # At the synthetic Kerr spin chi=.67, g^-2 widens the k=2 tooth.
    chi = Decimal(producer.DECLARED_REFERENCE_CHI)
    wide = producer.linewidth_fraction(Decimal(1), chi, 2, pi)
    narrow = producer.linewidth_fraction(Decimal(10), chi, 2, pi)
    assert 0.25 < float(wide) < 0.26
    assert 0.025 < float(narrow) < 0.026
    # Independent float path.
    for k in (2, 3, 5, 12):
        for a in (1, 4, 10):
            g = float(producer.g_of_chi(chi))
            expected = 64 * math.pi ** 2 * 2e-4 / (a * g * g * math.log(k))
            got = float(producer.linewidth_fraction(Decimal(a), chi, k, pi))
            assert abs(got - expected) < 1e-12 * expected


def test_reference_point_cross_check() -> None:
    receipt = _load_receipt()
    display = receipt["derived_for_display"]
    c = 299792458.0
    gm = 62.0 * 1.3271244e20
    chi = 0.67
    s = math.sqrt(1 - chi * chi)
    omega_h = c ** 3 * chi / (2 * gm * (1 + s))
    kappa = c ** 3 * s / (2 * gm * (1 + s))
    g = 2 * s / (1 + s)
    rot = 2 * omega_h / (2 * math.pi)
    base = c ** 3 * g / (16 * math.pi ** 2 * gm)

    def close(key: str, val: float) -> None:
        assert abs(display[key] - val) < 1e-12 * abs(val), key

    close("reference.omega_h_rad_per_s", omega_h)
    close("reference.kappa_per_s", kappa)
    close("reference.g_chi", g)
    close("reference.rotation_line_hz", rot)
    close("reference.base_spacing_hz_per_nat", base)
    for k in range(2, 13):
        close("reference.teeth.k%02d.delta_f_hz" % k, base * math.log(k))
        close("reference.teeth.k%02d.f_hz" % k, rot + base * math.log(k))


def test_mutation_guard_imported_template_constants() -> None:
    receipt = _load_receipt()
    law = receipt["imported_continuation_law"]
    assert law["ratio_law"] == (
        "(f_a - m*Omega_H/(2*pi)) / (f_b - m*Omega_H/(2*pi)) "
        "= ln(k_a)/ln(k_b), integers k >= 2"
    )
    assert law["universal_coordinate"].endswith("x_k = ln(k)/(8*pi)")
    assert law["kms_weight"] == "(k-1)/k"
    assert "not a normalized cross-k" in law["kms_scope"]
    assert law["linewidth_scope"].startswith("Gamma/Delta_E_k")
    assert law["transition_approximation"].startswith("leading small-transition")
    assert law["linewidth_fraction"] == "64*pi^2*p_0/(a*g(chi)^2*ln(k))"
    assert law["transition_status"].startswith("imported continuation")
    assert law["signed_black_hole_entropy_change"].endswith("= -ln(k)")
    assert law["positive_entropy_loss"].endswith("= ln(k)")
    assert "finite display/submodel" in law["computed_k_scope"]
    assert receipt["schema"].endswith(".v3")
    assert law["tooth_offset"] == "Delta_f_k = c^3*g(chi)*ln(k)/(16*pi^2*G*M)"
    consts = receipt["constants_exact"]
    assert consts["c_m_per_s"] == "299792458"
    assert consts["gm_sun_nominal_m3_per_s2"] == "1.3271244E20"
    ladder = receipt["universal_ladder"]
    assert [row["k"] for row in ladder] == list(range(2, 13))
    for row in ladder:
        assert row["x_exact"] == "ln(%d)/(8*pi)" % row["k"]
        assert row["ratio_to_k2_exact"] == "ln(%d)/ln(2)" % row["k"]
        assert row["kms_weight_exact"] == "%d/%d" % (row["k"] - 1, row["k"])
    # 17-significant-digit guard on the k = 2 universal position: any
    # mutation of the frozen ln(k)/(8*pi) normalization changes this.
    assert ladder[0]["x_sig40"].startswith("2.757945001908144")
    assert ladder[0]["x_sig40"].endswith("E-2")
    # The k = 4 ratio renders as exactly two at forty digits.
    k4 = next(row for row in ladder if row["k"] == 4)
    assert k4["ratio_to_k2_sig40"] == "2." + "0" * 39 + "E+0"
    ref = receipt["reference_point_synthetic"]
    assert ref["mass_nominal_solar"] == "62"
    assert ref["chi"] == "0.67"
    assert ref["m_azimuthal"] == 2
    assert "synthetic" in ref["label"]
    assert receipt["frame_contract"]["detector_mass"] == (
        "M_det = (1+z)*M_source"
    )
    assert "posterior samples alone" in receipt["comparison_contract_boundary"]


def test_receipt_determinism_and_hygiene() -> None:
    with open(RECEIPT, "rb") as fh:
        on_disk = fh.read()
    rebuilt = producer.canonical_bytes(producer.build_receipt())
    assert on_disk == rebuilt
    text = on_disk.decode("ascii")
    assert "/Users/" not in text
    assert "/tmp/" not in text


def test_independent_verifier_agreement() -> None:
    proc = subprocess.run(
        [sys.executable, os.path.join(HERE, "verify_integer_k_comb_independent.py")],
        capture_output=True,
        text=True,
        cwd=HERE,
        timeout=120,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "VERIFIED" in proc.stdout
