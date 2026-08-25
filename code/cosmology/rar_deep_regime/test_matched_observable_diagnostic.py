"""Tests for the matched-observable dark-sector diagnostic.

Covers the estimator-correctness control (both channels must recover a known
constant on synthetic matched data with paired ratio one) and mutation
guards for the wrong subtraction sign, a missing or moved deep cut, and
wrong mass-to-light values, plus receipt replay and the independent
verifier.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import matched_observable_diagnostic as mod  # noqa: E402
import rar_deep_regime_diagnostic as base  # noqa: E402

A0_TRUE = 1.0e-10
MASSES_MSUN = (1.0e9, 3.0e9, 1.0e10)


@pytest.fixture(scope="module")
def result() -> dict:
    return mod.run()


def synthetic_tables(
    masses_msun=MASSES_MSUN, a0: float = A0_TRUE, r_max_kpc: float = 60.0
) -> tuple[list[dict], list[dict]]:
    """Exact point-mass realizations of the deep laws with zero scatter.

    v_bar^2 = G M / r is encoded through the disk component at the committed
    mass-to-light value, and v_obs^2 = v_bar^2 + sqrt(G M a0) exactly, so
    both channels must recover a0 and their paired ratio must be one.
    """
    t1, t2 = [], []
    for k, m_sun in enumerate(masses_msun):
        name = f"SYN{k:02d}"
        mass = m_sun * base.M_SUN
        t1.append(
            {
                "name": name,
                "quality": 1,
                "inclination_deg": 60.0,
                "e_inclination_deg": 2.0,
                "L36_gsun": (mass / base.M_SUN) / base.UPSILON_DISK / 1e9,
                "MHI_gsun": 0.0,
                "vflat_kms": 0.0,
                "e_vflat_kms": 0.0,
            }
        )
        v_a2 = math.sqrt(base.G_SI * mass * a0)
        for r_kpc in np.geomspace(0.5, r_max_kpc, 25):
            r_m = float(r_kpc) * base.KPC_M
            vbar2 = base.G_SI * mass / r_m
            vobs = math.sqrt(vbar2 + v_a2)
            t2.append(
                {
                    "name": name,
                    "rad_kpc": float(r_kpc),
                    "vobs_kms": vobs / base.KM,
                    "e_vobs_kms": 0.05 * vobs / base.KM,
                    "vgas_kms": 0.0,
                    "vdisk_kms": math.sqrt(2.0 * vbar2) / base.KM,
                    "vbul_kms": 0.0,
                }
            )
    return t1, t2


def test_committed_receipt_matches_code(result: dict) -> None:
    receipt = json.loads(mod.RECEIPT_PATH.read_text(encoding="utf-8"))
    assert receipt == result
    assert receipt["schema"] == mod.SCHEMA
    assert receipt["physical_claim"] is False
    assert receipt["source_derived_output"] is False
    assert receipt["seen_data_postdiction"] is True


def test_channel_a_reuses_committed_corrected_fit(result: dict) -> None:
    committed = json.loads(
        (HERE / "receipts" / "sparc_deep_regime_diagnostic.json").read_text(
            encoding="utf-8"
        )
    )
    primary = result["primary_fraction_0p1"]["channel_a_rar"]
    assert primary["a0_m_s2"] == committed["deep_radial_acceleration"][
        "a0_fixed_exponent_m_s2"
    ]
    assert primary["n_points"] == committed["deep_radial_acceleration"]["n_points"]
    sweep_0p3 = committed["deep_radial_acceleration_fraction_sweep"][0]
    assert sweep_0p3["deep_fraction_of_reference_a0"] == 0.3
    assert result["sensitivity_fraction_0p3"]["channel_a_rar"]["a0_m_s2"] == (
        sweep_0p3["a0_fixed_exponent_m_s2"]
    )


def test_anomalous_subtraction_sign_and_mass_to_light() -> None:
    # Hand value: 1e10 - (30e3^2 + 0.5*40e3^2 + 0.7*10e3^2) = 8.23e9 m^2/s^2.
    # A flipped subtraction sign, an added baryonic term, or any other
    # mass-to-light value changes this number.
    v_a2, vbar2 = mod.anomalous_squared_speed_m2s2(100.0, 30.0, 40.0, 10.0)
    assert vbar2 == pytest.approx(1.77e9, rel=1e-12)
    assert v_a2 == pytest.approx(8.23e9, rel=1e-12)
    # Signed gas convention: a negative gas component subtracts.
    v_a2_neg, vbar2_neg = mod.anomalous_squared_speed_m2s2(100.0, -30.0, 40.0, 10.0)
    assert vbar2_neg == pytest.approx(-3.0e7, rel=1e-12)
    assert v_a2_neg == pytest.approx(1.003e10, rel=1e-12)


def test_deep_cut_fixed_absolute_and_shared(result: dict) -> None:
    primary = result["primary_fraction_0p1"]
    assert primary["deep_cut_g_bar_max_m_s2"] == 0.1 * 1.2e-10
    assert primary["channel_b_matched"][
        "deep_outermost_requirement_g_bar_max_m_s2"
    ] == primary["deep_cut_g_bar_max_m_s2"]
    sens = result["sensitivity_fraction_0p3"]
    assert sens["deep_cut_g_bar_max_m_s2"] == 0.3 * 1.2e-10
    assert result["deep_cut"]["reference_a0_m_s2"] == 1.2e-10
    # The cut is live on the committed snapshot: it excludes galaxies.
    assert primary["channel_b_matched"]["n_outermost_not_deep"] > 0
    assert (
        primary["channel_b_matched"]["n_used"]
        < primary["channel_b_matched"]["n_candidate_galaxies"]
    )


def test_synthetic_both_channels_recover_true_a0() -> None:
    t1, t2 = synthetic_tables()
    fraction = 0.3
    b_res = mod.channel_b(t1, t2, fraction)
    summary = b_res["summary"]
    assert summary["n_used"] == len(MASSES_MSUN)
    assert summary["n_outermost_not_deep"] == 0
    assert summary["a0_unweighted_log_mean_m_s2"] == pytest.approx(A0_TRUE, rel=1e-9)
    assert summary["a0_same_radius_mass_variant_m_s2"] == pytest.approx(
        A0_TRUE, rel=1e-9
    )
    acc = base.accelerations(t1, t2)
    block = mod.fraction_block(t1, t2, acc, fraction)
    assert block["channel_a_rar"]["a0_m_s2"] == pytest.approx(A0_TRUE, rel=1e-6)
    paired = block["paired_common_set"]
    assert paired["n_common_galaxies"] == len(MASSES_MSUN)
    assert abs(paired["point_log10_ratio_b_over_a"]) < 1e-6
    lo, hi = paired["log10_ratio_95pct"]
    assert abs(lo) < 1e-6 and abs(hi) < 1e-6


def test_unsubtracted_proxy_is_biased_upward() -> None:
    # The audited proxy defect: the finite-radius total speed without
    # baryonic subtraction overestimates a0 in every synthetic galaxy.
    t1, t2 = synthetic_tables()
    rows = mod.outermost_retained_points(t1, t2)
    meta = {r["name"]: r for r in t1}
    for name, r in rows.items():
        vo = r["vobs_kms"] * base.KM
        mass = (
            base.UPSILON_DISK * meta[name]["L36_gsun"] * 1e9
            + base.GAS_HELIUM * meta[name]["MHI_gsun"] * 1e9
        ) * base.M_SUN
        proxy = (vo * vo) ** 2 / (base.G_SI * mass)
        matched = r["v_a2"] * r["v_a2"] / (base.G_SI * mass)
        assert proxy > matched
        assert matched == pytest.approx(A0_TRUE, rel=1e-9)
        assert proxy > A0_TRUE


def test_missing_deep_cut_is_detected() -> None:
    # A galaxy whose outermost point is not deep must be excluded at the
    # committed cut and admitted only when the cut is moved, so deleting or
    # loosening the cut changes counted sample sizes.
    t1, t2 = synthetic_tables()
    t1_extra, t2_extra = synthetic_tables(
        masses_msun=(1.0e10,), r_max_kpc=3.0
    )
    t1_extra[0]["name"] = "SHALLOW"
    for row in t2_extra:
        row["name"] = "SHALLOW"
    b_res = mod.channel_b(t1 + t1_extra, t2 + t2_extra, 0.3)
    assert b_res["summary"]["n_outermost_not_deep"] == 1
    assert b_res["summary"]["n_used"] == len(MASSES_MSUN)
    assert "SHALLOW" not in b_res["names"].tolist()
    loose = mod.channel_b(t1 + t1_extra, t2 + t2_extra, 5.0)
    assert loose["summary"]["n_used"] == len(MASSES_MSUN) + 1
    assert "SHALLOW" in loose["names"].tolist()


def test_nonpositive_anomalous_speed_is_excluded_and_counted() -> None:
    t1 = [
        {
            "name": "NEG",
            "quality": 1,
            "inclination_deg": 60.0,
            "e_inclination_deg": 2.0,
            "L36_gsun": 10.0,
            "MHI_gsun": 1.0,
        }
    ]
    # Deep outermost point (g_bar ~ 4e-13) with v_bar above v_obs.
    t2 = [
        {
            "name": "NEG",
            "rad_kpc": 100.0,
            "vobs_kms": 40.0,
            "e_vobs_kms": 2.0,
            "vgas_kms": 0.0,
            "vdisk_kms": math.sqrt(2.0) * 50.0,
            "vbul_kms": 0.0,
        }
    ]
    b_res = mod.channel_b(t1, t2, 0.3)
    assert b_res["summary"]["n_nonpositive_anomalous_speed2"] == 1
    assert b_res["summary"]["n_used"] == 0
    assert b_res["summary"]["a0_unweighted_log_mean_m_s2"] is None
    assert "empty_selection" in b_res["summary"]


def test_zero_tallies_reported_with_plus_one_fractions(result: dict) -> None:
    for key in ("primary_fraction_0p1", "sensitivity_fraction_0p3"):
        paired = result[key]["paired_common_set"]
        n = paired["bootstrap_replicates"]
        for side in ("le_zero", "ge_zero"):
            count = paired[f"count_log10_ratio_{side}"]
            frac = paired[f"{side}_plus_one_fraction"]
            assert frac == pytest.approx((count + 1) / (n + 1))
            assert frac > 0.0


def test_verdict_rule_is_direction_neutral(result: dict) -> None:
    labels = {
        "TENSION_MATCHED_CHANNEL_B_ABOVE_CHANNEL_A",
        "TENSION_MATCHED_CHANNEL_B_BELOW_CHANNEL_A",
        "CONSISTENT_INTERVAL_CONTAINS_ZERO",
    }
    assert result["diagnostic_verdict"] in labels
    for label in labels:
        assert label in result["verdict_rule"]
    assert "not evidence for OPH" in result["verdict_rule"]
    assert mod._verdict(0.1, 0.2) == "TENSION_MATCHED_CHANNEL_B_ABOVE_CHANNEL_A"
    assert mod._verdict(-0.2, -0.1) == "TENSION_MATCHED_CHANNEL_B_BELOW_CHANNEL_A"
    assert mod._verdict(-0.1, 0.1) == "CONSISTENT_INTERVAL_CONTAINS_ZERO"
    assert mod._verdict(0.0, 0.1) == "CONSISTENT_INTERVAL_CONTAINS_ZERO"
    assert mod._verdict(-0.1, 0.0) == "CONSISTENT_INTERVAL_CONTAINS_ZERO"


def test_supersession_block_quantifies_both_diagnostics(result: dict) -> None:
    block = result["superseded_mixed_proxy_diagnostic"]
    committed = json.loads(
        (HERE / "receipts" / "sparc_deep_regime_diagnostic.json").read_text(
            encoding="utf-8"
        )
    )
    paired = committed["constant_consistency"]["paired_parent_galaxy_bootstrap"]
    assert block["mixed_proxy_log10_ratio_95pct"] == (
        paired["log10_ratio_btfr_over_rar_95pct"]
    )
    assert block["mixed_proxy_point_log10_ratio"] == (
        committed["constant_consistency"]["log10_ratio_btfr_over_deep"]
    )
    assert block["displacement_removed_dex"] == pytest.approx(
        block["mixed_proxy_point_log10_ratio"] - block["matched_point_log10_ratio"]
    )
    assert block["resolution_statement"]


def test_independent_verifier_passes() -> None:
    proc = subprocess.run(
        [sys.executable, str(HERE / "verify_matched_observable_independent.py")],
        capture_output=True,
        text=True,
        cwd=HERE,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert proc.stdout.startswith("PASS")
