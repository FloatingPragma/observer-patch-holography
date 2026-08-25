"""Tests for the joint rotation-curve and BTFR profile likelihood.

Covers a synthetic fixture with known a0 and known noise (recovery inside
the stated intervals plus interval-coverage sanity), mutation guards for a
flipped baryonic subtraction sign, a shifted deep cut, a wrong mass-to-light
convention, and a broken covariance block, and byte verification of the
committed receipt through the producer replay and the independent verifier.
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

import joint_rar_likelihood as mod  # noqa: E402

A0_TRUE = 1.0e-10
MASSES_MSUN = (1.0e9, 2.0e9, 3.0e9, 5.0e9, 1.0e10)

TEST_CONFIG = mod.Config(
    log10_a0_min=-10.25,
    log10_a0_max=-9.75,
    n_a0=201,
    n_upsilon=5,
    n_distance=3,
    n_inclination=3,
    rho_grid=(0.0, 0.3),
    fractions=(0.3,),
    bootstrap_replicates=200,
    seed=20260825,
)

CHANNEL_A_ONLY = mod.Config(
    log10_a0_min=-10.25,
    log10_a0_max=-9.75,
    n_a0=201,
    n_upsilon=5,
    n_distance=3,
    n_inclination=3,
    rho_grid=(0.0,),
    fractions=(),
    bootstrap_replicates=200,
    seed=20260825,
)


@pytest.fixture(scope="module")
def result() -> dict:
    return mod.run()


def synthetic_tables(
    masses_msun=MASSES_MSUN,
    a0: float = A0_TRUE,
    noise_kms: float = 0.0,
    e_vobs_kms: float = 4.0,
    rng: np.random.Generator | None = None,
    r_min_kpc: float = 2.0,
    r_max_kpc: float = 60.0,
    n_radii: int = 15,
    disk_encoding_upsilon: float = mod.UPSILON_DISK_CENTER,
) -> tuple[list[dict], list[dict]]:
    """Point-mass realizations of the deep laws with declared Gaussian noise.

    The baryonic squared speed v_bar^2 = G M / r is encoded through the disk
    component so that the fitted mass-to-light value
    ``disk_encoding_upsilon`` reproduces it exactly, and
    v_obs^2 = v_bar^2 + sqrt(v_bar^2 a0 r) exactly before noise.  Priors are
    centered at the truth: catalogue distance, inclination, and the
    committed disk mass-to-light 0.5 are the generating values whenever
    ``disk_encoding_upsilon`` is 0.5.
    """
    t1, t2 = [], []
    for k, m_sun in enumerate(masses_msun):
        name = f"SYN{k:02d}"
        mass = m_sun * mod.M_SUN
        t1.append(
            {
                "name": name,
                "dist_mpc": 10.0,
                "e_dist_mpc": 0.3,
                "inclination_deg": 60.0,
                "e_inclination_deg": 2.0,
                "L36_gsun": (mass / mod.M_SUN) / mod.UPSILON_DISK_CENTER / 1e9,
                "MHI_gsun": 0.0,
                "quality": 1,
            }
        )
        for r_kpc in np.geomspace(r_min_kpc, r_max_kpc, n_radii):
            r_m = float(r_kpc) * mod.KPC_M
            vbar2 = mod.G_SI * mass / r_m
            v_a2 = math.sqrt(vbar2 * a0 * r_m)
            vobs = math.sqrt(vbar2 + v_a2)
            if rng is not None and noise_kms > 0.0:
                vobs = vobs + float(rng.normal(0.0, noise_kms * mod.KM))
            t2.append(
                {
                    "name": name,
                    "rad_kpc": float(r_kpc),
                    "vobs_kms": vobs / mod.KM,
                    "e_vobs_kms": e_vobs_kms,
                    "vgas_kms": 0.0,
                    "vdisk_kms": math.sqrt(vbar2 / disk_encoding_upsilon)
                    / mod.KM,
                    "vbul_kms": 0.0,
                }
            )
    return t1, t2


# ---------------------------------------------------------------------------
# Receipt byte verification.


def test_committed_receipt_matches_code(result: dict) -> None:
    committed_bytes = mod.RECEIPT_PATH.read_bytes()
    assert mod.canonical_json(result).encode("utf-8") == committed_bytes
    receipt = json.loads(committed_bytes)
    assert receipt["schema"] == mod.SCHEMA
    assert receipt["physical_claim"] is False
    assert receipt["source_derived_output"] is False
    assert receipt["seen_data_postdiction"] is True
    assert receipt["seed"] == mod.SEED


def test_independent_verifier_passes() -> None:
    proc = subprocess.run(
        [sys.executable, str(HERE / "verify_joint_likelihood_independent.py")],
        capture_output=True,
        text=True,
        cwd=HERE,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert proc.stdout.startswith("PASS")


def test_receipt_pins_the_snapshot_digest(result: dict) -> None:
    # The digest pins the snapshot bytes the numbers were computed from; it
    # does not authenticate provenance or custody of the snapshot itself.
    import hashlib

    data = HERE.parent / "data"
    expected = hashlib.sha256(
        (data / "table1.dat").read_bytes() + (data / "table2.dat").read_bytes()
    ).hexdigest()
    assert result["sample"]["data_sha256"] == expected


def test_ml_inside_intervals_and_flags(result: dict) -> None:
    for block in result["subset_results"].values():
        for row in block["per_rho"]:
            ml = row["log10_a0_ml"]
            lo68, hi68 = row["interval68_log10"]
            lo95, hi95 = row["interval95_log10"]
            assert lo68 <= ml <= hi68
            assert lo95 <= lo68 and hi68 <= hi95
            assert row["interval95_pinned_at_grid_boundary"] == [False, False]
            if "paired_btfr" in row:
                p = row["paired_btfr"]
                n = p["bootstrap_replicates"]
                for side in ("le_zero", "ge_zero"):
                    count = p[f"count_log10_ratio_{side}"]
                    frac = p[f"{side}_plus_one_fraction"]
                    assert frac == pytest.approx((count + 1) / (n + 1))
                    assert frac > 0.0


# ---------------------------------------------------------------------------
# Synthetic recovery and coverage.


def test_noiseless_synthetic_recovers_a0_in_both_channels() -> None:
    t1, t2 = synthetic_tables()
    out = mod.run(t1, t2, TEST_CONFIG, include_provenance=False)
    for block in out["subset_results"].values():
        for row in block["per_rho"]:
            assert abs(row["log10_a0_ml"] - math.log10(A0_TRUE)) <= 0.005
    deep = out["subset_results"]["deep_f_0p3"]["per_rho"][0]
    paired = deep["paired_btfr"]
    assert paired["n_used"] == len(MASSES_MSUN)
    assert paired["a0_b_unweighted_log_mean_m_s2"] == pytest.approx(
        A0_TRUE, rel=1e-6
    )
    assert abs(paired["point_log10_ratio_b_over_a"]) <= 0.005
    assert paired["verdict"] == "CONSISTENT_INTERVAL_CONTAINS_ZERO"


def test_noisy_synthetic_interval_contains_truth() -> None:
    rng = np.random.default_rng(20260825)
    t1, t2 = synthetic_tables(noise_kms=4.0, rng=rng)
    out = mod.run(t1, t2, CHANNEL_A_ONLY, include_provenance=False)
    row = out["subset_results"]["full"]["per_rho"][0]
    lo95, hi95 = row["interval95_log10"]
    assert lo95 <= math.log10(A0_TRUE) <= hi95
    # Matched declared noise: the reduced chi-square sits near one, unlike
    # the committed-snapshot rows where the declared family underestimates
    # the scatter.
    assert 0.4 <= row["reduced_chi2_dof_no_nuisance"] <= 1.5


def test_interval_coverage_sanity() -> None:
    """68 percent interval coverage on repeated known-noise realizations.

    25 replicates with matched declared noise; the binomial band
    [10, 24] around 0.683 * 25 = 17 catches gross undercoverage (intervals
    far too narrow) and gross overcoverage (intervals far too wide).
    """
    hits = 0
    replicates = 25
    for k in range(replicates):
        rng = np.random.default_rng(1000 + k)
        t1, t2 = synthetic_tables(noise_kms=4.0, rng=rng)
        out = mod.run(t1, t2, CHANNEL_A_ONLY, include_provenance=False)
        row = out["subset_results"]["full"]["per_rho"][0]
        lo68, hi68 = row["interval68_log10"]
        if lo68 <= math.log10(A0_TRUE) <= hi68:
            hits += 1
    assert 10 <= hits <= 24, f"68 percent interval covered truth {hits}/25"


# ---------------------------------------------------------------------------
# Mutation guards.


def test_flipped_baryonic_subtraction_sign_fails_loudly() -> None:
    """The channel-B observable subtracts v_bar^2; addition is detected.

    On noiseless synthetic galaxies the subtracted estimator returns the
    true a0 to relative 1e-6 while the sign-flipped variant is biased above
    it by the factor (1 + 2 v_bar^2/v_A^2 + ...) for every galaxy, at least
    1.03 on this fixture, so a flipped subtraction cannot pass the recovery
    tolerance.
    """
    t1, t2 = synthetic_tables()
    galaxies, _ = mod.build_galaxies(t1, t2, TEST_CONFIG)
    assert len(galaxies) == len(MASSES_MSUN)
    for gal in galaxies:
        k = gal["idx_out"]
        vobs = gal["vobs_kms"][k] * mod.KM
        vbar2 = (
            gal["gas_term"][k]
            + mod.UPSILON_DISK_CENTER * gal["disk_term"][k]
            + gal["bul_term"][k]
        )
        mass = (
            mod.UPSILON_DISK_CENTER * gal["meta"]["L36_gsun"] * 1e9 * mod.M_SUN
        )
        correct = (vobs * vobs - vbar2) ** 2 / (mod.G_SI * mass)
        flipped = (vobs * vobs + vbar2) ** 2 / (mod.G_SI * mass)
        assert correct == pytest.approx(A0_TRUE, rel=1e-6)
        assert flipped > 1.03 * A0_TRUE
        assert flipped != pytest.approx(A0_TRUE, rel=1e-6)


def test_channel_b_estimate_hand_value() -> None:
    # One outermost point with round numbers: v_obs = 100 km/s, disk-only
    # v_bar^2 encoded for mass-to-light 0.5, at the central nuisance combo.
    # Hand value: v_A^2 = 1e10 - 0.5 * (2e9 * 2) = 8e9 m^2/s^2 (disk term
    # vdisk|vdisk| = 4e9 at vdisk = 63.2456 km/s squared);
    # M_b = (0.5 * 4 + 0) * 1e9 M_sun; a0_B = (8e9)^2 / (G M_b).
    gal = {
        "name": "HAND",
        "meta": {
            "name": "HAND",
            "dist_mpc": 10.0,
            "e_dist_mpc": 0.5,
            "inclination_deg": 60.0,
            "e_inclination_deg": 2.0,
            "L36_gsun": 4.0,
            "MHI_gsun": 0.0,
            "quality": 1,
        },
        "rad_kpc": np.array([10.0]),
        "vobs_kms": np.array([100.0]),
        "e_vobs_kms": np.array([2.0]),
        "gas_term": np.array([0.0]),
        "disk_term": np.array([4.0e9]),
        "bul_term": np.array([0.0]),
        "g_bar_cat": np.array([2.0e9 / (10.0 * mod.KPC_M)]),
        "idx_out": 0,
    }
    grids = mod.nuisance_grids(gal["meta"], TEST_CONFIG)
    center_flat = int(
        np.ravel_multi_index(
            (TEST_CONFIG.n_upsilon // 2, 1, 1),
            (TEST_CONFIG.n_upsilon, 3, 3),
        )
    )
    est = mod.channel_b_estimate(gal, grids, center_flat)
    assert est["upsilon"] == pytest.approx(0.5, rel=1e-12)
    assert est["distance_factor"] == pytest.approx(1.0, rel=1e-12)
    assert est["inclination_deg"] == pytest.approx(60.0, rel=1e-12)
    assert est["v_a2"] == pytest.approx(8.0e9, rel=1e-12)
    mass = 0.5 * 4.0 * 1e9 * mod.M_SUN
    assert est["mass"] == pytest.approx(mass, rel=1e-12)
    assert est["a0_b"] == pytest.approx(6.4e19 / (mod.G_SI * mass), rel=1e-12)


def test_shifted_deep_cut_is_detected(result: dict) -> None:
    """Deleting or moving the deep cut changes counted sample sizes.

    A galaxy with no point below the committed cut is absent from the deep
    subset and from the channel-B candidacy, and enters only when the cut is
    loosened; on the committed snapshot the cut is live (one f=0.3 candidate
    has a non-deep outermost point, and the deep subsets are proper subsets
    of the retained sample).
    """
    t1, t2 = synthetic_tables()
    t1_extra, t2_extra = synthetic_tables(
        masses_msun=(1.0e10,), r_min_kpc=0.5, r_max_kpc=2.5
    )
    t1_extra[0]["name"] = "SHALLOW"
    for row in t2_extra:
        row["name"] = "SHALLOW"
    out = mod.run(t1 + t1_extra, t2 + t2_extra, TEST_CONFIG, False)
    assert out["subset_results"]["full"]["n_galaxies"] == len(MASSES_MSUN) + 1
    deep = out["subset_results"]["deep_f_0p3"]
    assert deep["deep_cut_g_bar_max_m_s2"] == 0.3 * 1.2e-10
    assert deep["n_galaxies"] == len(MASSES_MSUN)
    paired = deep["per_rho"][0]["paired_btfr"]
    assert paired["n_candidate_galaxies"] == len(MASSES_MSUN)
    assert paired["n_used"] == len(MASSES_MSUN)
    loose = mod.Config(
        **{
            **{f: getattr(TEST_CONFIG, f) for f in TEST_CONFIG.__dataclass_fields__},
            "fractions": (5.0,),
        }
    )
    out_loose = mod.run(t1 + t1_extra, t2 + t2_extra, loose, False)
    paired_loose = out_loose["subset_results"]["deep_f_5p0"]["per_rho"][0][
        "paired_btfr"
    ]
    assert paired_loose["n_used"] == len(MASSES_MSUN) + 1
    # Committed-snapshot rows: the cut is live there as well.
    committed_0p3 = result["subset_results"]["deep_f_0p3"]
    assert committed_0p3["n_points"] < result["subset_results"]["full"]["n_points"]
    assert (
        committed_0p3["per_rho"][0]["paired_btfr"]["n_outermost_not_deep"] == 1
    )


def test_wrong_mass_to_light_shifts_the_estimate() -> None:
    """Data generated outside the declared mass-to-light convention misfit.

    The twin fixture encodes the same baryonic content at the committed 0.5
    and recovers a0 on the grid point of the truth; the mismatched fixture
    encodes it at 1.0, beyond the reach of the declared prior grid, and the
    recovered a0 is displaced, so a wrong mass-to-light convention cannot
    reproduce the committed receipt.
    """
    t1, t2 = synthetic_tables(disk_encoding_upsilon=0.5)
    out = mod.run(t1, t2, CHANNEL_A_ONLY, include_provenance=False)
    row = out["subset_results"]["full"]["per_rho"][0]
    assert abs(row["log10_a0_ml"] - math.log10(A0_TRUE)) <= 0.005

    t1_bad, t2_bad = synthetic_tables(disk_encoding_upsilon=1.0)
    out_bad = mod.run(t1_bad, t2_bad, CHANNEL_A_ONLY, include_provenance=False)
    row_bad = out_bad["subset_results"]["full"]["per_rho"][0]
    assert abs(row_bad["log10_a0_ml"] - math.log10(A0_TRUE)) > 0.04


def test_broken_covariance_block_fails_loudly() -> None:
    """Closed form against the explicit matrix inverse, plus validation."""
    rng = np.random.default_rng(7)
    n = 6
    u = rng.normal(size=n)
    for rho in (0.0, 0.2, 0.4, 0.6):
        R = (1.0 - rho) * np.eye(n) + rho * np.ones((n, n))
        explicit = float(u @ np.linalg.solve(R, u))
        assert mod.equal_correlation_chi2(u, rho) == pytest.approx(
            explicit, rel=1e-12
        )
        if rho > 0.0:
            # A block missing the 1/(1+(n-1)rho) correction is wrong.
            mutated = (np.sum(u * u) - rho * np.sum(u) ** 2) / (1.0 - rho)
            assert mutated != pytest.approx(explicit, rel=1e-6)
    for bad in (-0.1, 1.0, 1.5):
        with pytest.raises(ValueError):
            mod.equal_correlation_chi2(u, bad)
    with pytest.raises(ValueError):
        mod.profiled_curve(
            np.zeros((2, 1, 1, 1)),
            np.zeros((2, 1, 1, 1)),
            np.zeros(1),
            np.zeros((1, 1, 1)),
            3,
            1.0,
        )


def test_profiled_curve_matches_explicit_inverse_at_nonzero_rho() -> None:
    """The inline quadratic-form algebra of the production path, checked
    against the explicit per-point matrix inverse on nonzero residuals at
    nonzero rho (the guard for the production copy of the w algebra, in
    addition to the closed-form helper above)."""
    rng = np.random.default_rng(11)
    n = 5
    u = rng.normal(loc=0.7, size=n)
    s1 = float(np.sum(u * u))
    s2 = float(np.sum(u)) ** 2
    for rho in (0.2, 0.4, 0.6):
        R = (1.0 - rho) * np.eye(n) + rho * np.ones((n, n))
        explicit = float(u @ np.linalg.solve(R, u))
        out = mod.profiled_curve(
            np.full((1, 1, 1, 1), s1),
            np.full((1, 1, 1, 1), s2),
            np.zeros(1),
            np.zeros((1, 1, 1)),
            n,
            rho,
        )
        assert float(out["chi2_data"][0]) == pytest.approx(explicit, rel=1e-12)
        assert float(out["curve"][0]) == pytest.approx(explicit, rel=1e-12)
        assert mod.equal_correlation_chi2(u, rho) == pytest.approx(
            explicit, rel=1e-12
        )


def test_covariance_widens_or_shifts_reported_intervals(result: dict) -> None:
    # The rho scan is live on the committed snapshot: at least one deep
    # subset row changes its maximum-likelihood point or interval across
    # rho values, so silently dropping the scan is detectable.
    rows = result["subset_results"]["deep_f_0p1"]["per_rho"]
    assert len(rows) == len(result["error_model"]["rho_grid"])
    values = {json.dumps(r["interval95_log10"]) for r in rows}
    assert len(values) > 1


def test_verdict_rule_is_direction_neutral(result: dict) -> None:
    labels = {
        "TENSION_JOINT_CHANNEL_B_ABOVE_CHANNEL_A",
        "TENSION_JOINT_CHANNEL_B_BELOW_CHANNEL_A",
        "CONSISTENT_INTERVAL_CONTAINS_ZERO",
    }
    for label in labels:
        assert label in result["verdict_rule"]
    assert "not evidence for OPH" in result["verdict_rule"]
    assert mod._verdict(0.1, 0.2) == "TENSION_JOINT_CHANNEL_B_ABOVE_CHANNEL_A"
    assert mod._verdict(-0.2, -0.1) == "TENSION_JOINT_CHANNEL_B_BELOW_CHANNEL_A"
    assert mod._verdict(-0.1, 0.1) == "CONSISTENT_INTERVAL_CONTAINS_ZERO"
    assert mod._verdict(0.0, 0.1) == "CONSISTENT_INTERVAL_CONTAINS_ZERO"
    assert mod._verdict(-0.1, 0.0) == "CONSISTENT_INTERVAL_CONTAINS_ZERO"
    for f_key in ("deep_f_0p3", "deep_f_0p1"):
        for row in result["subset_results"][f_key]["per_rho"]:
            assert row["paired_btfr"]["verdict"] in labels
