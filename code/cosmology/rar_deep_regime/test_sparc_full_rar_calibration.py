"""Tests for the source-bound retrospective full-RAR calibration."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import sparc_full_rar_calibration as calibration  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_source_snapshot_and_cut_census() -> None:
    source = calibration.verify_source_snapshot()
    assert source["canonical_source_binding_receipt"] is True
    assert {row["path"] for row in source["files"]} == {
        "code/cosmology/rar_deep_regime/data/PROVENANCE.md",
        "code/cosmology/rar_deep_regime/data/ReadMe",
        "code/cosmology/rar_deep_regime/data/table1.dat",
        "code/cosmology/rar_deep_regime/data/table2.dat",
    }

    table1 = calibration.read_table1()
    table2 = calibration.read_table2()
    _, _, _, _, census = calibration.collect(
        table1,
        table2,
        disk_mass_to_light=0.5,
    )
    assert census == {
        "catalogue_galaxies": 175,
        "rotation_curve_points": 3391,
        "galaxies_excluded_quality_q_equals_3": 12,
        "galaxies_excluded_inclination_lt_30": 12,
        "galaxies_excluded_by_both_catalogue_cuts": 2,
        "galaxies_passing_catalogue_cuts": 153,
        "galaxies_with_retained_points": 147,
        "galaxies_passing_catalogue_cuts_with_zero_retained_points": 6,
        "catalogue_galaxies_at_inclination_eq_30": 5,
        "selected_galaxies_at_inclination_eq_30": 4,
        "contributing_galaxies_at_inclination_eq_30": 3,
        "retained_points_at_inclination_eq_30": 15,
        "points_after_catalogue_galaxy_cuts": 3168,
        "points_rejected_at_strict_relative_error_eq_0_1": 4,
        "galaxies_lost_at_strict_relative_error_eq_0_1": 2,
        "retained_points": 2696,
    }


def test_selection_boundaries_are_strict_and_shared() -> None:
    assert calibration.passes_galaxy_selection(30.0, 1)
    assert not calibration.passes_galaxy_selection(29.999, 1)
    assert not calibration.passes_galaxy_selection(30.0, 3)

    assert calibration.passes_point_selection(1.0, 100.0, 9.999)
    assert not calibration.passes_point_selection(1.0, 100.0, 10.0)
    assert not calibration.passes_point_selection(0.0, 100.0, 1.0)
    assert not calibration.passes_point_selection(1.0, 0.0, 1.0)


def test_committed_receipt_is_current_unit_safe_and_nonpromoting() -> None:
    receipt = json.loads(calibration.RECEIPT_PATH.read_text(encoding="utf-8"))
    fresh = calibration.build_report()
    assert calibration.render_report(fresh) == calibration.RECEIPT_PATH.read_bytes()
    assert receipt == fresh
    assert receipt["producer"]["script_sha256"] == sha256(
        Path(calibration.__file__)
    )

    gates = receipt["receipts"]
    assert gates["canonical_public_source_binding_receipt"] is True
    assert gates["deterministic_calibration_reproduction_receipt"] is True
    assert gates["independent_held_out_test_receipt"] is False
    assert gates["oph_specific_empirical_evidence_receipt"] is False
    assert gates["oph_validation_or_falsification_receipt"] is False
    assert receipt["interpretation"]["oph_candidate_values_consumed_by_fit"] is False
    assert receipt["interpretation"]["source_derived_a0"] is False
    assert receipt["interpretation"]["verdict"] == "NOT_ELIGIBLE_AS_OPH_EVIDENCE"

    for fit in receipt["fits"].values():
        assert fit["a0_in_1e_minus_10_m_per_s2"] == pytest.approx(
            fit["a0_si_m_per_s2"] / 1.0e-10,
            rel=1.0e-15,
        )
        assert fit[
            "galaxy_bootstrap_std_in_1e_minus_10_m_per_s2"
        ] == pytest.approx(
            fit["galaxy_bootstrap_std_si_m_per_s2"] / 1.0e-10,
            rel=1.0e-15,
        )
    standard = receipt["fits"]["disk_mass_to_light_0.5"]
    assert standard["a0_in_1e_minus_10_m_per_s2"] == pytest.approx(
        1.1612569877336687
    )
    assert standard["galaxy_bootstrap_std_in_1e_minus_10_m_per_s2"] == (
        pytest.approx(0.08020947274163046)
    )

    sensitivity = receipt["weighting_sensitivity"]["rows"][
        "disk_mass_to_light_0.5"
    ]
    equal_galaxy = sensitivity["equal_parent_galaxy_total_weight"]
    velocity_only = sensitivity["inverse_velocity_fractional_variance_only"]
    assert -10.0 < equal_galaxy["relative_to_primary_point_weighted_percent"] < -9.0
    assert 15.0 < velocity_only["relative_to_primary_point_weighted_percent"] < 16.0
    assert equal_galaxy["eligible_as_inference_or_oph_evidence"] is False
    assert velocity_only["eligible_as_inference_or_oph_evidence"] is False


def test_synthetic_full_rar_recovers_input() -> None:
    rng = np.random.default_rng(7)
    expected_a0 = 1.3e-10
    gbar = 10.0 ** rng.uniform(-13.0, -9.0, 1000)
    gobs = calibration.rar(gbar, expected_a0)
    fitted_a0 = calibration.fit(gbar, gobs)[0]
    assert fitted_a0 == pytest.approx(expected_a0, rel=1.0e-5)


def test_source_mutation_and_output_override_fail_closed(tmp_path: Path) -> None:
    source = tmp_path / "table1.dat"
    source.write_bytes(b"authentic fixture")
    expected = hashlib.sha256(source.read_bytes()).hexdigest()
    source.write_bytes(source.read_bytes() + b" tampered")
    with pytest.raises(
        calibration.SourceIntegrityError,
        match="source_sha256_mismatch:table1.dat",
    ):
        calibration._verify_pinned_file(source, expected)

    with pytest.raises(
        calibration.SourceIntegrityError,
        match="canonical_receipt_output_path_mismatch",
    ):
        calibration.main(["--output", str(tmp_path / "substitute.json")])
