"""Check calibrated receipt projection and outward reported bounds."""
from fractions import Fraction
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "native_postdiction_builder", ROOT / "code/particles/scripts/build_postdiction_ledger.py"
)
LEDGER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(LEDGER)


def test_receipt_projection_preserves_inputs_and_conditional_records():
    result = LEDGER._native_current_comparisons()
    assert len(result["source_sha256"]) == 10
    member = result["four_input_axial_member"]
    assert not member["independent_prediction"]
    for row in member["three_mass_fits"]:
        lower, upper = map(Fraction, row["same_source_energy_GeV_interval"])
        assert lower <= Fraction(row["target_GeV"]) <= upper
    neutral = result["calibrated_neutrinos"]
    assert not neutral["independent_prediction"]
    common = result["common_finite_source"]
    assert not common["physical_nucleon_pole_mass_claimed"]
    for row in common["energy_window_records"]:
        assert Fraction(row["original_inside_probability_lower"]) < Fraction(
            row["conditional_inside_probability_lower"]
        )
        assert Fraction(row["preparation_failure_probability"]) > 0
    rendered = "\n".join(LEDGER._render_native_current_comparisons(result))
    assert "The displayed agreements are calibrations" in rendered
    assert "momentum-sign parity contrast, distinct from group velocity" in rendered
    assert "zero lepton spatial kinetic coefficient and contains no photons" in rendered
    assert "absolute atmospheric phase at most 64*pi" in rendered


def test_decimal_bounds_are_outward():
    for rational in (str(Fraction(1) + Fraction(1, 10**200)), "1/3", "-1/3", "1137795909967161201100123216079496402779807388177/278759314981632789269196478408104518824755200000000000000"):
        for places in (0, 8, 10, 14):
            assert Fraction(LEDGER._directed_decimal(rational, places)) <= Fraction(rational)
            assert Fraction(LEDGER._directed_decimal(rational, places, upper=True)) >= Fraction(rational)


def test_stale_owner_source_aborts_projection(tmp_path, monkeypatch):
    target = tmp_path / "model/openwave_axial_state/results.json"
    target.parent.mkdir(parents=True)
    target.write_text(json.dumps({"source_sha256": {"source.py": "0" * 64}}))
    (tmp_path / "source.py").write_text("changed source")
    monkeypatch.setattr(LEDGER, "REPO", tmp_path)
    with pytest.raises(ValueError, match="native comparison source drift"):
        LEDGER._native_current_comparisons()


def test_mass_catalogue_projection_retains_calibration_and_failures():
    value = LEDGER._mass_and_decay_catalogue()
    assert len(value['Dirac_masses']) == 6
    assert len(value['complete_calibrated_comparison_catalogue']) == 17
    assert not value['source_joins']['all_sources_identified_as_one_interacting_theory']
    more = value['strong_masses_and_weak_widths']
    assert not more['tau_leptonic_width_ratio_comparison']['used_in_parameter_selection']
    assert more['nuclear_mass_calibration']['original_preparation_success'] == '1/17000'
    text = '\n'.join(LEDGER._render_mass_and_decay_catalogue(value))
    assert 'not total tau lifetimes' in text
    assert 'running quark parameters' in text


def test_hyperfine_projection_keeps_unused_failures_and_full_preparations():
    value = LEDGER._hyperfine_and_charge_profiles()
    assert len(value['five_mass_calibrations']) == 5
    assert value['all_unassigned_bound_nuclear_channels'] == []
    assert value['radius_comparison']['no_single_length_matches_both_radii_in_comparison_box']
    assert not value['radius_comparison']['physical_radius_target_used_to_fit_source_or_length']
    instrument = value['profile_instrument']
    assert Fraction(instrument['filter_success_exact']) + Fraction(instrument['filter_failure_exact']) == 1
    assert instrument['filter_success_exact'] == '1/120'
    text = '\n'.join(LEDGER._render_hyperfine_and_charge_profiles(value))
    assert 'seven empirical inputs' in text
    assert 'cannot both be fitted by any' in text
    assert 'not independent predictions or resonance poles' in text
    for row in value['profiles']:
        line = next(line for line in text.splitlines()
                    if line.startswith(f"| {row['assigned_channel']} | ["))
        endpoints = line.split('[', 1)[1].split(']', 1)[0].split(', ')
        lower, upper = map(Fraction, row['mean_square_charge_radius_fm2_interval'])
        assert Fraction(endpoints[0]) <= lower <= upper <= Fraction(endpoints[1])


def test_hyperfine_projection_rejects_changed_ancestor(tmp_path, monkeypatch):
    target = tmp_path / 'model/openwave_common_hyperfine_spectrum/results.json'
    target.parent.mkdir(parents=True)
    target.write_text(json.dumps({'verdict': 'PASS', 'source_sha256': {'source.py': '0' * 64}}))
    (tmp_path / 'source.py').write_text('different source')
    monkeypatch.setattr(LEDGER, 'REPO', tmp_path)
    with pytest.raises(ValueError, match='hyperfine source drift'):
        LEDGER._hyperfine_and_charge_profiles()
