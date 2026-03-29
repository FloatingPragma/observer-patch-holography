#!/usr/bin/env python3
"""Validate the D10 electroweak exactness audit artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
FAMILY_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_observable_family.py"
SOURCE_PAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_pair.py"
READOUT_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_readout.py"
POPULATION_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_population_evaluator.py"
EXACT_CLOSURE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_exact_closure_beyond_current_carrier.py"
FIBERWISE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.py"
OBSTRUCTION_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_tau2_current_carrier_obstruction.py"
EXACT_WZ_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_exact_wz_coordinate_beyond_single_tree_identity.py"
EXACT_CHART_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_exact_mass_pair_chart_current_carrier.py"
REPAIR_BRANCH_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_repair_branch_beyond_current_carrier.py"
FACTORIZATION_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_w_anchor_neutral_shear_factorization.py"
AUDIT_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_exactness_audit.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exactness_audit.json"


def test_d10_exactness_audit_records_mass_ratio_identity_obstruction() -> None:
    subprocess.run([sys.executable, str(FAMILY_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_PAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(READOUT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(POPULATION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EXACT_CLOSURE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FIBERWISE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(OBSTRUCTION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EXACT_WZ_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EXACT_CHART_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(REPAIR_BRANCH_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FACTORIZATION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(READOUT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(AUDIT_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d10_ew_exactness_audit"
    identity = payload["tree_level_identity_audit"]
    fixed_eta = payload["fixed_eta_single_sigma_audit"]
    assert abs(identity["mass_ratio_sin2_from_reference_WZ"] - payload["reference_wz_audit_slice"]["sin2w_eff"]) < 1.0e-12
    assert abs(identity["identity_residual_on_mixed_reference_surface"]) > 1.0e-3
    assert fixed_eta["verdict"]["mass_pair_nearly_coherent"] is True
    assert fixed_eta["verdict"]["running_alpha_conflicts_with_mass_pair"] is True
    assert fixed_eta["verdict"]["running_sin2_forces_family_escape"] is True
    assert payload["exact_closure_beyond_current_carrier"]["status"] == "closed"
    assert payload["fiberwise_population_tree_law_beneath_single_tree_identity"]["status"] == "closed_smaller_primitive"
    assert payload["tau2_current_carrier_obstruction"]["status"] == "closed_smaller_primitive"
    assert payload["exact_mass_pair_chart_current_carrier"]["status"] == "closed_smaller_primitive"
    current_pair = payload["current_carrier_closure_summary"]["current_carrier_exact_mass_pair"]
    assert abs(current_pair["MW_pole"] - 80.38629169244275) < 1.0e-12
    assert abs(current_pair["MZ_pole"] - 91.18290444674243) < 1.0e-12
    assert payload["current_carrier_closure_summary"]["current_carrier_builder_local_frontier"] == "EWExactMassPairSelector_D10"
    freeze_once = payload["freeze_once_coherent_repair_summary"]
    assert freeze_once["status"] == "closed_freeze_once_coherent_repair_law"
    assert abs(freeze_once["public_mass_pair"]["MW_pole"] - 80.377) < 1.0e-12
    assert abs(freeze_once["public_mass_pair"]["MZ_pole"] - 91.18797809193725) < 1.0e-12
    assert payload["active_builder_smallest_missing_object"] == "EWTargetFreeRepairValueLaw_D10"
    assert payload["broader_honest_repair_frontier"] == "EWTargetFreeRepairValueLaw_D10"
    assert payload["exact_pdg_wz_frontier"] == "EWTargetFreeRepairValueLaw_D10"
    assert payload["smallest_constructive_missing_object"] == "EWTargetFreeRepairValueLaw_D10"
    assert payload["d10_repair_branch_beyond_current_carrier"]["object_id"] == "D10RepairBranchBeyondCurrentCarrier"
    assert payload["d10_repair_branch_beyond_current_carrier"]["replaces_builder_local_frontier"] == "EWExactMassPairSelector_D10"
    assert payload["d10_repair_branch_beyond_current_carrier"]["stronger_residual_object"] == "EWSinglePostTransportTreeIdentity_D10"
    assert (
        payload["target_free_source_only_underdetermination"]["unconditional_theorem"]["name"]
        == "current_corpus_underdetermination_of_forward_d10_repair_law"
    )
    assert payload["target_free_source_only_candidate"]["object_id"] == "EWTargetEmitter_D10"
