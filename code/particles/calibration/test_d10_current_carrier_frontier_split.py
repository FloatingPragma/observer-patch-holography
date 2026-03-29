#!/usr/bin/env python3
"""Guard the split between builder-local and honest D10 repair frontiers."""

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
READOUT_OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_readout.json"
REPAIR_OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_repair_branch_beyond_current_carrier.json"


def test_d10_current_carrier_frontier_split_is_explicit() -> None:
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

    readout = json.loads(READOUT_OUTPUT.read_text(encoding="utf-8"))
    repair = json.loads(REPAIR_OUTPUT.read_text(encoding="utf-8"))

    assert readout["active_builder_smallest_missing_object"] == "EWTargetFreeRepairValueLaw_D10"
    assert readout["current_carrier_builder_local_frontier"] == "EWExactMassPairSelector_D10"
    assert readout["smallest_predictive_missing_object"] == "EWTargetFreeRepairValueLaw_D10"
    assert readout["exact_pdg_wz_frontier"] == "EWTargetFreeRepairValueLaw_D10"
    assert repair["object_id"] == "D10RepairBranchBeyondCurrentCarrier"
    assert repair["replaces_builder_local_frontier"] == "EWExactMassPairSelector_D10"
    assert repair["required_closure_kind"] == "single_family_single_P_no_mixed_readout"
    assert repair["operative_primitive"]["object_id"] == "EWSinglePostTransportTreeIdentity_D10"
