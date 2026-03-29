#!/usr/bin/env python3
"""Validate the D12 quark overlap transport law artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
D12_BRANCH_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_mass_branch_and_ckm_residual.py"
LAW_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_overlap_transport_law.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_overlap_transport_law.json"


def test_quark_d12_overlap_transport_law_collapses_tau_pair_to_one_scalar() -> None:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(D12_BRANCH_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(LAW_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_quark_d12_overlap_transport_law"
    assert payload["next_single_residual_object"] == "Delta_ud_overlap"
    assert payload["strictly_smaller_than"] == "source_readback_u_log_per_side_and_source_readback_d_log_per_side"
    candidate = payload["candidate_branch_from_t1_over_5"]
    best = payload["best_same_family_mass_point"]
    for branch in (candidate, best):
        assert abs(branch["tau_sum_half_delta_identity"]) < 1.0e-15
        assert abs(branch["tau_ratio_minus_sigma_ratio"]) < 1.0e-15
        assert abs(branch["lambda_minus_sigma_u_tau_u"]) < 1.0e-15
        assert abs(branch["lambda_minus_sigma_d_tau_d"]) < 1.0e-15
