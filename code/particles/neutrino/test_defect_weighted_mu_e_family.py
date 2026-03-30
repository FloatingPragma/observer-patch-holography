#!/usr/bin/env python3
"""Guard the defect-weighted mu_e family against readback drift."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "neutrino" / "derive_defect_weighted_mu_e_family.py"
OUTPUT = ROOT / "runs" / "neutrino" / "defect_weighted_mu_e_family.json"


def test_defect_weighted_mu_e_family_tracks_live_readback_status() -> None:
    subprocess.run(["python3", str(SCRIPT)], check=True)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_defect_weighted_majorana_edge_weight_family"
    assert payload["realized_same_label_gap_defect_readback_status"] == "complete_from_live_flavor_artifacts"
    assert payload["smallest_constructive_missing_object"] is None
    assert payload["strict_repo_missing_object"] is None
    assert payload["upstream_exact_clause"] == "same_label_overlap_nonzero_on_realized_refinement_arrows"
