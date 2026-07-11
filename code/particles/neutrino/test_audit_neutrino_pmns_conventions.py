#!/usr/bin/env python3
"""Regression checks for the exhaustive PMNS convention audit."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "neutrino" / "audit_neutrino_pmns_conventions.py"


def test_convention_audit_finds_no_normal_ordering_rescue() -> None:
    tables = pathlib.Path("/tmp/oph-nufit61")
    tb_off = tables / "v61.release-TBoff-NO.txt.xz"
    tb_yes = tables / "v61.release-TByes-NO.txt.xz"
    if not tb_off.exists() or not tb_yes.exists():
        return
    with tempfile.TemporaryDirectory(prefix="oph_pmns_convention_audit_") as tmpdir:
        output = pathlib.Path(tmpdir) / "audit.json"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--tb-off-no",
                str(tb_off),
                "--tb-yes-no",
                str(tb_yes),
                "--output",
                str(output),
            ],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(output.read_text(encoding="utf-8"))
        assert payload["enumeration"]["total_assignments"] == 72
        assert payload["decision"]["any_normal_ordering_consistent_assignment_passes_both_3sigma_gates"] is False
        assert payload["decision"]["declared_assignment_is_best_by_worst_dataset_score"] is True
        assert payload["decision"]["convention_error_found"] is False
        declared = payload["declared_assignment"]
        assert declared["scores"]["TByes-NO"]["T23/DCP_delta_chi2"] > 20.11
