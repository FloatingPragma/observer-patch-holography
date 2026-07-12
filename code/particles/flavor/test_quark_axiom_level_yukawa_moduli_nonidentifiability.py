#!/usr/bin/env python3
"""Guard the no-extra-axiom quark-Yukawa non-definability theorem."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = (
    ROOT
    / "particles"
    / "flavor"
    / "derive_quark_axiom_level_yukawa_moduli_nonidentifiability.py"
)
OUTPUT = (
    ROOT
    / "particles"
    / "runs"
    / "flavor"
    / "quark_axiom_level_yukawa_moduli_nonidentifiability.json"
)


def test_axioms_1_5_admit_physically_distinct_equal_mar_score_quark_spectra() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_axiom_level_yukawa_moduli_nonidentifiability"
    assert payload["proof_status"] == "closed_axiom_level_nondefinability_theorem"
    assert payload["additional_axioms_used"] is False
    assert payload["source_only_numeric_quark_spectrum_emitted"] is False
    assert payload["public_numeric_quark_rows_allowed"] is False

    family = payload["counterfamily"]
    assert family["parameter_space"] == "(lambda_u,lambda_d) in (R_{>0})^2"
    assert family["CKM_matrix_preserved"] is True
    assert family["CP_capability_preserved"] is True
    assert family["normalized_determinants_preserved"].endswith("= 1")
    assert family["quark_mass_singular_values_changed"] is True

    mar = payload["axiom_invariance_audit"]["Axiom_5_MAR"]
    assert mar["complexity_vector"] == "(chi_cpl,N_nonab,N_c,N_g)"
    assert mar["Yukawa_eigenvalues_are_components_of_complexity_vector"] is False
    assert mar["counterfamily_members_have_equal_MAR_score"] is True
    assert mar["counterfamily_members_remain_distinct_MAR_minima"] is True

    policy = payload["reference_data_policy"]
    assert policy["direct_input_artifacts"] == []
    assert policy["quark_reference_values_consumed"] is False
    assert policy["current_family_targets_consumed"] is False
    assert policy["fitted_spreads_consumed"] is False
    assert policy["numerical_flavor_template_consumed"] is False
    assert policy["no_target_leak_by_construction"] is True

    assert payload["corollaries"]["unique_source_map_P_to_six_quark_masses_exists"] is False
    assert payload["corollaries"]["MAR_breaks_counterfamily_under_current_definition"] is False
