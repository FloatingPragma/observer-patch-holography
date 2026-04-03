#!/usr/bin/env python3
"""Validate the machine-readable maximal theorem-emitted quark package boundary."""

from __future__ import annotations

import json
import pathlib


ROOT = pathlib.Path(__file__).resolve().parents[2]
MAXIMAL = ROOT / "particles" / "runs" / "flavor" / "quark_maximal_theorem_emitted_package.json"
NO_GO = ROOT / "particles" / "runs" / "flavor" / "quark_current_premise_no_go_theorem.json"
MINIMAL = ROOT / "particles" / "runs" / "flavor" / "quark_minimal_extension_closure_theorem.json"


def test_quark_maximal_theorem_emitted_boundary_is_explicit() -> None:
    maximal = json.loads(MAXIMAL.read_text(encoding="utf-8"))
    assert maximal["artifact"] == "oph_quark_maximal_theorem_emitted_package"
    assert maximal["status"] == "theorem_grade_boundary_emitted"
    assert maximal["emitted_package"]["same_label_left_handed_selector"]["sigma_id"] == "sigma_ref"
    assert maximal["emitted_package"]["shared_norm_binding"]["g_ch"] == 0.9231656602589082
    assert maximal["emitted_package"]["current_family_affine_mean_package"]["g_u"] == 0.7797392875757557
    assert maximal["emitted_package"]["current_family_affine_mean_package"]["g_d"] == 0.12172551081512113
    assert maximal["stronger_objects_not_emitted"] == [
        "Theta_ud^mass",
        "Theta_ud^phys",
        "Theta_ud^abs",
    ]


def test_quark_no_go_and_minimal_extension_reference_maximal_package() -> None:
    no_go = json.loads(NO_GO.read_text(encoding="utf-8"))
    minimal = json.loads(MINIMAL.read_text(encoding="utf-8"))
    assert no_go["maximal_theorem_emitted_package_artifact"] == "oph_quark_maximal_theorem_emitted_package"
    assert no_go["exact_minimal_extension_triple"] == {
        "H_mass": "ell_ud := log(c_d / c_u)",
        "H_phys": "s_ud^phys : M_ud^{CR,phys} -> Sigma_ud^phys",
        "H_abs": "A_q^phys : Sigma_ud^phys -> R",
    }
    assert minimal["base_theorem_emitted_package_artifact"] == "oph_quark_maximal_theorem_emitted_package"
    assert minimal["exact_minimal_extension_triple"] == no_go["exact_minimal_extension_triple"]
