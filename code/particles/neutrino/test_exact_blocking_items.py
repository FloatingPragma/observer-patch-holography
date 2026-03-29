#!/usr/bin/env python3
"""Validate the exact neutrino blocker audit."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "neutrino" / "derive_exact_blocking_items.py"


def test_exact_blocking_items_reports_isotropy_and_live_missing_objects() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_neutrino_blockers_") as tmpdir:
        tmp = pathlib.Path(tmpdir)
        forward = tmp / "forward.json"
        certificate = tmp / "certificate.json"
        pmns = tmp / "pmns.json"
        charged_left = tmp / "charged_left.json"
        eta_demo = tmp / "eta_demo.json"
        intrinsic = tmp / "intrinsic.json"
        repair = tmp / "repair.json"
        exact_out = tmp / "exact.json"
        summary_out = tmp / "summary.json"

        forward.write_text(
            json.dumps(
                {
                    "masses_gev_sorted": [
                        2.3986448447627196e-12,
                        2.3986448447627196e-12,
                        2.590074050773907e-12,
                    ],
                    "delta_m21_sq_gev2": 0.0,
                    "delta_m31_sq_gev2": 9.549864971855843e-25,
                    "ordering_phase_certified": "normal_like_collective_dominance",
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        certificate.write_text(
            json.dumps({"artifact": "oph_neutrino_same_label_scalar_certificate", "sufficient_for_intrinsic_mass_eigenstates": False}, indent=2)
            + "\n",
            encoding="utf-8",
        )
        pmns.write_text(json.dumps({"status": "open"}, indent=2) + "\n", encoding="utf-8")
        charged_left.write_text(json.dumps({"status": "open"}, indent=2) + "\n", encoding="utf-8")
        eta_demo.write_text(
            json.dumps({"eta_e": {"psi12": 0.1, "psi23": -0.2, "psi31": 0.1}}, indent=2) + "\n",
            encoding="utf-8",
        )
        intrinsic.write_text(
            json.dumps(
                {
                    "collective_vector_actual_aligned": True,
                    "solar_split_actual_gev2": 5.690121167370743e-26,
                    "delta_m31_actual_gev2": 9.798023388600230e-25,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        repair.write_text(json.dumps({}, indent=2) + "\n", encoding="utf-8")

        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--forward",
                str(forward),
                "--certificate",
                str(certificate),
                "--pmns",
                str(pmns),
                "--charged-left",
                str(charged_left),
                "--eta-demo",
                str(eta_demo),
                "--intrinsic-validation",
                str(intrinsic),
                "--repair",
                str(repair),
                "--exact-output",
                str(exact_out),
                "--summary-output",
                str(summary_out),
            ],
            check=True,
            cwd=ROOT,
        )

        exact_payload = json.loads(exact_out.read_text(encoding="utf-8"))
        summary_payload = json.loads(summary_out.read_text(encoding="utf-8"))
        assert exact_payload["artifact"] == "oph_exact_neutrino_blocker_audit_v8"
        assert exact_payload["neutrino_only_isotropy_obstruction"]["closed"] is True
        assert exact_payload["exact_blocker_counts"]["same_label_proof_facing_continuous_dof_mod_common_scale"] == 5
        assert exact_payload["exact_blocker_counts"]["same_label_builder_facing_centered_eta_dof"] == 2
        assert exact_payload["exact_blocker_counts"]["charged_left_basis_artifact_dof_before_phase_quotients"] == 9
        assert [item["name"] for item in exact_payload["exact_blockers"]] == [
            "live_same_label_scalar_certificate",
            "shared_charged_lepton_left_basis",
        ]
        assert summary_payload["exact_remaining_blockers"] == [
            "live_same_label_scalar_certificate",
            "shared_charged_lepton_left_basis",
        ]


def test_exact_blocking_items_close_when_certificate_basis_and_pmns_are_live() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_neutrino_blockers_closed_") as tmpdir:
        tmp = pathlib.Path(tmpdir)
        forward = tmp / "forward.json"
        certificate = tmp / "certificate.json"
        pmns = tmp / "pmns.json"
        charged_left = tmp / "charged_left.json"
        eta_demo = tmp / "eta_demo.json"
        intrinsic = tmp / "intrinsic.json"
        repair = tmp / "repair.json"
        exact_out = tmp / "exact.json"
        summary_out = tmp / "summary.json"

        forward.write_text(
            json.dumps(
                {
                    "masses_gev_sorted": [2.38e-12, 2.42e-12, 2.58e-12],
                    "delta_m21_sq_gev2": 1.7e-25,
                    "delta_m31_sq_gev2": 1.0e-24,
                    "ordering_phase_certified": "normal_like_collective_dominance",
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        certificate.write_text(
            json.dumps({"artifact": "oph_neutrino_same_label_scalar_certificate", "sufficient_for_intrinsic_mass_eigenstates": True}, indent=2)
            + "\n",
            encoding="utf-8",
        )
        pmns.write_text(json.dumps({"status": "closed"}, indent=2) + "\n", encoding="utf-8")
        charged_left.write_text(json.dumps({"status": "closed"}, indent=2) + "\n", encoding="utf-8")
        eta_demo.write_text(
            json.dumps({"eta_e": {"psi12": 0.1, "psi23": -0.2, "psi31": 0.1}}, indent=2) + "\n",
            encoding="utf-8",
        )
        intrinsic.write_text(
            json.dumps(
                {
                    "collective_vector_actual_aligned": True,
                    "solar_split_actual_gev2": 5.690121167370743e-26,
                    "delta_m31_actual_gev2": 9.798023388600230e-25,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        repair.write_text(json.dumps({}, indent=2) + "\n", encoding="utf-8")

        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--forward",
                str(forward),
                "--certificate",
                str(certificate),
                "--pmns",
                str(pmns),
                "--charged-left",
                str(charged_left),
                "--eta-demo",
                str(eta_demo),
                "--intrinsic-validation",
                str(intrinsic),
                "--repair",
                str(repair),
                "--exact-output",
                str(exact_out),
                "--summary-output",
                str(summary_out),
            ],
            check=True,
            cwd=ROOT,
        )

        exact_payload = json.loads(exact_out.read_text(encoding="utf-8"))
        summary_payload = json.loads(summary_out.read_text(encoding="utf-8"))
        assert exact_payload["artifact"] == "oph_exact_neutrino_blocker_audit_v8"
        assert exact_payload["fully_completed"] is False
        assert exact_payload["live_continuation_branch_status"]["status"] == "numerically_closed_but_quantitatively_wrong_branch"
        assert [item["name"] for item in exact_payload["exact_blockers"]] == [
            "physical_neutrino_branch_repair"
        ]
        assert summary_payload["exact_remaining_blockers"] == [
            "physical_neutrino_branch_repair"
        ]


def test_exact_blocking_items_reduce_to_one_absolute_normalization_after_repair() -> None:
    with tempfile.TemporaryDirectory(prefix="oph_neutrino_blockers_repaired_") as tmpdir:
        tmp = pathlib.Path(tmpdir)
        forward = tmp / "forward.json"
        certificate = tmp / "certificate.json"
        pmns = tmp / "pmns.json"
        charged_left = tmp / "charged_left.json"
        eta_demo = tmp / "eta_demo.json"
        intrinsic = tmp / "intrinsic.json"
        repair = tmp / "repair.json"
        exact_out = tmp / "exact.json"
        summary_out = tmp / "summary.json"

        forward.write_text(
            json.dumps(
                {
                    "masses_gev_sorted": [2.38e-12, 2.42e-12, 2.58e-12],
                    "delta_m21_sq_gev2": 1.7e-25,
                    "delta_m31_sq_gev2": 1.0e-24,
                    "ordering_phase_certified": "normal_like_collective_dominance",
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        certificate.write_text(
            json.dumps({"artifact": "oph_neutrino_same_label_scalar_certificate", "sufficient_for_intrinsic_mass_eigenstates": True}, indent=2)
            + "\n",
            encoding="utf-8",
        )
        pmns.write_text(json.dumps({"status": "closed"}, indent=2) + "\n", encoding="utf-8")
        charged_left.write_text(json.dumps({"status": "closed"}, indent=2) + "\n", encoding="utf-8")
        eta_demo.write_text(
            json.dumps({"eta_e": {"psi12": 0.1, "psi23": -0.2, "psi31": 0.1}}, indent=2) + "\n",
            encoding="utf-8",
        )
        intrinsic.write_text(
            json.dumps(
                {
                    "collective_vector_actual_aligned": True,
                    "solar_split_actual_gev2": 5.690121167370743e-26,
                    "delta_m31_actual_gev2": 9.798023388600230e-25,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        repair.write_text(
            json.dumps(
                {
                    "artifact": "oph_neutrino_weighted_cycle_repair",
                    "physical_window_status": "pmns_and_hierarchy_repaired",
                    "absolute_normalization_status": "open_one_positive_scale",
                    "pmns_observables": {"theta12_deg": 33.97},
                    "dimensionless_dm2": {"21": 1.0, "32": 2.0},
                    "compare_only_atmospheric_anchor": {"delta_m21_sq_eV2": 7.7e-5},
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--forward",
                str(forward),
                "--certificate",
                str(certificate),
                "--pmns",
                str(pmns),
                "--charged-left",
                str(charged_left),
                "--eta-demo",
                str(eta_demo),
                "--intrinsic-validation",
                str(intrinsic),
                "--repair",
                str(repair),
                "--exact-output",
                str(exact_out),
                "--summary-output",
                str(summary_out),
            ],
            check=True,
            cwd=ROOT,
        )

        exact_payload = json.loads(exact_out.read_text(encoding="utf-8"))
        summary_payload = json.loads(summary_out.read_text(encoding="utf-8"))
        assert exact_payload["artifact"] == "oph_exact_neutrino_blocker_audit_v8"
        assert exact_payload["live_continuation_branch_status"]["status"] == "physically_repaired_up_to_one_positive_scale"
        assert exact_payload["live_continuation_branch_status"]["absolute_scale_no_go"]["theorem"] == "neutrino_weighted_cycle_absolute_scale_no_go"
        assert exact_payload["live_continuation_branch_status"]["absolute_scale_no_go"]["proof_obstruction"] == "positive_rescaling_nonidentifiability"
        assert [item["name"] for item in exact_payload["exact_blockers"]] == [
            "one_positive_neutrino_mass_normalization_scalar"
        ]
        assert summary_payload["exact_remaining_blockers"] == [
            "one_positive_neutrino_mass_normalization_scalar"
        ]
