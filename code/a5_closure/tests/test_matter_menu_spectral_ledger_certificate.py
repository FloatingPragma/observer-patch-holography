#!/usr/bin/env python3
"""Regression and adversarial tests for the #612/#609 spectral-ledger certificate."""

from __future__ import annotations

import copy
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import matter_menu_spectral_ledger_certificate as cert  # noqa: E402
import super_tannakian_matter_lift_certificate as m314  # noqa: E402


class MatterMenuSpectralLedgerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_path = (
            MODULE_DIR / "manifests" / "matter_menu_spectral_ledger_reference.json"
        )
        cls.stored = cert.load_json(cls.manifest_path)
        cls.payload = cert.certificate_payload()
        cls.pinned = cert.load_pinned_matter_lift()

    def test_manifest_is_deterministic_and_matches_stored(self) -> None:
        recomputed = cert.certificate_payload()
        self.assertEqual(self.payload, recomputed)
        self.assertEqual(self.stored, self.payload)
        cert.verify_manifest(self.stored)
        self.assertEqual(
            cert.sha256_json(self.stored), cert.sha256_json(self.payload)
        )

    def test_schema_and_issues(self) -> None:
        self.assertEqual(
            self.payload["schema"], "oph.matter_menu_spectral_ledger_certificate.v1"
        )
        self.assertEqual(self.payload["issues"], [612, 609])

    def test_pinned_matter_lift_values(self) -> None:
        pinned_block = self.payload["pinned_matter_lift"]
        manifest = cert.load_json(MODULE_DIR / "manifests" / "super_tannakian_matter_reference.json")
        receipt = cert.load_json(MODULE_DIR / "receipts" / "super_tannakian_matter_reference.receipt.json")
        self.assertEqual(pinned_block["manifest_sha256"], cert.sha256_json(manifest))
        self.assertEqual(pinned_block["receipt_sha256"], cert.sha256_json(receipt))
        self.assertEqual(receipt["selection"]["projector_rank"], 15)
        self.assertEqual(self.pinned["a"], -2)
        self.assertEqual(self.pinned["b"], 3)
        algebra = self.payload["declared_algebra"]
        self.assertEqual(
            algebra["derived_block_charges"], {"a": -2, "b": 3, "normalization": 6}
        )
        self.assertEqual(
            algebra["charge_pair_y"], {"color_block": "-1/3", "weak_block": "1/2"}
        )

    def test_menu_is_complete_with_twelve_summands(self) -> None:
        menu = self.payload["menu"]
        self.assertEqual(menu["summand_count"], 12)
        self.assertEqual(menu["nontrivial_summand_count"], 10)
        self.assertEqual(len(menu["summands"]), 12)
        self.assertEqual(sum(row["dimension"] for row in menu["summands"]), 32)
        self.assertTrue(menu["projector_completeness"]["sum_equals_identity_on_32"])
        self.assertEqual(menu["projector_completeness"]["orthogonal_idempotents"], 12)
        self.assertEqual(menu["isotypic_completeness"]["verdict"], "exact")
        charges = [
            row["charge_q6"]
            for row in menu["summands"]
            if row["component_index"] is not None
        ]
        self.assertEqual(len(charges), 10)
        self.assertEqual(len(set(charges)), 10)
        self.assertNotIn(0, charges)
        self.assertEqual(
            sorted(charges), sorted([-2, 3, -4, 1, 6, 4, -1, -6, 2, -3])
        )
        gauss = [row["summand"] for row in menu["summands"] if row["gauss_vacuum_class"]]
        self.assertEqual(gauss, ["vacuum_line", "top_line"])

    def test_menu_agrees_with_scan_table_and_lean_indexing(self) -> None:
        components = m314._scan_components(-2, 3)
        by_index = {
            row["component_index"]: row
            for row in self.payload["menu"]["summands"]
            if row["component_index"] is not None
        }
        for i, component in enumerate(components):
            self.assertEqual(by_index[i]["charge_q6"], component["q"])
            self.assertEqual(by_index[i]["color_factor"], component["color"])
            self.assertEqual(by_index[i]["weak_dimension"], component["weak"])
            parity = "odd" if component["parity"] == 1 else "even"
            self.assertEqual(by_index[i]["fermionic_parity"], parity)

    def test_field_labels_match_pinned_dimensions(self) -> None:
        by_label = {
            row["field_label"]: row["dimension"]
            for row in self.payload["menu"]["summands"]
            if row["field_label"] is not None
        }
        for name, dimension in {"Q": 6, "u_c": 3, "e_c": 1, "d_c": 3, "L": 2}.items():
            self.assertEqual(by_label[name], dimension)
            self.assertEqual(by_label[name + "_conjugate"], dimension)

    def test_subset_classification_tally(self) -> None:
        classification = self.payload["subset_classification"]
        self.assertEqual(classification["subsets_enumerated"], 1024)
        self.assertEqual(
            classification["tally"],
            {
                "empty": 1,
                "vectorlike_nonchiral": 781,
                "anomalous": 240,
                "survivor": 2,
            },
        )
        self.assertEqual(sum(classification["tally"].values()), 1024)
        breakdown = classification["anomalous_breakdown_by_failing_traces"]
        self.assertEqual(sum(breakdown.values()), 240)
        self.assertEqual(breakdown["grav+su3+su2+u1_cubed"], 118)
        self.assertEqual(classification["survivor_count"], 2)
        self.assertEqual(
            [row["mask"] for row in classification["survivors"]], [227, 796]
        )
        for row in classification["survivors"]:
            self.assertEqual(row["dimension"], 15)
            self.assertEqual(row["weak_doublet_slots"], 4)
            self.assertTrue(row["chiral"])
            self.assertEqual(
                row["anomaly_traces"],
                {"grav": 0, "su3": 0, "su2": 0, "u1_cubed": 0},
            )

    def test_class_string_indexing(self) -> None:
        classification = self.payload["subset_classification"]
        class_string = classification["class_string"]
        self.assertEqual(len(class_string), 1024)
        self.assertEqual(class_string[0], "E")
        self.assertEqual(class_string[227], "S")
        self.assertEqual(class_string[796], "S")
        # Component 0 with its conjugate 8: anomaly free and vectorlike.
        self.assertEqual(class_string[(1 << 0) | (1 << 8)], "V")
        self.assertEqual(class_string.count("S"), 2)
        self.assertEqual(class_string.count("E"), 1)
        self.assertEqual(class_string.count("V"), 781)
        self.assertEqual(class_string.count("A"), 240)

    def test_lean_cross_reference(self) -> None:
        lean = self.payload["subset_classification"]["lean_cross_reference"]
        self.assertEqual(lean["file"], "Lean/Screen/ExteriorSelection.lean")
        self.assertEqual(lean["even_mask"], 796)
        self.assertEqual(lean["odd_mask"], 227)
        self.assertIn("selection_unique", lean["theorems"])
        self.assertTrue(lean["agreement"])
        lean_source = (
            MODULE_DIR.parents[1] / "Lean" / "Screen" / "ExteriorSelection.lean"
        ).read_text(encoding="utf-8")
        self.assertIn("def evenMask : Nat := 796", lean_source)
        self.assertIn("def oddMask : Nat := 227", lean_source)
        self.assertIn("theorem selection_unique", lean_source)

    def test_off_menu_vectorlike_pair(self) -> None:
        row = self.payload["off_menu_controls"]["vectorlike_pair"]
        self.assertEqual(row["components"], [0, 8])
        self.assertFalse(row["chiral"])
        self.assertEqual(
            row["anomaly_traces"], {"grav": 0, "su3": 0, "su2": 0, "u1_cubed": 0}
        )
        self.assertEqual(row["verdict"], "anomaly_free_but_fails_chirality_clause")

    def test_off_menu_higher_representation(self) -> None:
        row = self.payload["off_menu_controls"]["higher_representation"]
        self.assertEqual(row["menu_color_factors"], ["1", "3", "3bar"])
        self.assertEqual(row["menu_weak_dimensions"], [1, 2])
        self.assertEqual(
            row["verdict"], "not_a_summand_of_the_declared_exterior_algebra"
        )

    def test_sterile_countermodel_is_source_invisible(self) -> None:
        row = self.payload["off_menu_controls"]["neutral_singlet_sterile"]
        self.assertTrue(row["all_current_observable_couplings_zero"])
        couplings = row["couplings"]
        self.assertEqual(couplings["u1_charge_readback"], 0)
        self.assertEqual(couplings["triality_class"], 0)
        self.assertEqual(couplings["duality_class"], 0)
        self.assertEqual(couplings["weak_casimir_indicator"], 0)
        self.assertEqual(
            set(couplings["anomaly_form_contributions"].values()), {0}
        )
        self.assertEqual(len(row["kernel_residue_phases"]), 6)
        self.assertTrue(
            all(item["phase_integral"] for item in row["kernel_residue_phases"])
        )
        self.assertEqual(row["verdict"], "sterile_countermodel_source_invisible")
        self.assertEqual(
            row["scalar_channel_arithmetic"]["admissible_pairing"],
            ["L", "S", "sterile_line"],
        )

    def test_deferred_and_covered_rows(self) -> None:
        controls = self.payload["off_menu_controls"]
        self.assertEqual(
            controls["scalar_duplicates_inert_doublets"]["verdict"],
            "owned_by_issue_616",
        )
        self.assertFalse(
            controls["scalar_duplicates_inert_doublets"]["analyzed_here"]
        )
        self.assertEqual(
            controls["source_invisible_direct_sums"]["verdict"],
            "covered_by_sterile_countermodel",
        )

    def test_light_heavy_threshold_declaration(self) -> None:
        threshold = self.payload["light_heavy_threshold"]
        self.assertEqual(threshold["status"], "declared_threshold")
        self.assertEqual(threshold["retained_light_sector"]["masks"], [227, 796])
        self.assertEqual(threshold["retained_light_sector"]["rank"], 15)
        self.assertEqual(
            threshold["heavy_topological_side"]["lines"],
            ["vacuum_line", "top_line"],
        )
        self.assertEqual(
            threshold["physical_decoupling_interface"]["status"],
            "separate_open_physical_interface",
        )

    def test_verdicts(self) -> None:
        verdicts = self.payload["verdicts"]
        self.assertEqual(
            verdicts["menu_completeness_inside_declared_algebra"], "exact"
        )
        self.assertEqual(
            verdicts["every_admissible_in_algebra_object_factors_through_menu"],
            "exact",
        )
        self.assertEqual(verdicts["beyond_declared_algebra"], "independence_limited")
        self.assertIn("sterile", verdicts["beyond_declared_algebra_witness"])
        self.assertEqual(
            verdicts["light_heavy_threshold"],
            "declared_with_exact_in_algebra_clauses",
        )
        self.assertEqual(
            verdicts["continuum_spin_statistics_and_laboratory_identification"],
            "separate_physical_interfaces",
        )

    def test_fail_closed_control_rows(self) -> None:
        rows = {row["name"]: row for row in self.payload["fail_closed_controls"]}
        self.assertEqual(len(rows), 4)
        self.assertEqual(
            rows["claimed_third_scan_survivor"]["actual_error"], "SURVIVOR_COUNT"
        )
        self.assertEqual(
            rows["tampered_charge_assignment"]["actual_error"],
            "PINNED_ANOMALY_TALLY",
        )
        self.assertEqual(
            rows["sterile_summand_claimed_source_visible"]["actual_error"],
            "STERILE_SOURCE_INVISIBLE",
        )
        self.assertEqual(
            rows["dropped_menu_summand"]["actual_error"], "PROJECTOR_COMPLETENESS"
        )
        self.assertTrue(all(row["passed"] for row in rows.values()))

    def test_mutation_gates_fire_directly(self) -> None:
        components = m314._scan_components(-2, 3)
        with self.assertRaises(cert.CertificateError) as caught:
            cert.classify_subsets(components, claimed_extra_survivor=273)
        self.assertEqual(caught.exception.code, "SURVIVOR_COUNT")
        tampered = m314._scan_components(-2, 3)
        tampered[3] = dict(tampered[3], q=tampered[3]["q"] + 1)
        with self.assertRaises(cert.CertificateError) as caught:
            cert.classify_subsets(tampered)
        self.assertEqual(caught.exception.code, "PINNED_ANOMALY_TALLY")
        with self.assertRaises(cert.CertificateError) as caught:
            cert.sterile_countermodel(self.pinned, claim_source_visible=True)
        self.assertEqual(caught.exception.code, "STERILE_SOURCE_INVISIBLE")
        with self.assertRaises(cert.CertificateError) as caught:
            cert.menu_ledger(self.pinned, drop_summand_index=4)
        self.assertEqual(caught.exception.code, "PROJECTOR_COMPLETENESS")

    def test_tampered_component_table_breaks_menu_lane(self) -> None:
        # A tampered charge also fails closed in the menu lane through the
        # bidegree charge check, before any projector is formed.
        tampered_pinned = dict(self.pinned, a=-2, b=4)
        with self.assertRaises(cert.CertificateError) as caught:
            cert.menu_ledger(tampered_pinned)
        self.assertEqual(caught.exception.code, "MENU_CHARGE_MISMATCH")

    def test_tampered_manifest_is_rejected(self) -> None:
        mutant = copy.deepcopy(self.stored)
        mutant["subset_classification"]["tally"]["survivor"] = 3
        with self.assertRaises(cert.CertificateError) as caught:
            cert.verify_manifest(mutant)
        self.assertEqual(caught.exception.code, "MANIFEST_MISMATCH")

    def test_pinned_manifest_hash_drift_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            (base / "manifests").mkdir()
            (base / "receipts").mkdir()
            manifest = cert.load_json(
                MODULE_DIR / "manifests" / "super_tannakian_matter_reference.json"
            )
            receipt = cert.load_json(
                MODULE_DIR / "receipts" / "super_tannakian_matter_reference.receipt.json"
            )
            manifest = copy.deepcopy(manifest)
            manifest["exterior_matter_contract"]["block_trace_charges"]["weak_block"] = "1/3"
            cert.write_json(
                base / "manifests" / "super_tannakian_matter_reference.json", manifest
            )
            cert.write_json(
                base / "receipts" / "super_tannakian_matter_reference.receipt.json",
                receipt,
            )
            with self.assertRaises(cert.CertificateError) as caught:
                cert.load_pinned_matter_lift(base)
            self.assertEqual(caught.exception.code, "PINNED_HASH_LINK")

    def test_cli_emit_and_verify(self) -> None:
        script = MODULE_DIR / "matter_menu_spectral_ledger_certificate.py"
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "ledger.json"
            subprocess.run(
                [sys.executable, str(script), "emit", "--output", str(output)],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(cert.load_json(output), self.payload)
            subprocess.run(
                [sys.executable, str(script), "verify", "--manifest", str(output)],
                check=True,
                capture_output=True,
                text=True,
            )


if __name__ == "__main__":
    unittest.main()
