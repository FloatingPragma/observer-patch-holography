#!/usr/bin/env python3
"""Regression and boundary tests for the isolated issue #566 stage-one packet."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PRODUCER = HERE / "certify.py"
VERIFIER = HERE / "verify.py"
BASIS = HERE / "a5_alternating_bracket_reynolds_basis.json"
RECEIPT = HERE / "a5_alternating_bracket_space_stage1.receipt.json"
CARRIER = REPO / "code/a5_closure/manifests/echosahedral_federation_reference.json"


def canonical_digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True) + "\n", encoding="utf-8")


def rehash_receipt(receipt: dict[str, object]) -> None:
    receipt.pop("receipt_sha256", None)
    receipt["receipt_sha256"] = canonical_digest(receipt)


class BracketSpaceStage1Tests(unittest.TestCase):
    def copied_packet(self, root: Path) -> Path:
        packet = root / "code/a5_closure/issue_566_bracket_space_stage1"
        manifest_dir = root / "code/a5_closure/manifests"
        manifest_dir.mkdir(parents=True)
        shutil.copytree(HERE, packet)
        shutil.copy2(CARRIER, manifest_dir / CARRIER.name)
        return packet

    def run_copied_verifier(self, root: Path, packet: Path) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(packet / "verify.py")],
            cwd=root,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_generated_artifacts_are_current(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(PRODUCER), "--check"],
            cwd=REPO,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads(completed.stdout)
        self.assertEqual(summary["dimension"], 14)
        self.assertEqual(summary["basis_rank"], 14)

    def test_independent_verifier(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(VERIFIER)],
            cwd=REPO,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads(completed.stdout)
        self.assertTrue(summary["verified"])
        self.assertEqual(summary["proper_action_order"], 60)
        self.assertEqual(summary["dimension"], 14)
        self.assertEqual(summary["basis_rank"], 14)
        self.assertEqual(summary["exact_alternation_checks"], 24192)
        self.assertTrue(all(summary["independent_mutations"].values()))

    def test_basis_is_complete_reynolds_packet(self) -> None:
        basis = json.loads(BASIS.read_text(encoding="utf-8"))
        rows = basis["basis"]
        self.assertEqual(len(rows), 14)
        self.assertEqual(sum(row["orbit_size"] for row in rows), 792)
        self.assertEqual(sorted(row["orbit_size"] for row in rows), [12] + [60] * 13)
        for row in rows:
            self.assertEqual(len(row["entries"]), row["orbit_size"])
            self.assertTrue(all(abs(entry[3]) == 1 for entry in row["entries"]))
            self.assertTrue(all(entry[4] == row["orbit_size"] for entry in row["entries"]))

    def test_claim_boundary_stays_stage_one(self) -> None:
        receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
        self.assertEqual(receipt["status"], "EXACT_TARGET_FREE_SEARCH_SPACE_STAGE1_ONLY")
        self.assertTrue(receipt["target_firewall"]["enabled"])
        self.assertEqual(receipt["target_firewall"]["semantic_input_count"], 1)
        self.assertEqual(len(receipt["proper_port_action"]["permutation_rows"]), 60)
        self.assertFalse(any(receipt["later_gates"].values()))
        self.assertTrue(all(row["passed"] for row in receipt["mutation_tests"]))

    def test_independent_verifier_rejects_rehashed_receipt_semantic_tampering(self) -> None:
        def set_status(receipt: dict[str, object]) -> None:
            receipt["status"] = "ISSUE_566_CLOSED_AND_STANDARD_MODEL_FORCED"

        def set_boundary(receipt: dict[str, object]) -> None:
            receipt["claim_boundary"] = "The selected physical Standard Model bracket is proved."

        def set_environment_read(receipt: dict[str, object]) -> None:
            receipt["target_firewall"]["environment_target_read"] = True  # type: ignore[index]

        def clear_isomorphism(receipt: dict[str, object]) -> None:
            receipt["proper_port_action"]["a5_isomorphism"]["isomorphic"] = False  # type: ignore[index]

        def clear_orientation(receipt: dict[str, object]) -> None:
            receipt["proper_port_action"]["oriented_faces_preserved"] = False  # type: ignore[index]

        def change_dimension(receipt: dict[str, object]) -> None:
            receipt["representation"]["hom_ambient_dimension"] = 1  # type: ignore[index]

        def clear_covariance(receipt: dict[str, object]) -> None:
            receipt["reynolds_basis"]["exact_covariance_passed"] = False  # type: ignore[index]

        def clear_support_partition(receipt: dict[str, object]) -> None:
            receipt["reynolds_basis"]["supports_pairwise_disjoint"] = False  # type: ignore[index]

        def remove_later_gate(receipt: dict[str, object]) -> None:
            receipt["later_gates"].pop("physical_current_identified")  # type: ignore[union-attr]

        def inject_claim(receipt: dict[str, object]) -> None:
            receipt["physical_target_selected"] = True

        cases = {
            "false_closed_status": set_status,
            "false_physical_boundary": set_boundary,
            "environment_read_flag": set_environment_read,
            "a5_verdict": clear_isomorphism,
            "orientation_verdict": clear_orientation,
            "ambient_dimension": change_dimension,
            "covariance_verdict": clear_covariance,
            "support_partition": clear_support_partition,
            "missing_later_gate": remove_later_gate,
            "injected_root_claim": inject_claim,
        }
        for name, mutate in cases.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory(prefix="oph-bracket-receipt-") as raw_root:
                root = Path(raw_root)
                packet = self.copied_packet(root)
                receipt_path = packet / RECEIPT.name
                receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
                mutate(receipt)
                rehash_receipt(receipt)
                write_json(receipt_path, receipt)
                completed = self.run_copied_verifier(root, packet)
                self.assertNotEqual(completed.returncode, 0, completed.stdout)

    def test_independent_verifier_rejects_rehashed_basis_metadata_tampering(self) -> None:
        def change_issue(basis: dict[str, object]) -> None:
            basis["issue"] = 999

        def change_field(basis: dict[str, object]) -> None:
            basis["field"] = "measurement-fitted"

        def change_normalization(basis: dict[str, object]) -> None:
            basis["normalization"] = "coefficients selected from the Standard Model target"

        def change_stabilizer(basis: dict[str, object]) -> None:
            basis["basis"][0]["stabilizer_size"] = 999  # type: ignore[index]

        def inject_claim(basis: dict[str, object]) -> None:
            basis["physical_target_selected"] = True

        cases = {
            "issue": change_issue,
            "field": change_field,
            "normalization": change_normalization,
            "stabilizer": change_stabilizer,
            "injected_claim": inject_claim,
        }
        for name, mutate in cases.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory(prefix="oph-bracket-basis-") as raw_root:
                root = Path(raw_root)
                packet = self.copied_packet(root)
                basis_path = packet / BASIS.name
                receipt_path = packet / RECEIPT.name
                basis = json.loads(basis_path.read_text(encoding="utf-8"))
                receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
                mutate(basis)
                write_json(basis_path, basis)
                receipt["reynolds_basis"]["canonical_json_sha256"] = canonical_digest(basis)
                rehash_receipt(receipt)
                write_json(receipt_path, receipt)
                completed = self.run_copied_verifier(root, packet)
                self.assertNotEqual(completed.returncode, 0, completed.stdout)


if __name__ == "__main__":
    unittest.main()
