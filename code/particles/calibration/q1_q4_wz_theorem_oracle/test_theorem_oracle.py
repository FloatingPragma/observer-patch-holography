from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path

from verify_symbolic_q1_q4_wz import run_checks


HERE = Path(__file__).resolve().parent
RER_ROOT = HERE.parents[3]
WORKSPACE_ROOT = RER_ROOT.parent


class TheoremOracleTests(unittest.TestCase):
    def test_symbolic_oracle(self) -> None:
        result = run_checks()
        self.assertTrue(result["all_pass"])
        self.assertNotIn("promotion", result)

    def test_status_is_fail_closed(self) -> None:
        status = json.loads((HERE / "theorem_status.json").read_text())
        self.assertEqual(status["classification"], "theorem_oracle_only")
        self.assertFalse(status["physical_promotion"])
        self.assertFalse(status["oph_native_qft_producer_complete"])
        self.assertFalse(status["oph_native_wz_poles"])

    def test_namespaced_parallel_dag(self) -> None:
        dag = json.loads((HERE / "dependency_dag_v2.json").read_text())
        nodes = set(dag["nodes"])
        edges = {tuple(edge) for edge in dag["edges"]}
        self.assertTrue(all(node.startswith("SM_QFT_") or node in {
            "NONPERTURBATIVE_OBSERVABLE_TOWER", "ANALYTIC_SHEET_PACKET"
        } for node in nodes))
        self.assertIn(("SM_QFT_Q1_CLASSICAL_LOCAL", "SM_QFT_Q2E_EXACT_FINITE_MEASURE"), edges)
        self.assertIn(("SM_QFT_Q1_CLASSICAL_LOCAL", "SM_QFT_Q3_BV_FORMAL"), edges)
        self.assertNotIn(("SM_QFT_Q2E_EXACT_FINITE_MEASURE", "SM_QFT_Q3_BV_FORMAL"), edges)
        self.assertNotIn(("SM_QFT_Q3_WZ_STRICT_FIXED_PARAMETER", "SM_QFT_Q4_OS_OBSERVABLE_SECTOR"), edges)
        self.assertNotIn(("SOURCE_LAW", "SM_QFT_Q3_WZ_STRICT_FIXED_PARAMETER"), edges)
        self.assertIn(("NONPERTURBATIVE_OBSERVABLE_TOWER", "SM_QFT_Q4_OS_OBSERVABLE_SECTOR"), edges)

    def test_content_addressed_upstream(self) -> None:
        upstream = json.loads((HERE / "upstream_import.json").read_text())
        archive = WORKSPACE_ROOT / upstream["archive"]
        self.assertTrue(archive.is_file())
        digest = hashlib.sha256(archive.read_bytes()).hexdigest()
        self.assertEqual(digest, upstream["archive_sha256"])
        self.assertEqual(upstream["import_policy"], "content_addressed_theorem_algebra_only_no_promotion_flags")


if __name__ == "__main__":
    unittest.main()
