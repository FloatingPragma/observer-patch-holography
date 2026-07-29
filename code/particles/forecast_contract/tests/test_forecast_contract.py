"""Fail-closed tests for the issue-639 prediction-contract machinery."""

import json
import sys
import unittest
from pathlib import Path

PACKET_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKET_DIR))

import forecast_contract as fc  # noqa: E402


def _valid_contract() -> dict:
    payload = b"sealed comparison payload"
    return {
        "schema": "oph.prediction_contract.v1",
        "contract_id": "contract-test-0001",
        "candidate_id": "c2_family_multiplicity_three",
        "output_type": "discrete_verdict",
        "allowed_ancestry": [
            {
                "artifact": "band_decomposition",
                "sha256": "sha256:" + "a" * 64,
                "role": "exact_spectral_input",
            }
        ],
        "forbidden_input_classes": ["measured_particle_masses"],
        "generator": {"module": "gen.py", "sha256": "sha256:" + "b" * 64},
        "independent_checker": {
            "module": "check.py",
            "sha256": "sha256:" + "c" * 64,
        },
        "branch_rule": "declared branch rule under the frozen interface",
        "verdict_rule": "discrete equality verdict on the frozen output",
        "stop_rule": "one evaluation, no reruns after unsealing",
        "sealed_comparison": fc.seal_comparison(payload, "held out of tree"),
        "freeze": {"frozen_by_commit": True, "unsealing_after_commit_only": True},
        "promotion_conditions": ["none before scoring"],
    }, payload


class ForecastContractTests(unittest.TestCase):
    def test_valid_contract_passes_schema(self) -> None:
        contract, _ = _valid_contract()
        fc.validate_contract(contract)

    def test_missing_field_refused(self) -> None:
        contract, _ = _valid_contract()
        del contract["stop_rule"]
        with self.assertRaises(fc.ContractError) as ctx:
            fc.validate_contract(contract)
        self.assertEqual(ctx.exception.code, "CONTRACT_MISSING_FIELDS")

    def test_unknown_field_refused(self) -> None:
        contract, _ = _valid_contract()
        contract["target_hint"] = "injected"
        with self.assertRaises(fc.ContractError) as ctx:
            fc.validate_contract(contract)
        self.assertEqual(ctx.exception.code, "CONTRACT_UNKNOWN_FIELDS")

    def test_target_injection_flagged_by_ancestry_audit(self) -> None:
        inventory = fc.load_json(fc.INVENTORY_PATH)
        poisoned = {
            "allowed_ancestry": [
                "the measured lepton triple as calibration input"
            ]
        }
        verdict = fc.check_target_ancestry(
            poisoned, inventory["forbidden_input_classes"]
        )
        self.assertFalse(verdict["ancestry_clean"])
        self.assertEqual(
            verdict["ancestry_hits"][0]["forbidden_class"],
            "measured_particle_masses",
        )

    def test_ledger_candidates_are_ancestry_clean(self) -> None:
        state = fc.build_state()
        for name, row in state["candidate_ledger"].items():
            self.assertTrue(row["ancestry_clean"], name)
        self.assertEqual(
            state["contract_freeze_status"], "NO_CANDIDATE_ELIGIBLE_YET"
        )
        self.assertEqual(state["eligible_candidates"], [])

    def test_post_freeze_mutation_refused(self) -> None:
        contract, payload = _valid_contract()
        frozen = fc.sha256_bytes(fc.canonical_json(contract).encode("utf-8"))
        contract["branch_rule"] = "silently replaced branch rule"
        with self.assertRaises(fc.ContractError) as ctx:
            fc.unseal_and_score(contract, frozen, payload,
                                contract["branch_rule"], None)
        self.assertEqual(
            ctx.exception.code, "CONTRACT_POST_FREEZE_MUTATION"
        )

    def test_branch_replacement_refused(self) -> None:
        contract, payload = _valid_contract()
        frozen = fc.sha256_bytes(fc.canonical_json(contract).encode("utf-8"))
        with self.assertRaises(fc.ContractError) as ctx:
            fc.unseal_and_score(contract, frozen, payload,
                                "a different rule applied at scoring", None)
        self.assertEqual(ctx.exception.code, "BRANCH_RULE_REPLACED")

    def test_comparison_leakage_refused(self) -> None:
        contract, _ = _valid_contract()
        frozen = fc.sha256_bytes(fc.canonical_json(contract).encode("utf-8"))
        with self.assertRaises(fc.ContractError) as ctx:
            fc.unseal_and_score(contract, frozen, b"leaked other payload",
                                contract["branch_rule"], None)
        self.assertEqual(
            ctx.exception.code, "COMPARISON_LEAKAGE_OR_TAMPER"
        )

    def test_double_unsealing_refused(self) -> None:
        contract, payload = _valid_contract()
        frozen = fc.sha256_bytes(fc.canonical_json(contract).encode("utf-8"))
        record = fc.unseal_and_score(contract, frozen, payload,
                                     contract["branch_rule"], None)
        self.assertTrue(record["unsealed"])
        with self.assertRaises(fc.ContractError) as ctx:
            fc.unseal_and_score(contract, frozen, payload,
                                contract["branch_rule"], record)
        self.assertEqual(ctx.exception.code, "COMPARISON_ALREADY_UNSEALED")

    def test_stored_state_matches_rebuild(self) -> None:
        stored = json.loads(fc.OUT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(stored, fc.build_state())


if __name__ == "__main__":
    unittest.main()
