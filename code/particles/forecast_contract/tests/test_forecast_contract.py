"""Adversarial tests for the bounded issue-639 governance packet."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Callable

PACKET_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKET_DIR))

import forecast_contract as fc  # noqa: E402


FORBIDDEN_INPUT_CLASSES = [
    "comparison_payloads_before_unsealing",
    "consumer_policy_p_packet",
    "measured_atomic_and_si_constants",
    "measured_cosmology_and_horizon_scale",
    "measured_mixing_angles_and_widths",
    "measured_particle_masses",
    "prescribed_target_matched_simulator_charts",
    "target_residuals_and_fitted_exactifiers",
]


def _valid_contract_shape() -> dict:
    return {
        "schema": "oph.prediction_contract.v2",
        "contract_id": "contract-test-0001",
        "candidate_id": "c1_direct_public_record_capacity",
        "output_type": "discrete_verdict",
        "freeze_policy": {
            "artifact": "forecast_freeze_policy_v1.json",
            "sha256": "sha256:" + "8" * 64,
            "role": "immutable_issue_639_freeze_policy",
        },
        "candidate_inventory": {
            "artifact": "candidate_inventory.json",
            "sha256": "sha256:" + "9" * 64,
            "role": "frozen_exhaustive_candidate_inventory",
        },
        "live_fz_register": {
            "artifact": "claims/frozen_prediction_register.json",
            "sha256": "sha256:" + "7" * 64,
            "role": "live_v3_prediction_register",
        },
        "allowed_ancestry": [
            {
                "artifact": "band_decomposition.json",
                "sha256": "sha256:" + "a" * 64,
                "role": "exact_spectral_input",
            }
        ],
        "forbidden_input_classes": FORBIDDEN_INPUT_CLASSES,
        "discriminator_receipt": {
            "artifact": "discriminator.json",
            "sha256": "sha256:" + "d" * 64,
            "role": "source_visibility_and_target_ancestry",
        },
        "generator": {
            "module": "gen.py",
            "sha256": "sha256:" + "b" * 64,
        },
        "independent_checker": {
            "module": "check.py",
            "sha256": "sha256:" + "c" * 64,
        },
        "branch_rule": "declared branch rule under the frozen interface",
        "verdict_rule": "discrete equality verdict on the frozen output",
        "stop_rule": "one evaluation with no model changes after exposure",
        "sealed_comparison": {
            "payload_sha256": "sha256:" + "e" * 64,
            "byte_count": 25,
            "storage_note": "held by an independent custodian",
            "seal_scope": "integrity_and_external_custody",
            "custody_record": {
                "artifact": "custody.json",
                "sha256": "sha256:" + "f" * 64,
                "role": "comparison_access_custody",
            },
            "producer_access_before_freeze": False,
        },
        "freeze": {
            "repository_state_commit_sha": "1" * 40,
            "canonical_payload_sha256": "sha256:" + "2" * 64,
            "digest_scope": (
                "canonical_contract_without_freeze_"
                "canonical_payload_sha256"
            ),
            "external_digest_input_required": True,
            "unsealing_after_commit_only": True,
        },
        "readiness_status": (
            "DECLARED_PRE_FREEZE_CONTROLS_ONLY__NOT_SCORING_READY"
        ),
        "promotion_conditions": [
            "executable generator replay receipt",
            "checker implementation independence receipt",
            "comparison quarantine and access-control receipt",
            "durable single-use independent scoring receipt",
        ],
    }


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(
        json.dumps(payload, sort_keys=True, indent=1) + "\n",
        encoding="utf-8",
    )


def _pin(root: Path, relative: str, role: str) -> dict:
    return {
        "artifact": relative,
        "sha256": fc.sha256_bytes((root / relative).read_bytes()),
        "role": role,
    }


def _module_pin(root: Path, relative: str) -> dict:
    return {
        "module": relative,
        "sha256": fc.sha256_bytes((root / relative).read_bytes()),
    }


def _git(root: Path, *arguments: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(root), *arguments],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def _freeze_fixture(
    root: Path,
    *,
    selected_candidate: str = "c1_direct_public_record_capacity",
    inventory_mutator: Callable[[dict], None] | None = None,
) -> tuple[dict, dict, str]:
    contract = _valid_contract_shape()
    contract["candidate_id"] = selected_candidate
    (root / "band_decomposition.json").write_text(
        '{"rank":3}\n', encoding="utf-8"
    )
    (root / "gen.py").write_text("print(3)\n", encoding="utf-8")
    (root / "check.py").write_text(
        "assert int('3') == 3\n", encoding="utf-8"
    )

    ancestry_pin = _pin(
        root, "band_decomposition.json", "exact_spectral_input"
    )
    policy = fc.load_json(fc.POLICY_PATH)
    inventory = fc.load_json(fc.INVENTORY_PATH)
    inventory["completeness"] = {
        "status": "FROZEN_EXHAUSTIVE",
        "closure_criterion_frozen": True,
        **policy["closure_policy"],
    }
    inventory["selection_rule"]["status"] = "FROZEN_DETERMINISTIC"
    selected_row = inventory["candidates"][selected_candidate]
    selected_row["eligibility"] = "ELIGIBLE_SOURCE_VISIBLE"
    selected_row["description"] = "fixture candidate"
    selected_row["blocking_condition"] = "none"
    selected_row["allowed_ancestry"] = ["fixture source"]
    selected_row["allowed_ancestry_artifacts"] = [ancestry_pin]
    selected_row["public_knowledge_caveat"] = (
        "The fixture exercises governance mechanics only."
    )
    if inventory_mutator is not None:
        inventory_mutator(inventory)

    (root / "forecast_freeze_policy_v1.json").write_bytes(
        fc.POLICY_PATH.read_bytes()
    )
    live_register = fc.load_json(fc.FZ_REGISTER_PATH)
    (root / "claims").mkdir()
    (
        root / "claims" / "frozen_prediction_register.json"
    ).write_bytes(
        fc.FZ_REGISTER_PATH.read_bytes()
    )
    policy_sha256 = fc.sha256_bytes(
        (root / "forecast_freeze_policy_v1.json").read_bytes()
    )
    inventory_payload_sha256 = (
        fc.canonical_inventory_audit_payload_sha256(inventory)
    )
    _write_json(
        root / "omission_audit.json",
        {
            "schema": "oph.forecast_inventory_omission_audit.v1",
            "issue": 639,
            "inventory_payload_sha256": inventory_payload_sha256,
            "inventory_digest_scope": (
                "canonical_inventory_without_completeness_"
                "adversarial_omission_audit"
            ),
            "scope_id": policy["inventory_scope"]["scope_id"],
            "closure_criterion_id": policy["closure_policy"][
                "closure_criterion_id"
            ],
            "reviewed_candidate_ids": sorted(inventory["candidates"]),
            "reviewed_static_surface_ids": sorted(
                row["surface"]
                for row in inventory["known_surface_crosswalk"]
            ),
            "reviewed_live_register_surface_ids": sorted(
                row["surface"] for row in fc.build_fz_crosswalk(live_register)
            ),
            "freeze_policy_sha256": policy_sha256,
            "verdict": "NO_OMITTED_IN_SCOPE_CANDIDATE",
            "adversarial_review_complete": True,
            "reviewer_id": "fixture-adversarial-reviewer",
            "reviewer_attestation": (
                "I reviewed the policy-required candidates and crosswalk rows."
            ),
        },
    )
    inventory["completeness"]["adversarial_omission_audit"] = _pin(
        root,
        "omission_audit.json",
        "adversarial_inventory_omission_audit",
    )
    _write_json(root / "candidate_inventory.json", inventory)

    selected_row = inventory["candidates"][selected_candidate]
    selected_ancestry = selected_row.get(
        "allowed_ancestry_artifacts", [ancestry_pin]
    )
    reviewed_source_ids = sorted(
        pin["artifact"] for pin in selected_ancestry
    )
    _write_json(
        root / "ancestry_evidence.json",
        {
            "schema": "oph.forecast_target_ancestry_evidence.v1",
            "candidate_id": selected_candidate,
            "reviewed_source_ids": reviewed_source_ids,
            "reviewed_source_pins_sha256": fc.canonical_pin_list_sha256(
                selected_ancestry
            ),
            "forbidden_input_classes_sha256": fc.sha256_bytes(
                fc.canonical_json(
                    inventory["forbidden_input_classes"]
                ).encode("utf-8")
            ),
            "review_method": "bounded_declared_artifact_ancestry_review_v1",
            "verdict": "NO_DECLARED_FORBIDDEN_INPUT_FOUND",
            "semantic_input_closure_proved": False,
            "reviewer_id": "fixture-ancestry-reviewer",
            "reviewer_attestation": (
                "I reviewed the declared source identifiers and pinned bytes."
            ),
        },
    )
    evidence_pin = _pin(
        root,
        "ancestry_evidence.json",
        "bounded_declared_ancestry_review_evidence",
    )
    inventory_sha256 = fc.sha256_bytes(
        (root / "candidate_inventory.json").read_bytes()
    )
    live_register_sha256 = fc.sha256_bytes(
        (
            root / "claims" / "frozen_prediction_register.json"
        ).read_bytes()
    )
    _write_json(
        root / "discriminator.json",
        {
            "schema": "oph.forecast_discriminator_receipt.v1",
            "candidate_id": selected_candidate,
            "freeze_policy_sha256": policy_sha256,
            "inventory_sha256": inventory_sha256,
            "live_fz_register_sha256": live_register_sha256,
            "candidate_row_sha256": fc.canonical_candidate_row_sha256(
                selected_row
            ),
            "allowed_ancestry_artifacts_sha256": (
                fc.canonical_pin_list_sha256(selected_ancestry)
            ),
            "reviewed_source_ids": reviewed_source_ids,
            "evidence_pins": [evidence_pin],
            "verdict": (
                "ELIGIBLE_SOURCE_VISIBLE_WITH_BOUNDED_"
                "DECLARED_ANCESTRY_REVIEW"
            ),
            "source_visible": True,
            "declared_ancestry_review_complete": True,
            "semantic_input_closure_proved": False,
            "review_scope": (
                "bounded_declared_artifacts_not_transitive_semantic_closure"
            ),
        },
    )
    sealed = contract["sealed_comparison"]
    _write_json(
        root / "custody.json",
        {
            "schema": "oph.forecast_comparison_custody.v1",
            "contract_id": contract["contract_id"],
            "payload_sha256": sealed["payload_sha256"],
            "byte_count": sealed["byte_count"],
            "independent_custodian_id": "fixture-independent-custodian",
            "independence_attestation": (
                "The custodian is separate from producer and checker roles."
            ),
            "attestation_id": "fixture-custody-attestation-001",
            "storage_locator": "fixture://sealed-comparison",
            "access_history": [],
            "access_history_complete": True,
            "producer_access_before_freeze": False,
            "comparison_disclosed_before_freeze": False,
        },
    )

    _git(root, "init", "-q")
    _git(root, "config", "user.email", "forecast-test@example.invalid")
    _git(root, "config", "user.name", "Forecast Test")
    _git(root, "add", ".")
    _git(root, "commit", "-qm", "freeze fixture")
    commit_sha = _git(root, "rev-parse", "HEAD")

    contract["freeze_policy"] = _pin(
        root,
        "forecast_freeze_policy_v1.json",
        "immutable_issue_639_freeze_policy",
    )
    contract["candidate_inventory"] = _pin(
        root,
        "candidate_inventory.json",
        "frozen_exhaustive_candidate_inventory",
    )
    contract["live_fz_register"] = _pin(
        root,
        "claims/frozen_prediction_register.json",
        "live_v3_prediction_register",
    )
    contract["allowed_ancestry"] = selected_ancestry
    contract["forbidden_input_classes"] = list(
        inventory["forbidden_input_classes"]
    )
    contract["discriminator_receipt"] = _pin(
        root,
        "discriminator.json",
        "source_visibility_and_target_ancestry",
    )
    contract["generator"] = _module_pin(root, "gen.py")
    contract["independent_checker"] = _module_pin(root, "check.py")
    contract["sealed_comparison"]["custody_record"] = _pin(
        root,
        "custody.json",
        "comparison_access_custody",
    )
    contract["freeze"]["repository_state_commit_sha"] = commit_sha
    digest = fc.canonical_contract_payload_sha256(contract)
    contract["freeze"]["canonical_payload_sha256"] = digest
    return contract, inventory, digest


def _advance_freeze_commit(root: Path, contract: dict, message: str) -> str:
    _git(root, "add", ".")
    _git(root, "commit", "-qm", message)
    contract["freeze"]["repository_state_commit_sha"] = _git(
        root, "rev-parse", "HEAD"
    )
    digest = fc.canonical_contract_payload_sha256(contract)
    contract["freeze"]["canonical_payload_sha256"] = digest
    return digest


class ForecastContractTests(unittest.TestCase):
    def test_valid_contract_shape_passes_complete_schema(self) -> None:
        fc.validate_contract(_valid_contract_shape())

    def test_malformed_hash_is_refused(self) -> None:
        contract = _valid_contract_shape()
        contract["generator"]["sha256"] = "sha256:x"
        with self.assertRaises(fc.ContractError) as ctx:
            fc.validate_contract(contract)
        self.assertEqual(ctx.exception.code, "CONTRACT_SCHEMA_INVALID")

    def test_empty_forbidden_input_set_is_refused(self) -> None:
        contract = _valid_contract_shape()
        contract["forbidden_input_classes"] = []
        with self.assertRaises(fc.ContractError) as ctx:
            fc.validate_contract(contract)
        self.assertEqual(ctx.exception.code, "CONTRACT_SCHEMA_INVALID")

    def test_short_rules_and_identifiers_are_refused(self) -> None:
        for field in ("contract_id", "candidate_id", "branch_rule"):
            with self.subTest(field=field):
                contract = _valid_contract_shape()
                contract[field] = "x"
                with self.assertRaises(fc.ContractError) as ctx:
                    fc.validate_contract(contract)
                self.assertEqual(
                    ctx.exception.code, "CONTRACT_SCHEMA_INVALID"
                )

    def test_generator_cannot_be_its_own_checker(self) -> None:
        contract = _valid_contract_shape()
        contract["independent_checker"] = dict(contract["generator"])
        with self.assertRaises(fc.ContractError) as ctx:
            fc.validate_contract(contract)
        self.assertEqual(ctx.exception.code, "CHECKER_NOT_INDEPENDENT")

    def test_unimplemented_promotion_conditions_cannot_be_rewritten(self) -> None:
        contract = _valid_contract_shape()
        contract["promotion_conditions"][-1] = (
            "caller declares the contract scoring-ready"
        )
        with self.assertRaises(fc.ContractError) as ctx:
            fc.validate_contract(contract)
        self.assertEqual(
            ctx.exception.code,
            "PROMOTION_CONDITIONS_INCOMPLETE",
        )

    def test_negation_words_do_not_bypass_ancestry_warning(self) -> None:
        inventory = fc.load_json(fc.INVENTORY_PATH)
        for text in (
            "measured lepton triple without disclosure",
            "pole mass only if hidden",
            "comparison payload absent according to the author",
        ):
            with self.subTest(text=text):
                candidate = {
                    "allowed_ancestry": [text],
                    "description": "test",
                    "blocking_condition": "test",
                }
                verdict = fc.check_target_ancestry(
                    candidate, inventory["forbidden_input_classes"]
                )
                self.assertFalse(verdict["declared_vocabulary_hit_free"])
                self.assertTrue(verdict["ancestry_hits"])
                self.assertFalse(verdict["semantic_input_closure_proved"])

    def test_target_ancestry_row_is_not_reported_clean(self) -> None:
        state = fc.build_state()
        row = state["candidate_ledger"][
            "c1_koide_balanced_circulant_ratio"
        ]
        self.assertFalse(row["declared_vocabulary_hit_free"])
        self.assertTrue(row["ancestry_hits"])

    def test_fully_pinned_frozen_fixture_reaches_freeze_validation(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, external_digest = _freeze_fixture(root)
            fc.validate_contract_for_freeze(
                contract,
                repo_root=root,
                inventory=inventory,
                external_contract_digest_input=external_digest,
            )

    def test_provisional_inventory_cannot_reach_freeze(self) -> None:
        def mutate(inventory: dict) -> None:
            inventory["completeness"]["status"] = (
                "PROVISIONAL_NOT_EXHAUSTIVE"
            )

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, external_digest = _freeze_fixture(
                root, inventory_mutator=mutate
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(
            ctx.exception.code, "INVENTORY_NOT_FROZEN_EXHAUSTIVE"
        )

    def test_inventory_for_another_issue_cannot_reach_freeze(self) -> None:
        def mutate(inventory: dict) -> None:
            inventory["issue"] = 123

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, external_digest = _freeze_fixture(
                root, inventory_mutator=mutate
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(ctx.exception.code, "INVENTORY_ISSUE_MISMATCH")

    def test_forbidden_declared_text_cannot_be_marked_eligible(self) -> None:
        def mutate(inventory: dict) -> None:
            candidate = inventory["candidates"][
                "c1_direct_public_record_capacity"
            ]
            candidate["description"] = (
                "fixture candidate using a measured lepton mass"
            )

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, external_digest = _freeze_fixture(
                root, inventory_mutator=mutate
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(
            ctx.exception.code,
            "SELECTED_CANDIDATE_DECLARED_ANCESTRY_HIT",
        )

    def test_caller_cannot_weaken_anchored_forbidden_policy(self) -> None:
        def mutate(inventory: dict) -> None:
            inventory["forbidden_input_classes"] = {
                "harmless_placeholder": {
                    "description": "an attacker-selected empty policy",
                    "match_fragments": ["never appears in the candidate"],
                }
            }
            inventory["candidates"][
                "c1_direct_public_record_capacity"
            ]["description"] = (
                "fixture candidate fed the measured lepton mass"
            )

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, external_digest = _freeze_fixture(
                root, inventory_mutator=mutate
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(
            ctx.exception.code,
            "FORBIDDEN_INPUT_POLICY_MISMATCH",
        )

    def test_caller_cannot_narrow_anchored_inventory_scope(self) -> None:
        def mutate(inventory: dict) -> None:
            selected = inventory["candidates"][
                "c1_direct_public_record_capacity"
            ]
            selected["selection_priority"] = 999
            inventory["inventory_scope"] = {
                "scope_id": "attacker_selected_single_row",
                "inclusion_rule": "Only the caller-selected row is reviewed.",
                "candidate_class_legend": {"2": "caller selected"},
                "priority_policy": "The caller-selected row has priority.",
                "live_register_crosswalk": "No live surfaces are in scope.",
                "structural_known_value_policy": "No static surfaces are in scope.",
                "resource_policy": "Everything else is excluded.",
            }
            inventory["candidates"] = {
                "c1_direct_public_record_capacity": selected
            }
            inventory["known_surface_crosswalk"] = []

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, external_digest = _freeze_fixture(
                root, inventory_mutator=mutate
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(
            ctx.exception.code,
            "INVENTORY_SCOPE_POLICY_MISMATCH",
        )

    def test_contract_cannot_pin_a_rewritten_issue_policy(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, _ = _freeze_fixture(root)
            policy = fc.load_json(root / "forecast_freeze_policy_v1.json")
            policy["required_candidates"] = {
                contract["candidate_id"]: {
                    "candidate_class": 2,
                    "selection_priority": 999,
                }
            }
            _write_json(root / "forecast_freeze_policy_v1.json", policy)
            contract["freeze_policy"] = _pin(
                root,
                "forecast_freeze_policy_v1.json",
                "attacker_rewritten_issue_policy",
            )
            external_digest = _advance_freeze_commit(
                root,
                contract,
                "rewrite caller policy",
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(
            ctx.exception.code,
            "FREEZE_POLICY_NOT_CANONICAL",
        )

    def test_contract_cannot_pin_a_narrowed_live_register(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, _ = _freeze_fixture(root)
            register = fc.load_json(
                root / "claims" / "frozen_prediction_register.json"
            )
            register["rows"] = register["rows"][:1]
            register["retrospective_results"] = []
            _write_json(
                root / "claims" / "frozen_prediction_register.json",
                register,
            )
            contract["live_fz_register"] = _pin(
                root,
                "claims/frozen_prediction_register.json",
                "caller_narrowed_live_register",
            )
            external_digest = _advance_freeze_commit(
                root,
                contract,
                "narrow live register",
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(
            ctx.exception.code,
            "LIVE_REGISTER_NOT_CANONICAL",
        )

    def test_common_load_candidate_cannot_be_promoted_for_forecast(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, external_digest = _freeze_fixture(
                root,
                selected_candidate="c1_common_load_global_closure",
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(
            ctx.exception.code,
            "PERMANENTLY_INELIGIBLE_CANDIDATE_PROMOTED",
        )

    def test_archived_cap_crosswalk_cannot_be_reclassified(
        self,
    ) -> None:
        def mutate(inventory: dict) -> None:
            cap_row = next(
                row
                for row in inventory["known_surface_crosswalk"]
                if row["surface"]
                == "archived CAP-L/P/K 190-row capacity-map lattice"
            )
            cap_row["classification"] = "ELIGIBLE_NOVEL_FORECAST"
            cap_row["forecast_use"] = (
                "Select the archived aggregate as a successful forecast."
            )

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, external_digest = _freeze_fixture(
                root,
                inventory_mutator=mutate,
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(
            ctx.exception.code,
            "STATIC_CROSSWALK_POLICY_MISMATCH",
        )

    def test_unfrozen_selection_rule_cannot_reach_freeze(self) -> None:
        def mutate(inventory: dict) -> None:
            inventory["selection_rule"]["status"] = "DRAFT_UNFROZEN"

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, external_digest = _freeze_fixture(
                root, inventory_mutator=mutate
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(ctx.exception.code, "SELECTION_RULE_NOT_FROZEN")

    def test_exhaustive_label_requires_positive_pinned_omission_audit(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, _ = _freeze_fixture(root)
            _write_json(
                root / "omission_audit.json",
                {
                    "schema": "oph.forecast_inventory_omission_audit.v1",
                    "issue": 639,
                    "inventory_payload_sha256": (
                        fc.canonical_inventory_audit_payload_sha256(inventory)
                    ),
                    "inventory_digest_scope": (
                        "canonical_inventory_without_completeness_"
                        "adversarial_omission_audit"
                    ),
                    "scope_id": "test_scope",
                    "closure_criterion_id": "test-criterion-v1",
                    "reviewed_candidate_ids": [],
                    "verdict": "NO_OMITTED_IN_SCOPE_CANDIDATE",
                    "adversarial_review_complete": True,
                    "reviewer_id": "fixture-adversarial-reviewer",
                    "reviewer_attestation": (
                        "I claim to have reviewed an empty candidate list."
                    ),
                },
            )
            inventory["completeness"]["adversarial_omission_audit"] = _pin(
                root,
                "omission_audit.json",
                "adversarial_inventory_omission_audit",
            )
            _write_json(root / "candidate_inventory.json", inventory)
            _git(root, "add", ".")
            _git(root, "commit", "-qm", "record negative omission audit")
            commit_sha = _git(root, "rev-parse", "HEAD")
            contract["candidate_inventory"] = _pin(
                root,
                "candidate_inventory.json",
                "frozen_exhaustive_candidate_inventory",
            )
            contract["freeze"]["repository_state_commit_sha"] = commit_sha
            external_digest = fc.canonical_contract_payload_sha256(contract)
            contract["freeze"]["canonical_payload_sha256"] = external_digest
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(
            ctx.exception.code, "INVENTORY_OMISSION_AUDIT_NOT_POSITIVE"
        )

    def test_later_eligible_candidate_cannot_bypass_first_eligible_rule(
        self,
    ) -> None:
        def mutate(inventory: dict) -> None:
            inventory["candidates"][
                "c1_direct_public_record_capacity"
            ]["eligibility"] = "ELIGIBLE_SOURCE_VISIBLE"

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, external_digest = _freeze_fixture(
                root,
                selected_candidate="c1_screen_band_response_ratio",
                inventory_mutator=mutate,
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(ctx.exception.code, "CANDIDATE_SELECTION_DRIFT")

    def test_ineligible_selected_candidate_cannot_bypass_inventory(
        self,
    ) -> None:
        def mutate(inventory: dict) -> None:
            selected = inventory["candidates"][
                "c1_screen_band_response_ratio"
            ]
            selected["eligibility"] = "BLOCKED_SOURCE_INVISIBLE"
            inventory["candidates"][
                "c1_direct_public_record_capacity"
            ]["eligibility"] = "ELIGIBLE_SOURCE_VISIBLE"

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, external_digest = _freeze_fixture(
                root,
                selected_candidate="c1_screen_band_response_ratio",
                inventory_mutator=mutate,
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(ctx.exception.code, "CANDIDATE_NOT_ELIGIBLE")

    def test_boolean_only_discriminator_receipt_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, _ = _freeze_fixture(root)
            _write_json(
                root / "discriminator.json",
                {
                    "schema": "oph.forecast_discriminator_receipt.v1",
                    "candidate_id": contract["candidate_id"],
                    "verdict": "ELIGIBLE_SOURCE_VISIBLE",
                    "source_visible": True,
                },
            )
            contract["discriminator_receipt"] = _pin(
                root,
                "discriminator.json",
                "source_visibility_and_target_ancestry",
            )
            external_digest = _advance_freeze_commit(
                root, contract, "replace discriminator with booleans"
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(ctx.exception.code, "DISCRIMINATOR_SCHEMA")

    def test_discriminator_must_bind_inventory_and_candidate_row(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, _ = _freeze_fixture(root)
            discriminator = fc.load_json(root / "discriminator.json")
            discriminator["inventory_sha256"] = "sha256:" + "0" * 64
            _write_json(root / "discriminator.json", discriminator)
            contract["discriminator_receipt"] = _pin(
                root,
                "discriminator.json",
                "source_visibility_and_target_ancestry",
            )
            external_digest = _advance_freeze_commit(
                root, contract, "break discriminator inventory binding"
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(
            ctx.exception.code, "DISCRIMINATOR_BINDING_MISMATCH"
        )

    def test_discriminator_evidence_requires_exact_bound_schema(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, _ = _freeze_fixture(root)
            _write_json(
                root / "ancestry_evidence.json",
                {
                    "schema": "oph.forecast_target_ancestry_evidence.v1",
                    "verdict": "NO_DECLARED_FORBIDDEN_INPUT_FOUND",
                },
            )
            discriminator = fc.load_json(root / "discriminator.json")
            discriminator["evidence_pins"] = [
                _pin(
                    root,
                    "ancestry_evidence.json",
                    "bounded_declared_ancestry_review_evidence",
                )
            ]
            _write_json(root / "discriminator.json", discriminator)
            contract["discriminator_receipt"] = _pin(
                root,
                "discriminator.json",
                "source_visibility_and_target_ancestry",
            )
            external_digest = _advance_freeze_commit(
                root, contract, "replace ancestry evidence with assertion"
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(
            ctx.exception.code, "DISCRIMINATOR_EVIDENCE_NOT_POSITIVE"
        )

    def test_one_field_custody_record_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, _ = _freeze_fixture(root)
            _write_json(
                root / "custody.json",
                {"producer_access_before_freeze": False},
            )
            contract["sealed_comparison"]["custody_record"] = _pin(
                root,
                "custody.json",
                "comparison_access_custody",
            )
            external_digest = _advance_freeze_commit(
                root, contract, "replace custody record with one field"
            )
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(ctx.exception.code, "CUSTODY_RECORD_SCHEMA")

    def test_custody_record_must_bind_comparison_and_deny_access(
        self,
    ) -> None:
        for mutation, expected_code in (
            (
                {"payload_sha256": "sha256:" + "0" * 64},
                "CUSTODY_COMPARISON_BINDING_MISMATCH",
            ),
            (
                {"access_history": [{"actor": "producer"}]},
                "CUSTODY_PRE_FREEZE_ACCESS_NOT_DENIED",
            ),
        ):
            with self.subTest(mutation=mutation):
                with tempfile.TemporaryDirectory() as directory:
                    root = Path(directory)
                    contract, inventory, _ = _freeze_fixture(root)
                    custody = fc.load_json(root / "custody.json")
                    custody.update(mutation)
                    _write_json(root / "custody.json", custody)
                    contract["sealed_comparison"]["custody_record"] = _pin(
                        root,
                        "custody.json",
                        "comparison_access_custody",
                    )
                    external_digest = _advance_freeze_commit(
                        root, contract, "break custody binding"
                    )
                    with self.assertRaises(fc.ContractError) as ctx:
                        fc.validate_contract_for_freeze(
                            contract,
                            repo_root=root,
                            inventory=inventory,
                            external_contract_digest_input=external_digest,
                        )
                self.assertEqual(ctx.exception.code, expected_code)

    def test_integrity_digest_does_not_claim_access_quarantine(self) -> None:
        seal = fc.seal_comparison(b"comparison", "untrusted local note")
        self.assertEqual(seal["seal_scope"], "integrity_only")
        self.assertNotIn("producer_access_before_freeze", seal)

    def test_contract_payload_mutation_is_detected(self) -> None:
        contract = _valid_contract_shape()
        frozen = fc.canonical_contract_payload_sha256(contract)
        contract["freeze"]["canonical_payload_sha256"] = frozen
        contract["branch_rule"] = "silently replaced branch rule"
        with self.assertRaises(fc.ContractError) as ctx:
            fc.verify_frozen_contract(contract, frozen)
        self.assertEqual(
            ctx.exception.code, "CONTRACT_PAYLOAD_DIGEST_MISMATCH"
        )

    def test_contract_digest_is_not_self_referential(self) -> None:
        contract = _valid_contract_shape()
        first = fc.canonical_contract_payload_sha256(contract)
        contract["freeze"]["canonical_payload_sha256"] = first
        second = fc.canonical_contract_payload_sha256(contract)
        self.assertEqual(first, second)

    def test_rewriting_embedded_digest_cannot_bypass_external_record(
        self,
    ) -> None:
        contract = _valid_contract_shape()
        frozen = fc.canonical_contract_payload_sha256(contract)
        contract["freeze"]["canonical_payload_sha256"] = frozen
        contract["branch_rule"] = "attacker rewrote rule and embedded digest"
        contract["freeze"]["canonical_payload_sha256"] = (
            fc.canonical_contract_payload_sha256(contract)
        )
        with self.assertRaises(fc.ContractError) as ctx:
            fc.verify_frozen_contract(contract, frozen)
        self.assertEqual(
            ctx.exception.code, "EXTERNAL_FREEZE_DIGEST_INPUT_MISMATCH"
        )

    def test_external_freeze_digest_is_mandatory(self) -> None:
        contract = _valid_contract_shape()
        contract["freeze"]["canonical_payload_sha256"] = (
            fc.canonical_contract_payload_sha256(contract)
        )
        with self.assertRaises(fc.ContractError) as ctx:
            fc.verify_frozen_contract(contract, None)
        self.assertEqual(
            ctx.exception.code, "EXTERNAL_FREEZE_DIGEST_INPUT_REQUIRED"
        )

    def test_syntactic_commit_hash_must_name_a_real_commit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, _ = _freeze_fixture(root)
            contract["freeze"]["repository_state_commit_sha"] = "a" * 40
            external_digest = fc.canonical_contract_payload_sha256(contract)
            contract["freeze"]["canonical_payload_sha256"] = external_digest
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(ctx.exception.code, "FREEZE_COMMIT_NOT_FOUND")

    def test_worktree_repin_cannot_rewrite_frozen_commit_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, _ = _freeze_fixture(root)
            (root / "gen.py").write_text("print(4)\n", encoding="utf-8")
            contract["generator"] = _module_pin(root, "gen.py")
            external_digest = fc.canonical_contract_payload_sha256(contract)
            contract["freeze"]["canonical_payload_sha256"] = external_digest
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(ctx.exception.code, "COMMIT_PIN_HASH_MISMATCH")

    def test_inventory_argument_cannot_differ_from_pinned_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract, inventory, external_digest = _freeze_fixture(root)
            inventory["candidates"] = {}
            with self.assertRaises(fc.ContractError) as ctx:
                fc.validate_contract_for_freeze(
                    contract,
                    repo_root=root,
                    inventory=inventory,
                    external_contract_digest_input=external_digest,
                )
        self.assertEqual(ctx.exception.code, "INVENTORY_ARGUMENT_DRIFT")

    def test_scoring_interface_fails_closed(self) -> None:
        with self.assertRaises(fc.ContractError) as ctx:
            fc.unseal_and_score(_valid_contract_shape())
        self.assertEqual(ctx.exception.code, "SCORER_NOT_IMPLEMENTED")

    def test_state_records_unimplemented_controls_and_incomplete_inventory(
        self,
    ) -> None:
        state = fc.build_state()
        self.assertEqual(
            state["contract_freeze_status"],
            "DRAFT_GOVERNANCE_PACKET__NO_ELIGIBLE_CANDIDATE",
        )
        self.assertEqual(state["eligible_candidates"], [])
        self.assertEqual(
            state["inventory_completeness"]["status"],
            "PROVISIONAL_NOT_EXHAUSTIVE",
        )
        self.assertEqual(
            state["controls"]["executable_generator_validation"],
            "NOT_IMPLEMENTED",
        )
        self.assertEqual(
            state["controls"]["checker_true_independence_validation"],
            "NOT_IMPLEMENTED",
        )
        self.assertEqual(
            state["controls"]["comparison_access_quarantine"],
            "NOT_IMPLEMENTED",
        )
        self.assertEqual(
            state["controls"]["durable_single_use_unsealing_record"],
            "NOT_IMPLEMENTED",
        )
        self.assertEqual(
            state["controls"][
                "freeze_time_deterministic_first_eligible_selection"
            ],
            "IMPLEMENTED",
        )
        self.assertEqual(
            state["controls"][
                "non_self_referential_canonical_contract_digest"
            ],
            "IMPLEMENTED",
        )
        self.assertEqual(
            state["controls"]["durable_external_contract_digest_custody"],
            "NOT_IMPLEMENTED",
        )

    def test_known_surface_crosswalk_covers_frozen_register(self) -> None:
        state = fc.build_state()
        register = fc.load_json(fc.FZ_REGISTER_PATH)
        crosswalk = {
            row["surface"]: row
            for row in state["known_surface_crosswalk"]
        }
        for source_row in register["rows"]:
            projected = crosswalk[source_row["id"]]
            self.assertEqual(projected["registry_kind"], "prospective")
            self.assertEqual(
                projected["registry_status"], source_row["status"]
            )
        for source_row in register["retrospective_results"]:
            projected = crosswalk[source_row["id"]]
            self.assertEqual(projected["registry_kind"], "retrospective")
            self.assertEqual(
                projected["registry_status"], source_row["status"]
            )
        self.assertNotIn("FZ-04", crosswalk)
        self.assertIn("RR-506-ALPHA-HVP", crosswalk)
        self.assertIn("direct correctable-record N fixed point", crosswalk)
        self.assertIn(
            "common-load N closure and Lambda comparison", crosswalk
        )
        self.assertIn(
            "archived CAP-L/P/K 190-row capacity-map lattice", crosswalk
        )

    def test_common_load_row_is_not_promoted_as_target_clean(self) -> None:
        state = fc.build_state()
        candidate = state["candidate_ledger"][
            "c1_common_load_global_closure"
        ]
        self.assertEqual(
            candidate["eligibility"],
            (
                "INELIGIBLE_NON_BLIND_TARGET_ANCESTRY_"
                "AND_OPEN_PHYSICAL_IDENTITY"
            ),
        )
        crosswalk = next(
            row
            for row in state["known_surface_crosswalk"]
            if row["surface"] == "common-load N closure and Lambda comparison"
        )
        self.assertNotIn("target-clean screen", crosswalk["forecast_use"])
        self.assertIn(
            "not a blind target-clean forecast",
            crosswalk["forecast_use"],
        )

    def test_direct_source_n_precedes_common_load_in_class_one(self) -> None:
        inventory = fc.load_json(fc.INVENTORY_PATH)
        direct = inventory["candidates"][
            "c1_direct_public_record_capacity"
        ]
        common = inventory["candidates"]["c1_common_load_global_closure"]
        self.assertEqual(direct["candidate_class"], 1)
        self.assertEqual(common["candidate_class"], 1)
        self.assertLess(
            direct["selection_priority"], common["selection_priority"]
        )
        self.assertIn(
            "direct source-derived public-record N producer has first priority",
            inventory["inventory_scope"]["priority_policy"],
        )

    def test_potential_upgrade_routes_are_derived_from_inventory(self) -> None:
        state = fc.build_state()
        inventory = json.loads(fc.INVENTORY_PATH.read_text(encoding="utf-8"))
        expected = sorted(
            {
                issue
                for candidate in inventory["candidates"].values()
                for issue in candidate.get("potential_upgrade_issues", [])
            }
        )
        self.assertEqual(
            state["potential_upgrade_routes"]["issues"],
            expected,
        )

    def test_stored_state_matches_rebuild(self) -> None:
        stored = json.loads(fc.OUT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(stored, fc.build_state())


if __name__ == "__main__":
    unittest.main()
