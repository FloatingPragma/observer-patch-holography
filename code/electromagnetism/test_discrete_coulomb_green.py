"""Pytest for the discrete Coulomb Green producer, verifier, and guards.

Runs the producer build, the byte-level producer verification, the separate
independent replay, and mutation guards: a perturbed Green entry fails both
the producer identity check and the independent replay, and a wrong uniform
repair weight fails the pair-average identity in both implementations.
"""

from fractions import Fraction
from pathlib import Path
import json
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import discrete_coulomb_green as producer  # noqa: E402
import verify_discrete_coulomb_green_independent as verifier  # noqa: E402


def test_producer_build_matches_committed_receipt() -> None:
    committed = json.loads(producer.RECEIPT_PATH.read_bytes())
    assert committed == producer.build_receipt()
    assert producer.verify_committed_receipt()["issue"] == 733


def test_independent_replay_passes() -> None:
    assert verifier.main([]) == 0


def test_handoff_interface_is_design_only() -> None:
    committed = json.loads(producer.RECEIPT_PATH.read_bytes())
    handoff = committed["handoff_interface"]
    assert handoff["design_only"] is True
    assert handoff["frozen"] is False
    template = handoff["decision_rule_template"]
    assert template["template_only"] is True
    assert template["not_a_freeze"] is True
    assert template["verdict_labels"] == ["REPLICATED", "FAILED", "INCONCLUSIVE"]


def test_tree_solution_matches_committed_section() -> None:
    carrier_text = producer.load_pinned_source(dict(producer.SOURCE_CONTRACT[0]))
    left, right = producer.seam_table(carrier_text)
    for load_pairs in ({0: 1, 1: -1}, {3: 1, 8: -1}, {0: 2, 5: -1, 9: -1}):
        load = [Fraction(0)] * producer.PORTS
        for port, value in load_pairs.items():
            load[port] = Fraction(value)
        assert producer.tree_solution(left, right, load) == \
            verifier.committed_section(load)


def test_mutated_green_entry_fails_producer_identity() -> None:
    carrier_text = producer.load_pinned_source(dict(producer.SOURCE_CONTRACT[0]))
    left, right = producer.seam_table(carrier_text)
    lap = producer.laplacian_from_incidence(left, right)
    green = producer.green_matrix(lap)
    producer.verify_green_identities(lap, green)
    green[0][0] += Fraction(1, 180)
    with pytest.raises(producer.CoulombGreenError):
        producer.verify_green_identities(lap, green)


def test_mutated_green_entry_fails_independent_replay(tmp_path: Path) -> None:
    committed = json.loads(producer.RECEIPT_PATH.read_bytes())
    committed["green_matrix_times_180"][0][0] += 1
    # Re-sign the mutated body so only semantic recomputation can reject it.
    body = {k: v for k, v in committed.items() if k != "receipt_sha256"}
    committed["receipt_sha256"] = producer.tagged_sha256(
        producer.canonical_json_bytes(body)
    )
    mutated = tmp_path / "mutated_receipt.json"
    mutated.write_bytes(producer.canonical_json_bytes(committed))
    with pytest.raises(verifier.VerificationError):
        verifier.verify(
            verifier.DEFAULT_CARRIER, verifier.DEFAULT_GREEN_LEAN, mutated
        )


def test_wrong_repair_weight_fails_both_implementations() -> None:
    carrier_text = producer.load_pinned_source(dict(producer.SOURCE_CONTRACT[0]))
    left, right = producer.seam_table(carrier_text)
    lap = producer.laplacian_from_incidence(left, right)
    producer.pair_average_identity(left, right, lap, Fraction(1, 60))
    with pytest.raises(producer.CoulombGreenError):
        producer.pair_average_identity(left, right, lap, Fraction(1, 59))
    assert verifier.repair_identity(left, right, lap, Fraction(1, 60))
    assert not verifier.repair_identity(left, right, lap, Fraction(1, 59))


def test_pin_drift_fails_closed(tmp_path: Path) -> None:
    tampered = tmp_path / "SeamCurrentCarrierQuotient.lean"
    tampered.write_bytes(
        producer.CARRIER_LEAN_PATH.read_bytes().replace(b"![0, 0", b"![1, 0", 1)
    )
    contract = dict(producer.SOURCE_CONTRACT[0])
    contract["path"] = tampered
    with pytest.raises(producer.CoulombGreenError):
        producer.load_pinned_source(contract)
    with pytest.raises(verifier.VerificationError):
        verifier.verify(
            tampered, verifier.DEFAULT_GREEN_LEAN, verifier.DEFAULT_RECEIPT
        )
