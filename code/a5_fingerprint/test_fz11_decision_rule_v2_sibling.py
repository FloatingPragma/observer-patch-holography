from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent


def test_v2_sibling_matches_register_kill_band() -> None:
    """The in-repo v2 decision-rule sibling is the operative rule.

    The frozen prediction receipt embeds the superseded v1 rule and may
    not be mutated; this test binds the append-only sibling to the
    register's operative kill band so the two surfaces cannot drift.
    """

    sibling = json.loads(
        (HERE / "runtime" /
         "spin_six_primitive_port_prediction_decision_rule_v2.json"
         ).read_text()
    )
    register = json.loads(
        (REPO / "claims" / "frozen_prediction_register.json").read_text()
    )
    fz11 = next(r for r in register["rows"] if r["id"] == "FZ-11")
    kill = fz11["kill_band"]

    assert sibling["schema"] == "oph.fz11.decision_rule_erratum.v1"
    assert sibling["prediction_bytes_changed"] is False
    assert sibling["comparison_data_read"] is False
    # the four operative FAIL clauses appear in the register kill band
    assert "positive" in kill and "C4" in kill
    for j_clause in ("one", "five"):
        assert j_clause in kill
    for ratio in ("10/21", "32/315"):
        assert ratio in sibling["fail_at_five_sigma"][3]
        assert ratio in kill or ratio in fz11["content"]
    # the receipt the sibling corrects is the committed frozen receipt
    receipt = (HERE / "runtime" /
               "spin_six_primitive_port_prediction_receipt.json").read_bytes()
    import hashlib

    assert (
        hashlib.sha256(receipt).hexdigest()
        == sibling["prediction_receipt_sha256"]
    )
    # scope of failure stays branch-level pending issue #655
    assert "#655" in sibling["scope_of_failure"]
