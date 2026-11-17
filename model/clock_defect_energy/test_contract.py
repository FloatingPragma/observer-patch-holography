"""Custody checks for the common diagnostic domain and exact receipts."""
import json

import pytest

from . import contract


def test_complete_predecessor_sources_match():
    value = contract.load_contract()
    assert set(value["source_sha256"]) == contract.REQUIRED_SOURCES
    assert not value["predecessor"]["merged_into_main_claimed"]


@pytest.mark.parametrize("edit", ("empty_sources", "boolean_count", "accepted", "merged"))
def test_contract_rejects_changed_scope_or_custody(tmp_path, edit):
    value = contract.load_contract()
    if edit == "empty_sources":
        value["source_sha256"] = {}
    elif edit == "changed_source":
        value["source_sha256"][next(iter(value["source_sha256"]))] = "0"*64
    elif edit == "boolean_count":
        value["physical_comparison_count_added"] = False
    elif edit == "accepted":
        value["integrated_gate_accepted"] = True
    else:
        value["predecessor"]["merged_into_main_claimed"] = True
    path = tmp_path/"campaign.json"
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError):
        contract.load_contract(path)


@pytest.mark.parametrize("text", ('{"count":false}', '{"count":0.0}', '{"count":1,"count":0}'))
def test_receipt_rejects_type_changes_and_duplicate_keys(tmp_path, text):
    path = tmp_path/"receipt.json"
    path.write_text(text)
    with pytest.raises(ValueError):
        contract.check_receipt(path, {"count": 0})


def test_receipt_allows_only_formatting_changes(tmp_path):
    path = tmp_path/"receipt.json"
    path.write_text('{ "count": 0, "accepted": false }')
    contract.check_receipt(path, {"accepted": False, "count": 0})
