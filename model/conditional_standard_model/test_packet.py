"""Cold calculation, complete-output and inherited-custody checks."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope="module")
def fresh():
    return packet.report()


def test_complete_retained_receipt_matches_fresh_algebra_and_sources(fresh):
    stored = strict_json.load(packet.OUTPUT)
    assert packet.canonical(stored) == packet.canonical(fresh)
    assert len(fresh["catalogue"]["fermions"]) == 12
    assert len(fresh["catalogue"]["bosons"]["rows"]) == 6
    contract = fresh["theorem_contract"]
    for key in ("theorem_chain", "conditional_interfaces"):
        for row in contract[key]:
            assert all(ref["path"] in fresh["source_bindings"] for ref in row["references"])


@pytest.mark.parametrize("change", ["mass", "defect_charge", "physical_promotion", "omit_family"])
def test_retained_mass_omission_and_scientific_promotion_are_rejected(fresh, monkeypatch, change):
    value = deepcopy(fresh)
    if change == "mass":
        value["catalogue"]["fermions"][0]["tree_running_mass_exact"] = "1"
    elif change == "defect_charge":
        value["additive_defect_charge_constraints"]["general_integer_solution"]["muon"] = "2*E"
    elif change == "physical_promotion":
        value["source_selection_derived"] = True
    else:
        value["catalogue"]["fermions"].pop()
    monkeypatch.setattr(packet, "report", lambda: fresh)
    with pytest.raises(ValueError, match="full fresh calculation"):
        packet.validate(value)


def test_stale_inherited_proof_source_is_rejected(monkeypatch):
    original = strict_json.load

    def damaged(path):
        value = original(path)
        if str(path).endswith("native_matter_attachment_results.json"):
            value["source_bindings"]["model/common_source/native_matter_attachment.py"] = "0"*64
        return value

    monkeypatch.setattr(strict_json, "load", damaged)
    with pytest.raises(ValueError, match="stale inherited source"):
        packet.source_bindings()
