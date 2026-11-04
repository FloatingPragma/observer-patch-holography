"""Only one retained documentation-custody transition may cross frozen pins."""
import hashlib
import json
import shutil

import pytest

from tools import strict_json
from . import contract


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def renew_control_reference(root):
    path = root / contract._PREDECESSOR
    predecessor = strict_json.load(path)
    predecessor["retained_controls"][0]["sha256"] = digest(root / contract._CONTROL)
    write(path, predecessor)


@pytest.fixture
def source_tree(tmp_path, monkeypatch):
    original = contract.ROOT
    control = strict_json.load(original / contract._CONTROL)
    paths = (set(contract.REQUIRED_SOURCES) | set(control["source_sha256"])
             | {contract._CONTROL})
    for name in paths:
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original / name, target)
    monkeypatch.setattr(contract, "ROOT", tmp_path)
    return tmp_path


def test_current_metadata_transition_preserves_frozen_contract(source_tree):
    path = contract.HERE / "campaign.json"
    before = path.read_bytes()
    value = contract.load_contract()
    assert value["source_sha256"][contract._PREDECESSOR] == contract._PREDECESSOR_SHA256
    assert before == path.read_bytes()
    assert not value["integrated_gate_accepted"]
    assert value["physical_comparison_count_added"] == 0


def test_documentation_custody_can_change_with_both_references_validated(source_tree):
    (source_tree / contract._DOCUMENT).write_text("A revised research explanation.\n")
    path = source_tree / contract._CONTROL
    control = strict_json.load(path)
    control["source_sha256"][contract._DOCUMENT] = digest(source_tree / contract._DOCUMENT)
    write(path, control)
    renew_control_reference(source_tree)
    contract.load_contract()


def test_json_formatting_is_not_scientific_change(source_tree):
    control_path = source_tree / contract._CONTROL
    control_path.write_text(json.dumps(strict_json.load(control_path), sort_keys=True))
    renew_control_reference(source_tree)
    predecessor_path = source_tree / contract._PREDECESSOR
    predecessor_path.write_text(json.dumps(strict_json.load(predecessor_path), sort_keys=True))
    contract.load_contract()


@pytest.mark.parametrize("mutation", (
    "scientific_value", "executable_dependency", "source_inventory", "missing_document",
    "unrefreshed_document", "unrefreshed_control", "predecessor_scope",
    "other_retained_control", "retained_control_path", "other_predecessor_source",
))
def test_any_other_change_is_rejected(source_tree, mutation):
    control_path = source_tree / contract._CONTROL
    predecessor_path = source_tree / contract._PREDECESSOR
    control = strict_json.load(control_path)
    predecessor = strict_json.load(predecessor_path)
    if mutation == "scientific_value":
        control["scientific_mutation"] = True
    elif mutation == "executable_dependency":
        name = "model/common_source/native_defect_decision.py"
        path = source_tree / name
        path.write_bytes(path.read_bytes() + b"\n# changed executable source\n")
        control["source_sha256"][name] = digest(path)
    elif mutation == "source_inventory":
        name = "extra.py"
        (source_tree / name).write_text("pass\n")
        control["source_sha256"][name] = digest(source_tree / name)
    elif mutation == "missing_document":
        del control["source_sha256"][contract._DOCUMENT]
    elif mutation == "unrefreshed_document":
        (source_tree / contract._DOCUMENT).write_text("Unreferenced revision.\n")
    elif mutation == "unrefreshed_control":
        predecessor["retained_controls"][0]["sha256"] = "0" * 64
    elif mutation == "predecessor_scope":
        predecessor["whole_program_frozen"] = True
    elif mutation == "other_retained_control":
        predecessor["retained_controls"][1]["sha256"] = "0" * 64
    elif mutation == "retained_control_path":
        predecessor["retained_controls"][0]["path"] = contract._DOCUMENT
    elif mutation == "other_predecessor_source":
        path = source_tree / "model/completion/observer.py"
        path.write_bytes(path.read_bytes() + b"\n# changed predecessor\n")
    write(control_path, control)
    write(predecessor_path, predecessor)
    if mutation not in ("unrefreshed_control", "retained_control_path"):
        renew_control_reference(source_tree)
    with pytest.raises(ValueError, match="predecessor source changed"):
        contract.load_contract()


def test_changed_frozen_predecessor_pin_is_not_a_renewal(source_tree, tmp_path):
    value = strict_json.load(contract.HERE / "campaign.json")
    value["source_sha256"][contract._PREDECESSOR] = "0" * 64
    path = tmp_path / "forged_campaign.json"
    write(path, value)
    with pytest.raises(ValueError, match="predecessor source changed"):
        contract.load_contract(path)
