"""Mutation tests for the append-only audit register."""

from __future__ import annotations

import copy
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_audit_register as audit


def _data() -> dict:
    return copy.deepcopy(audit.load_json(audit.REGISTER_PATH))


def test_committed_register_and_surface_are_current() -> None:
    data = _data()
    records = audit.validate(data)
    assert audit.render(data, records).encode() == audit.SURFACE_PATH.read_bytes()


def test_unknown_reviewed_row_fails() -> None:
    data = _data()
    data["records"][0]["reviewed_observation_rows"].append("OL-Z9")
    with pytest.raises(SystemExit, match="outside the ledger"):
        audit.validate(data)


def test_attained_row_must_be_inside_reviewed_scope() -> None:
    data = _data()
    data["records"][0]["attained_rows_reviewed"].append("OL-A1")
    data["records"][0]["reviewed_observation_rows"].remove("OL-A1")
    with pytest.raises(SystemExit, match="outside the ledger/scope"):
        audit.validate(data)


def test_duplicate_finding_id_fails() -> None:
    data = _data()
    data["records"][0]["findings"][1]["id"] = data["records"][0]["findings"][0]["id"]
    with pytest.raises(SystemExit, match="finding ids"):
        audit.validate(data)
