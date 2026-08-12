"""Regression tests for duplicate-key-rejecting JSON loads."""

from __future__ import annotations

import pytest

import strict_json


def test_nested_duplicate_key_is_rejected() -> None:
    with pytest.raises(strict_json.DuplicateKeyError, match="duplicate JSON key 'id'"):
        strict_json.loads('{"row": {"id": "first", "id": "second"}}')


def test_distinct_keys_decode_normally() -> None:
    assert strict_json.loads('{"row": {"id": "only"}}') == {
        "row": {"id": "only"}
    }
