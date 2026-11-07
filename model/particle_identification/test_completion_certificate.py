"""Replay the exact finite interfaces and retain their physical scope."""
import pytest

from tools import strict_json
from .completion_certificate import OUTPUT, _same_json, report


def test_complete_interface_receipt_replays():
    assert _same_json(strict_json.load(OUTPUT), report())


@pytest.mark.parametrize("changed", [False, 0.0])
def test_receipt_equality_preserves_scalar_types(changed):
    assert not _same_json({"physical_predictions": 0},
                          {"physical_predictions": changed})


def test_receipt_equality_rejects_nonfinite_values():
    with pytest.raises(ValueError):
        _same_json({"error": float("nan")}, {"error": float("nan")})
