"""Replay the exact finite interfaces and retain their physical scope."""
import json

from .completion_certificate import OUTPUT, report


def test_complete_interface_receipt_replays():
    assert json.loads(OUTPUT.read_text()) == report()
