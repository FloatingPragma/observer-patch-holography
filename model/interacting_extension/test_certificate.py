"""The complete record receipt is reconstructed from its source operators."""
import json

from .certificate import OUTPUT,report


def test_complete_record_receipt_replays():
    assert json.loads(OUTPUT.read_text())==report()
