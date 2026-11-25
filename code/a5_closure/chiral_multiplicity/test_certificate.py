import json
from pathlib import Path

from .certificate import run


def test_exterior_multiplicity_receipt():
    actual = json.loads(json.dumps(run(), default=int))
    archived = json.loads(Path(__file__).with_name("results.json").read_text())
    assert actual == archived
    assert actual["replicas"]["3"]["index"] == 3
    assert actual["two_generations_plus_mirror"]["index"] == 1
    assert actual["three_plus_vectorlike_charged_singlet"]["dimensions"] == 47
