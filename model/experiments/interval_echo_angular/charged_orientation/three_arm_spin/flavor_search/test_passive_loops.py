import json
from pathlib import Path
from .passive_loops import run


def test_full_loop_identity_and_flux_comparator():
    record = run()
    assert record == json.loads(
        Path(__file__).with_name("passive_loops.json").read_text()
    )
    assert record["local_dimension"] == 99
    assert [r["bright_x"] for r in record["supplied_flux_control"]][::2] == [1, 0]
