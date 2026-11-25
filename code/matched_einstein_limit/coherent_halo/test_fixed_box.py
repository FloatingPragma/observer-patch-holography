import json
from pathlib import Path

from model.experiments.interval_echo_einstein_station.collar_pairing.measure import (
    compare,
)
from .fixed_box import run


def test_fixed_coupling_rescaling_puncture_and_probability():
    rows = run()
    compare(rows, json.loads(Path(__file__).with_name("fixed_box.json").read_text()))
    assert all(1.94 < r["mass_times_side"] < 2.01 for r in rows)
    assert all(
        b["norm_at_deleted_site"] < 0.53 * a["norm_at_deleted_site"]
        for a, b in zip(rows, rows[1:])
    )
    for r in rows:
        p = list(r["enclosed_probabilities"].values())
        assert 0 < p[0] < p[1] < p[2] < 1
