import json
from pathlib import Path

import numpy as np

from .measure import run


def test_crossed_material_and_pointer_instruments():
    record = run()
    receipt = json.loads(Path(__file__).with_name("results.json").read_text())
    for stage, archived in [
        (record, receipt),
        (record["filtered"], receipt["filtered"]),
    ]:
        assert stage["forecasts_overlapping"] == 16
        for row, frozen in zip(
            stage["experiments"], archived["experiments"], strict=True
        ):
            assert sum(row["counts"]) == stage["virtual_starts"]
            # Integer distributions come from sparse propagation/eigendecomposition.
            # Replays on a different BLAS may round a large bank differently.
            for key in ("remote_probability", "match_probability"):
                np.testing.assert_allclose(row[key], frozen[key], rtol=1e-6, atol=1e-17)
                lo, hi = frozen[key.replace("probability", "interval99")]
                assert lo <= row[key] <= hi
        matched = [
            r["match_interval99"]
            for r in stage["experiments"]
            if r["emitted"] == r["analyzed"]
        ]
        crossed = [
            r["match_interval99"]
            for r in stage["experiments"]
            if r["emitted"] != r["analyzed"]
        ]
        assert min(lo for lo, hi in matched) > max(hi for lo, hi in crossed)
    mixed = record["dephased"]
    assert max(row["match_interval99"][0] for row in mixed) < min(
        row["match_interval99"][1] for row in mixed
    )
