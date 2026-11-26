import json
from collections import Counter
from copy import deepcopy
from itertools import product

import numpy as np
import pytest

from . import apparatus as app
from . import measure as m


def test_complete_native_support_and_blind_marginals():
    for level in range(3):
        layout = app.board_layout(level)
        n = len(layout["owners"])
        source = np.array(list(product((0, 1), repeat=2 * n)), dtype=np.int8).reshape(
            -1, n, 2
        )
        for sharing in ("shared", "split"):
            reads = m.tapes(layout, source, sharing, "repair")
            np.testing.assert_array_equal(reads[:, 0], reads[:, 1])
            probability = m.oracle(n, "half", sharing, "repair")
            assert Counter(reads[:, 0]) == {
                w: p * len(source) for w, p in probability.items()
            }
            # Deleting any one flag makes the remaining n-1 joint bits uniform.
            for held in range(n):
                observed = Counter(reads[:, 0] & ~(1 << held))
                assert set(observed.values()) == {len(source) // (1 << (n - 1))}
        for mode in ("held_parked", "unchanged_loads"):
            reads = m.tapes(layout, source, "shared", mode)
            assert Counter(reads[:, 0]) == {
                w: p * len(source)
                for w, p in m.oracle(n, "half", "shared", mode).items()
            }


def test_all_start_replay_and_missing_vs_rejection():
    saved = json.loads((m.HERE / "records.json").read_text())
    assert json.loads(json.dumps(m.acquire())) == saved
    result, radius = m.score(saved)
    for i in range(3):
        shared = result[i, "half", "shared", "repair"]
        split = result[i, "half", "split", "repair"]
        assert shared[0] / m.SHOTS - radius > split[0] / m.SHOTS + radius
    for mutation, exc in [
        ("missing", "Missing panels"),
        ("mismatch", "Endpoint meters disagree"),
        ("odd", "forbidden tape"),
        ("exposure", "Missing starts"),
    ]:
        bad = deepcopy(saved)
        held = bad["panels"][0]["counts"][0]
        if mutation == "missing":
            bad["panels"].pop()
        if mutation == "mismatch":
            held[1] ^= 1
        if mutation == "odd":
            held[0] ^= 1
            held[1] ^= 1
        if mutation == "exposure":
            held[-1] -= 1
        with pytest.raises(ValueError, match=exc):
            m.score(bad)


def test_meter_routes_and_spare_stocks_are_local():
    for level in range(3):
        layout = app.board_layout(level)
        for s in layout["stations"]:
            slots = [*s["ports"], *s["spares"], s["standard"], *sum(s["coins"], ())]
            assert len(set(slots)) == 9
            for r in s["readers"]:
                assert r["input"] not in r["route"]
                assert r["comparison"] in app.federation(level).neighbours[r["input"]]
                assert (r["route"][0], r["route"][-1]) == (
                    s["standard"],
                    r["comparison"],
                )
            assert app.load_table(s, True).tolist() == [[0, 0], [0, 0], [1, 1], [1, 1]]
            assert app.load_table(s, False).tolist() == [[0, 0], [0, 1], [1, 0], [1, 1]]
