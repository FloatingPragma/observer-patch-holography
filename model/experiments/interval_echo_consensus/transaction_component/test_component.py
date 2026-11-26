"""Independent public decoder, flat-state transition checks and support audit."""

import gzip
from itertools import combinations, product
import json
from pathlib import Path

from .contract import batch, launch, validate
from .measure import run


def saved():
    return json.loads(
        gzip.decompress(Path(__file__).with_name("records.json.gz").read_bytes())
    )


def test_measurement_archive_and_every_commit_from_raw_echo_pairs():
    r = saved()
    expected_edges = {(1, 3), (2, 4)}
    for trace in [*r["orders"].values(), r["hidden_control"], *r["controls"].values()]:
        raw = trace["raw"]
        cursor = 0
        for _, side, station, word in trace["reads"]:
            a = r["allocation"][side]
            b, g, x, y, z, p, s = map(int, word)
            addresses = [a["p"], a["b"], a["g"], *a[f"bank{p}"], a["s"]]
            expected = [
                [station, slot[0], 2 - bit]
                for slot, bit in zip(addresses, [p, b, g, x, y, z, s], strict=True)
            ]
            while raw[cursor : cursor + 7] != expected and cursor < len(raw):
                cursor += 1
            assert raw[cursor : cursor + 7] == expected
            cursor += 7
        for event, side, members, _, cited, old, new in trace["events"]:
            assert old[:2] == new[:2]
            if new[5:] == [1, 0]:
                assert event == side + 1 and cited == [] and new[2:5] == [old[0]] * 3
                mismatches = {i for i in range(3) if old[2 + i] != old[0]}
                assert set(members) == {
                    f"{i}{j}"
                    for i, j in combinations(range(3), 2)
                    if {i, j} & mismatches
                }
            else:
                assert old[5:] == [1, 0] and new == old[:-1] + [1] and event == side + 3
                assert cited == [side + 1] and members == ["seal"]
        if len(trace["events"]) == 4:
            assert sorted(e[0] for e in trace["events"]) == [1, 2, 3, 4]
            assert {
                (parent, e[0]) for e in trace["events"] for parent in e[4]
            } == expected_edges


def test_all_words_and_source_batch_orders_against_direct_substitution():
    words = list(product((0, 1), repeat=3))
    tested = 0
    for a, b, x, y in product((0, 1), (0, 1), words, words):
        if x == (a,) * 3 or y == (b,) * 3:
            continue
        for data_first, seal_first in product((0, 1), repeat=2):
            states = launch(a, b, x, y)
            for stage, first in enumerate((data_first, seal_first)):
                prepared = {c.side: c for c in batch(states)}
                assert set(prepared) == {0, 1}
                for side in (first, 1 - first):
                    p = prepared[side]
                    old = states[side]
                    want = (old[0], 0, *([old[0]] * 3), 1, stage)
                    assert validate(p, old) == "ok" and p.result() == want
                    states[side] = want
                assert [u[6] for u in states] == [stage, stage]
            assert not batch(states)
            tested += 1
    assert tested == 784


def test_controls_and_semantic_dependencies():
    r = saved()
    controls = r["controls"]
    assert controls["interrupted"]["attempts"] == [[1, "interrupted"], [1, "published"]]
    assert controls["replay"]["attempts"] == [[1, "published"], [1, "stale"]]
    before = controls["boundary_version"]["reads"][-2:]
    assert before[0][3] == before[1][3] == "1100000"
    assert not controls["boundary_version"]["events"]
    assert "1010010" in [row[3] for row in controls["early_publication"]["reads"]]
    assert not controls["split"]["events"] and not controls["disagreement"]["events"]
    support = r["census"]["supports"]
    assert set(support["seal_ready"]) == set("bxyzps")
    assert set(support["pair_01"]) == set("bxyps")
    assert set(support["writers"]) == set("gps")
    assert support["boundary"] == ["b", "g"]
    for case in r["orders"].values():
        for who, side, _, word in case["reads"]:
            if who != "far" or word[5] != "0":
                continue
            assert list(map(int, word)) == r["source"][side]


def test_replay_including_private_bank_control():
    r = saved()
    current = run()
    assert json.loads(json.dumps(current)) == r
    assert r["orders"]["00"]["reads"] == r["hidden_control"]["reads"]


def test_counterfeit_parent_value_and_write_are_refused():
    from copy import deepcopy
    import pytest
    from .measure import authenticate

    events = saved()["orders"]["00"]["events"]
    authenticate(events)
    for change in ("value", "support", "version"):
        other = deepcopy(events)
        if change == "value":
            other[0][6][3] ^= 1
        elif change == "support":
            other[0][3].remove("y")
        else:
            other[2][5][5] = 0
        with pytest.raises(ValueError, match="unauthenticated"):
            authenticate(other)
