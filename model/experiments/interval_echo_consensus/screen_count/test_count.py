"""Decode the position-resolved rail records without the acquisition simulator."""

import gzip
from itertools import combinations
import json
from pathlib import Path
import numpy as np
import pytest

from .readback import HERE, summarize


def receipt(stem):
    with gzip.open(HERE / (stem + ".json.gz"), "rt") as f:
        return json.load(f)


def digit(raw):
    return {1: 0, 0: 1, 2: 2}[raw]


def word(raw, n):
    assert raw[0] == raw[1] and len(raw[1]) == 12 * n + 1
    trits = list(map(digit, raw[0]))
    return sum((sum(trits[j * 12 : j * 12 + 12]) % 2) << j for j in range(n)), trits[-1]


def test_bin_and_receiver_replay():
    archive = receipt("records")
    for row in archive:
        cfg = row["config"]
        owners = cfg["owners"]
        n = len(owners)
        bonds = cfg["inside"]
        root = owners.index(cfg["boundary"][0] // 12)
        masks = np.arange(1 << len(bonds))
        transfer = np.zeros(len(masks), dtype=int)
        for j, (a, b) in enumerate(bonds):
            bit = (masks >> j) & 1
            transfer ^= np.where(
                bit, 1 << owners.index(b // 12), 1 << owners.index(a // 12)
            )
        for sta in row["stations"]:
            assert (
                sta["sharp"] == [1, 0, 2] and sta["shuttered"][0] == sta["shuttered"][1]
            )
            for rec in sta["paths"]:
                assert rec["path"][0] == rec["port"] and rec["path"][-1] == sta["head"]
        for panel in row["panels"].values():
            observed = np.array(panel["rail_counts"])
            inferred = np.zeros((2, 1 << n), dtype=int)
            for dep in (0, 1):
                np.add.at(inferred[1 - dep], transfer ^ (dep << root), observed[dep])
            np.testing.assert_array_equal(inferred.ravel(), panel["joint"])
        if n == 6:
            positive = set()
            for r in row["targets"]:
                if r["raw"] is not None:
                    positive.add(r["target"])
                    assert word(r["raw"], n) == (r["target"], 1 - r["deposit"])
            assert (
                len(positive) == 64
                and sum(r["raw"] is None for r in row["targets"]) == 64
            )
        if row["cycle"]:
            x, y = row["cycle"]
            assert x != y and word(x, n) == word(y, n)
    assert summarize(archive) == json.loads((HERE / "summary.json").read_text())
    # Reassign a count inside the physically forbidden charge component.
    panel = archive[0]["panels"]["parked"]["joint"]
    full = next(i for i, v in enumerate(panel) if v > 0)
    panel[full] -= 1
    panel[full ^ 1] += 1
    with pytest.raises(ValueError, match="Out-of-support"):
        summarize(archive)


def test_hidden_checkpoint_capacity(monkeypatch):
    r = receipt("checkpoint")
    bank = [word(raw, 8)[0] for raw in r["raw"]]
    assert len(bank) == len(set(bank)) == 256
    hamming = {0, 1, 2, 4, 8, 16, 32, 64}
    for src, columns in enumerate(r["outcomes"]):
        assert {src ^ y for y in columns["unadjusted"]} == hamming
        assert all((z >> 7) == src >> 7 for z in columns["unadjusted"])
    for dep in (0, 1):
        reachable = set(range(128 * dep, 128 * dep + 128))
        messages = {r["outcomes"][x]["adjusted"][0] for x in reachable}
        balls = [set(r["outcomes"][x]["unadjusted"]) for x in messages]
        assert len(balls) == 16 and all(len(s) == 8 for s in balls)
        assert set.union(*balls) == reachable and sum(map(len, balls)) == len(reachable)
        assert all(len(a & b) == 0 for a, b in combinations(balls, 2))
    control = r["controls"]
    assert bank[control["diagonal"]] != bank[control["ordinary"]]
    assert (
        r["outcomes"][control["other_parent"]]["adjusted"][control["other_error"]]
        == control["two_errors"]
    )
    # The finite evaluator receives one joint stochastic channel: the edge label
    # is hidden. Distinct possible errors for two inputs must be allowed.
    monkeypatch.syspath_prepend(
        str(Path(__file__).resolve().parents[4] / "code" / "capacity_readback")
    )
    from correctable_public_record_capacity import compound_confusability_graph

    labels = [str(w) for w in bank]
    rows = {
        labels[x]: {labels[y]: 0.125 for y in r["outcomes"][x]["unadjusted"]}
        for x in range(256)
    }
    graph = compound_confusability_graph(
        labels, [{"name": "unread edge", "rows": rows}]
    )
    for x in range(256):
        expected = {
            labels[y]
            for y in range(256)
            if x >> 7 == y >> 7 and 1 <= (x ^ y).bit_count() <= 2
        }
        assert graph[labels[x]] == expected
    for panel in r["banks"]:
        terminal = [0] * 256
        for i, ct in enumerate(panel["starts"]):
            y = r["outcomes"][i][panel["source"]][panel["continuation"]]
            terminal[bank[y]] += ct
        assert terminal == panel["terminal_words"] and sum(terminal) == 2**17
        multiplicity = (2 if panel["boundary"] == "half" else 1) * (
            16 if panel["source"] == "adjusted" else 128
        )
        assert sum(k > 0 for k in terminal) == multiplicity
