import json
from itertools import combinations, product

import numpy as np
import pytest

from . import measure as m
from . import transfer


def brute_force(table):
    candidates = []
    for a,b,c in product(range(8), repeat=3):
        forward = table[a][b] * table[b][c] * table[c][a]
        backward = table[c][b] * table[b][a] * table[a][c]
        candidates.append((abs(forward - backward), (a,b,c), forward - backward))
    maximum = max(x[0] for x in candidates)
    if maximum == 0:
        return None, 0
    _, triple, gap = min(row for row in candidates if row[0] == maximum)
    return list(triple), gap


def test_public_logs_recount_the_source_and_all_selectors():
    data = json.loads((m.HERE / "results.json").read_bytes())
    pairs = list(combinations(range(18), 2))
    for history in data["histories"]:
        windows = []
        for row in history["episodes"]:
            stock = list(map(int, row["preparation"]))
            clicks = []
            for code in bytes.fromhex(row["offer_indices_hex"]):
                i,j = pairs[code]
                write = {stock[i], stock[j]} == {0,2}
                clicks.append(int(write))
                if write:
                    stock[i] = stock[j] = 1
            assert row["native"] == "".join(map(str, clicks))
            assert row["disabled"] == "0" * 60
            windows.extend(min(7, sum(clicks[i:i+20])) for i in (0,20,40))
        table = np.zeros((8,8), dtype=int)
        for a,b in zip(windows,windows[1:]):
            table[a,b] += 1
        assert table.tolist() == history["tables"]["forward"]
        assert table.T.tolist() == history["tables"]["backward"]
        for name,t in history["tables"].items():
            result = history["selectors"][name]
            assert (result["triple"], result["gap"]) == brute_force(t)
        assert history["selectors"]["within_episode"]["sign"] == -history["selectors"]["forward"]["sign"]


def test_lexical_ties_zero_and_bad_counts():
    t = [[0] * 8 for _ in range(8)]
    assert m.orientation(t)["sign"] == 0
    for a,b,c in [(0,1,2),(3,4,5)]:
        t[a][b] = t[b][c] = t[c][a] = 4
    assert m.orientation(t) == {"triple":[0,1,2], "gap":64, "sign":1}
    assert m.orientation(np.transpose(t).tolist()) == {"triple":[0,1,2], "gap":-64, "sign":-1}
    t[0][0] = -1
    with pytest.raises(ValueError):
        m.orientation(t)


def test_encoded_commands_and_held_readings():
    data = json.loads(transfer.RECEIPT.read_bytes())
    bank_count = 0
    refusal_count = 0
    for row in data["records"]:
        sign = row["selection"]["sign"]
        assert transfer.rail.decode(row["sent_clicks"]) == transfer.WORDS[sign]
        predicted_code = None if row["bridge_cut"] else transfer.WORDS[sign]
        assert transfer.rail.decode(row["received_clicks"]) == row["received_word"] == predicted_code
        if row["bridge_cut"] or sign == 0:
            assert row["no_read_per_angle"] == [8192] * 4
            assert "banks" not in row
            refusal_count += 4
            continue
        assert len(row["stock_conversion"]["writes"]) == 2
        for b in row["banks"]:
            target = (1 + sign * np.sin(np.pi * float(m.Fraction(b["theta_over_pi"])))) / 2
            assert sum(b["raw"]) == 8192 and len(b["raw"]) == 17 and sum(b["raw"][2:]) == 0
            assert abs(b["raw"][0]/8192 - target) < data["simultaneous_99_radius"]
            assert abs(b["raw"][0]/8192 - (1-target)) > data["simultaneous_99_radius"]
            bank_count += 1
    assert (bank_count, refusal_count) == (32,64)
