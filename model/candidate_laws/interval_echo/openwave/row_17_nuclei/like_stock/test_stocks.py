import json
from fractions import Fraction
from itertools import combinations

from ...row_14_profile.measure import EDGES, native_geometry
from . import measure, two


def test_all_stocks_and_all_partitions():
    raw = json.loads((measure.HERE / "counts.json").read_text())
    decoded = measure.decode(raw)
    assert decoded["energy_checks"] == 42
    wins = {
        k: sum(r["k"] == k and r["separation"][0] > 0 for r in decoded["rows"])
        for k in measure.CONTRASTS
    }
    assert wins == {8: 9, 4: 7, 0: 0}
    assert len(list(measure.partitions(6))) == 11
    assert all(r["partitions"] == 10 for r in decoded["rows"] if r["n"] == 6)


def test_rational_ground_certificate_and_native_clock():
    raw = json.loads((measure.HERE / "two_counts.json").read_text())
    mask, _ = native_geometry()
    for row in raw["rows"]:
        for n, name in ((1, "one"), (2, "two")):
            words = list(map(frozenset, combinations(range(12), n)))
            x = dict(zip(words, row["positive_" + name], strict=True))
            ratios = []
            for word in words:
                edges = [
                    {int(a), int(b)} for a, b in EDGES if (a in word) != (b in word)
                ]
                v = row["k"] * (len(edges) + 3 * sum(1 for i in word if mask[i]))
                twice = (2 * v + len(edges)) * x[word] - sum(x[word ^ e] for e in edges)
                ratios.append(Fraction(twice, 2 * x[word]))
            assert [min(ratios), max(ratios)] == list(map(Fraction, row[name]))
    native, half, zero = two.decode(raw)["rows"]
    assert 0.07648605 < native["separation"][0] < native["separation"][1] < 0.07648611
    for r in (native, half):
        assert r["limits"][0][0] > r["limits"][1][1]
    assert zero["separation"] == [0, 0]
