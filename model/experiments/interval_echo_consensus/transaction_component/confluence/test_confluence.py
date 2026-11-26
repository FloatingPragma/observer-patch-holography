"""Checks of the registered detector stream and the complete repair relation."""

from collections import Counter
import gzip
import json
from pathlib import Path

from .measure import run


def saved():
    return json.loads(
        gzip.decompress(Path(__file__).with_name("records.json.gz").read_bytes())
    )


def test_every_far_and_local_read_from_original_material_clicks():
    receipt = saved()
    for panel in receipt["panels"].values():
        for purpose, side, station, start, stop, bits, data, count in panel[
            "observations"
        ]:
            rows = panel["raw"][start:stop]
            allocation = receipt["allocation"][side]
            assert [row[1] for row in rows[:4]] == allocation["selector"]
            markers = [1 - row[2] for row in rows[:4]]
            assert markers.count(1) == 1 and markers.index(1) == count
            reads = rows[4:]
            pairs = [allocation["boundary"], *allocation["banks"][count]]
            assert all(row[0] == station for row in rows)
            assert [row[1] for row in reads] == [p for pair in pairs for p in pair]
            words = [tuple(r[2] for r in reads[i : i + 2]) for i in range(0, 8, 2)]
            assert all(v in {(0, 1), (1, 0)} for v in words)
            assert [v[1] for v in words] == [bits, *data]
            assert sum(z != bits for z in data) == count
        for side, mask, near, far, old, new, status, route in panel["offers"]:
            assert near == far and new == old - sum(near)
            if status == "published":
                assert new < old
                assert len(route) == (3 if (old, new) == (3, 0) else 1)
    for key, word in [("left_first", [1, 1, 0]), ("right_first", [0, 1, 1])]:
        reads = [
            z[-2]
            for z in receipt["panels"][key]["observations"]
            if z[:3] == ["published", 0, 15]
        ]
        assert reads == [word, [1, 1, 1]]


def test_all_starts_paths_and_fork_counts_without_producer():
    data = saved()["graph"]
    # Six independent error flags, with the least significant position in
    # each block denoting component zero. Enumerate every useful offered mask.
    masks = [v << offset for offset in (0, 3) for v in (3, 5, 6, 7)]

    def tails(error):
        if not error:
            return Counter({0: 1})
        count = Counter()
        for mask in masks:
            if error & mask:
                for k, v in tails(error & ~mask).items():
                    count[k + 1] += v
        return count

    graph = Counter()
    for b, words, rank, enabled, hist in data["rows"]:
        error = sum(
            int(b[i] != bit) << (3 * i + j)
            for i, word in enumerate(words)
            for j, bit in enumerate(word)
        )
        assert error.bit_count() == rank
        hits = [i for i, m in enumerate(masks) if error & m]
        assert enabled == hits
        counts = tails(error)
        assert hist == [counts.get(t, 0) for t in range(7)]
        for idx, left in enumerate(hits):
            for right in hits[idx + 1 :]:
                from_left = error & ~masks[left]
                from_right = error & ~masks[right]
                union = error & ~(masks[left] | masks[right])
                assert from_left & ~masks[right] == union == from_right & ~masks[left]
                graph["peaks"] += 1
                if (
                    from_left != from_right
                    and union not in (from_left, from_right)
                    and left // 4 == right // 4
                ):
                    graph["overlapping_two_proper_legs"] += 1
    assert graph == {"peaks": 4612, "overlapping_two_proper_legs": 384}
    assert max(i for row in data["rows"] for i, num in enumerate(row[-1]) if num) == 4


def test_nonseparating_read_and_previous_scheduler_boundary():
    r = saved()
    p = r["panels"]
    a, b = p["coarse_zero"], p["coarse_one"]
    assert a["initial"][1] == b["initial"][1] == a["end"][1] == b["end"][1] == 0
    assert [x[1] for x in a["end"][0]] == [[0] * 3] * 2
    assert [x[1] for x in b["end"][0]] == [[1] * 3] * 2
    assert p["parked"]["end"] == p["parked"]["initial"]
    assert p["hidden_ones"]["observations"] == p["left_first"]["observations"]
    assert p["consistent"]["initial"] == p["consistent"]["end"]
    old = r["prior_program"]
    assert old["orders"] == {
        s: 196 for s in ("1234", "1243", "1324", "2134", "2143", "2413")
    }


def test_strict_replay():
    assert json.loads(json.dumps(run())) == saved()
