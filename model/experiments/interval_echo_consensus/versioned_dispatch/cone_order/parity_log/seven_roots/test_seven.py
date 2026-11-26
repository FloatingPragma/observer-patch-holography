"""Independent rail decoder, Boolean forcing and archived order controls."""

from functools import cache
import gzip
import hashlib
import json

from .measure import HERE


@cache
def receipt():
    return json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))


def read(tape):
    assert len(tape) == 9
    codes = [divmod(int(n, 16), 4) for n in tape]
    assert all(pair in ((1, 2), (2, 1)) for pair in codes)
    return sum((pair[0] == 2) << i for i, pair in enumerate(codes))


def test_registry_layout_and_authentic_ancestry_in_every_sealed_history():
    r = receipt()
    mount = r["mount"]
    assert (mount["cells"], mount["capacity"]) == (320, 1914)
    assert mount["heads"][0] // 12 != mount["heads"][1] // 12
    words = mount["words"]
    occupied = {p for shelf in words for pair in shelf for p in pair}
    assert len(occupied) == 2 * 9 * 135
    assert occupied.isdisjoint(mount["heads"] + sum(mount["references"], []))
    for log in r["native"] + r["altered"]:
        assert log["final"][0] == log["final"][1]
        assert len(log["events"]) == len(log["final"][0]) == 135
        memory = {}
        incidence = []
        for reg, sources, commit in log["events"]:
            parsed = {}
            for where, first, second in sources:
                assert where < 7 and where not in parsed
                assert (where, read(first)) == memory[read(second) >> 1]
                assert first == second == log["final"][0][where]
                parsed[where] = read(first) & 1
            published = read(log["final"][0][reg])
            if reg >= 7:
                assert published & 1 == sum(parsed.values()) % 2
                incidence.append(sum(1 << i for i in parsed))
            else:
                assert parsed == {} and published & 1 == (log["word"] >> reg & 1)
            assert read(commit["before"]) == 0
            assert (
                commit["peers"] == [commit["after"]] * 2 == [log["final"][0][reg]] * 2
            )
            assert (published >> 1) == reg + 1 and (reg + 1) not in memory
            memory[reg + 1] = (reg, published)
            changed = {
                p
                for i, pair in enumerate(words[reg])
                if published >> i & 1
                for p in pair
            }
            assert set(commit["changed_ports"]) == changed
        expected = [
            33 if log["mode"] == "altered" and m == 65 else m for m in range(128)
        ]
        assert incidence == expected
        assert sum(len(log["events"][7 + m][1]) for m in range(128)) == 448


def test_every_pivotal_set_by_joint_variations_of_unselected_sources():
    r = receipt()
    table = {e["word"]: [read(s) & 1 for s in e["final"][1][7:]] for e in r["native"]}
    assert set(table) == set(range(128))
    positive = 0
    for mask in range(128):
        # Invariance on each complete joint fiber, including simultaneous omitted flips.
        fiber = {}
        for w in range(128):
            v = table[w][mask]
            assert fiber.setdefault(w & mask, v) == v
            for i in range(7):
                changes = v != table[w ^ (1 << i)][mask]
                assert changes == bool(mask & (1 << i))
                positive += changes
        assert len(fiber) == 2 ** mask.bit_count()
    assert 57344 == positive == r["analysis"]["native"]["verified_positive_flips"]
    assert r["analysis"]["native"]["pivotal"] == list(range(128))
    altered = {e["word"]: read(e["final"][1][7 + 65]) & 1 for e in r["altered"]}
    assert sorted(w for w, bit in altered.items() if bit != table[w][65]) == [32, 64]
    assert r["analysis"]["altered"]["pivotal"][65] == 33


def test_read_shutter_and_lost_id_preserve_every_other_event():
    r = receipt()
    baseline = next(e for e in r["native"] if e["word"] == 65)
    for row in r["controls"]:
        assert row["events"][:-1] == baseline["events"][:-1]
        _, samples, publication = row["events"][-1]
        assert publication is None
        assert (
            row["final"][0] == row["final"][1] == baseline["final"][0][:-1] + ["6" * 9]
        )
        assert samples[1:] == baseline["events"][-1][1][1:]
        _, a, b = samples[0]
        if row["mode"] == "shutter":
            assert read(a) == 3 and b == "0" * 9
        else:
            assert row["mode"] == "stamp" and read(a) == read(b) == 1
    for name, sha in r["sources"].items():
        assert hashlib.sha256((HERE.parents[6] / name).read_bytes()).hexdigest() == sha
