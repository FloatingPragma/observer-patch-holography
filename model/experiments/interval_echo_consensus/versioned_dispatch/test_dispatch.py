"""Independent decoding, subcube forcing, version checks and an apparatus replay."""

import gzip
import hashlib
import json
from copy import deepcopy
from itertools import combinations, product

import pytest

from .apparatus import Board
from .measure import HERE, ROOT, lineage


def archive():
    return json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))


def bits(raw):
    assert set(raw) <= {"1", "2"}
    return [int(d == "1") for d in raw]


def number(raw):
    return sum(b << i for i, b in enumerate(raw))


def read_out(row):
    a, b = row["dispatch_raw"]
    assert a == b
    s = bits(a)
    f, g = row["continuation_raw"]
    assert f == g
    t = bits(f)
    if not any(s):
        assert not any(t)
        return None
    assert s[0] + s[2] == s[4] + s[5] == 1
    assert s[0] >= s[1] and s[2] >= s[3]
    assert row["continuation_read"] == a
    result = (s[2], s[1] + s[3], s[5])
    assert t[:2] == [1, result[1] ^ result[2]]
    assert (number(s[6:]), number(t[2:])) == (8, 9)
    return result


def test_physical_addresses_stamps_and_all_parent_edges():
    data = archive()
    for path, digest in data["sources"].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest
    for bank in data["banks"].values():
        for row in bank["runs"]:
            gate, chan, left, right, turn, spare = row["source"]
            if row["mode"] == "hidden":
                turn ^= spare
            assert read_out(row) == (
                (chan, [left, right][chan], turn) if gate else None
            )
            tape = bits(row["dispatcher_raw"])
            tags = [7 if row["overwrite"] == i else i + 1 for i in range(6)]
            assert tape[::4] == row["source"]
            assert [number(tape[4 * i + 1 : 4 * i + 4]) for i in range(6)] == tags
            expected = [[tags[r], 8, r] for r in (0, 1, 2 + chan, 4)] + [[8, 9, 6]]
            assert row["reported_edges"] == lineage(row) == (expected if gate else [])
    parked = data["parked"]
    assert set(parked["dispatcher_raw"]) == {"0"}
    assert read_out(parked) is None and parked["declared_offer_support"] is None


def test_all_completions_for_each_subset():
    # Directly enumerate completions in a different ordering from the archive producer.
    for bank in archive()["banks"].values():
        observations = {tuple(r["source"]): read_out(r) for r in bank["runs"]}
        for k in range(7):
            for subset in combinations(range(6), k):
                groups = {}
                for word, res in observations.items():
                    key = tuple(word[i] for i in subset)
                    groups.setdefault(key, set()).add(res)
                fixed = sum(2**i for i in subset)
                for word in product((0, 1), repeat=6):
                    observed = len(groups[tuple(word[i] for i in subset)]) == 1
                    certificate = bank["certificates_by_input_mask"][number(word)]
                    assert (fixed in certificate["sufficient"]) == observed
        assert all(s == [63] for s in bank["diagnostic_least_by_input_mask"])


def test_receipt_rejects_tampered_version_and_frame():
    row = next(
        r
        for r in archive()["banks"]["native_2"]["runs"]
        if r["source"] == [1, 0, 0, 1, 1, 0]
    )
    changed = deepcopy(row)
    changed["certified_versions"]["8"][2][-1] = 3
    with pytest.raises(ValueError, match="different value or writer"):
        lineage(changed)
    changed = deepcopy(row)
    snap = changed["audit_snapshots"][-1]
    changed["audit_snapshots"][-1] = ("2" if snap[0] == "1" else "1") + snap[1:]
    with pytest.raises(ValueError, match="frame or stamp"):
        lineage(changed)
    replay = Board().run(tuple(row["source"]), 2)
    assert lineage(replay) == row["reported_edges"]
    del row["reported_edges"]
    assert json.loads(json.dumps(replay)) == row


def test_correlated_enrollment_has_two_minimal_certificates():
    # Restrict the measured nominal bank to equal payloads, retaining every other bit.
    rows = archive()["banks"]["native_None"]["runs"]
    restricted = {
        tuple(r["source"]): read_out(r)
        for r in rows
        if r["source"][2] == r["source"][3]
    }
    word = (1, 0, 1, 1, 0, 0)
    sufficient = []
    for size in range(7):
        for indices in combinations(range(6), size):
            matching = [
                delivery
                for other, delivery in restricted.items()
                if all(other[i] == word[i] for i in indices)
            ]
            if all(delivery == restricted[word] for delivery in matching):
                sufficient.append(set(indices))
    minima = [s for s in sufficient if not any(t < s for t in sufficient)]
    assert minima == [{0, 1, 2, 4}, {0, 1, 3, 4}]
