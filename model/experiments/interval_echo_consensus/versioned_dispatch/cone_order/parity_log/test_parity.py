"""Decode both standards, the material write frames and every single-log ancestry."""

import gzip
import hashlib
import json
from copy import deepcopy
from functools import cache

import pytest

from .measure import HERE, authenticate


@cache
def receipt():
    return json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))


def peer_word(tape):
    # Each rail receives a zero and a double standard; exactly one flags.
    value = 0
    for place, nibble in enumerate(tape):
        a, b = divmod(int(nibble, 16), 4)
        assert (a, b) in {(1, 2), (2, 1)}
        value += (a == 2) << place
    return value


def test_physical_stock_and_raw_reads_reconstruct_each_complete_nominal_run():
    data = receipt()
    mount = data["mount"]
    addresses = mount["words"]
    all_ports = [p for shelf in addresses for pair in shelf for p in pair]
    assert len(set(all_ports)) == 2 * 37 * 7
    assert set(all_ports).isdisjoint(mount["heads"] + sum(mount["references"], []))
    assert mount["heads"][0] // 12 != mount["heads"][1] // 12
    for mode, bank in data["banks"].items():
        for launch in bank:
            state = [0] * 37
            for event in launch["events"]:
                put = event["publication"]
                reg = event["target"]
                assert peer_word(put["before"]) == state[reg] == 0
                for root, left, right in event["reads"]:
                    assert state[root] == peer_word(left) == peer_word(right)
                    assert state[root] >> 1 == root + 1
                new = peer_word(put["after"])
                assert [peer_word(raw) for raw in put["peers"]] == [new, new]
                assert new >> 1 == event["id"] == reg + 1
                toggle = {
                    port
                    for k, pair in enumerate(addresses[reg])
                    if new & (1 << k)
                    for port in pair
                }
                assert set(put["changed_ports"]) == toggle
                state[reg] = new
            assert [[peer_word(raw) for raw in head] for head in launch["final"]] == [
                state,
                state,
            ]
            for mask in range(32):
                actual = 3 if mode == "altered" and mask == 5 else mask
                assert state[5 + mask] & 1 == (launch["input"] & actual).bit_count() % 2


def test_sufficiency_and_pivotal_flips_from_independent_remote_values():
    data = receipt()
    different = 0
    for mode, bank in data["banks"].items():
        table = [[peer_word(x) & 1 for x in launch["final"][1][5:]] for launch in bank]
        for word, analysis in enumerate(data["analysis"][mode]):
            assert analysis["heights"] == [0] * 6 + [1] * 31
            assert sum(map(len, analysis["parents"])) == 80
            for mask in range(32):
                support = sum(
                    1 << i
                    for i in range(5)
                    if table[word ^ (1 << i)][mask] != table[word][mask]
                )
                expected = 3 if mode == "altered" and mask == 5 else mask
                assert support == analysis["least"][mask] == expected
                assert all(
                    row[mask] == table[word][mask]
                    for other, row in enumerate(table)
                    if (other ^ word) & support == 0
                )
                upper = sum(1 << c for c in range(32) if c & support == support)
                assert int(analysis["forcing_hex"][mask], 16) == upper
                expected_past = [root + 1 for root in range(5) if support >> root & 1]
                assert (
                    analysis["ancestors"][5 + mask]
                    == analysis["parents"][5 + mask]
                    == expected_past
                )
            if mode == "native":
                assert len({tuple(past) for past in analysis["ancestors"][5:]}) == 32
            else:
                original = (
                    peer_word(data["banks"]["native"][word]["final"][1][5 + 5]) & 1
                )
                different += original != table[word][5]
    assert different == 16


def test_read_shutter_and_wrong_writer_with_unchanged_value():
    nominal = receipt()["banks"]["native"][21]
    controls = {e["mode"]: e for e in receipt()["controls"]}
    for mode, launch in controls.items():
        assert launch["events"][:-1] == nominal["events"][:-1]
        event = launch["events"][-1]
        assert event["publication"] is None
        assert (
            launch["final"][0]
            == launch["final"][1]
            == nominal["final"][0][:-1] + ["6666666"]
        )
        root, a, b = event["reads"][0]
        assert root == 0
        assert event["reads"][1:] == nominal["events"][-1]["reads"][1:]
        if mode == "shutter":
            assert peer_word(a) == 3 and b == "0000000"
        else:
            assert peer_word(a) == peer_word(b) == 1
        assert peer_word(nominal["events"][-1]["reads"][0][1]) == 3


def test_source_pins_and_independent_writer_authentication():
    data = receipt()
    for path, expected in data["sources"].items():
        assert (
            hashlib.sha256((HERE.parents[5] / path).read_bytes()).hexdigest()
            == expected
        )
    nominal = data["banks"]["native"][21]
    for fake in ("6996666", "6666666", "6966666"):
        row = deepcopy(nominal)
        # Both remote tapes now claim the same version; the actual published parent
        # or its value must still distinguish the falsification.
        row["events"][-1]["reads"][0][1:] = [fake, fake]
        with pytest.raises(ValueError, match="authentication"):
            authenticate(row, data["mount"], list(range(32)))
    altered = data["banks"]["altered"][21]
    with pytest.raises(ValueError, match="Pivotal source set"):
        authenticate(altered, data["mount"], data["declared_masks"])
