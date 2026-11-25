"""Decode deposited word counts independently from the physical contact graph."""

import json

import numpy as np

from .. import walls
from . import measure


def test_distal_and_primary_flags_from_port_records():
    archive = json.loads(measure.PATH.read_text())
    addresses = archive["ring"] + sum(archive["wings"], [])
    assert len(set(addresses)) == 20
    ring = archive["ring"]
    anchors = [ring[i] for i in archive["anchors"]]
    wings = archive["wings"]
    expected = []
    for value in range(1024):
        bits = [int(x) for x in format(value, "010b")]
        loads = dict(zip(addresses, walls.loads(bits), strict=True))
        ordered = all(loads[ring[i]] == loads[ring[(i + 1) % 8]] for i in (1, 3, 5, 7))
        clear = all(
            loads[a] == loads[w[0]] for a, w in zip(anchors, wings, strict=True)
        )
        mismatches = [[loads[w[i]] != loads[w[i + 1]] for i in (1, 3)] for w in wings]
        primary = ordered and clear and all(sum(side) == 1 for side in mismatches)
        expected.append((ordered, primary, primary and all(s[-1] for s in mismatches)))
    for case in archive["cases"]:
        for frame in case["readings"]:
            for preparation in ("sham", "walls"):
                r = frame[preparation]
                counts = np.array(r["counts_by_word"])
                assert (
                    len(counts) == 1024
                    and counts.sum() == measure.SHOTS
                    and min(counts) >= 0
                )
                assert (counts @ expected).tolist() == [
                    r["ordered"],
                    r["primary"],
                    r["distal"],
                ]
