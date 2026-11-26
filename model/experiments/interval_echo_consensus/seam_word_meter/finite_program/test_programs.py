import gzip
import hashlib
import itertools
import json
from fractions import Fraction

import numpy as np

from .machine import Bench, READ_CODES, WIDTH, decode_word, midpoint
from .measure import HERE, check


def receipt():
    return json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))


def test_every_digit_word_and_native_inverse():
    from ..apparatus import LAW

    u = LAW.echo()
    effects = np.tile([a for a, _ in LAW.SRC], 1 + LAW.N_LAB)
    for a, b in itertools.product((0, 2), repeat=2):
        source = np.eye(99)[:, LAW.SRC.index((a, b))]
        state = u @ source
        for effect in (effects, np.tile([j for _, j in LAW.SRC], 1 + LAW.N_LAB)):
            measured = state * (effect == (a + b) // 2)
            np.testing.assert_allclose(measured, state, atol=1e-14, rtol=0)
            np.testing.assert_allclose(u @ measured, source, atol=1e-14, rtol=0)
    for a, b in itertools.product(range(1 << WIDTH), repeat=2):
        trits = [((a >> k) & 1) + ((b >> k) & 1) for k in range(WIDTH)]
        assert midpoint(trits) == (None if (a + b) % 2 else (a + b) // 2)


def test_every_recorded_route_restores_all_intermediates():
    r = receipt()
    b = Bench()
    for u, v, path in r["layout"]["paths"]:
        initial = list(range(len(b.graph)))
        moved = initial.copy()
        for a, d in path:
            assert b.graph.has_edge(a, d)
            assert set((a, d)) != set(r["layout"]["writing_seam"])
            moved[a], moved[d] = moved[d], moved[a]
        initial[u], initial[v] = initial[v], initial[u]
        assert moved == initial


def test_raw_material_events_and_source_identity():
    r = receipt()
    check(r)
    for name, digest in r["sources"].items():
        assert (
            hashlib.sha256((HERE.parents[4] / name).read_bytes()).hexdigest() == digest
        )

    def all_tapes(value):
        if isinstance(value, dict) and "tape" in value:
            yield value["tape"]
        elif isinstance(value, (dict, list)):
            entries = value.values() if isinstance(value, dict) else value
            for child in entries:
                yield from all_tapes(child)

    count = 0
    for tape in all_tapes(r):
        count += 1
        for i, event in enumerate(tape):
            op = event["op"]
            if op == "store":
                assert decode_word(event["before"]) is not None
                assert decode_word(event["after"]) == event["command"]
            if op == "mean" and i + 4 < len(tape) and tape[i + 3]["op"] == "store":
                assert event["raw"][0] == event["raw"][1]
                left, right = (decode_word(tape[i + j]["raw"]) for j in (1, 2))
                assert [READ_CODES.index(k) for k in event["raw"][0]] == [
                    ((left >> k) & 1) + ((right >> k) & 1) for k in range(WIDTH)
                ]
                for out in tape[i + 3 : i + 5]:
                    assert 2 * decode_word(out["after"]) == left + right
    assert count == 40


def test_adaptive_semantics_and_discriminating_controls():
    r = receipt()
    for row in r["adaptive"]:
        a, b = row["inputs"]
        c = (a + b) // 2
        expect = (
            [c, c, 32, 96]
            if c == 64
            else [(c + 32) // 2, c, (c + 32) // 2, 96]
            if c < 64
            else [c, (c + 96) // 2, 32, (c + 96) // 2]
        )
        assert row["final"] == [expect] * 2
    for row in r["real_menu"]:
        x, y = row["inputs"]
        coeff = [(a + b) // 2 for a, b in zip(x, y)]
        # Positive alpha=sqrt(2) lies in (7/5,10/7). This independent enclosure
        # gives a sign on every nonzero expression in this finite menu.
        u, v = coeff
        interval = (v - Fraction(10, 7) * u, v - Fraction(7, 5) * u)
        route = None if interval == (0, 0) else "high" if interval[0] > 0 else "low"
        assert interval == (0, 0) or interval[0] > 0 or interval[1] < 0
        assert row["received"] == coeff and row["chosen"] == route
    assert all(row["chosen"] is None for row in r["controls"]["fold"])
    assert r["controls"]["force_low"]["final"][0] == [56, 80, 56, 96]
    assert r["controls"]["equal_weights"]["chosen"] is None
    assert r["controls"]["no_reset"]["stages"][-1][0] == [6, 0, 0, 18]


def test_checkpoint_collisions_and_atomic_failed_real_register():
    d = json.loads(gzip.decompress((HERE / "checkpoint.json.gz").read_bytes()))
    for name, digest in d["sources"].items():
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest
    for run in d["runs"]:
        record = run["returned"]
        if run["first"] == [1, 2]:
            assert record["prefix"] == [[64]] * 4
            assert record["classification"] is None and run["final"] == [[64] * 3] * 2
        else:
            root = (run["source"][0] + run["source"][1]) // 2
            assert record["prefix"] == [[64], [root]] and record["classification"] == (
                "low" if root == 48 else "high"
            )
    np.testing.assert_allclose(d["calendar_stay_cross"], [0.5, 0.5], rtol=0, atol=1e-13)
    b = Bench()
    b.start([0, 0, 2, 1] + [0] * 8)
    before = b.loads.copy()
    assert b.mean_axes([(0, 2), (1, 3)]) == "off_grid"
    assert b.loads == before and not any(e["op"] == "store" for e in b.tape)
    small = Bench(words=[b.words[0][:2], b.words[1][:2]])
    small.start([1, 3])
    assert small.mean(0, 1) == "accepted" and small.read(0, 1) == 2
