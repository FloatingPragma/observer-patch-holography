from collections import Counter
from itertools import product
import json
from math import log, sqrt
import re

import networkx as nx
import numpy as np
from scipy.linalg import expm

from . import measure as m


def test_literal_alphabets_supports_and_joint_records():
    form = m.formal_source()
    _, supports, shared, zones, mapping, work = m.owner_layout()
    receipt = m.report()
    assert len(mapping) == 260
    assert len(shared) == 28
    assert len(supports[86] & supports[88]) == 25
    assert supports[86].isdisjoint(supports[247])
    assert len(supports[88] & supports[247]) == 3
    lean = (m.ROOT / "Lean/QFT/SourceCorrelationCapstone.lean").read_text()
    for o in (86, 247):
        v = re.search(rf"def support{o}Full : Fin 96 → ℕ :=\s*!\[([^\]]+)\]", lean)
        assert supports[o] == set(map(int, v[1].split(",")))
    for i, j in ((86, 88), (86, 247), (88, 247)):
        assert set(zones[i]).isdisjoint(zones[j])
        assert not {mapping[n] for n in shared} & (set(zones[i]) | set(zones[j]))
    for o, (alphabet, path) in form.items():
        assert len(path) == 32 and set(path) == set(range(len(alphabet)))
        assert work[o][0] in zones[o]
        assert len(set(alphabet)) == len(alphabet)
    assert receipt["pairs"][88]["joint_walk"] == list(
        zip(form[86][1], form[88][1], strict=True)
    )
    assert receipt["pairs"][247]["joint_walk"] == list(
        zip(form[86][1], form[247][1], strict=True)
    )
    assert receipt["pairs"][88]["steps"] == [9, 4, 4, 14]
    assert receipt["pairs"][247]["steps"] == [9, 4, 5, 13]
    assert receipt["triple_dimension"] == 2366


def test_physical_routed_copy_on_long_separated_inputs():
    board, supports, shared, zones, mapping, work = m.owner_layout()
    graph = nx.Graph(board.seams)
    h, v = m.candidate.native_law.local_term_parts()
    full = h + m.candidate.KAPPA * v
    perm = [
        m.candidate.native_law.SRC.index((b, a)) for a, b in m.candidate.native_law.SRC
    ]
    np.testing.assert_allclose(
        expm(-1j * np.pi * full)[:, :9], np.eye(99)[:, perm], atol=2e-13
    )
    for o in m.IDS:
        base, registers = work[o]
        w = 12 * base
        g = graph.subgraph(12 * c + j for c in zones[o] for j in range(12))
        lengths = nx.single_source_shortest_path_length(g, w)
        source = 12 * registers[0] + 1
        farthest = max(
            (12 * mapping[n] for n in supports[o] - set(shared)), key=lengths.get
        )
        ports = [(source, w), (farthest, w + 6), (farthest + 1, w + 10)]
        hops = []
        for a, b in ports:
            path = nx.shortest_path(g, a, b)
            go = list(zip(path, path[1:]))
            hops.extend(go + go[:-1][::-1])
        for value, bit in product((1, 2), (0, 1)):
            # Distinct spectator symbols catch failure to put any intermediate
            # register back. Only the four tested ports enter the local word.
            before = np.arange(board.n) + 3
            before[[source, farthest, farthest + 1, w + 7]] = [
                value,
                2 - 2 * bit,
                2 * bit,
                1,
            ]
            state = before.copy()
            for a, b in hops:
                state[a], state[b] = state[b], state[a]
            assert tuple(state[[w, w + 6, w + 10, w + 7]]) == (
                value,
                2 - 2 * bit,
                2 * bit,
                1,
            )
            if state[w] != 1:
                state[w + 6], state[w + 10] = state[w + 10], state[w + 6]
            for a, b in reversed(hops):
                state[a], state[b] = state[b], state[a]
            if value != 1:
                before[farthest], before[farthest + 1] = (
                    before[farthest + 1],
                    before[farthest],
                )
            np.testing.assert_array_equal(state, before)


def test_all_code_columns_field_modes_and_region_units():
    progs, shared, _ = m.setup()
    for o, p in progs.items():
        for i, label in enumerate(p["labels"]):
            physical = m.encode(p, i)
            assert np.all(physical.sum(axis=-1) == 2)
            field_words = m.from_ports(physical)
            for row, node in zip(field_words, p["supports"], strict=True):
                assert tuple(row) == (m.ZERO if node in shared else label)
            by_field = tuple(
                min((-n, x) for x, n in Counter(v).items())[1] for v in field_words.T
            )
            assert by_field == label
            rebased = physical.copy()
            rebased[0, 0, 0] += 1
            # The constant family field has a material balance reader.
            assert m.from_ports(rebased)[0, 0] == 1
        # One independently acquired position projection and swap per matrix
        # unit. Column enumeration checks the projected operator on both axes.
        n = len(p["labels"])
        eye = np.eye(n, dtype=int)
        for i, j in product(range(n), repeat=2):
            perm = np.arange(n)
            perm[[i, j]] = perm[[j, i]]
            unit = np.diag(eye[i]) @ eye[:, perm] @ np.diag(eye[j])
            np.testing.assert_array_equal(unit, np.outer(eye[i], eye[j]))
    for right in (88, 247):
        left_dim, right_dim = len(progs[86]["labels"]), len(progs[right]["labels"])
        # Every product of one left and one right unit selects a distinct
        # joint row/column; this is the full matrix basis.
        assert (
            len(
                {
                    ((i, j), (k, b))
                    for i, k, j, b in product(
                        range(left_dim),
                        range(left_dim),
                        range(right_dim),
                        range(right_dim),
                    )
                }
            )
            == (left_dim * right_dim) ** 2
        )


def test_reading_reference_and_refusal_controls():
    d = json.loads((m.HERE / "records.json").read_text())
    N = d["shots"]
    eps = sqrt(log(2 * 80 / 0.01) / (2 * N))
    for r in d["interferometers"]:
        a = np.asarray(r["counts"])
        assert a.sum() == N
        if "observer" in r:
            o = r["observer"]
            if r["arm"] == "coherent":
                assert a[-2] == N
            elif r["arm"] == "published":
                assert abs(a[0] / N - 0.5) < eps and a[0] + a[-2] == N
            else:
                assert a[-1 if o == 88 else 0] == N
        else:
            assert a[:, -2].sum() == N
            if r["publish_left"]:
                assert abs(a[0, -2] / N - 0.5) < eps and a[0, -2] + a[-2, -2] == N
            else:
                assert a[-2, -2] == N
    for o, rows in d["quorum_controls"].items():
        row = rows[0]
        assert row["zero_overlap"] == row["last_overlap"]
        assert rows[3]["zero_overlap"] != rows[3]["last_overlap"]
        assert rows[4]["zero_overlap"] == rows[4]["last_overlap"]


def test_regional_expectations_and_common_tag_control():
    r = json.loads((m.HERE / "records.json").read_text())
    paths = m.formal_source()
    samples = r["shots"]
    eps = sqrt(log(2 * 2048 / 0.01) / (2 * samples))
    for bank in r["expectations"]:
        right = bank["right"]
        start = np.zeros((13, len(paths[right][0])))
        for i, j in zip(paths[86][1], paths[right][1], strict=True):
            start[i, j] += 1 / 32
        joint = np.array(bank["joint"])
        assert joint.sum() == sum(bank["source_times"]) == samples
        for misses in bank["refusals"].values():
            assert 0 < misses[0] < samples and misses[-1] == 0
            assert all(x >= y for x, y in zip(misses, misses[1:]))
        arm = bank["arm"]
        if arm == "shared":
            expected = sum(np.roll(start, (k, k), axis=(0, 1)) for k in range(13)) / 13
        else:
            expected = start
            if arm in ("left", "both"):
                expected = np.tile(expected.mean(axis=0), (13, 1))
            if arm in ("right", "both"):
                expected = np.broadcast_to(
                    expected.mean(axis=1, keepdims=True), start.shape
                )
        assert np.max(abs(joint / samples - expected)) < eps
        if right == 247:
            target = np.trace(expected)
            assert abs(np.trace(joint) / samples - target) < eps
            if arm == "shared":
                assert target > 0.3
            elif arm == "both":
                assert target < 0.08
