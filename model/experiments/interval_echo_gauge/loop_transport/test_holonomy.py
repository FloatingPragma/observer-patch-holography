"""Independent small configuration-space oracles and native phase/readout checks."""

from collections import deque
import json
from math import factorial
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import cell
from .holonomy import HERE, analyze_bright, ears, report
from .holonomy_records import folded_click, records, rotation


def board_from_paths(size, *paths):
    neighbours = [set() for _ in range(size)]
    for path in paths:
        for a, b in zip(path, path[1:]):
            neighbours[a].add(b)
            neighbours[b].add(a)
    return SimpleNamespace(n=size, neighbours=[sorted(x) for x in neighbours])


def enumerate_at_root(board):
    """Retain every reachable labelled board together with step count modulo four."""
    identity = tuple(range(board.n))
    todo = deque([(identity, 0)])
    seen = {(identity, 0)}
    while todo:
        labels, phase = todo.popleft()
        vacancy = labels.index(0)
        for neighbour in board.neighbours[vacancy]:
            moved = list(labels)
            moved[neighbour], moved[vacancy] = moved[vacancy], moved[neighbour]
            state = (tuple(moved), (phase + 1) % 4)
            if state not in seen:
                seen.add(state)
                todo.append(state)
    return {(labels, phase) for labels, phase in seen if labels[0] == 0}


def test_small_board_global_group_and_scalar_kernel():
    board = board_from_paths(6, [0, 1, 2, 0], [1, 3, 4, 2], [1, 5, 4])
    rooted = enumerate_at_root(board)
    perms = {x for x, _ in rooted}
    assert len(perms) == factorial(5)
    assert len(rooted) == 2 * len(perms)
    for perm, phase in rooted:
        inversions = sum(x > y for j, x in enumerate(perm) for y in perm[j + 1 :])
        assert phase % 2 == inversions % 2
    assert {phase for perm, phase in rooted if perm == tuple(range(board.n))} == {0, 2}
    triangle, extra = ears(board)
    assert len(triangle) + sum(len(ear) - 2 for ear in extra) == 6


def test_articulation_reduces_image():
    cut = board_from_paths(5, [0, 1, 2, 0], [0, 3, 4, 0])
    assert len({x for x, _ in enumerate_at_root(cut)}) == 4
    with pytest.raises(ValueError, match="just one"):
        ears(cut)


def test_changed_native_contrast_separates_coherent_paths():
    """Three possible initial mark sites, complete closed triangle pulse."""
    triangle = [0, 1, 2, 0]
    for kappa, wanted in [(8, [1, 1, 1]), (7.5, [1j, -1, -1j])]:
        gate = expm(-1j * np.pi * law.seam_term(kappa))
        phases = []
        for pos in (1, 2, 3):
            loads = [1] * 4
            loads[0] = 0
            loads[pos] = 2
            amplitude = 1
            for a, b in zip(triangle, triangle[1:]):
                first = law.SRC.index((loads[a], loads[b]))
                final = law.SRC.index((loads[b], loads[a]))
                amplitude *= gate[final, first]
                loads[a], loads[b] = loads[b], loads[a]
            phases.append(amplitude)
        np.testing.assert_allclose(phases, wanted, atol=2e-13)


def test_bright_meter_and_native_givens():
    root = 0
    board = cell()
    click, rotations = analyze_bright(board, root)
    v = np.ones(board.n, dtype=complex)
    v[root] = 0
    v /= np.linalg.norm(v)
    for _, _, a, b in rotations:
        theta = np.arctan2(np.sqrt(b), np.sqrt(a))
        c, s = np.cos(theta), np.sin(theta)
        np.testing.assert_allclose(rotation(theta), [[c, s], [-s, c]], atol=2e-14)
    assert abs(folded_click(v, click, rotations) - 1) < 1e-14
    v[:] = 0
    v[2:4] = [1 / np.sqrt(2), -1 / np.sqrt(2)]
    assert folded_click(v, click, rotations) < 1e-28


def test_registered_boards_and_counts():
    out = report()
    assert json.loads(json.dumps(out)) == json.loads(
        (HERE / "holonomy.json").read_text()
    )
    assert [v["permutation_degree"] for v in out["boards"].values()] == [
        11,
        239,
        959,
        3839,
    ]
    assert [
        out["boards"][key]["transport"]["contact_pulses"] for key in ("20", "80", "320")
    ] == [224, 210, 210]
    r = records()
    stored = json.loads((HERE / "holonomy_counts.json").read_text())
    for actual, expected in zip(r["loop_orders"], stored["loop_orders"]):
        assert actual["counts"] == expected["counts"]
    for item, p in zip(r["loop_orders"], [0.5 + np.sqrt(3) / 4, 0.5 - np.sqrt(3) / 4]):
        assert abs(item["probability"] - p) < 2e-14
        assert sum(item["counts"]) == r["shots"]
        assert abs(item["counts"][0] / r["shots"] - p) < 0.006
    np.testing.assert_allclose(
        [
            r["doublet_entanglement_fidelity"],
            r["dephased_doublet_entanglement_fidelity"],
        ],
        [1, 1 / 3],
        atol=2e-14,
    )
    np.testing.assert_allclose(
        r["ternary_probabilities"], np.array([23, 28, 13]) / 64, atol=2e-14
    )
    assert (
        r["ternary_counts"] == stored["ternary_counts"]
        and sum(r["ternary_counts"]) == r["shots"]
    )
    assert r["direct_coarse_counts"] == stored["direct_coarse_counts"]
    assert abs(r["direct_coarse_probability"] - 51 / 64) < 2e-14
