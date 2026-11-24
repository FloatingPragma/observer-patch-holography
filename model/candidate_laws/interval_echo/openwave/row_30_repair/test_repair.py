"""Check the rested wave against stored repairs and a dense local propagation."""

from collections import deque
import json

import numpy as np
from scipy.linalg import expm

from model.executable_model import federation
from model.candidate_laws.interval_echo.openwave.row_28_acceleration.experiment import (
    HERE as RANGE,
)
from .experiment import HERE, SHOTS, check_signed_collars, interval


def test_old_writes_independently_predict_home_return():
    board = federation(0)
    with np.load(RANGE / "records.npz") as source:
        sign = source["sign_280210"]
    neighbors = [[] for _ in range(board.n)]
    touched = np.zeros(board.n, int)
    for (a, b), marker in zip(board.seams, sign, strict=True):
        if marker:
            touched[a] += 1
            touched[b] += 1
        else:
            neighbors[a].append(b)
            neighbors[b].append(a)
    report = json.loads((HERE / "results.json").read_text())
    panel = next(
        row for row in report["panels"] if row["cells"] == 20 and row["seed"] == 280210
    )
    visited = {panel["home"]}
    todo = deque(visited)
    while todo:
        for nxt in neighbors[todo.popleft()]:
            if nxt not in visited:
                visited.add(nxt)
                todo.append(nxt)
    sites = sorted(visited)
    assert len(sites) == panel["ports"]
    assert touched[sites].sum() == panel["ones"]
    op = np.zeros((len(sites), len(sites)))
    for i, a in enumerate(sites):
        op[i, i] = len(neighbors[a]) / 2 + 8 * (6 - touched[a])
        for b in neighbors[a]:
            op[i, sites.index(b)] -= 0.5
    wave = np.ones(len(sites)) / np.sqrt(len(sites))
    for row in panel["measurements"][:2]:
        flight = expm(-1j * row["nominal_exposure"] * op) @ wave
        miss = 1 - abs(wave @ flight) ** 2
        assert row["measured_miss"][0] < miss < row["measured_miss"][1]
        assert sum(row["hand_return_counts"]) == SHOTS


def test_sharp_line_and_changed_law_leave_distinguishable_home_readings():
    check_signed_collars()
    panels = json.loads((HERE / "results.json").read_text())["panels"]
    assert len(panels) == 27
    assert sum(p["total_components"] for p in panels) == 940
    assert all(not p["uniform_components"] for p in panels)
    for panel in panels:
        var = 64 * panel["zeros"] * panel["ones"] / panel["ports"] ** 2
        for m in panel["measurements"]:
            assert m["compatible"]
            outcomes = m["hand_return_counts"]
            assert sum(outcomes) == SHOTS
            n = outcomes[1] + outcomes[3]
            assert list(interval(n)) == m["measured_miss"]
            if m["law"] == "registered":
                assert abs(var - m["variance"]) < 1e-9
                assert m["measured_miss"][0] > 0.035
            else:
                assert n == 0


def test_joint_width_bound_allows_a_prepared_superposition():
    # A state with a nonzero spread of incident write counts need not have a
    # narrow spatial or temporal spectrum. Verify the general inequality.
    w = np.array([0.0, 1.0, 0.0, 1.0])
    g = np.diag([1.0, 2.0, 2.0, 1.0]) - np.diag(np.ones(3), 1) - np.diag(np.ones(3), -1)
    psi = np.array([1, 3j, 4, 1j]) / np.sqrt(27)

    def spread(op):
        v = op @ psi
        return np.linalg.norm(v - (psi.conj() @ v) * psi)

    assert 8 * spread(np.diag(w)) <= spread(g / 2 - 8 * np.diag(w)) + spread(g) / 2
