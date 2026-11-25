"""Source, geometry, complete observer records and independent 27-word replay."""

import json
from itertools import combinations, product

import numpy as np
from scipy.linalg import expm
import sympy as sp

from model.executable_model import federation
from model.experiments.interval_echo_angular.charged_orientation.conical_cycle import measure as m


def test_native_face_source_and_doublet():
    board = federation(0)
    face = m.private_face()
    assert face == (0, 5, 11)
    assert face in board.faces() and all(board.degree[i] == 5 for i in face)
    for sign in (-1, 1):
        edges, words = m.terms(sign)
        h = sum(edges.values())
        np.testing.assert_allclose(h, (2*m.KAPPA + 1.5) * np.eye(3) - np.ones((3, 3))/2)
        np.testing.assert_array_equal(words.sum(axis=1) - 3, sign*np.ones(3))
        odd = np.array([0, 1, -1]) / np.sqrt(2)
        assert abs(abs(np.vdot(odd, m.prepare(edges, words))) - 1) < 1e-13


def test_first_order_cone_symbolically():
    c, s = sp.symbols("c s", real=True)
    diag = sp.diag(c, (-c + sp.sqrt(3)*s)/2, (-c - sp.sqrt(3)*s)/2)
    even = sp.Matrix([2, -1, -1])/sp.sqrt(6)
    odd = sp.Matrix([0, 1, -1])/sp.sqrt(2)
    e = sp.Matrix.hstack(even, odd)
    assert sp.simplify(e.T*diag*e) == sp.Matrix([[c, -s], [-s, -c]]) / 2


def test_full_27_word_replay_at_finite_rate():
    """Use the analytic native swap on every word, with all common phases present."""
    basis = list(product(range(3), repeat=3))
    identity = np.eye(len(basis))
    h = np.zeros_like(identity)
    words = np.asarray(basis)
    for a, b in combinations(range(3), 2):
        swaps = words.copy()
        swaps[:, [a, b]] = words[:, [b, a]]
        perm = [basis.index(tuple(word)) for word in swaps]
        h += (identity - identity[:, perm])/2 + m.KAPPA*np.diag((words[:, a] - words[:, b])**2)
    step = expm(-1j * m.STEP * h)
    directions = np.pi * (np.arange(64) + 0.5) / 64
    for sign in (-1, 1):
        edges, reduced = m.terms(sign)
        ports = [basis.index(tuple(word)) for word in reduced]
        source = np.zeros(len(basis), dtype=complex)
        source[ports] = m.prepare(edges, reduced)
        x, y = np.cos(m.ANGLES), np.sin(m.ANGLES)
        arms = []
        for direction in (1, -1):
            arm = source.copy()
            for t in directions:
                dwell = m.STEP*m.STRENGTH*(3 + x*np.cos(t) + direction*y*np.sin(t))/m.KAPPA
                kick = np.exp(-0.5j * m.KAPPA * (words**2 @ dwell))
                arm = kick*(step @ (kick*arm))
            np.testing.assert_allclose(np.delete(arm, ports), 0, atol=1e-12)
            arms.append(arm[ports])
        joint = m.probabilities(*arms)
        receipt = json.loads((m.HERE / "records.json").read_text())
        np.testing.assert_allclose(joint, receipt["runs"][str(sign)]["enclosing"]["64"]["probabilities"],
                                   atol=1e-11, rtol=0)


def test_deposited_terminal_banks_and_controls():
    archive = json.loads((m.HERE / "records.json").read_text())
    replay = m.acquire()
    assert archive["face"] == replay["face"]
    # 52 banks, six bins on each driven bank and 24 on spatial comparisons.
    # Conservatively bound all observed bins and the dark marginals together.
    radius = np.sqrt(np.log(2*500/0.01) / (2*m.SHOTS))

    def check(actual, stored):
        p = np.asarray(actual["probabilities"])
        np.testing.assert_allclose(p, stored["probabilities"], atol=2e-12, rtol=0)
        for bank in (actual, stored):
            n = np.asarray(bank["counts"])
            assert np.all(n >= 0) and n.sum() == m.SHOTS
            assert np.max(abs(n/m.SHOTS - p)) < radius
            assert abs(n[1].sum()/m.SHOTS - p[1].sum()) < radius

    for sign in ("-1", "1"):
        for key in ("runs", "changed_drive_contrast_zero"):
            for setting in m.SETTINGS:
                for size in m.SLICES:
                    a = replay[key][sign][setting][str(size)]
                    b = archive[key][sign][setting][str(size)]
                    check(a, b)
                    if key == "changed_drive_contrast_zero":
                        assert sum(b["counts"][1]) == 0
        for setting in ("full_turn", "half_turn_orders"):
            a = replay["runs"][sign]["material_spatial_reads"][setting]
            b = archive["runs"][sign]["material_spatial_reads"][setting]
            check(a, b)
            assert sum(b["counts"][1]) == 0
        def dark(setting):
            return sum(archive["runs"][sign][setting]["1024"]["counts"][1])/m.SHOTS
        assert dark("enclosing") - radius > 0.99
        assert dark("displaced") + radius < 0.01
        assert dark("twice") + radius < 0.02
        assert dark("parked") == 0
