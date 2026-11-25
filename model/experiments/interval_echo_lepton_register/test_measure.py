import json
from fractions import Fraction
from itertools import combinations

import numpy as np
from scipy.linalg import null_space
from sympy import Matrix

from model.executable_model.carrier import federation
from . import measure as m


def test_ideal_response_and_replayed_complete_port_banks():
    actual = m.measure()
    archive = json.loads((m.HERE / "results.json").read_text())
    angle = Fraction(1, 12) + Fraction(1, 6) - 2 * Fraction(1, 72)
    np.testing.assert_allclose(
        actual["branch_probabilities"], [8 / 12, 2 / 12, 1 / 12, 1 / 12]
    )
    np.testing.assert_allclose(
        [actual["independent_xor"], actual["shared_address_xor"]], [float(angle), 1 / 4]
    )
    assert actual["paired_port_counts"] == archive["paired_port_counts"]
    assert int(np.sum(actual["paired_port_counts"])) == m.SHOTS
    roots = (np.sqrt(2) + 2 * np.cos(float(angle) + 2 * np.pi * np.arange(3) / 3)) / 6
    assert min(roots) > 0
    for edge in actual["all_internal_seams"]:
        np.testing.assert_allclose(
            edge["character_probabilities"], roots**2, atol=3e-13
        )
    assert len(actual["all_internal_seams"]) == 30
    for bank, saved in zip(actual["banks"], archive["banks"], strict=True):
        np.testing.assert_allclose(
            bank["all_port_probabilities"], saved["all_port_probabilities"], atol=2e-12
        )
        assert bank["counts"] == saved["counts"] and sum(bank["counts"]) == m.SHOTS
    base, zero, flipped, incoherent, changed, measured = actual["banks"]
    raw = np.array(actual["paired_port_counts"])
    a, b, plus, _ = actual["ports"]
    observed = (
        raw[plus].sum() + raw[:, [a, b]].sum() - 2 * raw[plus, [a, b]].sum()
    ) / m.SHOTS
    assert observed == actual["measured_xor"] == measured["tip_phase_radians"]
    empirical = (np.sqrt(2) + 2 * np.cos(observed + 2 * np.pi * np.arange(3) / 3)) / 6
    np.testing.assert_allclose(measured["character_probabilities"], empirical**2)
    np.testing.assert_allclose(base["intensity_quotient"], 2 / 3)
    np.testing.assert_allclose(
        flipped["character_probabilities"], roots[[0, 2, 1]] ** 2
    )
    np.testing.assert_allclose(
        zero["character_probabilities"][1:], zero["character_probabilities"][1]
    )
    np.testing.assert_allclose(incoherent["character_probabilities"], [1 / 12] * 3)
    assert abs(changed["intensity_quotient"] - 2 / 3) > 0.15


def test_leakage_from_exact_adjacency_and_independent_generator():
    board = federation(0)
    for edge in m.CELL_SEAMS:
        a, b, plus, minus = m.addresses(edge)
        outside = sorted(set(range(m.N)) - {a, b, plus, minus})
        rows = [
            [
                sum(q in board.neighbours[p] for q in col)
                for col in [(a, b), (plus,), (minus,)]
            ]
            for p in outside
        ]
        # An integer unit minor certifies full rank without a numerical cutoff.
        assert any(
            abs(Matrix([rows[i] for i in ix]).det()) == 1
            for ix in combinations(range(len(rows)), 3)
        )
    a, b, plus, minus = m.addresses(m.CELL_SEAMS[0])
    cols = np.zeros((board.n, 3))
    cols[[a, b], 0] = 1 / np.sqrt(2)
    cols[plus, 1] = cols[minus, 2] = 1
    h = np.diag(m.KAPPA * (2 * board.degree - 6) + board.degree / 2)
    for u, v in board.seams:
        h[u, v] = h[v, u] = -0.5
    outside_modes = null_space(cols.T)
    leakage = outside_modes.T @ h @ cols
    archive = json.loads((m.HERE / "results.json").read_text())["release"]
    np.testing.assert_allclose(
        np.linalg.eigvalsh(leakage.T @ leakage), archive["leakage_squared_eigenvalues"]
    )
    freq, vec = np.linalg.eigh(h)
    analyzer = m.reader((a, b, plus, minus), 2 / 9, 1)
    source = np.zeros((board.n, 3), complex)
    source[: m.N] = analyzer.conj().T[:, [a, minus, plus]]
    for run in archive["records"]:
        state = vec @ (np.exp(-1j * freq * run["time"])[:, None] * (vec.T @ source))
        prob = abs(analyzer @ state[: m.N]) ** 2
        complete = np.column_stack((prob.T, 1 - prob.sum(axis=0)))
        counts = np.array(run["home_and_no_vacancy_counts"])
        assert np.all(counts.sum(axis=1) == m.SHOTS)
        np.testing.assert_allclose(counts / m.SHOTS, complete, atol=0.0015)
        np.testing.assert_allclose(prob.T, run["home_port_probabilities"], atol=2e-12)
        np.testing.assert_allclose(
            1 - prob.sum(axis=0), run["outside_home"], atol=2e-12
        )
