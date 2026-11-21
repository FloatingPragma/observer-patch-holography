from itertools import combinations
import json

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from .experiment import ANGLES, HERE, WEIGHTS, hardware
from .interval import SHOTS


def test_native_analyzer_and_gate_response_frequency():
    exchange, contrast = native_law.local_term_parts()
    h = exchange + KAPPA * contrast
    analyzer = expm(-.5j * np.pi * h)
    ids = [native_law.SRC.index(q) for q in ((1, 0), (0, 1))]
    phased = np.eye(len(h), dtype=complex)[:, ids] @ np.diag([1, -1j])
    minus = (analyzer @ phased)[ids[1]]
    np.testing.assert_allclose(np.outer(minus.conj(), minus), np.array([[1, -1], [-1, 1]]) / 2, atol=1e-12)
    for frequency in range(-9, 10):
        np.testing.assert_allclose(WEIGHTS @ np.exp(1j * frequency * ANGLES), 1j * frequency, atol=1e-12)


def test_balance_across_the_full_pair_sector():
    sector, _, _, meter, cells, bonds = hardware()
    h = sum(cells) + sum(b.h for b in bonds)
    region = cells[0]
    def bracket(other):
        return other @ region - region @ other
    residual = bracket(h) - sum(bracket(b.h) for b in bonds)
    assert residual.nnz == 0
    assert bracket(bonds[1].h).nnz > 0
    for a, b in combinations(cells, 2):
        assert (a @ b - b @ a).nnz == 0
    rng = np.random.default_rng(42)
    psi = rng.normal(size=sector.P) + 1j * rng.normal(size=sector.P)
    psi /= np.linalg.norm(psi)
    np.testing.assert_allclose(meter.read(psi), np.vdot(psi, region @ psi), atol=1e-12)
    for bond in bonds:
        read = [meter.read(bond.expose(psi, angle)) for angle in ANGLES]
        np.testing.assert_allclose(WEIGHTS @ read, 1j * np.vdot(psi, bracket(bond.h) @ psi), atol=1e-12)


def test_finite_counts_and_withheld_gate():
    result = json.loads((HERE / "interval.json").read_text())
    with np.load(HERE / "counts.npz") as archive:
        words = archive["diagonal_word_counts"]
        assert words.shape == (230, 79) and np.all(words.sum(axis=1) == SHOTS)
        clicks = archive["antisymmetric_click_counts"]
        assert clicks.shape == (230, 30) and np.all((clicks >= 0) & (clicks <= SHOTS))
    assert abs(result["expected_residual_6"]) < 1e-6
    assert result["endpoint_change"]["excludes_zero_at_4se"]
    assert all(g["excludes_zero_at_4se"] for g in result["gateways"])
    assert not result["endpoint_minus_both"]["excludes_zero_at_4se"]
    assert result["endpoint_minus_first"]["excludes_zero_at_4se"]
