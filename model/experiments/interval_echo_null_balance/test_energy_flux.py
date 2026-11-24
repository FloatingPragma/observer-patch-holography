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


def test_signed_cut_and_refusing_shared_storage_on_independent_single_load():
    from scipy import sparse
    from model.executable_model import federation

    board = federation(0)
    flight = board.laplacian()/2

    def term(a, b):
        return sparse.csr_matrix(([KAPPA + .5, KAPPA + .5, -.5, -.5],
                                  ([a, b, a, b], [a, b, b, a])), shape=flight.shape)

    def bracket(a, b):
        return a @ b - b @ a

    inside = [sparse.diags(KAPPA * (6 - board.degree) * (np.arange(board.n)//12 == cell))
              for cell in range(board.cells)]
    links = {}
    for a, b in board.seams:
        l, r = a//12, b//12
        if l == r:
            inside[l] = inside[l] + term(a, b)
        else:
            links[(l, r)] = term(a, b)

    owners = (0, 6)
    core = sum(inside[c] for c in owners)
    allocated, incoming = core.copy(), sparse.csr_matrix(flight.shape)
    for (l, r), bond in links.items():
        allocated += ((l in owners) + (r in owners)) * .5 * bond
        orientation = int(l in owners) - int(r in owners)
        incoming += .5 * orientation * bracket(bond, inside[l] - inside[r])
    assert (bracket(flight, allocated) - incoming).count_nonzero() == 0
    assert (bracket(flight, allocated) + incoming).count_nonzero() > 0
    assert (bracket(flight, core) - incoming).count_nonzero() > 0
    missing = np.isin(np.arange(board.n)//12, owners) * (6-board.degree) * KAPPA
    assert (bracket(flight, allocated-sparse.diags(missing)) - incoming).count_nonzero() > 0


def test_hidden_seam_currents_are_divergence_free():
    import sympy as sp
    from . import hidden_current
    from model.executable_model import cell

    modes = hidden_current.exact_modes(cell(), [1, 3, 10])
    balance = [sp.S.Zero] * 12
    for seam in hidden_current.seam_laws(cell(), *modes):
        a, b = seam['ports']
        flux = sp.sympify(seam['kinetic_current'])
        if (a, b) == (4, 5):
            assert abs(flux) == sp.sqrt(3)/4
        balance[a] -= flux
        balance[b] += flux
    assert all(sp.simplify(value) == 0 for value in balance)
