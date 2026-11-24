import sympy as sp

from model.executable_model import cell, federation
from .active_current import EDGE, exact_modes, leading_coefficients


def test_exterior_coefficients_against_full_cell_and_bridge_paths():
    board = federation(0)
    bridges = [(a, b) for a, b in board.external if a < 12]
    u, v = exact_modes(cell(), [a for a, _ in bridges])
    coefficients = leading_coefficients(board, bridges, u, v)
    # Use the full native kinetic operator on the cell and its first exterior receivers.
    # Higher exterior paths cannot enter its second power between these addresses.
    sites = [*range(12), *(b for _, b in bridges)]
    lap = sp.Matrix(board.laplacian()[sites, :][:, sites].toarray().astype(int)) / 2
    even_phase = (-1 - sp.I * sp.sqrt(3)) / 2
    odd_phase = -1
    amplitudes = []
    for sign in (1, -1):
        wave = (u + sign * sp.I * v) / sp.sqrt(2)
        left, right = (wave[p] for p in EDGE)
        even = (left + right) * even_phase / 2
        odd = (left - right) * odd_phase / 2
        wave[EDGE[0]], wave[EDGE[1]] = even + odd, even - odd
        initial = wave.col_join(sp.zeros(len(bridges), 1))
        assert (lap * initial)[12:, :] == sp.zeros(len(bridges), 1)
        amplitudes.append((-lap**2 * initial / 2)[12:, :])
    signed = [sp.simplify(abs(a)**2 - abs(b)**2) for a, b in zip(*amplitudes, strict=True)]
    assert signed == coefficients == [-sp.Rational(1, 128), sp.Rational(1, 128), 0]
