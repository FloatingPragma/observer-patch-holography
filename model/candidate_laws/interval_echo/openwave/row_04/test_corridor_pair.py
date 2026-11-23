import numpy as np
import pytest
from scipy.linalg import expm

from model.executable_model import federation
from model.executable_model.carrier import HERE as CARRIER_HOME
from ..submission.row_04_check import screen
from . import corridor_pair as corridor


def test_collar_preparation_and_original_patch_outcomes():
    setup = corridor.exact_flight()
    cut, home = corridor.geometry()
    assert set(setup["neutral_edges"]) == set(cut) and len(cut) == 12
    assert len(setup["live_swaps"]) == 4 and setup["final_pair"] == (156, 157)
    edges, index, _, adjacent = screen(CARRIER_HOME/"carriers"/"FEDERATION_LEVEL0.npz")
    addresses = [tuple(pair) for pair in edges if pair[0]//12 == pair[1]//12 == 13]
    assert addresses == list(home)
    for column, state in enumerate(corridor.BASIS):
        bins = corridor.receiver_bins(np.eye(3)[:, column])
        expected = addresses.index(state) if adjacent[index[state]] else 30
        np.testing.assert_array_equal(bins, np.eye(31)[expected])
    assert setup["dark_terminal"][0] == setup["dark_terminal"][7] == .5
    assert sum(setup["dark_terminal"]) == 1


@pytest.mark.parametrize("kappa", [-3, 0, .5, 8])
def test_complete_native_columns_before_restriction(kappa):
    # Assemble from full 99-state seam columns, with arbitrary residual load
    # pairs kept in the output dictionary. This detects escapes discarded by a
    # premature compression to three input states.
    board = federation(0)
    law = corridor.LAW
    blank = np.eye(law.DIM)[:, :law.N_SRC]
    h0, h1 = law.local_term_parts()
    h = h0+kappa*h1
    cut = set(corridor.geometry()[0])
    full = {}
    for col, pair in enumerate(corridor.BASIS):
        for edge in board.seams:
            word = tuple(int(v in pair) for v in edge)
            word_index = law.SRC.index(word)
            if edge in cut:
                np.testing.assert_allclose(h @ corridor.neutral_column(word), 0, atol=1e-13, rtol=0)
                continue
            out = h @ blank[:, word_index]
            np.testing.assert_allclose(out[law.N_SRC:], 0, atol=1e-13, rtol=0)
            for j in np.flatnonzero(np.abs(out[:law.N_SRC]) > 1e-13):
                destination = list(v for v in pair if v not in edge)
                for v, multiplicity in zip(edge, law.SRC[j], strict=True):
                    destination.extend([v]*multiplicity)
                key = tuple(sorted(destination)), col
                full[key] = full.get(key, 0.)+out[j]
        key = pair, col
        full[key] = full.get(key, 0.) + kappa*sum(6-board.degree[v] for v in pair)
    assert all(pair in corridor.BASIS for pair, _ in full)
    reduced = [[full.get((pair, i), 0.) for i in range(3)] for pair in corridor.BASIS]
    np.testing.assert_allclose(reduced, corridor.corridor_generator(kappa), atol=1e-12, rtol=0)
    dark = np.array([1, 0, -1])/np.sqrt(2)
    at = 2.31
    phase = np.exp(-1j*(.5+2*kappa)*at)
    np.testing.assert_allclose(expm(-1j*at*np.array(reduced)) @ dark, phase*dark, atol=5e-13)


def test_broken_cut_and_phase_control(monkeypatch):
    cut, home = corridor.geometry()
    exchange, _ = corridor.corridor_parts()
    monkeypatch.setattr(corridor, "geometry", lambda: (cut[1:], home))
    with pytest.raises(ArithmeticError, match="outside the path"):
        corridor.corridor_parts.__wrapped__()
    bright = np.array([1, 0, 1])/np.sqrt(2)
    at = 2*np.pi/3
    # At zero coupling the bright state reaches quiet probability 8/9.
    probability = abs((expm(-1j*at*exchange) @ bright)[1])**2
    assert abs(probability-8/9) < 2e-14


@pytest.mark.parametrize("invalid", [[-1., 1., 1.], [float("nan"), 0., 0.], [1., 1., 0.], [1.]])
def test_terminal_law_rejects_invalid_weights(invalid):
    with pytest.raises(ValueError, match="material weights"):
        corridor.receiver_bins(invalid)
