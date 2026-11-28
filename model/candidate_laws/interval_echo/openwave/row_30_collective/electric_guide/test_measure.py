from itertools import combinations

import numpy as np
import pytest
from scipy.linalg import expm

from . import connector_guard, local_signal, measure, release


def test_physical_dimer_words_produce_cross_site_signs():
    board = measure.federation(2)
    _, path, _ = measure.electric.guide(board)
    matrix, _, good = release.apparatus(board, path)
    physical = matrix[np.ix_(good, good)]
    words = [sum(1<<i for i in bits) for bits in combinations(range(8), 4)]
    words = list(np.array(words)[good])
    signs = []
    for site in range(4):
        op = np.zeros(physical.shape, dtype=complex)
        for col, word in enumerate(words):
            flipped = word ^ ((1 << (2*site+2)) - 1)
            y = 1j if word & (1 << 2*site) else -1j
            op[words.index(flipped), col] = (-1)**site*y
        signs.append(op)
    spectra, vectors = measure.modes(4, *measure.electric.couplings(8))
    for t in (.125, .54, 2.25, 16):
        prop = expm(-1j*t*physical)
        initial = prop @ signs[0] @ prop.conj().T
        measured = np.array([np.trace(op @ initial).real/16 for op in signs])
        np.testing.assert_allclose(measured, vectors @ (np.cos(spectra*t)*vectors[0]),
                                   atol=2e-11, rtol=0)


def test_spatial_fit_holds_out_third_addresses():
    x = np.arange(1, 45)
    for q in (.075, .51, 1.5, 2.9):
        prof = -.23*np.sin(q*x)
        fitted, residual = measure.wavelength(prof)
        assert abs(fitted-q) < 1e-7 and residual < 1e-6
        prof[x % 3 == 0] *= -1
        challenged, residual = measure.wavelength(prof)
        assert abs(challenged-fitted) < 1e-12 and residual > .05


def test_released_stock_changes_the_actual_signal():
    board = measure.federation(2)
    _, path, _ = measure.electric.guide(board)
    p = release.probabilities(board,path)
    assert np.allclose(p[:,:,1:].sum(axis=2), 1)
    assert np.allclose(p[0,:,0],1)
    assert .99 < p[1,-1,0] < .995
    assert abs((p[0,-1,3]-p[0,-1,1])-(p[1,-1,3]-p[1,-1,1])) > .18
    with np.load(measure.HERE/"counts.npz", allow_pickle=False) as source:
        bank = source["release"]
    assert np.max(abs(bank/measure.SHOTS-p))<.003
    bank[0,0,1] += 1
    with pytest.raises(ValueError, match="Incomplete guard"):
        release.read(bank, measure.SHOTS)


def test_coherent_connector_memory_retains_the_blank_leakage():
    from model.candidate_laws.contrast_memory_spherical_v1 import law

    assert connector_guard.exact_columns()["coefficient_entries_checked"] == 8712
    memory = np.arange(1,12)*np.exp(.7j*np.arange(11))
    memory /= np.linalg.norm(memory)
    p = abs(memory[0])**2
    h0, h1 = law.local_term_parts()
    H = h0+8*h1
    t = 1.37
    unitary = expm(-1j*t*H)
    mean = []
    for word in [(0,0), (2,0)]:
        col = np.kron(memory, np.eye(9)[:,law.SRC.index(word)])
        mean.append(np.vdot(col,H@col).real)
        if word==(2,0):
            evolved = (unitary@col).reshape(11,9)
            swapped = evolved[:,law.SRC.index((0,2))]
            assert np.vdot(swapped,swapped).real == pytest.approx(p*np.sin(t/2)**2, abs=1e-14)
    assert (mean[1]-mean[0])/2 == pytest.approx(16.25*p,abs=1e-12)


def test_local_dimer_return_from_physical_words():
    board = measure.federation(2)
    _, path, _ = measure.electric.guide(board)
    full, _, good = release.apparatus(board,path)
    matrix = full[np.ix_(good, good)]
    words = np.array([sum(1<<i for i in bits) for bits in combinations(range(8),4)])[good].tolist()
    x = np.zeros((4,16,16))
    for site in range(4):
        for col, word in enumerate(words):
            x[site,words.index(word ^ (3 << (2*site))),col] = 1
    j, h = measure.electric.couplings(8)
    formula = local_signal.response(4,1,j,h,"native")
    for t in (.25, 2, 16):
        prop = expm(-1j*t*matrix)
        state = prop@x[1]@prop.conj().T
        expected = [np.trace(a@state).real/16 for a in x]
        np.testing.assert_allclose(formula[local_signal.HOLDS.index(t)], expected, atol=2e-11)
