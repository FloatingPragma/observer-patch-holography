"""Independent dense density-matrix checks for the readout and dressed-band calculations."""
import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model import carrier

from .__main__ import bound_project, evaluate_arm


def test_cell_pointer_and_spectral_diagnostics_match_dense_ensemble():
    sector = PairSector(carrier.cell())
    h = sector.hamiltonian(8)
    eigenvalues, eigenvectors = eigh(h.toarray())
    band = eigenvectors[:, eigenvalues < -5]
    exact = band @ band.T
    approximate, error = bound_project(h, np.eye(sector.P))
    assert error < 1e-10
    assert np.linalg.norm(approximate - exact, ord=2) < error + 1e-13

    def unitary(t):
        return (eigenvectors * np.exp(-1j * t * eigenvalues)) @ eigenvectors.T.conj()

    internal = np.flatnonzero(sector.adjacent)
    target = internal[0]
    starts = np.eye(sector.P)[:, ~sector.adjacent]
    rho = (starts @ starts.T) / starts.shape[1]
    early = unitary(1) @ rho @ unitary(1).conj().T
    later = unitary(4)
    p = early[target, target].real
    quiet = early.copy()
    quiet[target, :] = 0
    quiet[:, target] = 0
    sharp = quiet.copy()
    sharp[target, target] = p
    repeat = np.sum(np.abs(later[internal, target])**2)
    result = evaluate_arm(sector, 8, starts, target, internal,
                          {"first_positive_per_trial": p, "repeat_positive_given_first": repeat})
    for name, ensemble in (("bypass", early), ("unread_sharp", sharp), ("half_serviced", (sharp + early)/2)):
        expected = np.real(np.diag(later @ ensemble @ later.conj().T))[internal]
        observed = result["complete_patch_distributions"][name]
        np.testing.assert_allclose(observed, np.append(expected, 1 - expected.sum()), atol=1e-12, rtol=0)
    diagnostic = result["global_dressed_band_diagnostic"]
    np.testing.assert_allclose(diagnostic["source_and_bypass"], np.trace(exact @ early).real, atol=1e-12, rtol=0)
    np.testing.assert_allclose(diagnostic["unread_sharp"], np.trace(exact @ sharp).real, atol=1e-12, rtol=0)
    np.testing.assert_allclose(diagnostic["conditional_sharp_quiet"],
                               np.trace(exact @ quiet).real / (1 - p), atol=1e-12, rtol=0)
