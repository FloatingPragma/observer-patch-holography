import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell
from .__main__ import evolve, perturbed_amplitudes, verify_local_reader


def test_release_matches_full_cell_exponential_at_late_times():
    sector = PairSector(cell())
    source = np.zeros(sector.P)
    source[sector.index(0, 1)] = 1
    for kappa in (0., candidate.KAPPA):
        h = sector.hamiltonian(kappa)
        times = (0, 7, 37)
        for t, v in zip(times, evolve(h, kappa, source, times), strict=True):
            np.testing.assert_allclose(v, expm(-1j * t * h.toarray()) @ source, atol=1e-11)


def test_actual_seam_measurement_and_preparation():
    verify_local_reader()
    a, b = perturbed_amplitudes()
    np.testing.assert_allclose(b / a, 1j * np.sqrt(.04/.96), atol=1e-12)
