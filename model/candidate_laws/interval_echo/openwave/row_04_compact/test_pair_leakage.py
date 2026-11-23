from fractions import Fraction

import numpy as np
import pytest
import sympy as sp
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model import cell, federation
from .pair_leakage import certify_pair_leakage, check_internal_spectrum, exit_block


def test_all_homes_and_the_sharp_m103_launch():
    assert certify_pair_leakage() == ((20, 80, 320), Fraction(1, 4), Fraction(5, 4))


def test_the_certificate_rejects_a_spurious_dark_pair():
    local = cell()
    block = sp.Matrix(list(exit_block(local, local.seams).values()))
    assert check_internal_spectrum(block) == Fraction(1, 4)
    block[:, 0] = sp.zeros(30, 1)
    with pytest.raises(ArithmeticError, match="spectrum"):
        check_internal_spectrum(block)


def test_quiet_detector_reads_the_same_initial_curvature_for_either_sign():
    board = federation(0)
    sector = PairSector(board)
    quiet = np.ones(sector.P, dtype=bool)
    quiet[[sector.index(a, b) for a, b in board.seams if a//12 == b//12 == 13]] = False
    state = np.zeros(sector.P)
    state[sector.index(156, 157)] = 1
    t = 0.0005
    for kappa in (-8, 0, 8):
        flight = expm_multiply((-1j*t)*sector.hamiltonian(kappa), state)
        loss = np.linalg.norm(flight[quiet])**2
        assert abs(loss/t**2 - 1.25) < 0.0001
