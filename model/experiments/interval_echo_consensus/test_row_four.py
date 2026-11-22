from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import MultisetSector, PairSector
from model.executable_model import cell
from . import row_four
from .bounded_clock_flight import DetectorEvolution, coefficients


def test_contact_selection_precedes_release_and_exposes_doubled_load():
    costs = row_four.source_costs()
    assert costs == {"square": (0,8,32), "unary": (0,8,16), "word": (0,8,8)}
    exchange, contrast = law.local_term_parts()
    assert row_four.cold_contact(costs, exchange + 8*contrast) == {1:[8,9], 2:[32,33]}
    for wrong in (exchange, exchange-8*contrast, exchange+16*contrast):
        with pytest.raises(AssertionError):
            row_four.cold_contact(costs, wrong)
    for metric in ("unary", "word"):
        with pytest.raises(AssertionError):
            row_four.cold_contact({"square": costs[metric]}, exchange+8*contrast)


def test_inertia_against_every_eigenvalue_of_the_full_native_two_load_space():
    board = cell()
    native = MultisetSector(board, 2)
    costs = row_four.source_costs()
    bounds = row_four.spectral(board, costs)
    occupied = np.array([kind == (2,) for kind in native.histogram])
    adjacent = np.array([sum(s[a]*s[b] for a,b in board.seams) == 1 for s in native.states])
    for k, name in ((8,"attractive"),(-8,"reversed"),(0,"exchange")):
        h = native.hamiltonian(k)-12*k*np.eye(len(native.states))
        eig, vec = np.linalg.eigh(h)
        assert np.count_nonzero(eig < -1e-9) == bounds["negative_levels_"+name]
        if name == "exchange":
            continue
        which = (eig < -1e-9) if k > 0 else (eig > bounds["exchange_upper"] + 1e-9)
        assert which.sum() == len(board.seams)
        lo, hi = bounds["attractive_adjacent_band" if k > 0 else "reversed_upper_band"]
        assert lo-1e-9 <= eig[which].min() <= eig[which].max() <= hi+1e-9
        floor = bounds["attractive_eigenstate_adjacency_lower" if k>0 else "reversed_upper_eigenstate_adjacency_lower"]
        assert (abs(vec[adjacent][:,which])**2).sum(axis=0).min() >= float(Fraction(floor))-1e-10
        np.testing.assert_allclose(h[occupied][:,~occupied], 0, atol=0)
        if k < 0:
            np.testing.assert_allclose(abs(vec[~occupied][:,eig<-1e-9]), 0, atol=1e-10)


def test_fixed_point_bessel_and_flight_enclosure():
    with mp.workdps(70):
        for x in (1,16,592):
            c, error = coefficients(x)
            for n in (0,1,len(c)//3,len(c)-1):
                want = mp.besselj(n,x)*(1 if n==0 else 2)*(-1j)**(n%4)
                assert abs(want-complex(c[n])) < error
    sector = PairSector(cell())
    detector = sector.seam_states
    for k in (8,-8,0):
        h = sector.hamiltonian(k)
        solver = DetectorEvolution(h, row_four.TIMES)
        values, vectors = np.linalg.eigh(h.toarray())
        initial = detector[0]
        ports = abs((np.exp(-1j*np.outer(row_four.TIMES,values))*vectors[initial]) @ vectors[detector].T)**2
        expected = np.column_stack((ports,1-ports.sum(axis=1)))
        actual = solver.read(initial, detector)
        assert np.max(abs(actual-expected).sum(axis=1)) < solver.probability_error < 1e-5


def test_census_retains_a_failing_launch_and_rejects_missing_outcomes():
    p = np.zeros((3,20,30,6,31))
    p[0,...,0] = p[2,...,0] = 1.
    p[1,...,-1] = 1.
    p[0,17,3,5,[0,30]] = [0.,1.]
    rows = row_four.gates(p, 1e-6)["menus"]
    assert rows[0]["failing_launches"] == [[17,3]] and not rows[1]["failing_launches"]
    for invalid in (p[...,:-1], p*0, p*np.nan):
        with pytest.raises(ValueError):
            row_four.gates(invalid, 1e-6)
