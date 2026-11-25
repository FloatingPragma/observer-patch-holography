import json
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.experiments.interval_echo_einstein_station.collar_pairing.ranging.certificate import (
    certificate, integer_flight,
)
from model.experiments.interval_echo_einstein_station.collar_pairing.ranging.measure import (
    clocks, first_hits, formation, generator, kernels, ordinary_write, range_panel,
)


def test_ordinary_energy_and_clock_parities():
    ordinary_write()
    h = formation.emission.LAW.seam_term(8)
    ids = [formation.emission.LAW.SRC.index(x) for x in ((2,1),(1,2))]
    for n, recorded in enumerate(clocks(.25,24)[0]):
        u = expm(-.25j*n*h[np.ix_(ids,ids)])
        np.testing.assert_allclose([abs(u[1,0])**2,
                                    abs((u @ (np.array([1,1j])/np.sqrt(2)))[1])**2],
                                   recorded, atol=3e-14)


def test_integer_bound_and_hamiltonian_normalization():
    c = certificate()
    assert c["state_vector_error_upper"] < 3.4e-11
    board = formation.federation(1)
    route = formation.region(board)
    exact = integer_flight(board,route)
    delta = exact-2*formation.hamiltonian(board,route,1,8)
    assert delta.nnz == board.n*4
    np.testing.assert_allclose(delta.diagonal(), delta[0,0], atol=1e-13)


def test_no_click_history_against_explicit_measurements():
    board = formation.federation(1)
    route = formation.region(board)
    h = generator(board,route,1,8,True)
    u = kernels(h,route,.25,4)
    p,q=0,1
    expected = first_hits(u,p,q)
    surviving = np.zeros((h.shape[0],4), complex)
    surviving[4*route[p]:4*route[p]+4] = np.eye(4)
    target = slice(4*route[q],4*route[q]+4)
    for n in (1,2,3,4):
        surviving = expm_multiply(-.25j*h,surviving,traceA=0)
        np.testing.assert_allclose(surviving[target], expected[n], atol=3e-13)
        surviving[target] = 0


def test_static_reciprocity_and_censoring_invariants():
    archive = json.loads(Path(__file__).with_name("results.json").read_text())
    ratios = [r["triangles_in_lexical_order"]["formed"][0] for r in archive["clock_comparisons"]
              if r["step"] == .25 and r["horizon"] == 6 and not r["left_odd"]]
    assert ratios == pytest.approx([.74915630246,1.30429335216,1.],abs=1e-9)
    for run in archive["runs"]:
        for panel in run["panels"]:
            for pair in panel["pairs"]:
                for data in pair["arms"].values():
                    for row in data["horizons"]:
                        assert row["complete"]+row["waiting_for_return"]+row["missed_outward"] == pytest.approx(1.)
                        assert 0 < row["range"] <= row["horizon"]/2+1e-13
