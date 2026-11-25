import json
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.distance import cdist

from model.experiments.interval_echo_einstein_station.collar_pairing.ranging.pool import (
    correction_at_sixth_order, distance_check, generator_difference, leading_reciprocity, pooled_ranges,
)
from model.experiments.interval_echo_einstein_station.collar_pairing.ranging.measure import (
    formation, generator, kernels, prepared_states,
)


def test_free_gains_and_tetrahedron():
    xyz = np.array([[0,0,0],[2,0,0],[0,3,0],[0,0,4]])
    gain = np.array([.5,1.,1.,2.])
    d = cdist(xyz,xyz)
    launch = gain[:,None]*d
    check = distance_check(launch)
    np.testing.assert_allclose(check["cycle_ratios"],1,atol=1e-13)
    np.testing.assert_allclose(check["gains"],gain)
    np.testing.assert_allclose(check["projected"],d)
    assert min(check["centered_eigenvalues"])>0
    launch[1,2] *= 1.03
    assert distance_check(launch)["necessary_log_error_floor"] > .004
    cycle = np.array([[0,1,2,1],[1,0,1,2],[2,1,0,1],[1,2,1,0.]])
    assert min(distance_check(cycle)["centered_eigenvalues"]) < 0


def test_species_erasures_before_range_calculation():
    archive = json.loads(Path(__file__).with_name("results.json").read_text())
    sources = [r for r in archive["runs"] if r["contrast"] == 8 and not r["left_odd"]]
    for step,panels in enumerate(zip(*(r["panels"] for r in sources),strict=True)):
        for k,horizon in enumerate((2,4,6)):
            for arm in ("formed","independent"):
                distances = pooled_ranges(panels,arm,horizon)
                for x,y in zip(panels[0]["pairs"],panels[1]["pairs"],strict=True):
                    average = sum(p["arms"][arm]["horizons"][k]["range"] for p in (x,y))/2
                    assert distances[x["p"],x["q"]] == pytest.approx(average,abs=1e-13)
                if step==0 and horizon==6 and arm=="formed":
                    assert distance_check(distances)["cycle_ratios"][2] == pytest.approx(.89650832734)


def test_leading_cancellation_and_sixth_order_reversal():
    assert leading_reciprocity()==12
    assert generator_difference()==32
    correction_at_sixth_order()
    board = formation.federation(1)
    route = formation.region(board)
    rho = prepared_states(board,route,1)[0]
    polarization = -np.trace(np.diag([1,-1,1,-1])@rho).real
    t = 1/256
    difference = 0.
    for q in (-1,1):
        u = kernels(generator(board,route,q,8,False),route,t,1)[1]
        for p,r,sign in ((2,3,1),(3,2,-1)):
            turn = u[p,r]@u[r,p]
            difference += sign*np.trace(turn.conj().T@turn@rho).real/2
    assert difference/t**6 == pytest.approx(2048*polarization/3,rel=.002)
