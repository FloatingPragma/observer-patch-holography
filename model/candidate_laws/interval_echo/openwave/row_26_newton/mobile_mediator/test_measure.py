"""Independent checks of the clamped native columns, moments and observable bands."""

from fractions import Fraction
import json
from itertools import product

import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import eigvalsh

from . import matched_energy, measure


@pytest.fixture(scope="module")
def saved():
    return json.loads((measure.HERE/"results.json").read_text())


def test_direct_full_board_native_action():
    board = measure.federation(0)
    _, _, pins = next(x for x in measure.sources(board) if x[0] == 2)
    law = measure.candidate.native_law
    allowed = np.setdiff1d(np.arange(board.n), pins)
    seam = law.seam_term(8)
    edges = np.array(board.seams)
    def full_energy(load):
        ids = [law.SRC.index((load[a],load[b])) for a,b in edges]
        return sum(seam[ids,ids]) + 8*np.dot(6-board.degree,load**2)
    blank = np.ones(board.n, int)
    for _, _, sites in measure.sources(board):
        for orientation, expected in ((1,133),(-1,101)):
            target = blank.copy()
            target[list(sites[2:])] += [-orientation, orientation]
            assert full_energy(target) - full_energy(blank) == expected
    for q, reverse in product((-1,1),(False,True)):
        graph, h, va, vb = measure.operators(board,pins,8,q,reverse)
        for i,j in measure.SETTINGS:
            host = np.ones(board.n, int)
            host[list(pins)] += np.array([-i, i, j*(2*reverse-1), j*(1-2*reverse)])
            e0 = full_energy(host)
            for loc in range(len(allowed)):
                test = host.copy()
                test[allowed[loc]] += q
                assert 2*(full_energy(test)-e0) == h[loc,loc] + i*va[loc] + j*vb[loc]
            # The exact local swap produces all unconstrained off-diagonals.
            swap = seam[law.SRC.index((1,1+q)), law.SRC.index((1+q,1))]
            assert np.array_equal(h - np.diag(h.diagonal()), 2*swap*graph)


def test_return_order_from_independent_native_powers_modulo_two_to_64(saved):
    board = measure.federation(0)
    for panel in saved["panels"]:
        if panel["kappa"] != 8 or panel["mobile_charge"] != 1:
            continue
        graph, h, va, vb = measure.operators(board,panel["pins"],8,1,panel["reversed_receiver"])
        exact = measure.first_mixed_moment(graph,va,vb)
        assert exact == panel["first_mixed"]
        order = exact["order"]
        sums = [0]*(order+1)
        # Integer matrix operations in Z/(2^64) give exact trace residues
        # even when individual powers exceed the machine word.
        for sign,(i,j) in zip(measure.SIGNS, measure.SETTINGS, strict=True):
            op = sparse.csr_matrix((h + np.diag(i*va+j*vb)).astype(np.uint64))
            power = np.eye(len(va),dtype=np.uint64)
            for t in range(1, order+1):
                power = op @ power
                sums[t] += int(sign)*int(np.trace(power))
        modulus = 2**64
        assert all(s % modulus == 0 for s in sums[:order])
        expected = Fraction(exact["trace_moment"]) * 2**order
        assert expected.denominator == 1 and sums[order] % modulus == int(expected) % modulus


def test_far_thermal_signs_ground_controls_and_sampled_returns(saved):
    radius = saved["joint_99pct_connected_ramsey_radius"]
    assert len(saved["panels"]) == 40
    for q, reverse in product((-1,1),(False,True)):
        readings = [p for p in saved["panels"] if p["kappa"] == 8 and p["distance"] >= 2
                    and p["mobile_charge"] == q and p["reversed_receiver"] == reverse]
        assert [p["first_mixed"]["order"] for p in readings] == [6,12,18,22]
        potentials = [p["potential_per_activity_beta_eighth_z1e_minus55"] for p in readings]
        if reverse:
            assert all(a[0] > b[1] > 0 for a,b in zip(potentials,potentials[1:]))
        else:
            assert all(a[1] < b[0] < 0 for a,b in zip(potentials,potentials[1:]))
    for p in saved["panels"]:
        for read in p["ramsey"]:
            assert all(0 <= k <= measure.SHOTS for k in read["clicks"])
            assert 2*(measure.SIGNS @ read["clicks"])/measure.SHOTS == read["measured"]
            assert abs(read["predicted"] - read["measured"]) < radius
        if p["kappa"] == 0:
            assert p["heat_trace_beta_eighth"][0] > 0
    near = next(p for p in saved["panels"] if p["kappa"]==8 and p["distance"]==2
                and p["mobile_charge"]==1 and not p["reversed_receiver"])
    x = next(r for r in near["ramsey"] if r["wait"]==8 and r["component"]=="real")
    assert x["measured"] + radius < -.008
    assert near["ground_error_radius"] < 1e-8 and near["frozen_connected_ground"] == 17


def test_thermal_replay_and_independent_trace_at_nearest_noncontact_shell(saved):
    board = measure.federation(0)
    for p in saved["panels"]:
        if p["kappa"] != 8 or p["mobile_charge"] != 1:
            continue
        graph,h,va,vb = measure.operators(board,p["pins"],8,1,p["reversed_receiver"])
        box = measure.thermal_trace(h,va,vb)
        a,b = p["heat_trace_beta_eighth"]
        assert max(a,box[0]) <= min(b,box[1])
        assert np.all(va*vb == 0)
        if p["distance"] == 2:
            # Direct subtraction is numerically resolved on this shell.
            levels = np.array([eigvalsh((h+np.diag(i*va+j*vb))/2) for i,j in measure.SETTINGS])
            raw = measure.SIGNS @ np.exp(-levels/8).sum(axis=1)
            assert a-2e-13 < raw < b+2e-13


def test_equal_static_energy_and_opposite_interaction_across_four_separations():
    board = measure.federation(0)
    archive = json.loads((measure.HERE/"matched_energy.json").read_text())
    assert len(archive["panels"]) == 20
    for q, reverse in product((-1,1),(False,True)):
        selected = [p for p in archive["panels"] if p["distance"] >= 2
                    and p["mobile_charge"] == q and p["reversed_receiver"] == reverse]
        assert [p["first_mixed"]["order"] for p in selected] == [10,16,22,26]
        bounds = [p["potential_per_activity_beta_eighth_z1e_minus72"] for p in selected]
        for p in selected:
            assert list(board.degree[p["pins"]]) == [5]*4
            for a,b in (p["pins"][:2],p["pins"][2:]):
                assert b in board.neighbours[a]
                # Each adjacent degree-five dipole has eight unit-contrast
                # seams and one contrast-two seam plus two counting-reference units.
                assert matched_energy.source_energy(board,(a,b),False) == 8*12 + 9/2 + 8*2
                assert matched_energy.source_energy(board,(a,b),True) == 116.5
            assert p["first_source_energy"] == p["second_source_energy"] == 116.5
            graph,h,va,vb = measure.operators(board,p["pins"],8,q,reverse)
            assert measure.first_mixed_moment(graph,va,vb) == p["first_mixed"]
            if p["distance"] == 5 and q == 1:
                lo, hi = measure.thermal_trace(h,va,vb)
                a,b = p["heat_trace_beta_eighth"]
                assert max(lo,a) <= min(b,hi)
        if reverse:
            assert all(a[0] > b[1] > 0 for a,b in zip(bounds,bounds[1:]))
        else:
            assert all(a[1] < b[0] < 0 for a,b in zip(bounds,bounds[1:]))
