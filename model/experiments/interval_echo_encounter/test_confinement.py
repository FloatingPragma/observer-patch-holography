"""Independent global-sector tests of local closure, sign read and controls."""
from copy import deepcopy
from itertools import product
import json

import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model.carrier import federation
from . import confinement, conjugation, sign_read, transport


def test_exact_spaces_and_cell_counts():
    actual = confinement.certificate()
    assert actual == json.loads(confinement.RECEIPT.read_text())
    assert actual["distinct_cell_product_dimensions"] == {"0": 760, "1": 26500, "2": 451420}
    assert actual["quiet_quadratic_coefficients"] == {"home": "1/4", "seams": "5/4"}
    for p in actual["profiles"]:
        for window in ("home", "seams"):
            assert [case["dimension"] for case in p["windows"][window].values()] == [0]*3
    # An isolated cell has a wholly invariant home space. A solver always
    # returning zero would pass every connected two-load negative case.
    _, matrix, outward = confinement.blocks((), "home", 8)
    _, basis = confinement.maximal_space(matrix, outward.values())
    np.testing.assert_array_equal(basis, np.eye(66))


def test_global_sector_includes_every_exit_and_diagonal():
    car = federation(0)
    pair = PairSector(car)
    h = 2 * pair.hamiltonian(8)
    cell = 13
    gates = confinement.gate_profiles(car)[cell]
    outside = [next(n for n in car.neighbours[12*cell+g] if n//12 != cell) for g in gates]
    address = list(range(12*cell, 12*cell+12)) + outside
    for window in ["home", "seams"]:
        states, inside, exits = confinement.blocks(gates, window, 8)
        index = [pair.index(address[a], address[b]) for a,b in states]
        columns = h[:, index]
        np.testing.assert_array_equal(columns[index].toarray(), inside)
        for (a,b), row in exits.items():
            np.testing.assert_array_equal(columns[pair.index(address[a],address[b])].toarray()[0], row)
        assert columns.count_nonzero() == np.count_nonzero(inside) + np.count_nonzero(list(exits.values()))


def test_maximality_by_independent_modular_observability_gramian():
    # A nonzero minor modulo a prime also survives over the rationals.
    # The rational invariant witnesses supply the opposite dimension bound.
    prime = 1009
    for profile in json.loads(confinement.RECEIPT.read_text())["profiles"]:
        for window, cases in profile["windows"].items():
            for strength, case in cases.items():
                _, a, exits = confinement.blocks(profile["gates"], window, int(strength))
                rows = np.array(list(exits.values()))
                gram = np.zeros_like(a)
                for _ in range(len(a)):
                    gram = (gram + rows.T @ rows) % prime
                    rows = (rows @ a) % prime
                rank = 0
                for col in range(len(a)):
                    pivots = np.flatnonzero(gram[rank:,col])
                    if not pivots.size:
                        continue
                    pivot = rank + pivots[0]
                    gram[[rank,pivot]] = gram[[pivot,rank]]
                    factor = pow(int(gram[rank,col]), -1, prime)
                    gram[rank] = (gram[rank]*factor) % prime
                    gram[rank+1:] = (gram[rank+1:] -
                                     gram[rank+1:,col,None]*gram[rank]) % prime
                    rank += 1
                assert rank == len(a) - case["dimension"]


def test_distinct_neighbor_cell_product_in_global_hamiltonian():
    car = federation(0)
    pair = PairSector(car)
    i,j = (p//12 for p in car.external[0])
    profiles = confinement.gate_profiles(car)
    modes = []
    for cell in (i,j):
        _, a, out = confinement.blocks(profiles[cell], "single", 0)
        _, v = confinement.maximal_space(a, out.values())
        modes.append(v.astype(float))
    factors = np.kron(*modes)
    indices = [pair.index(12*i+a,12*j+b) for a in range(12) for b in range(12)]
    columns = np.zeros((pair.P, factors.shape[1]))
    columns[indices] = factors
    for strength in (-8,0,8):
        np.testing.assert_array_equal(pair.hamiltonian(strength) @ columns, 6*columns)


def test_local_pulses_on_full_native_collar_and_registered_sign_read():
    k, v = law.local_term_parts()
    binary = [law.SRC.index((a,b)) for a in range(2) for b in range(2)]
    for strength in (8,0,-8):
        native = expm(-0.5j*np.pi*(k+strength*v))
        reduced = expm(-0.5j*np.pi*law.seam_term(strength))
        np.testing.assert_allclose(native[:, binary], np.eye(law.DIM,9) @ reduced[:,binary], atol=2e-13)
    result = sign_read.campaign()
    expected = json.loads(sign_read.RECEIPT.read_text())
    sign_read.verify(result, expected)
    for strength, arm in result["arms"].items():
        assert arm["slopes"] == [int(strength)/2, -1/8]
    def readings(strength):
        return np.array(result["arms"][strength]["readings"]["0.0625"])[:,0]
    coherent, dephased = readings("8") - readings("-8")
    assert coherent > 0.10
    assert abs(dephased) < 0.001
    corrupted = deepcopy(expected)
    corrupted["arms"]["8"]["readings"]["0.0625"][0][-1] += 0.01
    with pytest.raises(ValueError):
        sign_read.verify(result, corrupted)


def test_published_release_and_complementary_outcomes():
    data = json.loads(transport.RECEIPT.read_text())
    # Sums of 30 rounded probabilities in the M10.3 deposited receipt.
    published = {"8": [0.990982505,0.989653268,0.980271946],
                 "0": [0.450001679,0.275768631,0.202905801]}
    for strength, known in published.items():
        hits = [sum(data["arms"][strength][t][:30]) for t in ("1","7","37")]
        np.testing.assert_allclose(hits, known, rtol=0, atol=2e-8)
    for readings in data["arms"].values():
        for bins in readings.values():
            assert len(bins) == 34 and min(bins) >= 0 and abs(sum(bins)-1) < 2e-11
    corrupt = deepcopy(data)
    corrupt["arms"]["-8"]["148"].pop()
    with pytest.raises(ValueError):
        transport.verify(data, corrupt)


def test_load_conjugation_on_all_qutrit_states_and_terminal_records():
    flip = np.eye(9)[::-1]
    seam = law.seam_term(8)
    np.testing.assert_array_equal(flip @ seam @ flip, seam)
    for ref, (a,b) in product((0,1), law.SRC):
        # One missing contact at a, none at b. The identity extends by summation.
        difference = 8*((2-a-ref)**2-(a-ref)**2)
        assert difference == 32*(1-ref)*(1-a)
    result = conjugation.campaign()
    conjugation.verify(result, json.loads(conjugation.RECEIPT.read_text()))
    for board in result.values():
        port = board["ports"][0] % 12
        native, exchange, centered = board["arms"]
        assert native["excess"]["slope"] - native["vacancy"]["slope"] == 8
        assert (native["excess"]["readings"]["0.0625"][port] -
                native["vacancy"]["readings"]["0.0625"][port]) > 0.10
        for law_control in (exchange, centered):
            assert law_control["excess"] == law_control["vacancy"]
    corrupted = deepcopy(result)
    corrupted["2"]["arms"][-1]["vacancy"]["readings"]["0.125"][-1] += 0.001
    with pytest.raises(ValueError):
        conjugation.verify(result, corrupted)
