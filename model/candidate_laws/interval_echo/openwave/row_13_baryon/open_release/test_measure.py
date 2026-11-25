from copy import deepcopy
from fractions import Fraction
from itertools import combinations, product
import json

import numpy as np
import pytest
from scipy import linalg

from . import measure


@pytest.fixture(scope="module")
def replay():
    result = measure.run()
    measure.compare(result, json.loads(measure.ARCHIVE.read_text()))
    # Exercise the same primitive-type contract as both CLI replay modes.
    measure.compare(result, json.loads(json.dumps(result, allow_nan=False)))
    return result


def test_same_prepared_source_and_all_complete_board_certificates(replay):
    assert replay["source"]["sharp_face"] == [0, 1, 5]
    assert replay["source"]["target_dimension"] == 4
    assert .1995 < replay["source"]["finite_filter_admission"] < .1996
    assert replay["source"]["finite_vs_ideal_trace_distance_upper"] < .006
    audit = replay["source"]["spectral_audit"]
    assert Fraction(audit["selected_energy_upper"]) < Fraction(201, 2)
    assert Fraction(audit["target_cluster_width_upper"]) < Fraction(1, 10**6)
    assert [(r["boundary"]["connected_triples"], r["boundary"]["fragmented_boundary_triples"])
            for r in replay["boards"]] == [(1900, 5416), (7600, 21648), (30400, 86608)]
    for board in replay["boards"]:
        assert board["boundary"]["max_cross_departures"] == 12
        assert board["boundary"]["max_integer_Gram_row_sum"] == 28
    bound = replay["single_terminal_read_bound_uniform_in_exposure"]
    assert bound["finite_filter_connected_lower"] > .976
    assert bound["original_start_joint_lower"] > .194
    assert min(r["conditional_witness_lower"] for r in replay["all_board_local_colour_witness_bounds"]) > .005


def test_local_connected_predicate_and_boundary_gram_by_full_diagonalization():
    board = measure.cell()
    sector = measure.BinarySector(board, 3)
    edges = set(board.seams)
    counted = np.array([sum(tuple(sorted(e)) in edges for e in combinations(word, 2))
                        for word in sector.states])
    connected = counted >= 2
    certificate = measure.connected_boundary(board)
    assert certificate["connected_triples"] == int(connected.sum())
    # Independent full-space matrices distinguish connected internal motion
    # from every disconnected three-unit configuration on this small board.
    m = -2*sector.K[np.ix_(~connected, connected)]
    assert np.array_equal(m, np.rint(m))
    assert certificate["max_integer_Gram_row_sum"] == int((m.T@m).sum(axis=1).max())
    h = sector.hamiltonian(8, True)
    assert linalg.eigvalsh(h[np.ix_(~connected, ~connected)])[0] >= -16-1e-12


def test_native_energy_meter_and_same_source_changed_law(replay):
    for board in replay["boards"]:
        record = board["complete_seam_energy_records"]
        source = next(r for r in board["source_moment_diagnostics"] if r["source"] == "actual_finite_filter")
        assert abs(source["shifted_energy"]-record["same_source_shifted_energy_kappa8"]) < 1e-11
        assert record["same_source_shifted_energy_kappa8"] < -42
        assert record["same_source_energy_exchange_only"] > 0
        for row in record["nontrivial_edges"]:
            values = row["odd_both_other"]
            assert min(values) >= 0 and abs(sum(values)-1) < 2e-12


def test_local_triplet_witness_has_the_stated_spectral_diameter():
    saved = json.loads(measure.FRAMES.read_text())
    groups = saved["groups"]
    sector = measure.BinarySector(measure.cell(), 3)
    frames = []
    for frame in saved["sharp_energy_filtered"]["frames"]:
        u, _, vt = linalg.svd(np.array(frame), full_matrices=False)
        frames.append(u@vt)
    transform = np.zeros((27, len(sector.states)))
    for row, colour in enumerate(product(range(3), repeat=3)):
        for local in product(*[range(len(g)) for g in groups]):
            word = tuple(sorted(groups[k][local[k]] for k in range(3)))
            transform[row, sector.index[word]] = np.prod([frames[k][local[k], colour[k]] for k in range(3)])
    singlet = measure.EPS.ravel()@transform
    witness = np.outer(singlet, singlet)-transform.T@transform/3
    energies = linalg.eigvalsh(witness)
    assert abs(energies[0]+1/3) < 2e-14 and abs(energies[-1]-2/3) < 2e-14
    assert abs(energies[-1]-energies[0]-1) < 3e-14
    np.testing.assert_allclose(linalg.svdvals(measure.EPS.reshape(3, 9))**2,
                               np.ones(3)/3, rtol=0, atol=2e-15)


def test_checker_rejects_changed_binding_certificate_and_colour_prediction(replay):
    original = measure.ARCHIVE.read_bytes()
    for label in ("boundary", "witness"):
        wrong = deepcopy(replay)
        if label == "boundary":
            wrong["boards"][0]["boundary"]["max_integer_Gram_row_sum"] += 1
        else:
            wrong["all_board_local_colour_witness_bounds"][1]["conditional_witness_lower"] = 0.
        with pytest.raises(ArithmeticError):
            measure.compare(wrong, replay)
    assert measure.ARCHIVE.read_bytes() == original
