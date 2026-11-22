import json
from itertools import combinations

import numpy as np
import pytest

from model.executable_model import Carrier, cell, federation
from . import adjacency_bound, record_action, record_clock, record_readiness, record_receipt
from . import record_release, registered_forecast, sign_graph


def test_primitive_record_queries_and_complete_operator():
    result = record_action.run()
    assert result["checked_word_cases"] == 6750
    assert [r["counts"] for r in result["readers"].values()] == [[8,8], [8,16], [8,32]]
    assert result["untransported_word"] == 9
    e = result["relative_rate_controls"]
    np.testing.assert_allclose([e["0.5"], e["1.0"], e["2.0"]], [16, 0, 32], atol=1e-12, rtol=0)
    bank = sign_graph.witness_masks()
    assert bank["sign_patterns"] == 256 and bank["rejected_masks"] == [128,64]
    assert bank["recovered_single_rates"] == [1]*8 and bank["recovered_zero_joint_rates"] == 247
    assert bank["equal_and_unequal_occupations"] == {"00":0,"01":8,"10":8,"11":0}
    assert bank["rival_equal_load_phase_rate"] == 16
    assert bank["frequency_register_bits"] == 8 and bank["maximum_forbidden_read_probability"] < 1e-23
    board = cell()
    bad = record_action.signed_table().copy()
    bad[2,0] = 0
    with pytest.raises(ArithmeticError):
        record_action.incidence_check(board, bad)
    with pytest.raises(ValueError):
        list(record_action.flags([4], "square"))


def test_marker_preparation_probes_and_static_lifts():
    assert record_clock.prepare_marker() < 1e-13
    assert record_clock.parity_clock() < 1e-13
    assert record_clock.literal_lift()["microscopic_states"] == 768
    probes = record_readiness.check()
    assert probes["odd"][1] > 1-1e-14 and probes["even"][1] < 1e-14
    bright, dark = record_readiness.parity_instrument()
    rho = np.eye(2)/2
    np.testing.assert_allclose(bright@rho@bright.T + dark@rho@dark.T, rho, rtol=0, atol=1e-15)
    assert sign_graph.microscopic_path()["max_multiplicity"] == 2
    assert sign_graph.microscopic_path(True)["max_multiplicity"] == 1


def test_lazy_private_graph_reversible_and_exact_on_every_small_basis_column():
    board = Carrier("open path", 3, [(0,1),(1,2)])
    walk = sign_graph.PrivateRecordWalk(board)
    pairs = list(combinations(range(board.n),2))
    width = len(walk.offers)
    banks = 2**width
    assert banks == 16
    matrix = np.zeros((len(pairs)*banks,)*2, dtype=int)
    for j, pair in enumerate(pairs):
        for s in range(banks):
            start = j*banks+s
            ends = list(walk.neighbors(pair,s))
            assert len(ends) == len(set(ends)) and (pair,s) not in ends
            for q, t in ends:
                matrix[pairs.index(q)*banks+t, start] = -1
            matrix[start,start] = len(ends)
    np.testing.assert_array_equal(matrix, matrix.T)
    np.testing.assert_array_equal(matrix.sum(axis=0), np.zeros(len(matrix), dtype=int))
    kinetic = np.array([[1,-1,0],[-1,2,-1],[0,-1,1]])
    for mask in range(banks):
        char = [(-1)**(s & mask).bit_count() for s in range(banks)]
        v = np.kron(np.eye(len(pairs), dtype=int), np.array(char)[:,None])
        energies = [sum(bit for j, bit in enumerate(sign_graph.defects(q,walk.offers)) if mask & (1<<j))
                    for q in pairs]
        np.testing.assert_array_equal(matrix@v, v@(kinetic+2*np.diag(energies)))
    with pytest.raises(ValueError):
        list(walk.neighbors((0,0), 0))
    with pytest.raises(ValueError):
        list(walk.neighbors((0,1), banks))


def test_proof_hypotheses_and_complete_projected_matrix():
    for board, bound in ((cell(),"121/137"), *((federation(j),"5/6") for j in range(3))):
        assert adjacency_bound.certificate(board)["all_time_adjacent_lower"] == bound
    board = federation(0)
    a, b = board.external[0]
    occupied = set(x for e in board.external for x in e)
    other = tuple(next(v for v in range(12*(x//12),12*(x//12+1)) if v not in occupied) for x in (a,b))
    with pytest.raises(ValueError):
        adjacency_bound.certificate(Carrier("duplicate bridge", board.n, [*board.seams, other]))
    for arm, resources in (("primitive",6480),("assembled",8040)):
        result = sign_graph.full_projection(arm)
        assert result["private_signs"] == resources and result["material_states"] == 28680


def test_receipt_rejects_lost_outcomes_and_nonfinite_numbers():
    measured = {"reads": [1, [0.01, 0.99]], "passed": True}
    for field in ({"reads": [1,[0.01]], "passed": True},
                  {"reads": [1,[0.01,0.99]], "passed": 1},
                  {"reads": [1,[float("nan"),0.99]], "passed": True}, {}):
        with pytest.raises(ValueError):
            record_receipt.verify(measured, field)
    with pytest.raises(ValueError):
        record_release.release("off", 0)
    with pytest.raises(ValueError):
        record_release.source("missing")


def test_held_launch_predicts_all_bins_and_distinguishes_assembled():
    saved = json.loads(record_receipt.RECEIPT.read_text())
    for actual, row in zip((registered_forecast.compare(a) for a in record_receipt.ARMS[:3]),
                           saved["externally_enrolled"], strict=True):
        record_receipt.verify(actual, row)
        assert actual["launch"] == [87,92]
        assert actual["family_rejected"] == (actual["arm"] == "assembled")
    assert len(saved["release"]) == 8
    for row in saved["release"].values():
        assert set(row["appointments"]) == {"1","7","37"}
        for p in row["appointments"].values():
            assert len(p) == 31 and min(p) >= 0 and abs(sum(p)-1) < 1e-9
    assert saved["release"]["primitive/16"]["continuous_floor"] > .9632
    assert saved["release"]["assembled/16"]["continuous_floor"] > .9599
    assert not saved["release"]["assembled/8"]["threshold_met"]
    assert 1-saved["release"]["off/16"]["appointments"]["37"][-1] < .203
