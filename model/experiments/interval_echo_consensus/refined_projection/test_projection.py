from fractions import Fraction
from math import log

from .measure import agreement, push, run, stratum_optimum


def test_native_refinement_and_changed_wiring():
    results = run()["levels"]
    for depth, row in enumerate(results, 1):
        assert row["ancestor_equality"] == str(Fraction(5, 8)**depth)
        for record in row["marked_trials_by_menu"].values():
            assert record["complete_histories"] == 32**depth
            assert record["distinct_ordered_pair_readbacks"] == 16**depth
            assert record["half_contacts"] == 5*depth*32**depth
    assert results[0]["cover_weights_coarse_to_fine"] == ["1"]
    assert results[1]["cover_weights_coarse_to_fine"] == ["1/2", "1/2"]
    expected_excess = (3/8) * ((5/8)*log(5/2) + (3/8)*log(1/2))
    assert abs(results[1]["changed_wiring_KL_excess_nats"] - expected_excess) < 1e-12


def test_exact_projection_on_deeper_record_alphabet():
    p = stratum_optimum(3)
    assert agreement(p) == Fraction(125, 512)
    assert push(p) == stratum_optimum(2)
