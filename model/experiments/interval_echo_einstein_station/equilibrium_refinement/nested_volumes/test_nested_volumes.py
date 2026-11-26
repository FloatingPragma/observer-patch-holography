"""Reachable hardware, independent geometry and charge, different admission laws."""

from fractions import Fraction
from functools import lru_cache
from math import comb, log
import json
import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix
from model.candidate_laws.interval_echo import candidate
from ..tracer_flux import measure as tracer
from . import flight, measure


@lru_cache(maxsize=1)
def acquired():
    return measure.run()


def saved():
    return json.loads(
        measure.Path(measure.__file__).with_name("results.json").read_text()
    )


def test_native_generators_propagation_and_numerical_error():
    h, states, _ = flight.hardware(1)
    older, order, _, _ = tracer.source(candidate.KAPPA)
    assert order == [(tag[0], w) for w, tag in states]
    np.testing.assert_array_equal(h.toarray(), older)
    toy = np.array([[15, -2, 1], [-2, -3, -4], [1, -4, 19]], dtype=float)
    for z, t in [(Fraction(1, 8), 0.125), (Fraction(1, 4), 0.25)]:
        (p,), proof = flight.propagators(csr_matrix(toy), np.array([0, 2]), [z])
        exact = abs(expm(-1j * t * toy / 2)[:, [0, 2]]) ** 2
        assert (
            np.max(abs(p - exact).sum(0))
            < proof["column_probability_error_upper"]
            < 1e-7
        )


def test_two_nested_radii_ordinary_holdout_and_no_matching_control():
    for receipt in (acquired(), saved()):
        interior = receipt["independent_interior"]["measured_held_capacity_nats"]
        assert interior[0] < 0.5 * log(2 / 5) < interior[1] < -0.45
        for j, bank in enumerate(receipt["tagged_volumeters"], 1):
            assert bank["states"] == (660, 2970)[j - 1]
            assert bank["dynamics_bound"]["column_probability_error_upper"] < 1e-8
            work = (0, log(5) / 8)[j - 1]
            for frame in bank["times"]:
                assert (
                    frame["entirely_counted_prediction_nats"][0]
                    < work
                    < frame["entirely_counted_prediction_nats"][1]
                )
                for name, stock in measure.STOCKS.items():
                    arms = frame["arms"][name]
                    for rad, arm in enumerate(arms, 1):
                        sites = [
                            i for i, d in enumerate(receipt["clock_radii"]) if d <= rad
                        ]
                        count = sum(
                            comb(len(sites) - sum(t in sites for t in w), j)
                            for w in stock
                        ) / len(stock)
                        assert (
                            arm["independently_counted_99"][0]
                            < count
                            < arm["independently_counted_99"][1]
                        )
                        assert abs(arm["recovered_placements"] - count) < 1e-10
                        for d in ("forward", "reverse"):
                            assert sum(arm[d]["outcomes"]) == measure.SHOTS
                if j == 2:
                    lo, hi = frame["finite_shell_secant_prediction_nats"]
                    assert 0 < lo < log(5) * (np.sqrt(85) - 9) / 2 < hi
                    assert frame["arms"]["split"][1]["without_source_match_99"][0] > 36
        for probes in range(1, 6):
            held = (comb(4, probes) + comb(6, probes)) / 2 - comb(5, probes)
            ordinary = comb(5, probes) - comb(4, probes)
            assert held / ordinary == (probes - 1) / (2 * (6 - probes))


def test_material_closures_and_gluing():
    for receipt in (acquired(), saved()):
        shutter = receipt["shutters"]
        assert len(shutter["closed_word_histograms"]) == 55
        assert shutter["native_gates"]["echo_read_undo"] == 21690
        for row in shutter["closed_word_histograms"]:
            k, stop = row["k"], row["stop"]
            counts = row["permutations_by_inserted_prefix"]
            assert sum(counts) == 120
            for steps in range(stop):
                assert Fraction(
                    sum(counts[steps + 1 :]), sum(counts[steps:])
                ) == Fraction(5 - steps - k, 5 - steps)
            assert counts[-1] * comb(5, k) == 120
        a, b = shutter["measured_held_nats"]
        assert a < 0.5 * log(2 / 5) < b < -0.45
        a, b = shutter["arms"]["spread"]["mean_work_nats"]
        assert (
            a < log(10) / 2 < b
            and shutter["arms"]["spread"]["negative_log_pooled_success"] < 0.6
        )
