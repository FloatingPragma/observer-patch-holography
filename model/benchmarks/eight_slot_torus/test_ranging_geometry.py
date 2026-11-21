from fractions import Fraction

import numpy as np
import pytest

from .ranging_geometry import PAIRS, _euclidean_obstruction, _octahedron, classify, distance_bands
from .station_ranging import tag_outcomes


def test_tag_timeouts_and_false_echoes_are_never_graded_as_round_trips():
    data = np.zeros((13,13), dtype=np.int64)
    data[0,0], data[0,6], data[7,0], data[7,1], data[1,6], data[6,7] = (3,5,7,11,13,17)
    result = tag_outcomes(data, 0)
    assert result == {"both_missed": 3, "receiver_missed_sender_clicked": 5,
                      "receiver_clicked_sender_missed": 7, "both_clicked": 41,
                      "wrong_receiver_letter": 18, "wrong_sender_letter": 17,
                      "correct_both_letters": 13, "returned_receiver_letter": 13,
                      "changed_receiver_letter": 28, "sender_correct_without_receiver": 5}
    assert tag_outcomes(data, 1)["correct_both_letters"] == 0


def test_orientation_is_conditioned_before_averaging_and_misses_are_retained():
    n = 2**32
    h = np.zeros((2,6,6,13,13), dtype=np.int64)
    for i,j in PAIRS:
        h[:,i,j,1,3] = n
        h[:,j,i,0,1] = n//2  # False echo.
        h[:,j,i,5,0] = n//4  # Capture followed by timeout.
        h[:,j,i,9,11] = n//4
    bands = distance_bands(h, trials=n, tv=Fraction(0))
    assert all(Fraction(599,100) < lo < 6 < hi < Fraction(601,100)
               for block in bands for lo,hi in block.values())
    h[:,1,0,9,11] = 0
    h[:,1,0,0,0] = n//4
    assert distance_bands(h, trials=n)[0][0,1] is None
    h[1,0,0,0,0] = n
    with pytest.raises(ValueError, match="self-pair"):
        distance_bands(h, trials=n)


def test_rare_arrivals_cannot_beat_the_declared_source_bias():
    n = 2**34
    h = np.zeros((2,6,6,13,13), dtype=np.int64)
    for i,j in PAIRS:
        for a,b in ((i,j),(j,i)):
            h[:,a,b,0,0], h[:,a,b,1,1] = n-1, 1
    assert all(b is None for b in distance_bands(h, trials=n)[0].values())
    bands = distance_bands(h, trials=n, tv=Fraction(0))[1]
    assert bands[1,4] == (1,12)
    verdict = classify(bands, ((0,1), (2,3), (4,5)))
    assert all(row["six_station_verdict"] == "unresolved" for row in verdict["thresholds"])


@pytest.mark.parametrize("tick", [1, 12])
def test_observations_at_either_deadline_feed_the_exact_certificate(tick):
    n = 2**34
    h = np.zeros((2,6,6,13,13), dtype=np.int64)
    for i,j in PAIRS:
        h[:,i,j,tick,tick] = h[:,j,i,tick,tick] = n
    bands = distance_bands(h)[0]
    assert all(isinstance(x, Fraction) for band in bands.values() for x in band)
    assert all(lo <= tick <= hi for lo,hi in bands.values())
    assert classify(bands, ((0,1), (2,3), (4,5)))["thresholds"][0]["global_3d_fit_excluded"]


def test_negative_type_witness_and_single_global_fit_are_distinct():
    opposite = ((0,1), (2,3), (4,5))
    bad = {p: (Fraction(2),Fraction(2)) if p in opposite else (1,1) for p in PAIRS}
    x = _euclidean_obstruction(bad, Fraction(1,20))
    assert x and sum(x["zero_sum_coefficients"]) == 0 and Fraction(x["quadratic_form_lower"]) > 0
    good = {p: (2,2) if p in opposite else (Fraction(1414,1000),Fraction(1415,1000)) for p in PAIRS}
    eps = Fraction(1,100)
    witness = _octahedron(good, eps, opposite)
    radius2 = Fraction(witness["radius_squared"])
    for p,(lo,hi) in good.items():
        d2 = (4 if p in opposite else 2)*radius2
        assert ((1-eps)*hi)**2 <= d2 <= ((1+eps)*lo)**2
    assert _euclidean_obstruction(good, eps) is None
    good[0,1] = None
    assert _octahedron(good, eps, opposite) is None
    with pytest.raises(ValueError):
        classify(good, ((0,1),(0,2),(4,5)))
