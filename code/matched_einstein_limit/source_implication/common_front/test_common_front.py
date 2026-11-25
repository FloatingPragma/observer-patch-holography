from math import factorial

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .certificate import certify
from .measure import Ruler, native, optics, run


def test_continuous_time_arrivals_approach_the_common_clock_coefficient():
    board = optics.pairing.formation.federation(1)
    route = optics.pairing.formation.region(board)
    ruler = Ruler(board, route, 3)
    state = optics.source()["state"]
    p, q = route[0], route[3]
    kind = ruler.index[ruler.ids[(p, q)]]
    leading = ruler.coefficients(state)[kind] / (2**6 * factorial(3) ** 2)
    # The two charges and directions propagate through their actual Hamiltonians.
    times = np.array([0.04, 0.02, 0.01, 0.005])
    observations = []
    for charge in (1, -1):
        pooled = np.zeros(len(times))
        for parity in (0, 1):
            h = native.generator(board, route, charge, parity)
            seeds = np.zeros((h.shape[0], 8))
            seeds[4 * p : 4 * p + 4, :4] = np.eye(4)
            seeds[4 * q : 4 * q + 4, 4:] = np.eye(4)
            for i, t in enumerate(times):
                pulses = expm_multiply(-0.5j * t * h, seeds)
                f = pulses[4 * q : 4 * q + 4, :4]
                r = pulses[4 * p : 4 * p + 4, 4:]
                effect = (f.conj().T @ f + r.conj().T @ r) / 2
                pooled[i] += np.trace(effect @ state[parity]).real
        rescaled = pooled / times**6
        errors = np.abs(rescaled / leading - 1)
        assert np.all(errors[1:] < errors[:-1] / 3)
        assert errors[-1] < 5e-4
        observations.append(rescaled)
    assert abs(observations[0][0] - observations[1][0]) > 0.001


def test_clock_reference_and_own_volume_measurement():
    p = run()["panels"]
    assert [row["edges"] for row in p] == [2520, 3600, 4704, 6672]
    assert [row["blank_reference"] for row in p] == [
        [2, 3],
        [12, 44],
        [13, 792],
        [21, 792],
    ]
    assert np.max(np.abs(p[0]["slopes"]["unpair"])) < 1e-12
    assert p[2]["slopes"]["cap_six_matched_ordinary"][0][3] > 0.06
    assert p[2]["slopes"]["cap_six_matched_unpair"][0][3] > 0.03
    own = p[3]["own_volume"]
    assert 83 < own["ordinary"]["individual_gain"] < 84
    for name in ("ordinary", "pairing"):
        m = own[name]
        observed = m["finite"][-1]["centered_area_slope"]
        assert abs(observed - m["area_slope"]) < 1e-8


def test_independent_source_enclosures_and_held_calibration():
    data = certify()
    three, four = data["3"], data["4"]
    assert three["direct"][1] < three["two_link_lower"]
    assert three["convex_line_order_gap"] > 0.56
    assert three["four_gram"] == "negative_direction"
    assert four["four_gram"] == "positive"
    assert four["five_gram"] == "negative_direction"
    assert four["fixed_area"]["ordinary"][1] < -0.0048
    assert four["fixed_area"]["unpair"][0] > 0.000015
    assert -0.00588 < four["held_fraction"][0] < four["held_fraction"][1] < -0.00573
