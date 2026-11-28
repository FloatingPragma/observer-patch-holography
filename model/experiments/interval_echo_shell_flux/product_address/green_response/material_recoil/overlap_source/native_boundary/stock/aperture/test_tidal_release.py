from collections import defaultdict
from fractions import Fraction as Q
from itertools import permutations

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    tidal_release as tidal,
)


def test_closed_pair_census_from_local_currents_against_every_four_station_history():
    ends = (0, 3)
    layers = tidal.flight.matching_layers(4, [(0, 1), (1, 2), (2, 3)])
    groups = tidal.category_weights(*ends, layers)
    pop = (1, 2, 1)
    labels = (2, 9)
    dist = {(pop, ((), (), labels)): Q(1)}
    words = set(permutations((0, 1, 1, 2)))
    half = tidal.lab.gates()[1]
    w = np.array([[0, 7, 3, 8], [7, 0, 9, 4], [3, 9, 0, 5], [8, 4, 5, 0]])
    t = 1
    measured = []
    for launch in [(0, 0), (t, 0), (0, t), (t, t)]:
        sharp = defaultdict(Q)
        for word in words:
            for m in labels:
                sharp[((launch[0], *word[1:3], launch[1]), m)] += Q(
                    1, len(words) * len(labels)
                )

        def score():
            return sum(
                p
                * sum(
                    int(w[i, j])
                    for i in range(4)
                    for j in range(i)
                    if x[i] == x[j] == t
                )
                for (x, _), p in sharp.items()
            )

        readings = [score()]
        for sequence in (layers, layers[::-1]):
            service = [ends] + [
                (i, int(j))
                for partner in sequence
                for i, j in enumerate(partner)
                if i < j
            ]
            for i, j in service:
                after = defaultdict(Q)
                for (x, m), p in sharp.items():
                    col = 9 * (m if (i, j) == ends else 0) + 3 * x[i] + x[j]
                    for row in np.flatnonzero(half[:, col]):
                        q = list(x)
                        q[i], q[j] = divmod(int(row) % 9, 3)
                        after[tuple(q), int(row) // 9 if (i, j) == ends else m] += (
                            p * Q(int(half[row, col]), 8)
                        )
                sharp = after
            readings.append(score())
        measured.append(np.array(readings, object))
    _, C = tidal.native_source(dist, t, ends, layers, groups)
    exact = tidal.pair_read(w, ends, layers, C[3] - C[1] - C[2] + C[0])
    np.testing.assert_array_equal(
        exact, measured[3] - measured[1] - measured[2] + measured[0]
    )


def test_second_reader_on_the_same_contact_cannot_count_a_new_pair():
    w = np.array([[0, 3], [3, 0]])
    assert np.array_equal(tidal.pair_back(w, [np.array([1, 0])]), 4 * w)
