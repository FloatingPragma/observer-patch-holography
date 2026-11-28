from collections import defaultdict
from fractions import Fraction as Q
from itertools import permutations

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    material_release as release,
)


def test_two_visits_preserve_the_joint_ancestry_and_resident_head():
    # Exhaust the sharp outcomes on a four-station cycle. Unlike the large
    # board reduction, this needs no assumptions about ancestor categories.
    ends = (0, 3)
    layers = release.matching_layers(4, [(0, 1), (1, 2), (2, 3)])
    ancestors = release.ancestry(ends, layers)
    original = set(permutations((0, 0, 1, 2)))
    labels = (2, 9)
    dist = {((2, 1, 1), ((), (), labels)): Q(1)}
    half = release.lab.gates()[1]
    for launch in ((0, 0), (1, 0), (0, 2)):
        paths = defaultdict(Q)
        for word in original:
            for mark in labels:
                paths[((launch[0], *word[1:3], launch[1]), mark)] += Q(
                    1, len(original) * len(labels)
                )
        for phase in (layers, layers[::-1]):
            service = [ends]
            for partners in phase:
                service += [(i, int(j)) for i, j in enumerate(partners) if i < j]
            for i, j in service:
                after = defaultdict(Q)
                for (word, mark), p in paths.items():
                    col = 9 * (mark if (i, j) == ends else 0) + 3 * word[i] + word[j]
                    for row in np.flatnonzero(half[:, col]):
                        dest = list(word)
                        dest[i], dest[j] = divmod(int(row) % 9, 3)
                        key = (tuple(dest), int(row) // 9 if (i, j) == ends else mark)
                        after[key] += p * Q(int(half[row, col]), 8)
                paths = after
        exact = np.full((4, 3), Q(0), object)
        for (word, _), p in paths.items():
            for i, t in enumerate(word):
                exact[i, t] += p
        prediction = release.trial(dist, ends, layers, ancestors, launch)[0]
        np.testing.assert_array_equal(exact, prediction)


def test_two_ends_on_one_blank_share_one_swap_coin():
    # Independent one-body walkers would incorrectly put both loads on one port.
    partners = release.matching_layers(3, [(0, 1)])
    assert release.ancestry((0, 1), partners) == {(0, 1): Q(1, 2), (1, 0): Q(1, 2)}
