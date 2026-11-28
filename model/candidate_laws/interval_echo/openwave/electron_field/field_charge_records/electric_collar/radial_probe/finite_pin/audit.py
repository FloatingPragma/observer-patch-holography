"""Independently enclose selected finite-pin energies using exact row ratios."""

from fractions import Fraction
import json
from math import inf, nextafter
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu

from .....fan_path.fixture import FanPath, archived_fixture
from ..apparatus import RadialCollar


SETTINGS = ((0, 0), (1, 0), (0, 1), (1, 1))
WEIGHTS = (1, -1, -1, 1)


def prepared_matrix(collar, rate, pins, radius, background, setting):
    terms = collar.terms
    matrix = terms["fan_rings"] + terms["outer_ring"] + 2 * rate * terms["connector"]
    if setting[0]:
        matrix += terms["source_hop"] + terms["source_degree"]
    if setting[1]:
        matrix += terms["probe_hop"] + terms["probe_degree"]
    target = collar.radial[radius]
    initial = [
        collar.ids[collar.hub, target, w | (background << collar.inside_bits)]
        for w in collar.start_words
    ]
    reached = {initial[0]}
    agenda = [initial[0]]
    while agenda:
        row = matrix.getrow(agenda.pop())
        for neighbor, value in zip(row.indices, row.data, strict=True):
            if value < 0 and neighbor not in reached:
                reached.add(neighbor)
                agenda.append(neighbor)
    assert initial[1] in reached
    scope = sorted(reached)
    pin = []
    for x in scope:
        source, probe, _ = collar.basis[x]
        assert setting[0] or source == collar.hub
        assert setting[1] or probe == target
        pin.append(pins[0] * (source != collar.hub) + pins[1] * (probe != target))
    return matrix[scope][:, scope].toarray() + np.diag(pin)


def exact_row_bounds(matrix):
    estimate = eigvalsh(matrix, subset_by_index=(0, 0))[0]
    inverse = splu(csc_matrix(matrix - (estimate - 2**-12) * np.eye(len(matrix))))
    trial = np.ones(len(matrix))
    for _ in range(6):
        trial = inverse.solve(trial)
        trial /= max(trial)
    assert all(trial > 0)
    # The proposal can use floating point. Once it is positive, the
    # exact ratios bound the ground energy by Perron positivity.
    vector = [Fraction(float(x)) for x in trial]
    ratios = [
        sum(Fraction(float(row[j])) * vector[j] for j in np.flatnonzero(row))
        / vector[i]
        for i, row in enumerate(matrix)
    ]
    return min(ratios), max(ratios)


def combine(intervals, signs):
    return tuple(
        sum(
            sign * interval[end if sign == 1 else 1 - end]
            for interval, sign in zip(intervals, signs, strict=True)
        )
        for end in (0, 1)
    )


def bounds(collar, rate, pins, radius, word):
    return combine(
        [
            exact_row_bounds(prepared_matrix(collar, rate, pins, radius, word, s))
            for s in SETTINGS
        ],
        WEIGHTS,
    )


def read_record(arm, radius, word, setting):
    return combine(
        [
            (Fraction(e) - Fraction(err), Fraction(e) + Fraction(err))
            for e, err, *_ in arm["presence_raw_records"][radius][word][setting]
        ],
        WEIGHTS,
    )


def main():
    data = json.loads((Path(__file__).with_name("results.json")).read_text())
    positive, negative = data["cases"][:2]
    with archived_fixture() as archive:
        path = FanPath(archive, 4, 1, asymmetric=True)
        collar = RadialCollar(path, 1, [334, 365])
        opposite = RadialCollar(path, -1, [334, 365])
        for arm, j in ((0, 0), (0, 3), (1, 0), (2, 0)):
            bank = positive["connectors"][arm]
            rate = bank["connector_rate"]
            pins = data["source_and_probe_pin_settings"][j]
            exact = [bounds(collar, rate, pins, r, b) for r, b in enumerate((45, 47))]
            saved = [read_record(bank, r, b, j) for r, b in enumerate((45, 47))]
            for outer, inner in zip(saved, exact, strict=True):
                assert outer[0] <= inner[0] <= inner[1] <= outer[1]
            radial = combine(exact, (-1, 1))
            assert radial[1] - radial[0] < Fraction(1, 10**10)
            if rate == 0:
                assert radial[0] <= 0 <= radial[1]
            else:
                assert radial[1] < 0
            low = nextafter(float(radial[0]), -inf)
            high = nextafter(float(radial[1]), inf)
            print(f"Pin {pins}, ring {rate:g}, ΔU in [{low}, {high}]")
        reverse = bounds(opposite, 0.5, [64, 64], 0, 45)
        reverse_saved = read_record(negative["connectors"][0], 0, 45, 0)
        assert reverse_saved[0] <= reverse[0] <= reverse[1] <= reverse_saved[1]
        same_pin = combine((bounds(collar, 0.5, [64, 64], 0, 45), reverse), (1, 1))
        assert same_pin[0] > 0
        print("Verified exact row bounds for the matched pin and polarity controls.")


if __name__ == "__main__":
    main()
