"""Direct addressed-memory probability spectra, including every vacant outcome."""

from fractions import Fraction
import json
from itertools import permutations

import numpy as np
import pytest
from model.executable_model import federation
from ... import shared_variance, volume
from . import word_entropy


def word_spectrum(law, geometry, per_ancestor):
    # Sum over all ordered source placements, tracing away outside stations.
    # Keeping distinct admitted positions separate matters for local tags.
    N = 80
    dist = len(geometry)
    head_law = law.sum(1)
    product = np.outer(head_law, head_law)
    slots = {}
    for a, b in permutations(range(N), 2):
        first = a if a < dist else -1
        second = b if b < dist else -1
        address = (first, second)
        if min(address) >= 0:
            outcome = (
                law
                if not per_ancestor or geometry[first] == geometry[second]
                else product
            )
        elif max(address) < 0:
            outcome = np.ones(1)
        else:
            outcome = head_law
        slots[address] = slots.get(address, 0.0) + outcome.ravel()
    values = np.concatenate([slots[a] for a in sorted(slots)]) / (N * (N - 1))
    assert values.min() > 0 and values.sum() == pytest.approx(1.0, abs=5e-15)
    return values


@pytest.mark.parametrize("per_ancestor", [False, True])
def test_full_regional_entropy(per_ancestor):
    stored = json.loads((word_entropy.HERE / "word_entropy.json").read_text())
    sectors = json.loads((word_entropy.HERE / "ancestor_limit.json").read_text())[
        "midpoint_surveys"
    ][0]
    geometry = np.repeat(np.arange(20), sectors["ancestor_populations"])
    _, _, _, prep = volume.radar_source_and_hand(federation(0))
    laws, _ = shared_variance.source_words(prep)
    rev, fresh, shared = [word_spectrum(law, geometry, per_ancestor) for law in laws]
    diff = shared - fresh
    ell = Fraction(sectors["radius"], sectors["eccentricity"])
    measured = float(-diff @ np.log(rev)) / float(ell**4)
    name = "ancestor" if per_ancestor else "global"
    lo, hi = stored["finite_resolutions"][0][name + "_entropy_slope_per_fourth"]
    assert (lo + hi) / 2 == pytest.approx(measured, rel=6e-12)
    assert abs(diff @ np.log(fresh)) < 2e-15

    def h(p):
        return -p @ np.log(p)

    assert h(shared) < h(fresh)
    for t, positive in [(0.5, True), (0.66, True), (0.68, False)]:
        ds = h((1 - t) * rev + t * shared) - h((1 - t) * rev + t * fresh)
        assert bool(ds > 0) == positive


def test_receipt_and_bracket():
    saved = json.loads((word_entropy.HERE / "word_entropy.json").read_text())
    assert word_entropy.run() == saved
    a, b = map(
        Fraction, saved["zero_of_shared_minus_fresh_entropy_at_common_admixture"]
    )
    sums = word_entropy.probabilities()
    assert 0.668 < a < b < 0.669 and b - a < Fraction(1, 10**10)
    assert word_entropy.contrast(sums, a).a > 0 and word_entropy.contrast(sums, b).b < 0
