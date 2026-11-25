"""Independent finite spin checks of the mode assignment and two-wing readout."""

from itertools import product
import json

import numpy as np
import pytest

from . import measure


@pytest.mark.parametrize("cut", [True, False])
def test_modes_exhaust_the_direct_material_spectrum(cut):
    for n, boundary, kappa in product((1, 2, 3, 5), (0, 2), (0, .125, 8)):
        states = list(product((0, 1), repeat=n))
        h = np.zeros((len(states), len(states)))
        for i, word in enumerate(states):
            ports = np.array([load for x in word for load in (2*x, 2-2*x)])
            boundaries = np.r_[boundary, ports]
            if cut:
                boundaries = ports
            for a, b in zip(boundaries, boundaries[1:]):
                h[i,i] += kappa*(a-b)**2 + .5*(a != b)
            for s in range(n):
                changed = list(word)
                changed[s] ^= 1
                h[states.index(tuple(changed)),i] = -.5
        energies = np.linalg.eigvalsh(h)
        predicted, _ = measure.energy_table(n, kappa, cut)
        np.testing.assert_allclose(energies-energies.min(), np.sort(predicted), atol=5e-12, rtol=0)


def test_calibration_and_dephasing_on_a_known_two_wing_state():
    left, right = [measure.calibrated_modes(b) for b in (0,2)]
    bell = (np.outer(left[:,0], right[:,1]) + np.outer(left[:,1], right[:,0]))/np.sqrt(2)
    np.testing.assert_allclose(measure.outcomes(bell, left, right, "XX", False)[:2,:2],
                               np.eye(2)/2, atol=1e-14)
    np.testing.assert_allclose(measure.outcomes(bell, left, right, "XX", True)[:2,:2],
                               np.ones((2,2))/4, atol=1e-14)


def test_deposited_counts_and_coherence():
    receipt = json.loads(measure.ARCHIVE.read_text())
    for panel in receipt["coherence"]:
        for setting in panel["reads"].values():
            assert len(setting["all_19_counts"]) == 19
            assert sum(setting["all_19_counts"]) == receipt["shots"]
    wall = [p for p in receipt["coherence"] if (p["contrast"],p["cut"],p["source"]) == (8,False,"wall")]
    assert wall[0]["simultaneous_99_witness_lower"] > 1.95
    assert wall[1]["simultaneous_99_witness_lower"] < 1
    archived = [p for p in receipt["energy_counts"] if p["contrast"] == 8 and not p["cut"]
                and p["wall"] and p["exposure"] == 32][0]
    assert sum(sum(row) for row in archived["one_each_mode_counts"]) == 2064278
    assert sum(sum(row) for row in archived["mode_number_counts"])+archived["ring_miss"] == 2**26
