import numpy as np

from .__main__ import Patch, derivatives, hardware, reduced


def test_modular_slope_against_two_rail_spectrum():
    patch = Patch()
    slope = 24 / 35 * np.log(4 / 11)
    biased = derivatives(patch, 3 / 4)
    np.testing.assert_allclose(biased['reference_spectrum'], [0.2, 0.25, 0.55], atol=1e-12)
    assert abs(biased['regional_modular_slope'] - slope) < 1e-12
    assert abs(derivatives(patch, 0.5)['regional_modular_slope']) < 1e-12
    assert biased['maximum_diagonal_derivative'] < 1e-12


def test_splitter_hides_launch_from_cut_but_preserves_patch_coherence():
    patch = Patch()
    for control in (False, True):
        difference = reduced(patch.state(1, splitter=control) - patch.state(0, splitter=control), [0, 2])
        distance = sum(abs(np.linalg.eigvalsh(difference))) / 2
        assert abs(distance - (0 if control else 0.8)) < 1e-12
    assert abs(patch.region(1)[1, 2] - patch.region(0)[1, 2]) > 0.6
    for level in range(3):
        board = hardware(level)
        inside = set(a for a, _ in board['boundary_seams'])
        a, b, c, d = board['rails_ABCD']
        assert a in inside and b not in inside and d not in inside and c >= 12
