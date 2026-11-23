import numpy as np
from scipy.sparse import csr_matrix

from .measure import TICKS, aperture, balance
from model.experiments.contact import radar


def test_analytic_crossing_with_unequal_launch_volumes():
    # Exactly two ports exchange. A third outside port is a spectator.
    h = csr_matrix([[0., -0.5, 0.], [-0.5, 0., 0.], [0., 0., 3.]])
    hand = radar.DT*np.asarray(TICKS)
    result = aperture(h, np.array([0]), hand, dense=True, held=False)
    expected = np.sin(hand/2)**2
    np.testing.assert_allclose(result["outward"], expected)
    np.testing.assert_allclose(result["inward"], expected/2)
    np.testing.assert_allclose(result["inward"]/(result["outward"]+result["inward"]), 1/3)
    np.testing.assert_allclose(result["area"], 4*expected/hand**2)
    assert result["leading_area"] == 1


def test_known_area_entropy_cancellation_and_independent_counterdirection():
    entropy = np.array([0.25, -0.5, 0.25, 0.5])
    match = balance(-entropy/3, entropy)
    assert np.isclose(match["coefficient"],3)
    assert match["relative_residual"] < 1e-13
    mismatch = balance(np.array([1., 0., -1., 0.]), entropy)
    assert np.isclose(mismatch["relative_residual"],1.)
