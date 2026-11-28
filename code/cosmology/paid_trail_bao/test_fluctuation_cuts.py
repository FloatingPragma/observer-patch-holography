import numpy as np

from . import fluctuation_cuts as cuts
from .fluctuation_lensing import envelope


def test_upper_envelope_and_cut_radius_integral():
    p = np.linspace(0, 1, 1001)
    top = envelope(p)
    assert min(top - p**2 * (1 - p) ** 2) >= -1e-15
    assert np.max(np.diff(top, n=2)) < 1e-14
    assert np.allclose(
        cuts.integrated(
            np.array([1, 3]), np.array([0.2, 0.6]), np.array([0.5, 1, 2, 3, 5])
        ),
        [0, 0, 0.1, 0.4 / 3, 0.32],
    )


def test_free_outer_mass_and_convex_chord():
    tab = dict(
        radius_kpc=[2, 3, 4],
        baryons_kms=[0] * 3,
        error_kms=[0] * 3,
        observed_kms=np.sqrt([4 / 2, 8 / 3, 9 / 4]),
    )
    result = cuts.convex_test(tab, True)
    assert not result["feasible"] and result["witness"]["kind"] == "convex"
    tab["observed_kms"] = np.sqrt([4 / 2, 6 / 3, 9 / 4])
    assert cuts.convex_test(tab, True)["feasible"]
