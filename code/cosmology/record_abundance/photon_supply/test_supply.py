import numpy as np
from cosmology.record_abundance.photon_supply import measure


def test_photon_energy_and_scale_factor_weighting():
    photons = measure.photon_density(measure.T0)
    assert 2.46e-5 < photons < 2.48e-5
    q = 1037 / 202
    shares = [0.1, 0.9]
    a = measure.history(q, measure.BC, [2000, 1e5], shares, photons)
    R = sum(
        f
        * measure.history(q, measure.BC, [z], [1], photons)[
            "lost_over_final_comoving_photon_energy"
        ]
        for f, z in zip(shares, (2000, 1e5))
    )
    np.testing.assert_allclose(a["neff_after"], measure.NSTAR * (1 + R))
    bounds = measure.run()["boards"][1]["calibration_bounds"]
    assert 48000 < bounds["baryon_ratio"]["minimum_single_burst_redshift"] < 49000
