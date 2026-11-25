"""An omitted amplitude can interfere at first order in a complete detector."""

import numpy as np
from .measure import enclosure


def test_unfiltered_enclosure_contains_dense_flight_with_nonzero_cross_terms():
    rng = np.random.default_rng(109)
    eig, _ = np.linalg.qr(rng.normal(size=(4, 4)))
    energy = np.array([0.0, 0.13, 0.91, 2.0])
    x = np.array([1.0, 0.42, -0.8, 0.06])
    g = np.array([0.8, -0.6])
    dwell, strength = 3.25, 0.05
    result = enclosure(
        eig[:, :2], energy[:2], g, 1.0, 0.0, 0.0, 0.0, x, x, strength, dwell
    )
    full, retained, omitted = [], [], []
    for sign in (-1, 1):
        y = np.exp(sign * 1j * strength * x) * (eig[:, :2] @ g)
        evolved = np.exp(-1j * dwell * energy) * (eig.T @ y)
        full.append(float(x @ (abs(eig @ evolved) ** 2)))
        retained.append(float(x @ (abs(eig[:, :2] @ evolved[:2]) ** 2)))
        omitted.append(float(np.linalg.norm(evolved[2:]) ** 2))
    difference = full[1] - full[0]
    assert result["displacement"][0] <= difference <= result["displacement"][1]
    assert abs(difference - (retained[1] - retained[0])) > sum(omitted)
