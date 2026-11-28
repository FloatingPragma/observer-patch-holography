import numpy as np
from mpmath import iv, mp

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    occupied_equilibrium as heat,
    occupied_release as flight,
)


def test_source_cancellation_with_initial_material_current():
    k = np.array([[2.0, -0.5, 0], [-0.5, -1.0, -0.5], [0, -0.5, 1.0]])
    marks = np.array([[1.0, 1.0, 0], [0.0, 1.0, 1.0]])
    inside = np.array([True, False, False])
    o = np.diag(inside.astype(float))
    measured = np.array([1, 1j, -0.5j])
    measured /= np.linalg.norm(measured)
    base, terms = flight.material_ops(k, marks, inside)
    second = []
    for a, b in [(0, 0), (0, 33), (32, 0), (32, 33)]:
        h = k + np.diag(a * marks[0] + b * marks[1])
        direct = (h @ o - o @ h) @ h - h @ (h @ o - o @ h)
        np.testing.assert_array_equal(direct, base + a * terms[0] + b * terms[1])
        second.append((measured.conj() @ direct @ measured).real)
    assert abs(flight.MIX @ second) < 1e-13


def test_heat_enclosure_against_high_precision_exponential():
    twice = np.array(
        [[24, -1, 0, 0], [-1, 5, -1, 0], [0, -1, -8, -1], [0, 0, -1, 3]], dtype=int
    )
    before = iv.dps
    try:
        iv.dps = 65
        values, energy = heat.heat_shells(twice, np.arange(4))
        with mp.workdps(80):
            h = mp.matrix(twice.tolist()) / 2
            gibbs = mp.expm(-mp.log(mp.mpf(8) / 7) * h)
            z = sum(gibbs[i, i] for i in range(4))
            exact = [gibbs[i, i] / z for i in range(4)] + [mp.mpf(0)] * 3
            exact.append(sum((h * gibbs)[i, i] for i in range(4)) / z)
            for result, val in zip([*values, energy], exact, strict=True):
                assert (
                    mp.make_mpf(result._mpi_[0]) <= val <= mp.make_mpf(result._mpi_[1])
                )
    finally:
        iv.dps = before
