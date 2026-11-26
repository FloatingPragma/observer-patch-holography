from itertools import combinations

import mpmath as mp
import numpy as np
import pytest
from mpmath.ctx_iv import MPIntervalContext

from . import measure


@pytest.fixture(scope="module")
def data():
    return measure.run()


def test_native_contacts_routes_and_source(data):
    assert data["ring_size"] == 246 and data["half_filling"] == 123
    assert data["source_orbital_maximum_error"] < 1e-12
    assert data["adjacent_source_rotations"] == 246 * 245 // 2
    controls = data["contact_controls"]
    assert controls["native"]["pulse_operator_error"] < 1e-12
    assert controls["exchange_removed"]["transfer_probability"] < 1e-25
    assert controls["contrast_halved"]["pulse_operator_error"] > 1
    passive = data["all_open"]
    assert passive["sea_energy_lower_bound"] == 3900
    assert passive["same_filling_comparator_energy"] == 1008
    assert passive["outside_receiver_quadratic_coefficient"] == 1 / 8


def test_overlapping_regions_and_signed_local_inversion(data):
    assert len(data["intervals"]) == 8 and len(data["local_stress"]) == 24
    for row in data["intervals"]:
        assert row["commutator_norm"] < 1e-11
        assert row["uniform_control_commutator_norm"] > 0.27
        assert max(x["relative_difference"] for x in row["readings"].values()) < 0.002
    for row in data["local_stress"]:
        assert abs(row["bond_energy"] - row["energy_inversion"]) < 1e-13
        assert row["entropy_relative_difference"] < 0.0075
        if row["source"] == "signed_rotation":
            assert (row["right_bond"] == 51) == (row["entropy_inversion"] > 0)


def test_four_site_rescaling_obstruction_with_directed_rounding():
    iv = MPIntervalContext()
    iv.dps = 40
    n = 246
    angle = iv.pi / n
    c1 = 1 / (n * iv.sin(angle))
    c3 = -1 / (n * iv.sin(3 * angle))
    gap = iv.sqrt(((c3 - c1) / 2) ** 2 + c1**2)
    centre = (1 + c3 + c1) / 2
    p0, p1 = centre - gap, centre + gap
    d1, d2 = [
        -n * iv.sin(j * angle) * iv.sin((4 - j) * angle) / iv.sin(4 * angle)
        for j in (1, 2)
    ]
    gap = iv.sqrt(d2**2 + 4 * d1**2)
    r0 = iv.log((1 - p0) / p0) / ((d2 + gap) / 2)
    r1 = iv.log((1 - p1) / p1) / ((d2 - gap) / 2)
    assert iv.mpf("0.987") < r0 < iv.mpf("0.988")
    assert iv.mpf("1.014") < r1 < iv.mpf("1.016")


def test_rare_mode_entropy_from_direct_high_precision_covariance(data):
    n, length = data["ring_size"], 24
    with mp.workdps(45):
        def covariance(i, j):
            d = i - j
            return mp.mpf("0.5") if d == 0 else mp.sin(mp.pi * d / 2) / (n * mp.sin(mp.pi * d / n))

        full_covariance = mp.matrix([[covariance(i, j) for j in range(length)] for i in range(length)])
        eigenvalues, basis = mp.eigsy(full_covariance)
        smallest = float(eigenvalues[0])
        holes = [mp.exp(1j * i * (mp.pi / 2 - mp.pi / n)) / mp.sqrt(n) for i in range(length)]
        particles = [mp.exp(1j * i * (mp.pi / 2 + mp.pi / n)) / mp.sqrt(n) for i in range(length)]
        projected_hole, projected_particle = (basis.T * mp.matrix(v) for v in (holes, particles))
        entropy = float(sum(
            mp.log((1 - p) / p) * (abs(projected_particle[i])**2 - abs(projected_hole[i])**2)
            for i, p in enumerate(eigenvalues)
        ))
    row = data["intervals"][-1]
    np.testing.assert_allclose(smallest, row["minimum_mode_probability"], rtol=3e-6, atol=0)
    np.testing.assert_allclose(entropy, row["readings"]["plane_wave"]["entropy"], atol=1e-10, rtol=0)


def gaussian_density(c):
    """All occupation-word matrix elements, used on a four-station interval."""
    n = len(c)
    one = np.eye(n)
    odds = c @ np.linalg.inv(one - c)
    rho = np.zeros((2**n, 2**n), dtype=complex)
    for size in range(n + 1):
        words = list(combinations(range(n), size))
        for x in words:
            for y in words:
                i, j = (sum(1 << k for k in w) for w in (x, y))
                rho[i, j] = np.linalg.det(odds[np.ix_(x, y)])
    return np.linalg.det(one - c).real * rho


def test_spin_interval_entropy_and_synthetic_confidence(data):
    momenta, phase, filled = measure.ring_modes(data["ring_size"])
    kernel = phase[46:50, filled] @ phase[46:50, filled].conj().T
    hole, particle = (v[46:50] for v in measure.tangents(momenta, phase, filled)["coherent_pair"])
    excited = kernel + np.outer(particle, particle.conj()) - np.outer(hole, hole.conj())
    vacuum = gaussian_density(kernel)
    alternate = gaussian_density(excited)
    np.testing.assert_allclose([vacuum.trace(), alternate.trace()], [1, 1], atol=1e-12)
    levels, modes = np.linalg.eigh(vacuum)
    derivative = -np.log(levels) @ np.diag(modes.conj().T @ (alternate - vacuum) @ modes).real
    readings = data["intervals"][0]["readings"]["coherent_pair"]
    np.testing.assert_allclose(derivative, readings["entropy"], rtol=0, atol=2e-12)
    bank = data["count_bank"]
    assert sum(map(sum, bank["reference_bonds"])) == 3 * bank["shots_per_setting"]
    for label in ("entropy", "weighted_energy"):
        lo, hi = bank[f"{label}_interval"]
        assert 0 < lo <= readings[label] <= hi
