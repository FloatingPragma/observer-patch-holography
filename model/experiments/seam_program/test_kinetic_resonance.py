import json
from hashlib import sha256

import mpmath as mp
import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model.carrier import Carrier
from . import cubic_resonance as cube
from . import guide_scaling as guide


def sector(n, dim):
    sites = np.arange(n**dim).reshape((n,) * dim)
    edges = [
        (i, j)
        for axis in range(dim)
        for i, j in zip(
            sites.ravel(), np.roll(sites, -1, axis=axis).ravel(), strict=True
        )
    ]
    return PairSector(Carrier("periodic guide", sites.size, edges))


def test_cubic_relative_generator_is_native_on_coherent_states():
    n = 4
    physical = sector(n, 3)
    coordinates = np.array(np.unravel_index(np.arange(n**3), (n,) * 3)).T
    r = (coordinates[physical.bb] - coordinates[physical.aa]) % n
    indices = np.ravel_multi_index(r.T, (n,) * 3) - 1
    opposite = np.ravel_multi_index((-coordinates[1:] % n).T, (n,) * 3) - 1
    rng = np.random.default_rng(9701)
    state = rng.normal(size=n**3 - 1) + 1j * rng.normal(size=n**3 - 1)
    state = (state + state[opposite]) / 2
    kinetic, shell, _ = cube.punctured_cube(n)
    for kappa in (0, 0.968, 8):
        relative = kinetic @ state - 2 * kappa * shell * state
        np.testing.assert_allclose(
            physical.hamiltonian(kappa) @ state[indices], relative[indices], atol=1e-12
        )


def test_ring_reduction_against_full_native_release_and_boundary():
    n = 16
    physical = sector(n, 1)
    distances = np.minimum(physical.bb - physical.aa, n - (physical.bb - physical.aa))
    initial = np.full(physical.P, 1 / np.sqrt(physical.P))
    for kappa in (8.0, -1 / (4 * n), 0.0):
        data = guide.ring(n, kappa)
        row = data["readouts"][0]
        state = expm_multiply(
            -1j * row["time"] * n**2 * physical.hamiltonian(kappa), initial
        )
        observed = np.bincount(distances, weights=abs(state) ** 2)[1:]
        np.testing.assert_allclose(
            observed, row["separation_probabilities"], atol=1e-12
        )
        p1, p2 = 1.3, -0.8
        s = guide.exchange_amplitude(p1, p2, 1 / n, kappa)
        z1, z2 = np.exp(1j * np.array([p1, p2]) / n)

        def psi(x, y):
            return z1**x * z2**y + s * z2**x * z1**y

        energy = 2 - np.cos(p1 / n) - np.cos(p2 / n)
        np.testing.assert_allclose(
            energy * psi(0, 1),
            (1 - 2 * kappa) * psi(0, 1) - (psi(-1, 1) + psi(0, 2)) / 2,
            atol=3e-13,
        )


def test_cubic_green_function_and_archive():
    x = json.loads((cube.HERE / "cubic_resonance.json").read_text())
    assert (
        x["protocol_sha256"]
        == sha256(
            (cube.HERE / "cubic_resonance_protocol.json").read_bytes()
        ).hexdigest()
    )
    np.testing.assert_allclose(x["green_at_zero"], cube.green(0.0), rtol=1e-12)
    # Independent finite momentum sum at negative energy, with fast finite-volume convergence.
    k = 2 * np.pi * np.arange(96) / 96
    e = 2 * (
        3
        - np.cos(k[:, None, None])
        - np.cos(k[None, :, None])
        - np.cos(k[None, None, :])
    )
    np.testing.assert_allclose(cube.green(1.0), np.mean(1 / (e + 1)), rtol=1e-11)
    assert cube.independent_shell_check(8, x["kappa_star"], 1) < 1e-12
    pulse = x["resonance_counterpulse"]
    assert pulse["maximum_column_error"] < 1e-10
    expected = np.exp(2j * (x["kappa_star"] - cube.KAPPA) / 64)
    np.testing.assert_allclose(
        complex(pulse["diagonal_real"][-1], pulse["diagonal_imag"][-1]),
        expected,
        atol=1e-10,
    )
    for row in x["kinetic_energy_ratio_at_double_momentum"]:
        h = 1 / row["n"]
        np.testing.assert_allclose(row["ratio"], (np.sin(h) / np.sin(h / 2)) ** 2)
    for r in x["scattering"]:
        h = 1 / r["n"]
        for j, m in enumerate(x["physical_inverse_lengths"]):
            f = cube.scattering(h, r["kappa_eff"], -((h * m) ** 2), cube.green(h * m))
            np.testing.assert_allclose(
                f, r["imaginary_momentum_response"][j], rtol=1e-10, atol=1e-12
            )
            rad = 4 * np.sin(h * m / 2) ** 2
            z = cube.scattering(h, r["kappa_eff"], rad, cube.outgoing_green(rad))
            # Scattering is a complex amplitude. At the pole its real component
            # approaches zero, while roundoff in the resolvent perturbs both components.
            archived = complex(*r["real_momentum_response"][j])
            np.testing.assert_allclose(z, archived, rtol=1e-10, atol=1e-12)
            if r["n"] == 256 and r["arm"] == "critical":
                assert abs(z - 1j / m) < 0.004 and abs(f - 1 / m) < 0.004
    for r in x["packets"]:
        for reading in r["readouts"]:
            np.testing.assert_allclose(sum(reading["separation_bins"]), 1, atol=1e-9)
        if r["n"] == 8:
            fresh = cube.packet(8, r["arm"], x["kappa_star"])
            for a, b in zip(fresh["readouts"], r["readouts"], strict=True):
                np.testing.assert_allclose(
                    a["separation_bins"], b["separation_bins"], atol=1e-9
                )
                np.testing.assert_allclose(
                    a["uniform_return"], b["uniform_return"], atol=1e-9
                )


def test_outgoing_green_at_the_finest_resonance_against_high_precision():
    # Resolve the logarithmic saddle with Carlson's symmetric integral and
    # independent high-precision quadrature on either side of the saddle.
    with mp.workdps(40):
        energy = 4 * mp.sin(mp.mpf(1) / 1024) ** 2
        saddle = mp.sqrt(4 - energy)

        def real_segment(sign, width):
            def integrand(t):
                if t == 0:
                    return 0
                offset = sign * width * t**2
                v = saddle + offset
                complementary = (offset * (2 * saddle + offset) / 4) ** 2
                return 2 * width * t * mp.elliprf(0, complementary, 1) / (
                    mp.pi**2 * mp.sqrt(4 + v**2)
                )
            return mp.quad(integrand, [0, 1])

        real = real_segment(-1, saddle) + real_segment(1, mp.sqrt(8 - energy) - saddle)
        imaginary = mp.quad(
            lambda v: mp.elliprf(0, ((energy - v**2 - 4) / 4) ** 2, 1)
            / (mp.pi**2 * mp.sqrt(4 - v**2)),
            [0, mp.sqrt(energy)],
        )
        value = cube.outgoing_green(float(energy))
        assert abs(value.real - float(real)) < 2e-14
        assert abs(value.imag - float(imaginary)) < 2e-16


def test_guide_archive():
    saved = json.loads((guide.HERE / "guide_scaling.json").read_text())
    fresh = guide.scan()
    assert fresh["protocol_sha256"] == saved["protocol_sha256"]
    for a, b in zip(fresh["rows"], saved["rows"], strict=True):
        assert (a["n"], a["arm"]) == (b["n"], b["arm"])
        for m, k in zip(a["exchange"], b["exchange"], strict=True):
            np.testing.assert_allclose(list(m.values()), list(k.values()), atol=1e-12)
        for m, k in zip(a["readouts"], b["readouts"], strict=True):
            np.testing.assert_allclose(
                m["separation_probabilities"], k["separation_probabilities"], atol=1e-9
            )
