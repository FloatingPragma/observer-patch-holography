"""Exactness and regression checks of the simulator.

Each check returns a short summary and raises on failure. ``python -m simulator
check`` prints them; pytest runs them as tests.
"""
from __future__ import annotations

import numpy as np

from .experiments import clock_journey, continuum_run
from .law import LAWS, dispersion, step, step_pair, walk_matrix
from .sector import sector_step
from .spectral import Packet

RNG_SEED = 20261120


def _random_state(shape, rng):
    psi = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    return psi / np.linalg.norm(psi)


def check_unitarity():
    rng = np.random.default_rng(RNG_SEED)
    worst = 0.0
    for law in LAWS.values():
        psi = _random_state((5, 4, 3, 6, 8), rng)
        for _ in range(7):
            psi = step(psi, law)
        worst = max(worst, abs(np.linalg.norm(psi) - 1))
    assert worst < 1e-12, worst
    return f"norm preserved by all {len(LAWS)} members, worst drift {worst:.1e}"


def check_sector_reduction():
    rng = np.random.default_rng(RNG_SEED)
    law = LAWS["reflected-routing"]
    size, private, order = 6, 8, 3
    chi = _random_state((size, size, size, 8), rng)
    phase = np.exp(2j * np.pi * order * np.arange(private) / private)[None, None, None, :, None]
    full = chi[:, :, :, None, :] * phase
    for _ in range(9):
        chi = sector_step(chi, law, 2 * np.pi * order / private)
        full = step(full, law)
    error = float(np.max(np.abs(full - chi[:, :, :, None, :] * phase)))
    assert error < 1e-12, error
    return f"private-sector update equals the literal update, max deviation {error:.1e}"


def check_field_identity():
    rng = np.random.default_rng(RNG_SEED)
    equal = np.ones(8) / np.sqrt(8)
    worst = 0.0
    for _ in range(20):
        angles = rng.uniform(-np.pi, np.pi, 4)
        energy, walk, state = dispersion(angles), walk_matrix(angles), equal.astype(complex)
        for n in range(1, 16):
            state = walk @ state
            worst = max(worst, abs(np.vdot(equal, state) - np.cos(n * energy)))
    assert worst < 1e-12, worst
    return f"<e|W^n e> = cos(nE) on 20 random wave vectors, max deviation {worst:.1e}"


def check_encounter():
    rng = np.random.default_rng(RNG_SEED)
    carrier = (4, 3, 1, 2)
    a = _random_state(carrier + (8,), rng)
    b = _random_state(carrier + (8,), rng)
    pair = np.einsum("ijklm,nopqr->ijklmnopqr", a, b)
    free = LAWS["reflected-routing"]
    met = LAWS["reflected-routing-encounter"]
    product, interacting = pair.copy(), pair.copy()
    for _ in range(5):
        a, b = step(a, free), step(b, free)
        product = step_pair(product, free)
        interacting = step_pair(interacting, met)
    independent = float(np.max(np.abs(product - np.einsum("ijklm,nopqr->ijklmnopqr", a, b))))
    differs = float(np.max(np.abs(interacting - product)))
    norm = abs(np.linalg.norm(interacting) - 1)
    assert independent < 1e-12 and differs > 1e-3 and norm < 1e-12
    return (f"zero encounter phase evolves two excitations independently ({independent:.1e}); "
            f"phase pi/2 changes the pair state by up to {differs:.3f} and keeps its norm")


def check_spectral_engine():
    """Port of the step-3 gate audit: spectral readings against the literal update."""
    law = LAWS["reflected-routing"]
    size, origin, sites = 8, (2, 4, 3), [(2, 4, 3), (5, 6, 7)]
    worst = 0.0
    for kind in ("optical", "clock", "reference"):
        packet = Packet(kind, (1, -1), size=size, private_size=size, carrier_order=1,
                        private_order=1, width=0.8, support=1, transverse_sign=-1, origin=origin)
        fourier = np.zeros((size,) * 4 + (8,), complex)
        fourier[np.ix_(*(o % size for o in packet.orders))] = packet.spinors
        psi = np.roll(np.fft.ifftn(fourier, axes=(0, 1, 2, 3), norm="ortho"), origin, axis=(0, 1, 2))
        for t in range(6):
            density = np.sum(np.abs(psi) ** 2, axis=(3, 4))
            for axis in range(3):
                projected = density.sum(axis=tuple(a for a in range(3) if a != axis))
                worst = max(worst, float(np.max(np.abs(projected - packet.marginal(t, axis)))))
            if kind == "clock":
                for site in sites:
                    local = np.sum(np.abs(psi[site]) ** 2, axis=-1)[:: size // 4]
                    worst = max(worst, float(np.max(np.abs(local - packet.ramsey(t, site)))))
            psi = step(psi, law)
    assert worst < 1e-12, worst
    return f"spectral marginals and four-aperture readings equal the literal update, max deviation {worst:.1e}"


def check_continuum():
    coarse = continuum_run(16, 2.0, keep_frames=False)["max_error"]
    fine = continuum_run(32, 2.0, keep_frames=False)["max_error"]
    ratio = coarse / fine
    assert 3.6 < ratio < 4.4 and fine < 0.012, (coarse, fine)
    return (f"largest deviation from Klein-Gordon {coarse:.4f} at N=16 and {fine:.4f} at N=32, "
            f"ratio {ratio:.2f} (second order)")


def check_axial_journey():
    j = clock_journey((1, 0))
    assert j["reunion"] == 1459 and abs(j["forecast"] - 0.71128) < 2e-4
    assert abs(j["beat_ratio"] - j["forecast"]) < 0.035 and abs(j["optical_ratio"] - j["forecast"]) < 0.035
    return (f"axial step-3 journey: reunion {j['reunion']}, forecast {j['forecast']:.5f}, "
            f"beat {j['beat_ratio']:.5f}, optical {j['optical_ratio']:.5f}")


CHECKS = [
    ("unitarity", check_unitarity),
    ("sector reduction", check_sector_reduction),
    ("field identity", check_field_identity),
    ("encounter phase", check_encounter),
    ("spectral engine", check_spectral_engine),
    ("continuum limit", check_continuum),
    ("two-clock journey", check_axial_journey),
]


def test_unitarity():
    check_unitarity()


def test_sector_reduction():
    check_sector_reduction()


def test_field_identity():
    check_field_identity()


def test_encounter():
    check_encounter()


def test_spectral_engine():
    check_spectral_engine()


def test_continuum():
    check_continuum()


def test_axial_journey():
    check_axial_journey()
