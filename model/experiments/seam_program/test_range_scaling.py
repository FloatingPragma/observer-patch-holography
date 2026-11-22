"""Independent microscopic and Fourier checks of the interaction range experiment."""

import json
from hashlib import sha256

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from .range_scaling import (
    ARMS,
    G,
    HERE,
    aperture,
    dicke_control,
    free_step,
    held_law,
    simulate,
)
from .range_contact import native_check


def test_native_occupation_phase_and_work_rail_restoration():
    basis = [native_law.SRC.index(x) for x in ((0, 0), (0, 1), (1, 0), (1, 1))]
    embed = np.eye(native_law.DIM)[:, basis]
    h0, h1 = native_law.local_term_parts()
    sign = np.diag([1, 1, -1, -1])
    for tau in (1 / 16, 0.39, -0.6):
        gate = expm(-1j * (tau / 2) * (h0 + KAPPA * h1)) @ embed
        u = embed.T @ gate
        np.testing.assert_allclose(gate, embed @ u, atol=1e-13)
        expected = np.diag(np.exp(-1j * tau * np.array([0, 8.5, 8.5, 0])))
        np.testing.assert_allclose(sign @ u @ sign @ u, expected, atol=1e-13)
    for active in (True, False):
        assert native_check(1 / 16, active)["maximum_column_error"] < 1e-10


def walk_symbol(k):
    coin = np.ones((8, 8)) / 4 - np.eye(8)
    directions = np.repeat(np.r_[k, 0.0], 2) * np.tile([1, -1], 4)
    swap = np.eye(8)[np.arange(8) ^ 1]
    return swap @ np.diag(np.exp(-1j * directions)) @ coin


def test_relative_shift_against_independent_one_body_symbols():
    rng = np.random.default_rng(921)
    state = rng.normal(size=(4, 4, 4, 8, 8)) + 1j * rng.normal(size=(4, 4, 4, 8, 8))
    modes = np.fft.fftn(state, axes=(0, 1, 2))
    for address in np.ndindex(4, 4, 4):
        k = 2 * np.pi * np.array(address) / 4
        modes[address] = walk_symbol(-k) @ modes[address] @ walk_symbol(k).T
    np.testing.assert_allclose(
        free_step(state), np.fft.ifftn(modes, axes=(0, 1, 2)), atol=1e-13
    )


def test_all_range_and_dose_controls_at_common_momenta_and_times():
    saved = json.loads((HERE / "range_scaling.json").read_text())
    assert (
        saved["protocol_sha256"]
        == sha256((HERE / "range_scaling_protocol.json").read_bytes()).hexdigest()
    )
    assert {(r["n"], r["arm"]) for r in saved["refinements"]} == {
        (n, a) for n in (4, 8, 16, 32) for a in ARMS
    }
    for row in saved["refinements"]:
        n, arm = row["n"], row["arm"]
        _, window, strength = aperture(n, arm)
        fraction = window.mean()
        prob = 1 - 4 * fraction * (1 - fraction) * np.sin(strength / 2) ** 2
        np.testing.assert_allclose(row["held_at_one"][0], prob, atol=3e-13)
        np.testing.assert_allclose(sum(row["clocked_at_one"]), 1, atol=1e-10)
        for record in row["fixed"].values():
            assert min(record) >= 0
            np.testing.assert_allclose(sum(record), 1, atol=1e-10)
        if n <= 8:
            fresh = simulate(n, arm)
            for t, record in row["fixed"].items():
                np.testing.assert_allclose(fresh["fixed"][t], record, atol=1e-11)
            np.testing.assert_allclose(
                fresh["clocked_at_one"], row["clocked_at_one"], atol=1e-11
            )
    for row in saved["symmetric_density"]:
        assert row == dicke_control(row["n"], row["filling"])
    assert 0 < held_law(aperture(8, "half_width")[1], G, 1)[1] < 1
