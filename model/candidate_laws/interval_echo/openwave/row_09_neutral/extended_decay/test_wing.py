"""Independent checks of the guide energy and the two currents on the boundary."""

import json

import numpy as np
from scipy.linalg import expm

from ..joined_decay.measure import core_pulse
from . import wing


def test_direct_flights_and_isolated_calibration():
    receipt = json.loads((wing.HERE / "wing.json").read_text())
    wing.decode(receipt)
    _, _, settings = wing.installation()
    x = np.array([[0.0, 1.0], [1.0, 0.0]])
    exchange = -(np.kron(x, np.eye(2)) + np.kron(np.eye(2), x)) / 2
    source = exchange + np.diag([0.0, 33.0, 32.0, 33.0])
    local = exchange + np.diag([0.0, 0.0, 0.0, -1.0])
    gaps, frame = np.linalg.eigh(local)
    np.testing.assert_allclose(
        gaps - gaps[0], receipt["isolated_wing_gaps"], atol=1e-13
    )
    for value, line in zip(
        (gaps - gaps[0])[1:], receipt["clock_simulation"]["fringes"], strict=True
    ):
        assert abs(value - line["gap"]) < line["radius"]
    assert np.ptp(gaps) < receipt["parent_gap"] / 10

    # Closing the C,D triangle factors the physical generator across this cut.
    cut = settings["boundary_cut"]
    np.testing.assert_array_equal(
        cut, np.kron(cut[::4, ::4], np.eye(4)) + np.kron(np.eye(8), local)
    )
    i = np.arange(32)
    np.testing.assert_array_equal(
        cut - settings["native"], np.diag((1 - (i >> 2 & 1)) * (i >> 1 & 1))
    )
    _, source_frame = np.linalg.eigh(source)
    _, vectors = np.linalg.eigh(settings["native"])
    cold = vectors[:, 0]
    for name, h in settings.items():
        stocks = {
            "idle": cold,
            "pulse": np.kron(core_pulse(source_frame), np.eye(8)) @ cold,
        }
        for arm, stock in stocks.items():
            for wait, panel in zip(
                receipt["waits"], receipt["panels"][name][arm], strict=True
            ):
                # Direct matrix exponential and simultaneous projectors, independent
                # of the spectral evolution and compiled readout used to acquire.
                state = (expm(-1j * wait * h) @ stock).reshape(4, 2, 4)
                amplitudes = np.einsum("ai,de,acd->ice", source_frame, frame, state)
                np.testing.assert_allclose(
                    (abs(amplitudes) ** 2).sum(axis=1),
                    panel["core_by_wing"],
                    atol=2e-12,
                    rtol=0,
                )
                rails = abs(state.reshape(4, 2, 2, 2)) ** 2
                np.testing.assert_allclose(
                    rails.sum(axis=(0, 3)), panel["crossing_c_by_d"], atol=2e-12, rtol=0
                )
        if name != "native":
            for wait in receipt["waits"]:
                levels, energies = wing.contrast(receipt, name, wait)
                np.testing.assert_allclose(levels, 0, atol=1e-11)
                np.testing.assert_allclose(energies, 0, atol=1e-11)


def test_power_laws_for_short_and_long_clocks():
    _, _, settings = wing.installation()
    h = settings["native"]
    source = np.zeros(32, complex)
    source[8:12] = [1.0, -1j, 1j, -1.0]
    source /= np.linalg.norm(source)
    y = np.array([[0.0, -1j], [1j, 0.0]])
    projectors = ((np.eye(2) - y) / 2, (np.eye(2) + y) / 2)
    for stop, bins in ((0.1, 3), (32.0, 8)):
        panels = wing.boundary_power(h, source, stop=stop, bins=bins)
        for name, panel in panels.items():
            negative, positive = [
                np.kron(np.diag([1, 0]), eig)
                if name == "bare"
                else np.kron(eig, np.diag([0, 1]))
                for eig in projectors
            ]
            counts = np.zeros(3)
            for t in (np.arange(bins) + 0.5) * stop / bins:
                psi = expm(-1j * t * h) @ source
                for j, read in enumerate(
                    (negative, np.eye(4) - positive - negative, positive)
                ):
                    operator = np.kron(np.eye(4), np.kron(read, np.eye(2)))
                    counts[j] += np.vdot(psi, operator @ psi).real / bins
            np.testing.assert_allclose(
                counts, panel["negative_zero_positive"], atol=2e-12, rtol=0
            )
