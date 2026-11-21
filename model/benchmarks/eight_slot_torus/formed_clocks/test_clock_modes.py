"""Independent source, polarization and finite-window checks for the formed-pair read."""
import json

import numpy as np
import pytest

from . import clock_modes as modes
from .. import law, pair_persistence, spectral


def test_fibre_separation_endpoints_and_both_poles():
    for k, ranks in ((np.zeros(4), (1, 3, 4)), (np.full(4, np.pi), (1, 4, 3)),
                     (np.array([0.2, -0.7, 0.9, 1.4]), (2, 3, 3)),
                     (np.array([0., 0., 1e-10, 0.]), (2, 3, 3))):
        routing, bright, (even, odd) = modes.contact_fibre(k)
        projectors = (bright, even, odd)
        np.testing.assert_allclose(sum(projectors), np.eye(8), atol=2e-14)
        for projector, rank in zip(projectors, ranks):
            np.testing.assert_allclose(projector @ projector, projector, atol=2e-14)
            np.testing.assert_allclose(projector.conj().T, projector, atol=2e-14)
            assert np.trace(projector).real == pytest.approx(rank)
        np.testing.assert_allclose((even + odd) @ law.EQUAL, 0, atol=2e-14)
        np.testing.assert_allclose(routing @ even, even, atol=2e-14)
        np.testing.assert_allclose(routing @ odd, -odd, atol=2e-14)
        for rail, phi in enumerate(law.LAWS["torus"].rail_contact_phases):
            # Independent source matrix uses the literal simulator's free coin/routing.
            op = np.exp(0.5j * phi) * law.walk_matrix(k, np.exp(-1j * phi))
            E, minus, plus = modes.mature_born_poles(k, rail)
            np.testing.assert_allclose(op @ minus, np.exp(-1j * E) * minus, atol=2e-14)
            np.testing.assert_allclose(op @ plus, np.exp(1j * E) * plus, atol=2e-14)
            np.testing.assert_allclose(plus @ minus, 0, atol=2e-14)
            np.testing.assert_allclose(op @ even, -np.exp(0.5j * phi) * even, atol=2e-14)
            np.testing.assert_allclose(op @ odd, np.exp(0.5j * phi) * odd, atol=2e-14)
            for age in range(128):
                theta = law.hop_angle(age, law.LAWS["torus"].maturation)
                hop = np.exp(-1j * theta) * (np.cos(theta) * np.eye(8) + 1j * np.sin(theta) * routing)
                child = hop @ law.contact_coin(np.eye(8), phi) @ bright
                np.testing.assert_allclose((even + odd) @ child, 0, atol=2e-14)


@pytest.mark.parametrize("size", (17, 32))
@pytest.mark.parametrize("rail", (0, 1))
def test_sharp_birth_aperture_against_literal_each_age(size, rail):
    state = np.zeros((1, 1, size, 1, law.SLOTS), complex)
    state[0, 0, 0, 0] = law.EQUAL
    for t in range(194):
        if t in (0, 1, 63, 127, 128, 193):
            h = np.vdot(law.EQUAL, state[0, 0, :2, 0].sum(axis=0)) / np.sqrt(2)
            assert abs(modes.born_aperture(size, t, rail) - h) <= 2e-12
        state = law.step_maturing(state, law.LAWS["torus"], t, rail)


def test_held_port_against_literal_opening_and_scalar_recurrence():
    # The small real-space carrier contains the entire light cone through 129 steps.
    # Obtain both seeds of the mature scalar recurrence by applying the literal law.
    protocol = json.loads(modes.PROTOCOL.read_text(encoding="utf-8"))
    waits = [n + shift for n in protocol["mature_waits"] for shift in (0, protocol["offset"])]
    offsets = np.arange(-129, 130)
    for rail, phi in enumerate(law.LAWS["torus"].rail_contact_phases):
        state = np.zeros((1, 1, len(offsets), 1, 8), complex)
        state[0, 0, 129, 0] = law.EQUAL
        for age in range(128):
            state = law.step_maturing(state, law.LAWS["torus"], age, rail)
        next_state = law.step_bound(state, phi)
        seeds = [field[0, 0, :, 0] @ law.EQUAL for field in (state, next_state)]
        for n in waits:
            size = 2 * (128 + n) + 3
            K = 2 * np.pi * np.fft.fftfreq(size)
            pad = np.zeros((2, size), complex)
            pad[:, offsets % size] = seeds
            q, uq = np.fft.fft(pad, axis=1)
            E = np.arccos((3 + np.cos(K)) * np.cos(phi / 2) / 4)
            evolved = (np.sin(n * E) * uq - np.sin((n - 1) * E) * q) / np.sin(E)
            actual = np.fft.ifft(evolved)[:2].sum() / np.sqrt(2)
            assert abs(actual - modes.born_aperture(size, 128 + n, rail)) <= 2e-12


def test_compact_plaquette_signed_control_and_bright_capture_exclusion():
    even = pair_persistence.circulation()
    x, y = np.indices(even.shape[:2])
    odd = even * ((-1.) ** (x + y))[:, :, None, None, None]
    site = (3, 3, 0, 0)
    delta = np.diff(law.LAWS["torus"].rail_contact_phases)[0] / 2
    for field, parity in ((even, 1), (odd, -1)):
        np.testing.assert_allclose(law.route(field), parity * field)
        assert pair_persistence.capture_projection(field) < 2e-14
    for second, flip in ((even, 0), (odd, np.pi)):
        # Both parities can interfere at this one corner although they are globally orthogonal.
        before = spectral.rail_phase(np.array([even[site], second[site]]))
        after = spectral.rail_phase(np.array([law.step_maturing(field, law.LAWS["torus"], 128, rail)[site]
                                              for rail, field in enumerate((even, second))]))
        assert abs(np.angle(np.exp(1j * (after - before - flip + delta)))) < 2e-14
    assert abs(np.vdot(even, odd)) < 2e-14
    expected = -delta
    for item in modes.flat_control()["one_update"]:
        assert abs(item["hand_increment"] - expected) < 2e-14
        assert sum(item["rail_ports"]) == pytest.approx(0.25)
        assert item["miss"] == pytest.approx(.75)
    with pytest.raises(ValueError, match="pre-wrap"):
        modes.born_aperture(32, 8192, 0)



def test_original_start_records_agree_with_held_replay():
    expected = json.loads(modes.PROTOCOL.with_name("clock_modes_records.json").read_text(encoding="utf-8"))
    actual = json.loads(json.dumps(modes.station_records(), allow_nan=False))

    def compare(observed, sealed, field=None):
        assert type(observed) is type(sealed)
        if isinstance(observed, dict):
            assert observed.keys() == sealed.keys()
            for key in sealed:
                compare(observed[key], sealed[key], key)
        elif isinstance(observed, list):
            assert len(observed) == len(sealed)
            for left, right in zip(observed, sealed, strict=True):
                compare(left, right, field)
        elif isinstance(observed, float):
            if field not in ("relative_errors", "final_advance_error"):
                assert observed == pytest.approx(sealed, abs=2e-13, rel=1e-10)
        else:
            assert observed == sealed

    compare(actual, expected)
    for record in (actual, expected):
        final_advance = record["pairs"][-1]["advance"]
        difference = abs(final_advance - record["advance_forecast"])
        assert abs(record["final_advance_error"] - difference) < 3e-16
        assert record["final_phase_gate"] == (final_advance > 0 and difference <= 0.03)
        for run in record["pairs"]:
            for obs in run["appointments"]:
                bins = [obs["birth_refusal"], obs["born_port_miss"], *obs["rail_ports"]]
                assert min(bins) >= 0 and abs(sum(bins) - 1) < 1e-14
                # Check each ratio against its own full-precision components.
                # A small absolute perturbation of a complex amplitude is
                # amplified on division by a low-admission forecast.
                ratios = [abs(complex(*h)-complex(*l))/abs(complex(*l))
                          for h,l in zip(obs["conditional_rail_amplitudes"],
                                         obs["leading_amplitudes"], strict=True)]
                np.testing.assert_allclose(obs["relative_errors"], ratios, atol=3e-16, rtol=2e-14)
                assert obs["complex_gate"] == (max(ratios) <= run["tolerance"])
