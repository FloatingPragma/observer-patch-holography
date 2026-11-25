import numpy as np
import pytest

from model.experiments.interval_echo_lepton_register.pole_meter.measure import setup

from .aperture import compact_instrument


def _cell_projector():
    _, numerator, denominator, *_ = setup()
    return numerator[72:84, 72:84] / denominator


def test_donor_reset_and_repeatability_on_all_port_inputs():
    p = _cell_projector()
    result = compact_instrument(p, 6)
    yes, no = result["accepted"], result["refused"]
    identity = np.eye(12)
    assert result["flag_ports"] == [72, 73]
    assert [route[0] for route in result["routes"]] == result["flag_ports"]
    assert len({route[-1] for route in result["routes"]}) == 1
    assert result["rotation_count"] > 0
    assert result["max_donor_residue"] < 1e-12
    for actual, expected in [
        (yes, p),
        (no, identity - p),
        (yes.conj().T @ yes + no.conj().T @ no, identity),
        (yes @ yes, yes),
        (no @ no, no),
        (yes @ no, np.zeros_like(p)),
    ]:
        np.testing.assert_allclose(actual, expected, atol=2e-12, rtol=0)


def test_one_donor_read_preserves_the_entangled_reference():
    # This second geometry has nonzero support on every possible port.
    angles = 2 * np.pi * np.arange(12) / 12
    basis = np.column_stack((np.cos(angles), np.sin(angles))) / np.sqrt(6)
    p = basis @ basis.T
    coherent = compact_instrument(p, 6)
    resolved = compact_instrument(p, 6, resolve_address=True)
    assert resolved["accepted"] is None
    assert len(resolved["accepted_kraus"]) == 2
    bell = basis / np.sqrt(2)
    for instrument, expected_fidelity in [(coherent, 1), (resolved, 0.5)]:
        branches = [k @ bell for k in instrument["accepted_kraus"]]
        np.testing.assert_allclose(sum(np.vdot(b, b) for b in branches), 1, atol=1e-12)
        np.testing.assert_allclose(
            sum(abs(np.vdot(bell, b)) ** 2 for b in branches),
            expected_fidelity,
            atol=1e-12,
        )
        refused_bell = np.column_stack(
            (np.cos(2 * angles), np.sin(2 * angles))
        ) / np.sqrt(12)
        np.testing.assert_allclose(
            instrument["refused"] @ refused_bell, refused_bell, atol=1e-12
        )


def test_projector_and_layout_contract():
    p = _cell_projector()
    with pytest.raises(ValueError, match="workstation"):
        compact_instrument(p, 0)
    with pytest.raises(ValueError, match="rank two"):
        compact_instrument(np.eye(12), 6)
    unitary = np.eye(12, dtype=complex)
    pivot = np.unravel_index(np.argmax(abs(p - np.diag(p.diagonal()))), p.shape)[0]
    unitary[pivot, pivot] = 1j
    with pytest.raises(ValueError, match="real"):
        compact_instrument(unitary @ p @ unitary.conj().T, 6)
