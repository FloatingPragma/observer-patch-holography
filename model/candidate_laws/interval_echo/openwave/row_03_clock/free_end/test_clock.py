"""Independent checks use dense native wing matrices and the committed integer records."""

import numpy as np
import pytest
from scipy.linalg import eigh, expm

from ...row_08_annihilation import experiment
from ...row_08_annihilation.daughter_modes.measure import frequencies
from . import measure as clock


@pytest.mark.parametrize("length", [1, 3, 5])
@pytest.mark.parametrize("pin", [0, 2])
@pytest.mark.parametrize("contrast", [4, 8])
def test_free_tip_clifford_correlation(length, pin, contrast):
    h = experiment.wing(contrast, length, pin)
    tip = np.kron(np.eye(2**(length-1)), [[0,-1j],[1j,0]])
    k = np.pi*np.arange(1,length+1)/(length+1)
    weights = 2*np.sin(k)**2/(length+1)
    lines = frequencies(length, contrast, False)
    energy, frame = eigh(h)
    one = frame[:, 1:1+length]
    # Single-occupancy states lie below two-occupancy threshold in these fixtures.
    np.testing.assert_allclose(energy[1:length+1]-energy[0], lines, atol=2e-12)
    np.testing.assert_allclose(one.T @ tip @ one, 0, atol=3e-11)
    for time in (0, .032, 1.3, 22.45):
        u = expm(-1j*time*(h-np.eye(len(h))*np.mean(h.diagonal())))
        yt = u.conj().T @ tip @ u
        expected = weights @ np.cos(lines*time)
        np.testing.assert_allclose((tip@yt+yt@tip)/2, np.eye(len(h))*expected, atol=2e-11)
        np.testing.assert_allclose(one.T @ yt @ one, 0, atol=4e-11)


def test_signed_load_read_pulse_and_return():
    x = np.kron(np.eye(4), [[0,1],[1,0]])
    z = np.kron(np.eye(4), [[1,0],[0,-1]])
    pulse = expm(1j*np.pi*x/4)
    # Orientation -Z after the positive quarter-turn gives the declared Y.
    np.testing.assert_allclose(pulse.conj().T @ (-z) @ pulse, clock.CURRENT, atol=1e-14)
    for s, p in zip(clock.SIGNS, clock.PROJECTORS, strict=True):
        np.testing.assert_allclose(pulse.conj().T @ ((clock.IDENTITY-s*z)/2) @ pulse,
                                   p, atol=1e-14)


def test_complete_records_and_held_predictions():
    with np.load(clock.HERE/"counts.npz", allow_pickle=False) as receipts:
        clock.validate(clock.decode(receipts))


def test_foreign_reference_can_only_supply_its_own_fast_phase():
    # Coherent superposition across two single-quantum energies. Coupling only at readout.
    e, v = eigh(experiment.wing(8, 3, 0))
    target = (v[:,1]+1j*v[:,2])/np.sqrt(2)
    target_dipole = np.outer(v[:,0],v[:,1])+np.outer(v[:,1],v[:,0])
    reference = np.array([1,1])/np.sqrt(2)
    x = np.array([[0,1],[1,0]])
    for t in (0,.4,1.7):
        target_now = v @ (np.exp(-1j*(e-e[0])*t)*(v.T @ target))
        reference_now = np.array([1,np.exp(-31.5j*t)])*reference
        joint = np.kron(target_now,reference_now)
        cross = np.kron(target_dipole,x)
        np.testing.assert_allclose(np.vdot(joint,cross @ joint), 0, atol=4e-13)
        np.testing.assert_allclose(np.vdot(joint,np.kron(clock.IDENTITY,x) @ joint),
                                   np.cos(31.5*t),atol=4e-13)
