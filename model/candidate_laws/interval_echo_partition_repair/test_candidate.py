"""Native contrast coupling and equilibrium with a retained reference system."""
import numpy as np

from . import candidate


def test_contrast_reads_and_unwrites_on_entangled_inputs():
    rng = np.random.default_rng(806)
    shape = (candidate.material.DIM, 5, 2)
    amplitudes = rng.normal(size=shape) + 1j*rng.normal(size=shape)
    amplitudes /= np.linalg.norm(amplitudes)
    after = candidate.contrast_aperture(amplitudes)
    np.testing.assert_allclose(candidate.contrast_aperture(after), amplitudes, rtol=0, atol=2e-15)
    old = amplitudes.reshape(-1, 2)
    new = after.reshape(-1, 2)
    np.testing.assert_allclose(new.conj().T @ new, old.conj().T @ old, rtol=0, atol=2e-15)
    rows = candidate.native_readout()["contrast_aperture"]
    assert len(rows) == 9 and sum(row["open"] for row in rows) == 6


def test_record_service_and_native_reservoir_act_locally():
    for volume, weights in [(False, np.full(5, 0.2)), (True, np.array([2/3, 1/3]))]:
        n = len(weights)
        remote = np.array([[0.5, 0.1j], [-0.1j, 0.5]])
        eq = np.einsum("ij,ab->iajb", np.diag(weights), remote)
        np.testing.assert_allclose(candidate.repair_record(eq, volume=volume), eq)
        psi = np.eye(n, 2)/np.sqrt(2)
        density = np.einsum("ia,jb->iajb", psi, psi.conj())
        serviced = candidate.repair_record(density, volume=volume)
        np.testing.assert_allclose(np.trace(serviced, axis1=0, axis2=2), np.eye(2)/2)
        expected = np.einsum("ij,ab->iajb", np.eye(n)/n, np.eye(2)/2)
        np.testing.assert_allclose(candidate.prepare_reference_payload(density), expected)


def test_native_gates_preserve_a_slowest_repair_mode():
    # Check the sharp witness with the full native contrast gate on all
    # seam/word columns, including coherent states of the seam memory.
    size = candidate.material.DIM
    mode = np.cos((2*np.pi/5) * (np.arange(5) - 0.5))
    mode /= max(abs(mode))
    weights = (1 + 0.5*mode)/5
    diagonal = np.tile(weights/size, size)
    joint = np.eye(size*5).reshape(size, 5, size*5)
    gate = candidate.contrast_aperture(joint).reshape(size*5, size*5)
    np.testing.assert_allclose(diagonal[:, None]*gate, gate*diagonal, atol=2e-17, rtol=0)

    after = candidate.repair_record(np.diag(weights).reshape(5, 1, 5, 1))
    deviation = np.diag(after[:, 0, :, 0]) - 1/5
    np.testing.assert_allclose(deviation, (3+np.sqrt(5))/8 * (weights-1/5),
                               atol=1e-16, rtol=0)


def test_complete_choi_matrices_and_selective_remote_states():
    # An unnormalized maximally entangled input tests the channel on every
    # matrix unit at once. The remote identity is its initial marginal.
    for size, volume in ((5, False), (2, True)):
        eye = np.eye(size)
        entangled = np.einsum("ia,jb->iajb", eye, eye)
        if volume:
            probabilities = np.array([[3, 1], [2, 2]])/4
        else:
            probabilities = (2*eye + np.roll(eye, 1, axis=1) + np.roll(eye, -1, axis=1))/4
        serviced = candidate.repair_record(entangled, volume=volume)
        joint = np.diag(probabilities.T.ravel()).reshape(entangled.shape)
        np.testing.assert_array_equal(serviced, joint)
        np.testing.assert_array_equal(np.trace(serviced, axis1=0, axis2=2), eye)
        selected = serviced[0, :, 0, :]
        assert not np.allclose(selected/np.trace(selected), eye/size)
        reset = candidate.prepare_reference_payload(entangled)
        np.testing.assert_array_equal(reset, np.eye(size**2).reshape(entangled.shape)/size)
