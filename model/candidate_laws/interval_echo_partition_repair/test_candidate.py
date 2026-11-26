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
