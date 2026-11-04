"""Full-Fock Lindblad controls at a fixed nonzero repair rate.

These dimer checks test the source algebra and limiting mechanism; they
are not evaluations at the conservative cubic theorem threshold.
"""
import numpy as np
import pytest
from scipy.linalg import expm

from . import algebra as a


def source_generator(c, U, rate, law):
    H = U/2*np.diag(a.defects(2))+c/2*a.hopping_twice(2, ((0, 1, -1),))
    identity = np.eye(16)
    repair = sum(np.kron(K, K)/2 for K in a.branch_matrices(2, law, (0, 1)))
    generator = -1j*(np.kron(identity, H)-np.kron(H.T, identity))+rate*(repair-np.eye(256))
    return H, generator


def spin_hamiltonian():
    X = np.array([[0., 1], [1, 0]])
    Y = np.array([[0., -1j], [1j, 0]])
    Z = np.diag([1., -1])
    return (np.kron(X, X)+np.kron(Y, Y)+np.kron(Z, Z)-np.eye(4))/4


@pytest.mark.parametrize('law', ['C', 'E'])
def test_full_source_at_fixed_repair_rate_converges_to_nontrivial_spin_motion(law):
    initial = np.zeros((16, 16), complex)
    initial[9, 9] = 1
    spin = np.zeros((4, 4), complex)
    spin[2, 2] = 1
    Hspin = spin_hamiltonian()
    time, rate = .7, .8
    expected = expm(-1j*time*Hspin)@spin@expm(1j*time*Hspin)
    readout = np.diag([1., -1, 1, -1])
    expected_value = np.trace(expected@readout).real
    assert abs(expected_value-np.trace(spin@readout).real) > .2
    errors, charges = [], []
    for c in (4, 8, 16, 32):
        H, generator = source_generator(c, c*c, rate, law)
        evolved = (expm(time*generator)@initial.reshape(-1, order='F')).reshape((16, 16), order='F')
        assert abs(np.trace(evolved)-1) < 1e-10
        assert np.linalg.eigvalsh(evolved).min() > -1e-11
        errors.append(abs(np.trace(evolved@a.spin_embedding(2, readout)).real-expected_value))
        charges.append(np.trace(evolved@np.diag(a.defects(2))).real/2)
    assert errors[-1] < errors[0]/10
    assert charges[-1] < charges[0]/10
    assert errors[-1] < .01


@pytest.mark.parametrize('law', ['C', 'E'])
def test_energy_derivative_and_charge_conservation_follow_actual_generator(law):
    rng = np.random.default_rng(118)
    M = rng.normal(size=(16, 16))+1j*rng.normal(size=(16, 16))
    density = M@M.conj().T
    density /= np.trace(density)
    H, generator = source_generator(3, 17, .7, law)
    derivative = (generator@density.reshape(-1, order='F')).reshape((16, 16), order='F')
    N = np.diag(a.defects(2))
    repair_energy = a.twice_unread_dual(H, 2, law, (0, 1))/2-H
    assert abs(np.trace(derivative@H)-.7*np.trace(density@repair_energy)) < 1e-10
    assert np.array_equal(a.twice_unread_dual(N, 2, law, (0, 1)), 2*N)
    assert np.linalg.norm(generator@np.eye(16).reshape(-1, order='F')) < 1e-12


def test_faithful_full_source_input_keeps_a_positive_floor_under_actual_process():
    c, U = 8, 64
    delta = (2*c/U)**2
    density = delta*np.eye(16)/16
    density[9, 9] += 1-delta
    H, generator = source_generator(c, U, 1.3, 'C')
    assert abs(np.trace(density@H)/2-c*c/U) < 1e-12
    evolved = (expm(.4*generator)@density.reshape(-1, order='F')).reshape((16, 16), order='F')
    assert np.linalg.eigvalsh(evolved).min() >= delta/16-1e-12
    assert np.linalg.norm(evolved-density) > .01
