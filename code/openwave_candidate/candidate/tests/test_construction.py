import numpy as np
from oph_openwave.construction import Construction


def test_one_generator_cost_and_relativistic_square():
    c = Construction(3, ((0, 1), (1, 2), (2, 0)), 0.17)
    b, h = c.incidence, c.hamiltonian
    expected = np.block([[b @ b.T + c.mass**2 * np.eye(3), np.zeros((3, 3))],
                         [np.zeros((3, 3)), b.T @ b + c.mass**2 * np.eye(3)]])
    assert np.allclose(h @ h, expected)
    n = np.array([2, -3, 7])
    assert n @ b @ b.T @ n == sum((n[a] - n[d])**2 for a, d in c.seams)


def test_independent_loads_have_additive_eigenphase_and_no_signalling():
    c = Construction(2, ((0, 1),), 0.31)
    energies, modes = np.linalg.eigh(c.hamiltonian)
    a, d = modes[:, 0], modes[:, -1]
    t = 0.73
    ea, ed = c.evolve_independent((a, d), t)
    assert np.allclose(np.kron(ea, ed), np.exp(-1j*t*(energies[0]+energies[-1]))*np.kron(a,d))
    assert np.allclose(c.evolve_independent((a, a), t)[0], ea)


def test_literal_signed_records_preserve_exact_gauss_and_mean_law():
    c = Construction(3, ((0, 1), (1, 2), (2, 0)))
    n, e = (10**30 + 1, -7, 4), (0, 0, 0)
    initial = c.gauss_register(n, e)
    for step in range(100):
        n, e, _ = c.committed_event(n, e, step % 3, step % 2)
        assert c.gauss_register(n, e) == initial
    n = (3, -2, 7)
    outcomes = [c.committed_event(n, (0, 0, 0), seam, coin)[0]
                for seam in range(3) for coin in (0, 1)]
    assert np.allclose(np.mean(outcomes, axis=0), (np.eye(3)-c.incidence@c.incidence.T/6)@n)


def test_nonfinite_states_and_noninteger_graph_are_rejected():
    import pytest
    for ports, seams, mass in [(True, (), 0), (2.5, (), 0), (2, ((0, 1.0),), 0),
                               (2, ((False, 1),), 0), (2, ((0, 1),), float('nan')),
                               (2, ((0, 1),), True)]:
        with pytest.raises(ValueError):
            Construction(ports, seams, mass)
    c = Construction(2, ((0, 1),))
    with pytest.raises(ValueError):
        c.evolve([float('nan'), 0, 0], .1)
