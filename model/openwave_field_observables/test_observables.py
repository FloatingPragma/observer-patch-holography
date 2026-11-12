"""Independent small-matrix oracles and adversarial source changes."""
from itertools import product
import numpy as np
import pytest
from scipy.linalg import expm
from . import lattice as l, waves, electrostatic, mass_shell


def incidence(n):
    sites = tuple(product(range(n), repeat=3)); lookup = {x: i for i, x in enumerate(sites)}
    d = np.zeros((3*n**3, n**3))
    c = np.zeros((3*n**3, 3*n**3))
    for v, x in enumerate(sites):
        for j in range(3):
            y = tuple((a+int(k == j)) % n for k, a in enumerate(x))
            d[3*v+j, v] -= 1; d[3*v+j, lookup[y]] += 1
        for f, (i, j) in enumerate(l.PAIRS):
            xi = tuple((a+int(k == i)) % n for k, a in enumerate(x))
            xj = tuple((a+int(k == j)) % n for k, a in enumerate(x))
            c[3*v+f, 3*v+i] += 1
            c[3*v+f, 3*lookup[xi]+j] += 1
            c[3*v+f, 3*lookup[xj]+i] -= 1
            c[3*v+f, 3*v+j] -= 1
    return d, c


def test_real_space_evolution_against_independent_dense_generator():
    n = 3; d, c = incidence(n)
    assert np.array_equal(c@d, np.zeros((81, 27)))
    rng = np.random.default_rng(234)
    e = c.T@rng.normal(size=81)
    b = c@rng.normal(size=81)
    generator = np.block([[np.zeros((81, 81)), c.T], [-c, np.zeros((81, 81))]])
    reference = expm(.2*generator)@np.r_[e, b]
    got_e, got_b = l.verlet(e.reshape(n,n,n,3), b.reshape(n,n,n,3), 1, .2, 128)
    assert np.linalg.norm(np.r_[got_e.ravel(), got_b.ravel()]-reference) < 5e-5
    assert np.linalg.norm(d.T@got_e.ravel()) < 1e-11


def test_source_green_against_dense_pseudoinverse():
    d, _ = incidence(3)
    lap = d.T@d
    expected = np.linalg.pinv(lap, hermitian=True)[:, 0].reshape(3,3,3)
    assert np.max(np.abs(l.green(3)-expected)) < 1e-13


def test_orientation_mutation_breaks_chain_identity():
    d, c = incidence(3)
    c[0, 0] *= -1
    assert np.linalg.norm(c@d) > 1


@pytest.mark.parametrize('cap', [1, 2, 3])
def test_complete_quantum_boundary(cap):
    value = waves.cutoff_theorem(cap)
    assert value['complete_dimension'] == (cap+1)*(cap+2)//2
    assert value['omitted_boundary_completion_detected']


def test_unstable_step_is_rejected():
    e = np.zeros((3,3,3,3))
    with pytest.raises(ValueError):
        l.verlet(e, e, 1, 1, 1)


def test_gauss_failure_cannot_be_hidden_by_curl():
    phi = np.random.default_rng(221).normal(size=(4,4,4))
    e = l.gradient(phi)
    assert np.linalg.norm(l.curl(e)) < 1e-12
    assert np.linalg.norm(l.divergence_adjoint(e)) > 1


def test_invalid_coupling_rejected():
    with pytest.raises(ValueError):
        electrostatic.pair_energy(l.green(4), (1,0,0), 0)


def test_unreported_energy_shift_breaks_complete_wave_identity():
    h = mass_shell.hamiltonian((0,0,0), 1)
    assert np.max(np.abs(np.sort(np.linalg.eigvalsh(h))-np.repeat([-5,-1,1,5], 4))) < 1e-12
    shifted = h+.1*np.eye(16)
    assert np.linalg.norm(shifted@shifted-h@h) > 1


def test_parity_reduction_intertwines_complete_small_source():
    from scipy import sparse
    from . import orbitals
    n = 5
    lap, potential, a = orbitals.dirichlet_source(n, 12.)
    eye = sparse.eye(n)
    full = (sparse.kron(sparse.kron(lap, eye), eye)
            +sparse.kron(sparse.kron(eye, lap), eye)
            +sparse.kron(sparse.kron(eye, eye), lap))/(2*a*a)
    full = full+sparse.diags(potential.ravel())
    for parity in ((False,False,False),(True,False,False)):
        frames = [orbitals.parity_frame(n, odd)[0] for odd in parity]
        isometry = sparse.kron(sparse.kron(frames[0],frames[1]),frames[2])
        reduced, _ = orbitals.reduced_hamiltonian(n, 12., parity)
        assert np.max(np.abs((full@isometry-isometry@reduced).toarray())) < 1e-13
        defect = sparse.diags(np.arange(n**3)/n**3)
        assert np.linalg.norm(((full+defect)@isometry-isometry@reduced).toarray()) > .1


def test_invalid_parity_grid_rejected():
    from . import orbitals
    with pytest.raises(ValueError):
        orbitals.parity_frame(4)
