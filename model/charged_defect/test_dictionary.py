"""Same-candidate energy, complete full-Fock instrument and leakage controls."""
from itertools import product
import math

import numpy as np
import pytest
from scipy import sparse

from ..particle_transport import algebra
from . import dictionary as d


def independent_local_generator(axis):
    # Tensor Jordan--Wigner construction, independent of the field builder.
    one, z, low = np.eye(2), np.diag([1, -1]), np.array([[0, 1], [0, 0]])
    operators = []
    for mode in range(4):
        value = np.ones((1, 1))
        for j in reversed(range(4)):
            value = np.kron(value, low if j == mode else z if j < mode else one)
        operators.append(value)
    sigma = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
             np.diag([1, -1]))[axis]
    result = np.zeros((16, 16), complex)
    for a, b in product(range(2), repeat=2):
        term = (1j/16)*sigma[a, b]*operators[a].conj().T@operators[2+b]
        result += term+term.conj().T
    n = [x.conj().T@x for x in operators]
    return result+n[0]@n[1]+n[2]@n[3]-sum(n)/2


@pytest.mark.parametrize('kind,charge', [('hole', -1), ('particle', 1)])
@pytest.mark.parametrize('axis', range(3))
def test_identical_gap_and_charge_on_every_pair_spin_column(kind, charge, axis):
    k = independent_local_generator(axis)
    v, g = algebra.candidate_columns(kind, axis), algebra.ground_column(axis)
    delta = math.sqrt(17)/8-1/16
    ground_energy = -(4+math.sqrt(17))/8
    np.testing.assert_allclose(k@g, ground_energy*g, atol=1e-15, rtol=0)
    np.testing.assert_allclose(k@v, (ground_energy+delta)*v, atol=1e-15, rtol=0)
    embedding = d.embedded_columns(3, kind, axis)
    np.testing.assert_allclose(embedding.conj().T@embedding, np.eye(12), atol=1e-15, rtol=0)
    generator = d.tensor_generator(3, axis)
    np.testing.assert_allclose(generator@embedding,
                               (3*ground_energy+2*delta)*embedding, atol=1e-15, rtol=0)
    number = np.array([i.bit_count() for i in range(4096)])
    np.testing.assert_array_equal(embedding[number != 6+2*charge], 0)
    # Every omitted block is exactly in g, including reference entanglement
    # among the two internal doublets.  First pair is blocks 0 and 1.
    rng = np.random.default_rng(1143)
    spin_reference = rng.normal(size=(4, 3))+1j*rng.normal(size=(4, 3))
    spin_reference /= np.linalg.norm(spin_reference)
    state = (embedding[:, :4]@spin_reference).reshape(16, 256, 3)
    exterior = np.einsum('air,bir->ab', state, state.conj())
    np.testing.assert_allclose(exterior, np.outer(g, g.conj()), atol=1e-15, rtol=0)


@pytest.mark.parametrize('kind', ('hole', 'particle'))
@pytest.mark.parametrize('law', ('C', 'E'))
@pytest.mark.parametrize('seam', ((0, 1), (0, 2), (1, 2)))
def test_full_fock_extension_retains_outside_sector_and_energy(kind, law, seam):
    outcomes = d.full_extension(kind=kind, law=law, seam=seam)
    generator = d.tensor_generator(3)
    complete = sum((x.conj().T@x for _, x in outcomes),
                   start=sparse.csr_matrix((4096, 4096)))
    assert sparse.linalg.norm(complete-sparse.eye(4096)) < 1e-12
    assert max(sparse.linalg.norm(generator@x-x@generator) for _, x in outcomes) < 1e-12
    embedding = sparse.csr_matrix(d.embedded_columns(3, kind))
    outside = outcomes[-1][1]
    assert outcomes[-1][0] == ('outside',)
    assert sparse.linalg.norm(outside@embedding) < 1e-13
    # Empty input is outside this charged-pair sector and still has a full
    # output with probability one.  Removing this outcome fails completeness.
    vacuum = np.zeros(4096); vacuum[0] = 1
    np.testing.assert_array_equal(outside@vacuum, vacuum)
    assert all(np.linalg.norm(x@vacuum) == 0 for _, x in outcomes[:-1])
    assert sparse.linalg.norm(complete-outside.conj().T@outside-sparse.eye(4096)) > 60
    # Every selected native branch stays in the same candidate carrier.
    for _, x in outcomes[:-1]:
        assert sparse.linalg.norm(outside@x@embedding) < 1e-13


def test_same_protected_triplet_is_degenerate_and_uniform_field_destroys_protection():
    receipt = d.report()
    assert receipt['two_defect_carrier_energy_spread'] == 0
    assert receipt['new_absolute_physical_masses_predicted'] == 0
    assert receipt['triplet_is_three_generations'] is False
    for row in receipt['joined_uniform_field_outside_sector_leakage'].values():
        assert row['exact_squared_Frobenius_norm'] == '3/256'
        assert row['Frobenius_norm'] == pytest.approx(math.sqrt(3)/16, abs=1e-15)
        assert row['triplet_normalized_squared_leakage'] == '1/512'
        assert row['nonzero_nonleaking_input_exists'] is False
    # Independent normalization of the four internal columns; the protected
    # triplet has strictly positive joining leakage, not just the singlet.
    from ..full_fock_field.field import FieldGeometry, FullFock
    from .branches import COUPLING
    field = FullFock(8).hamiltonian(FieldGeometry((4, 1, 1), speed=1/8, interaction=1.))
    coupling = COUPLING/np.sqrt([2., 1., 2., 1.])
    for kind in ('hole', 'particle'):
        embedding = d.embedded_columns(2, kind)
        image = field@embedding
        leak = image-embedding@(embedding.conj().T@image)
        gram = coupling.conj().T@leak.conj().T@leak@coupling
        np.testing.assert_allclose(gram, np.diag([3, 1, 1, 1])/512, atol=1e-17, rtol=0)


@pytest.mark.parametrize('kwargs', [{'blocks': 12}, {'blocks': True}, {'kind': 'pair_hole'},
                                  {'axis': True}])
def test_invalid_embedding_domains_reject(kwargs):
    with pytest.raises(ValueError):
        d.embedded_columns(**kwargs)


@pytest.mark.parametrize('bad', ['3', 3.0, True, None, 12])
def test_invalid_full_extension_block_count_rejects_cleanly(bad):
    with pytest.raises(ValueError):
        d.full_extension(blocks=bad)
