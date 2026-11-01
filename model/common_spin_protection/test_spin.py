"""Independent same-field symmetry, complete-repair and frame controls."""
from itertools import product

import numpy as np
import pytest
from scipy import sparse

from ..full_fock_field.field import FullFock
from ..field_symmetry.symmetry import hidden_spin_matrix
from . import spin


def independent_hidden_spin():
    field = FullFock(12)
    c = [field.annihilator(j) for j in range(12)]
    sigma = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
             np.diag([1, -1]))
    output = []
    for axis in range(3):
        value = sparse.csr_matrix((4096, 4096), dtype=complex)
        for site in range(6):
            sign = 1 if axis == 0 else (-1)**site
            for a, b in product(range(2), repeat=2):
                value += sign*sigma[axis][a, b]*(c[2*site+a].conj().T@c[2*site+b])
        output.append(value)
    return output


def test_full_field_hidden_spin_from_independent_car_and_literal_staggering():
    spins = independent_hidden_spin()
    field = FullFock(12).hamiltonian(spin.GEOMETRY)
    for calculated, observed in zip(spins, spin.operators().spin_twice, strict=True):
        assert sparse.linalg.norm(calculated-observed.matrix()) == 0
        assert sparse.linalg.norm(field@calculated-calculated@field) == 0
    # Uniform bare spin is not the hidden conserved triplet on this chain.
    c = [FullFock(12).annihilator(j) for j in range(12)]
    bare_z = sum((c[2*j].conj().T@c[2*j]-c[2*j+1].conj().T@c[2*j+1]
                  for j in range(6)), start=sparse.csr_matrix((4096, 4096)))
    assert sparse.linalg.norm(field@bare_z-bare_z@field) > 1


def test_exact_local_premises_and_all_complete_full_fock_outcomes():
    local = spin.local_certificate()
    assert local['exact_charged_projector_commutators'] == 6
    assert local['exact_complete_256_state_swap_commutators'] == 3
    complete = spin.full_certificate()
    assert complete['complete_repair_outcomes_checked'] == 156
    assert complete['outcome_spin_commutators_checked'] == 468
    assert complete['largest_repair_spin_Frobenius_residual'] < 1e-12
    assert complete['outside_outcome_retained'] is True
    assert complete['joined_energy_conserved_by_repair'] is False
    assert complete['joined_energy_pullback_Frobenius_change'] > .4


def test_general_alignment_and_even_side_three_dimensional_route():
    certificate = spin.alignment_certificate()
    assert [row['displacement_parity'] for row in certificate['exact_parity_cases']
            if row['aligned']] == [[0, 0, 0], [1, 1, 1]]
    for displacement in product(range(-3, 4), repeat=3):
        expected = all((displacement[a]+displacement[b]) % 2 == 0
                       for a, b in ((0, 1), (1, 2), (0, 2)))
        assert spin.aligned_translation(displacement) == expected
    # Compare actual site frames on an even-sided cube and translated copies.
    for base in ((0, 0, 0), (2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 1)):
        for point in product(range(2), repeat=3):
            for axis in range(3):
                translated = tuple(a+b for a, b in zip(base, point, strict=True))
                np.testing.assert_array_equal(hidden_spin_matrix(translated, axis),
                                              hidden_spin_matrix(point, axis))
    assert certificate['new_cube_candidate_or_source_execution_supplied'] is False


def test_cube_dimer_swap_obstruction_is_not_a_claim_of_no_possible_algebra():
    result = spin.frame_obstruction()
    assert result['exact_commutator_Gram'] == [[64, 0, 0], [0, 32, 0], [0, 0, 32]]
    assert result['surviving_global_spin_generator_linear_dimension'] == 0
    assert result['all_possible_noiseless_algebras_excluded'] is False
    # Independently check the full one-body cube generators under the two
    # undressed whole-dimer permutations, before any representation reduction.
    sites = tuple(product(range(2), repeat=3))
    site_index = {site: j for j, site in enumerate(sites)}
    block_bases = ((0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1))
    surviving = {0, 1, 2}
    for target in (1, 2):
        permutation = np.eye(16, dtype=int)
        for x in range(2):
            a = site_index[(x, *block_bases[0][1:])]
            b = site_index[(x, *block_bases[target][1:])]
            for component in range(2):
                permutation[[2*a+component, 2*b+component]] = permutation[[2*b+component, 2*a+component]]
        for axis in range(3):
            s = sparse.block_diag([sparse.csr_matrix(hidden_spin_matrix(site, axis)) for site in sites]).toarray()
            if np.any(permutation@s != s@permutation):
                surviving.discard(axis)
    assert not surviving


@pytest.mark.parametrize('bad', [(0, 0), (True, 0, 0), [0, 0, 0], (0., 0, 0)])
def test_translation_domain_rejects_noninteger_or_ambiguous_inputs(bad):
    with pytest.raises(ValueError):
        spin.aligned_translation(bad)
