import numpy as np
import pytest
import sympy as s

from model.cube_eta_protection import cube
from model.interacting_vacuum import ground as g


def test_source_negative_band_is_exact_and_full_spin_doublet_is_retained():
    a, b, even, odd, kr, ks = g.band_data()
    assert sorted(even+odd) == list(range(8))
    assert np.array_equal(a@a, 3*np.eye(8, dtype=np.int64))
    assert np.array_equal(b@b.T, 3*np.eye(4, dtype=np.int64))
    assert np.array_equal(kr.T@kr+3*ks.T@ks, 6*np.eye(4, dtype=np.int64))
    assert np.array_equal(a@kr, -3*ks)
    assert np.array_equal(a@ks, -kr)


def test_full_contact_rank_and_balanced_block_use_independent_rational_elimination():
    c = g.contact_matrix()
    balanced = [j for j, mask in enumerate(g.BASIS4) if (mask & 0x55).bit_count() == 2]
    # SymPy eliminates over Q; this independently confirms the prime minor.
    assert s.Matrix(c[:, balanced]).rank() == 35
    assert sum(row['rank_lower_bound'] for row in g.certificate()['spin_blocks']) == 65
    assert c.shape == (224, 70)


def test_each_retained_prime_minor_is_nonzero_as_an_exact_integer():
    c = g.contact_matrix()
    for record in g.certificate()['spin_blocks']:
        up = record['hidden_up_particles']
        columns = [j for j, mask in enumerate(g.BASIS4) if (mask & 0x55).bit_count() == up]
        block = c[:, columns]
        minor = block[np.ix_(record['minor_rows'], record['minor_columns_within_block'])]
        exact = int(s.Matrix(minor).det(method='domain-ge'))
        assert exact != 0
        assert exact % 101 == record['minor_determinant_mod101']


def test_omitting_odd_sublattice_contacts_enlarges_kernel():
    even = g.band_data()[2]
    reduced = g.contact_matrix(even)
    assert g.modular_rank(reduced)[0] == 54
    assert 70-g.modular_rank(reduced)[0] == 16
    assert g.modular_rank(g.contact_matrix())[0] == 65


def test_all_five_spin_two_weights_are_present_and_contact_free():
    vectors = g.multiplet()
    assert [int(v@v) for v in vectors] == [1, 4, 24, 144, 576]
    assert np.array_equal(np.column_stack(vectors).T@np.column_stack(vectors),
                          np.diag([1, 4, 24, 144, 576]))
    for k, vector in enumerate(vectors):
        assert not np.any(g.contact_matrix()@vector)
        assert all((mask & 0x55).bit_count() == 4-k
                   for mask, coefficient in zip(g.BASIS4, vector, strict=True) if coefficient)


def test_highest_weight_is_normalized_and_is_an_exact_physical_frame_eigenstate():
    a, b = g.reference_arrays()
    assert int(a@a+3*b@b) == 1296
    assert int(a@b) == 0
    ops = cube.operators()
    r = cube.restrict(ops.frame.real, 4)
    i = cube.restrict(ops.frame.imag, 4)
    ar, ai, br, bi = r.T@a, -i.T@a, r.T@b, -i.T@b
    hr, hi = cube.restrict(ops.H.real, 4), cube.restrict(ops.H.imag, 4)
    # Exact Q(i,sqrt3) equality 16H psi=-4sqrt3 psi.
    assert np.array_equal(hr@ar-hi@ai, -12*br)
    assert np.array_equal(hr@ai+hi@ar, -12*bi)
    assert np.array_equal(hr@br-hi@bi, -4*ar)
    assert np.array_equal(hr@bi+hi@br, -4*ai)
    assert int(ar@ar+ai@ai+3*(br@br+bi@bi)) == 1296


def test_reference_is_explicit_sparse_fixed_number_hidden_up_state():
    reference = g.reference()
    assert len(reference) == 62
    assert all(mask.bit_count() == 4 and mask & 0xAAAA == 0 for mask in reference)
    assert all(type(x) is int for pair in reference.values() for x in pair)
    reference.clear()
    assert len(g.reference()) == 62
    result = g.certificate()
    assert result['N4_ground_dimension'] == 5
    assert result['reference_prepared_by_native_source'] is False
    assert result['other_particle_number_ground_comparison_proved_here'] is False


@pytest.mark.parametrize('action', [
    lambda: g.modular_rank([[1.1]]),
    lambda: g.modular_rank([[1]], prime=2),
    lambda: g.contact_matrix((0, 0)),
    lambda: g.contact_matrix((8,)),
    lambda: g.spin_ladder(np.zeros(70), lower=True),
    lambda: g.spin_ladder(np.zeros(70, dtype=np.int64), lower='yes'),
])
def test_inexact_or_unadmitted_inputs_fail(action):
    with pytest.raises(ValueError): action()
