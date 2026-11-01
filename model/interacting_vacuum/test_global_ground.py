from fractions import Fraction as F

import numpy as np
import pytest
from scipy import sparse

from model.interacting_vacuum import global_ground as g


def test_all_four_full_balanced_blocks_certify_the_global_ground():
    result = g.certificate()
    assert result['global_ground_dimension'] == 5
    assert [x['balanced_dimension'] for x in result['Gram_blocks']] == [1568, 3136, 3920, 4900]
    assert all(x['residual_strict_upper'] == '1/512' for x in result['Gram_blocks'])
    assert result['all_N_not4_energy_gap_strict_lower'] == '7/1024'
    assert result['numerical_eigenvalue_admission'] is False


def test_exact_scalar_gram_does_not_require_cholesky_or_full_rank():
    h = sparse.csr_matrix(np.zeros((2, 2), dtype=np.int64))
    q = np.zeros((2, 2), dtype=np.int64)
    assert g.verify_gram(h, q, F(1)) == F(27, 64)
    q[0, 0] = 1 << 19
    assert g.verify_gram(h, q, F(1)) == F(27, 64)


def test_corrupted_factor_is_rejected_by_integer_residual():
    h = sparse.csr_matrix(np.eye(2, dtype=np.int64))
    q = g.propose_gram(h)
    assert g.verify_gram(h, q) < g.CEILING
    q[0, 0] += 10000
    with pytest.raises(ArithmeticError, match='residual'):
        g.verify_gram(h, q)


def test_duplicate_sparse_entries_are_canonicalized_before_exact_addition():
    h = sparse.csr_matrix((np.array([1, 2], dtype=np.int64), np.array([0, 0]), np.array([0, 2])), shape=(1, 1))
    canonical = sparse.csr_matrix([[3]], dtype=np.int64)
    q = g.propose_gram(canonical)
    assert g.verify_gram(h, q) == g.verify_gram(canonical, q)


def test_integer_overflow_and_nonhermitian_inputs_are_rejected():
    with pytest.raises(OverflowError):
        g.verify_gram(sparse.eye(2, dtype=np.int64), np.full((2, 2), 1 << 40, dtype=np.int64))
    with pytest.raises(ValueError, match='symmetric'):
        g.verify_gram(sparse.csr_matrix([[0, 1], [0, 0]], dtype=np.int64), np.zeros((2, 2), dtype=np.int64))
    with pytest.raises(ValueError):
        g.verify_gram(sparse.eye(2), np.zeros((2, 2), dtype=np.int64))


@pytest.mark.parametrize('value', [-1, 17, True, 4.0])
def test_particle_number_domain(value):
    with pytest.raises(ValueError): g.balanced_indices(value)
