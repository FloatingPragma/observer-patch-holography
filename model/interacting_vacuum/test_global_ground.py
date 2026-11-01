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


def test_public_certificate_mutations_cannot_poison_the_cached_proof():
    result = g.certificate()
    result['global_ground_dimension'] = 1
    result['Gram_blocks'][0]['exact_integer_Gram_verified'] = False
    result['Gram_blocks'][1]['residual_strict_upper'] = F(1, 2)
    result['Gram_blocks'][1]['residual_strict_upper']._denominator = -1
    result['Gram_blocks'].pop()
    fresh = g.certificate()
    assert fresh['global_ground_dimension'] == 5
    assert len(fresh['Gram_blocks']) == 4
    assert all(row['exact_integer_Gram_verified'] is True for row in fresh['Gram_blocks'])
    assert all(row['residual_strict_upper'] == '1/512' for row in fresh['Gram_blocks'])


def test_exact_scalar_gram_does_not_require_cholesky_or_full_rank():
    h = sparse.csr_matrix(np.zeros((2, 2), dtype=np.int64))
    q = np.zeros((2, 2), dtype=np.int64)
    assert g.verify_gram(h, q, F(1)) == F(27, 64)
    q[0, 0] = 1 << 19
    assert g.verify_gram(h, q, F(1)) == F(27, 64)


@pytest.mark.parametrize('base', [sparse.csr_matrix, sparse.coo_matrix])
def test_sparse_subclass_cannot_hide_original_negative_source_entries(base):
    class HiddenSource(base):
        def tocoo(self, copy=False):
            return sparse.coo_matrix((1, 1), dtype=np.int64)
    h = HiddenSource([[-16]], dtype=np.int64)
    q = g.propose_gram(sparse.csr_matrix([[0]], dtype=np.int64))
    assert int(h.toarray()[0, 0]) == -16
    assert -F(1)+F(27, 64) < 0  # The actual shifted source is not PSD.
    with pytest.raises(ValueError, match='standard.*sparse'):
        g.verify_gram(h, q)


def test_native_sparse_instance_cannot_shadow_its_conversion_method():
    h = sparse.csr_matrix([[-16]], dtype=np.int64)
    h.tocoo = lambda copy=False: sparse.coo_matrix((1, 1), dtype=np.int64)
    q = g.propose_gram(sparse.csr_matrix([[0]], dtype=np.int64))
    assert int(h.data[0]) == -16
    with pytest.raises(ValueError, match='native sparse storage'):
        g.verify_gram(h, q)


@pytest.mark.parametrize('form', ['coo', 'bsr'])
def test_sparse_stored_array_subclass_cannot_erase_original_entries(form):
    class HiddenEntries(np.ndarray):
        def copy(self, *args, **kwargs):
            return np.zeros(self.shape, dtype=np.int64)
    h = sparse.csr_matrix([[-16]], dtype=np.int64).asformat(form)
    h.data = h.data.view(HiddenEntries)
    assert int(np.asarray(h.data).flat[0]) == -16
    q = g.propose_gram(sparse.csr_matrix([[0]], dtype=np.int64))
    with pytest.raises(ValueError, match='native sparse storage'):
        g.verify_gram(h, q)


def test_dok_stored_dictionary_cannot_hide_original_entries():
    class HiddenEntries(dict):
        def items(self): return {}.items()
        def keys(self): return {}.keys()
        def values(self): return {}.values()
        def __len__(self): return 0
    h = sparse.dok_matrix([[-16]], dtype=np.int64)
    h._dict = HiddenEntries(h._dict)
    assert dict.__getitem__(h._dict, (0, 0)) == -16
    q = g.propose_gram(sparse.csr_matrix([[0]], dtype=np.int64))
    with pytest.raises(ValueError, match='native sparse storage'):
        g.verify_gram(h, q)


def test_dok_stored_integer_subclass_cannot_change_conversion_value():
    class HiddenInteger(int):
        def __int__(self): return 0
    h = sparse.dok_matrix([[-16]], dtype=np.int64)
    h._dict[(0, 0)] = HiddenInteger(-16)
    assert int.__int__(h._dict[(0, 0)]) == -16
    q = g.propose_gram(sparse.csr_matrix([[0]], dtype=np.int64))
    with pytest.raises(ValueError, match='native sparse storage'):
        g.verify_gram(h, q)


@pytest.mark.parametrize('numerator,denominator', [(1, -1), (1, 0), (True, 1), (1, True)])
def test_malformed_fraction_ceiling_cannot_bypass_exact_residual_admission(numerator, denominator):
    ceiling = F(1)
    ceiling._numerator, ceiling._denominator = numerator, denominator
    with pytest.raises(ValueError, match='positive exact rational'):
        g.verify_gram(sparse.csr_matrix([[0]], dtype=np.int64),
                      np.zeros((1, 1), dtype=np.int64), ceiling)


def test_positive_nonreduced_fraction_ceiling_is_detached_and_canonicalized():
    ceiling = F(1)
    ceiling._numerator, ceiling._denominator = 2, 4
    h, q = sparse.csr_matrix([[0]], dtype=np.int64), np.zeros((1, 1), dtype=np.int64)
    assert g.verify_gram(h, q, ceiling) == g.verify_gram(h, q, F(1, 2)) == F(27, 64)
    assert (ceiling.numerator, ceiling.denominator) == (2, 4)


@pytest.mark.parametrize('argument', ['numerator', 'factor'])
def test_masked_inputs_cannot_silently_become_complete_certificates(argument):
    h = np.zeros((1, 1), dtype=np.int64)
    q = np.zeros((1, 1), dtype=np.int64)
    if argument == 'numerator':
        h = np.ma.array(h, mask=True)
    else:
        q = np.ma.array(q, mask=True)
    with pytest.raises(ValueError, match='unmasked'):
        g.verify_gram(h, q, F(1))


@pytest.mark.parametrize('form', ['csr', 'csc', 'coo', 'bsr'])
def test_sparse_stored_masks_are_rejected_before_format_conversion(form):
    h = sparse.csr_matrix([[3]], dtype=np.int64).asformat(form)
    h.data = np.ma.array(h.data, mask=True)
    with pytest.raises(ValueError, match='unmasked'):
        g.verify_gram(h, np.zeros((1, 1), dtype=np.int64), F(1))
    assert np.ma.isMaskedArray(h.data) and np.all(h.data.mask)


@pytest.mark.parametrize('order', ['C', 'F'])
def test_corrupted_factor_is_rejected_by_integer_residual(order):
    h = sparse.csr_matrix(np.eye(2, dtype=np.int64))
    q = np.array(g.propose_gram(h), order=order)
    assert g.verify_gram(h, q) < g.CEILING
    q[0, 0] += 10000
    with pytest.raises(ArithmeticError, match='residual'):
        g.verify_gram(h, q)


@pytest.mark.parametrize('order', ['C', 'F'])
def test_integer_gram_layout_preserves_exact_scalar_residual_and_input(order):
    h = np.array([[0, 1, -2], [1, 3, 1], [-2, 1, -4]], dtype=np.int64)
    q = np.array([[1 << 19, 2, -3], [4, -(1 << 18), 6],
                  [7, 8, 1 << 17]], dtype=np.int64, order=order)
    before = q.copy()
    # A scalar rational oracle neither calls matmul nor reuses the verifier.
    expected = max(sum(abs(F(int(h[i, j]), 16)+(g.SHIFT if i == j else 0)-
                           sum(F(int(q[i, k])*int(q[j, k]), 1 << (2*g.BITS))
                               for k in range(3))) for j in range(3))
                   for i in range(3))
    assert g.verify_gram(sparse.csr_matrix(h), q, F(1)) == expected
    assert np.array_equal(q, before)
    assert q.flags.c_contiguous == (order == 'C')
    assert q.flags.f_contiguous == (order == 'F')


def test_duplicate_sparse_entries_are_canonicalized_before_exact_addition():
    h = sparse.csr_matrix((np.array([1, 2], dtype=np.int64), np.array([0, 0]), np.array([0, 2])), shape=(1, 1))
    canonical = sparse.csr_matrix([[3]], dtype=np.int64)
    q = g.propose_gram(canonical)
    assert g.verify_gram(h, q) == g.verify_gram(canonical, q)


def _duplicates(values, dtype, form):
    data = np.asarray(values, dtype=dtype)
    indices = np.zeros(len(data), dtype=np.int32)
    if form == 'csr':
        return sparse.csr_matrix((data, indices, np.array([0, len(data)], dtype=np.int32)), shape=(1, 1))
    return sparse.coo_matrix((data, (indices, indices)), shape=(1, 1))


@pytest.mark.parametrize('form', ['csr', 'coo'])
def test_duplicate_source_sum_cannot_wrap_to_a_false_ground_certificate(form):
    h = _duplicates([-(1 << 62)]*4, np.int64, form)
    q = g.propose_gram(sparse.csr_matrix([[0]], dtype=np.int64))
    # The mathematical numerator is -2**64; fixed-width coalescing made it0.
    with pytest.raises(OverflowError, match='source numerator sum'):
        g.verify_gram(h, q)
    assert h.nnz == 4 and all(int(value) == -(1 << 62) for value in h.data)


@pytest.mark.parametrize('form', ['csr', 'coo'])
@pytest.mark.parametrize('dtype,value', [(np.int8, -64), (np.uint8, 64)])
def test_narrow_integer_duplicates_are_summed_without_dtype_wrap(form, dtype, value):
    h = _duplicates([value]*4, dtype, form)
    q = np.zeros((1, 1), dtype=np.int64)
    expected = abs(F(4*value, 16)+g.SHIFT)
    assert g.verify_gram(h, q, F(20)) == expected
    assert h.dtype == np.dtype(dtype) and h.nnz == 4
    assert all(int(entry) == value for entry in h.data)


@pytest.mark.parametrize('form', ['csr', 'coo'])
def test_large_duplicate_cancellation_remains_exact_and_admissible(form):
    h = _duplicates([-(1 << 62)]*4+[(1 << 62)]*4, np.int64, form)
    assert g.verify_gram(h, np.zeros((1, 1), dtype=np.int64), F(1)) == g.SHIFT
    assert h.nnz == 8


@pytest.mark.parametrize('form', ['dense', 'csr', 'coo'])
def test_unsigned_source_values_cannot_silently_narrow(form):
    h = (np.array([[1 << 63]], dtype=np.uint64) if form == 'dense'
         else _duplicates([1 << 63], np.uint64, form))
    with pytest.raises(OverflowError, match='source numerator sum'):
        g.verify_gram(h, np.zeros((1, 1), dtype=np.int64))


@pytest.mark.parametrize('kind', ['matrix', 'array'])
@pytest.mark.parametrize('form', ['dense', 'csr', 'csc', 'coo', 'dia', 'bsr', 'lil', 'dok'])
def test_safe_unsigned_integer_formats_preserve_the_exact_numerator(form, kind):
    original = np.array([[3]], dtype=np.uint64)
    h = original if form == 'dense' else getattr(sparse, f'{form}_{kind}')(original)
    assert g.verify_gram(h, np.zeros((1, 1), dtype=np.int64), F(1)) == F(39, 64)


def test_source_range_is_checked_before_fixed_width_symmetry_subtraction():
    h = np.array([[0, 1 << 62], [-(1 << 62), 0]], dtype=np.int64)
    with pytest.raises(OverflowError, match='accumulation or residual'):
        g.verify_gram(h, np.zeros((2, 2), dtype=np.int64))


@pytest.mark.parametrize('order', ['C', 'F'])
def test_integer_overflow_and_nonhermitian_inputs_are_rejected(order):
    with pytest.raises(OverflowError):
        g.verify_gram(sparse.eye(2, dtype=np.int64),
                      np.full((2, 2), 1 << 40, dtype=np.int64, order=order))
    with pytest.raises(ValueError, match='symmetric'):
        g.verify_gram(sparse.csr_matrix([[0, 1], [0, 0]], dtype=np.int64), np.zeros((2, 2), dtype=np.int64))
    with pytest.raises(ValueError):
        g.verify_gram(sparse.eye(2), np.zeros((2, 2), dtype=np.int64))


@pytest.mark.parametrize('value', [-1, 17, True, 4.0])
def test_particle_number_domain(value):
    with pytest.raises(ValueError): g.balanced_indices(value)
