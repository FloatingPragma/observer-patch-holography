from fractions import Fraction as F

import numpy as np
import pytest
from scipy import sparse

from . import gap


def small_witness():
    # A separate complete two-dimensional integer Hamiltonian has exact
    # centered energies -33/8 and -4. It tests the certificate arithmetic.
    h = sparse.csr_matrix(np.diag([-2, 0]).astype(np.int64))
    p, q = gap.propose_witness(h)
    return h, p, q


def scalar_oracle(h, p, q):
    dense = h.toarray()
    n = len(p)
    norm = sum(int(x)**2 for x in p)
    energy = sum(F(int(p[i])*int(dense[i, j])*int(p[j]), 16*norm)
                 for i in range(n) for j in range(n))-4
    residual = max(sum(abs(
        F(int(dense[i, j]), 16)
        -(4+gap.PROPOSAL_CUTOFF if i == j else 0)
        +F(int(p[i])*int(p[j]), 1 << (2*gap.TRIAL_BITS))
        -sum(F(int(q[i, k])*int(q[j, k]), 1 << (2*gap.FACTOR_BITS))
             for k in range(n))) for j in range(n)) for i in range(n))
    return energy, residual


@pytest.mark.parametrize('order', ['C', 'F'])
def test_exact_rank_one_arithmetic_matches_independent_scalar_oracle(order):
    h, p, q = small_witness()
    q = np.array(q, order=order)
    before_p, before_q = p.copy(), q.copy()
    energy, residual = scalar_oracle(h, p, q)
    value = gap.verify_witness(h, p, q)
    assert energy == -F(33, 8) == value['trial_rayleigh']
    assert residual == value['gram_residual'] < F(1, 512)
    assert value['neutral_gap_strict_lower'] == F(1, 128)
    np.testing.assert_array_equal(p, before_p)
    np.testing.assert_array_equal(q, before_q)
    assert q.flags.c_contiguous == (order == 'C')
    assert q.flags.f_contiguous == (order == 'F')


def test_offdiagonal_source_terms_are_checked_by_independent_rational_oracle():
    h = sparse.csr_matrix(np.array([[-2, 1, 0], [1, 2, 1], [0, 1, 4]], dtype=np.int64))
    p, q = gap.propose_witness(h)
    energy, residual = scalar_oracle(h, p, q)
    value = gap.verify_witness(h, p, q)
    assert value['trial_rayleigh'] == energy
    assert value['gram_residual'] == residual
    changed = h.copy()
    changed[0, 1] = changed[1, 0] = -1
    with pytest.raises(ArithmeticError, match='Rayleigh|Gram residual'):
        gap.verify_witness(changed, p, q)


def test_wrong_relative_trial_sign_is_not_hidden_by_normalization():
    h = sparse.csr_matrix(np.array([[-2, 1], [1, 0]], dtype=np.int64))
    p, q = gap.propose_witness(h)
    gap.verify_witness(h, p, q)
    damaged = p.copy()
    damaged[1] *= -1
    assert sum(int(x)**2 for x in p) == sum(int(x)**2 for x in damaged)
    with pytest.raises(ArithmeticError, match='Rayleigh'):
        gap.verify_witness(h, damaged, q)


def test_wrong_gram_factor_cannot_be_accepted_from_a_good_trial():
    h, p, q = small_witness()
    q[0, 0] += 100000
    with pytest.raises(ArithmeticError, match='Gram residual'):
        gap.verify_witness(h, p, q)


def test_a_second_low_energy_direction_defeats_rank_one_certificate():
    h, p, q = small_witness()
    doubled = sparse.csr_matrix(np.diag([-2, -2]).astype(np.int64))
    with pytest.raises(ArithmeticError, match='Gram residual'):
        gap.verify_witness(doubled, p, q)


def test_sparse_subclass_cannot_hide_a_second_low_energy_direction():
    h, p, q = small_witness()
    class HiddenSecondGround(sparse.csr_matrix):
        def tocoo(self, copy=False):
            return h.tocoo(copy=copy)
    doubled = HiddenSecondGround(np.diag([-2, -2]).astype(np.int64))
    assert np.array_equal(doubled.toarray(), np.diag([-2, -2]))
    with pytest.raises(ValueError, match='standard.*sparse'):
        gap.verify_witness(doubled, p, q)


def test_returned_fraction_bound_cannot_relax_a_later_ground_admission():
    h, p, q = small_witness()
    value = gap.verify_witness(h, p, q)
    zero = sparse.csr_matrix(np.zeros((2, 2), dtype=np.int64))
    p0, q0 = gap.propose_witness(zero)
    original = (gap.GROUND_UPPER.numerator, gap.GROUND_UPPER.denominator)
    value['ground_energy_strict_upper']._numerator = 0
    value['ground_energy_strict_upper']._denominator = 1
    try:
        # This source has two exactly degenerate centered eigenvalues -4.
        with pytest.raises(ArithmeticError, match='Rayleigh'):
            gap.verify_witness(zero, p0, q0)
    finally:
        # Also leave a deliberately unfixed implementation clean after failure.
        gap.GROUND_UPPER._numerator, gap.GROUND_UPPER._denominator = original


def test_every_returned_fraction_is_detached_from_later_checks_and_proof_constants():
    h, p, q = small_witness()
    value = gap.verify_witness(h, p, q)
    expected = {name: F(bound.numerator, bound.denominator) for name, bound in value.items()}
    names = ('GROUND_UPPER', 'EXCITED_LOWER', 'GAP_LOWER')
    original = {name: (getattr(gap, name).numerator, getattr(gap, name).denominator) for name in names}
    try:
        for bound in value.values():
            bound._numerator = 17
            bound._denominator = 19
        assert gap.verify_witness(h, p, q) == expected
        assert all(getattr(gap, name) == F(*parts) for name, parts in original.items())
    finally:
        for name, parts in original.items():
            getattr(gap, name)._numerator, getattr(gap, name)._denominator = parts


def test_exact_Temple_and_fidelity_bounds_cover_a_mixed_trial():
    h, _, _ = small_witness()
    p = np.array([1 << gap.TRIAL_BITS, 1 << (gap.TRIAL_BITS-3)], dtype=np.int64)
    scale = 1 << gap.TRIAL_BITS
    shifted = h.toarray()/16+np.outer(p/scale, p/scale)
    shifted[np.diag_indices(2)] -= float(4+gap.PROPOSAL_CUTOFF)
    q = np.rint(np.linalg.cholesky(shifted)*(1 << gap.FACTOR_BITS)).astype(np.int64)
    value = gap.verify_witness(h, p, q)
    weights = [F(int(x)**2, sum(int(y)**2 for y in p)) for x in p]
    energies = [-F(33, 8), -F(4)]
    rho = sum(w*e for w, e in zip(weights, energies))
    variance = sum(w*(e-rho)**2 for w, e in zip(weights, energies))
    assert value['trial_rayleigh'] == rho
    assert value['trial_residual_variance'] == variance > 0
    assert value['Temple_ground_energy_lower'] <= energies[0]
    assert value['trial_ground_infidelity_upper'] >= weights[1]
    assert value['trial_ground_infidelity_upper'] > gap.INFIDELITY_CEILING
    assert value['Temple_ground_energy_lower'] < rho


def test_eigensolver_values_are_not_admission_premises(monkeypatch):
    h = sparse.csr_matrix(np.diag([-2, 0, 1]).astype(np.int64))
    monkeypatch.setattr(gap, 'eigsh', lambda *args, **kwargs:
                        (np.array([1e100]), np.array([[1.], [0.], [0.]])))
    p, q = gap.propose_witness(h)
    assert gap.verify_witness(h, p, q)['trial_rayleigh'] == -F(33, 8)


@pytest.mark.parametrize('argument', ['trial', 'factor'])
@pytest.mark.parametrize('dtype', [np.float64, np.int16, np.uint64, np.bool_])
def test_no_silent_integer_narrowing_or_coercion(argument, dtype):
    h, p, q = small_witness()
    if argument == 'trial':
        p = p.astype(dtype)
    else:
        q = q.astype(dtype)
    with pytest.raises(ValueError, match='int64'):
        gap.verify_witness(h, p, q)


@pytest.mark.parametrize('argument', ['source', 'stored_source', 'trial', 'factor'])
def test_masks_cannot_be_erased_before_certificate_validation(argument):
    h, p, q = small_witness()
    if argument == 'source':
        h = np.ma.array(h.toarray(), mask=True)
    elif argument == 'stored_source':
        h.data = np.ma.array(h.data, mask=True)
    elif argument == 'trial':
        p = np.ma.array(p, mask=True)
    else:
        q = np.ma.array(q, mask=True)
    with pytest.raises(ValueError, match='unmasked'):
        gap.verify_witness(h, p, q)


@pytest.mark.parametrize('argument', ['trial', 'factor'])
def test_wrong_shapes_and_zero_trial_are_refused(argument):
    h, p, q = small_witness()
    if argument == 'trial':
        p = p[:, None]
    else:
        q = q[:, :1]
    with pytest.raises(ValueError, match='shape'):
        gap.verify_witness(h, p, q)
    h, p, q = small_witness()
    with pytest.raises(ValueError, match='nonzero'):
        gap.verify_witness(h, np.zeros_like(p), q)


@pytest.mark.parametrize('order', ['C', 'F'])
def test_products_are_bounded_before_integer_gram_arithmetic(order):
    h, p, _ = small_witness()
    q = np.full((2, 2), 1 << 40, dtype=np.int64, order=order)
    with pytest.raises(OverflowError, match='accumulation'):
        gap.verify_witness(h, p, q)


def test_source_range_is_checked_before_symmetry_subtraction():
    h, p, q = small_witness()
    h = np.array([[0, 1 << 62], [-(1 << 62), 0]], dtype=np.int64)
    with pytest.raises(OverflowError, match='accumulation'):
        gap.verify_witness(h, p, q)


@pytest.mark.parametrize('form', ['csr', 'coo'])
def test_original_sparse_duplicates_cannot_wrap_to_a_forged_source(form):
    _, p, q = small_witness()
    data = np.array([-(1 << 62)]*4, dtype=np.int64)
    indices = np.zeros(4, dtype=np.int32)
    h = (sparse.csr_matrix((data, indices, [0, 4, 4]), shape=(2, 2)) if form == 'csr'
         else sparse.coo_matrix((data, (indices, indices)), shape=(2, 2)))
    with pytest.raises(OverflowError, match='source numerator sum'):
        gap.verify_witness(h, p, q)


def test_safe_duplicate_cancellation_has_the_same_exact_certificate():
    h, p, q = small_witness()
    values = np.array([-(1 << 62)]*4+[(1 << 62)]*4+[-1, -1], dtype=np.int64)
    indices = np.zeros(len(values), dtype=np.int32)
    duplicate = sparse.coo_matrix((values, (indices, indices)), shape=(2, 2))
    assert gap.verify_witness(duplicate, p, q) == gap.verify_witness(h, p, q)
    assert duplicate.nnz == 10


def spin_pair_trials():
    indices = gap.global_ground.balanced_indices(8)
    positions = {int(mask): j for j, mask in enumerate(indices)}
    # Three doublons plus two singly occupied sites: N8 and Sz0 exactly.
    doublons = sum(3 << (2*j) for j in (2, 3, 4))
    up_down, down_up = doublons | 1 | (2 << 2), doublons | 2 | (1 << 2)
    triplet = np.zeros(4900, dtype=np.int64)
    triplet[positions[up_down]] = triplet[positions[down_up]] = 1
    singlet = triplet.copy()
    singlet[positions[down_up]] = -1
    return indices, singlet, triplet


def test_complete_spin_raising_matches_independent_ordered_CAR_oracle():
    indices = gap.global_ground.balanced_indices(8)
    raising, targets = gap.spin_raising_matrix(indices)
    target_rows = {int(mask): row for row, mask in enumerate(targets)}
    expected = {}
    for column, state in enumerate(indices):
        for site in range(8):
            up, down = 2*site, 2*site+1
            mask = int(state)
            if not (mask >> down & 1) or mask >> up & 1:
                continue
            sign = (-1)**(mask & ((1 << down)-1)).bit_count()
            mask ^= 1 << down
            sign *= (-1)**(mask & ((1 << up)-1)).bit_count()
            mask |= 1 << up
            expected[target_rows[mask], column] = sign
    actual = raising.tocoo()
    assert raising.shape == (3136, 4900) and raising.nnz == 9800
    assert {(int(i), int(j)): int(x) for i, j, x in
            zip(actual.row, actual.col, actual.data, strict=True)} == expected
    assert all(x == 1 for x in expected.values())
    assert all(int(mask).bit_count() == 8 and (int(mask) & 0x5555).bit_count() == 5
               for mask in targets)


def test_exact_spin_expectation_separates_singlet_and_wrong_relative_sign_triplet():
    indices, singlet, triplet = spin_pair_trials()
    assert gap.trial_spin_expectation(singlet, indices) == 0
    assert gap.trial_spin_expectation(triplet, indices) == 2
    assert sum(int(x)**2 for x in singlet) == sum(int(x)**2 for x in triplet) == 2
    # Cross terms vanish between exact singlet and triplet eigenvectors.
    mixture = 3*singlet+triplet
    assert gap.trial_spin_expectation(mixture, indices) == F(1, 5)
    assert gap.trial_spin_expectation(mixture*12345, indices) == F(1, 5)
    assert gap.SPIN_EXPECTATION_CEILING < 2*(1-gap.INFIDELITY_CEILING)


@pytest.mark.parametrize('mutation', ['reorder', 'duplicate', 'wrong_weight'])
def test_ladder_rejects_changed_source_address_basis(mutation):
    indices, singlet, _ = spin_pair_trials()
    if mutation == 'reorder':
        indices[[0, 1]] = indices[[1, 0]]
    elif mutation == 'duplicate':
        indices[1] = indices[0]
    else:
        indices[0] = 0
    with pytest.raises(ValueError, match='canonical balanced source address order'):
        gap.trial_spin_expectation(singlet, indices)


@pytest.mark.parametrize('argument', ['trial', 'addresses'])
@pytest.mark.parametrize('dtype', [np.float64, np.int8, np.uint64, np.bool_])
def test_spin_ladder_does_not_narrow_or_coerce_inputs(argument, dtype):
    indices, singlet, _ = spin_pair_trials()
    if argument == 'trial':
        singlet = singlet.astype(dtype)
    else:
        indices = indices.astype(dtype)
    with pytest.raises(ValueError, match='int64'):
        gap.trial_spin_expectation(singlet, indices)


@pytest.mark.parametrize('argument', ['trial', 'addresses'])
def test_spin_ladder_rejects_masks_before_numpy_coercion(argument):
    indices, singlet, _ = spin_pair_trials()
    if argument == 'trial':
        singlet = np.ma.array(singlet, mask=True)
    else:
        indices = np.ma.array(indices, mask=True)
    with pytest.raises(ValueError, match='unmasked'):
        gap.trial_spin_expectation(singlet, indices)


@pytest.mark.parametrize('coefficient', [-(1 << 63), 1 << 62])
def test_spin_ladder_checks_accumulation_before_fixed_width_multiplication(coefficient):
    indices, singlet, _ = spin_pair_trials()
    singlet[np.flatnonzero(singlet)[0]] = coefficient
    with pytest.raises(OverflowError, match='spin-raising accumulation'):
        gap.trial_spin_expectation(singlet, indices)


def test_spin_ladder_requires_the_complete_nonzero_trial_and_returns_fresh_data():
    indices, singlet, _ = spin_pair_trials()
    with pytest.raises(ValueError, match='shape'):
        gap.trial_spin_expectation(singlet[:, None], indices)
    with pytest.raises(ValueError, match='nonzero'):
        gap.trial_spin_expectation(np.zeros_like(singlet), indices)
    raising, targets = gap.spin_raising_matrix(indices)
    raising.data[:] = 0
    targets[:] = 0
    assert gap.trial_spin_expectation(singlet, indices) == 0
    fresh, addresses = gap.spin_raising_matrix(indices)
    assert fresh.nnz == 9800 and np.all(fresh.data == 1) and np.all(addresses != 0)


@pytest.mark.parametrize('mutation', ['zero', 'relative_sign'])
def test_a_false_small_spin_norm_cannot_replace_the_actual_source_ladder(monkeypatch, mutation):
    indices = gap.global_ground.balanced_indices(8)
    h = gap.global_ground.source_block(8)
    original = gap.spin_raising_matrix
    def corrupt(addresses):
        raising, targets = original(addresses)
        if mutation == 'zero':
            raising.data[:] = 0
        else:
            raising.data[0] *= -1
        return raising, targets
    monkeypatch.setattr(gap, 'spin_raising_matrix', corrupt)
    with pytest.raises(ArithmeticError, match='different ordered addresses'):
        gap._verify_spin_source(h, indices)


def test_source_spin_binding_rejects_a_different_balanced_hamiltonian():
    indices = gap.global_ground.balanced_indices(8)
    h = gap.global_ground.source_block(8).copy()
    h.data[0] += 1
    with pytest.raises(ArithmeticError, match='different ordered addresses'):
        gap._verify_spin_source(h, indices)


def test_full_centered_source_has_a_certified_gap_in_every_number_sector(monkeypatch):
    def forbidden_import():
        raise AssertionError('Lieb theorem cannot be a new gap/uniqueness admission premise')
    monkeypatch.setattr(gap.cube, 'lieb_premises', forbidden_import)
    gap._certificate_items.cache_clear()
    value = gap.certificate()
    assert value['full_Fock_dimension'] == 65536
    assert value['balanced_N8_dimension'] == 4900
    assert value['ground_energy_strict_upper'] == '-1055/256'
    assert value['N8_excited_energy_strict_lower'] == '-1053/256'
    assert value['neutral_gap_strict_lower'] == '1/128'
    assert value['other_charge_sector_gap_strict_lower'] == '1/8'
    assert value['full_centered_gap_strict_lower'] == '1/128'
    assert value['ground_energy_enclosure_strict_lower'] == '-206203/50000'
    assert value['ground_energy_enclosure_strict_upper'] == '-82481/20000'
    assert value['normalized_trial_ground_infidelity_strict_upper'] == '1/1000000'
    assert value['normalized_trial_ground_half_trace_distance_strict_upper'] == '1/1000'
    assert value['bounded_observable_expectation_error_strict_factor'] == '1/500'
    assert value['normalized_trial_spin_squared_expectation_strict_upper'] == '1'
    assert value['spin_raising_target_dimension'] == 3136
    assert value['full_source_SU2_and_spin_frame_checked'] is True
    assert value['ground_hidden_spin'] == 0 and value['full_ground_multiplicity'] == 1
    assert value['imported_Lieb_theorem_used_for_gap_or_uniqueness'] is False
    assert value['native_ground_preparation_executed'] is False
    assert value['numerical_eigenvalue_or_residual_used_for_admission'] is False
    assert value['thermodynamic_gap_proved'] is False
    assert value['relativistic_mass_identified'] is False
    assert value['localized_topological_protection_proved'] is False
    value['full_centered_gap_strict_lower'] = '1'
    value['native_ground_preparation_executed'] = True
    again = gap.certificate()
    assert again['full_centered_gap_strict_lower'] == '1/128'
    assert again['native_ground_preparation_executed'] is False
