from fractions import Fraction as F
from itertools import product

import pytest

from .symmetry import (abundance_kernel_certificate, bias_bound, conjugate, local_certificate, local_step,
                       paired_history, selection_bound)


def test_every_complete_local_branch_conjugates():
    receipt = local_certificate()
    assert receipt['complete_local_branches'] == 2592
    assert 0 < receipt['accepted_nonzero_twists'] < 2592


def test_the_coin_transformation_is_necessary():
    state = (5, 0, 0, 0)
    a = local_step(state, 0)
    b = local_step(conjugate(state), 0)
    assert b[:4] != conjugate(a[:4])
    assert local_step(conjugate(state), 1)[:4] == conjugate(a[:4])


def test_rejections_and_no_operations_have_conjugate_records():
    rejected = local_step((2, 0, 0, 0), 0)
    no_op = local_step((3, 3, 0, 0), 1)
    assert rejected[-2:] == (1, False)
    assert no_op[-2:] == (0, False)


def test_native_history_and_single_sample_fluctuation():
    result = paired_history(sweeps=12)
    assert result['paired_mean_zero']
    assert result['total_integer_lift'] == 0
    assert result['maximum_single_history_odd_count'] > 0
    assert len(set(result['complete_record_sha256'])) == 2


def test_actual_population_kernel_has_the_formal_symmetry_premises():
    receipt = abundance_kernel_certificate()
    assert receipt['states'] == 1296
    assert receipt['complete_record_branches'] == 7776
    assert receipt['equivariant_kernel'] and receipt['even_preparation']


def test_bias_bound_against_all_two_state_grid_kernels():
    # Independent finite distributions and kernels, not the native checker.
    for a, b, p in product((F(i, 4) for i in range(5)), repeat=3):
        kernel = ((a, 1-a), (b, 1-b))
        defect = abs(a+b-1)
        law = (p, 1-p)
        initial = abs(2*p-1)
        for n in range(1, 5):
            law = tuple(sum(law[i]*kernel[i][j] for i in (0, 1)) for j in (0, 1))
            assert abs(law[0]-law[1]) <= bias_bound(initial, [defect]*n)


def test_symmetric_filter_and_asymmetric_postselection_control():
    assert selection_bound(F(1, 100), F(1, 4)) == F(1, 25)
    assert selection_bound(0, F(1, 100)) == 0
    assert selection_bound(F(1, 2), F(1, 4)) == 1
    symmetric_mean = F(1, 2)*1 + F(1, 2)*(-1)
    positive_only_mean = (F(1, 2)*1) / F(1, 2)
    assert symmetric_mean == 0 and positive_only_mean == 1


def test_bias_bound_rejects_nonprobabilities_and_inexact_inputs():
    for initial, defects in ((-1, []), (0.0, []), (0, [F(2)]), (0, [-1])):
        with pytest.raises(ValueError):
            bias_bound(initial, defects)
    with pytest.raises(ValueError):
        selection_bound(0, 0)
    assert bias_bound(F(1, 10), [0]*100) == F(1, 10)
    assert bias_bound(0, [F(1, 500)]*10) == 1 - F(499, 500)**10


def test_bias_bound_consumes_a_generator_once_without_losing_its_defects():
    assert bias_bound(0, (F(1, 2) for _ in range(2))) == F(3, 4)


@pytest.mark.parametrize('coin', [False, True, 0.0, 1.0])
def test_exact_native_coin_rejects_boolean_and_float_aliases(coin):
    with pytest.raises(ValueError):
        local_step((2, 0, 0, 0), coin)
