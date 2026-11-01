from fractions import Fraction as F

import numpy as np
import pytest

from model.retarded_matter_response import static as s


def test_exact_residual_admits_static_response_and_memory():
    result=s.certificate()
    assert result['C_operator_error_strict']=='1/10000000'
    assert result['static_memory_operator_error_strict']=='1/1000000'
    assert result['floating_point_solve_admits_the_result'] is False
    assert result['exact_equality_to_affine_adjacency_asserted'] is False


def test_corrupt_or_zero_proposal_is_not_admitted():
    with pytest.raises(ArithmeticError):s.verify(np.zeros((560,8),dtype=np.int64))
    value=s.propose();value[0,0]+=s.SCALE//100
    with pytest.raises(ArithmeticError):s.verify(value)
    with pytest.raises(ValueError):s.verify(np.zeros((560,8)))


def test_integer_square_root_bound_preserves_cancellation():
    value=s.quadratic_frobenius_bound(np.array([[-1732050808]],dtype=np.int64),
                                     np.array([[1000000000]],dtype=np.int64),10**9)
    assert F(0)<value<F(1,10**9)
    assert s.SQRT_LOWER**2<3*s.SQRT_SCALE**2<(s.SQRT_LOWER+1)**2
