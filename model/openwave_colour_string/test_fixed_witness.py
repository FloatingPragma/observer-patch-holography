"""A complete proof frame is checked, never trusted as a saved spectrum."""
from copy import deepcopy
from fractions import Fraction as F
from unittest.mock import patch
import numpy as np
import pytest
import sympy as sy
from .fixed_witness import spectral_enclosure


def witness():
    D = 2**40
    return {'dyadic_denominator': D,
            'dyadic_eigenframe': [[str(D), '0'], ['0', str(D)]],
            'dyadic_eigenvalues': [str(D), str(2*D)]}


def test_native_source_and_full_residual_are_recomputed_without_eigensolver():
    with patch('model.openwave_six_quark.certificate.eigh', side_effect=AssertionError('fresh frame')):
        exact = sy.diag(1, 2)
        result = spectral_enclosure(exact, witness())
        changed = spectral_enclosure(exact+sy.Matrix([[0, sy.Rational(1, 100)], [sy.Rational(1, 100), 0]]), witness())
    assert F(changed['all_eigenvalue_error_upper']) > F(result['all_eigenvalue_error_upper'])
    values = np.linalg.eigvalsh(np.array([[1, .01], [.01, 2]]))
    for value, (lo, hi) in zip(values, changed['all_eigenvalue_intervals']):
        assert float(F(lo)) <= value <= float(F(hi))


@pytest.mark.parametrize('case', ['missing_column', 'zero_column', 'boolean', 'overflow', 'unsorted', 'denominator'])
def test_malformed_or_incomplete_frame_is_rejected(case):
    data = deepcopy(witness())
    if case == 'missing_column':
        data['dyadic_eigenframe'][0].pop()
    elif case == 'zero_column':
        data['dyadic_eigenframe'][1][1] = '0'
    elif case == 'boolean':
        data['dyadic_eigenframe'][0][0] = True
    elif case == 'overflow':
        data['dyadic_eigenvalues'][1] = str(2**80)
    elif case == 'unsorted':
        data['dyadic_eigenvalues'].reverse()
    else:
        data['dyadic_denominator'] = True
    with pytest.raises(AssertionError):
        spectral_enclosure(sy.diag(1, 2), data)
