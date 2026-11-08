"""Exact constants in the source Gram lower bound and its finite compression."""
from fractions import Fraction as F


def gram_bounds(m_lower=F(1, 1000)):
    if type(m_lower) is not F or not 0 < m_lower <= F(1, 2):
        raise ValueError('a positive rational spin order lower bound is required')
    return {
        'weighted_field_gram_lower': m_lower*m_lower,
        'weighted_field_gram_upper': F(3),
        'gram_condition_number_upper': F(3)/(m_lower*m_lower),
        'inverse_square_root_norm_upper': 1/m_lower,
        'cyclic_translation_spectral_multiplicity': 1,
        'translation_spectral_support': 'all R3 up to Lebesgue null sets',
        'source_gram_whitening_implemented_as_native_measurement': False,
    }


def frame_perturbation(m_lower, gram_error):
    """An estimated finite Gram retains invertibility only below its true gap.

    Tests must be orthonormal in the weighted Fourier norm before this applies.
    An arbitrary probe family need not have this lower eigenvalue bound.
    """
    lower = gram_bounds(m_lower)['weighted_field_gram_lower']
    if type(gram_error) is not F or not 0 <= gram_error < lower:
        raise ValueError('a rational Gram error strictly below m_lower squared is required')
    return {'estimated_gram_lower': lower-gram_error,
            'estimated_gram_upper': F(3)+gram_error,
            'estimated_gram_condition_number_upper': (F(3)+gram_error)/(lower-gram_error)}
