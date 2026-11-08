"""Exact finite Gram admission and original-probability error bounds."""
from fractions import Fraction as F


def _nonnegative(name, value):
    if type(value) not in (int, F) or value < 0:
        raise ValueError(name+' must be an exact nonnegative number')
    return F(value)


def probability_error(gamma, gram_error, overlap_norm, overlap_error, instrument_error=F(), effect_error=F()):
    """Bound one original event weight; never divide by its success probability.

    gamma certifies the reference Gram lower bound. The perturbation must be
    strictly smaller. The source compiler's complete error is supplied and
    kept separate from Gram and overlap errors.
    """
    gamma = _nonnegative('Gram gap', gamma)
    dc = _nonnegative('Gram error', gram_error)
    norm = _nonnegative('overlap norm', overlap_norm)
    dz = _nonnegative('overlap error', overlap_error)
    detector = _nonnegative('instrument error', instrument_error)
    da = _nonnegative('coordinate effect error', effect_error)
    if gamma == 0 or dc >= gamma:
        raise ValueError('Gram invertibility requires error strictly below its positive gap')
    if detector > 1:
        raise ValueError('complete half-trace error must be at most one')
    alpha = gamma-dc
    raw = (2*norm+dz)*dz/alpha+norm*norm*dc/(alpha*alpha)+norm*norm*da/alpha
    return {'certified_common_gram_gap': alpha,
            'original_probability_error': raw,
            'coordinate_effect_approximation_error': da,
            'complete_implementation_error': detector,
            'normalized_input_event_error_upper_bound': min(F(1), raw+detector),
            'division_by_success_probability': False}


def bounded_probe_errors(modes, amplitude, raw_gram_error, raw_overlap_error):
    """Convert complete bounded-score errors using known LCU coefficients.

    A complex test is the sum of two sine probes, so ||A||<=2/|s|.
    Inputs bound complex normalized kernel errors, including preparation,
    control and sampling contributions. No observed success is used.
    """
    if type(modes) is not int or modes<1:
        raise ValueError('positive finite mode count required')
    if type(amplitude) not in (int,F) or amplitude==0:
        raise ValueError('exact nonzero sine amplitude required')
    eg=_nonnegative('raw Gram error',raw_gram_error)
    ez=_nonnegative('raw overlap error',raw_overlap_error)
    normalization=F(2)/abs(amplitude)
    return {'known_lcu_norm_upper_bound':normalization,
            'gram_operator_error':modes*normalization**2*eg,
            'overlap_vector_error':modes*normalization*ez,
            'division_by_observed_acceptance':False}
