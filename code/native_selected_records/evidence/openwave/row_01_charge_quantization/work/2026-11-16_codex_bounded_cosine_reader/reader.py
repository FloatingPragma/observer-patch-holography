"""Bounded first-moment reader; standard library, no file/process/network access."""
import math

def bounded_frequency(record, eta, probability_error=0.0):
    """eta is an explicit external bound on TOTAL unwanted-band population.

    The target band may contain both frequency signs. The result is unsigned
    on [0,pi]; neither signed population nor eta is inferred from this record.
    probability_error is a deterministic bound, not an estimated confidence.
    """
    if not math.isfinite(eta) or not 0 <= eta < 1:
        raise ValueError('eta must be a finite external certificate in [0,1)')
    if not math.isfinite(probability_error) or not 0 <= probability_error <= .5:
        raise ValueError('invalid deterministic probability-error bound')
    for axis in ('x', 'y'):
        pair = [record[axis+'_plus'], record[axis+'_minus']]
        if any(not math.isfinite(x) or not 0 <= x <= 1 for x in pair):
            raise ValueError('invalid binary outcome probability')
        if abs(sum(pair)-1) > 1e-12:
            raise ValueError('binary outcomes do not sum to one')
    real = record['x_plus']-record['x_minus']
    imaginary = record['y_plus']-record['y_minus']
    lower_cosine = max(-1., min(1., (real-2*probability_error-eta)/(1-eta)))
    upper_cosine = max(-1., min(1., (real+2*probability_error+eta)/(1-eta)))
    return {'frequency_interval': [math.acos(upper_cosine), math.acos(lower_cosine)],
            'chart': [0., math.pi], 'real_moment': real,
            'imaginary_moment_separate_sign_record': imaginary,
            'eta_external_certificate': eta,
            'deterministic_probability_error': probability_error,
            'finite_shot_confidence_claim': False}
