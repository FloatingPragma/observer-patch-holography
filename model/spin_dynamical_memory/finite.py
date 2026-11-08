"""Rational full-source memory and bounded-history certificates."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
from pathlib import Path
import numpy as np
from model.native_ordered_ground import certificate, source
from model.centered_particle_spectrum.certificate import decimal_bounds, sqrt_bounds
from .algebra import memory_interval

ROOT = Path(__file__).resolve().parents[2]
WITNESS = 'model/native_ordered_ground/witness.npz'


@lru_cache(maxsize=1)
def verified_source():
    with np.load(ROOT/WITNESS, allow_pickle=False) as archive:
        witness = {key: archive[key] for key in archive.files}
    proof = certificate.verify(witness)
    h, scale = source.matrix()
    ground = witness['vectors'][:, 0].astype(object)
    flip = np.arange(256)^255
    if not np.array_equal(h[np.ix_(flip, flip)], h):
        raise ArithmeticError('exact transverse parity symmetry changed')
    return witness, proof, h, scale, ground


def _interval(center, error, digits=12, nonnegative=False):
    low, high = center-error, center+error
    if nonnegative:
        low = max(F(), low)
    return decimal_bounds(low, high, digits)


def memory_certificate():
    witness, proof, h, scale, ground = verified_source()
    norm = sum(int(x)**2 for x in ground)
    z = np.array([1-2*(mask&1) for mask in range(256)], dtype=object)
    a = z*ground
    ha = h.astype(object)@a
    hg = h.astype(object)@ground
    hg_mean = F(int(ground@hg), scale*norm)
    ha_mean = F(int(a@ha), scale*norm)
    ha_second = F(int(ha@ha), scale*scale*norm)
    state_error = F(proof['normalized_ground_vector_error'])
    # The exact source norm is at most J|E|+BN=14. Comparing the two
    # energy expectations gives4*14*state_error. Variance of H on AΩ
    # is invariant under the reference energy subtraction.
    mean = ha_mean-hg_mean
    mean_error = 56*state_error
    variance = ha_second-ha_mean**2
    variance_error = 6*14**2*state_error
    mlo, mhi = mean-mean_error, mean+mean_error
    vlo, vhi = max(F(), variance-variance_error), variance+variance_error
    if vlo <= F(43, 100):
        raise ArithmeticError('source complement memory strength not certified')
    unit = 1 << certificate.BITS
    vectors = witness['vectors'].astype(object)
    energies = [F(int(x), unit) for x in witness['energies']]
    coefficients = vectors.T@a
    weights = [F(int(c)**2, unit**2*norm) for c in coefficients]
    e0lo, e0hi = map(F, proof['ground_energy_interval'])
    origin = energies[0]
    delta = F(proof['complete_polar_Weyl_radius'])+max(origin-e0lo, e0hi-origin)
    gram = F(proof['Gram_defect'])
    shape = gram*(2+gram)+2*state_error
    rows = []
    for eta in (F(1, 4), F(1), F(4)):
        response = sum(w/(e-origin+eta) for w, e in zip(weights, energies))
        error = delta/eta**2+shape/eta
        gl, gh = response-error, response+error
        sl, sh = memory_interval((mlo, mhi), (gl, gh), eta)
        if sl <= 0:
            raise ArithmeticError('source memory not strictly positive')
        rows.append({'eta': str(eta), 'normalized_resolvent_interval': decimal_bounds(gl, gh, 12),
            'complement_memory_interval': decimal_bounds(sl, sh, 12),
            'complete_bounded_score_interval': decimal_bounds(eta*gl, eta*gh, 12),
            'retained_complement_score_interval': decimal_bounds(1-eta*gh, 1-eta*gl, 12),
            'all_spectral_weights_count': len(weights)})
    return {'source_dimension': 256, 'bonds': 12, 'J': '1', 'selection_field': '1/4',
        'probe': 'local encoded Z0, a=Z0 Omega, norm(a)=1',
        'exact_vanishing_ground_mean_from_parity': True,
        'mean_excitation_energy_interval': decimal_bounds(mlo, mhi, 12),
        'complement_memory_strength_interval': decimal_bounds(vlo, vhi, 12),
        'rows': rows, 'source_ground_certificate_replayed': True,
        'finite_field_ground_substituted_for_infinite_ordered_ground': False}


def _phase_polynomial_apply(h, vector_real, vector_imag, denominator, sign, order):
    """Integer Taylor action for exp(sign*i*h/denominator)."""
    if sign not in (-1, 1) or type(order) is not int or order < 0:
        raise ValueError('signed nonnegative Taylor order required')
    current_r, current_i = vector_real.copy(), vector_imag.copy()
    out_r = np.zeros(len(vector_real), dtype=object)
    out_i = np.zeros(len(vector_real), dtype=object)
    common = denominator**order*factorial(order)
    for degree in range(order+1):
        coefficient = denominator**(order-degree)*(factorial(order)//factorial(degree))
        out_r += coefficient*current_r
        out_i += coefficient*current_i
        if degree != order:
            next_r, next_i = -sign*(h@current_i), sign*(h@current_r)
            current_r, current_i = next_r, next_i
    return out_r, out_i, common


def bounded_history_certificate(order=18, time=F(1, 16), site=1):
    """Exact polynomial certificate for V=Z0(t) Zsite Z0(t) Zsite."""
    if type(order) is not int or not 4 <= order <= 32 or type(time) is not F or time <= 0:
        raise ValueError('order4..32 and positive rational source time required')
    if type(site) is not int or not 1 <= site < 8:
        raise ValueError('distinct cube site required')
    _, proof, h, scale, ground = verified_source()
    numerator, denominator = time.numerator, time.denominator*scale
    h = h.astype(object)*numerator
    diagonal = lambda x: np.array([1-2*((mask>>x)&1) for mask in range(256)], dtype=object)
    z0, zy = diagonal(0), diagonal(site)
    real, imag = ground.copy(), np.zeros(256, dtype=object)
    divisor = 1
    # Right-to-left exact target e^(itH)Z0e^(-itH)Zy e^(itH)Z0e^(-itH)Zy.
    for _ in range(2):
        real, imag = zy*real, zy*imag
        real, imag, d = _phase_polynomial_apply(h, real, imag, denominator, -1, order)
        divisor *= d
        real, imag = z0*real, z0*imag
        real, imag, d = _phase_polynomial_apply(h, real, imag, denominator, 1, order)
        divisor *= d
    norm = sum(int(x)**2 for x in ground)
    center_r = F(int(ground@real), divisor*norm)
    center_i = F(int(ground@imag), divisor*norm)
    remainder = (14*time)**(order+1)/factorial(order+1)
    product_error = 4*remainder*(1+remainder)**3
    state_error = F(proof['normalized_ground_vector_error'])
    error = product_error+2*state_error
    abslo, abshi = sqrt_bounds(center_r**2+center_i**2)
    lower = max(F(), 1-(abshi+error)**2)
    upper = min(F(1), 1-max(F(), abslo-error)**2)
    if order >= 18 and time == F(1, 16) and site == 1 and lower <= F(6, 1000):
        raise ArithmeticError('noncentral source history not certified')
    return {'time': str(time), 'site_pair': [0, site], 'Taylor_order': order,
        'source_clock': 'exp(i t H_B), J=1, B=1/4',
        'history': 'V=Z0(t) Zsite Z0(t) Zsite',
        'real_expectation_interval': _interval(center_r, error),
        'imaginary_expectation_interval': _interval(center_i, error),
        'noncentral_history_variance_interval': decimal_bounds(lower, upper, 12),
        'single_unitary_Taylor_remainder': str(remainder),
        'four_evolution_operator_error': str(product_error),
        'ground_vector_error': str(state_error),
        'complete_real_interferometer_probabilities': {
            'plus': _interval((1+center_r)/2, error/2),
            'minus': _interval((1-center_r)/2, error/2)},
        'complete_imaginary_interferometer_probabilities': {
            'plus': _interval((1+center_i)/2, error/2),
            'minus': _interval((1-center_i)/2, error/2)},
        'all_input_reference_safe_controlled_history_trace_error': str(2*product_error+product_error**2),
        'native_controlled_history_tape_executed': False,
        'nonzero_continuum_limit_or_scattering_inferred': False}
