"""Exact character-frame algebra and source-certified finite localization."""
from fractions import Fraction as F
from pathlib import Path
import numpy as np

from model.native_ordered_ground import certificate, source
from model.native_joint_spin_spectroscopy import algebra
from model.centered_particle_spectrum.certificate import sqrt_bounds, decimal_bounds

ROOT = Path(__file__).resolve().parents[2]


def character(k, x):
    if any(type(a) is not int or not 0 <= a < 8 for a in (k, x)):
        raise ValueError('literal three-bit character and group element required')
    return (-1)**((k & x).bit_count())


def frame_gram(support):
    """Gram of the eight detector vectors for an admitted cyclic support.

    In character coordinates eta_x=(chi_k(x)/sqrt(8))_(k in support).
    Its synthesis operator obeys sum_x|eta_x><eta_x|=I_support even when
    rank<8. The eight effects then form a POVM, not an orthogonal PVM.
    """
    if (type(support) is not tuple or any(type(k) is not int or not 0 <= k < 8 for k in support)
            or len(set(support)) != len(support)):
        raise ValueError('distinct supported characters required')
    gram = tuple(tuple(F(sum(character(k, x)*character(k, y) for k in support), 8)
                       for y in range(8)) for x in range(8))
    for k in support:
        for l in support:
            if F(sum(character(k, x)*character(l, x) for x in range(8)), 8) != int(k == l):
                raise ArithmeticError('complete finite frame lost its identity resolution')
    return gram


def square_interval(lo, hi):
    if lo > hi:
        raise ValueError('ordered real interval required')
    return (F() if lo <= 0 <= hi else min(lo*lo, hi*hi), max(lo*lo, hi*hi))


def finite_fourier_uncertainty():
    """Coordinate/character readout commutators in the full eight-mode code."""
    ident = np.eye(8, dtype=object)
    rows = []
    for axis in range(3):
        x = np.diag([F((n >> axis) & 1) for n in range(8)])
        translation = ident[[n ^ (1 << axis) for n in range(8)]]
        r = F(1, 2)*(ident-translation)
        commutator = x@r-r@x
        if (not np.array_equal(x@x, x) or not np.array_equal(r@r, r)
                or np.trace(commutator) != 0
                or not np.array_equal(commutator@commutator, -F(1, 4)*ident)):
            raise ArithmeticError('finite coordinate/character commutator changed')
        rows.append({'axis': axis, 'commutator_trace': '0',
                     'commutator_squared': '-I/4',
                     'half_commutator_modulus_norm': '1/4',
                     'constant_Heisenberg_lower_bound_inferred': False})
    return rows


def source_dictionary():
    with np.load(ROOT/'model/native_ordered_ground/witness.npz', allow_pickle=False) as archive:
        witness = {k: archive[k] for k in archive.files}
    checked = certificate.verify(witness)
    h, scale = source.matrix()
    group = algebra.audit(h, scale)
    ground = witness['vectors'][:, 0].astype(object)
    norm = sum(int(x)**2 for x in ground)
    probe = ground*np.array([1-2*(mask & 1) for mask in range(256)], dtype=object)
    columns = tuple(algebra.character_projection(probe, k) for k in range(8))
    norms = tuple(sum(int(x)**2 for x in column) for column in columns)
    weights = tuple(F(n, 64*norm) for n in norms)
    if sum(weights) != 1:
        raise ArithmeticError('the complete character weight lost probe normalization')
    for k, a in enumerate(columns):
        if norms[k] <= 0:
            raise ArithmeticError('source character disappeared; use the incomplete-frame construction')
        for l, b in enumerate(columns):
            if sum(int(x)*int(y) for x, y in zip(a, b, strict=True)) != (norms[k] if k == l else 0):
                raise ArithmeticError('source character vectors ceased to be orthogonal')
        for shift, permutation in enumerate(algebra.reflections()):
            if not np.array_equal(a[list(permutation)], character(k, shift)*a):
                raise ArithmeticError('source character covariance failed')
    state_error = F(checked['normalized_ground_vector_error'])
    minimum_weight = min(weights)-2*state_error
    if minimum_weight <= F(1, 100):
        raise ArithmeticError('all-character positive source weight was not certified')
    # ||u/||u||-v/||v|||| <=2||u-v||/||u||.  sqrt(w)>=w for w<=1.
    mode_error = 2*state_error/minimum_weight
    # Eight mode columns give an isometry error <=sqrt(8)*mode_error<3d.
    # Readout rows and the complementary projection give dilation error
    # <=sqrt(5)*3d<9d; this bounds complete half-diamond distance.
    instrument_error = 9*mode_error
    return {'checked': checked, 'h': h, 'scale': scale, 'group': group,
            'columns': columns, 'norms': norms, 'weights': weights,
            'state_error': state_error, 'minimum_weight': minimum_weight,
            'mode_error': mode_error, 'instrument_error': instrument_error}


def report():
    data = source_dictionary()
    columns, norms, weights = data['columns'], data['norms'], data['weights']
    if frame_gram(tuple(range(8))) != tuple(tuple(F(x == y) for y in range(8)) for x in range(8)):
        raise ArithmeticError('full cyclic localization is not a complete orthogonal code')
    roots = [sqrt_bounds(w, 100) for w in weights]
    probabilities = []
    # Ideal detector built from the actual ground, on the actual probed ground.
    # Proposal-mode replacement and ground approximation are both included.
    error = data['instrument_error']+data['state_error']
    for position in range(8):
        lo = hi = F()
        for k, (a, b) in enumerate(roots):
            if character(k, position) > 0:
                lo += a; hi += b
            else:
                lo -= b; hi -= a
        lower, upper = square_interval(lo, hi)
        probabilities.append({'position_bits': format(position, '03b'),
            'original_probability_interval': decimal_bounds(
                max(F(), lower/8-error), min(F(1), upper/8+error), 11)})
    h, scale = data['h'].astype(object), data['scale']
    leakage = []
    for k, column in enumerate(columns):
        acted = h@column
        first = F(sum(int(a)*int(b) for a, b in zip(column, acted, strict=True)), scale*norms[k])
        second = F(sum(int(x)**2 for x in acted), scale*scale*norms[k])
        variance = second-first*first
        # For normalized vectors distance d and ||H||<=14, the variance
        # difference is at most6*14^2*d. The scalar energy origin cancels.
        delta = 6*14**2*data['mode_error']
        lower, upper = max(F(), variance-delta), variance+delta
        leakage.append({'character_bits': format(k, '03b'),
            'squared_source_leakage_interval': decimal_bounds(lower, upper, 11)})
    if F(leakage[7]['squared_source_leakage_interval'][0]) <= F(1, 100):
        raise ArithmeticError('finite dynamics leakage control lost its strict positive witness')
    return {'source': 'selected open eight-site spin cube, J=1, B=1/4; local encoded Z probe',
            'source_group': data['group'], 'full_source_dimension': 256,
            'cyclic_probe_rank': 8, 'dark_complement_rank': 248,
            'all_true_character_weights_strictly_above': '1/100',
            'character_weight_intervals': [decimal_bounds(
                max(F(), w-2*data['state_error']), min(F(1), w+2*data['state_error']), 11)
                for w in weights],
            'complete_detector_frame_halfdiamond_error_upper': str(data['instrument_error']),
            'position_probabilities': probabilities,
            'actual_cyclic_probe_dark_probability': '0',
            'implemented_proposal_dark_probability_upper': str(data['state_error']**2),
            'source_dynamics_leakage': leakage,
            'character111_squared_leakage_strict_lower': '1/100',
            'full_effects_identity_resolution_exact': True,
            'covariance_under_all8_native_group_elements_exact': True,
            'all8_positions_and_dark_complement_retained': True,
            'nonconstant_structure_factor_requires_whitening': max(weights) > 10*min(weights),
            'coordinate_character_uncertainty': finite_fourier_uncertainty(),
            'native_detector_tape_executed': False,
            'cyclic_probe_space_invariant_under_source_H': False,
            'detector_locality_proved': False,
            'graph_reflection_labels_are_continuous_momentum': False,
            'particle_position_mass_or_lorentz_covariance_identified': False}
