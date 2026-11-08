"""Finite source controls of polar compression and complete dark readout.

Rational Gram data come from the unchanged complete eight-site source
witness. Floating matrix functions below are independent diagnostics of
finite-section identities; they do not certify an infinite limit.
"""
from fractions import Fraction as F
import numpy as np
from model.spin_covariant_localization import localization


def _inverse_root(matrix):
    values, vectors = np.linalg.eigh(matrix)
    if values[0] <= 0:
        raise ValueError('positive Gram matrix required')
    return (vectors/np.sqrt(values))@vectors.conj().T


def controls():
    source = localization.source_dictionary()
    weights = source['weights']
    # Character and finite-coordinate Fourier bases both come from this
    # same cube. A second orthogonal basis mixes all distinct weights.
    fourier = np.array([[localization.character(k, x) for x in range(8)]
                        for k in range(8)], dtype=float)/np.sqrt(8)
    t = np.diag(np.sqrt([float(w) for w in weights]))@fourier
    full = t@_inverse_root(t.T@t)
    errors=[]
    for n in (1, 2, 4, 8):
        partial = t[:, :n]
        v = partial@_inverse_root(partial.T@partial)
        pi = v@v.T
        # Fixed coordinate partition in the input Fourier space.
        effects=[v@(fourier[:, :n].T@np.diag(np.eye(8)[x])@fourier[:, :n])@v.T
                 for x in range(8)]
        dark=np.eye(8)-pi
        complete=sum(effects)+dark
        gram_residual=np.max(np.abs(v.T@v-np.eye(n)))
        complete_residual=np.max(np.abs(complete-np.eye(8)))
        if max(gram_residual, complete_residual)>1e-11:
            raise ArithmeticError('finite-source Gram compression diagnostic failed')
        for e in effects+[dark]:
            if np.linalg.eigvalsh(e)[0]<-1e-11:
                raise ArithmeticError('finite source detector lost positivity')
        target=full@fourier.T@np.diag(np.eye(8)[0])@fourier@full.T
        errors.append({'rank':n,'complete_dark_instrument':True,
                       'orthogonality_diagnostic_below_1e_minus_11':True,
                       'completeness_diagnostic_below_1e_minus_11':True,
                       'full_rank_agrees_with_original_position': bool(n!=8 or np.linalg.norm(effects[0]-target)<1e-11)})
    # Actual retained source weights are unequal. Compressing first and
    # inverting later differs from compressing the full inverse.
    c=t.T@t
    inverse_then_compress=_inverse_root(c)[:2,:2]
    compress_then_inverse=_inverse_root(c[:2,:2])
    if np.linalg.norm(inverse_then_compress-compress_then_inverse)<1e-3:
        raise ArithmeticError('noncommuting source finite-section countercontrol vanished')
    return {'actual_source_character_weights':[str(w) for w in weights],
            'certified_positive_character_gap':str(source['minimum_weight']),
            'source_ground_certificate_rechecked':True,
            'finite_polar_compressions':errors,
            'compression_does_not_commute_with_inverse_square_root':True,
            'matrix_function_checks_are_floating_diagnostics':True,
            'finite_controls_are_not_proof_of_continuum_convergence':True}
