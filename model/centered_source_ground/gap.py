"""Exact finite centered-cube gap from an untrusted rank-one Gram proposal.

The proposal uses floating eigensolvers and Cholesky. Admission uses only
integer matrix arithmetic, rational Rayleigh quotients and rank-one minmax.
No proposal array, eigensolver residual or observed gap is a retained premise.
"""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np
from scipy import linalg
from scipy.sparse.linalg import eigsh

from ..interacting_vacuum import global_ground
from . import cube


TRIAL_BITS = 23
FACTOR_BITS = 23
RESIDUAL_CEILING = F(1, 512)
PROPOSAL_CUTOFF = -F(2105, 512)
EXCITED_LOWER = PROPOSAL_CUTOFF-RESIDUAL_CEILING
GROUND_UPPER = -F(1055, 256)
GAP_LOWER = EXCITED_LOWER-GROUND_UPPER
INFIDELITY_CEILING = F(1, 10**6)
GROUND_ENCLOSURE_LOWER = -F(206203, 50000)
GROUND_ENCLOSURE_UPPER = -F(82481, 20000)
INTMAX = int(np.iinfo(np.int64).max)


def _integer_array(value, shape, label):
    if np.ma.isMaskedArray(value):
        raise ValueError('unmasked integer '+label+' required')
    result = np.asarray(value)
    if result.shape != shape or result.dtype != np.dtype(np.int64):
        raise ValueError('int64 '+label+' of the complete declared shape required')
    return result


def _source(numerator):
    # This shared helper sums original sparse duplicates over Python integers
    # and rejects masked and overflowing values before any fixed-width sum.
    h = global_ground._integer_numerator(numerator)
    if not h.shape[0]:
        raise ValueError('nonempty integer source matrix required')
    return h


def _maximum(array):
    return max(abs(int(array.min())), abs(int(array.max())))


def propose_witness(numerator):
    """Return fresh integer trial/factor arrays; neither is trusted evidence."""
    h = _source(numerator)
    n = h.shape[0]
    # This proposer is used on the exact source. Its symmetry check has its
    # own prebound because even a preparatory fixed-width subtraction can wrap.
    maxh = max((abs(int(x)) for x in h.data), default=0)
    if 2*maxh > INTMAX:
        raise OverflowError('source symmetry subtraction could overflow')
    if (h-h.T).nnz:
        raise ValueError('symmetric source Hamiltonian required')
    floating = h.astype(float)/16
    if n <= 2:
        _, vectors = linalg.eigh(floating.toarray())
        vector = vectors[:, 0]
    else:
        initial = np.arange(1, n+1, dtype=float)
        initial /= np.linalg.norm(initial)
        _, vectors = eigsh(floating, k=1, which='SA', tol=1e-11,
                          v0=initial, maxiter=10000)
        vector = vectors[:, 0]
    rounded = np.rint(vector*(1 << TRIAL_BITS))
    if not np.isfinite(rounded).all() or np.max(np.abs(rounded)) >= 1 << 30:
        raise ArithmeticError('trial proposal exceeds the exact integer domain')
    trial = rounded.astype(np.int64)
    # Avoid integer outer-product overflow even in this untrusted proposer.
    trial_float = trial.astype(float)/(1 << TRIAL_BITS)
    shifted = floating.toarray()+np.outer(trial_float, trial_float)
    shifted[np.diag_indices(n)] -= float(4+PROPOSAL_CUTOFF)
    rounded = np.rint(linalg.cholesky(shifted, lower=True)*(1 << FACTOR_BITS))
    if not np.isfinite(rounded).all() or np.max(np.abs(rounded)) >= 1 << 30:
        raise ArithmeticError('Gram proposal exceeds the exact integer domain')
    return trial, np.ascontiguousarray(rounded, dtype=np.int64)


def verify_witness(numerator, trial, factor):
    """Check a complete finite witness; return exact rational diagnostics.

    Here numerator denotes 16H and K=H-4I on N=8. The rank-one
    correction is pp^T/2^(2 TRIAL_BITS), with no normalization assumption.
    The factor need not be triangular or nonsingular.
    """
    h = _source(numerator)
    n = h.shape[0]
    p = _integer_array(trial, (n,), 'trial')
    q = _integer_array(factor, (n, n), 'Gram factor')
    maxp, maxq = _maximum(p), _maximum(q)
    maxh = max((abs(int(x)) for x in h.data), default=0)
    scale = 1 << (2*FACTOR_BITS)
    source_factor = scale//16
    rank_factor = 1 << (2*(FACTOR_BITS-TRIAL_BITS))
    diagonal = -(4+PROPOSAL_CUTOFF)*scale
    if diagonal.denominator != 1:
        raise ArithmeticError('nonintegral fixed Gram scaling')
    diagonal = diagonal.numerator
    accumulation = n*maxq**2+rank_factor*maxp**2+source_factor*maxh+abs(diagonal)
    if accumulation > INTMAX or n*maxh*maxp > INTMAX:
        raise OverflowError('integer Gram or Rayleigh accumulation could overflow')
    if (h-h.T).nnz:
        raise ValueError('symmetric source Hamiltonian required')
    norm = sum(int(x)**2 for x in p)
    if not norm:
        raise ValueError('nonzero trial required')
    hp = h@p
    rayleigh = F(sum(int(x)*int(y) for x, y in zip(p, hp, strict=True)), 16*norm)-4
    if rayleigh >= GROUND_UPPER:
        raise ArithmeticError('exact trial Rayleigh quotient misses the ground upper bound')
    q = np.ascontiguousarray(q)
    maximum = 0
    # Row blocks retain the exact full matrix and keep the verifier's scratch
    # memory small. Every source entry and every factor column participates.
    for start in range(0, n, 128):
        stop = min(start+128, n)
        residual = rank_factor*np.outer(p[start:stop], p)-q[start:stop]@q.T
        residual += h[start:stop].toarray()*source_factor
        residual[np.arange(stop-start), np.arange(start, stop)] += diagonal
        if n*_maximum(residual) > INTMAX:
            raise OverflowError('integer residual row-sum reduction could overflow')
        maximum = max(maximum, int(np.max(np.abs(residual).sum(axis=1))))
    error = F(maximum, scale)
    if error >= RESIDUAL_CEILING:
        raise ArithmeticError('exact rank-one Gram residual exceeds the fixed ceiling')
    variance = F(sum(int(x)**2 for x in hp), 256*norm)-(rayleigh+4)**2
    if variance < 0:
        raise ArithmeticError('negative exact residual variance')
    separation = EXCITED_LOWER-rayleigh
    infidelity = variance/separation**2
    temple = rayleigh-variance/separation
    return {'trial_rayleigh': rayleigh, 'gram_residual': error,
            'trial_residual_variance': variance,
            'trial_ground_infidelity_upper': infidelity,
            'Temple_ground_energy_lower': temple,
            'ground_energy_strict_upper': GROUND_UPPER,
            'excited_energy_strict_lower': EXCITED_LOWER,
            'neutral_gap_strict_lower': GAP_LOWER}


@lru_cache(maxsize=1)
def _certificate_items():
    """Cache immutable scalar claims after the complete source check."""
    original = global_ground.certificate()
    symmetry = cube.source_symmetry()
    lieb = cube.lieb_premises()
    sectors = cube.sector_bounds()
    if (original['full_Fock_dimension'] != 65536 or
            original['global_ground_energy'] != '-sqrt(3)/4' or
            symmetry['exact_K_transform'] != 'PH(H-N/2)=H-N/2' or
            lieb['ground_hidden_spin'] != 0 or lieb['ground_multiplicity'] != 1 or
            sectors['centered_global_ground_number'] != 8 or
            F(sectors['charged_sector_gap_strict_rational_floor']) != F(1, 8)):
        raise ArithmeticError('full-source sector or spin premises changed')
    h = global_ground.source_block(8)
    if h.shape != (4900, 4900):
        raise ArithmeticError('complete balanced N8 source block changed')
    trial, factor = propose_witness(h)
    checked = verify_witness(h, trial, factor)
    if (checked['trial_ground_infidelity_upper'] >= INFIDELITY_CEILING or
            checked['Temple_ground_energy_lower'] <= GROUND_ENCLOSURE_LOWER or
            checked['trial_rayleigh'] >= GROUND_ENCLOSURE_UPPER):
        raise ArithmeticError('exact ground-vector approximation misses its fixed bounds')
    if not 0 < GAP_LOWER < F(1, 8):
        raise ArithmeticError('neutral and charged gap comparison failed')
    result = {
        'full_Fock_dimension': 65536,
        'balanced_N8_dimension': 4900,
        'centered_action': 'K=H-N/2; c=1/8,U=1, open2x2x2 cube',
        'trial_dyadic_bits': TRIAL_BITS,
        'Gram_factor_dyadic_bits': FACTOR_BITS,
        'Gram_residual_strict_upper': str(RESIDUAL_CEILING),
        'rank_one_shifted_Gram_identity':
            'K+2105/512 I+pp^T/2^46=QQ^T/2^46+R; ||R||<=||R||infinity<1/512',
        'ground_energy_strict_upper': str(GROUND_UPPER),
        'N8_excited_energy_strict_lower': str(EXCITED_LOWER),
        'neutral_gap_strict_lower': str(GAP_LOWER),
        'other_charge_sector_gap_strict_lower': '1/8',
        'full_centered_gap_strict_lower': str(GAP_LOWER),
        'ground_energy_enclosure_strict_lower': str(GROUND_ENCLOSURE_LOWER),
        'ground_energy_enclosure_strict_upper': str(GROUND_ENCLOSURE_UPPER),
        'normalized_trial_ground_infidelity_strict_upper': str(INFIDELITY_CEILING),
        'normalized_trial_ground_half_trace_distance_strict_upper': '1/1000',
        'bounded_observable_expectation_error_strict_factor': '1/500',
        'ground_approximation_proof':
            'For rho=<p,Kp>/<p,p>, sigma^2=||(K-rho)p||^2/<p,p> and B=-1053/256, '
            'excited weight<=sigma^2/(B-rho)^2 and E0>=rho-sigma^2/(B-rho). '
            'The expectation error is less than ||A||/500 for every bounded A.',
        'native_ground_preparation_executed': False,
        'balanced_sector_completeness':
            'At N8 every hidden-spin irrep has integer spin and contains Sz=0. '
            'The imported Lieb theorem fixes the unique ground spin to zero; '
            'all other full N8 energies therefore occur in the balanced block above its ground.',
        'rank_one_minmax':
            'A positive rank-one correction makes K+1053/256 I strictly positive; '
            'at most one balanced eigenvalue lies at or below -1053/256.',
        'numerical_eigenvalue_or_residual_used_for_admission': False,
        'integer_proposals_recomputed_and_verified': True,
        'full_source_couplings_changed': False,
        'thermodynamic_gap_proved': False,
        'relativistic_mass_identified': False,
        'localized_topological_protection_proved': False,
    }
    return tuple(sorted(result.items()))


def certificate():
    """Check the actual source; return detached fixed rational claims only."""
    return dict(_certificate_items())
