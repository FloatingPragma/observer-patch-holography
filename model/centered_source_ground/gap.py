"""Exact finite centered-cube gap from an untrusted rank-one Gram proposal.

The proposal uses floating eigensolvers and Cholesky. Admission uses only
integer matrix arithmetic, rational Rayleigh quotients and rank-one minmax.
No proposal array, eigensolver residual or observed gap is a retained premise.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import eigsh

from ..interacting_vacuum import global_ground
from ..field_stress.energy import GaussianInteger as GI
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
DIMENSION = 1 << 16
SPIN_EXPECTATION_CEILING = F(1)


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


def _spin_indices(up):
    # The source orders fixed-number addresses by occupied-mode combinations,
    # which differs from increasing bitmask order.
    return np.asarray([sum(1 << j for j in occupied) for occupied in combinations(range(16), 8)
                       if sum(j % 2 == 0 for j in occupied) == up],
                      dtype=np.int64)


def spin_raising_matrix(indices):
    """Exact S+ from the complete ordered N8, Sz0 basis to Sz1.

    In the interleaved hidden CAR ordering, annihilating the occupied down
    mode and creating the adjacent empty up mode gives the same parity
    twice, hence coefficient +1. No spin projection of the trial is made.
    """
    columns = _integer_array(indices, (4900,), 'balanced source addresses')
    if not np.array_equal(columns, _spin_indices(4)):
        raise ValueError('complete canonical balanced source address order required')
    targets = _spin_indices(5)
    row_of = {int(mask): row for row, mask in enumerate(targets)}
    rows, cols = [], []
    for col, value in enumerate(columns):
        mask = int(value)
        for site in range(8):
            if (mask >> (2*site)) & 3 == 2:
                rows.append(row_of[mask ^ (3 << (2*site))])
                cols.append(col)
    matrix = sparse.csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)),
                               shape=(3136, 4900))
    return matrix, targets


def trial_spin_expectation(trial, indices=None):
    """Exact <S^2>=||S+ p||^2/||p||^2 on the full balanced basis."""
    if indices is None:
        indices = global_ground.balanced_indices(8)
    raising, _ = spin_raising_matrix(indices)
    p = _integer_array(trial, (4900,), 'balanced trial')
    # At most eight unit entries can contribute to any output address.
    # Norms are accumulated over Python integers, never fixed-width squares.
    if 8*_maximum(p) > INTMAX:
        raise OverflowError('integer spin-raising accumulation could overflow')
    norm = sum(int(x)**2 for x in p)
    if not norm:
        raise ValueError('nonzero balanced trial required')
    raised = raising@p
    return F(sum(int(x)**2 for x in raised), norm)


def _verify_spin_source(numerator, indices):
    """Bind the CAR ladder to the exact full source, frame and SU(2)."""
    restricted, targets = spin_raising_matrix(indices)
    rows, cols = [], []
    for mask in range(DIMENSION):
        for site in range(8):
            if (mask >> (2*site)) & 3 == 2:
                rows.append(mask ^ (3 << (2*site)))
                cols.append(mask)
    plus = sparse.csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)),
                            shape=(DIMENSION, DIMENSION))
    twice_z = sparse.diags([2*(mask & 0x5555).bit_count()-mask.bit_count()
                           for mask in range(DIMENSION)], dtype=np.int64, format='csr')
    ops = global_ground.cube.operators()
    if ops.scale != 16 or ops.hubbard_H.imag.nnz:
        raise ArithmeticError('exact real hidden-frame numerator16H required')
    h = _source(ops.hubbard_H.real)
    if h.shape != (DIMENSION, DIMENSION):
        raise ArithmeticError('complete full-CAR source shape changed')
    maxh = max((abs(int(x)) for x in h.data), default=0)
    if 2*DIMENSION*max(16, maxh)**2 > INTMAX:
        raise OverflowError('full-source spin commutator could overflow')
    if ((plus@plus.T-plus.T@plus-twice_z).nnz or
            (twice_z@plus-plus@twice_z-2*plus).nnz or
            (h@plus-plus@h).nnz or (h@twice_z-twice_z@h).nnz):
        raise ArithmeticError('exact full-source hidden SU2 identity failed')
    block = _source(numerator)
    if ((plus[targets][:, indices]-restricted).nnz or block.shape != restricted.shape[1:]*2 or
            (h[indices][:, indices] != block).nnz):
        raise ArithmeticError('spin ladder and source block use different ordered addresses')
    hidden = (GI.from_exact(plus+plus.T), GI.from_exact(plus.T-plus).times_i(),
              GI.from_exact(twice_z))
    identity = GI.from_exact(sparse.eye(DIMENSION, dtype=np.int64, format='csr'))
    if (not (ops.frame@ops.frame.adjoint()-identity).is_zero() or
            not (ops.frame@ops.H@ops.frame.adjoint()-ops.hubbard_H).is_zero()):
        raise ArithmeticError('full source hidden frame changed')
    for actual, expected in zip(ops.spin_twice, hidden, strict=True):
        if not (ops.frame@actual@ops.frame.adjoint()-expected).is_zero():
            raise ArithmeticError('CAR ladder differs from the actual source spin frame')


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
    checked = {'trial_rayleigh': rayleigh, 'gram_residual': error,
            'trial_residual_variance': variance,
            'trial_ground_infidelity_upper': infidelity,
            'Temple_ground_energy_lower': temple,
            'ground_energy_strict_upper': GROUND_UPPER,
            'excited_energy_strict_lower': EXCITED_LOWER,
            'neutral_gap_strict_lower': GAP_LOWER}
    # Fraction's internal slots are writable. No returned diagnostic or bound
    # may alias the proof constants or another invocation's values.
    return {name: F(value.numerator, value.denominator) for name, value in checked.items()}


@lru_cache(maxsize=1)
def _certificate_items():
    """Cache immutable scalar claims after the complete source check."""
    original = global_ground.certificate()
    symmetry = cube.source_symmetry()
    sectors = cube.sector_bounds()
    if (original['full_Fock_dimension'] != 65536 or
            original['global_ground_energy'] != '-sqrt(3)/4' or
            symmetry['exact_K_transform'] != 'PH(H-N/2)=H-N/2' or
            sectors['centered_global_ground_number'] != 8 or
            F(sectors['charged_sector_gap_strict_rational_floor']) != F(1, 8)):
        raise ArithmeticError('full-source sector or spin premises changed')
    h = global_ground.source_block(8)
    if h.shape != (4900, 4900):
        raise ArithmeticError('complete balanced N8 source block changed')
    indices = global_ground.balanced_indices(8)
    _verify_spin_source(h, indices)
    trial, factor = propose_witness(h)
    checked = verify_witness(h, trial, factor)
    spin_expectation = trial_spin_expectation(trial, indices)
    if (checked['trial_ground_infidelity_upper'] >= INFIDELITY_CEILING or
            checked['Temple_ground_energy_lower'] <= GROUND_ENCLOSURE_LOWER or
            checked['trial_rayleigh'] >= GROUND_ENCLOSURE_UPPER):
        raise ArithmeticError('exact ground-vector approximation misses its fixed bounds')
    if (spin_expectation >= SPIN_EXPECTATION_CEILING or
            not SPIN_EXPECTATION_CEILING < 2*(1-INFIDELITY_CEILING)):
        raise ArithmeticError('exact trial spin expectation does not exclude a nonsinglet ground')
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
        'normalized_trial_spin_squared_expectation_strict_upper': str(SPIN_EXPECTATION_CEILING),
        'spin_raising_target_dimension': 3136,
        'full_source_SU2_and_spin_frame_checked': True,
        'ground_hidden_spin': 0,
        'full_ground_multiplicity': 1,
        'imported_Lieb_theorem_used_for_gap_or_uniqueness': False,
        'singlet_ground_proof':
            'The unique balanced ground is an S^2 eigenvector. At N8 its spin is integer. '
            'If S>=1, trial overlap>1-1/1000000 implies <S^2>>2(1-1/1000000)>1, '
            'contradicting the exact ||S+ p||^2/||p||^2<1 on the same unprojected trial.',
        'ground_approximation_proof':
            'For rho=<p,Kp>/<p,p>, sigma^2=||(K-rho)p||^2/<p,p> and B=-1053/256, '
            'excited weight<=sigma^2/(B-rho)^2 and E0>=rho-sigma^2/(B-rho). '
            'The expectation error is less than ||A||/500 for every bounded A.',
        'native_ground_preparation_executed': False,
        'balanced_sector_completeness':
            'At N8 every hidden-spin irrep has integer spin and contains Sz=0. '
            'Rank-one minmax and the exact trial spin bound prove the unique ground is a singlet; '
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
