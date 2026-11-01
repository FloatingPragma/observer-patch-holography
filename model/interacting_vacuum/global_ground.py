"""Global ground comparison with integer Gram witnesses, not eigensolver admission.

Numerical Cholesky proposes Q.  Only the exact identity involving Q Q^T
and its integer residual proves the lower bound.  No floating eigenvalue,
Cholesky success flag, or nonzero-pivot assumption enters that proof.
"""
from functools import lru_cache
from fractions import Fraction as F

import numpy as np
from scipy import linalg, sparse

from ..cube_eta_protection import cube
from ..field_stress.energy import commutator
from . import ground


BITS = 20
SHIFT = F(27, 64)
CEILING = F(1, 512)
INTMAX = np.iinfo(np.int64).max


def balanced_indices(number):
    if type(number) is not int or not 0 <= number <= 16:
        raise ValueError('integer particle number0..16 required')
    up = (number+1)//2
    return np.asarray([int(mask) for mask in cube.weights(number)
                       if (int(mask) & 0x5555).bit_count() == up], dtype=np.int64)


def source_block(number):
    """Balanced hidden-spin weight of the unchanged full-Fock numerator16H."""
    indices = balanced_indices(number)
    value = cube.operators().hubbard_H
    if value.imag.nnz or cube.operators().scale != 16:
        raise ArithmeticError('source is no longer the exact real numerator16H')
    return value.real[indices][:, indices].tocsr()


def propose_gram(numerator, bits=BITS):
    """Untrusted numerical proposer; its output is always checked over Z."""
    if bits != BITS:
        raise ValueError('this certificate uses20 dyadic bits')
    a = numerator.toarray().astype(float)/16
    a[np.diag_indices_from(a)] += float(SHIFT)
    lower = linalg.cholesky(a, lower=True, overwrite_a=True, check_finite=True)
    q = np.rint(lower*(1 << bits))
    if not np.isfinite(q).all() or np.max(np.abs(q)) >= 1 << 30:
        raise ArithmeticError('untrusted Gram proposal exceeds the exact integer domain')
    return q.astype(np.int64)


def verify_gram(numerator, q, ceiling=CEILING):
    """Exact PSD+residual decomposition; returns the rational residual norm.

    Every product and accumulation is bounded before int64 matrix multiply.
    After forming the exact residual, its absolute row sums are separately
    bounded before their integer reduction.  The matrix need not have full
    rank and Q need not be triangular.
    """
    h = sparse.csr_matrix(numerator, copy=True)
    h.sum_duplicates()
    h.eliminate_zeros()
    q = np.asarray(q)
    if (h.shape[0] != h.shape[1] or q.shape != h.shape or
            h.dtype.kind not in 'iu' or q.dtype != np.dtype(np.int64)):
        raise ValueError('square integer H numerator and int64 Gram factor required')
    if not isinstance(ceiling, F) or ceiling <= 0:
        raise ValueError('positive exact rational residual ceiling required')
    if (h-h.T).nnz:
        raise ValueError('symmetric source Hamiltonian required')
    n = len(q)
    maxq = max((abs(int(x)) for x in (q.min(), q.max())), default=0)
    maxh = max((abs(int(x)) for x in h.data), default=0)
    factor, shift = 1 << (2*BITS-4), 27*(1 << (2*BITS-6))
    if n*maxq**2+factor*maxh+shift > INTMAX:
        raise OverflowError('integer Gram accumulation or residual could overflow')
    residual = q@q.T
    residual *= -1
    residual[np.diag_indices(n)] += shift
    coo = h.tocoo()
    residual[coo.row, coo.col] += coo.data*factor
    maximum = max(abs(int(residual.min())), abs(int(residual.max())))
    if n*maximum > INTMAX:
        raise OverflowError('integer residual row-sum reduction could overflow')
    row_bound = int(np.max(np.abs(residual).sum(axis=1)))
    result = F(row_bound, 1 << (2*BITS))
    if result >= ceiling:
        raise ArithmeticError('exact Gram residual exceeds the declared lower-bound domain')
    return result


@lru_cache(maxsize=1)
def certificate():
    """Full65536-state ground proof; only four nontrivial Gram blocks needed."""
    ground.certificate()
    ops = cube.operators()
    # Conserved exact hidden SU2 makes a balanced weight meet every spin irrep.
    for spin in ops.spin_twice:
        if not commutator(ops.H, spin).is_zero():
            raise ArithmeticError('the source hidden-spin symmetry changed')
    blocks = []
    for number in range(5, 9):
        numerator = source_block(number)
        q = propose_gram(numerator)
        verify_gram(numerator, q)
        # The receipt intentionally binds a fixed rational ceiling, not a
        # machine-dependent untrusted Cholesky proposal or rounded residual.
        blocks.append({'particle_number': number, 'balanced_dimension': len(q),
                       'dyadic_factor_bits': BITS, 'exact_integer_Gram_verified': True,
                       'residual_strict_upper': str(CEILING),
                       'sector_energy_lower': str(-SHIFT-CEILING)})
    # 7/1024 is a conservative certified full charged-sector gap.
    if F(3, 16) <= F(441, 1024)**2:
        raise ArithmeticError('rational square-root separation failed')
    return {'original_action': 'H=dGamma(h)+D; c=1/8,U=1, open2x2x2 cube',
            'full_Fock_dimension': 65536, 'energy_shift_or_chemical_potential': None,
            'global_ground_energy': '-sqrt(3)/4', 'global_ground_dimension': 5,
            'ground_particle_number': 4, 'ground_hidden_spin': '2',
            'all_N_not4_energy_gap_strict_lower': '7/1024',
            'low_number_proof': 'N0..3: D>=0 and dGamma(h)>=-N sqrt3/16',
            'high_number_proof': 'N9..16: D>=N-8 and dGamma(h)>=-N sqrt3/16; lower bound>0',
            'balanced_weight_reduction': 'H commutes with exact hidden SU2. Every integer-spin irrep contains S_z=0 and every half-integer irrep contains S_z=1/2, so each fixed-N eigenvalue occurs in this balanced block.',
            'Gram_blocks': blocks,
            'Gram_identity': 'H=2^(-40) Q Q^T-27/64 I+R; ||R||2<=||R||infinity<1/512',
            'nontrivial_sector_lower_bound': '-217/512',
            'full_ground_order': '-sqrt3/4 < -441/1024 < -217/512',
            'numerical_factor_used_only_as_untrusted_proposal': True,
            'numerical_eigenvalue_admission': False,
            'all_ground_states_are_passive_for_this_finite_H': True,
            'unique_ground_claimed': False,
            'thermodynamic_or_relativistic_vacuum_claimed': False}
