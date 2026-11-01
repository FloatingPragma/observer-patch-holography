"""Exact charged greater/removal measures above an actual interacting ground.

All vectors use the source hidden frame; at the origin its onsite rotation
is the identity, so the four local probes are original physical CAR modes.
No complete pole resolution or native preparation is claimed here.
"""
from functools import lru_cache
from fractions import Fraction as F

import numpy as np

from ..cube_eta_protection import cube, repair
from . import ground


def reference_vector():
    a, b = np.zeros(65536, dtype=np.int64), np.zeros(65536, dtype=np.int64)
    for mask, (x, y) in ground.reference().items():
        a[mask], b[mask] = x, y
    return a, b


def car(vector, mode, create):
    """Exact CAR operation on Q(sqrt3) numerator arrays; no renormalization."""
    if type(mode) is not int or not 0 <= mode < 16 or type(create) is not bool:
        raise ValueError('one of16 original CAR modes and a creation flag required')
    a, b = vector
    if a.shape != (65536,) or b.shape != a.shape or a.dtype != np.int64 or b.dtype != np.int64:
        raise ValueError('two complete-Fock int64 numerator arrays required')
    result = (np.zeros_like(a), np.zeros_like(b))
    for mask in np.flatnonzero((a != 0) | (b != 0)):
        mask = int(mask)
        if bool(mask & (1 << mode)) == create:
            continue
        target, sign = mask ^ (1 << mode), (-1)**((mask & ((1 << mode)-1)).bit_count())
        result[0][target], result[1][target] = sign*a[mask], sign*b[mask]
    return result


def spin_ladder(vector, *, lower):
    result = (np.zeros(65536, dtype=np.int64), np.zeros(65536, dtype=np.int64))
    for site in range(8):
        term = car(car(vector, 2*site+(0 if lower else 1), False),
                   2*site+(1 if lower else 0), True)
        for dst, src in zip(result, term, strict=True): dst += src
    return result


def eta_ladder(vector, *, lower):
    matrix = cube.operators().eta_minus.real if lower else cube.operators().eta_plus.real
    return tuple(matrix@x for x in vector)


def inner(left, right):
    a, b = left; c, d = right
    return int(a@c+3*(b@d)), int(a@d+b@c)


def gap_numerator(vector):
    """16(H-E_G) times a vector: H16+4sqrt3 on the full Fock space."""
    a, b = vector; h = cube.operators().hubbard_H.real
    return h@a+12*b, h@b+4*a


def _equal(left, right):
    return all(np.array_equal(a, b) for a, b in zip(left, right, strict=True))


@lru_cache(maxsize=1)
def spectral_certificate():
    ground.certificate()
    g = reference_vector()
    if inner(g, g) != (1296, 0): raise ArithmeticError('ground normalization failed')
    up_remove, up_add = car(g, 0, False), car(g, 0, True)
    down_remove, v = car(g, 1, False), car(g, 1, True)
    for state in (up_remove, up_add):
        if inner(state, state) != (648, 0): raise ArithmeticError('up spectral sum rule failed')
        a, b = state
        if not _equal(gap_numerator(state), (3*b, a)):
            raise ArithmeticError('up addition/removal line is not exactly sqrt3/16')
    if inner(down_remove, down_remove) != (0, 0) or inner(v, v) != (1296, 0):
        raise ArithmeticError('down spectral sum rule failed')
    s = spin_ladder(spin_ladder(v, lower=False), lower=True)
    e = eta_ladder(eta_ladder(v, lower=True), lower=False)
    # s/5 and e/5 are the exact orthogonal spin5/2 and eta5/2 components.
    if (inner(s, s) != (3240, 0) or inner(e, e) != (3240, 0) or
            inner(s, e) != (0, 0) or inner(v, s) != (648, 0) or inner(v, e) != (648, 0)):
        raise ArithmeticError('exact spin/eta charged residues failed')
    if not _equal(gap_numerator(s), (3*s[1], s[0])):
        raise ArithmeticError('spin pole has the wrong source energy')
    if not _equal(gap_numerator(e), (16*e[0]+3*e[1], 16*e[1]+e[0])):
        raise ArithmeticError('eta pole has the wrong source energy')
    r = tuple(5*x-y-z for x, y, z in zip(v, s, e, strict=True))
    if inner(r, r) != (25920, 0) or inner(r, s) != (0, 0) or inner(r, e) != (0, 0):
        raise ArithmeticError('charged residual measure is not an orthogonal weight4/5')
    gap = gap_numerator(v)
    mean, second = inner(v, gap), inner(gap, gap)
    if mean != (10368, 0) or second != (169776, 0):
        raise ArithmeticError('source charged first/second moment failed')
    # Independent exact moments of the residual probability measure.
    rmean, rsecond = inner(r, gap_numerator(r)), inner(gap_numerator(r), gap_numerator(r))
    if (tuple(F(x, 16*25920) for x in rmean) != (F(1, 2), -F(1, 64)) or
            tuple(F(x, 256*25920) for x in rsecond) != (F(131, 256), -F(1, 64))):
        raise ArithmeticError('residual normalized moments failed')
    return {'reference': 'normalized highest-weight global ground G; N4,S2,M2; E_G=-sqrt3/4',
            'probe_site': [0, 0, 0], 'probe_frame': 'original physical CAR equals hidden frame at this site',
            'lambda': 'sqrt3/16',
            'up_removal': {'particle_number': 3, 'total_weight': '1/2', 'gap': 'sqrt3/16'},
            'up_addition': {'particle_number': 5, 'total_weight': '1/2', 'gap': 'sqrt3/16'},
            'down_removal': {'particle_number': 3, 'total_weight': '0'},
            'down_addition': {'particle_number': 5, 'total_weight': '1',
                              'measure': '(1/10)delta_lambda+(1/10)delta_(1+lambda)+(4/5)mu_residual',
                              'resolved_symmetry_components': [{'gap': 'sqrt3/16', 'weight': '1/10', 'selection': 'S=5/2'},
                                                     {'gap': '1+sqrt3/16', 'weight': '1/10', 'selection': 'eta=5/2'}],
                              'residual_weight': '4/5', 'residual_measure_normalized': True,
                              'positive_support_strict_lower': '7/1024',
                              'gap_mean': '1/2', 'gap_second_moment': '131/256', 'gap_variance': '67/256',
                              'normalized_residual_gap_mean': '1/2-sqrt3/64',
                              'normalized_residual_gap_second_moment': '131/256-sqrt3/64',
                              'normalized_residual_gap_variance': '1069/4096',
                              'all_residual_atoms_resolved': False},
            'sum_rules': '||c_up G||^2+||c_up^dagger G||^2=1; ||c_down G||^2+||c_down^dagger G||^2=1',
            'charge_relative_to_G': 'addition+1; removal-1 under conserved model N',
            'CAR_parity_relative_to_G': 'odd',
            'positive_support_uses_global_ground_certificate': True,
            'finite_gaps_are_relativistic_rest_masses': False,
            'native_reference_preparation_or_detector_executed': False}


@lru_cache(maxsize=1)
def repair_certificate():
    """Passive outside-service protection for all N3/N4/N5, including coherences."""
    a, b = repair.exact_embedding()
    support = set(map(int, a.nonzero()[0])) | set(map(int, b.nonzero()[0]))
    if any(mask.bit_count() != 6 for mask in support):
        raise ArithmeticError('the admitted repair support is not wholly N6')
    count = 0
    for law in ('C', 'E'):
        for seam in repair.SEAMS:
            matrices = repair.matrices(law, seam)
            if not np.array_equal(sum(m.T@m for _, m in matrices), 2*np.eye(24, dtype=np.int64)):
                raise ArithmeticError('an admitted repair lost complete outcome accounting')
            count += len(matrices)+1
    if count != 130: raise ArithmeticError('complete repair outcome count changed')
    return {'source_active_support_particle_number': 6,
            'passive_protected_numbers': [3, 4, 5],
            'passive_protected_dimension': 560+1820+4368,
            'complete_outcomes_including_zero_and_outside': count,
            'exact_statement': 'V^dagger P_(N3+N4+N5)=0, all native Kraus words vanish there and QP=P; every cross-number coherence and arbitrary reference is preserved.',
            'Hamiltonian_preserves_this_sum': True,
            'arbitrary_finite_interleaving': 'H evolves the state normally; each admitted repair acts as identity on this sum, with the outside record certain.',
            'protection_kind': 'passive outside-service identity of the declared full-Fock completion',
            'active_error_correction_or_arbitrary_repair_protection': False}
