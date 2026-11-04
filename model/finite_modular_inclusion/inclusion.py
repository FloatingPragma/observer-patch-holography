"""Complete CAR boundary witnesses and rational positive-time leakage.

The algebra of a region contains every local CAR operator, with all number
sectors retained. Any such operator commutes with an exterior number
projection. The exact commutator witnesses below therefore separate an
ambient modular derivative from the entire regional algebra.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

import numpy as np
from scipy import sparse

from ..field_stress.energy import GaussianInteger as GI, commutator
from ..field_symmetry.symmetry import operators
from ..full_fock_field.field import FieldGeometry


BETA = F(1, 32)
TIME = F(1, 32)
NORM_BOUND = F(3)
DIMENSION = 256
MODES = 8


def checked_region(sites):
    if (type(sites) is not tuple or not 0 < len(sites) < 4
            or any(type(x) is not int or not 0 <= x < 4 for x in sites)
            or tuple(sorted(set(sites))) != sites):
        raise ValueError('a sorted nonempty proper subset of the four sites is required')
    return tuple(2*x+spin for x in sites for spin in (0, 1))


@lru_cache(maxsize=2)
def source(interaction=1):
    if type(interaction) is not int or interaction not in (0, 1):
        raise ValueError('the retained interaction or its zero control is required')
    data = operators(FieldGeometry((2, 2, 1), interaction=float(interaction)))
    identity = GI.from_exact(sparse.eye(DIMENSION))
    centered = data.H-GI.from_exact(data.scale*sparse.eye(DIMENSION))
    row_bound = int(np.max((abs(centered.real)+abs(centered.imag)).sum(axis=1)))
    if interaction == 1 and F(row_bound, data.scale) != NORM_BOUND:
        raise ArithmeticError('the full-source norm certificate changed')
    if not (data.H-data.H.adjoint()).is_zero():
        raise ArithmeticError('the complete source is not self-adjoint')
    return data.H, data.scale, identity


def number(mode):
    if type(mode) is not int or not 0 <= mode < MODES:
        raise ValueError('an actual source mode is required')
    return GI.from_exact(sparse.diags([(mask >> mode) & 1 for mask in range(DIMENSION)],
                                     dtype=np.int64))


def crossing_modes(sites, interaction=1):
    kept = checked_region(sites)
    h, scale, _ = source(interaction)
    pairs = []
    for p in kept:
        for q in range(MODES):
            if q in kept:
                continue
            real, imag = int(h.real[1 << p, 1 << q]), int(h.imag[1 << p, 1 << q])
            if real or imag:
                if real*real+imag*imag != (scale//16)**2:
                    raise ArithmeticError('crossing one-particle hopping changed')
                pairs.append((p, q, real, imag))
    return tuple(pairs)


def leakage_bound(beta=BETA, time=TIME, generator_error=F(0)):
    """Distance to every regional operator, for one norm-one occupation probe.

    The Gibbs derivative has one forbidden entry of magnitude beta/16.
    Twice the squared centered-H norm bounds its second-derivative remainder.
    The optional full modular-generator error is transferred by Duhamel.
    """
    if any(type(x) not in (int, F) for x in (beta, time, generator_error)):
        raise ValueError('exact rational parameters required')
    beta, time, generator_error = map(F, (beta, time, generator_error))
    if min(beta, time, generator_error) < 0:
        raise ValueError('nonnegative parameters required')
    return beta*time/16-2*NORM_BOUND**2*beta**2*time**2-2*time*generator_error


def deleted_boundary(sites, interaction=1):
    """Delete only hopping entries that change the total regional occupation."""
    kept = checked_region(sites)
    h, _, _ = source(interaction)
    counts = np.array([sum((mask >> p) & 1 for p in kept) for mask in range(DIMENSION)])
    parts = []
    for part in (h.real, h.imag):
        coo = part.tocoo()
        admitted = counts[coo.row] == counts[coo.col]
        parts.append(sparse.coo_matrix((coo.data[admitted],
                                       (coo.row[admitted], coo.col[admitted])),
                                      shape=(DIMENSION, DIMENSION)).tocsr())
    return GI(*parts)


def tensor_sum_certificate(matrix, sites):
    """Exact signed CAR check H=H_A tensor I+I tensor H_B up to scalars."""
    kept = checked_region(sites)
    outside = tuple(p for p in range(MODES) if p not in kept)
    da, db = 1 << len(kept), 1 << len(outside)
    order, signs = [], []
    for mask in range(DIMENSION):
        a = sum(((mask >> p) & 1) << j for j, p in enumerate(kept))
        b = sum(((mask >> p) & 1) << j for j, p in enumerate(outside))
        inversions = sum(((mask >> p) & 1)*((mask >> q) & 1)
                         for p in kept for q in outside if q < p)
        order.append(a*db+b)
        signs.append((-1)**inversions)
    for part in (matrix.real, matrix.imag):
        reordered = np.zeros((DIMENSION, DIMENSION), dtype=np.int64)
        entries = part.tocoo()
        for i, j, value in zip(entries.row, entries.col, entries.data, strict=True):
            reordered[order[i], order[j]] = int(value)*signs[i]*signs[j]
        tensor = reordered.reshape(da, db, da, db)
        trace_b = np.trace(tensor, axis1=1, axis2=3)
        trace_a = np.trace(tensor, axis1=0, axis2=2)
        residual = (DIMENSION*reordered-da*np.kron(trace_b, np.eye(db, dtype=np.int64))
                    -db*np.kron(np.eye(da, dtype=np.int64), trace_a)
                    +np.trace(reordered)*np.eye(DIMENSION, dtype=np.int64))
        if np.any(residual):
            return False
    return True


def region_certificate(sites):
    kept = checked_region(sites)
    h, scale, _ = source()
    pairs = crossing_modes(sites)
    if not pairs:
        raise ArithmeticError('a proper connected-source region has no boundary witness')
    p, q, real, imag = pairs[0]
    derivative = commutator(h, number(p))
    separated = commutator(number(q), derivative)
    # The full matrix, not a fixed-particle truncation, defines both witnesses.
    if separated.is_zero():
        raise ArithmeticError('the modular derivative failed to leave the local CAR algebra')
    if (int(derivative.real[1 << p, 1 << q]),
            int(derivative.imag[1 << p, 1 << q])) != (-real, -imag):
        raise ArithmeticError('the exact derivative matrix element changed')
    if crossing_modes(sites, 0) != pairs:
        raise ArithmeticError('the interaction-deletion control changed a kinetic witness')
    if tensor_sum_certificate(h, sites):
        raise ArithmeticError('crossing hopping unexpectedly factorized')
    if not tensor_sum_certificate(deleted_boundary(sites), sites):
        raise ArithmeticError('boundary deletion did not give an exact tensor sum')
    # Removing all hopping retains the onsite interaction. It is local to
    # either side and hence preserves both complete regional CAR algebras.
    onsite = GI(sparse.diags(h.real.diagonal(), format='csr', dtype=np.int64),
                sparse.csr_matrix((DIMENSION, DIMENSION), dtype=np.int64))
    if not tensor_sum_certificate(onsite, sites):
        raise ArithmeticError('the zero-speed onsite control failed to factorize')
    norm_squared = sum(int(x)**2 for part in (separated.real, separated.imag) for x in part.data)
    return {'sites': list(sites), 'kept_modes': list(kept),
            'regional_algebra_complex_dimension': 4**len(kept),
            'crossing_mode_pairs': [[a, b] for a, b, _, _ in pairs],
            'probe_mode': p, 'exterior_mode': q,
            'hopping_entry': {'real': str(F(real, scale)), 'imag': str(F(imag, scale))},
            'double_commutator_Hilbert_Schmidt_squared': str(F(norm_squared, scale**2)),
            'gibbs_distance_from_entire_regional_algebra_lower': str(leakage_bound()),
            'noninteracting_control_has_same_witness': True,
            'deleted_boundary_is_exact_tensor_sum': True,
            'zero_speed_onsite_interaction_is_exact_tensor_sum': True}


def certificate():
    rows = [region_certificate(sites) for size in (1, 2, 3)
            for sites in combinations(range(4), size)]
    if leakage_bound() != F(23, 524288):
        raise ArithmeticError('the positive-time bound changed')
    return {'all_proper_site_regions': rows, 'full_Fock_dimension': DIMENSION,
            'beta': str(BETA), 'positive_modular_time': str(TIME),
            'centered_H_norm_upper': str(NORM_BOUND),
            'zero_beta_flow_is_identity': True,
            'empty_and_full_algebras_are_invariant': True,
            'finite_dimensional_half_sided_invariance_implies_two_sided': True,
            'geometric_null_translation_constructed': False,
            'continuum_half_sided_inclusion_excluded': False,
            'Einstein_or_physical_particle_gate_accepted': False}
