"""Exact N4 ground space and its hidden-frame highest-weight Slater vector.

The original cube at c=1/8,U=1 has one-body eigenvalues +/-sqrt(3)/16.
On N4 its lower kinetic edge is -sqrt(3)/4.  Equality for the repulsive
Hamiltonian requires four negative-band orbitals and zero onsite pairs.
The latter kernel is certified over integers, with an explicit spin-two
multiplet providing the matching rank upper bound.  Ground comparison
against other particle numbers is a separate theorem.
"""
from functools import lru_cache
from itertools import combinations, permutations

import numpy as np

from ..cube_eta_protection import cube


PRIME = 101
BASIS4 = tuple(sum(1 << j for j in sites) for sites in combinations(range(8), 4))
BASIS2 = tuple(sum(1 << j for j in sites) for sites in combinations(range(8), 2))


@lru_cache(maxsize=1)
def band_data():
    """Source-bound integer A and B, with K=[sqrt3 I;-B^T], K*K=6I."""
    ops = cube.operators()
    h = cube.restrict(ops.hubbard_H.real, 1).toarray()
    if ops.scale != 16 or cube.restrict(ops.hubbard_H.imag, 1).nnz:
        raise ArithmeticError('the exact hidden-frame one-body action changed')
    a = h[::2, ::2]
    if (not np.array_equal(h[1::2, 1::2], a) or
            np.any(h[::2, 1::2]) or np.any(h[1::2, ::2])):
        raise ArithmeticError('the source is not the same scalar spin doublet')
    even = tuple(j for j, site in enumerate(cube.GEOMETRY.sites) if sum(site) % 2 == 0)
    odd = tuple(j for j in range(8) if j not in even)
    b = a[np.ix_(even, odd)]
    if (len(even) != 4 or len(odd) != 4 or np.any(a[np.ix_(even, even)]) or
            np.any(a[np.ix_(odd, odd)]) or not np.array_equal(a, a.T) or
            not np.array_equal(b@b.T, 3*np.eye(4, dtype=np.int64))):
        raise ArithmeticError('source bipartite negative-band identity failed')
    # K=kr+sqrt3*ks, in physical hidden-frame site order.
    kr, ks = np.zeros((8, 4), dtype=np.int64), np.zeros((8, 4), dtype=np.int64)
    ks[np.asarray(even)] = np.eye(4, dtype=np.int64)
    kr[np.asarray(odd)] = -b.T
    if (not np.array_equal(a@kr, -3*ks) or not np.array_equal(a@ks, -kr) or
            not np.array_equal(kr.T@kr+3*ks.T@ks, 6*np.eye(4, dtype=np.int64)) or
            np.any(kr.T@ks+ks.T@kr)):
        raise ArithmeticError('negative orbitals lack the claimed eigenvalue or Gram matrix')
    for value in (a, b, kr, ks): value.setflags(write=False)
    return a, b, even, odd, kr, ks


def _annihilate(mask, mode):
    if not mask & (1 << mode): return None
    return mask ^ (1 << mode), (-1)**((mask & ((1 << mode)-1)).bit_count())


def _create(mask, mode):
    if mask & (1 << mode): return None
    return mask | (1 << mode), (-1)**((mask & ((1 << mode)-1)).bit_count())


def contact_matrix(sites=tuple(range(8))):
    """6 times every local c_down*c_up on the orthonormal negative band.

    Rows are (physical site, remaining two band orbitals); columns are the
    full70 four-particle configurations of eight spin/orbital modes.
    Products of K entries are integers because each site lies entirely in
    the rational or sqrt3 part.  No projection after contact is assumed:
    annihilation leaves the remaining orbitals inside the negative band.
    """
    if (type(sites) is not tuple or len(set(sites)) != len(sites) or
            any(type(x) is not int or not 0 <= x < 8 for x in sites)):
        raise ValueError('distinct physical site indices required')
    _, _, _, _, kr, ks = band_data()
    lookup = {mask: j for j, mask in enumerate(BASIS2)}
    result = np.zeros((len(sites)*len(BASIS2), len(BASIS4)), dtype=np.int64)
    for block, x in enumerate(sites):
        products = np.outer(kr[x], kr[x])+3*np.outer(ks[x], ks[x])
        if np.any(np.outer(kr[x], ks[x])+np.outer(ks[x], kr[x])):
            raise ArithmeticError('local pair coefficients are unexpectedly irrational')
        for col, mask in enumerate(BASIS4):
            for a, b in zip(*np.nonzero(products), strict=True):
                first = _annihilate(mask, 2*int(a))
                if first is None: continue
                second = _annihilate(first[0], 2*int(b)+1)
                if second is None: continue
                result[block*len(BASIS2)+lookup[second[0]], col] += int(products[a, b])*first[1]*second[1]
    return result


def spin_ladder(vector, *, lower=True):
    """Exact second-quantized hidden-spin raising/lowering on the70 states."""
    value = np.asarray(vector)
    if value.shape != (70,) or value.dtype.kind not in 'iu' or type(lower) is not bool:
        raise ValueError('70 exact integer coefficients required')
    lookup = {mask: j for j, mask in enumerate(BASIS4)}
    result = np.zeros(70, dtype=np.int64)
    for col, mask in enumerate(BASIS4):
        for a in range(4):
            first = _annihilate(mask, 2*a+(0 if lower else 1))
            if first is None: continue
            second = _create(first[0], 2*a+(1 if lower else 0))
            if second is not None:
                result[lookup[second[0]]] += int(value[col])*first[1]*second[1]
    return result


def multiplet():
    """Five unnormalized exact S=2 vectors, in descending S_z order."""
    highest = np.zeros(70, dtype=np.int64)
    highest[BASIS4.index(sum(1 << (2*a) for a in range(4)))] = 1
    vectors = [highest]
    for _ in range(4): vectors.append(spin_ladder(vectors[-1]))
    if np.any(spin_ladder(vectors[0], lower=False)) or np.any(spin_ladder(vectors[-1])):
        raise ArithmeticError('spin-two highest/lowest weight condition failed')
    for k in range(1, 5):
        if not np.array_equal(spin_ladder(vectors[k], lower=False), k*(5-k)*vectors[k-1]):
            raise ArithmeticError('exact spin-two ladder coefficient failed')
    if [int(v@v) for v in vectors] != [1, 4, 24, 144, 576]:
        raise ArithmeticError('spin multiplet norm failed')
    return tuple(vectors)


def modular_rank(matrix, prime=PRIME):
    """Positive integer-minor witness; no rank upper bound from prime failure."""
    if prime != PRIME: raise ValueError('the retained rank certificate uses prime101')
    value = np.asarray(matrix)
    if value.ndim != 2 or value.dtype.kind not in 'iu': raise ValueError('integer matrix required')
    a = np.array(value % prime, dtype=np.int64, copy=True)
    rows, columns, rank = list(range(len(a))), [], 0
    for col in range(a.shape[1]):
        found = next((j for j in range(rank, len(a)) if a[j, col]), None)
        if found is None: continue
        a[[rank, found]] = a[[found, rank]]
        rows[rank], rows[found] = rows[found], rows[rank]
        a[rank] = a[rank]*pow(int(a[rank, col]), -1, prime) % prime
        for j in range(rank+1, len(a)):
            a[j] = (a[j]-a[j, col]*a[rank]) % prime
        columns.append(col); rank += 1
        if rank == len(a): break
    return rank, rows[:rank], columns


def _det_mod(matrix):
    a = np.array(matrix, dtype=np.int64, copy=True) % PRIME
    if a.shape[0] != a.shape[1]: raise ValueError('square minor required')
    result = 1
    for col in range(len(a)):
        found = next((j for j in range(col, len(a)) if a[j, col]), None)
        if found is None: return 0
        if found != col:
            a[[col, found]] = a[[found, col]]; result = -result
        pivot = int(a[col, col]); result = result*pivot % PRIME
        for j in range(col+1, len(a)):
            a[j] = (a[j]-a[j, col]*pow(pivot, -1, PRIME)*a[col]) % PRIME
    return result % PRIME


def _mul(a, b):
    return a[0]*b[0]+3*a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def _determinant(rows):
    result = [0, 0]
    for permutation in permutations(range(4)):
        inversions = sum(a > b for j, a in enumerate(permutation) for b in permutation[j+1:])
        value = ((-1)**inversions, 0)
        for j, k in enumerate(permutation): value = _mul(value, rows[j][k])
        result[0] += value[0]; result[1] += value[1]
    return tuple(result)


@lru_cache(maxsize=1)
def _reference_items():
    _, _, _, _, kr, ks = band_data()
    result = []
    for sites in combinations(range(8), 4):
        coefficient = _determinant([[(int(kr[x, a]), int(ks[x, a])) for a in range(4)] for x in sites])
        if coefficient != (0, 0): result.append((sum(1 << (2*x) for x in sites), coefficient))
    return tuple(result)


def reference():
    """Sparse hidden-frame N4 state: mask -> (A,B), amplitude (A+B sqrt3)/36.

    This is the all-hidden-up highest-weight state, not a source-prepared
    field state.  The existing physical frame is cube.operators().frame^*.
    """
    return dict(_reference_items())


def reference_arrays():
    """Integer numerator arrays in cube.weights(4) order, denominator36."""
    lookup = {int(mask): j for j, mask in enumerate(cube.weights(4))}
    a, b = np.zeros(1820, dtype=np.int64), np.zeros(1820, dtype=np.int64)
    for mask, (x, y) in _reference_items(): a[lookup[mask]], b[lookup[mask]] = x, y
    return a, b


@lru_cache(maxsize=1)
def certificate():
    _, b, even, odd, _, _ = band_data()
    c = contact_matrix(); vectors = multiplet(); blocks = []
    if any(np.any(c@v) for v in vectors):
        raise ArithmeticError('a claimed ground multiplet contains contact pairs')
    for up in range(4, -1, -1):
        columns = [j for j, mask in enumerate(BASIS4) if (mask & 0x55).bit_count() == up]
        matrix = c[:, columns]
        rank, rows, cols = modular_rank(matrix)
        if rank != len(columns)-1 or not np.any(vectors[4-up][columns]):
            raise ArithmeticError('positive rank witness and explicit kernel do not match')
        if (np.any(matrix@vectors[4-up][columns]) or
                any(vectors[4-up][j] for j in range(70) if j not in columns)):
            raise ArithmeticError('spin-block null vector is not a valid rank upper bound')
        determinant = _det_mod(matrix[np.ix_(rows, cols)])
        if not determinant: raise ArithmeticError('the retained integer minor vanished modulo101')
        blocks.append({'hidden_up_particles': up, 'dimension': len(columns),
                       'rank_lower_bound': rank, 'rank_upper_bound_from_explicit_kernel': len(columns)-1,
                       'kernel_dimension': 1, 'minor_rows': rows, 'minor_columns_within_block': cols,
                       'minor_determinant_mod101': determinant})
    a, d = reference_arrays()
    h = cube.restrict(cube.operators().hubbard_H.real, 4)
    if (not np.array_equal(h@a, -12*d) or not np.array_equal(h@d, -4*a) or
            int(a@a+3*(d@d)) != 36**2 or int(2*(a@d)) != 0):
        raise ArithmeticError('normalized original-action highest-weight ground reference failed')
    if any(any((mask >> (2*x)) & 3 == 3 for x in range(8)) for mask in reference()):
        raise ArithmeticError('highest-weight reference contains a doublon')
    return {'particle_number': 4, 'source_c': '1/8', 'source_U': '1',
            'one_body_eigenvalues': ['-sqrt(3)/16', 'sqrt(3)/16'],
            'negative_one_body_dimension': 8, 'negative_band_N4_dimension': 70,
            'even_sites': list(even), 'odd_sites': list(odd), 'B': b.tolist(),
            'negative_orbital_matrix': 'K=[sqrt3 I_even;-B^T_odd], K^*K=6I4',
            'contact_matrix_shape': list(c.shape),
            'contact_map_normalization': 'matrix equals6 times the stacked local c_down*c_up maps on the normalized band',
            'rank_prime': PRIME, 'spin_blocks': blocks,
            'N4_ground_energy': '-sqrt(3)/4', 'N4_ground_dimension': 5,
            'ground_hidden_spin': 'S=2', 'spin_ladder_squared_norms': [1, 4, 24, 144, 576],
            'ground_proof': 'Kinetic energy is bounded below by4*(-sqrt3/16), with equality exactly in the negative-band wedge. U D>=0 and U>0 force all local pair annihilators to vanish. Nonzero integer minors give rank lower bounds; the five explicit spin-two vectors give matching upper bounds.',
            'highest_weight_reference_frame': 'hidden Hubbard frame; physical state is the existing Fock frame adjoint applied to it',
            'reference_amplitude_denominator': 36,
            'reference_sparse_numerators': [{'mask': mask, 'rational': x, 'sqrt3': y} for mask, (x, y) in _reference_items()],
            'reference_prepared_by_native_source': False,
            'other_particle_number_ground_comparison_proved_here': False}
