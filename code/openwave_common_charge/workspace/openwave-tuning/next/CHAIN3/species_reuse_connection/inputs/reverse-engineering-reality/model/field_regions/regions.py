"""CAR restrictions with every local occupancy sector and reference retained.

The maps come from the signed full-Fock split F(K) tensor F(K-complement).
The input N=2 sector is invariant for the field, but its regional image has
N_local=0,1,2. No conditional normalization or projection onto N_local=2 is
performed. Reference indices are always fastest in a density-matrix index.
"""
from __future__ import annotations
from itertools import combinations
from math import comb
import numpy as np
from scipy import sparse


def _modes(values, *, ambient=None):
    if not isinstance(values, (tuple, list)):
        raise ValueError('a sorted list or tuple of distinct integer modes is required')
    values = tuple(values)
    if any(type(p) is not int or p < 0 for p in values) or tuple(sorted(set(values))) != values:
        raise ValueError('modes must be distinct nonnegative integers in ascending order')
    if ambient is not None and not set(values).issubset(ambient):
        raise ValueError('retained modes must belong to the declared ambient modes')
    return values


def regional_basis(modes):
    """Vacuum, ascending singletons, lexicographic pairs, in that order."""
    modes = _modes(modes)
    return ((),) + tuple((p,) for p in modes) + tuple(combinations(modes, 2))


def _matrix(value, dimension):
    if np.ma.isMaskedArray(value):
        raise ValueError('masked matrices are not admitted')
    value = np.asarray(value)
    if value.shape != (dimension, dimension) or value.dtype.kind not in 'ficu' or not np.isfinite(value).all():
        raise ValueError('a finite matrix on the full declared input is required')
    return np.asarray(value, complex)


def _reference_dimension(value):
    if type(value) is not int or value < 1:
        raise ValueError('reference dimension must be a positive integer')
    return value


def _groups(input_basis, kept_modes):
    """Signed full-Fock split restricted to the explicitly supplied basis.

    Moving all retained creations before all exterior creations contributes
    one minus sign for each exterior mode preceding a retained one.
    """
    kept = set(kept_modes)
    local_index = {occupied: j for j, occupied in enumerate(regional_basis(kept_modes))}
    groups = {}
    for j, occupied in enumerate(input_basis):
        inside = tuple(p for p in occupied if p in kept)
        outside = tuple(p for p in occupied if p not in kept)
        sign = (-1)**sum(q < p for q in outside for p in inside)
        groups.setdefault(outside, []).append((j, local_index[inside], sign))
    return groups


def _trace_matrix(density, input_basis, kept_modes, reference_dimension):
    ref = _reference_dimension(reference_dimension)
    density = _matrix(density, len(input_basis)*ref)
    dim = len(regional_basis(kept_modes))*ref
    result = np.zeros((dim, dim), complex)
    for group in _groups(input_basis, kept_modes).values():
        src = np.array([j*ref+r for j, _, _ in group for r in range(ref)])
        dst = np.array([j*ref+r for _, j, _ in group for r in range(ref)])
        sign = np.repeat([sign for _, _, sign in group], ref)
        result[np.ix_(dst, dst)] += density[np.ix_(src, src)]*sign[:, None]*sign[None, :]
    return result


def restrict_density(space, density, kept_modes, reference_dimension=1):
    """Trace exterior modes of an N=2 density, preserving its actual trace.

    The function is linear on arbitrary finite matrices; positivity is a
    property of the Kraus map, not repaired by clipping numerical eigenvalues.
    With a reference, input order is (N=2 basis, reference), output order is
    (regional_basis(kept_modes), reference).
    """
    kept = _modes(kept_modes, ambient=range(space.modes))
    return _trace_matrix(density, tuple(combinations(range(space.modes), 2)), kept,
                         reference_dimension)


def restrict_pure(space, state, kept_modes):
    """Efficient restriction of a pure N=2 state with arbitrary reference.

    A shape (space.dimension, ref) array represents one joint pure vector;
    the reference is retained, not treated as a mixture of input states.
    """
    kept = _modes(kept_modes, ambient=range(space.modes))
    state = space.state(state)
    if state.ndim == 1:
        state = state[:, None]
    ref = state.shape[1]
    if ref < 1:
        raise ValueError('empty reference space is not admitted')
    dim = len(regional_basis(kept))*ref
    result = np.zeros((dim, dim), complex)
    basis = tuple(combinations(range(space.modes), 2))
    for group in _groups(basis, kept).values():
        src = [j for j, _, _ in group]
        dst = np.array([j*ref+r for _, j, _ in group for r in range(ref)])
        sign = np.array([sign for _, _, sign in group])
        vector = (state[src]*sign[:, None]).reshape(-1)
        result[np.ix_(dst, dst)] += np.outer(vector, vector.conj())
    return result


def restrict_regional(density, ambient_modes, kept_modes, reference_dimension=1):
    """Restrict a previously reduced local N<=2 density, without postselection.

    Coherences between different total-number sectors, if supplied, are also
    handled by the full-Fock split. The field's fixed-N input has none after
    tracing its exterior, but this map does not assume that simplification.
    """
    ambient = _modes(ambient_modes)
    kept = _modes(kept_modes, ambient=ambient)
    return _trace_matrix(density, regional_basis(ambient), kept, reference_dimension)


def lift_number_observable(space, observable, kept_modes):
    """Heisenberg inclusion of a local number-preserving observable on N=2.

    This inclusion is a unital star homomorphism on the number-preserving
    regional algebra. Arbitrary local pair-creation operators do not preserve
    the global N=2 sector and are deliberately rejected here.
    """
    kept = _modes(kept_modes, ambient=range(space.modes))
    basis = regional_basis(kept)
    observable = _matrix(observable, len(basis))
    numbers = np.array([len(x) for x in basis])
    if np.any(observable[numbers[:, None] != numbers[None, :]] != 0):
        raise ValueError('the N=2 observable lift requires exact local number preservation')
    rows, cols, values = [], [], []
    input_basis = tuple(combinations(range(space.modes), 2))
    for group in _groups(input_basis, kept).values():
        src = np.array([j for j, _, _ in group])
        dst = np.array([j for _, j, _ in group])
        sign = np.array([sign for _, _, sign in group])
        block = observable[np.ix_(dst, dst)]*sign[:, None]*sign[None, :]
        i, j = np.nonzero(block)
        rows.extend(src[i]); cols.extend(src[j]); values.extend(block[i, j])
    return sparse.coo_matrix((values, (rows, cols)), shape=(space.dimension, space.dimension),
                             dtype=complex).tocsr()


def conditional_expectation(space, matrix, kept_modes):
    """Trace-preserving CP projection onto the regional algebra inside N=2.

    This is the ambient-algebra expectation from CPRestrictionNet, distinct
    from the trace-preserving state restriction, whose codomain is smaller.
    Each exterior sector is replaced by its normalized identity. All exterior
    multiplicities are retained. The map is not multiplicative.
    """
    kept = _modes(kept_modes, ambient=range(space.modes))
    reduced = restrict_density(space, matrix, kept)
    outside = space.modes-len(kept)
    basis = regional_basis(kept)
    for number in range(3):
        indices = [j for j, occupied in enumerate(basis) if len(occupied) == number]
        multiplicity = comb(outside, 2-number) if outside >= 2-number else 0
        if multiplicity:
            reduced[np.ix_(indices, indices)] /= multiplicity
    return lift_number_observable(space, reduced, kept).toarray()


def embed_two_particles(source_space, target_space, mode_map, state):
    """Exterior-square embedding of an injective one-mode map, with CAR sign."""
    if not isinstance(mode_map, (tuple, list)) or len(mode_map) != source_space.modes:
        raise ValueError('one target mode per source mode is required')
    if any(type(p) is not int or not 0 <= p < target_space.modes for p in mode_map) or len(set(mode_map)) != len(mode_map):
        raise ValueError('the mode map must be an injection into the target')
    state = source_space.state(state)
    out = np.zeros((target_space.dimension,)+state.shape[1:], complex)
    for j, (p, q) in enumerate(zip(source_space.first, source_space.second, strict=True)):
        a, b = mode_map[p], mode_map[q]
        sign = 1 if a < b else -1
        out[target_space.index(min(a, b), max(a, b))] = sign*state[j]
    return out


def cube_mode_map(smaller, larger):
    """Same coordinate/spin modes in nested declared Weyl cubes."""
    if smaller.speed != larger.speed or smaller.interaction != larger.interaction:
        raise ValueError('a field refinement must retain the same speed and interaction')
    try:
        return tuple(2*larger.site_index[x]+spin for x in smaller.sites for spin in range(2))
    except KeyError as exc:
        raise ValueError('the target cube must contain every source site') from exc


def two_particle_embedding(source_space, target_space, mode_map):
    """Sparse exterior-square isometry, also usable for operator comparisons."""
    # Reuse the public validator and exact index convention without allocating
    # a dense matrix with one column for every source basis state.
    probe = embed_two_particles(source_space, target_space, mode_map,
                                np.zeros(source_space.dimension, complex))
    del probe
    rows, signs = [], []
    for p, q in zip(source_space.first, source_space.second, strict=True):
        a, b = mode_map[p], mode_map[q]
        rows.append(target_space.index(min(a, b), max(a, b)))
        signs.append(1 if a < b else -1)
    return sparse.coo_matrix((signs, (rows, range(source_space.dimension))),
                             shape=(target_space.dimension, source_space.dimension)).tocsr()
