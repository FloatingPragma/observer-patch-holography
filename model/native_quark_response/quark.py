"""Literal native CAR maps into the complete selected quark image."""
from collections import defaultdict
from functools import lru_cache
from math import factorial
from itertools import permutations, product

import sympy as s

from model.coherent_defect_wave.chiral_joint.local_reference.algebra import (
    components, one_body_fock, quark_clifford, transitions,
)
from model.coherent_defect_wave.chiral_joint.local_reference.reduction import (
    coordinate_kernel, occupations, quark_current_matrices, quark_module,
    restrict, symmetric_six, vector_five,
)


@lru_cache(maxsize=None)
def symmetric_power(degree):
    vector = vector_five()
    keys = occupations(degree, 5)
    positions = {key: i for i, key in enumerate(keys)}

    def lift(matrix):
        entries = defaultdict(lambda: s.Integer(0))
        for col, key in enumerate(keys):
            for (i, j), value in matrix.todok().items():
                if key[j]:
                    target = list(key)
                    target[j] -= 1
                    target[i] += 1
                    entries[positions[tuple(target)], col] += key[j] * value
        return s.ImmutableSparseMatrix(len(keys), len(keys), entries)

    diagonal = [s.prod(factorial(n) * vector['metric'][i, i]**n
                       for i, n in enumerate(key)) for key in keys]
    metric = s.ImmutableSparseMatrix(s.diag(*diagonal))
    generators = tuple(lift(t) for t in vector['generators'])
    return dict(keys=keys, metric=metric, generators=generators,
                transitions=tuple(lift(t) for t in vector['transitions']),
                currents=tuple(lift(t) for t in vector['currents']))


@lru_cache(maxsize=1)
def exterior_data():
    rows, f, gamma, rotations, _, _ = quark_clifford()
    generators = gamma + rotations
    two = tuple(mask for mask in range(16) if mask.bit_count() == 2)
    casimir = sum((one_body_fock(t)**2 for t in generators), s.zeros(16))
    vector, _ = coordinate_kernel(casimir.extract(two, two) - 4*s.eye(6))
    assert vector.H * vector == vector_five()['metric']
    return rows, f, generators, two, vector


def elementary_matrix(mode, create):
    entries = {}
    for mask in range(16):
        if bool(mask >> mode & 1) == create:
            continue
        target = mask ^ (1 << mode)
        entries[target, mask] = (-1)**((mask & ((1 << mode)-1)).bit_count())
    return s.ImmutableSparseMatrix(16, 16, entries)


@lru_cache(maxsize=1)
def copy_alignment():
    """Signed one-body permutations align every actual four-mode copy."""
    _, prototype, _, _, _, _ = quark_clifford()
    alignments = []
    for rows in (row for row in components() if len(row) == 4):
        actual = tuple(t.extract(rows, rows) for t in transitions())
        candidates = (s.ImmutableSparseMatrix(4, 4, {(p[j], j): signs[j] for j in range(4)})
                      for p in permutations(range(4)) for signs in product((1, -1), repeat=4))
        unitary = next(u for u in candidates if all(a*u == u*t
                       for a, t in zip(actual, prototype, strict=True)))
        assert unitary.H*unitary == s.eye(4)
        alignments.append((rows, unitary))
    assert len(alignments) == 6
    return tuple(alignments)


@lru_cache(maxsize=None)
def target_module(create):
    """Maps are rational; the target metric contains the split factor 1/6."""
    rows, f, generators, two, vector = exterior_data()
    power = symmetric_power(5)
    source = quark_module()
    source_keys = symmetric_six()['keys']
    positions = {key: i for i, key in enumerate(power['keys'])}
    odd = tuple(mask for mask in range(16) if mask.bit_count() == (3 if create else 1))
    assert len(odd) == 4

    def lifted(matrices, spinor_matrices):
        return tuple(s.ImmutableSparseMatrix(s.kronecker_product(t, s.eye(4))
                     + s.kronecker_product(s.eye(126), u))
                     for t, u in zip(matrices, spinor_matrices, strict=True))

    spinor_generators = tuple(one_body_fock(t).extract(odd, odd) for t in generators)
    full_generators = lifted(power['generators'], spinor_generators)
    casimir = s.ImmutableSparseMatrix(sum((t*t for t in full_generators), s.zeros(504)))
    full_transitions = lifted(power['transitions'],
                             tuple(one_body_fock(t).extract(odd, odd) for t in f))
    full_currents = list(lifted(power['currents'],
                               tuple(one_body_fock(t).extract(odd, odd)
                                     for t in quark_current_matrices())))
    full_currents[3] += (s.Rational(1, 6) if create else -s.Rational(1, 6))*s.eye(504)
    assert all(t == s.diag(*t.diagonal()) for t in full_currents[2:])
    full_metric = s.ImmutableSparseMatrix(s.kronecker_product(power['metric'], s.eye(4))/6)
    maps = []
    for mode in range(4):
        local = elementary_matrix(mode, create).extract(odd, two)*vector
        entries = defaultdict(lambda: s.Integer(0))
        for col, key in enumerate(source_keys):
            for (spinor, i), value in local.todok().items():
                if key[i]:
                    target = list(key)
                    target[i] -= 1
                    entries[4*positions[tuple(target)]+spinor, col] += key[i]*value
        raw = s.ImmutableSparseMatrix(504, 210, entries)*source['embedding']
        assert casimir*raw == s.Rational(95, 2)*raw
        maps.append(s.ImmutableSparseMatrix(raw))

    groups = defaultdict(list)
    for i in range(504):
        groups[(full_currents[2][i, i], full_currents[3][i, i])].append(i)
    entries = {}
    coordinates = []
    weights = []
    image_rank = 0
    for weight, indices in sorted(groups.items()):
        basis, pivots = coordinate_kernel(casimir.extract(indices, indices)
                                          - s.Rational(95, 2)*s.eye(len(indices)))
        if not basis.cols:
            assert all(raw.extract(indices, range(140)) == s.zeros(len(indices), 140)
                       for raw in maps)
            continue
        image = s.Matrix.hstack(*(raw.extract(indices, range(140)) for raw in maps))
        rank = image.to_DM().rank()
        assert rank == basis.cols
        image_rank += rank
        start = len(coordinates)
        entries.update({(indices[i], start+j): value
                        for (i, j), value in basis.todok().items()})
        coordinates.extend(indices[i] for i in pivots)
        weights.extend([weight]*basis.cols)
    basis = s.ImmutableSparseMatrix(504, len(coordinates), entries)
    assert basis.cols == image_rank == 224
    metric = s.ImmutableSparseMatrix(basis.H*full_metric*basis)
    reduced_maps = tuple(s.ImmutableSparseMatrix(raw.extract(coordinates, range(140)))
                         for raw in maps)
    assert all(basis*m == raw for m, raw in zip(reduced_maps, maps, strict=True))
    ts = tuple(restrict(t, basis, coordinates) for t in full_transitions)
    currents = tuple(restrict(t, basis, coordinates) for t in full_currents)
    assert metric*ts[2] == ts[0].H*metric
    assert metric*ts[3] == ts[1].H*metric
    assert all(metric*t == t.H*metric for t in currents)
    assert sum((m.H*metric*m for m in reduced_maps), s.zeros(140)) == 2*source['metric']
    # All ten exact CAR covariance identities are evaluated before compression.
    source_power = symmetric_power(6)
    source_generators = tuple(restrict(t, source['embedding'], source['coordinates'])
                              for t in source_power['generators'])
    for target_t, source_t, fundamental in zip(full_generators, source_generators,
                                               generators, strict=True):
        for alpha, raw in enumerate(maps):
            expected = sum(((fundamental[beta, alpha] if create else -fundamental[alpha, beta])
                            * maps[beta] for beta in range(4)), s.zeros(504, 140))
            assert target_t*raw - raw*source_t == expected
    return dict(create=create, number=13 if create else 11, basis=basis,
                coordinates=tuple(coordinates), metric=metric, transitions=ts,
                currents=currents, maps=reduced_maps,
                weights=tuple((int(2*w), int(6*y)) for w, y in weights),
                native_component=rows, source_dimension=140,
                ambient_dimension=504, image_dimension=224,
                Casimir=s.Rational(95, 2))


if __name__ == '__main__':
    copy_alignment()
    for creation in (False, True):
        data = target_module(creation)
        print(dict(create=creation, number=data['number'],
                   source_dimension=data['source_dimension'],
                   ambient_dimension=data['ambient_dimension'],
                   image_dimension=data['image_dimension'],
                   Casimir=str(data['Casimir'])), flush=True)
