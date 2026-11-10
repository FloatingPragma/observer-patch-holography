"""Proper cubic rotations on the same fermion, scalar and quantum-link carrier."""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product

import sympy as s

from ...chiral_source.model import PAULI
from ..model import clean
from ..periodic_family.model import Configuration, cells, cell_index, edges, kinetic_terms
from .model import validate_key, reverse_label


@lru_cache(maxsize=1)
def cubic_rotations():
    """One exact SU(2) lift of each of the 24 signed proper permutations."""
    quaternions = []
    for i in range(4):
        quaternions.append(tuple(s.Integer(j == i) for j in range(4)))
    quaternions.extend((s.Rational(1, 2), *(s.Rational(sign, 2) for sign in signs))
                       for signs in product((-1, 1), repeat=3))
    for first, second in combinations(range(4), 2):
        for sign in (-1, 1):
            q = [s.Integer(0)] * 4
            q[first], q[second] = 1 / s.sqrt(2), sign / s.sqrt(2)
            quaternions.append(tuple(q))
    result = []
    for q in quaternions:
        spin = q[0] * s.eye(2) - s.I * sum(
            (q[i + 1] * PAULI[i] for i in range(3)), s.zeros(2))
        spin = s.ImmutableMatrix(spin)
        rotated = [s.simplify(spin * sigma * spin.H) for sigma in PAULI]
        matrix = s.ImmutableMatrix(3, 3, lambda j, i:
            s.simplify(s.trace(PAULI[j] * rotated[i]) / 2))
        assert matrix.det() == 1 and matrix.T * matrix == s.eye(3)
        assert all(entry in (-1, 0, 1) for entry in matrix)
        result.append((matrix, spin))
    assert len({tuple(r) for r, _ in result}) == 24
    return tuple(sorted(result, key=lambda pair: tuple(pair[0])))


def rotate_cell(xyz, rotation, side):
    lifted = tuple(sum(int(rotation[i, j]) * xyz[j] for j in range(3)) for i in range(3))
    winding = sum(value // side for value in lifted)
    return tuple(value % side for value in lifted), s.Integer((-1) ** (winding % 2))


@lru_cache(maxsize=48)
def geometry(rotation, side):
    """Site permutation/AP phases and edge permutation/orientation signs."""
    site_map = tuple((cell_index(side, xyz), phase)
                     for xyz, phase in (rotate_cell(x, rotation, side) for x in cells(side)))
    edge_map = []
    for xyz, axis in edges(side):
        target_axis = next(j for j in range(3) if rotation[j, axis])
        sign = int(rotation[target_axis, axis])
        target, _ = rotate_cell(xyz, rotation, side)
        if sign == -1:
            shifted = list(target)
            shifted[target_axis] = (shifted[target_axis] - 1) % side
            target = tuple(shifted)
        edge_map.append((3 * cell_index(side, target) + target_axis, sign))
    assert sorted(index for index, _ in site_map) == list(range(side ** 3))
    assert sorted(index for index, _ in edge_map) == list(range(3 * side ** 3))
    return site_map, tuple(edge_map)


def hamiltonian_words(config=Configuration()):
    """Complete path words with the native internal indices left symbolic.

    A word is c_left^dagger sigma_axis W_path c_right. Native covariance
    and the full 119-state J identity justify expansion in all internal indices.
    """
    result = defaultdict(lambda: s.Integer(0))
    for term in kinetic_terms(config):
        path = tuple((edge, 1) for edge in term.path)
        key = term.left_cell, term.right_cell, term.axis, path
        adjoint = term.right_cell, term.left_cell, term.axis, tuple(
            (edge, -1) for edge in reversed(term.path))
        result[key] += term.coefficient
        result[adjoint] += s.conjugate(term.coefficient)
    return clean(result)


def rotate_words(words, rotation, side, *, retain_seam=True):
    sites, edge_map = geometry(rotation, side)
    out = defaultdict(lambda: s.Integer(0))
    for (left, right, axis, path), coefficient in words.items():
        target_axis = next(j for j in range(3) if rotation[j, axis])
        sign = int(rotation[target_axis, axis])
        a, phase_a = sites[left]
        b, phase_b = sites[right]
        transformed_path = tuple((edge_map[e][0], direction * edge_map[e][1])
                                 for e, direction in path)
        phase = phase_a * phase_b if retain_seam else 1
        out[a, b, target_axis, transformed_path] += phase * sign * coefficient
    return clean(out)


def axial_rotation_column(key, axis=2, config=Configuration()):
    """The exact pi rotation: spin -i sigma_axis, including every Fock sector.

    Pauli matrices are monomial, so this unitary acts by one signed permutation
    on the complete basis. The native R15 label is left unchanged.
    """
    if axis not in range(3):
        raise ValueError('three rotation axes')
    mask, links, scalars = validate_key(key, config)
    rotation = s.ImmutableMatrix.diag(*(1 if j == axis else -1 for j in range(3)))
    sites, edge_map = geometry(rotation, config.side)
    spin = -s.I * PAULI[axis]
    modes = []
    coefficient = s.Integer(1)
    for mode in range(30 * config.side ** 3):
        if not mask >> mode & 1:
            continue
        cell, local = divmod(mode, 30)
        spin_in, internal = divmod(local, 15)
        spin_out = next(j for j in range(2) if spin[j, spin_in])
        target_cell, phase = sites[cell]
        modes.append(30 * target_cell + 15 * spin_out + internal)
        coefficient *= phase * spin[spin_out, spin_in]
    inversions = sum(a > b for j, a in enumerate(modes) for b in modes[j + 1:])
    coefficient *= (-1) ** inversions
    target_links = [0] * len(links)
    for edge, value in enumerate(links):
        target_edge, sign = edge_map[edge]
        target_links[target_edge] = value if sign == 1 else reverse_label(value)
    target_scalars = [0] * len(scalars)
    for cell, value in enumerate(scalars):
        target_scalars[sites[cell][0]] = value
    target_mask = sum(1 << mode for mode in modes)
    return {(target_mask, tuple(target_links), tuple(target_scalars)): s.expand(coefficient)}
