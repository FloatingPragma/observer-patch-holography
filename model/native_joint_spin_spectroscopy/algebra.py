"""Exact graph-reflection actions, separately from periodic translations."""
from functools import lru_cache
from fractions import Fraction as F
import numpy as np


def reflections(dimension=3):
    if type(dimension) is not int or not 1 <= dimension <= 3:
        raise ValueError('one through three cube coordinates required')
    return _reflections(dimension)


@lru_cache(maxsize=3)
def _reflections(dimension):
    sites = 1 << dimension
    return tuple(tuple(sum(((mask >> x) & 1) << (x ^ shift)
                           for x in range(sites))
                       for mask in range(1 << sites))
                 for shift in range(sites))


def character_projection(vector, character, dimension=3):
    actions = reflections(dimension)
    if type(character) is not int or not 0 <= character < len(actions):
        raise ValueError('actual character required')
    vector = np.asarray(vector, dtype=object)
    if vector.shape != (len(actions[0]),):
        raise ValueError('complete spin vector required')
    # The returned numerator must be divided by the group order.
    return sum(((-1) ** (character & shift).bit_count()
                * vector[list(action)]
                for shift, action in enumerate(actions)),
               np.zeros_like(vector))


def audit(matrix, scale, dimension=3):
    actions = reflections(dimension)
    h = np.asarray(matrix, dtype=object)
    size = len(actions[0])
    if h.shape != (size, size) or type(scale) is not int or scale <= 0:
        raise ValueError('complete exact source Hamiltonian required')
    for a, pa in enumerate(actions):
        if not np.array_equal(h[np.ix_(pa, pa)], h):
            raise ValueError('graph reflection does not preserve Marshall source')
        for b, pb in enumerate(actions):
            if tuple(pa[pb[i]] for i in range(size)) != actions[a ^ b]:
                raise ValueError('reflection group law failed')
    # Marshall conjugation of reflection equals reflection times global Z
    # precisely when shift interchanges the two sublattices.
    odd = [x for x in range(1 << dimension) if x.bit_count() % 2]
    marshall = [(-1) ** sum((m >> x) & 1 for x in odd) for m in range(size)]
    for shift, action in enumerate(actions):
        for m in range(size):
            if marshall[m] * marshall[action[m]] != (-1) ** (m.bit_count() * (shift.bit_count() % 2)):
                raise ValueError('physical combined-reflection phase failed')
    return {'dimension': dimension, 'sites': 1 << dimension,
            'complete_spin_states': size, 'group_order': len(actions),
            'exact_group_and_source_checks': True,
            'open_cube_labels': 'characters of combined coordinate reflections; not low-momentum periodic dispersion'}


def dimer_common_kernel():
    swap = np.eye(4, dtype=object)[[0, 2, 1, 3]]
    exchange = F(1, 2) * (swap - np.eye(4, dtype=object))
    x0 = np.eye(4, dtype=object)[[1, 0, 3, 2]]
    x1 = np.eye(4, dtype=object)[[2, 3, 0, 1]]
    field = F(-1, 8) * (x0 - x1)
    basis = np.array([[1, 0], [0, 1], [0, 1], [1, 0]], dtype=object)
    if np.any(exchange @ basis) or np.any(field @ basis):
        raise ArithmeticError('transverse probe is not in the common source kernel')
    return {'basis_columns': basis.tolist(), 'exchange_annihilates': True,
            'selected_field_annihilates': True,
            'arbitrary_input_splitting_error_is_zero': False}
