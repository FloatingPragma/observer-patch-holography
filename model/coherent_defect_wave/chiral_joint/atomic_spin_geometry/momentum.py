"""Actual source translation labels of the rotated cube eigenstates."""
from functools import lru_cache
from itertools import product
import sympy as s
from .cube import full_operators, magnetization_block
from ..atomic_momentum_response.model import shift_spins


def embedding(basis):
    return s.ImmutableSparseMatrix(256, len(basis), {(mask, j): 1 for j, mask in enumerate(basis)})


def translation(axis):
    displacement = tuple(int(i == axis) for i in range(3))
    return s.ImmutableSparseMatrix(256, 256, {
        (shift_spins(mask, 2, displacement), mask): 1 for mask in range(256)})


@lru_cache(maxsize=1)
def certificate():
    plain, twisted, unitary = full_operators()
    operators = [translation(i) for i in range(3)]
    for t in operators:
        assert t.H * t == s.eye(256) and t*t == s.eye(256)
        assert twisted*t == t*twisted
    basis, block = magnetization_block(4)
    polynomial = block**3 + 20*block**2 - 16*block - 576*s.eye(70)
    kernel = polynomial.nullspace()
    assert len(kernel) == 3
    singlets = unitary * embedding(basis) * s.Matrix.hstack(*kernel)
    ground_characters = []
    for t in operators:
        if t*singlets == singlets:
            ground_characters.append(1)
        elif t*singlets == -singlets:
            ground_characters.append(-1)
        else:
            raise AssertionError('ground polynomial carries different translation characters')
    columns = []
    for down in (3, 4, 5):
        basis, block = magnetization_block(down)
        kernel = (block + 16*s.eye(len(basis))).nullspace()
        assert len(kernel) == 1
        columns.append(unitary*embedding(basis)*kernel[0])
    triplet = s.Matrix.hstack(*columns)
    gram = triplet.H*triplet
    assert gram.det() > 0
    reduced = []
    for t in operators:
        small = gram.inv()*triplet.H*t*triplet
        assert t*triplet == triplet*small
        assert small**2 == s.eye(3)
        reduced.append(small)
    records = []
    for signs in product((-1, 1), repeat=3):
        projector = s.eye(3)
        for sign, t in zip(signs, reduced):
            projector = projector*(s.eye(3)+sign*t)/2
        assert projector**2 == projector
        rank = s.trace(projector)
        assert rank in (0, 1, 2, 3)
        records.append({'translation_characters': list(signs),
                        'momenta_in_pi_units': [0 if sign == 1 else 1 for sign in signs],
                        'first_triplet_rank': int(rank)})
    assert sum(row['first_triplet_rank'] for row in records) == 3
    return {'same_atomic_source_translation_used': True,
            'ground_translation_characters': ground_characters,
            'complete_first_triplet_translation_records': records,
            'spin_axis_rotation_is_a_unitary_change_of_coordinates': True,
            'native_clock_or_source_translation_is_replaced': False,
            'relativistic_dispersion_or_named_particle_identification': False}
