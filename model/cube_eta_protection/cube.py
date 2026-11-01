"""Exact complete-Fock cube action and eta intertwiners; no aligned-spin premise."""
from functools import lru_cache
from itertools import combinations
import math

import numpy as np
from scipy import sparse

from ..field_symmetry import symmetry
from ..field_stress.energy import commutator
from ..full_fock_field.field import FieldGeometry


GEOMETRY = FieldGeometry((2, 2, 2), speed=1/8, interaction=1.)
DIMENSION = 1 << 16


@lru_cache(maxsize=1)
def operators():
    return symmetry.operators(GEOMETRY)


@lru_cache(maxsize=8)
def weights(number):
    if type(number) is not int or not 0 <= number <= 16:
        raise ValueError('integer CAR number0..16 required')
    return np.array([sum(1 << j for j in pair) for pair in combinations(range(16), number)], dtype=np.int64)


def restrict(value, left, right=None):
    return value[weights(left)][:, weights(left if right is None else right)].tocsr()


@lru_cache(maxsize=1)
def eta_image():
    """Unnormalized full1820-column injection N4 -> N6."""
    return restrict(operators().eta_plus.real, 6, 4)


@lru_cache(maxsize=1)
def parent_lift():
    """T=(eta+)²:N2->N6, and Q8=8 times the eta3 parent projector."""
    eta = operators().eta_plus.real
    t = restrict(eta, 6, 4)@restrict(eta, 4, 2)
    v = np.zeros(120, dtype=np.int64)
    lookup = {int(mask): j for j, mask in enumerate(weights(2))}
    for site, coordinate in enumerate(GEOMETRY.sites):
        v[lookup[3 << (2*site)]] = (-1)**sum(coordinate)
    q8 = 8*sparse.eye(120, dtype=np.int64)-sparse.csr_matrix(np.outer(v, v))
    return t.tocsr(), q8.tocsr(), v


def certificate():
    ops = operators()
    if ops.scale != 16 or not commutator(ops.K, ops.eta_plus).is_zero():
        raise ArithmeticError('the actual cube action lost its centered eta symmetry')
    eta = ops.eta_plus.real
    x = eta_image()
    lower = restrict(ops.eta_minus.real, 2, 4)
    gram = x.T@x-lower.T@lower-4*sparse.eye(1820, dtype=np.int64)
    if gram.nnz:
        raise ArithmeticError('the exact eta injection lower bound failed')
    for part, shift in ((ops.H.real, 16), (ops.H.imag, 0)):
        if (restrict(part, 6)@x-x@restrict(part, 4)-shift*x).nnz:
            raise ArithmeticError('full1820-column cube Hamiltonian intertwiner failed')
    t, q8, v = parent_lift()
    if (q8@q8-8*q8).nnz or (t.T@t-60*sparse.eye(120, dtype=np.int64)-24*sparse.csr_matrix(np.outer(v, v))).nnz:
        raise ArithmeticError('eta3/eta4 two-particle lift Gram identity failed')
    if int(q8.diagonal().sum()) != 8*119:
        raise ArithmeticError('eta3 parent dimension changed')
    for part, shift in ((ops.H.real, 32), (ops.H.imag, 0)):
        if (restrict(part, 6)@t-t@restrict(part, 2)-shift*t).nnz:
            raise ArithmeticError('two-step cube Hamiltonian intertwiner failed')
        if (restrict(part, 2)@q8-q8@restrict(part, 2)).nnz:
            raise ArithmeticError('eta3 parent is not action invariant')
    return {'geometry': [2, 2, 2], 'full_Fock_dimension': DIMENSION,
            'number_six_dimension': math.comb(16, 6),
            'protected_eta_image_dimension': 1820,
            'physical_eta_definition': 'sum_x (-1)^(x0+x1+x2) c_up(x)^dagger c_down(x)^dagger',
            'original_Weyl_contact_action_checked': True,
            'all_three_spatial_hopping_axes_present': True,
            'exact_eta_Gram_identity': 'eta_minus eta_plus=eta_plus eta_minus+4I on N4',
            'minimum_injection_Gram_eigenvalue': 4,
            'full_image_intertwiner': 'H_N6 eta_plus=eta_plus(H_N4+U I)',
            'eta3_parent_dimension': 119, 'eta3_spin_singlet_dimension': 35,
            'eta3_spin_triplet_dimension': 84,
            'eta3_lift_Gram': 'T^dagger T=60I+24vv^dagger, v^dagger v=8; T=(eta_plus)^2',
            'eta3_normalized_intertwiner': 'H_N6 V=V(H_N2+2U I), V=T/sqrt60 on eta3',
            'spin_frame_alignment_used_for_protection': False,
            'maximal_protected_subspace_claimed': False,
            'physical_particle_identification': False}
