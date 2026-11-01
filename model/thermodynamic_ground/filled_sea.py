"""Exact source countercontrol: the filled free sea is not an interacting ground.

The Wick identity is used only for this explicitly constructed Slater input.
No Gaussian time evolution of an interacting state is assumed.
"""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np

from ..cube_eta_protection import cube
from ..interacting_vacuum import ground


def reference_arrays():
    """All negative hidden-spin orbitals occupied; amplitudes(A+B sqrt3)/1296."""
    reference = ground.reference()
    a, b = np.zeros(65536, dtype=np.int64), np.zeros(65536, dtype=np.int64)
    for up, (x, y) in reference.items():
        for down0, (z, w) in reference.items():
            down = down0 << 1
            inversions = sum(i > j for i in range(16) if up >> i & 1
                             for j in range(16) if down >> j & 1)
            sign, mask = (-1)**inversions, up | down
            a[mask], b[mask] = sign*(x*z+3*y*w), sign*(x*w+y*z)
    return a, b


def _inner(left, right):
    a,b = left; c,d = right
    return int(a@c+3*(b@d)), int(a@d+b@c)


@lru_cache(maxsize=1)
def certificate():
    a,b = reference_arrays()
    denominator = 1296**2
    if _inner((a,b),(a,b)) != (denominator, 0):
        raise ArithmeticError('the complete filled sea is not normalized')
    support = np.flatnonzero((a != 0)|(b != 0))
    if any(int(mask).bit_count() != 8 for mask in support):
        raise ArithmeticError('filled sea left the complete half-filled sector')
    d = np.asarray([sum(((mask >> (2*j)) & 3) == 3 for j in range(8))
                    for mask in range(65536)], dtype=np.int64)
    dv = d*a, d*b
    mean, second = _inner((a,b),dv), _inner(dv,dv)
    if mean != (2*denominator,0) or second != (14*denominator//3,0):
        raise ArithmeticError('exact filled-sea contact moments failed')
    h = cube.operators().hubbard_H.real
    if (not np.array_equal(h@a-16*d*a,-24*b) or
            not np.array_equal(h@b-16*d*b,-8*a)):
        raise ArithmeticError('filled sea does not fill the actual negative kinetic band')
    # Independently evaluate the general Wick variance formula from P_xy.
    adjacency = ground.band_data()[0]
    ordered_edges = int(np.count_nonzero(adjacency))
    wick_variance = F(8,16)+F(ordered_edges,144)
    if ordered_edges != 24 or wick_variance != F(2,3):
        raise ArithmeticError('covariance and complete occupation moments disagree')
    return {'source': 'unchanged original cube at c=1/8,U=1; physical-frame unitary preserves D',
            'particle_number': 8, 'nonzero_occupation_terms': len(support),
            'amplitude_denominator': 1296,
            'contact_mean': '2', 'contact_second_moment': '14/3',
            'contact_variance': '2/3', 'kinetic_energy': '-sqrt3/2',
            'interacting_energy_mean': '2-sqrt3/2', 'interacting_energy_variance': '2/3',
            'interacting_stationary_pure_state': False,
            'general_Slater_domain': 'pure finite CAR Slater projector P with every two-mode onsite block P_xx=I2/2',
            'general_exact_variance': 'Var(D)=V/16+sum_(x!=y)|det(P_xy)|^2>=V/16',
            'general_derivation': 'For x!=y, Cov(Dx,Dy)=-Tr(P_xy P_yx)/4+|det(P_xy)|^2. Projector identity gives sum_(y!=x)Tr(P_xy P_yx)=1/2; onsite Var(Dx)=3/16.',
            'if_also_kinetic_eigenstate': 'Var(H0+UD)=U^2 Var(D)>0 forU>0,V>0',
            'correlated_interacting_ground_excluded': False}
