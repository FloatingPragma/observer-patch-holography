"""Exact original-law join to the half-filled Hubbard charge-gap theorem."""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np
from scipy import sparse

from ..field_symmetry import symmetry
from ..field_stress.energy import GaussianInteger as GI
from ..periodic_momentum import momentum


SPEED = F(1, 8)
REPULSION = F(1)
HOPPING = SPEED/2


def frame(site):
    """G_x=(-1)^x(i sigma_x)^x, a period-four onsite SU2 frame."""
    if type(site) is not int or site < 0: raise ValueError('nonnegative site required')
    return (-1)**(site % 2)*symmetry.onsite_frame((site % 4, 0, 0))


def one_body(length):
    if type(length) is not int or length < 4 or length % 4 or length > 64:
        raise ValueError('explicit source ring length4m between4 and64 required')
    geometry = momentum.PeriodicGeometry((length, 1, 1), SPEED, REPULSION)
    h = np.array((16*momentum.one_particle(geometry)).tolist(), dtype=complex)
    g = sparse.block_diag([sparse.csr_matrix(frame(x)) for x in range(length)], format='csr').toarray()
    target = np.zeros_like(h)
    for x in range(length):
        for spin in (0, 1):
            p, q = 2*x+spin, 2*((x+1) % length)+spin
            target[p, q] -= 1; target[q, p] -= 1
    exact_g, exact_h = GI.from_exact(g), GI.from_exact(h)
    if not (exact_g@exact_h@exact_g.adjoint()-GI.from_exact(target)).is_zero():
        raise ArithmeticError('unchanged periodic Weyl chain does not map to Hubbard hopping')
    return h, target, g


@lru_cache(maxsize=1)
def certificate():
    sigma = np.asarray(symmetry.SIGMA[0])
    for x in range(4):
        g = frame(x)
        if not np.array_equal(g@g.conj().T, np.eye(2)) or g[0, 0]*g[1, 1]-g[0, 1]*g[1, 0] != 1:
            raise ArithmeticError('onsite source frame is not exactly special unitary')
        if not np.array_equal(g@(1j*sigma)@frame(x+1).conj().T, -np.eye(2)):
            raise ArithmeticError('period-four local link theorem failed')
        # Exterior lift on vacuum, both one-particle states and doublon.
        local = np.eye(4, dtype=complex); local[1:3, 1:3] = g
        doublon = np.diag([0, 0, 0, 1])
        if not np.array_equal(local@doublon@local.conj().T, doublon):
            raise ArithmeticError('quartic contact was changed by the frame')
    for length in (4, 8, 12, 16):
        h, _, _ = one_body(length)
        zero_modes = np.zeros((2*length, 4), dtype=np.int64)
        for x in range(length):
            for spin in (0, 1):
                zero_modes[2*x+spin, spin] = 1
                zero_modes[2*x+spin, 2+spin] = (-1)**x
        if np.any(h@zero_modes): raise ArithmeticError('free source lost its four exact zero modes')
    # Independent full256-dimensional CAR witness, beyond the one-body map.
    h, target, g = one_body(4)
    physical = symmetry._lift(h, 8)
    standard = symmetry._lift(target, 8)
    rows, values = [], []
    for mask in range(256):
        occupied = [j for j in range(8) if mask >> j & 1]
        destinations = [int(np.flatnonzero(g[:, j])[0]) for j in occupied]
        target_mask = sum(1 << j for j in destinations)
        inversion = sum(a > b for i, a in enumerate(destinations) for b in destinations[i+1:])
        phase = (-1)**inversion
        for source, destination in zip(occupied, destinations, strict=True): phase *= g[destination, source]
        rows.append(target_mask); values.append(phase)
    transform = GI.from_exact(sparse.coo_matrix((values, (rows, range(256))), shape=(256, 256)))
    contact = GI.from_exact(sparse.diags([16*sum((mask >> (2*x)) & 3 == 3 for x in range(4))
                                        for mask in range(256)], format='csr', dtype=np.int64))
    if not (transform@(physical+contact)@transform.adjoint()-standard-contact).is_zero():
        raise ArithmeticError('full-Fock unchanged-law Hubbard join failed')
    return {'source': 'periodic_momentum.PeriodicGeometry((L,1,1),c=1/8,U=1)',
            'mathematical_lengths': 'L=4m, m>=1; the local rule and period-four frame prove every such ring',
            'explicit_one_body_lengths': [4, 8, 12, 16], 'complete_Fock_witness_dimension': 256,
            'onsite_frame': 'G_x=(-1)^x(i sigma_x)^x',
            'hopping_after_exact_frame': 't=c/2=1/16', 'repulsion_after_exact_frame': 'U=1',
            'quarter_site_frames_checked': 4, 'contact_changed': False,
            'Hubbard_operator': '-t sum(c_x^dagger c_(x+1)+h.c.)+U sum(n_up n_down)',
            'empirical_mass_input': False,
            'native_reference_log_statistics_used_as_new_couplings': False,
            'half_filled_canonical_ground_sector_required': True,
            'half_filling_derived_as_mu_zero_global_ground': False,
            'centered_hamiltonian_silently_substituted': False,
            'zero_interaction_half_filled_charge_gap': '0; chiral spectral pairing and four zero modes on L=4m rings'}
