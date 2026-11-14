"""Exact spin/eta symmetry of the unchanged full-Fock Weyl/contact action.

The onsite frame exposes real pi-flux hopping. This is an internal symmetry,
not spatial/Lorentz spin or a topological particle assignment. All explicit
checks use Gaussian-integer sparse matrices on every occupation sector.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
import math

import numpy as np
from scipy import sparse

from ..field_stress.energy import GaussianInteger as GI, commutator, sum_exact
from ..full_fock_field.field import FieldGeometry, FullFock


SIGMA = (((0, 1), (1, 0)), ((0, -1j), (1j, 0)), ((1, 0), (0, -1)))
IDENTITY = ((1, 0), (0, 1))
MAX_SITES = 8


def _frozen(value):
    value = np.asarray(value, dtype=complex)
    return np.frombuffer(value.tobytes(), dtype=value.dtype).reshape(value.shape)


def _site(site):
    if (type(site) is not tuple or len(site) != 3 or
            any(type(x) is not int or not 0 <= x <= 4096 for x in site)):
        raise ValueError('three bounded nonnegative integer coordinates required')
    return site


def _axis(axis):
    if type(axis) is not int or axis not in range(3):
        raise ValueError('integer axis 0..2 required')
    return axis


def _times(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def _adjoint(a):
    return tuple(tuple(complex(a[j][i]).conjugate() for j in range(2)) for i in range(2))


def _frame(site):
    out = IDENTITY
    for coordinate, sigma in zip(_site(site), SIGMA, strict=True):
        unitary = tuple(tuple(1j*z for z in row) for row in sigma)
        for _ in range(coordinate % 4):
            out = _times(out, unitary)
    return out


def onsite_frame(site):
    """R_x=(i sigma_0)^x (i sigma_1)^y (i sigma_2)^z, with d_x=R_x c_x."""
    return _frozen(_frame(site))


def link_sign(site, axis):
    """R_x (i sigma_j) R_(x+e_j)^* = link_sign(x,j) I."""
    return (-1)**sum(_site(site)[_axis(axis)+1:])


def hidden_spin_matrix(site, axis):
    frame = _frame(site)
    return _frozen(_times(_times(_adjoint(frame), SIGMA[_axis(axis)]), frame))


def _geometry(geometry, *, explicit=False):
    if type(geometry) is not FieldGeometry:
        raise ValueError('the declared open-box FieldGeometry required')
    if type(geometry.speed) is not float or type(geometry.interaction) is not float:
        raise ValueError('canonical finite field coefficients required')
    FieldGeometry(geometry.sides, geometry.speed, geometry.interaction)
    if explicit and math.prod(geometry.sides) > MAX_SITES:
        raise ValueError('explicit complete-Fock certificate admits at most eight sites')
    return tuple(product(*(range(length) for length in geometry.sides)))


def one_body_frame(geometry=FieldGeometry((2, 2, 1))):
    return sparse.block_diag([sparse.csr_matrix(onsite_frame(site)) for site in _geometry(geometry)], format='csr')


def _scale(operator, coefficient):
    if (type(coefficient) not in (int, complex) or
            not math.isfinite(coefficient.real) or not math.isfinite(coefficient.imag) or
            coefficient.real != int(coefficient.real) or coefficient.imag != int(coefficient.imag)):
        raise ValueError('Gaussian-integer coefficient required')
    if (abs(coefficient.real)+abs(coefficient.imag) > 2**20 or
            operator.maximum()*(abs(coefficient.real)+abs(coefficient.imag)) > 2**30):
        raise OverflowError('exact rescaling numerator bound exceeded')
    return GI.from_exact(operator.matrix()*coefficient)


def _identity(dimension):
    return GI.from_exact(sparse.eye(dimension, format='csr'))


def _annihilator(modes, mode):
    addresses = tuple(mask for mask in range(1 << modes) if mask >> mode & 1)
    signs = [(-1)**(mask & ((1 << mode)-1)).bit_count() for mask in addresses]
    return GI.from_exact(sparse.coo_matrix((signs, ([mask ^ (1 << mode) for mask in addresses], addresses)),
                         shape=(1 << modes, 1 << modes)))


def _lift(matrix, modes):
    """Independent occupation routing of a Gaussian-integer one-body matrix."""
    matrix = sparse.coo_matrix(matrix)
    rows, columns, values = [], [], []
    for p, q, coefficient in zip(matrix.row, matrix.col, matrix.data, strict=True):
        p, q = int(p), int(q)
        if not coefficient:
            continue
        for mask in range(1 << modes):
            if not mask >> q & 1:
                continue
            cleared = mask ^ (1 << q)
            if cleared >> p & 1:
                continue
            sign = (-1)**((mask & ((1 << q)-1)).bit_count()+
                          (cleared & ((1 << p)-1)).bit_count())
            rows.append(cleared | (1 << p)); columns.append(mask); values.append(sign*coefficient)
    return GI.from_exact(sparse.coo_matrix((values, (rows, columns)), shape=(1 << modes,)*2))


def _fock_frame(sites):
    """Full exterior-power frame; every onsite determinant is exactly one."""
    dimension = 1 << (2*len(sites))
    rows, phases = [], []
    frames = [onsite_frame(site) for site in sites]
    for mask in range(dimension):
        target, phase = mask, 1
        for x, frame in enumerate(frames):
            pair = (mask >> (2*x)) & 3
            if pair in (1, 2):
                spin = pair-1
                output = int(np.flatnonzero(frame[:, spin])[0])
                target = (target & ~(3 << (2*x))) | (1 << (2*x+output))
                phase *= frame[output, spin]
        rows.append(target); phases.append(phase)
    return GI.from_exact(sparse.coo_matrix((phases, (rows, range(dimension))), shape=(dimension,)*2))


@dataclass(frozen=True)
class Operators:
    geometry: FieldGeometry
    scale: int
    H: GI
    K: GI
    N: GI
    spin_twice: tuple[GI, GI, GI]
    eta_plus: GI
    eta_minus: GI
    eta_z_twice: GI
    bonds: tuple[GI, ...]
    frame: GI
    hubbard_H: GI


def operators(geometry=FieldGeometry((2, 2, 1))):
    """H/K/bond numerators have denominator scale; named generators do not.

    Binary-rational c,U are admitted with denominator at most 2**16. Exact
    sparse multiplication uses the predecessor's checked integer bounds.
    """
    sites = _geometry(geometry, explicit=True)
    modes, dimension = 2*len(sites), 1 << (2*len(sites))
    c, u = Fraction(geometry.speed), Fraction(geometry.interaction)
    scale = math.lcm((c/2).denominator, (u/2).denominator)
    if scale > 2**16:
        raise ValueError('exact coefficient denominator exceeds the certificate domain')
    if max(scale*c/2, scale*u) > 2**20:
        raise ValueError('exact field coefficient numerator exceeds the certificate domain')
    index = {site: k for k, site in enumerate(sites)}
    one = np.zeros((modes, modes), complex)
    scalar = np.zeros_like(one)
    bonds = []
    for x, site in enumerate(sites):
        for axis in range(3):
            neighbor = list(site); neighbor[axis] += 1
            y = index.get(tuple(neighbor))
            if y is None:
                continue
            block = 1j*int(scale*c/2)*np.asarray(SIGMA[axis])
            bond = np.zeros_like(one)
            bond[2*x:2*x+2, 2*y:2*y+2] = block
            bond[2*y:2*y+2, 2*x:2*x+2] = block.conj().T
            one += bond; bonds.append(_lift(bond, modes))
            scalar[2*x:2*x+2, 2*y:2*y+2] = int(scale*c/2)*link_sign(site, axis)*np.eye(2)
            scalar[2*y:2*y+2, 2*x:2*x+2] = int(scale*c/2)*link_sign(site, axis)*np.eye(2)
    if not np.array_equal(one, scale*geometry.single_h().toarray()):
        raise ArithmeticError('literal Weyl links differ from the actual field action')
    numbers = np.array([mask.bit_count() for mask in range(dimension)])
    contacts = np.array([sum((mask >> (2*x)) & 3 == 3 for x in range(len(sites)))
                         for mask in range(dimension)])
    N = GI.from_exact(sparse.diags(numbers, dtype=np.int64))
    contact = GI.from_exact(sparse.diags(int(scale*u)*contacts, dtype=np.int64))
    H = sum_exact(bonds, dimension)+contact
    actual = GI.from_exact(scale*FullFock(modes).hamiltonian(geometry))
    if not (H-actual).is_zero():
        raise ArithmeticError('independent full occupation routing differs from the field action')
    spins = tuple(_lift(sparse.block_diag([sparse.csr_matrix(hidden_spin_matrix(site, axis)) for site in sites]), modes)
                  for axis in range(3))
    rows, columns, values = [], [], []
    for x, site in enumerate(sites):
        for mask in range(dimension):
            if (mask >> (2*x)) & 3 == 0:
                rows.append(mask | (3 << (2*x))); columns.append(mask); values.append((-1)**sum(site))
    eta = GI.from_exact(sparse.coo_matrix((values, (rows, columns)), shape=(dimension,)*2))
    return Operators(geometry, scale, H, H-_scale(N, int(scale*u/2)), N,
                     spins, eta, eta.adjoint(), GI.from_exact(sparse.diags(numbers-len(sites), dtype=np.int64)),
                     tuple(bonds), _fock_frame(sites), _lift(scalar, modes)+contact)


def eta_product_state(geometry=FieldGeometry((2, 2, 1)), theta=math.pi/2, phase=0.):
    """Product of cos(theta/2)|00>+(-1)^|x| exp(i phase)sin(theta/2)|11>.

    Different even-number sectors are coherent; a phase/charge apparatus is
    supplied. This is an exact eta-tower state, not the interacting ground.
    """
    sites = _geometry(geometry, explicit=True)
    if any(type(value) not in (int, float) or not math.isfinite(value) for value in (theta, phase)):
        raise ValueError('finite real pair angle and phase required')
    empty, pair = math.cos(theta/2), np.exp(1j*phase)*math.sin(theta/2)
    result = np.array([1.], complex)
    for site in sites:
        result = np.kron([empty, 0, 0, (-1)**sum(site)*pair], result)
    return _frozen(result)


def irrep_multiplicities(sites):
    """Full kinematical SU(2)xSU(2) decomposition, not an energy classification."""
    if type(sites) is not int or not 1 <= sites <= 128:
        raise ValueError('integer site count in 1..128 required')
    def weight(spin, eta):
        a, b = sites+eta+spin, sites+eta-spin
        if a % 2 or b % 2 or not 0 <= a <= 2*sites or not 0 <= b <= 2*sites:
            return 0
        return math.comb(sites, a//2)*math.comb(sites, b//2)
    rows = []
    for spin in range(sites+1):
        for eta in range(sites+1):
            count = weight(spin, eta)-weight(spin+2, eta)-weight(spin, eta+2)+weight(spin+2, eta+2)
            if count < 0:
                raise ArithmeticError('negative highest-weight multiplicity')
            if count:
                if (spin+eta-sites) % 2:
                    raise ArithmeticError('spin/eta common-center parity failed')
                rows.append(dict(twice_spin=spin, twice_eta=eta, multiplicity=count))
    if sum((row['twice_spin']+1)*(row['twice_eta']+1)*row['multiplicity'] for row in rows) != 4**sites:
        raise ArithmeticError('full Fock representation dimension lost')
    return tuple(rows)


def frame_certificate():
    """Exact finite residue proof covering every integer-coordinate open box."""
    links = 0
    for site in product(range(4), repeat=3):
        frame = _frame(site)
        if _times(frame, _adjoint(frame)) != IDENTITY or frame[0][0]*frame[1][1]-frame[0][1]*frame[1][0] != 1:
            raise ArithmeticError('onsite special-unitary frame failed')
        for axis in range(3):
            neighbor = list(site); neighbor[axis] += 1
            u = tuple(tuple(1j*z for z in row) for row in SIGMA[axis])
            result = _times(_times(frame, u), _adjoint(_frame(tuple(neighbor))))
            expected = tuple(tuple(link_sign(site, axis)*z for z in row) for row in IDENTITY)
            if result != expected:
                raise ArithmeticError('Weyl-to-real-pi-flux link identity failed')
            expected_spin = (-1)**sum(site[j] for j in range(3) if j != axis)*np.asarray(SIGMA[axis])
            if not np.array_equal(hidden_spin_matrix(site, axis), expected_spin):
                raise ArithmeticError('independent staggered-spin formula failed')
            links += 1
    for a, b in ((0, 1), (0, 2), (1, 2)):
        ua = tuple(tuple(1j*z for z in row) for row in SIGMA[a])
        ub = tuple(tuple(1j*z for z in row) for row in SIGMA[b])
        loop = _times(_times(_times(ua, ub), _adjoint(ua)), _adjoint(ub))
        if loop != ((-1, 0), (0, -1)):
            raise ArithmeticError('central pi-flux plaquette identity failed')
    return dict(coordinate_period=4, residue_sites=64, exact_link_identities=links,
                all_integer_open_boxes=True, all_onsite_determinants_one=True,
                elementary_plaquette_holonomy='-I2',
                scalar_hopping_signs=['(-1)^(y+z)', '(-1)^z', '1'],
                untwisted_odd_periodic_boundaries_included=False)


def certify(geometry=FieldGeometry((2, 2, 1)), *, wards=True):
    if type(wards) is not bool:
        raise ValueError('explicit Boolean Ward-check selection required')
    ops = operators(geometry)
    sites = _geometry(geometry, explicit=True)
    dimension = 1 << geometry.modes
    def zero(value, label):
        if not value.is_zero():
            raise ArithmeticError(label)
    zero(ops.frame@ops.frame.adjoint()-_identity(dimension), 'full Fock frame is not unitary')
    zero(ops.frame@ops.H@ops.frame.adjoint()-ops.hubbard_H, 'full Fock pi-flux equivalence failed')
    for a, b, c in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        zero(commutator(ops.spin_twice[a], ops.spin_twice[b])-_scale(ops.spin_twice[c], 2j), 'spin Lie algebra failed')
    zero(commutator(ops.eta_plus, ops.eta_minus)-ops.eta_z_twice, 'eta Lie algebra failed')
    zero(commutator(ops.eta_z_twice, ops.eta_plus)-_scale(ops.eta_plus, 2), 'eta ladder charge failed')
    zero(commutator(ops.H, ops.eta_plus)-_scale(ops.eta_plus, int(Fraction(geometry.interaction)*ops.scale)), 'eta frequency failed')
    zero(commutator(ops.K, ops.eta_plus), 'grand generator does not preserve eta')
    for spin in ops.spin_twice:
        zero(commutator(ops.H, spin), 'hidden spin is not conserved')
        zero(commutator(spin, ops.eta_plus), 'spin and eta fail to commute')
        zero(commutator(spin, ops.eta_z_twice), 'spin and eta charge fail to commute')
    # Each complete spatial bond, not either component hop alone, commutes
    # with eta. As it kills vacuum, it kills every (eta_plus)^m vacuum.
    for bond in ops.bonds:
        zero(commutator(bond, ops.eta_plus), 'spatial bond does not preserve the eta tower')
        if bond.real[:, 0].nnz or bond.imag[:, 0].nnz:
            raise ArithmeticError('spatial bond fails to kill vacuum')
    checked_wards = 0
    if wards:
        for x, site in enumerate(sites):
            local = [_annihilator(geometry.modes, 2*x+s) for s in range(2)]
            frame = _frame(site)
            d = [sum_exact([_scale(local[j], complex(frame[s][j])) for j in range(2)], dimension)
                 for s in range(2)]
            for axis, spin in enumerate(ops.spin_twice):
                for s in range(2):
                    rhs = sum_exact([_scale(d[t].adjoint(), complex(SIGMA[axis][t][s])) for t in range(2)], dimension)
                    zero(commutator(spin, d[s].adjoint())-rhs, 'spinor creation Ward identity failed')
                    checked_wards += 1
            sign = (-1)**sum(site)
            for s, coefficient in ((0, -sign), (1, sign)):
                zero(commutator(ops.eta_plus, d[s])-_scale(d[1-s].adjoint(), coefficient), 'eta particle/hole Ward identity failed')
                checked_wards += 1
    return dict(sides=list(geometry.sides), sites=len(sites), field_modes=geometry.modes,
                complete_Fock_dimension=dimension, number_sector_dimensions=[math.comb(geometry.modes, n) for n in range(geometry.modes+1)],
                action_denominator=ops.scale, spatial_bonds=len(ops.bonds),
                all_spatial_bonds_annihilate_entire_eta_tower=True,
                complete_CAR_Ward_identities_checked=checked_wards,
                full_Fock_pi_flux_equivalence_exact=True, spin_eta_algebra_exact=True,
                eta_H_frequency=str(Fraction(geometry.interaction)), eta_K_frequency='0',
                central_product=(-1)**len(sites),
                linear_SO4_state_action=(len(sites) % 2 == 0),
                physical_Lorentz_spin_identified=False)


def certificate():
    return dict(schema='oph.full_field_spin_eta_symmetry.v1',
                frame=frame_certificate(),
                full_carrier_checks=[certify(FieldGeometry(sides)) for sides in ((2, 2, 1), (2, 2, 2))],
                full_carrier_irreps={str(v): list(irrep_multiplicities(v)) for v in (1, 2, 4, 8)},
                eta_pair_charge=2, CAR_tensor_operator_twice_spin_twice_eta=[1, 1],
                exact_relations='[H,eta_plus]=U eta_plus; [H-U*N/2,eta_plus]=0',
                odd_site_state_action_is_projective_SO4=True,
                full_field_observable_adjoint_action_descends_to_SO4=True,
                physical_particle_or_topological_defect_identified=False,
                interacting_bulk_gap_or_localized_spinor_proved=False)


report = certificate
