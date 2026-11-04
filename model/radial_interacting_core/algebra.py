"""Exact bounded-radial source symmetries, local frame and parity identities."""
from itertools import product

import sympy as s

from ..bounded_defect_index.algebra import onsite
from ..native_mass_channels.algebra import SIGMA


def source_matrix(radius=1, mass=64):
    if type(radius) is not int or radius != 1:
        raise ValueError('the exact displayed matrix uses the radius-one cube')
    sites = tuple(product(range(-radius, radius+1), repeat=3))
    h = s.zeros(2*len(sites))
    for address, point in enumerate(sites):
        h[2*address:2*address+2, 2*address:2*address+2] = onsite(point, mass)
        for axis in range(3):
            neighbor = tuple(n+(j == axis) for j, n in enumerate(point))
            if neighbor in sites:
                other = sites.index(neighbor)
                h[2*address:2*address+2, 2*other:2*other+2] = s.I*SIGMA[axis]/2
                h[2*other:2*other+2, 2*address:2*address+2] = -s.I*SIGMA[axis]/2
    return s.ImmutableMatrix(h), sites


def symmetries():
    h, sites = source_matrix()
    reflection = s.zeros(27)
    for index, point in enumerate(sites):
        reflection[index, sites.index(tuple(-n for n in point))] = 1
    stagger = s.diag(*((-1)**(sum(point) % 2) for point in sites))
    J = s.I*SIGMA[1]
    R = s.kronecker_product(reflection, s.eye(2))
    C = s.kronecker_product(stagger, J)
    Theta = R*C
    if (h != h.H or R*h*R != -h or C*h.T*C.H != -h
            or Theta*h.conjugate()*Theta.H != h
            or Theta*Theta.conjugate() != -s.eye(54)):
        raise ArithmeticError('the full bounded radial source symmetry failed')
    plain = s.kronecker_product(s.eye(27), J)
    if plain*h.conjugate()*plain.H == h or R*plain*h.conjugate()*plain.H*R == h:
        raise ArithmeticError('an invalid plain or inversion-only time reversal was accepted')
    return {'complete_source_one_body_dimension': 54,
            'mass_profile': 'B_x=m sum_j (-1)^x_j (x_j/|x|) sigma_j; B_0=0',
            'kinetic_source_block': 'i v sigma_j/(2a)',
            'inversion_anticommutes': True,
            'particle_hole_coefficient': 'D J, D_x=(-1)^sum(x), J=i sigma_2',
            'particle_hole_preserves_centered_contact': True,
            'commuting_Kramers_antiunitary': 'R D J complex_conjugation',
            'one_body_Kramers_square': -1,
            'Fock_Kramers_square': '(-1)^N',
            'plain_spin_time_reversal_commutes': False,
            'inversion_times_plain_spin_time_reversal_commutes': False,
            'physical_Lorentz_spin_identification': False}


def local_frame():
    a, b, c, d, r0, r1, r2, r3 = s.symbols('a b c d r0 r1 r2 r3', real=True)
    rho = a*a+b*b+c*c+d*d
    R = r0*s.eye(2)+sum((r*sig for r, sig in zip((r1, r2, r3), SIGMA)), s.zeros(2))
    J = s.I*SIGMA[1]
    for parity in (-1, 1):
        W = s.Matrix([[a+s.I*b, -parity*(c-s.I*d)],
                      [c+s.I*d, parity*(a-s.I*b)]])
        if ((W.H*W-rho*s.eye(2)).applyfunc(s.expand) != s.zeros(2)
                or (W*W.H-rho*s.eye(2)).applyfunc(s.expand) != s.zeros(2)
                or (parity*J*W.conjugate()-W*J).applyfunc(s.expand) != s.zeros(2)
                or s.expand(s.conjugate(W.det())*W.det()-rho*rho) != 0):
            raise ArithmeticError('local Kramers frame does not have scalar weight')
        partner = s.eye(2)-rho*s.eye(2)-R
        pair = (W.H*R*W+W.H*partner*W-(rho-rho*rho)*s.eye(2)).applyfunc(s.expand)
        if pair != s.zeros(2):
            raise ArithmeticError('paired filled-sea matrix identity failed')
    alpha = s.symbols('alpha', real=True)
    sea = (1-alpha)/2
    parity_probability = lambda x, y: 1-x-y+2*x*y
    single = s.expand(parity_probability(sea+alpha, sea))
    empty = s.expand(parity_probability(sea, sea))
    full = s.expand(parity_probability(sea+alpha, sea+alpha))
    if (single != (1-alpha**2)/2 or empty != (1+alpha**2)/2 or full != empty):
        raise ArithmeticError('actual center parity Gaussian formula failed')
    nu, nd, U = s.symbols('n_up n_down U')
    D0 = 1-nu-nd+2*nu*nd
    if s.expand(U*(nu-s.Rational(1, 2))*(nd-s.Rational(1, 2))-U*D0/2+U/4) != 0:
        raise ArithmeticError('actual contact is not the stated parity observable')
    return {'normalized_zero_frame': 'W=(I,-D_Q^-1 V) G^-1/2; W_0=sqrt(alpha) I_2',
            'zero_frame_is_inversion_even': True,
            'local_zero_frame_identity': 'W_x^* W_x=W_x W_x^*=rho_x I_2',
            'rho_normalization': 'sum_x rho_x=1',
            'effective_contact_IPR': 'I=sum_x rho_x^2',
            'paired_sea_identity': 'sum_x W_x^* R_-(x,x) W_x=(1-I) I_2/2',
            'exact_core_contact': 'C+U I (n_1 n_2-(n_1+n_2)/2)',
            'center_sea': 'R_-(0,0)=(1-alpha) I_2/2',
            'center_parity_operator': '(n_center-1)^2',
            'single_core_center_parity': '(1-alpha^2)/2',
            'empty_or_full_core_center_parity': '(1+alpha^2)/2',
            'ideal_charged_minus_ground_center_parity': 'alpha^2',
            'actual_center_contact': '(U/2) D0-U/4'}


def certificate():
    return {'source_symmetries': symmetries(), 'local_zero_frame_and_sea': local_frame()}
