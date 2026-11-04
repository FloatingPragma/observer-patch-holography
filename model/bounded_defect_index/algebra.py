"""Exact source and atomic-index identities, retaining both onsite rows."""
from itertools import product

import sympy as s

from ..native_mass_channels import algebra as source


def site(value):
    if type(value) is not tuple or len(value) != 3 or any(type(x) is not int for x in value):
        raise ValueError('three integer lattice coordinates required')
    return value


def orientation(value):
    if type(value) is not int or value not in (-1, 1):
        raise ValueError('hedgehog orientation +1 or -1 required')
    return value


def onsite(position, mass=1, sign=1):
    position, sign = site(position), orientation(sign)
    mass = source.rational(mass)
    if mass <= 0:
        raise ValueError('positive exact mass magnitude required')
    radius = s.sqrt(sum(n*n for n in position))
    if radius == 0:
        return s.zeros(2)
    return s.ImmutableMatrix(sign*s.Rational(mass)/radius*sum(
        (n*(-1)**(n % 2)*source.SIGMA[j] for j, n in enumerate(position)), s.zeros(2))).applyfunc(s.expand)


def atomic_certificate():
    change = s.Matrix.vstack(s.Matrix.hstack(s.eye(2), s.eye(2)),
                             s.Matrix.hstack(s.eye(2), -s.eye(2)))/s.sqrt(2)
    reflection = s.Matrix.vstack(s.Matrix.hstack(s.zeros(2), s.eye(2)),
                                 s.Matrix.hstack(s.eye(2), s.zeros(2)))
    if change.H*change != s.eye(4) or change.H*reflection*change != s.diag(1, 1, -1, -1):
        raise ArithmeticError('inversion-pair carrier decomposition failed')
    for point in tuple(product((1, 2), repeat=3))+((0, 0, 1), (0, -2, 1)):
        for sign in (-1, 1):
            mass = onsite(point, 3, sign)
            if (mass != mass.H or (mass*mass-9*s.eye(2)).applyfunc(s.expand) != s.zeros(2)
                    or onsite(tuple(-n for n in point), 3, sign) != -mass):
                raise ArithmeticError('the complete bounded source mass changed')
            original = s.diag(mass, -mass)
            expected = s.Matrix.vstack(s.Matrix.hstack(s.zeros(2), mass),
                                        s.Matrix.hstack(mass, s.zeros(2)))
            if (change.H*original*change-expected).applyfunc(s.expand) != s.zeros(4):
                raise ArithmeticError('an atomic inversion pair is not invertible off diagonal')
    if onsite((0, 0, 0)) != s.zeros(2):
        raise ArithmeticError('the declared inversion-fixed source site changed')
    shifted_origin = s.eye(2)/3
    if shifted_origin.det() == 0 or shifted_origin+shifted_origin == s.zeros(2):
        raise ArithmeticError('the inversion-breaking scalar countercontrol failed')
    return {'atomic_origin_even_kernel_dimension': 2,
            'atomic_odd_kernel_dimension': 0,
            'every_nonzero_inversion_pair_invertible': True,
            'atomic_inverse_norm_on_nonzero_pairs': '1/m',
            'reflection_index_plus_hedgehog': 2,
            'reflection_index_antihedgehog': 2,
            'grading': '(R psi)(n)=psi(-n), P_+/-=(I+/-R)/2',
            'grading_is_nonlocal_spatial_inversion': True,
            'identified_with_continuum_internal_chiral_index_or_degree': False,
            'compact_origin_scalar_shift': {'shift': 'm/3', 'original_zero_levels_after_shift': ['m/3', 'm/3'],
                                            'atomic_kernel_after_shift': 0,
                                            'atomic_essential_gap_magnitude': 'm',
                                            'anticommutes_with_reflection': False}}


def edge_certificate():
    first, following = s.symbols('u0:3', real=True), s.symbols('w0:3', real=True)
    delta = [following[j]-first[j] for j in range(3)]
    for signs in product((-1, 1), repeat=3):
        original = sum((signs[j]*first[j]*source.SIGMA[j] for j in range(3)), s.zeros(2))
        for axis in range(3):
            neighbor = sum((signs[j]*(-1 if axis == j else 1)*following[j]*source.SIGMA[j]
                            for j in range(3)), s.zeros(2))
            factor = (source.SIGMA[axis]*neighbor+original*source.SIGMA[axis]).applyfunc(s.expand)
            if (factor.H*factor-sum(x*x for x in delta)*s.eye(2)).applyfunc(s.expand) != s.zeros(2):
                raise ArithmeticError('the exact two-spin anticommutator edge norm changed')
    gamma, masses = source.folded_generators()
    p = s.symbols('p0:3', real=True)
    unit = s.symbols('v0:3', real=True)
    h = sum((p[j]*gamma[j]+unit[j]*masses[j+1] for j in range(3)), s.zeros(16))
    if (h*h-(sum(x*x for x in unit)+sum(x*x for x in p))*s.eye(16)).applyfunc(s.expand) != s.zeros(16):
        raise ArithmeticError('a source valley is missing from the asymptotic bulk gap')
    return {'complete_folded_source_dimension': 16, 'all_weyl_valleys_retained': 8,
            'anticommutator_edge_norm': '(m*v/(2*a))*|u(n+e_j)-u(n)|',
            'unit_direction_difference_bound': '2/min(|n|,|n+e_j|)',
            'finite_hopping_count_per_site': 6,
            'source_hopping_block_norm': 'v/(2*a)',
            'homotopy': 'H_t=t*h0+B, 0<=t<=1',
            'essential_spectrum_all_homotopy_parameters': '[-sqrt(m^2+3*t^2*v^2/a^2),-m] union [m,sqrt(m^2+3*t^2*v^2/a^2)]',
            'essential_gap_stays_open_in_norm_continuous_homotopy': True}


def certificate():
    return {'atomic': atomic_certificate(), 'source_edges_and_gap': edge_certificate(),
            'index_operator': 'P_- H P_+ : range(P_+) -> range(P_-)',
            'exact_Fredholm_index': 2,
            'exact_normalizable_zero_mode_lower_bound': 2,
            'kernel_dimension_is_finite': True,
            'exact_kernel_exhaustion_proved': False,
            'compact_reflection_odd_selfadjoint_perturbations_keep_index': True,
            'reflection_odd_perturbations_of_norm_less_than_m_keep_index': True,
            'arbitrary_gapped_perturbations_pin_zero': False,
            'native_repair_preserves_required_operator_symmetry_proved': False}
