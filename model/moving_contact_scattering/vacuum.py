"""Vacuum obstruction for the additive reflected Fock generator.

Negative energy per separated pair rules out an empty-Fock ground, not
stability of matter at fixed N. The fixed-particle lower bound is retained.
"""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from ..contact_continuum.continuum import rational
from ..moving_bound_pair.limit import moving_pole_certificate


def negative_packet_certificate(spacing, gamma):
    a = rational(spacing, 'spacing', positive=True)
    g = rational(gamma, 'inverse length', positive=True)
    if a*g > F(1, 4096):
        raise ValueError('negative packet certificate requires a*gamma<=1/4096')
    cap = g/4
    pole = moving_pole_certificate(a, g, cap)
    energy = 3*cap*cap/4+pole['center_of_mass_edge_error_upper']-pole['pole_offset_lower']
    if energy >= -g*g/8:
        raise ArithmeticError('uniformly negative bound packet was lost')
    norm_error = a*a*g*g/768
    error = 24*norm_error/(a*a)
    if energy+error >= -g*g/16:
        raise ArithmeticError('compact approximation lost the negative variational margin')
    return {'spacing': a, 'gamma': g, 'momentum_cap': cap,
            'bound_wavepacket_energy_upper': energy,
            'two_particle_operator_norm_upper': 12/(a*a),
            'sufficient_normalized_compact_vector_norm_error': norm_error,
            'compact_energy_approximation_error_upper': error,
            'compact_pair_energy_upper': -g*g/16,
            'm_disjoint_pairs_energy_upper': '-m*gamma^2/16',
            'compact_support_radius_effectively_computed': False,
            'compact_approximation_justification': 'finite-support vectors are dense and J_a on N=2 is bounded'}


def fixed_number_lower_bound(number, interaction):
    if type(number) is not int or number < 0:
        raise ValueError('nonnegative integer particle number required')
    u = rational(interaction, 'onsite attraction magnitude')
    return -u*(number//2)


def affine_shell_certificate():
    radius, alpha, offset, gamma, mu, q, q1, q0, speed, mass = s.symbols(
        'radius alpha offset gamma mu q q1 q0 speed mass', real=True)
    energy = alpha*(radius**2/4-gamma**2+2*mu)+offset
    residual = s.Poly(s.expand(energy**2-speed**2*(q*radius**2+q1*radius+q0)
                              -mass**2*speed**4), radius)
    if residual.coeff_monomial(radius**4) != alpha**2/16:
        raise ArithmeticError('the quadratic bound band lost its quartic Lorentz residual')
    return {'energy_map': 'E=alpha*(|P|^2/4-gamma^2+2*mu)+offset; alpha != 0',
            'momentum_map': 'p=L P+p0 with fixed invertible real L and fixed p0',
            'quartic_coefficient_on_P=r_e1': 'alpha^2/16 > 0',
            'exact_Lorentz_shell_on_any_open_momentum_set': False,
            'nonrelativistic_low_momentum_approximation_excluded': False,
            'arbitrary_nonlinear_dictionary_excluded': False}


def _move(mask, first, second):
    if not (mask >> second) & 1:
        return None
    sign = (-1)**((mask & ((1 << second)-1)).bit_count())
    target = mask ^ (1 << second)
    if (target >> first) & 1:
        return None
    sign *= (-1)**((target & ((1 << first)-1)).bit_count())
    return target | (1 << first), sign


@lru_cache(maxsize=2)
def finite_operator(sites=4):
    """Exact small reflected Weyl/contact source, including inter-core hops.

    The four-site oriented ring has two spatial cores with interleaved mode
    addresses. The two-site restriction has the single internal core link.
    This checks the CAR product identity; it is not a critical-family pole.
    """
    if type(sites) is not int or sites not in (2, 4):
        raise ValueError('the independent control uses two or four sites')
    links = ((0, 1),) if sites == 2 else ((0, 2), (2, 1), (1, 3), (3, 0))
    modes = 2*sites
    entries = {}
    for mask in range(1 << modes):
        doubles = sum((mask >> (2*x)) & 3 == 3 for x in range(sites))
        diagonal = 4*mask.bit_count()-10*doubles
        if diagonal:
            entries[mask, mask] = diagonal
        for x, y in links:
            for spin in (0, 1):
                first, second = 2*x+spin, 2*y+1-spin
                for p, q, coefficient in ((first, second, -s.I),
                                          (second, first, s.I)):
                    moved = _move(mask, p, q)
                    if moved is not None:
                        target, sign = moved
                        entries[target, mask] = entries.get((target, mask), 0)+sign*coefficient
    return s.ImmutableSparseMatrix(1 << modes, 1 << modes, entries)


def wedge_vector(left, right, *, graded=True):
    """Map two four-mode vectors to disjoint interleaved sets of eight modes."""
    core_a, core_b = (0, 1, 4, 5), (2, 3, 6, 7)
    entries = {}
    for a, x in left.items():
        for b, y in right.items():
            addresses = ([core_a[j] for j in range(4) if (a >> j) & 1]+
                         [core_b[j] for j in range(4) if (b >> j) & 1])
            inversions = sum(p > q for i, p in enumerate(addresses) for q in addresses[i+1:])
            mask = sum(1 << j for j in addresses)
            entries[mask, 0] = (-1)**(inversions if graded else 0)*x*y
    return s.ImmutableSparseMatrix(256, 1, entries)


def finite_additivity_certificate():
    left, right = {3: 1, 12: 2, 5: -s.I}, {3: 1, 12: 2, 5: -s.I}
    a = s.ImmutableSparseMatrix(16, 1, {(j, 0): v for j, v in left.items()})
    b = s.ImmutableSparseMatrix(16, 1, {(j, 0): v for j, v in right.items()})
    product = wedge_vector(left, right)
    wrong = wedge_vector(left, right, graded=False)
    norms = [(v.H*v)[0] for v in (a, b, product)]
    energies = [s.simplify((v.H*h*v)[0]/norm) for v, h, norm in
                ((a, finite_operator(2), norms[0]),
                 (b, finite_operator(2), norms[1]),
                 (product, finite_operator(4), norms[2]))]
    wrong_energy = s.simplify((wrong.H*finite_operator(4)*wrong)[0]/norms[2])
    if norms != [6, 6, 36] or energies[2] != energies[0]+energies[1]:
        raise ArithmeticError('disjoint fixed-number fermion cores failed additivity')
    if wrong_energy == energies[2]:
        raise ArithmeticError('ungraded tensor control did not change the energy')
    for mask in range(256):
        pairs = sum((mask >> (2*x)) & 3 == 3 for x in range(4))
        if pairs > mask.bit_count()//2:
            raise ArithmeticError('onsite pairs exceed half the particle number')
    return {'finite_sites': 4, 'full_Fock_dimension': 256,
            'source': 'reflected nearest-neighbor sigma_x Weyl hopping c=2, number shift 4, attraction 10',
            'core_particle_numbers': [2, 2], 'norms_squared': [str(x) for x in norms],
            'core_energies': [str(x) for x in energies[:2]],
            'graded_product_energy': str(energies[2]),
            'ungraded_product_energy': str(wrong_energy),
            'inter_core_hopping_expectation': '0',
            'onsite_pair_bound_exhausted_occupation_states': 256,
            'critical_family_bound_state_asserted_by_finite_control': False}


def certificate():
    return {'negative_pair': negative_packet_certificate(F(1, 16384), 1),
            'fixed_N_lower_bound': '-U_a floor(N/2)',
            'fixed_N_upper_bound': '6N/a^2',
            'finite_spacing_stability_constant': 'U_a/2 per particle',
            'fixed_N_self_adjoint': True,
            'Fock_direct_sum_self_adjoint': True,
            'Fock_infimum_for_a_gamma_at_most_1_over_4096': '-infinity',
            'empty_Fock_vacuum_is_ground_state': False,
            'necessary_finite_number_shift_limit_for_empty_ground': 'mu>=gamma^2/2',
            'necessity_proves_sufficiency': False,
            'uniform_continuum_many_body_stability_proved': False,
            'different_finite_density_representations_excluded': False,
            'Lorentz_shell': affine_shell_certificate(),
            'finite_fermion_product_control': finite_additivity_certificate()}
