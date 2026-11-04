"""Complete supplied-mass bulk branches and punctured hedgehog level counts."""
from functools import lru_cache
from itertools import product

import sympy as s

from ..native_mass_channels import algebra as source
from .algebra import generators


@lru_cache(maxsize=1)
def dirac_blocks():
    gamma, masses = source.folded_generators()
    scalar, spin = masses[0], masses[3]
    charge = scalar*spin
    rows = []
    for sign in (-1, 1):
        for by in (0, 1):
            indices = [2*source.SITES.index((0, by, bz))+sigma
                       for bz in (0, 1) for sigma in (0, 1)]
            inclusion = s.eye(16)[:, indices]
            block = s.ImmutableMatrix((s.eye(16)+sign*charge)*inclusion/s.sqrt(2))
            if block.H*block != s.eye(4) or charge*block != sign*block:
                raise ArithmeticError('mass-sector Dirac inclusion failed')
            small = tuple(s.ImmutableMatrix(block.H*g*block) for g in gamma+(spin,))
            if any(g*block != block*h for g, h in zip(gamma+(spin,), small, strict=True)):
                raise ArithmeticError('four-row source block is not invariant')
            words = []
            for mask in range(16):
                word = s.eye(4)
                for j, g in enumerate(small):
                    if (mask >> j) & 1:
                        word *= g
                words.append(s.Matrix(list(word)))
            if s.Matrix.hstack(*words).rank() != 16:
                raise ArithmeticError('source Dirac block is reducible')
            rows.append((sign, by, block))
    full = s.ImmutableMatrix.hstack(*(row[2] for row in rows))
    if full.H*full != s.eye(16) or full*full.H != s.eye(16):
        raise ArithmeticError('bulk decomposition omitted source rows')
    return tuple(rows)


def spin_rotation(direction):
    values = source.vector(direction, 3)
    norm = s.sqrt(sum(s.Rational(x)**2 for x in values))
    if norm == 0:
        raise ValueError('nonzero exact spin-mass direction required')
    beta = generators()[1]
    unit = tuple(s.Rational(x)/norm for x in values)
    target = sum((x*g for x, g in zip(unit, beta, strict=True)), s.zeros(16))
    if unit[2] == -1:
        result = beta[0]*beta[2]
    else:
        result = (s.eye(16)+target*beta[2])/s.sqrt(2*(1+unit[2]))
    return s.ImmutableMatrix(result.applyfunc(s.simplify))


def coulomb_level(principal):
    if type(principal) is not int or principal < 1:
        raise ValueError('positive integer Coulomb principal quantum number required')
    return {'principal_quantum_number': principal,
            'energy_squared_over_mass_squared': str(s.Rational(principal**2-1, principal**2)),
            'zero_multiplicity': 2 if principal == 1 else 0,
            'each_signed_energy_multiplicity': 0 if principal == 1 else 2*principal**2,
            'orbital_degeneracy_per_internal_copy': principal**2,
            'internal_copies': 2}


def angular_potential_certificate():
    _, _, _, b, _ = generators()
    number = sum((x.H*x for x in b), s.zeros(16))
    # Rational unit direction gives exact cross terms, not an axis-only check.
    direction = (s.Rational(2, 7), s.Rational(3, 7), s.Rational(6, 7))
    radial = sum((x*g for x, g in zip(direction, b, strict=True)), s.zeros(16))
    radial_number = (radial.H*radial).applyfunc(s.expand)
    zero = s.zeros(16)
    if ((radial_number*radial_number-radial_number).applyfunc(s.expand) != zero or
            (number*radial_number-radial_number*number).applyfunc(s.expand) != zero):
        raise ArithmeticError('radial fermion occupation is not a commuting projection')
    potential = s.ImmutableMatrix(number-s.eye(16)-radial_number)
    form_rows = []
    for form in range(4):
        projector = s.eye(16)
        for other in range(4):
            if other != form:
                projector = projector*(number-other*s.eye(16))/s.Rational(form-other)
        projector = s.ImmutableMatrix(projector.applyfunc(s.expand))
        rank = 2*s.binomial(3, form)
        if (projector*projector-projector).applyfunc(s.expand) != zero or s.trace(projector) != rank:
            raise ArithmeticError('an exterior-form source sector was omitted')
        if form == 0 and (potential*projector+projector).applyfunc(s.expand) != zero:
            raise ArithmeticError('Coulomb attractive scalar sign failed')
        if form == 3 and (potential*projector-projector).applyfunc(s.expand) != zero:
            raise ArithmeticError('Coulomb repulsive top-form sign failed')
        if form == 1 and ((potential+radial_number)*projector).applyfunc(s.expand) != zero:
            raise ArithmeticError('one-form radial attraction sign failed')
        if form == 2 and ((potential-s.eye(16)+radial_number)*projector).applyfunc(s.expand) != zero:
            raise ArithmeticError('two-form positivity identity failed')
        form_rows.append({'form_degree': form, 'rank': int(rank),
                          'angular_potential_in_units_2vm_over_r':
                              ('-1', '-n_radial', '1-n_radial', '+1')[form]})
    return {'full_source_dimension': 16, 'form_sectors': form_rows,
            'gradient_Hessian_term': '(m/r)*(N_f-n_radial)',
            'squared_H_angular_term': '(2*v*m/r)*(N_f-1-n_radial)',
            'two_and_three_form_subthreshold_states': 0,
            'punctured_constant_magnitude_core_profile_supplied': True,
            'smooth_core_excited_levels_identified_with_Coulomb_levels': False}


def certificate():
    rows = dirac_blocks()
    return {'bulk': {'source_dimension': 16, 'all_weyl_valleys': 8,
                     'four_row_irreducible_blocks': len(rows),
                     'blocks_per_scalar_spin_mass_sector': 2,
                     'mass_branches': ['abs(m0+|m|)', 'abs(m0-|m|)'],
                     'pure_spin_equal_mass_Dirac_blocks': 4,
                     'at_m0_equal_abs_m': {'massless_Dirac_blocks': 2, 'mass_2abs_m_Dirac_blocks': 2},
                     'mass_coefficients_supplied': True, 'physical_species_identified': False},
            'Coulomb_angular': angular_potential_certificate(),
            'Coulomb_levels': [coulomb_level(n) for n in range(1, 8)],
            'Coulomb_general_signed_multiplicity': '2*n^2 for integer n>=2',
            'Coulomb_zero_multiplicity': 2,
            'Coulomb_essential_spectrum': '(-infinity,-m] union [m,infinity)',
            'Coulomb_excited_eigenvalues': '+/-m*sqrt(1-1/n^2), n>=2',
            'generic_smooth_core_all_excited_levels_solved': False}
