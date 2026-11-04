"""Exact full-valley Clifford and oscillator algebra for supplied textures."""
from functools import lru_cache
from itertools import product
from math import comb

import sympy as s

from ..native_mass_channels import algebra as source


@lru_cache(maxsize=1)
def generators():
    gamma, masses = source.folded_generators()
    beta = masses[1:]
    rho = tuple(s.ImmutableMatrix(s.I*g*b) for g, b in zip(gamma, beta, strict=True))
    annihilators = tuple(s.ImmutableMatrix((b+s.I*g)/2) for g, b in zip(gamma, beta, strict=True))
    return gamma, beta, rho, annihilators, s.ImmutableMatrix(rho[0]*rho[1]*rho[2])


def zero_basis(orientation=1):
    if type(orientation) is not int or orientation not in (-1, 1):
        raise ValueError('texture orientation +1 or -1 required')
    result = s.zeros(16, 2)
    for address, valley in enumerate(source.SITES):
        value = (orientation*s.I)**sum(valley)/s.sqrt(8)
        result[2*address, 0] = value
        result[2*address+1, 1] = value
    return s.ImmutableMatrix(result)


def occupation_projector(occupation):
    if (type(occupation) is not tuple or len(occupation) != 3 or
            any(type(x) is not int or x not in (0, 1) for x in occupation)):
        raise ValueError('three binary occupation labels required')
    result = s.eye(16)
    for rho, value in zip(generators()[2], occupation, strict=True):
        result = result*(s.eye(16)+(1-2*value)*rho)/2
    return s.ImmutableMatrix(result)


@lru_cache(maxsize=1)
def internal_fock_basis():
    b = generators()[3]
    columns = []
    for occupation in product((0, 1), repeat=3):
        pair = zero_basis()
        for j in reversed(range(3)):
            if occupation[j]:
                pair = b[j].H*pair
        columns.extend(pair[:, j] for j in range(2))
    return s.ImmutableMatrix.hstack(*columns)


def clifford_certificate():
    gamma, beta, rho, b, chirality = generators()
    eye, zero = s.eye(16), s.zeros(16)
    for i in range(3):
        if rho[i] != rho[i].H or rho[i]*rho[i] != eye:
            raise ArithmeticError('source internal grading changed')
        for j in range(3):
            if (b[i]*b[j]+b[j]*b[i] != zero or
                    b[i]*b[j].H+b[j].H*b[i] != int(i == j)*eye):
                raise ArithmeticError('source Clifford CAR failed')
            if rho[i]*rho[j] != rho[j]*rho[i]:
                raise ArithmeticError('source internal gradings do not commute')
    expected = s.I*gamma[0]*gamma[1]*gamma[2]*beta[0]*beta[1]*beta[2]
    if chirality != expected or chirality*chirality != eye or chirality != chirality.H:
        raise ArithmeticError('source chiral orientation changed')
    if any(chirality*g+g*chirality != zero for g in gamma+beta):
        raise ArithmeticError('texture operator lost its chiral grading')
    rows = []
    for occupation in product((0, 1), repeat=3):
        projector = occupation_projector(occupation)
        if projector != projector.H or projector*projector != projector or s.trace(projector) != 2:
            raise ArithmeticError('an internal occupation sector was omitted')
        rows.append({'occupation': list(occupation), 'rank': 2, 'chirality': (-1)**sum(occupation)})
    change = internal_fock_basis()
    if change.H*change != eye or change*change.H != eye:
        raise ArithmeticError('source-to-internal-Fock map is not unitary')
    for orientation in (-1, 1):
        basis = zero_basis(orientation)
        maps = b if orientation == 1 else tuple(x.H for x in b)
        if (basis.H*basis != s.eye(2) or any(x*basis != s.zeros(16, 2) for x in maps)
                or chirality*basis != orientation*basis):
            raise ArithmeticError('normalized defect polarization columns changed')
    return {'source_folded_dimension': 16, 'weyl_valleys_retained': 8,
            'internal_CAR_modes': 3, 'internal_Fock_multiplicity': 2,
            'occupation_sectors': rows, 'full_internal_Fock_change_of_basis_unitary': True,
            'positive_texture_zero_polarization_rank': 2, 'negative_texture_zero_polarization_rank': 2,
            'positive_texture_chiral_index': 2, 'negative_texture_chiral_index': -2,
            'source_selected_order_parameter_proved': False}


def radial_hessian_certificate():
    x = s.Matrix(s.symbols('x0:3', real=True))
    mass, length = s.symbols('m ell', positive=True)
    radius = s.sqrt(x.dot(x)+length**2)
    potential = mass*(radius-length)
    hessian = s.hessian(potential, x)
    expected = mass/radius*(s.eye(3)-x*x.T/radius**2)
    if (hessian-expected).applyfunc(s.simplify) != s.zeros(3):
        raise ArithmeticError('radial Hessian differs from the potential derivative')
    if (hessian*x-mass*length**2/radius**3*x).applyfunc(s.simplify) != s.zeros(3, 1):
        raise ArithmeticError('radial Hessian eigenvalue changed')
    if s.simplify(hessian.det()-mass**3*length**2/radius**5) != 0:
        raise ArithmeticError('strict Hessian positivity determinant changed')
    return {'potential': 'm*(sqrt(|x|^2+ell^2)-ell), m>0, ell>0',
            'gradient': 'm*x/sqrt(|x|^2+ell^2)',
            'radial_Hessian_eigenvalue': 'm*ell^2/(|x|^2+ell^2)^(3/2)',
            'tangential_Hessian_eigenvalue': 'm/sqrt(|x|^2+ell^2)',
            'tangential_multiplicity': 2, 'strict_positive_Hessian_every_finite_point': True,
            'zero_mode_envelope': 'exp(-m*(sqrt(|x|^2+ell^2)-ell)/v)',
            'normalizable_zero_mode_dimension': 2,
            'other_subthreshold_discrete_levels_classified': False}


def oscillator_multiplicity(level):
    if type(level) is not int or level < 0:
        raise ValueError('nonnegative integer oscillator level required')
    reduced = sum(comb(3, f)*comb(level-f+2, 2) for f in range(min(3, level)+1))
    if reduced != (1 if level == 0 else 4*level*level+2):
        raise ArithmeticError('oscillator full-sector combinatorics failed')
    return {'level': level, 'energy_squared': f'2*v*kappa*{level}',
            'zero_multiplicity': 2 if level == 0 else 0,
            'each_signed_energy_multiplicity': 0 if level == 0 else reduced,
            'squared_energy_multiplicity': 2*reduced}


def oscillator_states(level):
    oscillator_multiplicity(level)
    rows = []
    for occupation in product((0, 1), repeat=3):
        total = level-sum(occupation)
        if total < 0:
            continue
        for nx in range(total+1):
            for ny in range(total-nx+1):
                rows.extend(((nx, ny, total-nx-ny), occupation, copy) for copy in (0, 1))
    return tuple(rows)


def oscillator_block(level):
    """H/sqrt(2v*kappa) on its complete invariant q-level space."""
    states = oscillator_states(level)
    index, entries = {row: j for j, row in enumerate(states)}, {}
    for column, (spatial, occupation, copy) in enumerate(states):
        for axis in range(3):
            target_n, target_f = list(spatial), list(occupation)
            sign = (-1)**sum(occupation[:axis])
            if occupation[axis]:
                target_n[axis] += 1
                target_f[axis] = 0
                coefficient = sign*s.sqrt(spatial[axis]+1)
            else:
                if spatial[axis] == 0:
                    continue
                target_n[axis] -= 1
                target_f[axis] = 1
                coefficient = sign*s.sqrt(spatial[axis])
            row = index[tuple(target_n), tuple(target_f), copy]
            entries[row, column] = entries.get((row, column), 0)+coefficient
    return s.ImmutableSparseMatrix(len(states), len(states), entries)


def certificate():
    return {'clifford': clifford_certificate(), 'smooth_radial_texture': radial_hessian_certificate(),
            'linear_texture_levels': [oscillator_multiplicity(n) for n in range(9)],
            'linear_texture_general_signed_multiplicity': '4*q^2+2 for integer q>=1',
            'linear_texture_zero_multiplicity': 2,
            'linear_texture_is_unbounded_at_spatial_infinity': True}
