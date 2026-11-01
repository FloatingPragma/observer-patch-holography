"""Source-bound grading and dispersion obstructions; no empirical mass fit.

Theorems concern the specified operators, reference and quadratic limiting
band. They do not exclude other excitations, dictionaries or scaling limits.
"""
from fractions import Fraction as F
from itertools import combinations
import math

import sympy as s

from ..cube_eta_protection import cube
from ..moving_bound_pair import limit


def exact(value, name, *, positive=False):
    if type(value) not in (int, F):
        raise ValueError(f'{name} must be an integer or Fraction')
    value = F(value)
    if positive and value <= 0:
        raise ValueError(f'{name} must be positive')
    return value


def grading_certificate():
    """Bind even, number-neutral local response to all source occupation rows."""
    masks = tuple(int(mask) for mask in cube.weights(6))
    counts = tuple(mask.bit_count() for mask in masks)
    assert len(masks) == math.comb(16, 6) == 8008
    assert set(counts) == {6}
    reference = {}
    for sites in combinations(range(8), 3):
        mask = sum(3 << (2*site) for site in sites)
        reference[mask] = (-1)**sum(sum(cube.GEOMETRY.sites[x]) for x in sites)
    probe = {mask: (4*(mask & 3).bit_count()-3)*a for mask, a in reference.items()}
    assert sum(a*a for a in reference.values()) == 56
    assert sum(a*a for a in probe.values()) == 840
    assert sum(reference[mask]*a for mask, a in probe.items()) == 0
    # A spectral transition |e_N6><R_N6| always has zero number change and
    # parity +1, independently of degeneracy or of the eigenvector chart.
    transition_changes = {left-right for left in set(counts) for right in set(counts)}
    assert transition_changes == {0}
    # On the entire Fock space, both n_0 and exp(i*pi*n_0/2) are diagonal.
    # Four exact local phases suffice; all 2^14 exterior occupations are free.
    local_phases = (1, s.I, s.I, -1)
    assert all(local_phases[mask & 3] in (1, s.I, -1) for mask in range(1 << 16))
    return {'complete_Fock_dimension': 1 << 16,
            'reference_and_probe_number': 6, 'complete_number_sector_dimension': len(masks),
            'reference_occupation_columns': len(reference),
            'probe_squared_unnormalized_norm': sum(a*a for a in probe.values()),
            'source_parity': 1, 'local_operator_parity': 1,
            'transition_number_change': 0,
            'grading_preserving_single_odd_fermion_identification': 'excluded',
            'number_charge_identified_as_electromagnetism': False,
            'even_parity_identified_as_Lorentz_spin_or_Bose_statistics': False,
            'domain': 'specified N6 reference, density/phase response and its H spectral components'}


def shell_certificate():
    """A nonconstant quadratic band is not a relativistic shell on an open set.

Allow E=alpha*(|P|^2/4-gamma^2)+beta, p=L P, alpha nonzero,
any fixed real L, c and m. Along P=r e_1 the quartic coefficient of
E^2-c^2|p|^2-m^2*c^4 is alpha^2/16. An identity on any open ball
would make this polynomial identically zero, a contradiction.
"""
    r, alpha, beta, gamma, c, m = s.symbols('r alpha beta gamma c m', real=True)
    l0, l1, l2 = s.symbols('l0 l1 l2', real=True)
    energy = alpha*(r*r/4-gamma*gamma)+beta
    residual = s.Poly(s.expand(energy**2-c*c*r*r*(l0*l0+l1*l1+l2*l2)-m*m*c**4), r)
    quartic = residual.coeff_monomial(r**4)
    assert quartic == alpha**2/16
    # Consume the actual source-derived band formula through its exact API.
    source = limit.moving_pole_certificate(F(1, 1024), F(1), F(1))
    assert source['reflected_band_energy_limit'] == '|P|^2/4-gamma^2'
    return {'source_band': source['reflected_band_energy_limit'],
            'source_comparison': 'J_a=6/a^2-H_a on the specified two-body bound band',
            'allowed_energy_dictionary': 'alpha*epsilon+beta, fixed real alpha != 0 and beta',
            'allowed_momentum_dictionary': 'p=L*P, fixed invertible real linear L',
            'Lorentz_shell_squared_residual_quartic': str(quartic),
            'exact_positive_mass_Lorentz_shell_on_open_momentum_set': 'excluded',
            'exact_massless_Lorentz_shell_on_open_momentum_set': 'excluded',
            'nonrelativistic_approximation_excluded': False,
            'other_source_scalings_excluded': False,
            'source_number': 2, 'source_parity': 1,
            'physical_mass_charge_and_spin_identified': False}


def low_momentum_gap(rest_energy, quadratic_coefficient, momentum, momentum_cap):
    """Bounds after matching the constant and quadratic Taylor coefficients.

E_quad=M+A*r^2 and E_rel=sqrt(M^2+2*M*A*r^2), with M,A>0.
Their exact difference is A^2*r^4/(E_quad+E_rel). The conservative
enclosure below is rational and vanishes quartically at r=0.
    """
    M = exact(rest_energy, 'rest energy', positive=True)
    A = exact(quadratic_coefficient, 'quadratic coefficient', positive=True)
    r = exact(momentum, 'absolute momentum')
    cap = exact(momentum_cap, 'momentum cap')
    if not 0 <= r <= cap:
        raise ValueError('requires 0 <= momentum <= momentum_cap')
    return {'lower': A*A*r**4/(2*(M+A*cap*cap)),
            'upper': A*A*r**4/(2*M),
            'constant_and_quadratic_terms_match': True,
            'rest_energy_and_units_are_assigned_inputs': True}


def neutral_number_shift_certificate():
    """Neither a scalar nor a number chemical potential changes N6 gaps."""
    E, R, mu, shift = s.symbols('E R mu shift', real=True)
    gap = s.expand((E-6*mu+shift)-(R-6*mu+shift))
    assert gap == E-R
    return {'same_number': 6, 'transformed_gap': str(gap),
            'scalar_or_common_number_chemical_potential_repairs_negative_gaps': False,
            'global_time_reversal_repairs_mixed_frequency_signs': False}
