"""Energy tests for the configured coherent phase and neutral attenuation.

All coefficients are read from the native coherent analyzer. No measured
mass, fitted phase, fitted attenuation or SI conversion enters this module.
Minimization over controller settings is an external variational hypothesis.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
import json

import sympy as s

from model.native_flavor_orbit import intervals as iv
from model.standard_model_source_map import koide_source
from model.particle_candidate_atlas.interacting_chiral import quantum_higgs_operators


SOURCE_INPUTS = (
    'model/standard_model_source_map/koide_source.py',
    'model/standard_model_source_map/KOIDE_SOURCE.md',
    'model/particle_candidate_atlas/attenuation.py',
    'model/particle_candidate_atlas/BARE_EFFECTIVE_CONTROLLER.md',
    'model/charged_lepton_response/relativistic_mass.py',
    'model/charged_lepton_response/RELATIVISTIC_MASS.md',
    'model/particle_candidate_atlas/interacting_chiral.py',
    'model/particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md',
    'model/interacting_extension/higgs_join/interacting_carrier.py',
    'model/charged_lepton_response/SOURCE_PHASE_SELECTION.md',
    'model/particle_source_selection/PROOF.md',
    'model/particle_candidate_atlas/CAPACITY_NON_SELECTION.md',
    'Lean/ObserverPatchHolography/KoideCirculant.lean',
)


def _rational(value):
    if type(value) not in (int, str, F):
        raise ValueError('an exact rational input is required')
    value = F(value)
    return s.Rational(value.numerator, value.denominator)


def _clean(value):
    return s.factor(s.radsimp(s.expand(value)))


def source_coefficients():
    weights = tuple(koide_source.WEIGHTS)
    if weights != tuple(map(s.Rational, ('10/72', '5/72', '5/72', '2/72', '50/72'))):
        raise ArithmeticError('the native coherent branch weights differ')
    return s.sqrt(weights[0]), s.sqrt(weights[1])


def roots_and_masses(real=1, imag=0, attenuation='1'):
    """Exact C eigenvalues and Y eigenvalues at a declared unit phase.

    The mass normalization is |v_H|=1. A common external energy conversion
    is absent. The formula also applies outside the positive-root chamber,
    where signed-root Koide and positive-root Koide differ.
    """
    x, y = map(s.sympify, (real, imag))
    if _clean(x*x+y*y-1) != 0:
        raise ValueError('the coherent phase must have unit modulus')
    r = _rational(attenuation)
    if not 0 <= r <= 1:
        raise ValueError('attenuation must lie in [0,1]')
    a, b = source_coefficients()
    roots = tuple(_clean(a+r*b*t) for t in (2*x, -x-s.sqrt(3)*y, -x+s.sqrt(3)*y))
    return roots, tuple(_clean(c*c/3) for c in roots)


def moment_certificate():
    """Laurent-polynomial identities on the literal native three-cycle."""
    a, b = source_coefficients()
    r = s.Symbol('r', real=True)
    z = s.Symbol('z', nonzero=True)
    rot = koide_source.P
    c = a*s.eye(3)+r*b*(z*rot+rot**2/z)
    y = c*c/3
    tr_y = 5*(1+r*r)/36
    tr_y2 = (s.Rational(25, 3888)+s.Rational(25, 648)*r*r
             +s.Rational(25, 2592)*r**4
             +s.Rational(25, 1944)*s.sqrt(2)*r**3*(z**3+z**-3)/2)
    if _clean(s.trace(y)-tr_y) != 0 or _clean(s.trace(y*y)-tr_y2) != 0:
        raise ArithmeticError('the source spectral moments differ')
    return {'trace_Y': '5(1+r^2)/36',
            'trace_Y_squared': '25/3888 + 25r^2/648 + 25r^4/2592 + 25sqrt(2)r^3 cos(3theta)/1944',
            'r_one_trace_Y_squared': '(425+100sqrt(2)cos(3theta))/7776',
            'literal_three_cycle_replayed': True,
            'phase_flat_free_zero_momentum_energy': '-2|v_H| Tr(Y)',
            'normalized_coherent_Koide_Q': '(1+r^2)/3'}


def convexity_certificate():
    """The exact single-root derivative used in the attenuation proof."""
    x = s.Symbol('x', real=True)
    p, h = s.symbols('p h', positive=True)
    f = s.sqrt(p*p+h*h*x**4/9)
    expected = (2*h*h*x*x/9)*(3*p*p+h*h*x**4/9)/(p*p+h*h*x**4/9)**s.Rational(3, 2)
    if s.simplify(s.diff(f, x, 2)-expected) != 0:
        raise ArithmeticError('the convexity identity differs')
    return {'positive_momentum_second_derivative':
            '2h^2 x^2 (3p^2+h^2x^4/9)/(9(p^2+h^2x^4/9)^(3/2))',
            'zero_momentum_function': '|h| x^2/3',
            'root_slopes_sum': '0',
            'free_energy_selected_attenuation_under_external_minimization': '1',
            'applies_to_every_positive_weighted_finite_momentum_set': True,
            'native_attenuation_controller_selected': False}


def _sqrt_box(value, digits=48):
    value = iv.box(value)
    if value.lo < 0:
        raise ValueError('nonnegative square-root interval required')
    return iv.I(iv.sqrt_rational(value.lo, digits).lo,
                iv.sqrt_rational(value.hi, digits).hi)


def free_energy_interval(momentum, real=1, imag=0, attenuation='1'):
    """Charged free-CAR sea energy per momentum, with both spin states.

    Phase-independent massless spectators are retained by the source field
    and cancel in the energy comparisons. Momentum denotes c|p| for the
    continuum symbol, or the same lattice-symbol kinetic norm.
    """
    p = F(momentum)
    if p < 0:
        raise ValueError('nonnegative kinetic norm required')
    _, masses = roots_and_masses(real, imag, attenuation)
    energy = -2*sum((_sqrt_box(p*p+iv.expression(m)**2) for m in masses), iv.I(0))
    return iv.rounded_out(energy, 30)


def finite_momentum_comparisons():
    """Rational enclosures at declared phases and momenta, without mass data."""
    phase = (s.Rational(63, 65), s.Rational(16, 65))
    endpoint = ((s.sqrt(6)+s.sqrt(2))/4, (s.sqrt(6)-s.sqrt(2))/4)
    rows = []
    for p in (F(1, 256), F(1, 64)):
        energies = {'zero': free_energy_interval(p),
                    'rational_control': free_energy_interval(p, *phase),
                    'positive_chamber_endpoint': free_energy_interval(p, *endpoint)}
        delta = energies['rational_control']-energies['zero']
        endpoint_delta = energies['rational_control']-energies['positive_chamber_endpoint']
        if p == F(1, 256) and not (delta.hi < 0 and endpoint_delta.hi < 0):
            raise ArithmeticError('interior-phase energy control failed')
        if p == F(1, 64) and not delta.lo > 0:
            raise ArithmeticError('momentum-dependent ordering control failed')
        rows.append({'kinetic_norm': str(p),
                     'energies': {k: v.strings() for k, v in energies.items()},
                     'control_minus_zero': delta.strings(),
                     'control_minus_endpoint': endpoint_delta.strings(),
                     'interior_minimizer_exists_in_positive_chamber': p == F(1, 256)})
    return {'rational_phase': '(63+16i)/65', 'attenuation': '1',
            'common_neutral_Higgs_norm': '1', 'spin_degeneracy': 2,
            'rows': rows,
            'global_minimizing_phase_value_computed': False,
            'source_selects_either_momentum_sample': False,
            'measured_inputs_used': []}


def _create(mask, mode):
    if mask >> mode & 1:
        return None
    return mask | (1 << mode), (-1)**((mask & ((1 << mode)-1)).bit_count())


def _annihilate(mask, mode):
    if not mask >> mode & 1:
        return None
    return mask ^ (1 << mode), (-1)**((mask & ((1 << mode)-1)).bit_count())


def _pair_terms(families):
    """Literal charged-Higgs spin singlet, including its CAR order."""
    for h in range(2):
        for f in range(families):
            for g in range(families):
                for left_spin, right_spin, sign in ((0, 1, 1), (1, 0, -1)):
                    yield h, f, g, 6*f+2*h+left_spin, 6*g+4+right_spin, sign


def vertex_action(y, scalar, mask):
    """Full finite-Higgs vertex on any occupation input; no sector cutoff."""
    y = s.ImmutableMatrix(y)
    families = y.rows
    if y.cols != families or not 1 <= families <= 3:
        raise ValueError('one to three family coefficients required')
    if type(scalar) is not int or not 0 <= scalar <= 2 or type(mask) is not int or not 0 <= mask < 1 << (6*families):
        raise ValueError('valid scalar and full CAR occupation required')
    result = {}
    for h, f, g, left, right, sign in _pair_terms(families):
        if scalar == 0:
            # c_left^dagger c_right^dagger acts from right to left.
            first = _create(mask, right)
            second = None if first is None else _create(first[0], left)
            coefficient, destination = y[f, g], h+1
        elif scalar == h+1:
            first = _annihilate(mask, left)
            second = None if first is None else _annihilate(first[0], right)
            coefficient, destination = s.conjugate(y[f, g]), 0
        else:
            continue
        if first is not None and second is not None and coefficient != 0:
            key = destination, second[0]
            result[key] = result.get(key, 0)+sign*first[1]*second[1]*coefficient
    return {k: _clean(v) for k, v in result.items() if _clean(v) != 0}


def zero_charge_basis(families=3):
    if type(families) is not int or not 1 <= families <= 3:
        raise ValueError('one to three families required')
    pairs = tuple((1 << i)+(1 << j) for i, j in combinations(range(6*families), 2))
    return ((0, 0), *((h, mask) for h in (1, 2) for mask in pairs))


def zero_charge_block(y, pair_cost=0):
    """The whole N_F-2N_H=0 block, including all dark occupation states.

    The optional pair cost is 2epsilon_F+epsilon_H for common number terms.
    Spatial hopping, link terms and scalar transport are excluded here.
    """
    y = s.ImmutableMatrix(y)
    basis = zero_charge_basis(y.rows)
    indices = {v: i for i, v in enumerate(basis)}
    cost = s.sympify(pair_cost)
    matrix = s.MutableSparseMatrix(len(basis), len(basis), {})
    for col, state in enumerate(basis):
        for destination, coefficient in vertex_action(y, *state).items():
            if destination not in indices:
                raise ArithmeticError('the claimed conserved block is not invariant')
            matrix[indices[destination], col] = coefficient
        if col:
            matrix[col, col] += cost
    return s.ImmutableSparseMatrix(matrix), basis


def quantum_block_certificate():
    """Compare the actual 192-state vertex, then certify the 307-state block."""
    actual = quantum_higgs_operators()
    for scalar in range(3):
        for mask in range(64):
            got = vertex_action(s.eye(1), scalar, mask)
            expected = {divmod(row, 64): coefficient
                        for (row, _), coefficient in actual['V'][:, scalar*64+mask].todok().items()}
            if got != expected:
                raise ArithmeticError('the literal source quantum-Higgs tensor differs')
    # A complex off-diagonal control checks every family-pair and adjoint sign.
    y = s.ImmutableMatrix([[2, s.I, 1], [-s.I, 3, 2*s.I], [1, -2*s.I, 4]])
    block, basis = zero_charge_block(y, 7)
    if block != block.H:
        raise ArithmeticError('the complete conserved block is not Hermitian')
    bright = block[:, 0]
    expected_norm2 = 4*s.trace(y.H*y)
    if (bright.H*bright)[0] != expected_norm2:
        raise ArithmeticError('the bright vector norm differs')
    vacuum = s.ImmutableSparseMatrix(len(basis), 1, {(0, 0): 1})
    if block*bright != expected_norm2*vacuum+7*bright:
        raise ArithmeticError('the vacuum-bright span is not invariant')
    if block[1:, 1:] != 7*s.eye(len(basis)-1):
        raise ArithmeticError('the complementary dark states were changed')
    return {'literal_source_full_CAR_states_compared': 192,
            'three_family_active_CAR_modes': 18,
            'full_local_dimension': 3*(1 << 18),
            'conserved_operator': 'N_F-2N_H',
            'complete_zero_charge_block_dimension': len(basis),
            'vacuum_bright_invariant_dimension': 2,
            'retained_dark_dimension': len(basis)-2,
            'bright_squared_norm': '4 Tr(Y^dagger Y)',
            'two_bright_energies': '(d +/- sqrt(d^2+16 Tr(Y^dagger Y)))/2',
            'dark_energy': 'd',
            'd_definition': '2epsilon_F+epsilon_H',
            'whole_sector_ground_is_lower_bright_energy': True,
            'other_conserved_sectors_retained_by_full_vertex': True,
            'full_local_ground_identified': False,
            'spatial_interacting_ground_identified': False,
            'external_minimization_on_positive_root_chamber': {'theta': '0', 'r': '1'},
            'native_controller_or_ground_preparation_selected': False}


def certificate():
    return {'source_inputs': list(SOURCE_INPUTS),
            'native_coefficients': list(map(str, source_coefficients())),
            'moments': moment_certificate(),
            'free_attenuation': convexity_certificate(),
            'finite_momentum_controls': finite_momentum_comparisons(),
            'finite_quantum_Higgs': quantum_block_certificate(),
            'source_selected_phase': None,
            'source_selected_attenuation': None,
            'measured_inputs_used': [],
            'possible_construction_routes': [
            'Compute the effective phase and attenuation from the action and complete correlations of a fixed declared program.',
            'Specify a dynamical phase/attenuation register, its coupling and mixing or preparation law, and derive its complete outputs on one energy reference.'],
            'architecture_uniqueness_required': False}


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2, sort_keys=True))
