"""Exact finite Wilson decomposition and volume-uniform local spectral bounds.

The spectral projector is that of the uncompressed coordinate symbol.
Its expectation is a complete finite-code effect. It is not a particle
mass projector, a locality estimate or a global kernel-gap probability.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import itertools
import json
from pathlib import Path

import sympy as s

from model.interacting_extension.chiral_lattice import source_data
from model.particle_local_gauge_state import state
from model.particle_observer_outputs import field_route


ROOT = Path(__file__).resolve().parents[2]
SOURCE_INPUTS = (
    'model/particle_local_gauge_state/state.py',
    'model/particle_local_gauge_state/PROOF.md',
    'model/particle_field_preparation/preparation.py',
    'model/particle_field_preparation/REFINEMENT.md',
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/interacting_extension/chiral_lattice/source_data.py',
    'model/interacting_extension/chiral_lattice/PROOF.md',
    'model/common_source/overlap_regulator.py',
    'model/observer_particles/results/representation_identification_operators.json',
)


def integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(name+' must be an integer at least '+str(minimum))
    return value


def rational(value, name):
    if type(value) not in (int, str, F):
        raise ValueError(name+' requires exact rational data')
    value = F(value)
    if value < 0:
        raise ValueError(name+' must be nonnegative')
    return value


def clean(matrix):
    return s.ImmutableSparseMatrix(matrix).applyfunc(s.simplify)


def comm(a, b):
    return a*b-b*a


@lru_cache(maxsize=1)
def gamma_data():
    native = source_data.load()
    raw, _ = native['gamma']
    gamma = tuple(clean(s.Matrix(g).applyfunc(lambda x: s.nsimplify(x, rational=True)))
                  for g in raw)
    eye = s.eye(4)
    for i, a in enumerate(gamma):
        if a != a.H:
            raise ArithmeticError('the source gamma is not Hermitian')
        for j, b in enumerate(gamma):
            if clean(a*b+b*a) != (2*eye if i == j else s.zeros(4)):
                raise ArithmeticError('the source Clifford relations differ')
    return gamma


def native_representation(fundamental):
    """Actual native even exterior basis Lambda²(C5) plus Lambda⁴(C5).

    This helper uses the frozen source basis and verifies its q6Y weights;
    it does not infer the representation from its dimension.
    """
    g = clean(fundamental)
    if g.shape != (5, 5) or clean(g.H*g) != s.eye(5) or s.simplify(g.det()) != 1:
        raise ValueError('a determinant-one unitary fundamental matrix is required')
    if any(g[i, j] != 0 for i in range(5) for j in range(5) if (i < 3) != (j < 3)):
        raise ValueError('the source fundamental must preserve its color/weak blocks')
    raw = json.loads((ROOT/SOURCE_INPUTS[-1]).read_text())
    basis = tuple(tuple(raw['basis_subsets'][i]) for i in raw['even_indices'])
    expected = tuple(sum(-2 if i < 3 else 3 for i in b) for b in basis)
    if expected != source_data.load()['charges']:
        raise ArithmeticError('the frozen native exterior basis has different weights')
    out = s.zeros(15)
    for i, left in enumerate(basis):
        for j, right in enumerate(basis):
            if len(left) == len(right):
                out[i, j] = g.extract(left, right).det()
    out = clean(out)
    if clean(out.H*out) != s.eye(15):
        raise ArithmeticError('the native exterior lift is not unitary')
    return out


def native_weak_link(pauli_index, cosine, sine):
    c, r = s.sympify(cosine), s.sympify(sine)
    if s.simplify(c*c+r*r) != 1:
        raise ValueError('exact unit-circle coefficients are required')
    pauli = (s.Matrix([[0, 1], [1, 0]]),
             s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    index = integer(pauli_index, 'Pauli index')
    if index > 2:
        raise ValueError('Pauli index must be zero, one or two')
    return native_representation(s.diag(s.eye(3), c*s.eye(2)+s.I*r*pauli[index]))


def decomposition(shifts):
    """The complete HJL identity for three exact finite unitary shifts.

    Each shift acts on spatial sites and the native internal representation;
    spin is tensored in without dropping any Clifford or internal columns.
    """
    if len(shifts) != 3:
        raise ValueError('three covariant shifts are required')
    shifts = tuple(clean(t) for t in shifts)
    size = shifts[0].rows
    eye = s.eye(size)
    if not size or any(t.shape != (size, size) or clean(t.H*t) != eye for t in shifts):
        raise ValueError('all shifts must be unitary matrices of one common dimension')
    gamma = gamma_data()
    D = tuple(t-eye for t in shifts)
    Db = tuple(eye-t.H for t in shifts)
    # A=I-aD_W at Wilson parameter s=0. The spatial source H_W=-beta A.
    A = s.eye(4*size)
    for mu in range(3):
        A += s.kronecker_product((s.eye(4)-gamma[mu])/2, D[mu])
        A -= s.kronecker_product((s.eye(4)+gamma[mu])/2, Db[mu])
    A = clean(A)
    kernel = clean(-s.kronecker_product(gamma[3], eye)*A)
    positive = s.zeros(4*size)
    error = s.zeros(4*size)
    ordered = []
    for mu, nu in itertools.permutations(range(3), 2):
        pair = D[nu]*D[mu]
        P = clean(pair.H*pair)
        R = clean(-Db[mu]*comm(D[mu], Db[nu]-D[nu]))
        B = clean(Db[mu]*D[mu]*Db[nu]*D[nu])
        if clean(B-P-R) != s.zeros(size):
            raise ArithmeticError('positive Wilson pair decomposition failed')
        sigma = s.I*comm(gamma[mu], gamma[nu])/2
        C = s.kronecker_product(s.I*sigma/2, comm(Db[mu]+D[mu], Db[nu]+D[nu]))
        Q = -s.kronecker_product(gamma[mu], comm(Db[mu]+D[mu], Db[nu]-D[nu]))
        positive += s.kronecker_product(s.eye(4), P)/4
        error += (s.kronecker_product(s.eye(4), R)+C+Q)/4
        ordered.append({'mu': mu, 'nu': nu,
                        'B_error_commutator_coefficient': 4,
                        'C_commutator_coefficient': 2,
                        'D_commutator_coefficient': 4})
    positive, error = clean(positive), clean(error)
    residual = clean(kernel*kernel-s.eye(4*size)-positive-error)
    if kernel != kernel.H or error != error.H or residual != s.zeros(4*size):
        raise ArithmeticError('the full spatial Wilson identity failed')
    return {'kernel': kernel, 'positive': positive, 'error': error,
            'identity_residual': residual, 'ordered_pairs': ordered,
            'spin_dimension': 4, 'site_internal_dimension': size}


def local_bound(magnetic_cost_density, cells=1, degree=128):
    """Translation-invariant coordinate marginal; full low-mode effect.

    Inputs are exact upper bounds for <S>/V, not just probabilities of a
    fixed-threshold rough event. The result includes arbitrary untouched
    references but does not assert a dynamical sign-response estimate.
    """
    density = rational(magnetic_cost_density, 'magnetic cost density')
    count = integer(cells, 'detector cells', 1)
    m = integer(degree, 'Legendre degree')
    # [400 |X| sqrt(6 <S>/V)]². Store squares rather than a rounded root.
    trace_square = 960000*count*count*density
    rank = 60*count
    epsilon = field_route.polynomial_bounds(m)['sign_error_on_admitted_gap']
    hs = {'rational_part': rank*epsilon*epsilon,
          'square_root_coefficient': (2+epsilon)**2,
          'square_root_argument': min(F(rank*rank), trace_square)}
    return {'detector_cells': count, 'one_body_local_rank': rank,
            'magnetic_cost_density_upper': density,
            'local_low_mode_trace_upper_squared': min(F(rank*rank), trace_square),
            'low_mode_threshold_absolute_HW': F(1, 2),
            'local_sign_minus_polynomial_HS_squared_upper': hs,
            'Legendre_degree': m, 'sign_approximation_error_on_good_spectrum': epsilon,
            'auxiliary_uniform_local_mode_probe_low_probability_upper_squared':
                min(F(1), trace_square/(rank*rank)),
            'translation_invariant_coordinate_marginal_required': True,
            'compressed_effect_not_truncated_kernel_projector': True,
            'global_gap_probability_claimed': False,
            'local_spectral_projector_decay_claimed': False,
            'interacting_Duhamel_record_bound_claimed': False}


def plan_bound(target_bits=20, degree=128, spacing='1', cells=1):
    plan = state.density_plan(target_bits, degree, spacing)
    density = plan['initial_magnetic_cost_per_cell_upper']+(
        plan['initial_electric_energy_per_cell_upper']
        +2*plan['full_Fock_kinetic_norm_per_cell_upper'])/plan['declared_plaquette_coefficient']
    output = local_bound(density, cells, degree)
    eta = plan['target_probability_upper']
    if output['local_low_mode_trace_upper_squared'] > 20*eta*cells*cells:
        raise ArithmeticError('the source density does not imply its advertised spectral bound')
    return {**output, 'requested_roughness_target': eta,
            'simpler_local_low_mode_trace_upper_squared': min(F((60*cells)**2), 20*eta*cells*cells),
            'same_plaquette_coefficient': plan['declared_plaquette_coefficient'],
            'same_per_link_word_degree': plan['single_link']['per_link_word_degree'],
            'all_time_on_same_prepared_polynomial_field': True,
            'bound_and_local_action_independent_of_total_volume': True,
            'complete_probe_effects': ['P_K Tr(P_X Pi_low P_X)/(60 |X|) P_K',
                                      'I_code minus low-mode effect'],
            'original_probability_sum': F(1),
            'conditioned_outcomes_used': False,
            'exact_detector_is_declared_primitive': True,
            'compiled_detector_error_included': False}


@lru_cache(maxsize=1)
def exact_controls():
    identity = s.eye(15)
    flat = decomposition((identity, identity, identity))
    if flat['kernel']**2 != s.eye(60):
        raise ArithmeticError('flat native Wilson symbol differs')
    first = native_weak_link(0, s.Rational(1, 2), s.sqrt(3)/2)
    second = native_weak_link(1, s.Rational(1, 2), s.sqrt(3)/2)
    rough = decomposition((first, second, identity))
    square = clean(rough['kernel']**2)
    low = clean((square-s.eye(60))*(square-3*s.eye(60))/3)
    if clean(low**2) != low or low.trace() != 16 or clean(square*low) != s.zeros(60):
        raise ArithmeticError('the complete native zero-mode projector differs')
    one = clean(-square*(square-3*s.eye(60))/2)
    three = clean(square*(square-s.eye(60))/6)
    if clean(low+one+three) != s.eye(60) or clean(square*(square-s.eye(60))*(square-3*s.eye(60))) != s.zeros(60):
        raise ArithmeticError('the full spectral menu is incomplete')
    spectrum = {0: low.trace(), 1: one.trace(), 3: three.trace()}
    plaquette = clean(first*second*first.H*second.H)
    cost = s.simplify(1-s.re(plaquette.trace())/15)
    probe = {str(k): s.Rational(v, 60) for k, v in spectrum.items()}
    return {'actual_native_basis_used': True, 'spatial_beta': 'gamma_tuple[3], not gamma5',
            'flat_symbol_dimension': 60, 'flat_squared_kernel_is_identity': True,
            'noncommuting_native_symbol': {
                'weak_links': ['(I+i sqrt3 sigma1)/2', '(I+i sqrt3 sigma2)/2', 'I'],
                'full_dimension': 60, 'ordered_pair_identity_checks': len(rough['ordered_pairs']),
                'H_W_squared_spectral_ranks': spectrum,
                'complete_auxiliary_symbol_probe_probabilities': probe,
                'positive_plaquette_cost_sum': cost,
                'low_mode_count': low.trace(),
                'universal_count_bound_squared': 960000*cost,
                'full_spectral_probability_sum': sum(probe.values())},
            'point_coordinate_normalizable_preparation_claimed': False,
            'native_weak_multiplicities_retained': True}


def serial(value):
    if isinstance(value, (F, s.Basic)):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(v) for v in value]
    return value


def certificate():
    return serial({'schema': 'oph.local_wilson_density.v1',
        'source_inputs': list(SOURCE_INPUTS),
        'primary_algebra': 'Hernandez-Jansen-Luscher, hep-lat/9808010v2, Appendix C',
        'trace_error_constant': 20, 'low_mode_count_squared_constant': 960000,
        'exact_native_controls': exact_controls(),
        'same_source_all_time_local_spectral_plan': plan_bound(),
        'all_time_bound_is_on_uncompressed_symbol_effect': True,
        'physical_mass_or_species_assignments': [],
        'remaining_response_premise': 'volume-uniform local interacting sign/polynomial Duhamel bound on the same quantum state'})
