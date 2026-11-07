"""One biased native edge with gauge transport and protected memory.

The matter/link registers and flux-dependent holds define an explicit law
extension. The memory charge is not identified with a particle's statistics.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from ..gauge_defect_transport import transport
from . import reference

TARGET_PAIR = (0, 5)
SOURCE_PAIR = (5, 11)
SEAM = (0, 11)
SOURCE_FILES = (
    'model/protected_reference/gauge_join.py',
    'model/protected_reference/test_gauge_join.py',
    'model/protected_reference/GAUGE_JOIN.md',
    *transport.SOURCE_FILES,
)


def native_edge():
    """Actual opposite moves, ordered parity and endpoint coin weights."""
    p = reference.parameters()
    x, y = (reference.exchange.INDEX[pair] for pair in (TARGET_PAIR, SOURCE_PAIR))
    rows = []
    for law in ('C', 'E'):
        for coin, origin, target in ((1, y, x), (0, x, y)):
            selected = [row for row in reference.protection.branches(law)
                        if row.seam == SEAM and row.coin == coin and origin in row.inputs
                        and row.outputs[row.inputs.index(origin)] == target]
            if len(selected) != 1 or selected[0].signs[selected[0].inputs.index(origin)] != -1:
                raise ArithmeticError('the native move lost its crossed spectator parity')
            row = selected[0]
            rows.append({'law': law, 'coin': coin, 'raw': row.raw, 'orientation_flip': 1})
    pi = reference.local.reference(p)
    forward = reference.local.coin_probability(p, SEAM, 1)
    reverse = reference.local.coin_probability(p, SEAM, 0)
    if forward+reverse != 1 or forward*pi[y] != reverse*pi[x]:
        raise ArithmeticError('the actual endpoint coin lost detailed balance')
    return {'target_index': x, 'source_index': y,
            'target_reference_weight': pi[x], 'source_reference_weight': pi[y],
            'forward_probability': forward, 'reverse_probability': reverse,
            'native_columns': rows}


@lru_cache(maxsize=5)
def active_maps(charge):
    """Active order: position, orientation, internal/link.

    The other 64 source configurations are retained. Both moved outcomes
    vanish there; both held outcomes are identity there.
    """
    value = transport.literal_transport(charge, completed=False)
    b = value['normalized']
    n = b.rows
    flip = s.ImmutableSparseMatrix([[0, 1], [1, 0]])
    dressed = s.ImmutableSparseMatrix(s.kronecker_product(flip, b))
    maps = dict(transport.complete_partial_transport(dressed))
    maps['charge'] = s.ImmutableSparseMatrix(s.kronecker_product(s.SparseMatrix.eye(2), flip, s.SparseMatrix.eye(n)))
    maps['internal_link_dimension'] = n
    maps['matter_dimension'] = value['dimension']
    maps['full_source_memory_internal_link_dimension'] = 132*n
    maps['complement_dimension'] = 128*n
    return maps


def active_reference(charge, plus_weight=F(1, 3)):
    """Unnormalized active block of the full normalized stationary state."""
    a = reference._probability(plus_weight)
    edge = native_edge()
    n = active_maps(charge)['internal_link_dimension']
    sigma = s.ImmutableSparseMatrix([[s.Rational(1, 2), s.Rational(2*a-1)/2],
                                    [s.Rational(2*a-1)/2, s.Rational(1, 2)]])
    source = s.ImmutableSparseMatrix(s.diag(s.Rational(edge['target_reference_weight']), s.Rational(edge['source_reference_weight'])))
    return s.ImmutableSparseMatrix(s.kronecker_product(source, sigma, s.SparseMatrix.eye(n)/n))


def instrument_on_active(charge, matrix):
    maps, edge = active_maps(charge), native_edge()
    result = s.SparseMatrix.zeros(maps['forward'].rows)
    for direction, probability in (('forward', edge['forward_probability']),
                                    ('reverse', edge['reverse_probability'])):
        for key in (direction, direction+'_hold'):
            k = maps[key]
            result += s.Rational(probability)*k*matrix*k.H
    return transport.link.clean(result)


def case_certificate(charge):
    maps = active_maps(charge)
    value = transport.literal_transport(charge, completed=False)
    edge = native_edge()
    size = maps['forward'].rows
    identity = s.ImmutableSparseMatrix(s.SparseMatrix.eye(size))
    effects = s.SparseMatrix.zeros(size)
    q = maps['charge']
    for direction, probability in (('forward', edge['forward_probability']),
                                    ('reverse', edge['reverse_probability'])):
        for key in (direction, direction+'_hold'):
            k = maps[key]
            effects += s.Rational(probability)*k.H*k
            if not transport.link.zero(q*k-k*q):
                raise ArithmeticError('gauge transport changed the protected native charge')
    if not transport.link.zero(effects-identity):
        raise ArithmeticError('the complete four-outcome gauge instrument is incomplete')
    for a in (F(), F(1, 3), F(1)):
        state = active_reference(charge, a)
        if not transport.link.zero(instrument_on_active(charge, state)-state):
            raise ArithmeticError('the same native reference is not stationary')
    # In units log(25/8), source K has active levels 1 and 2. Its memory
    # term commutes with all maps on every supported stationary mixture.
    level = s.ImmutableSparseMatrix(size, size, {(i, i): 1+i//(size//2) for i in range(size)})
    for key, gap in (('forward', -1), ('reverse', 1),
                     ('forward_hold', 0), ('reverse_hold', 0)):
        if not transport.link.zero(level*maps[key]-maps[key]*level-gap*maps[key]):
            raise ArithmeticError('the same-state modular frequency of a record changed')
    if value['endpoint_checks'] != 72:
        raise ArithmeticError('the endpoint gauge-current inventory changed')
    return {'primitive_charge': charge, 'matter_dimension': maps['matter_dimension'],
            'link_dimension': 60, 'active_dimension': size,
            'retained_inactive_dimension': maps['complement_dimension'],
            'full_dimension': maps['full_source_memory_internal_link_dimension'],
            'transport_rank_per_memory_character': int(s.trace(value['source_projector'])),
            'complete_outcomes': 4, 'complete_effect_identity': True,
            'all_maps_preserve_protected_charge': True,
            'stationary_reference_checked_including_memory_coherences': True,
            'same_state_modular_gap_multiples': {'forward': -1, 'reverse': 1,
                                                'forward_hold': 0, 'reverse_hold': 0},
            'literal_endpoint_Gauss_intertwiner_checks': value['endpoint_checks']}


def singlet_certificate(charge):
    """Stationary reducing corner with a declared dual reference at y.

    Matrix columns are amplitudes of one entangled vector, not classical
    alternatives. The source marginal is conditioned on the active pair.
    """
    value = transport.literal_transport(charge, completed=False)
    d, n, maps = value['dimension'], value['normalized'].rows, value['maps']
    initial = s.ImmutableSparseMatrix(2*n, d, {(n+60*b, b): 1/s.sqrt(d) for b in range(d)})
    target = maps['forward']*initial
    if (s.trace(initial.H*initial) != 1 or s.trace(target.H*target) != 1
            or s.trace(initial.H*target) != 0):
        raise ArithmeticError('the neutralized source and target are not orthonormal')
    for key, incoming, expected in (
            ('forward', initial, target), ('reverse', target, initial),
            ('forward_hold', target, target), ('reverse_hold', initial, initial),
            ('forward', target, s.SparseMatrix.zeros(2*n, d)),
            ('reverse', initial, s.SparseMatrix.zeros(2*n, d)),
            ('forward_hold', initial, s.SparseMatrix.zeros(2*n, d)),
            ('reverse_hold', target, s.SparseMatrix.zeros(2*n, d))):
        if not transport.link.zero(maps[key]*incoming-expected):
            raise ArithmeticError('the full singlet span is not a reducing transport corner')
    gauge = transport.neutralized_control(charge)
    if gauge['endpoint_singlet_checks'] != 48:
        raise ArithmeticError('the declared dual reference lost endpoint singlet support')
    edge = native_edge()
    mass = edge['target_reference_weight']+edge['source_reference_weight']
    target_weight = edge['target_reference_weight']/mass
    source_weight = edge['source_reference_weight']/mass
    if edge['forward_probability']*source_weight != edge['reverse_probability']*target_weight:
        raise ArithmeticError('the conditioned native weights lost singlet stationarity')
    logical = (s.Matrix([[0, 1], [0, 0]]), s.diag(1, 0),
               s.Matrix([[0, 0], [1, 0]]), s.diag(0, 1))
    weights = (edge['forward_probability'],)*2+(edge['reverse_probability'],)*2
    target_density = s.diag(s.Rational(target_weight), s.Rational(source_weight))
    for i in range(2):
        for j in range(2):
            matrix_unit = s.SparseMatrix(2, 2, {(i, j): 1})
            output = sum((s.Rational(w)*k*matrix_unit*k.H for w, k in zip(weights, logical, strict=True)), s.zeros(2))
            if output != (target_density if i == j else s.zeros(2)):
                raise ArithmeticError('the supported native attempt is not the claimed replacement channel')
    return {'primitive_charge': charge, 'dual_reference_dimension': d,
            'full_source_memory_link_matter_reference_dimension': 132*n*d,
            'endpoint_singlet_checks': gauge['endpoint_singlet_checks'],
            'two_orthonormal_singlets_form_reducing_corner': True,
            'target_probability': str(target_weight), 'source_probability': str(source_weight),
            'corner_state': '(25|psi_x><psi_x|+8|psi_y><psi_y|)/33 tensor sigma_orientation',
            'state_is_annihilated_by_endpoint_Gauss_generators': True,
            'full_native_source_marginal_replaced_by_active_conditional_marginal': True,
            'dual_reference_and_singlet_preparation_are_declared': True,
            'literal_singlet_preparation_pulses_executed': False,
            'one_unread_attempt_exact_preparation_in_each_sharp_charge_corner': True,
            'preparation_domain': 'arbitrary supported two-singlet state with untouched reference and one fixed Q character',
            'preparation_archive_traced_only_for_replacement_comparison': True,
            'clock_domain': 'supported singlet and memory-charge corner only',
            'forward_modular_energy_change': 'log(8/25)'}


def certificate():
    edge = native_edge()
    return {
        'schema': 'oph.protected_reference.gauge_join.v1',
        'source_pair': list(SOURCE_PAIR), 'target_pair': list(TARGET_PAIR), 'native_seam': list(SEAM),
        'crossed_spectator_port': 5, 'native_columns': edge['native_columns'],
        'forward_coin_probability': str(edge['forward_probability']),
        'reverse_coin_probability': str(edge['reverse_probability']),
        'target_reference_weight': str(edge['target_reference_weight']),
        'source_reference_weight': str(edge['source_reference_weight']),
        'forward_modular_energy_change': 'log(8/25)',
        'stationary_reference': 'rho_* tensor (a P_plus+(1-a)P_minus) tensor I_(60d)/(60d)',
        'full_source_dimension': 66, 'persistent_orientation_dimension': 2,
        'inactive_native_configurations_retained_as_identity': 64,
        'rows': [case_certificate(charge) for charge in transport.NAMES],
        'Gauss_singlet_stationary_corners': [singlet_certificate(charge) for charge in transport.NAMES],
        'scope': {'derived_native_coin_and_modular_energy_ratios': True,
                  'faithful_full_reference_is_Gauss_singlet_supported': False,
                  'matter_and_link_registers_are_declared': True,
                  'flux_dependent_gauge_extension_is_original_C_E_law': False,
                  'single_moving_internal_fiber_is_two_particle_colored_code': False,
                  'all_native_seams_have_a_common_gauge_network': False,
                  'stationary_state_is_unique_or_mixing_proved': False,
                  'same_stationary_reference_selects_coherent_hopping_Hamiltonian': False,
                  'particle_species_or_masses_identified': False},
    }
