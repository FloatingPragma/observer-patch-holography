"""Actual bulk transfer and spatial smearing bounds, without a momentum claim.

Infinite-volume physical algebras supply a comparison domain. Their state
is the actual finite Gauss state with invariant exterior blanks. No finite
periodic carrier or different preparation is substituted.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import product

import sympy as s

from model.particle_common_contract import contract
from model.particle_primary_field_algebra import algebra
from model.particle_primary_locality import locality
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS = (
    'model/particle_primary_spatial_algebra/spatial.py',
    'model/particle_primary_spatial_algebra/PROOF.md',
    'model/particle_common_contract/contract.py',
    'model/particle_common_contract/PROOF.md',
    'model/particle_primary_field_algebra/algebra.py',
    'model/particle_primary_field_algebra/PROOF.md',
    'model/particle_primary_locality/locality.py',
    'model/particle_primary_locality/PROOF.md',
    'model/particle_gauss_ground_state/PROOF.md',
    'model/particle_field_causality/PROOF.md',
    'model/particle_field_causality/bounds.py',
    'model/continuum_field_net/README.md',
    'Lean/Variational/TranslationInvariantComposedInstance.lean',
)


def _vector(value, name):
    if not isinstance(value, (tuple, list)) or len(value) != 3:
        raise ValueError(name+' requires three exact rational coordinates')
    return tuple(algebra._rational(x, name) for x in value)


def hat_integral(left, right, center='0', radius='1'):
    """Exact integral of the normalized triangular spatial density on one cell."""
    l = algebra._rational(left, 'left endpoint')
    r = algebra._rational(right, 'right endpoint')
    c = algebra._rational(center, 'hat center')
    R = ground._q(radius, 'positive hat radius', True)
    if r < l:
        raise ValueError('cell endpoints must be ordered')
    def primitive(x):
        y = (x-c)/R
        if y <= -1: return F(0)
        if y >= 1: return F(1)
        if y <= 0: return (y+1)**2/2
        return 1-(1-y)**2/2
    return primitive(r)-primitive(l)


def hat_weights(spacing='1/2', center=('0', '0', '0'), radius='1', max_cells=100000):
    """Finite exact spatial weights; large source plans do not allocate this menu."""
    a = ground._q(spacing, 'positive spacing', True)
    c = _vector(center, 'spatial center')
    R = ground._q(radius, 'positive support radius', True)
    cap = ground._i(max_cells, 'allocation cell cap', 1)
    intervals = [range(((x-R)/a).__floor__(), ((x+R)/a).__ceil__()) for x in c]
    count = 1
    for row in intervals: count *= len(row)
    if count > cap:
        raise ValueError('exact weight enumeration exceeds its explicit allocation cap')
    weights = {}
    for cell in product(*intervals):
        weight = F(1)
        for i, k in enumerate(cell):
            weight *= hat_integral(a*k, a*(k+1), c[i], R)
        if weight: weights[cell] = weight
    if sum(weights.values()) != 1 or any(w <= 0 for w in weights.values()):
        raise ArithmeticError('the exact positive spatial weights lost normalization')
    return weights


def rounded_translation(level, displacement):
    level = ground._i(level, 'nonnegative source level', 0)
    y = _vector(displacement, 'translation')
    a = F(1, 2**level)
    steps = tuple((x/a).__floor__() for x in y)
    actual = tuple(a*k for k in steps)
    return {'spacing': a, 'requested': y, 'lattice_steps': steps,
            'actual': actual, 'coordinate_error': tuple(x-z for x, z in zip(y, actual))}


def translation_group_defect(level, first, second, gradient_L1=('2', '2', '2')):
    g = tuple(ground._q(x, 'gradient L1 cap') for x in gradient_L1)
    if len(g) != 3: raise ValueError('three derivative bounds required')
    left = rounded_translation(level, first); right = rounded_translation(level, second)
    combined = rounded_translation(level, tuple(x+y for x,y in zip(left['requested'], right['requested'])))
    defect = tuple(x+y-z for x,y,z in zip(left['lattice_steps'], right['lattice_steps'], combined['lattice_steps']))
    if any(k not in (0, -1) for k in defect):
        raise ArithmeticError('rounded translation composition exceeds one lattice step')
    error = left['spacing']*sum(abs(k)*v for k,v in zip(defect,g))
    return {'lattice_group_defect': defect, 'smeared_operator_norm_difference_upper': error,
            'vanishing_upper': left['spacing']*sum(g),
            'finite_rounded_maps_are_exact_R3_group': False,
            'regular_quotient_maps_form_exact_R3_group': True}


def spatial_translation_bound(displacement, gradient_L1=('2', '2', '2'), norm_upper='1'):
    y = _vector(displacement, 'translation')
    g = tuple(ground._q(x, 'gradient L1 cap') for x in gradient_L1)
    if len(g) != 3: raise ValueError('three derivative bounds required')
    cap = ground._q(norm_upper, 'uniform norm cap')
    return min(2*cap, sum(abs(x)*v for x,v in zip(y,g)))


def crossing_incidence(side):
    """Independently enumerate every interaction crossing the cube's owner cut."""
    L = ground._i(side, 'box side', 2)
    def inside(x): return all(0 <= t < L for t in x)
    bonds = []; plaquettes = []
    for x in product(range(-1,L), repeat=3):
        for axis in range(3):
            y = tuple(k+int(i==axis) for i,k in enumerate(x))
            if inside(x) != inside(y): bonds.append((x,axis))
        for i,j in ((0,1),(0,2),(1,2)):
            owners = (x, tuple(k+int(q==i) for q,k in enumerate(x)),
                      tuple(k+int(q==j) for q,k in enumerate(x)))
            flags = [inside(c) for c in owners]
            if any(flags) and not all(flags): plaquettes.append((x,i,j))
    if len(bonds) != 6*L*L or len(plaquettes) != 12*L*L-3*L:
        raise ArithmeticError('open-boundary crossing incidence differs')
    return {'side': L, 'crossing_bonds': len(bonds),
            'crossing_plaquettes': len(plaquettes),
            'plaquette_count_upper': 12*L*L,
            'crossing_norm_times_spacing_upper': F(1272,7)*L*L,
            'dangling_owned_links_keep_their_electric_energy': True}


def bulk_transfer(member, support_radius='1', duration='1/1000', kind='plaquette'):
    """Bound finite/infinite dynamics on a normalized compact spatial smear.

    The state extension is physical, but no isometric finite physical
    algebra embedding is asserted. This is an ambient norm comparison.
    """
    contract.verify_member(member)
    if kind not in ('charge','plaquette'):
        raise ValueError('the declared smeared probes are charge and plaquette')
    source = locality.selected_source_bound(member)
    a = F(member['source']['spacing']); L = member['source']['side']
    R = ground._q(support_radius, 'positive physical support radius', True)
    time = ground._q(duration, 'source time extent')
    support = 1 if kind == 'charge' else 3
    collar = L//2-(R/a).__ceil__()-2
    cross = F(1272,7)*L*L/a
    W = source['family']['weighted_rate_times_spacing_upper']
    margin = F(collar,2)-2*W*time/a
    factor = 2*support*time*cross
    if time == 0 and collar >= 1:
        error = {'exact_zero': True, 'base': 2, 'negative_exponent': None}
    else:
        bits = -1
        if collar >= 1 and margin >= 0 and factor > 0:
            bits = max(-1, margin.__floor__()-locality.cone.ceil_log2(factor))
        error = {'exact_zero': False, 'base': 2, 'negative_exponent': bits}
    return contract.serial({
        'member_id': member['member_id'], **member['ids'],
        'kind': kind, 'support_radius': R, 'duration': time,
        'source_spacing': a, 'side': L, 'centered_collar_lower': collar,
        'physical_bulk_support_admitted': collar >= 1,
        'crossing_interaction_norm_upper': cross,
        'single_probe_support_cell_count': support,
        'spatial_smear_L1_norm_upper': F(1),
        'rational_exponential_margin': margin,
        'Duhamel_prefactor': factor,
        'operator_difference_norm_upper': error,
        'initial_state': 'exact finite Gauss state and conjugate reference, invariant exterior blanks',
        'time_group_comparison': 'homogeneous infinite Wilson time versus actual finite open source time',
        'finite_state_translation_averaged': False,
        'finite_physical_algebra_embedding_assumed_isometric': False,
        'infinite_state_assumed_stationary_before_limit': False,
        'full_spatial_matrices_executed': False,
    })


def noninvariant_state_control():
    """Time positivity and a continuous spatial action do not fix the state.

    C(S1), spatial rotation, time identity, evaluation at angle0. The four
    exact samples also supply a finite commuting-matrix countercontrol.
    """
    A = s.diag(1,0,-1,0); U = s.zeros(4)
    for i in range(4): U[(i+1)%4,i] = 1
    rho = s.diag(1,0,0,0); shifted = U*A*U.H
    Z = A-s.eye(4)
    return {'H': s.zeros(4), 'rho': rho, 'observable': A,
            'translated': shifted,
            'original_expectation': s.trace(rho*A),
            'translated_expectation': s.trace(rho*shifted),
            'null_vector_norm_squared': s.trace(rho*Z.H*Z),
            'translated_null_norm_squared': s.trace(rho*(U*Z*U.H).H*(U*Z*U.H)),
            'continuous_parent': 'C(S1), beta_y f(theta)=f(theta-y), omega(f)=f(0)',
            'canonical_spatial_GNS_implementation_exists': False,
            'control_is_primary_source': False}


def open_shift_counterexample():
    H = s.Matrix([[0,1,0],[1,0,1],[0,1,0]])
    shift = s.Matrix([[0,0,1],[1,0,0],[0,1,0]])
    defect = shift*H*shift.H-H
    return {'H_open': H, 'cyclic_permutation': shift, 'Hamiltonian_defect': defect,
            'squared_Hilbert_Schmidt_defect': s.trace(defect.H*defect),
            'cyclic_label_is_open_source_symmetry': False,
            'control_is_primary_source': False}


def common_payload(member):
    contract.verify_member(member)
    level = member['source']['level']
    return contract.serial({
        'member_id': member['member_id'], **member['ids'],
        'actual_bulk_transfer': bulk_transfer(member),
        'charge_seed': algebra.probe(member, kind='charge'),
        'plaquette_seed': algebra.probe(member, kind='plaquette'),
        'rounded_translation_example': rounded_translation(level, ('1/3','-2/5','1/7')),
        'group_defect_example': translation_group_defect(level, ('1/3','-2/5','1/7'), ('2/7','1/9','-1/6')),
        'regular_algebra_has_point_norm_continuous_R3_action': True,
        'regular_spatial_action_commutes_with_source_time': True,
        'limiting_regular_probe_word_moments_match_actual_finite_source': True,
        'stationary_limiting_source_state_and_nonnegative_time_energy': True,
        'norm_isomorphism_to_frozen_unsmeared_point_algebra_claimed': False,
        'selected_limiting_state_spatial_invariance_proved': False,
        'canonical_GNS_momentum_or_vacuum_momentum_proved': False,
        'finite_Fourier_probe_label_identified_as_momentum': False,
        'physical_zero_momentum_rest_energy_or_mass_proved': False,
        'nontrivial_spatial_field_or_species_proved': False,
        'extra_spatial_averaging_coin_inserted': False,
    })


def bind_member(member):
    result = common_payload(member)
    return {'schema':'oph.primary_spatial_regular_algebra.v1', 'result':result,
            'result_sha256':contract.digest(result)}


def verify_bound(member, receipt):
    if receipt != bind_member(member):
        raise ValueError('spatial algebra receipt differs from actual source or scope')


def certificate():
    counter = noninvariant_state_control(); opened = open_shift_counterexample()
    return contract.serial({
        'schema':'oph.primary_spatial_regular_algebra_certificate.v1',
        'selected_members':[bind_member(contract.member(n)) for n in (2,4,6)],
        'crossing_incidence_controls':[crossing_incidence(L) for L in (2,3,4)],
        'independent_spatial_weights':{'spacing':'1/2','radius':'1',
            'nonzero_cells':len(hat_weights()),'sum':sum(hat_weights().values())},
        'noninvariant_state_control':{
            'expectation_before':int(counter['original_expectation']),
            'expectation_after':int(counter['translated_expectation']),
            'null_norm_before':int(counter['null_vector_norm_squared']),
            'null_norm_after':int(counter['translated_null_norm_squared']),
            'is_primary':False},
        'open_shift_counterexample':{'squared_HS_defect':int(opened['squared_Hilbert_Schmidt_defect']),
                                    'is_primary':False},
        'actual_primary_momentum_mass_or_GeV_output':None,
    })
