"""An exact protected-label/interaction join on one declared finite program.

The controller is an architectural choice. It makes a source invariant
observable through a complete interacting scalar instrument, while keeping
that invariant distinct from scalar occupation and particle identity.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import json

import sympy as s

from model.defect_protection import protection
from model.exchange_sector import statistics as exchange
from model.common_source import slow_geometry, native_gauge_object
from model.particle_observer_outputs import higgs_records as field
from model.particle_observer_outputs import holonomy_transport as holonomy


SOURCE_INPUTS = (
    'model/defect_protection/protection.py',
    'model/defect_protection/README.md',
    'model/exchange_sector/statistics.py',
    'model/common_source/native_gauge_object.py',
    'model/protected_defect_catalog/PROOF.md',
    'model/gauge_defect_transport/PROOF.md',
    'model/standard_model_source_map/FERMION_ATTACHMENT.md',
    'model/standard_model_source_map/KINETIC_ATTACHMENT.md',
    'model/observer_particles/spin_matter_join.py',
    'model/particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md',
    'model/particle_observer_outputs/higgs_records.py',
    'model/particle_observer_outputs/HIGGS_RECORDS.md',
    'model/particle_observer_outputs/holonomy_transport.py',
    'model/particle_observer_outputs/HOLONOMY_TRANSPORT.md',
)


@lru_cache(maxsize=1)
def scalar_source():
    """The literal full 307-state quantum scalar, at declared pair cost zero."""
    return field.sector_data(s.diag(1, 2, 3), 0)


def strength(charge):
    if type(charge) is not int or charge not in (-1, 1):
        raise ValueError('a sharp source charge +1 or -1 is required')
    return (3-charge)//2


def common_time():
    return s.pi/(6*s.sqrt(scalar_source()['g_squared']))


@lru_cache(maxsize=2)
def scalar_tick(charge):
    """Same time, distinct declared charge-dependent Hamiltonian blocks."""
    base = scalar_source()
    angle = strength(charge)*s.pi/6
    active = base['vacuum_projector']+base['bright_projector']
    return field._matrix(base['dark_projector']+s.cos(angle)*active
        -s.I*s.sin(angle)*base['vertex']/s.sqrt(base['g_squared']))


@lru_cache(maxsize=1)
def controller():
    """All 614 charge/field columns, including both 305-state dark complements."""
    base = scalar_source()
    x = s.ImmutableMatrix([[0, 1], [1, 0]])
    sectors = {q: (s.eye(2)+q*x)/2 for q in (-1, 1)}
    h = s.ImmutableSparseMatrix(s.kronecker_product((3*s.eye(2)-x)/2, base['vertex']))
    u = s.ImmutableSparseMatrix(sum((s.kronecker_product(sectors[q], scalar_tick(q))
                                    for q in (-1, 1)), s.zeros(614)))
    transform = s.ImmutableSparseMatrix(s.kronecker_product(s.Matrix([[1, 1], [1, -1]])/s.sqrt(2), s.eye(307)))
    return {'H': h, 'U': u, 'Q': s.ImmutableSparseMatrix(s.kronecker_product(x, s.eye(307))),
            'charge_projectors': sectors, 'charge_basis_transform': transform,
            'native_source_dimension': 132, 'field_dimension': 307,
            'complete_native_field_dimension': 132*307,
            'tick_time': common_time(), 'field': base}


def complete_detector_maps():
    """Read Q and scalar occupation, preserving every quantum complement."""
    data = controller()
    return {(q, a): s.ImmutableSparseMatrix(s.kronecker_product(data['charge_projectors'][q], p)*data['U'])
            for q in (-1, 1) for a, p in field.measurement_projectors(data['field']).items()}


def sector_record_probabilities():
    base = scalar_source()
    result = {}
    for q in (-1, 1):
        evolved = scalar_tick(q)*base['vacuum']
        probabilities = {name: field._clean((evolved.H*p*evolved)[0])
                         for name, p in field.measurement_projectors(base).items()}
        if sum(probabilities.values()) != 1:
            raise ArithmeticError('a complete source sector lost detector weight')
        result[q] = probabilities
    return result


@lru_cache(maxsize=1)
def native_branch_join():
    """Actual all-column C/E charge decomposition; no dense 40524-square matrix.

    In the Q basis, every cover branch is its original positive/negative
    exchange-character branch. Identity on all 307 field columns and the
    controlled scalar block then give the full tensor-product intertwiner.
    """
    result = {}
    for law in ('C', 'E'):
        counts = [0]*132
        tested, nonzero, raw = 0, 0, []
        for row in protection.branches(law):
            mapping = dict(protection.lifted_pairs(row))
            signed = dict(zip(row.inputs, zip(row.outputs, row.signs, strict=True), strict=True))
            for column in range(132):
                flipped = (column+66) % 132
                image = mapping.get(column)
                if mapping.get(flipped) != (None if image is None else (image+66) % 132):
                    raise ArithmeticError('a literal source column changes its protected charge')
                counts[column] += image is not None
                tested += 1
                nonzero += image is not None
            for source in range(66):
                expected = signed.get(source)
                for q in (-1, 1):
                    actual = {}
                    for orientation, amplitude in ((0, 1), (1, q)):
                        output = mapping.get(66*orientation+source)
                        if output is not None:
                            actual[output] = amplitude
                    target = {} if expected is None else {
                        expected[0]: (1 if q == 1 else expected[1]),
                        66+expected[0]: q*(1 if q == 1 else expected[1])}
                    if actual != target:
                        raise ArithmeticError('charge-basis decoding changes an original exchange column')
            raw.append({'law': law, 'seam': list(row.seam), 'coin': row.coin,
                        'raw': row.raw, 'abort': row.abort,
                        'input_column_count': len(row.inputs),
                        'original_uniform_source_probability': F(len(row.inputs), 66*60)})
        if counts != [60]*132 or sum(r['original_uniform_source_probability'] for r in raw) != 1:
            raise ArithmeticError('the uniform seam/coin source is not complete')
        result[law] = {'all_input_columns_checked': tested,
                       'nonzero_partial_injection_columns': nonzero,
                       'all_raw_slots': raw, 'zero_raw_slots': sum(not r['input_column_count'] for r in raw),
                       'field_columns_per_native_column': 307,
                       'complete_joint_column_identities': tested*307,
                       'every_branch_commutes_with_Q_and_controlled_H': True}
    return result


def native_joint_records(law):
    """Uniform loads, orientation0, scalar empty; every raw/Q/scalar event."""
    source = native_branch_join()[law]
    probabilities = sector_record_probabilities()
    rows = []
    for raw in source['all_raw_slots']:
        for q in (-1, 1):
            for outcome, probability in probabilities[q].items():
                rows.append(dict(raw, charge=q, scalar_record=outcome,
                    joint_probability=field._clean(s.Rational(raw['original_uniform_source_probability'])*probability/2)))
    if sum(row['joint_probability'] for row in rows) != 1:
        raise ArithmeticError('joint native/field records were conditioned or dropped')
    return rows


def rotation_and_parity_boundaries():
    """Reuse actual port actions; do not call an exchange sign a spin lift."""
    rotations = tuple(tuple(map(int, p)) for p in slow_geometry.source().rotations())
    column_checks, order_two = 0, None
    for g in rotations:
        images, signs = exchange.rotation_action(g, -1)
        cover = tuple(66*(r ^ int(signs[i] == -1))+images[i]
                      for r in (0, 1) for i in range(66))
        for column in range(132):
            if cover[(column+66) % 132] != (cover[column]+66) % 132:
                raise ArithmeticError('native port rotation does not preserve the coupled charge')
            column_checks += 1
        if g != tuple(range(12)) and all(g[g[i]] == i for i in range(12)):
            if any(cover[cover[i]] != i for i in range(132)):
                raise ArithmeticError('the inherited source rotation became a projective spin lift')
            order_two = g
    if order_two is None:
        raise ArithmeticError('the native rotation group lost its order-two control')
    axes = native_gauge_object.axis_action(order_two)
    if any(axes[axes[i]] != i for i in range(6)):
        raise ArithmeticError('the protected six-object action changed its ordinary rotation law')
    parity = tuple((-1)**mask.bit_count() for _, mask in scalar_source()['basis'])
    if set(parity) != {1}:
        raise ArithmeticError('the complete neutral scalar sector is not CAR even')
    return {'native_rotation_Q_commutation_columns': column_checks,
            'source_order_two_rotation': list(order_two),
            'ordered_cover_rotation_square': '+I132',
            'protected_axis_object_rotation_square': '+I6',
            'spin_half_lift_of_same_half_turn_square': '-I2',
            'intertwiner_preserving_these_two_rotation_squares_must_vanish': True,
            'complete_neutral_field_CAR_parity': '+I307',
            'Q_minus_to_neutral_field_odd_parity_intertwiner_exists': False,
            'native_exchange_character_identified_with_spatial_spin': False,
            'field_current_and_spin_Ward_proofs_reused': True}


def occupation_obstruction():
    data = scalar_source()
    commutator = data['H']*data['occupation']-data['occupation']*data['H']
    norm_squared = field._clean(((commutator*data['vacuum']).H*(commutator*data['vacuum']))[0])
    if norm_squared != data['g_squared'] or sector_record_probabilities()[1]['occupied'] != s.Rational(1, 4):
        raise ArithmeticError('the actual source occupation ceased to be dynamical')
    return {'H_NH_commutator_on_vacuum_norm_squared': norm_squared,
            'empty_input_occupied_probability_at_common_tick': s.Rational(1, 4),
            'exactly_protected_charge_can_equal_dynamic_scalar_occupation': False,
            'no_go_domain': 'an intertwiner preserving the same evolution, sharp charge and the nonzero bright transition',
            'legitimate_preserved_label': 'Q becomes the conserved direct-sum interaction label, not NH'}


def geometry_join():
    """Use the actual two geometric work states and every retained field input."""
    path = holonomy.graphical_path()
    stages = path['stages']
    region = 1 << path['start_face']
    first_edge = path['dual_path'][0]['edge']
    initial = (0, 1)
    pair = (holonomy.edge_syndrome(first_edge), 0)
    if holonomy.work_exchange(*initial, first_edge) != pair:
        raise ArithmeticError('the reused local birth operation differs')
    # These two gauge-orbit basis vectors have equal work-plus-defect cost.
    if (holonomy.energy_units(*initial), holonomy.energy_units(*pair)) != (2, 2):
        raise ArithmeticError('the two geometric code states have unequal H_cost')
    if (initial[0] & region).bit_count() % 2 or (pair[0] & region).bit_count() % 2 != 1:
        raise ArithmeticError('the actual birth does not change the selected boundary parity')
    transport_checks = 0
    for stage in stages[2:-1]:
        if stage['boundary_holonomy_sign'] != -1 or not stage['update_preserves_region']:
            raise ArithmeticError('the supplied protected path crosses its own selected region')
        if stage['total_energy_in_J'] != 2:
            raise ArithmeticError('the protected path changes the original work cost')
        transport_checks += 1
    prior = {1: s.Rational(9, 25), -1: s.Rational(16, 25)}
    records = {(q, outcome): field._clean(prior[q]*p)
               for q, outcomes in sector_record_probabilities().items() for outcome, p in outcomes.items()}
    if sum(records.values()) != 1:
        raise ArithmeticError('geometric creation probabilities lost an outcome')
    return {'same_carrier_start_face': path['start_face'], 'other_face': path['separated_face'],
            'first_creation_edge': first_edge,
            'code_gauge_orbit_states': {'Q_plus': {'curvature_bits': initial[0], 'work': 1},
                                      'Q_minus': {'curvature_bits': pair[0], 'work': 0}},
            'isometry_field_columns_per_charge': 307, 'complete_join_columns': 614,
            'geometric_code_cost': '2J I614',
            'full_gauge_link_work_field_dimension': (1 << 20)*307,
            'full_ambient_link_work_field_dimension': (1 << 31)*307,
            'every_other_curvature_work_field_input_retained': True,
            'outside_Gauss_outcome_retained': True,
            'native_Q_and_face_Q_identified': False,
            'all_protected_path_steps_checked': transport_checks,
            'same_clock_complete_joint_records': [{'charge': q, 'scalar_record': outcome, 'probability': p}
                                                 for (q, outcome), p in records.items()],
            'fixed_region_is_part_of_apparatus': True,
            'rotation_covariance_requires_rotating_region_and_its_readout_together': True,
            'protected_moves_preserve_full_H_cost_plus_H_Q': True,
            'birth_and_annihilation_preserve_full_H_cost_plus_H_Q': False}


def geometry_rotation_covariance():
    """The actual sixty carrier rotations move the region and its apparatus."""
    geometry = holonomy.complex_data()
    edge_index = {frozenset(edge): i for i, edge in enumerate(geometry['edges'])}
    face_index = {frozenset(face): i for i, face in enumerate(geometry['faces'])}
    tested = 0
    for rotation in slow_geometry.source().rotations():
        edges = tuple(edge_index[frozenset(rotation[v] for v in edge)]
                      for edge in geometry['edges'])
        faces = tuple(face_index[frozenset(rotation[v] for v in face)]
                      for face in geometry['faces'])
        for face in range(20):
            for edge in range(30):
                before = (geometry['face_masks'][face] >> edge) & 1
                after = (geometry['face_masks'][faces[face]] >> edges[edge]) & 1
                if before != after:
                    raise ArithmeticError('a native rotation changes the actual face-edge incidence')
                tested += 1
    # Equality on each elementary face/edge extends by F2 bilinearity to all
    # regions and all link inputs, including their retained references.
    return {'actual_rotation_incidence_columns_checked': tested,
            'all_regions_and_all_link_inputs_by_bilinearity': True,
            'simultaneously_rotated_region_charge': 'R_g Q_R R_g^*=Q_{gR}',
            'fixed_region_global_rotation_invariance_claimed': False}


def preparation_energy():
    """Birth work for the coupled Hamiltonian, with measurement backaction explicit."""
    # In charge coordinates, the literal work-exchange pulse is cI-isX.
    base = scalar_source()
    x = s.Matrix([[0, 1], [1, 0]])
    birth = s.Rational(3, 5)*s.eye(2)-s.I*s.Rational(4, 5)*x
    strength_matrix = s.diag(1, 2)
    if birth.H*birth != s.eye(2):
        raise ArithmeticError('the inherited geometric birth pulse lost unitarity')
    comm = strength_matrix*birth-birth*strength_matrix
    if comm == s.zeros(2):
        raise ArithmeticError('the source controller incorrectly commutes with birth')
    initial = s.Matrix([1, 0]); prepared = birth*initial
    initial_square = (initial.H*strength_matrix**2*initial)[0]*base['g_squared']
    prepared_square = field._clean((prepared.H*strength_matrix**2*prepared)[0]*base['g_squared'])
    return {'scalar_mean_energy_before_and_after_birth': 0,
            'scalar_energy_second_moment_before_birth': initial_square,
            'scalar_energy_second_moment_after_birth': prepared_square,
            'birth_commutes_with_H_cost': True,
            'birth_commutes_with_controlled_scalar_H': False,
            'separate_two_time_energy_audit_work': [
                {'work': 0, 'probability': s.Rational(9, 25)},
                {'work': s.sqrt(base['g_squared']), 'probability': s.Rational(8, 25)},
                {'work': -s.sqrt(base['g_squared']), 'probability': s.Rational(8, 25)}],
            'two_time_energy_audit_dephases_initial_scalar': True,
            'two_time_audit_is_same_occupation_protocol': False,
            'preparation_service': 'pause or refocus scalar dynamics, apply birth pulse, restore stated H_Q; controller work is part of the declared program',
            'closed_autonomous_energy_conserving_birth_claimed': False}


def _serial(value):
    if isinstance(value, dict):
        return {str(k): _serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serial(v) for v in value]
    if isinstance(value, (s.Basic, F)):
        return str(value)
    return value


def certificate():
    base = scalar_source()
    return _serial({'schema': 'oph.finite_defect_field_attachment.v1',
        'source_inputs': list(SOURCE_INPUTS),
        'declared_coupled_H': 'H_Q=(3I-Q)V(Y)/2, Y=diag(1,2,3), d=0',
        'one_shared_time': common_time(),
        'sector_interaction_strengths': {'+1': 1, '-1': 2},
        'sector_bright_gaps': {q: 2*strength(q)*s.sqrt(base['g_squared']) for q in (-1, 1)},
        'same_time_scalar_records': sector_record_probabilities(),
        'full_native_branch_attachment': native_branch_join(),
        'native_complete_records': {law: native_joint_records(law) for law in ('C', 'E')},
        'geometry_attachment': geometry_join(),
        'geometry_rotation_covariance': geometry_rotation_covariance(),
        'preparation_energy': preparation_energy(),
        'rotation_and_statistics': rotation_and_parity_boundaries(),
        'occupation_obstruction': occupation_obstruction(),
        'architecture_choice': 'charge-controlled coupling is declared; it is not implied by protected topology',
        'native_repair_alone_selected_this_coupling': False,
        'same_source_state_clock_operator_and_original_weights': True,
        'physical_particle_assignments': [], 'physical_mass_predictions': [],
        'measured_particle_inputs': []})


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2, sort_keys=True))
