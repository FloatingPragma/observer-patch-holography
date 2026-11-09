"""Complete protected-source/field joins and controls against false identities."""
from pathlib import Path

import pytest
import sympy as s

from model.particle_defect_attachment import attachment as a


def zero(matrix):
    return all(s.simplify(value) == 0 for value in matrix.todok().values())


def test_literal_full_scalar_source_and_dark_complement():
    data = a.scalar_source()
    assert data['Y'] == s.diag(1, 2, 3)
    assert data['pair_cost'] == 0
    assert data['H'] == data['vertex']
    assert data['g_squared'] == 56
    assert len(data['basis']) == 307
    assert s.trace(data['dark_projector']) == 305
    assert zero(data['H']*data['dark_projector'])
    assert {mask.bit_count() % 2 for _, mask in data['basis']} == {0}


def test_all_614_controller_columns_decode_same_charge_and_hamiltonian():
    data = a.controller()
    h, u, q, f = (data[key] for key in ('H', 'U', 'Q', 'charge_basis_transform'))
    assert zero(h-h.H)
    assert zero(h*q-q*h)
    assert zero(u.H*u-s.eye(614))
    assert zero(f.H*f-s.eye(614))
    v = data['field']['vertex']
    assert zero(f.H*h*f-s.diag(v, 2*v))
    assert zero(f.H*u*f-s.diag(a.scalar_tick(1), a.scalar_tick(-1)))
    assert zero(f.H*q*f-s.diag(s.eye(307), -s.eye(307)))


def test_complete_detector_keeps_dark_inputs_and_an_entangled_reference():
    data = a.controller(); field = data['field']
    maps = a.complete_detector_maps()
    total = s.zeros(614, 614)
    for operator in maps.values():
        total += operator.H*operator
    assert zero(total-s.eye(614))
    dark = a.field.gauge_invariant_dark_probe(field)['vector']
    plus, minus = s.Matrix([1, 1])/s.sqrt(2), s.Matrix([1, -1])/s.sqrt(2)
    # Coefficient columns refer to two orthogonal, untouched reference states.
    psi = s.ImmutableSparseMatrix(s.Matrix.hstack(
        s.kronecker_product(plus, field['vacuum'])/s.sqrt(2),
        s.kronecker_product(minus, dark)/s.sqrt(2)))
    reference = s.zeros(2)
    probabilities = {}
    for label, operator in maps.items():
        output = operator*psi
        reference += output.T*output.conjugate()
        probabilities[label] = s.simplify(s.trace(output.H*output))
    assert zero(reference-s.eye(2)/2)
    assert probabilities[(-1, 'occupied')] == s.Rational(1, 2)
    assert probabilities[(-1, 'empty')] == 0
    # A detector that silently discarded D would lose half this original input.
    retained_bright = s.kronecker_product(s.eye(2), field['vacuum_projector']+field['bright_projector'])
    assert s.simplify(s.trace(psi.H*retained_bright*psi)) == s.Rational(1, 2)


def test_actual_native_branch_columns_original_weights_and_zero_slots():
    joined = a.native_branch_join()
    for law, raw_count, zero_count, columns in [('C', 240, 120, 31680), ('E', 480, 300, 63360)]:
        source = joined[law]
        assert source['all_input_columns_checked'] == columns
        assert source['complete_joint_column_identities'] == columns*307
        assert len(source['all_raw_slots']) == raw_count
        assert source['zero_raw_slots'] == zero_count
        records = a.native_joint_records(law)
        assert len(records) == 4*raw_count
        assert sum(row['joint_probability'] for row in records) == 1
        for index, original in enumerate(source['all_raw_slots']):
            assert sum(row['joint_probability'] for row in records[4*index:4*index+4]) == original['original_uniform_source_probability']
        assert any(row['abort'] for row in records)


def test_sharp_sector_outputs_require_the_declared_controller():
    records = a.sector_record_probabilities()
    assert records == {-1: {'empty': s.Rational(1, 4), 'occupied': s.Rational(3, 4)},
                       1: {'empty': s.Rational(3, 4), 'occupied': s.Rational(1, 4)}}
    assert a.common_time() == s.pi/(12*s.sqrt(14))
    # Equal-strength coupling cannot distinguish Q through this field detector.
    field = a.scalar_source(); u = a.scalar_tick(1)
    vacuum = field['vacuum']; probability = (vacuum.H*u.H*field['occupation']*u*vacuum)[0]
    assert s.simplify(probability) == s.Rational(1, 4)
    assert sum(abs(records[1][k]-records[-1][k]) for k in records[1])/2 == s.Rational(1, 2)
    for invalid in (0, 2, True, s.Integer(1)):
        with pytest.raises(ValueError):
            a.strength(invalid)


def test_actual_geometric_code_and_all_original_creation_records():
    joined = a.geometry_join()
    assert joined['complete_join_columns'] == 614
    assert joined['full_gauge_link_work_field_dimension'] == (1 << 20)*307
    assert joined['full_ambient_link_work_field_dimension'] == (1 << 31)*307
    assert joined['outside_Gauss_outcome_retained']
    assert not joined['native_Q_and_face_Q_identified']
    assert joined['all_protected_path_steps_checked'] > 0
    probabilities = {(row['charge'], row['scalar_record']): row['probability']
                     for row in joined['same_clock_complete_joint_records']}
    assert probabilities == {(1, 'empty'): s.Rational(27, 100), (1, 'occupied'): s.Rational(9, 100),
                             (-1, 'empty'): s.Rational(16, 100), (-1, 'occupied'): s.Rational(48, 100)}
    assert joined['protected_moves_preserve_full_H_cost_plus_H_Q']
    assert not joined['birth_and_annihilation_preserve_full_H_cost_plus_H_Q']


def test_protected_transport_commutes_with_cost_and_region_on_every_local_column():
    path = a.holonomy.graphical_path(); region = 1 << path['start_face']
    data = a.holonomy.complex_data()
    for edge in range(30):
        if not a.holonomy.preserves_region(1 << edge, region):
            continue
        f, g = data['edge_faces'][edge]
        spectator = next(i for i in range(20) if i not in (f, g))
        for first in range(2):
            for second in range(2):
                bits = (first << f) | (second << g) | ((first ^ second) << spectator)
                for work in range(2):
                    before = (bits, work)
                    after = a.holonomy.move_or_hold(*before, edge)
                    assert a.holonomy.energy_units(*before) == a.holonomy.energy_units(*after)
                    assert (before[0] & region).bit_count() % 2 == (after[0] & region).bit_count() % 2
    result = a.geometry_rotation_covariance()
    assert result['actual_rotation_incidence_columns_checked'] == 60*20*30


def test_preparation_work_changes_full_energy_and_audit_changes_detector_law():
    report = a.preparation_energy()
    assert report['scalar_energy_second_moment_before_birth'] == 56
    assert report['scalar_energy_second_moment_after_birth'] == s.Rational(73, 25)*56
    assert not report['birth_commutes_with_controlled_scalar_H']
    # Independent complete four-dimensional bright energy measurement.
    x = s.Matrix([[0, 1], [1, 0]])
    birth = s.kronecker_product(s.Rational(3, 5)*s.eye(2)-s.I*s.Rational(4, 5)*x, s.eye(2))
    h = s.kronecker_product(s.diag(1, 2), s.sqrt(56)*x)
    initial = s.Matrix([1, 0, 0, 0])
    projects = {(q, sign): s.kronecker_product(s.diag(1, 0) if q == 1 else s.diag(0, 1),
                                             (s.eye(2)+sign*x)/2)
                for q in (1, -1) for sign in (1, -1)}
    work_law = {}
    dephased = s.zeros(4)
    for (q0, sign0), p0 in projects.items():
        e0 = a.strength(q0)*sign0*s.sqrt(56)
        dephased += p0*initial*initial.H*p0
        for (q1, sign1), p1 in projects.items():
            e1 = a.strength(q1)*sign1*s.sqrt(56)
            vector = p1*birth*p0*initial
            probability = s.simplify((vector.H*vector)[0])
            if probability:
                work_law[e1-e0] = work_law.get(e1-e0, 0)+probability
    assert work_law == {0: s.Rational(9, 25), 2*s.sqrt(14): s.Rational(8, 25),
                       -2*s.sqrt(14): s.Rational(8, 25)}
    assert sum(work*weight for work, weight in work_law.items()) == 0
    assert h*birth != birth*h
    assert dephased != initial*initial.H
    # Dephased scalar stays maximally mixed on the bright block at any tick.
    final = birth*dephased*birth.H
    for q in (1, -1):
        charge = s.diag(1, 0) if q == 1 else s.diag(0, 1)
        probability = s.trace(s.kronecker_product(charge, s.diag(0, 1))*final)
        prior = s.trace(s.kronecker_product(charge, s.eye(2))*final)
        assert s.simplify(probability/prior) == s.Rational(1, 2)


def test_occupation_rotation_and_statistics_nonidentification_is_exact():
    assert a.occupation_obstruction()['H_NH_commutator_on_vacuum_norm_squared'] == 56
    result = a.rotation_and_parity_boundaries()
    assert result['native_rotation_Q_commutation_columns'] == 60*132
    assert result['intertwiner_preserving_these_two_rotation_squares_must_vanish']
    assert result['complete_neutral_field_CAR_parity'] == '+I307'
    assert not result['Q_minus_to_neutral_field_odd_parity_intertwiner_exists']


def test_primary_source_paths_resolve():
    root = Path(__file__).resolve().parents[2]
    assert all((root/path).is_file() for path in a.SOURCE_INPUTS)
