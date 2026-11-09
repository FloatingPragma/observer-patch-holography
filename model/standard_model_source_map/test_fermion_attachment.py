"""Independent full-Fock and source-weight checks of the existing CAR crosswalk."""
import itertools
import json

import pytest
import sympy as s

from . import fermion_attachment as f


def creators(modes):
    """Independent tensor-Pauli convention, little-endian occupation order."""
    z = s.diag(1, -1)
    raising = s.Matrix([[0, 0], [1, 0]])
    return tuple(s.kronecker_product(*[raising if j == i else z if j < i else s.eye(2)
                                       for j in reversed(range(modes))])
                 for i in range(modes))


def test_every_three_mode_CAR_column_and_primitive_current_identity():
    cs = creators(3)
    for i, j in itertools.product(range(3), repeat=2):
        assert cs[i]*cs[j]+cs[j]*cs[i] == s.zeros(8)
        assert cs[i].H*cs[j]+cs[j]*cs[i].H == int(i == j)*s.eye(8)
        for k in range(3):
            q = cs[i]*cs[j].H
            assert q*cs[k]-cs[k]*q == int(j == k)*cs[i]
    for i, a in enumerate(cs):
        for mask in range(8):
            column = a[:, mask]
            result = f.creation_action(i, mask, 3)
            if result is None:
                assert column == s.zeros(8, 1)
            else:
                assert column[result[0]] == result[1]
                assert sum(value != 0 for value in column) == 1
        assert a.H*a+a*a.H == s.eye(8)


def test_complex_bilinear_action_matches_independent_full_matrix():
    g = s.Matrix([[2, s.I, s.Rational(1, 3)], [-s.I, -1, 2*s.I],
                  [s.Rational(1, 3), -2*s.I, 3]])
    cs = creators(3)
    q = sum((g[i, j]*cs[i]*cs[j].H for i, j in itertools.product(range(3), repeat=2)), s.zeros(8))
    for mask in range(8):
        actual = f.one_body_action(g, mask)
        assert actual == {i: q[i, mask] for i in range(8) if q[i, mask] != 0}
        for mode in range(3):
            assert f.current_creation_residual(g, mode, mask) == {}


def test_signed_permutation_matches_all_exterior_minors():
    permutation = (2, 4, 0, 3, 1)
    p = s.zeros(5)
    for j, i in enumerate(permutation):
        p[i, j] = 1
    for mask in range(32):
        occupied = [i for i in range(5) if mask >> i & 1]
        image, sign = f.exterior_permutation(mask, permutation)
        targets = [i for i in range(5) if image >> i & 1]
        assert sign == p.extract(targets, occupied).det()
        assert image.bit_count() == mask.bit_count()
    assert f.wedge_census(iter(permutation))['full_Fock_dimension'] == 32


def test_unsigned_swap_passes_one_particle_but_fails_occupied_creation():
    permutation = (1, 0)
    assert f.exterior_permutation(1, permutation) == (2, 1)
    assert f.exterior_permutation(2, permutation) == (1, 1)
    assert f.exterior_permutation(3, permutation) == (3, -1)
    # a0^dagger|01_mode1> has sign+; swapping the resulting two particles
    # has sign-. The unsigned permutation gives the wrong intertwiner.
    created = f.creation_action(0, 2, 2)
    rhs = f.creation_action(1, 1, 2)
    assert created == (3, 1) and rhs == (3, -1)


def test_second_quantized_one_body_Casimir_is_not_the_total_Casimir():
    cs = creators(2)
    paulis = (s.Matrix([[0, 1], [1, 0]])/2,
              s.Matrix([[0, -s.I], [s.I, 0]])/2, s.diag(1, -1)/2)
    charges = [sum((t[i, j]*cs[i]*cs[j].H for i, j in itertools.product(range(2), repeat=2)), s.zeros(4))
               for t in paulis]
    number = sum((a*a.H for a in cs), s.zeros(4))  # creators a: a*a.H=n
    total = sum((q*q for q in charges), s.zeros(4))
    assert total[3, 3] == 0
    assert (s.Rational(3, 4)*number)[3, 3] == s.Rational(3, 2)
    for a in cs:
        assert sum((q*(q*a-a*q)-(q*a-a*q)*q for q in charges), s.zeros(4)) == s.Rational(3, 4)*a


def test_full_add_remove_probabilities_depend_on_the_supplied_input():
    a = creators(1)[0]
    rho = s.diag(s.Rational(1, 3), s.Rational(2, 3))
    assert s.trace(rho*a.H*a) == s.Rational(1, 3)
    assert s.trace(rho*a*a.H) == s.Rational(2, 3)
    assert a.H*a+a*a.H == s.eye(2)


@pytest.fixture(scope='module')
def certificate():
    return f.certificate()


def test_actual_source_intertwiner_and_complete_occupation_census(certificate):
    assert certificate['all_twelve_native_port_currents_intertwined']
    assert certificate['physical_CAR_modes'] == 90
    assert certificate['full_Fock_dimension'] == str(2**90)
    census = certificate['complete_internal_wedge_control']
    assert census['full_Fock_dimension'] == 32768
    assert census['creation_columns_including_zero_columns'] == 491520
    assert certificate['actual_current_occupation_column_controls'] == 720
    assert certificate['ninety_mode_creation_column_controls'] == 270
    assert sorted(certificate['native_to_action_one_particle_permutation']) == list(range(90))
    json.dumps(certificate, allow_nan=False)


def test_exact_target_component_records_cover_each_field_mode_once(certificate):
    rows = certificate['field_operator_charge_records']
    assert len(rows) == 15
    assert sorted(i for row in rows for i in row['native_field_modes']) == list(range(90))
    assert sorted(i for row in rows for i in row['action_field_modes']) == list(range(90))
    assert {row['occupation_mask'] for row in rows} == {3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30}
    assert {row['Q'] for row in rows} == {'-2/3', '1/3', '2/3', '-1/3', '0', '-1', '1'}
    assert all(row['auxiliary_occupation_parity'] == 1 and row['field_creation_parity'] == -1 for row in rows)
    assert all(row['physical_species_assignment'] is None for row in rows)
    for row in rows:
        assert s.Rational(row['Q']) == s.Rational(row['T3'])+s.Rational(row['Y'])


def test_scope_retains_the_supplied_architecture_and_missing_state(certificate):
    assert certificate['three_family_factor_is_declared'] is True
    assert certificate['Casimir_labels_use_creation_operator_double_commutators'] is True
    for flag in ('one_body_Casimir_identified_with_total_many_body_Casimir',
                 'annihilation_identified_with_a_created_antiparticle_without_reference',
                 'auxiliary_even_register_identified_with_odd_field_Fock_space',
                 'configured_CAR_architecture_selected_by_unread_repair',
                 'mass_or_family_eigenbasis_selected', 'default_two_flavour_interacting_ground_used',
                 'additional_native_pulse_tape_executed',
                 'interacting_relativistic_state_or_particle_pole_proved', 'empirical_mass_input_used'):
        assert certificate[flag] is False
    assert certificate['physical_species'] == certificate['physical_masses'] == []


def test_returned_record_mutation_cannot_corrupt_source_weights(certificate):
    rows = f.casimir_and_weight_records()
    rows[0]['color_weight'][0] = '1000'
    rows[0]['native_field_modes'][0] = 999
    assert f.casimir_and_weight_records() == certificate['field_operator_charge_records']


@pytest.mark.parametrize('args', [(True, 0, 0), (0, 15, 0), (0, 0, 3), (-1, 0, 0)])
def test_invalid_physical_mode_indices_rejected(args):
    with pytest.raises((TypeError, ValueError)):
        f.mode_index(*args)


@pytest.mark.parametrize('permutation', [(), (0, 0), (True, 1), (0, 2)])
def test_incomplete_or_inexact_permutation_rejected(permutation):
    with pytest.raises(ValueError):
        f.exterior_permutation(0, permutation)


def test_no_dense_ninety_mode_census_or_floating_generator_admitted():
    with pytest.raises(ValueError, match='sixteen'):
        f.wedge_census(range(90))
    with pytest.raises(ValueError, match='exact finite'):
        f.one_body_action(s.Matrix([[1.0]]), 0)
    with pytest.raises(ValueError, match='occupation mask'):
        f.creation_action(0, 8, 3)


def test_mutated_actual_current_fails_its_literal_fock_commutator_dictionary(monkeypatch):
    original = f.one_body_action
    def changed(generator, mask):
        out = original(generator, mask)
        if mask == 1:
            out[1] = out.get(1, 0)+1
        return out
    monkeypatch.setattr(f, 'one_body_action', changed)
    assert f.current_creation_residual(s.eye(3), 1, 1) != {}
