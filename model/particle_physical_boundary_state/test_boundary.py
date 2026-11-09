from copy import deepcopy
from fractions import Fraction as F
import json

import pytest
import sympy as s

from model.particle_common_contract import contract
from . import boundary


def test_full_word_charge_spectrum_has_constructive_endpoint_bounds():
    assert boundary.charge_spectrum(0) == (0,)
    assert boundary.charge_spectrum(1) == (-6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6)
    assert boundary.charge_spectrum(2) == (-12, *range(-10, 11), 12)
    for degree in range(1, 5):
        values = set(boundary.charge_spectrum(degree))
        assert set(boundary.charge_spectrum(degree-1)) <= values
        assert values == {-q for q in values}
        assert max(values) == 6*degree
        assert 6*degree-1 not in values


def test_entire_native_v1_retains_all_charge_multiplicities():
    control = boundary.native_degree_one_control()
    assert control['dimension'] == 119
    assert control['outcome_ranks'] == {
        '-6': 1, '-4': 9, '-3': 4, '-2': 9, '-1': 36, '0': 1,
        '1': 36, '2': 9, '3': 4, '4': 9, '6': 1}
    q = control['matrix']
    assert sum((p for _, p in control['projectors']), s.zeros(119)) == s.eye(119)
    assert sum((value*p for value, p in control['projectors']), s.zeros(119)) == q


def test_charge_functional_calculus_resolves_non_diagonal_source_and_all_records():
    rotation = s.Matrix([[s.Rational(3, 5), -s.Rational(4, 5)],
                         [s.Rational(4, 5), s.Rational(3, 5)]])
    q = rotation*s.diag(-3, 6)*rotation.T
    pvm = boundary.charge_projectors(q, (-3, 0, 6))
    assert pvm[1][1] == s.zeros(2)  # zero outcomes are retained
    assert pvm[0][1] == rotation[:, 0]*rotation[:, 0].T
    assert pvm[2][1] == rotation[:, 1]*rotation[:, 1].T
    record_isometry = s.Matrix.vstack(*(p for _, p in pvm))
    assert record_isometry.H*record_isometry == s.eye(2)


@pytest.mark.parametrize('matrix,spectrum', [
    (s.diag(-3, 6), (-3, 0)),
    (s.Matrix([[0, 1], [0, 0]]), (0, 1)),
    (s.diag(0, 1), (0, 0, 1)),
    (s.diag(0, s.Float(1)), (0, 1)),
])
def test_incomplete_or_nonexact_charge_measurement_is_rejected(matrix, spectrum):
    with pytest.raises(ValueError):
        boundary.charge_projectors(matrix, spectrum)


def test_native_character_blocks_have_scalar_commutants():
    rows = boundary.source_irreducibility_control()
    assert {(row['q6Y'], row['dimension']) for row in rows} == {
        (-4, 3), (-3, 2), (1, 6), (2, 3), (6, 1)}
    assert all(row['commutant_dimension'] == 1 for row in rows)


def test_loop_schmidt_contraction_and_eleven_distinct_physical_charge_tuples():
    rows = boundary.loop_witnesses()
    assert len(rows) == len({row['cut_charges'] for row in rows}) == 11
    for row in rows:
        d = row['irrep_dimension']
        # Haar-normalized path coefficients: chi(AB)=sum_ab A_ab B_ba.
        coefficients = s.zeros(d*d)
        for a in range(d):
            for b in range(d):
                coefficients[a*d+b, b*d+a] = s.Rational(1, d)
        marginal = coefficients*coefficients.H
        assert s.trace(marginal) == 1
        assert marginal == s.eye(d*d)*s.Rational(1, d*d)
        assert row['Schmidt_rank'] == marginal.rank()
        assert row['regional_Schmidt_eigenvalue'] == F(1, d*d)
        assert sum(row['cut_charges']) == 0


def test_selected_state_has_actual_boundary_weights_and_complete_menu():
    member = contract.member(0)
    payload = boundary.common_contract_payload(member)
    menu = payload['complete_central_charge_measurement']
    assert menu['outcome_count'] == 1331
    assert menu['effect_sum'] == 'identity on the entire regional carrier'
    assert menu['all_zero_probability_outcomes_retained']
    assert not menu['all_nonabelian_center_labels_resolved']
    assert payload['state']['conjugate_reference_idle_H'] == 0
    assert payload['state']['beta'] == member['state']['beta']
    assert payload['state']['dyadic_scale'] == member['state']['dyadic_scale']
    assert payload['finite_physical_preparation_including_boundary_charges_constructed']
    assert not payload['refinement_compatible_state_family_proved']
    assert not payload['full_spatial_density_or_gate_tape_executed']
    floor = payload['physical_density_floor_factored']['negative_exponent']
    assert floor > 0
    for row in payload['strictly_positive_selected_state_witnesses']:
        assert row['selected_state_sector_probability_lower']['negative_exponent'] == floor
        assert F(row['regional_Schmidt_span_floor']['coefficient']) == F(1, row['Schmidt_rank'])


def test_exact_consumer_rejects_forged_source_and_result():
    member = contract.member(0)
    receipt = boundary.bind_member(member)
    boundary.verify_bound(member, receipt)
    bad_member = deepcopy(member)
    bad_member['reference']['idle_H'] = '1'
    with pytest.raises(ValueError):
        boundary.common_contract_payload(bad_member)
    bad_receipt = deepcopy(receipt)
    bad_receipt['result']['refinement_compatible_state_family_proved'] = True
    bad_receipt['result_sha256'] = contract.digest(bad_receipt['result'])
    with pytest.raises(ValueError):
        boundary.verify_bound(member, bad_receipt)


def test_finite_next_member_keeps_boundary_domains_separate():
    result = boundary.common_contract_payload(contract.member(1))
    assert result['complete_central_charge_measurement']['outcome_count'] == 23**3
    assert not result['full_spatial_multiplicities_computed']
    assert not result['interacting_relativistic_limit_or_particle_identification']


def test_certificate_is_json_safe_without_dense_spatial_output():
    certificate = boundary.certificate()
    encoded = json.dumps(certificate, allow_nan=False)
    assert len(encoded) < 20000
    assert certificate['boundary_preparation_requirement_closed_at_finite_member']
    assert not certificate['selected_full_spatial_run_executed']


@pytest.mark.parametrize('degree', [-1, True, F(1, 2)])
def test_invalid_word_degree_rejected(degree):
    with pytest.raises(ValueError):
        boundary.charge_spectrum(degree)
