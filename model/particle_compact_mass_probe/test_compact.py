"""Exact blank maps, distinct failures, complete references and resource checks."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s

from model.particle_common_contract import contract
from model.particle_filtered_probe import filter as apparatus
from . import compact as c


@pytest.fixture(scope='module')
def control():
    return c.quantum_control()


@pytest.fixture(scope='module')
def member():
    return contract.member(1)


def array(value):
    return np.array(value.evalf(18).tolist(), dtype=complex)


def test_exact_direct_dilation_is_unitary_and_uses_the_same_blank(control):
    direct = control['compact']
    assert direct['success'] == control['original']['success']
    assert direct['unitary'].H*direct['unitary'] == s.eye(4)
    assert direct['kraus'][1] == s.Rational(4, 5)*s.eye(2)
    assert direct['selector_register_qubits'] == 0


@pytest.mark.parametrize('imaginary', [False, True])
def test_success_ramsey_maps_and_original_scores_equal_for_both_phases(control, imaginary):
    old = apparatus.ramsey_instrument(control['H'], control['original']['kraus'], '1/2', 1, imaginary)
    direct = c.compact_ramsey(control['H'], control['B'], '1/2', 1, imaginary)['rows']
    assert len(old) == 8 and len(direct) == 4
    assert [x['K'] for x in old[:2]] == [x['K'] for x in direct[:2]]
    probabilities = [s.simplify(s.trace(control['rho']*row['K'].H*row['K'])) for row in direct]
    assert sum(probabilities) == 1 and min(probabilities) >= 0
    assert sum(probabilities[:2]) == s.Rational(9, 25)
    assert probabilities[0]-probabilities[1] == (s.Rational(9, 50) if imaginary else 0)


def test_failure_channels_and_failure_ramsey_records_are_different(control):
    assert control['old_failure_on_zero'] == s.diag(0, s.Rational(16, 25))
    assert control['compact_failure_on_zero'] == s.diag(s.Rational(16, 25), 0)
    real = control['receipt']['phases'][0]
    assert real['compact_probabilities'][2:] == ['16/25', '0']
    assert real['original_probabilities'][4:6] == ['8/25', '8/25']


@pytest.mark.parametrize('imaginary', [False, True])
def test_full_four_outcome_state_retains_the_entangled_idle_reference(control, imaginary):
    psi = np.diag([np.sqrt(3)/2, 1/2]).reshape(4)
    rows = c.compact_ramsey(control['H'], control['B'], '1/2', 1, imaginary)['rows']
    outputs = [np.kron(array(row['K']), np.eye(2))@psi for row in rows]
    assert abs(sum(np.vdot(x, x).real for x in outputs)-1) < 1e-14
    reference = sum(x.reshape(2, 2).T@x.reshape(2, 2).conj() for x in outputs)
    assert np.linalg.norm(reference-np.diag([3/4, 1/4])) < 1e-14
    old = apparatus.ramsey_instrument(control['H'], control['original']['kraus'], '1/2', 1, imaginary)
    for index in (0, 1):
        old_output = np.kron(array(old[index]['K']), np.eye(2))@psi
        assert np.linalg.norm(old_output-outputs[index]) < 1e-14


def test_zero_probe_and_involution_keep_all_four_records():
    h = s.diag(0, 1)
    for b in (s.zeros(2), s.Matrix([[0, 1], [1, 0]])):
        direct = c.compact_ramsey(h, b)
        assert len(direct['rows']) == 4
        total = sum((row['K'].H*row['K'] for row in direct['rows']), s.zeros(2))
        assert total.applyfunc(s.simplify) == s.eye(2)
    zero = c.compact_instrument(s.zeros(2))
    assert zero['kraus'] == (s.zeros(2), s.eye(2))


def test_complete_operator_error_covers_a_perturbed_entangled_instrument():
    h = s.diag(0, 4); x = s.Matrix([[0, 1], [1, 0]])
    old = c.compact_ramsey(h, 3*x/5, '1/2', 1)['rows']
    perturbed = c.compact_ramsey(h, 61*x/100, '1/2', 1)['rows']
    psi = np.array([1, 0, 0, 1], dtype=complex)/np.sqrt(2)
    distance = 0
    for left, right in zip(old, perturbed, strict=True):
        u = np.kron(array(left['K']), np.eye(2))@psi
        v = np.kron(array(right['K']), np.eye(2))@psi
        difference = np.outer(u, u.conj())-np.outer(v, v.conj())
        distance += np.abs(np.linalg.eigvalsh(difference)).sum()/2
    bound = c.complete_error_bound('1/100')
    assert distance <= float(bound['complete_halftrace_error_upper'])
    assert bound['complete_halftrace_error_upper'] < 1
    assert c.complete_error_bound()['complete_halftrace_error_upper'] == 0


def test_factored_four_outcome_compiler_matches_inherited_finite_formula():
    direct = c.conditional_compiler(3, F(1, 2), 4)
    inherited = apparatus.factored_compiler_plan(3, 4, F(1, 512))
    for key in ('padded_source_bits', 'retained_record_bits', 'padded_unitary_bits',
                'pulse_bound_exact_expression', 'pulse_bound_factored_upper',
                'port_fraction_bits', 'rank_one_tree_angle_bits'):
        assert direct[key] == inherited[key]
    sigma = direct['complete_halfdiamond_allowance']
    assert F(sigma['coefficient'])/2**sigma['negative_exponent'] == F(1, 512)


def test_primary_plan_reuses_success_weight_and_source_with_small_quantum_carrier(member):
    value = c.common_contract_bound(member)
    c.verify_bound(member, value)
    result = value['result']; budget = result['carrier_budget']
    assert result['same_blank_Ramsey_Kraus_operators_on_all_inputs'] is True
    assert result['complete_failure_channel_equal_to_selector_apparatus'] is False
    assert budget['padded_source_reference_and_apparatus_qubit_upper'] == 116738
    assert budget['one_shot_joint_capacity_upper']['exact_log_terms'] == [
        {'base': 2, 'coefficient': 61442}, {'base': 451, 'coefficient': 6144}]
    assert budget['total_self_contained_executable_capacity_proved'] is False
    assert budget['conditional_compiler_extra_service_flags_included'] is False
    assert result['quadrature_nodes_in_operator_definition']['exponent'] > 10**12
    assert result['conditional_compiler']['retained_record_bits'] == 2
    assert result['conditional_compiler']['port_fraction_bits'] > 10**12
    bad = deepcopy(value)
    bad['result']['complete_failure_channel_equal_to_selector_apparatus'] = True
    with pytest.raises(ValueError):
        c.verify_bound(member, bad)


@pytest.mark.parametrize('probe', [s.Matrix([[0, 1], [0, 0]]), s.diag(2, 0)])
def test_invalid_direct_dilation_targets_are_rejected(probe):
    with pytest.raises(ValueError):
        c.compact_instrument(probe)
