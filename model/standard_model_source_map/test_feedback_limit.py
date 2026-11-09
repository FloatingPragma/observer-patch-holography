from fractions import Fraction as F
from itertools import product
import json

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from model.common_source import native_family_record_memory as memory
from model.common_source import native_gauge_family_records as gates
from model.common_source import native_record_yukawa_field as field
from model.common_source import native_repair_mass as repair
from . import feedback_limit as limit


@pytest.mark.parametrize('policy,trace', [('always', F(1)), ('C', F(1, 12)),
                                        ('E', F(1, 6)), ('CE', F(1, 72))])
def test_actual_raw_records_give_the_unconditional_mean(policy, trace):
    data = memory.architecture()
    expected = s.zeros(3)
    total = [s.Integer(0)]*12
    for old in range(12):
        for cc, cr, ec, er, c, weight in repair.native_records():
            e = ec == 0 and er in (4, 5)
            dest = memory._destination(old, c, e, 2, data['gC'], data['gE'])
            total[old] += weight
            if limit.enabled(policy, c, e):
                expected += data['stationary'][old]*weight*data['projectors'][dest]
    observed = limit.record_model(F(1, 12), F(1, 6), policy)
    assert total == [1]*12
    assert all(s.simplify(x) == 0 for x in expected-observed['mean'])
    assert s.simplify(s.trace(expected)-s.Rational(trace.numerator, trace.denominator)) == 0
    assert len(observed['rows']) == 48
    assert observed['transition']*observed['stationary'] == observed['stationary']


def test_C_duty_is_not_renormalized_or_multiplied_mean():
    p, q = F(1, 12), F(1, 6)
    a = limit.record_model(p, q, 'always')
    c = limit.record_model(p, q, 'C')
    assert any(s.simplify(x) != 0 for x in c['mean']-s.Rational(1, 12)*a['mean'])
    assert all(s.simplify(x) == 0 for x in c['mean']-a['mean']+s.Rational(11, 12)*a['projectors'][2])


@pytest.mark.parametrize('p,q', [(F(0), F(0)), (F(0), F(1)), (F(1, 3), F(0)), (F(1, 3), F(1))])
def test_endpoint_zero_records_remain_complete(p, q):
    for duty in limit.POLICIES:
        model = limit.record_model(p, q, duty)
        assert any(row['probability'] == 0 for row in model['rows'])
        for old in range(12):
            assert sum(row['probability'] for row in model['rows'] if row['old'] == old) == 1
        paths = limit.path_instrument(p, q, duty, steps=2, time=F(1, 32), initial_memory=0)
        assert len(paths) == 16 and any(row['probability'] == 0 for row in paths)
        gram = sum((row['kraus'].conj().T@row['kraus'] for row in paths), np.zeros((6, 6), complex))
        np.testing.assert_allclose(gram, np.eye(6), atol=2e-14, rtol=0)


@pytest.mark.parametrize('duty', limit.POLICIES)
def test_complete_reference_entangled_paths_with_noncommuting_common_clock(duty):
    h0 = np.diag([20, -7, 3, 5, -11, 9]).astype(complex)
    h0[0, 4] = 2j
    h0[4, 0] = -2j
    psi = np.arange(1, 13)+1j*np.arange(12, 0, -1)
    psi = psi/np.linalg.norm(psi)
    time, p, q = F(1, 16), F(1, 12), F(1, 6)
    bound = limit.averaging_bound(3, time, p)
    actual = np.zeros((12, 12), complex)
    ideal = np.zeros_like(actual)
    archive = 0.
    for old, weight in ((0, .4), (5, .6)):
        local = psi if old == 0 else np.roll(psi, 3)
        paths = limit.path_instrument(p, q, duty, steps=3, time=time,
                                      initial_memory=old, common_hamiltonian=h0)
        assert len(paths) == 64
        assert sum(r['probability'] for r in paths) == 1
        gram = sum((r['kraus'].conj().T@r['kraus'] for r in paths), np.zeros((6, 6), complex))
        np.testing.assert_allclose(gram, np.eye(6), atol=2e-14, rtol=0)
        for row in paths:
            v = np.kron(row['kraus'], np.eye(2))@local
            w = np.kron(row['ideal_kraus'], np.eye(2))@local
            a, b = np.outer(v, v.conj()), np.outer(w, w.conj())
            actual += weight*a
            ideal += weight*b
            archive += weight*np.abs(np.linalg.eigvalsh(a-b)).sum()
    reduced = np.abs(np.linalg.eigvalsh(actual-ideal)).sum()
    assert 1e-8 < archive <= float(bound['archive_distance_upper'])+1e-12
    assert reduced <= float(bound['reduced_distance_upper'])+1e-12
    assert reduced <= archive+1e-12


def test_stationary_product_reduces_only_the_reduced_bound():
    arbitrary = limit.averaging_bound(400, 1, F(1, 12))
    stationary = limit.averaging_bound(400, 1, F(1, 12), stationary_product=True)
    assert stationary['reduced_distance_upper'] < arbitrary['reduced_distance_upper']
    assert stationary['archive_distance_upper'] == arbitrary['archive_distance_upper']
    assert limit.averaging_bound(10, 0, F(1, 12))['archive_distance_upper'] == 0
    assert limit.averaging_bound(10, 1, F(1, 12), feedback_norm=0)['archive_distance_upper'] == 0


def test_noncommuting_reset_pulses_are_old_memory_independent():
    for duty in limit.POLICIES:
        data = limit.record_model(F(1, 12), F(1, 6), duty)
        for e in (False, True):
            rows = [r for r in data['rows'] if not r['C'] and r['E'] == e]
            assert len({(r['destination'], r['enabled']) for r in rows}) == 1
    # A controller using the old memory on reset violates the proof premise.
    P = [np.asarray(x, complex) for x in data['projectors']]
    eta = np.zeros((3, 3), complex)
    eta[0, 0] = 1
    assert np.linalg.norm(expm(-.1j*P[2])@eta@expm(.1j*P[2])
                          -expm(-.1j*P[3])@eta@expm(.1j*P[3])) > 1e-4


def test_exact_archive_budget_and_outward_square_root():
    for p, time, norm, error in product((F(0), F(1, 12), F(9, 10)),
                                       (F(0), F(1, 3)), (F(1), F(14)), (F(1, 16), F(1, 3))):
        n = limit.steps_for_archive(error, time, p, feedback_norm=norm)
        bound = limit.averaging_bound(n, time, p, feedback_norm=norm)
        assert bound['archive_distance_squared_upper'] <= error*error
        assert bound['archive_distance_upper']**2 >= bound['archive_distance_squared_upper']
    assert limit.sqrt_upper(F(2))**2 >= 2
    assert limit.sqrt_upper(F(9, 16)) == F(3, 4)


def _creator(modes, site):
    a = np.zeros((1 << modes, 1 << modes), complex)
    for word in range(1 << modes):
        if not word & (1 << site):
            a[word | (1 << site), word] = (-1)**((word & ((1 << site)-1)).bit_count())
    return a


def test_pairing_norm_on_all_four_mode_Fock_rows():
    K = np.array([[1+2j, 3-1j], [3-1j, -2+.5j]])
    delta = np.kron([[0, 1], [-1, 0]], K.conj())
    creators = [_creator(4, i) for i in range(4)]
    V = sum((delta[i, j]*creators[i]@creators[j] for i in range(4) for j in range(i+1, 4)),
            np.zeros((16, 16), complex))
    V += V.conj().T
    assert abs(np.linalg.norm(V, 2)-np.linalg.svd(K, compute_uv=False).sum()) < 2e-14
    assert np.linalg.norm(V, 2) > np.linalg.norm(K, 2)


def test_actual_Higgs_map_lifts_all_three_source_means_and_rank_norms():
    P = np.asarray(memory.architecture()['projectors'][2], complex)
    norms = []
    for gate in gates.higgs_supports():
        channel = {'Y': P, 'Higgs_components': np.asarray(gate['components'], complex),
                   'left_projector': np.asarray(gate['left_projector'], complex),
                   'right_projector': np.asarray(gate['right_projector'], complex)}
        K = field.assemble_kernel({'sector': channel})
        lift = limit.pairing_lift(K)
        assert np.max(abs(lift['creation_kernel']+lift['creation_kernel'].T)) < 1e-14
        norm = 2*gate['right_dimension']
        assert abs(lift['Fock_operator_norm_diagnostic']-norm) < 2e-13
        norms.append(norm)
        for duty in limit.POLICIES:
            data = limit.record_model(F(1, 12), F(1, 6), duty)
            direct = field.assemble_kernel({'sector': dict(channel, Y=np.asarray(data['mean'], complex))})
            mean = np.zeros((45, 45), complex)
            for row in data['rows']:
                if row['enabled']:
                    coeff = float(data['stationary'][row['old']])*float(row['probability'])
                    pp = np.asarray(data['projectors'][row['destination']], complex)
                    mean += coeff*field.assemble_kernel({'sector': dict(channel, Y=pp)})
            np.testing.assert_allclose(mean, direct, atol=2e-14, rtol=0)
    assert sorted(norms) == [2, 6, 6] and sum(norms) == 14


def test_finite_volume_and_attempt_cost_are_not_silently_uniform():
    one = limit.physical_budget(10000, 1, (F(1, 12),)*3)
    two = limit.physical_budget(10000, 1, (F(1, 12),)*3, cells=2)
    assert one['pairing_norm_upper'] == 14
    assert one['joint_C_probability'] == 1-F(11, 12)**3
    assert two['pairing_norm_upper'] == 28
    assert two['uniform_fixed_memory_bound'] > one['uniform_fixed_memory_bound']
    assert two['source_volume_uniform'] is False


@pytest.mark.parametrize('args', [(0, 1, F(1, 12)), (True, 1, F(1, 12)),
                                  (1, -1, 0), (1, 1., 0), (1, 1, 1), (1, 1, -1)])
def test_invalid_bounds_are_rejected(args):
    with pytest.raises(ValueError):
        limit.averaging_bound(*args)


@pytest.mark.parametrize('p,q,policy', [(1, 0, 'C'), (0, 2, 'C'), (0., 0, 'C'),
                                      (0, 0, 'selected_success')])
def test_invalid_event_laws_and_postselected_policy_are_rejected(p, q, policy):
    with pytest.raises(ValueError):
        limit.record_model(p, q, policy)


def test_certificate_retains_scopes_and_is_detached_compact_json():
    value = limit.certificate()
    encoded = json.dumps(value, sort_keys=True)
    assert len(encoded) < 30000
    assert value['finite_step_channel_is_exact_mean_unitary'] is False
    assert value['physical_feedback_tape_executed'] is False
    assert value['flavor_singular_values_modified_by_limit'] is False
    histories = value['complete_finite_history_control']
    assert histories['total_complete_coarse_histories'] == 256
    assert len(histories['original_path_probabilities_exact']) == 64
    assert sum(map(F, histories['original_path_probabilities_exact'])) == 1
    assert histories['native_fine_tape_executed'] is False
    value['all_four_duties'][0]['mean_exact'][0][0] = 'corrupted'
    assert json.dumps(limit.certificate(), sort_keys=True) == encoded
    argument = F(1, 12)
    result = limit.averaging_bound(4096, 1, argument)
    result['uniform_fixed_memory_bound']._numerator = 1
    assert argument == F(1, 12)
