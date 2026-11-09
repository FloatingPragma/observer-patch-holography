"""Exact source, inference and scope controls for short-record reconstruction."""
from fractions import Fraction as F

import pytest
import sympy as s

from model.native_record_response import short_records as r, response, transfer
from model.observer import native


@pytest.fixture(scope='module')
def model():
    return r._model()


def test_full_raw_minimality_and_two_record_inverse(model):
    _, features, X, A, _, G, Gi = model
    assert r.one_record_rank() == 42
    assert X*A == A*X == G*Gi == s.eye(66)
    assert [sum(f[0] == kind for f in features) for kind in ('constant', 'single', 'pair')] == [1, 11, 54]
    assert 1/r._norminf(Gi) == r.GRAM_FLOOR
    assert r._norm1(A) == F(7955, 396)
    assert r._norminf(A) == F(464, 11)


def test_two_record_features_against_independent_literal_endpoints(model):
    d, features, X, *_ = model
    weights = (25, 25)+(8,)*10
    for column, feature in enumerate(features):
        for i, state in enumerate(d['states']):
            if feature[0] == 'constant':
                expected = 1
            else:
                e = response.local.SEAMS[feature[1]]
                contrast = state[e[0]]-state[e[1]]
                if feature[0] == 'single':
                    expected = contrast
                else:
                    f = response.local.SEAMS[feature[2]]
                    expected = F()
                    for coin in (0, 1):
                        following, _ = native.native(state, e, coin)
                        weight = F(weights[e[0] if coin else e[1]], weights[e[0]]+weights[e[1]])
                        expected += weight*contrast*(following[f[0]]-following[f[1]])
            assert X[i, column] == expected
    # A repeated passive product omits the first measurement's disturbance.
    assert any(X[:, j] != s.Matrix([d['D'][f[1], i]*d['D'][f[2], i] for i in range(66)])
               for j, f in enumerate(features) if f[0] == 'pair')


def test_reflected_five_event_entries_from_raw_source_insertions(model):
    d, features, _, _, _, _, _ = model
    chosen = next(row for row in d['branches'] if row.law == 'E' and row.current == 1)
    label = transfer.label(chosen)
    matrix = r.moment_matrix(label)
    def weighted(vector, feature, reverse=False):
        for e in (feature[1:] if not reverse else tuple(reversed(feature[1:]))):
            out = s.zeros(66, 1)
            for row in d['branches']:
                raw = transfer.reverse_label(transfer.label(row)) if reverse else transfer.label(row)
                score_row = r._branch(raw)
                score = score_row.feature(response.local.SEAMS[e])
                if score:
                    for i, j in row.pairs:
                        out[j] += row.weight*score*vector[i]
            vector = out
        return vector
    for i, j in ((12, 17), (36, 65), (0, 62), (1, 4)):
        state = weighted(d['pi'], features[i], reverse=True)
        state = r._backward(chosen).T*state
        state = weighted(state, features[j])
        assert sum(state) == matrix[i, j]


def test_all_native_selective_maps_and_complete_held_out_law():
    c = r.certificate()
    assert c['exact_all_selective_map_inverse_checks'] == 721
    assert c['zero_raw_slots_retained'] == 421
    probabilities = list(map(F, c['complete_seven_event_sign_law'].values()))
    assert len(probabilities) == 128 and min(probabilities) > 0 and sum(probabilities) == 1
    raw = [F(row['probability']) for row in c['held_out_original_word_probabilities']]
    assert raw[0] > 0 and raw[1] > 0 and raw[2] == 0


def test_held_out_raw_words_against_literal_interpreter(model):
    d = model[0]
    for word in r._positive_words():
        program = response.source()
        density = {(i, 0, i, 0): F(p) for i, p in enumerate(d['pi'])}
        memory = archive = ()
        for label in word:
            outputs = response.interpreter._outcomes(program, program['commands'][0], density, memory, archive)
            assert len(outputs) == 721
            selected = [out for out in outputs if (out.record.kind, out.record.law,
                        out.record.seam, out.record.coin, out.record.raw) == label]
            outcome, = selected
            density, memory, archive = outcome.density, outcome.memory, outcome.archive
        assert response.interpreter.trace(density) == r.word_probability(word)


def test_data_mutation_changes_recovered_original_probability(model):
    row = next(row for row in model[0]['branches'] if row.pairs)
    moment = s.MutableDenseMatrix(r.moment_matrix(transfer.label(row)))
    good = r.recover_insertion(moment)
    moment[0, 0] += s.Rational(1, 1000)
    changed = r.recover_insertion(moment)
    G = r.gram(); e0 = s.eye(66)[:, 0]
    assert (e0.T*G*(changed-good)*e0)[0] == s.Rational(1, 1000)


def test_zero_flags_and_unconditioned_scores(model):
    labels = tuple(transfer.label(row) for row in model[0]['branches'])
    empty = next(transfer.label(row) for row in model[0]['branches'] if not row.pairs)
    for feature in r.feature_labels():
        assert abs(r.score(feature, (labels[1], labels[2]))) <= r.SCORE_BOUND
        assert r.score(feature, (empty, empty)) == int(feature[0] == 'constant')
    assert r.word_probability((empty,)) == 0
    assert r.word_probability(()) == 1


def test_complete_prediction_budget_has_original_probability_and_preparation(model):
    epsilon = F(1, 2**60)
    budget = r.noise_budget(epsilon, 16, preparation_error=F(1, 10**6))
    C = F(1248314867975, 264)
    assert F(budget['complete_kernel_TV_lipschitz_factor']) == C
    assert F(budget['complete_history_TV_upper']) == F(1, 10**6)+16*C*epsilon
    assert r.noise_budget(F(1), 1)['complete_history_TV_upper'] == '1'
    assert r.noise_budget(F(0), 2)['complete_history_TV_upper'] == '0'


def test_Gram_and_log_error_are_rational_admissions():
    budget = r.transfer_energy_error(F(1, 2**40))
    delta = F(66, 2**40)
    eta = 2*delta/(r.GRAM_FLOOR-delta)
    assert F(budget['transfer_eigenvalue_error_upper']) == eta
    assert F(budget['transfer_log_energy_error_upper']) == eta/(transfer.LOWER-eta)
    with pytest.raises(ValueError, match='positive Gram'):
        r.transfer_energy_error(r.GRAM_FLOOR/66)
    with pytest.raises(ValueError, match='positive transfer'):
        r.transfer_energy_error(r.GRAM_FLOOR/132)


def test_sampler_plan_is_sufficient_and_not_executed():
    e, alpha = F(1, 2**40), F(1, 100)
    plan = r.sample_budget(e, alpha)
    N = plan['independent_stationary_trials_per_menu_sufficient']
    q = plan['simultaneous_scalar_quantities']
    L = 0
    while 2**L < 2*q/alpha:
        L += 1
    assert N*e**2 >= 2*r.SCORE_BOUND**4*L
    assert N > 10**40 and plan['sampled_trials_executed'] is False


@pytest.mark.parametrize('call', [lambda:r.noise_budget(-1, 2), lambda:r.noise_budget(F(0), True),
    lambda:r.noise_budget(.1, 2), lambda:r.noise_budget(0, 2, preparation_error=2),
    lambda:r.sample_budget(0, F(1, 2)), lambda:r.sample_budget(F(1, 10), 1),
    lambda:r.transfer_energy_error(True), lambda:r.coarse_history(9),
    lambda:r.gram(2), lambda:r.word_probability([]), lambda:r.moment_matrix(('unknown',)),
    lambda:r.recover_insertion(s.eye(2)), lambda:r.recover_insertion(s.eye(66)*.1)])
def test_input_domains_are_checked(call):
    with pytest.raises(ValueError):
        call()


def test_public_certificate_is_detached_and_retains_scope():
    value = r.certificate()
    value['feature_basis'][0][0] = 'changed'
    value['scope']['full_quantum_instrument_tomography'] = True
    fresh = r.certificate()
    assert fresh['feature_basis'][0] == ['constant']
    assert fresh['scope']['full_quantum_instrument_tomography'] is False
    assert fresh['scope']['physical_clock_or_particle_energy_identified'] is False
    assert fresh['scope']['coherent_CAR_gauge_Higgs_action_selected'] is False


@pytest.mark.parametrize('mutation', ('physical_clock_or_particle_energy_identified',
    'full_quantum_instrument_tomography', 'coherent_CAR_gauge_Higgs_action_selected',
    'empirical_sampling_executed'))
def test_owner_rejects_scope_forgery(tmp_path, monkeypatch, mutation):
    import json
    actual = {'short_record_population_realization': r.certificate()}
    forged = json.loads(json.dumps(actual))
    forged['short_record_population_realization']['scope'][mutation] = True
    path = tmp_path/'forged.json'
    path.write_text(json.dumps(forged))
    monkeypatch.setattr(response, 'report', lambda: actual)
    with pytest.raises(ValueError, match='receipt or source bindings differ'):
        response.check(path)


def test_owner_binds_all_short_record_proof_code_and_Lean():
    bindings = response.source_bindings()
    assert all(path in bindings for path in r.SOURCE_FILES)


def test_actual_mixing_martingale_and_program_cap():
    b = r.continuing_calibration_budget()
    assert F(b['score_conditional_mean_bias_upper']) <= F(1, 2**41)
    assert b['burn_in_attempts_per_trial'] == 8192
    assert b['maximum_complete_trials_in_one_admitted_program'] == (2**32)//8197
    assert b['raw_attempts_for_two_menus_upper'] > 2**32
    assert b['required_sampler_fits_one_admitted_program'] is False
    assert b['stationary_independent_resets_assumed'] is False
    assert b['full_source_and_new_burn_in_archive_product_reset_claimed'] is False
    assert b['martingale_trials_per_menu_sufficient'] == r.sample_budget(F(1, 2**41), F(1, 100))['independent_stationary_trials_per_menu_sufficient']
    with pytest.raises(ValueError, match='fixed native burn-in'):
        r.continuing_calibration_budget(F(1, 2**120), F(1, 100))


def test_actual_complete_interpreter_resource_stop():
    result = r.resource_stop_control()
    assert result['attempts_executed'] == 1
    assert result['retained_terminal_blocks'] == 721
    assert result['total_original_probability'] == '1'
    assert result['completed_calibration_trials'] == 0
