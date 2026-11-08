"""Independent spectral, full-record, intervention and admission controls."""
from fractions import Fraction as F
from itertools import combinations, product
import json
from math import factorial

import numpy as np
import pytest

from . import response as r, packet as p


@pytest.fixture(scope='module')
def source():
    family, first, second, h = p.local_case()
    states = tuple(combinations(range(12), 2)); index = {s: i for i, s in enumerate(states)}
    transfer = np.zeros((66, 66))
    for j, state in enumerate(states):
        mask = sum(1 << i for i in state)
        for edge, coin in product(family.seams, (0, 1)):
            target, _, _ = r.exchange.sector.source.image(mask, edge, coin, 12)
            pair = tuple(i for i in range(12) if target >> i & 1)
            transfer[index[pair], j] += 1/60
    eigenvalues, vectors = np.linalg.eigh(transfer)
    assert eigenvalues.min() >= 2/3-1e-14
    probe = np.array([h.get(s, 0) for s in states], dtype=float)
    return family, first, second, h, transfer, eigenvalues, vectors, probe


@pytest.mark.parametrize('time', (F(0), F(1, 64), F(-1, 64), F(1, 8)))
def test_exact_coefficient_polynomial_matches_native_ground_response(source, time):
    family, _, _, h, transfer, vals, vectors, probe = source
    row = r.continuation(30, 5, 2, time, 12)
    lags = r.exact_native_lags(family, h, 2, row['lag_attempts'])
    actual = r.reconstruct(row['power_coefficients'], lags)
    estimate = complex(float(actual[0]), float(actual[1]))
    exact = ((vectors.T@probe)**2/4) @ np.exp(1j*60*float(time)*np.log(vals))
    assert abs(estimate-exact) <= float(row['operator_error_upper'])+1e-12
    for lag, response in zip(row['lag_attempts'], lags, strict=True):
        assert abs(float(response)-probe@np.linalg.matrix_power(transfer, lag)@probe/4) < 1e-13
    assert r.reconstruct(row['power_coefficients'], (F(1),)*13) == (F(1), F(0))


def test_blocking_spectral_floor_and_continuation_bound_on_independent_scalar_grid():
    for seams in (30, 630, 2520, 10080):
        degree = 5 if seams == 30 else 6
        row = r.continuation(seams, degree, 2, F(1, 64), 12)
        b = row['block_attempts']
        assert b >= F(seams, 8*degree)
        assert abs(row['phase_parameter']) <= 16*degree/F(64)
        coefficients = [complex(float(a), float(c)) for a, c in row['power_coefficients']]
        for probability in np.linspace(1-2*degree/seams, 1, 101):
            q = probability**b
            assert q >= .5-1e-14
            target = np.exp(1j*(2*seams/64)*np.log(probability))
            assert abs(np.polynomial.polynomial.polyval(q, coefficients)-target) < float(row['operator_error_upper'])+1e-11


def test_signed_time_conjugation_and_zero_time_are_exact():
    a = r.continuation(630, 6, 2, F(1, 64), 12)
    b = r.continuation(630, 6, 2, F(-1, 64), 12)
    assert b['power_coefficients'] == tuple((x, -y) for x, y in a['power_coefficients'])
    zero = r.continuation(630, 6, 2, F(), 12)
    assert zero['operator_error_upper'] == 0
    assert zero['coefficient_l1_upper'] == 1
    assert zero['power_coefficients'] == ((F(1), F()),)+((F(), F()),)*12


def test_positive_energy_matrix_kernel_and_reference_distinction(source):
    _, _, _, _, transfer, vals, vectors, probe = source
    hamiltonian = (vectors*(-60*np.log(vals)))@vectors.T
    omega = np.ones(66)/np.sqrt(66)
    assert np.linalg.norm(hamiltonian@omega) < 1e-12
    assert np.linalg.norm(np.outer(omega, omega)-np.eye(66)/66) > .9
    other = np.roll(probe, 1); columns = []
    for t, v in ((-.03, probe), (0, probe), (.02, other), (.05, probe)):
        columns.append(vectors@(np.exp(1j*60*t*np.log(vals))*(vectors.T@v/2)))
    gram = np.column_stack(columns).conj().T@np.column_stack(columns)
    assert np.linalg.eigvalsh(gram).min() > -1e-13
    assert np.linalg.norm(transfer@omega-omega) < 1e-13


def test_same_primitive_response_at_distinct_occupation_without_state_identification(source):
    family, _, _, h, _, _, _, _ = source
    for q in (2, 3, 4):
        lifted = r.exchange.sector.inclusion(h, 12, q, 2)
        norm = r.exchange.sector.dot(lifted, lifted)
        current = lifted
        for n in range(4):
            expected = r.exact_native_lags(family, h, 2, (n,))[0]
            assert F(r.exchange.sector.dot(lifted, current), norm*60**n) == expected
            current = r.exchange.sector.counts_step(family, current, q)


def test_exact_taylor_moments_and_finite_intervention_have_correct_sign(source):
    family, _, _, h, _, vals, vectors, probe = source
    t = F(1, 64); epsilon = F(1, 4096)
    witness = r.response_witness(family, h, time=t, kick=epsilon)
    assert witness['transfer_finite_kick_score_difference_lower'] == F(4551, 1476395008)
    g = 60*(1-vals); hm = -60*np.log(vals); weights = (vectors.T@probe)**2/4
    actual_g = -2*np.imag(weights@np.exp(-1j*float(t)*g))
    actual_h = -2*np.imag(weights@np.exp(-1j*float(t)*hm))
    assert actual_g >= float(witness['local_exchange_normalized_susceptibility_lower'])
    assert actual_h >= float(witness['transfer_normalized_susceptibility_lower'])
    omega = np.ones(66)/np.sqrt(66)
    kicked = np.exp(1j*float(epsilon)*probe)*omega
    evolved = vectors@(np.exp(-1j*float(t)*hm)*(vectors.T@kicked))
    difference = np.real(np.vdot(evolved, probe*evolved))
    assert difference >= float(witness['transfer_finite_kick_score_difference_lower'])
    wrong_kick = np.exp(-1j*float(epsilon)*probe)*omega
    wrong_evolved = vectors@(np.exp(-1j*float(t)*hm)*(vectors.T@wrong_kick))
    assert np.real(np.vdot(wrong_evolved, probe*wrong_evolved)) < 0
    moments = r.exchange.moments(family, h, 2, 8)
    a, b = r.unitary_taylor(moments, t)
    assert abs(complex(float(a), float(b))-weights@np.exp(-1j*float(t)*g)) <= float((20*t)**9/factorial(9))


def test_all_raw_records_and_reference_are_phase_blind_but_quantum_output_is_not():
    row = p.phase_blindness()
    assert row['complete_raw_slots'] == 721
    assert row['terminal_work_changed_slots'] > 0
    assert row['complete_record_reference_marginals_equal']
    plus = np.ones(2)/np.sqrt(2)
    phase = np.diag([1, 1j])
    assert abs(np.vdot(plus, phase@plus))**2 == pytest.approx(.5)
    assert abs(np.vdot(plus, plus))**2 == pytest.approx(1)


def test_mutated_native_partial_injection_is_rejected(monkeypatch):
    source = r.exchange.sector.source.interpreter
    original = source.local.raw_maps
    def bad(*args, **kwargs):
        for coin, label, pairs in original(*args, **kwargs):
            if len(pairs) > 1:
                pairs = list(pairs); pairs[1] = (pairs[1][0], pairs[0][1])
            yield coin, label, pairs
    source._local_mapping.cache_clear()
    try:
        with monkeypatch.context() as changes:
            changes.setattr(source.local, 'raw_maps', bad)
            with pytest.raises(ArithmeticError, match='partial injection|diagonal-phase sensitivity'):
                p.phase_blindness()
    finally:
        source._local_mapping.cache_clear()


def test_actual_native_program_budgets_keep_nonadmission_and_sampling_cost():
    rows = [r.source_budget(level, 80*4**level) for level in (0, 1, 2)]
    assert [x['all_setting_programs_admitted'] for x in rows] == [True, False, False]
    assert len(rows[0]['program_digests']) == 13
    assert rows[1]['program_digests'] == rows[2]['program_digests'] == ()
    assert rows[0]['maximum_attempts_per_experiment'] == 116424316
    for row in rows:
        assert row['sampling']['total_independent_repetitions'] > 10**17
        assert row['transfer_vacuum_systematic_error_upper'] < F(1, 4000)
        assert row['maximum_attempts_per_experiment'] > row['preparation']['warmup_attempts']


def test_sampling_union_bound_and_small_occupation_cost_are_explicit():
    budget = r.sampling_budget(F(10), F(1, 5), F(1, 100), F(1, 20), 4)
    n = budget['repetitions_per_setting']
    assert F(4*10**2, n)*25*10000 <= F(1, 20)
    high = r.sampling_budget(F(10), r.preparation.variance(240, 80), F(1, 100), F(1, 20), 4)
    low = r.sampling_budget(F(10), r.preparation.variance(240, 2), F(1, 100), F(1, 20), 4)
    assert low['total_independent_repetitions'] > high['total_independent_repetitions']*100000


@pytest.mark.parametrize('bad', (True, False, 1.0, '1'))
def test_typed_exact_domain_rejects_aliases(bad):
    with pytest.raises(ValueError):
        r.continuation(630, 6, 2, bad, 12)
    with pytest.raises(ValueError):
        r.continuation(630, 6, 2, F(1, 64), bad)
    with pytest.raises(ValueError):
        r.source_budget(bad, 80)


def test_resource_and_empty_normalization_domains_fail_closed(source):
    for args in ((10, 6, 2, F(1), 12), (630, 6, 2, F(10), 12), (630, 6, 2, F(0), 97)):
        with pytest.raises(ValueError): r.continuation(*args)
    with pytest.raises(ValueError): r.source_budget(3, 80)
    with pytest.raises(ValueError): r.source_budget(0, 1)
    with pytest.raises(ValueError): r.reconstruct(((F(1), F()),), ())
    with pytest.raises(ValueError): r.sampling_budget(F(1), F(0), F(1), F(1, 2), 1)
    family, _, _, h, *_ = source
    with pytest.raises(ArithmeticError): r.response_witness(family, h, kick=F(1))
    with pytest.raises(ValueError): r.exact_native_lags(family, {(0, 1): 1}, 2, (1,))


@pytest.fixture(scope='module')
def produced():
    return p.produce()


def test_retained_receipt_and_complete_import_source_bindings(produced):
    assert p.canonical(json.loads(p.OUTPUT.read_text())) == p.canonical(produced)
    assert all(name in produced['source_bindings'] for name in p.SOURCE_FILES)
    assert 'model/native_program/interpreter.py' in produced['source_bindings']
    assert 'model/native_correlation_preparation/PROOF.md' in produced['source_bindings']
    assert not produced['scope']['unitary_real_time_executed_by_native_repair']


@pytest.mark.parametrize('mutation', ('bool', 'witness', 'missing-source'))
def test_receipt_mutations_rejected_with_typed_canonical_comparison(produced, monkeypatch, tmp_path, mutation):
    altered = json.loads(p.canonical(produced))
    if mutation == 'bool': altered['local_source']['occupation'] = True
    elif mutation == 'witness': altered['local_source']['finite_intervention_witness']['native_bare_record_difference'] = '1'
    else: altered['source_bindings'].pop('model/native_program/interpreter.py')
    target = tmp_path/'bad.json'; target.write_text(p.canonical(altered))
    monkeypatch.setattr(p, 'produce', lambda: produced)
    with pytest.raises(ValueError): p.verify(target)


def test_source_race_guard_cannot_be_bypassed(monkeypatch):
    original = p.source_bindings; count = 0
    def changed():
        nonlocal count
        count += 1; value = original()
        if count > 1: value['model/native_program/interpreter.py'] = '0'*64
        return value
    monkeypatch.setattr(p, 'source_bindings', changed)
    with pytest.raises(ArithmeticError, match='source bytes changed'):
        p.produce()
