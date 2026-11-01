"""Full-input source, complete-output and phase-sensitive mutation controls."""
from dataclasses import replace
from fractions import Fraction
from types import MappingProxyType

import numpy as np
import pytest

from . import many_body_service as s
from ..defect_field_bridge import full_domain


@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('character', [1, -1])
def test_every_source_column_matches_actual_native_seam_on_all_occupations(law, character):
    bits = s._domain(law, character)
    target = s.target_columns(law, character)
    for mask in range(16):
        native = sum(((mask >> i) & 1) << p for i, p in enumerate(s.NATIVE_PORTS))
        for coin in (0, 1):
            raw, moved, sign, _ = full_domain.branch_action(
                native, law, coin, sites=12, flavors=1, seam=s.NATIVE_SEAM)
            packed = sum(((moved >> p) & 1) << i for i, p in enumerate(s.NATIVE_PORTS))
            assert target[packed+(coin << 4)+(raw << 5), mask] == (
                (sign if character == -1 else 1)/np.sqrt(2))
        assert np.count_nonzero(target[:, mask]) == 2
    word = s.repair_word(law, character)
    initial = np.eye(1 << bits, dtype=complex)[:, :16]
    run = s.recovery.execute_complete(word, bits, initial)
    assert s.admit(run, target, word=word, bits=bits, initial=initial) < float(s.SOURCE_CEILING)
    assert run.failure_weights.shape == (len(run.word), 4)
    # Actual complement outputs remain available, including zero outcomes.
    for event in (0, len(run.word)//2, len(run.word)-1):
        failed = run.failure(event)
        assert len(failed) == (2 if run.word[event][0] == 'J' else 4)
        assert sum(np.vdot(x, x).real for x in failed) == pytest.approx(
            sum(run.failure_weights[event]), abs=1e-28)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_fixed_number_observer_distinguishes_lifts_and_which_path_erases_readout(law):
    clean = [s.experiment(law, c) for c in (1, -1)]
    marked = [s.experiment(law, c, which_path=True) for c in (1, -1)]
    assert clean[0]['ideal_nonzero_terminal_probabilities'] != clean[1]['ideal_nonzero_terminal_probabilities']
    assert marked[0]['ideal_nonzero_terminal_probabilities'] == marked[1]['ideal_nonzero_terminal_probabilities']
    for result in clean+marked:
        assert sum(map(Fraction, result['ideal_nonzero_terminal_probabilities'].values())) == 1
        assert result['total_input_occupation'] == 2
        assert result['resources']['proper_code_service_events'] == 2*result['resources']['native_pulses']
        assert result['resources']['first_failure_Kraus_slots'] > 0


def test_wrong_but_complete_positive_lift_does_not_certify_negative_lift():
    word, initial = s.repair_word('C', 1), np.eye(128, dtype=complex)[:, :16]
    run = s.recovery.execute_complete(word, 7, initial)
    assert np.linalg.norm(run.success.conj().T@run.success-np.eye(16)) < 1e-10
    with pytest.raises(ArithmeticError, match='admission ceiling'):
        s.admit(run, s.target_columns('C', -1), word=word, bits=7, initial=initial)


@pytest.mark.parametrize('angle, accepted', [(1e-12, True), (1e-4, False)])
def test_complete_native_export_leakage_is_not_replaced_by_exact_zero(monkeypatch, angle, accepted):
    original = s.compiler.native_pulse
    def perturbed(kind, theta):
        value = original(kind, theta)
        if kind == 'Y' and theta == -np.pi/4:
            rotation = np.eye(6, dtype=complex)
            rotation[0, 0] = rotation[5, 5] = np.cos(angle)
            rotation[0, 5], rotation[5, 0] = -np.sin(angle), np.sin(angle)
            return rotation@value
        return value
    monkeypatch.setattr(s.compiler, 'native_pulse', perturbed)
    word, initial = s.repair_word('C', 1), np.eye(128, dtype=complex)[:, :16]
    run = s.recovery.execute_complete(word, 7, initial)
    assert run.failure_probability > 0
    if accepted:
        s.admit(run, s.target_columns('C', 1), word=word, bits=7, initial=initial)
    else:
        with pytest.raises(ArithmeticError, match='admission ceiling'):
            s.admit(run, s.target_columns('C', 1), word=word, bits=7, initial=initial)


@pytest.mark.parametrize('kind', ['integer', 'subclass', 'nonfinite', 'indefinite_gram'])
def test_malformed_outputs_cannot_hide_wrong_source_or_failure_mass(kind):
    target = s.target_columns('C', 1)
    word, initial = s.repair_word('C', 1), np.eye(128, dtype=complex)[:, :16]
    run = s.recovery.execute_complete(word, 7, initial)
    actual = target.copy()
    gram = np.zeros((16, 16), complex)
    weights = np.zeros((len(word), 4), float)
    if kind == 'integer':
        actual = np.zeros_like(target, dtype=np.int64)
        actual[0, 0] = 2**32
    elif kind == 'subclass':
        class Liar(np.ndarray):
            def __sub__(self, other):
                return np.zeros(self.shape)
        actual = actual.view(Liar)
        actual[0, 0] = 1
    elif kind == 'nonfinite':
        weights[0, 0] = np.nan
    else:
        gram[0, 0], gram[1, 1] = 1e-10, -1e-10
    with pytest.raises((ValueError, ArithmeticError)):
        s.admit(replace(run, success=actual, failure_gram=gram, failure_weights=weights),
                target, word=word, bits=7, initial=initial)


def test_public_experiment_rejects_dropped_failure_archive_event(monkeypatch):
    original = s.recovery.execute_complete
    def truncated(*args, **kwargs):
        run = original(*args, **kwargs)
        return replace(run, word=run.word[:-1])
    monkeypatch.setattr(s.recovery, 'execute_complete', truncated)
    with pytest.raises(ValueError, match='exact requested source word'):
        s.experiment('C', -1)


def test_indefinite_failure_gram_is_rejected_even_when_norm_equals_positive_trace():
    word, initial = s.repair_word('C', -1), np.eye(128, dtype=complex)[:, :16]
    run = s.recovery.execute_complete(word, 7, initial)
    gram = np.zeros((16, 16), complex)
    gram[0, 0] = gram[1, 1] = 2e-20
    gram[2, 2] = -1e-20
    weights = np.zeros((len(word), 4), float)
    weights[0, 0] = 3e-20
    forged = replace(run, failure_gram=gram, failure_weights=weights)
    with pytest.raises(ArithmeticError, match='literal source prefixes'):
        s.admit(forged, s.target_columns('C', -1), word=word, bits=7, initial=initial)


@pytest.mark.parametrize('mutation', (
    'word', 'bits', 'initial', 'weight_row', 'pair_padding', 'primitive', 'missing_primitive',
))
def test_complete_archive_is_bound_to_requested_execution(mutation):
    word, initial = s.repair_word('E', -1), np.eye(256, dtype=complex)[:, :16]
    run = s.recovery.execute_complete(word, 8, initial)
    if mutation == 'word':
        forged = replace(run, word=tuple(reversed(run.word)))
    elif mutation == 'bits':
        forged = replace(run, bits=True)
    elif mutation == 'initial':
        other = run.initial.copy(); other[0, 0] *= -1
        forged = replace(run, initial=other)
    elif mutation == 'weight_row':
        forged = replace(run, failure_weights=run.failure_weights[:-1])
    elif mutation == 'pair_padding':
        weights = run.failure_weights.copy()
        event = next(i for i, (kind, _, _) in enumerate(run.word) if kind == 'J')
        weights[event, 3] = 1e-30
        forged = replace(run, failure_weights=weights)
    else:
        primitives = dict(run.primitives)
        key = next(iter(primitives))
        if mutation == 'primitive':
            primitives[key] = -primitives[key]
        else:
            del primitives[key]
        forged = replace(run, primitives=MappingProxyType(primitives))
    with pytest.raises((ValueError, ArithmeticError)):
        s.admit(forged, s.target_columns('E', -1), word=word, bits=8, initial=initial)


def _literal_step(state, full, targets, bits):
    """Independent bit-address routing of every complement amplitude."""
    size = 1 << len(targets)
    others = [p for p in range(bits) if p not in targets]
    success = np.zeros_like(state)
    failures = [np.zeros((1 << len(others), state.shape[1]), complex) for _ in range(6-size)]
    for address, vector in enumerate(state):
        local = sum(((address >> p) & 1) << (len(targets)-1-i) for i, p in enumerate(targets))
        remainder = sum(((address >> p) & 1) << i for i, p in enumerate(others))
        blank = address & ~sum(1 << p for p in targets)
        for output in range(size):
            moved = blank | sum(((output >> (len(targets)-1-i)) & 1) << p
                                for i, p in enumerate(targets))
            success[moved] += full[output, local]*vector
        for slot, failure in enumerate(failures):
            failure[remainder] += full[size+slot, local]*vector
    return success, failures


def test_nonzero_failure_vectors_match_independent_address_prefixes_with_reference(monkeypatch):
    original = s.compiler.native_pulse
    def leaking(kind, angle):
        full = original(kind, angle)
        rotation = np.eye(6, dtype=complex)
        theta = 1e-10
        rotation[0, 0] = rotation[5, 5] = np.cos(theta)
        rotation[0, 5], rotation[5, 0] = -np.sin(theta), np.sin(theta)
        return rotation@full
    monkeypatch.setattr(s.compiler, 'native_pulse', leaking)
    word = s.repair_word('C', -1)
    rng = np.random.default_rng(139)
    initial = rng.normal(size=(128, 3))+1j*rng.normal(size=(128, 3))
    initial /= np.linalg.norm(initial)
    run = s.recovery.execute_complete(word, 7, initial)
    checked = {0, len(word)-1,
               next(i for i, (kind, _, targets) in enumerate(run.word)
                    if kind == 'J' and abs(targets[0]-targets[1]) > 1)}
    state = initial.copy()
    for event, (kind, angle, targets) in enumerate(run.word):
        state, failures = _literal_step(state, run.primitives[kind, angle], targets, 7)
        if event in checked:
            retained = run.failure(event)
            assert len(retained) == len(failures)
            for actual, expected in zip(retained, failures, strict=True):
                assert np.linalg.norm(actual-expected) <= (
                    64*np.finfo(float).eps*(event+1)*np.linalg.norm(expected))
    assert run.failure_probability > 0


@pytest.mark.parametrize('name,value', (
    ('NATIVE_PORTS', (0, 6, 5, 1)), ('NATIVE_PORTS', (0, True, 5, 6)),
    ('NATIVE_SEAM', (0, 1)), ('SEAM', (0, 3)),
))
def test_public_observer_rejects_wrong_native_embedding(monkeypatch, name, value):
    monkeypatch.setattr(s, name, value)
    with pytest.raises(ValueError, match='active-port chart'):
        s.experiment('C', -1)


def test_source_certificate_consumes_native_embedding_not_only_local_target(monkeypatch):
    original = s.many_body.column
    def wrong_sign(mask, law, seam, coin, character=-1):
        raw, output, sign, current = original(mask, law, seam, coin, character)
        return raw, output, -sign if mask == 1 and current else sign, current
    monkeypatch.setattr(s.many_body, 'column', wrong_sign)
    with pytest.raises(ArithmeticError, match='embedded native instrument'):
        s.certificate()


@pytest.mark.parametrize('law, character', [('X', 1), ('C', True), ('E', 0), ('C', -1.)])
def test_source_domain_is_explicit(law, character):
    with pytest.raises(ValueError):
        s.repair_word(law, character)
