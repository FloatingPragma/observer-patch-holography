"""Sparse target, complete native source, relative phases and negative controls."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy import sparse

from ..common_spin_protection import spin
from ..observer.source_scalar import compiler
from . import circuits as c, exact, execute, experiment as e


def blank(bits=13):
    out = np.zeros((1 << bits, 1), complex); out[0, 0] = 1
    return out


def test_sparse_exact_preparation_has_all_four_signed_amplitudes():
    value = exact.preparation_certificate()
    assert value['Hubbard_occupation_numerators'] == [[53, 1], [197, -1], [773, 1], [3077, -1]]
    assert value['generic_amplitude_tree_used'] is False
    assert value['charge_postselection_used'] is False
    state = exact.target()
    assert np.linalg.norm(state) == 1
    for law in ('C', 'E'):
        for seam in ((0, 1), (0, 2), (1, 2)):
            for label, op in spin.repair_outcomes('hole', law, seam):
                expected = state if label == ('outside',) else np.zeros_like(state)
                assert np.linalg.norm(op@state-expected) < 1e-15


def test_literal_twelve_bit_preparation_and_pointer_reference_columns():
    word = c.preparation_word(pointer=False)
    run = execute.complete(word, 12, blank(12))
    assert len(word) == 292
    assert np.linalg.norm(run.success[:, 0]-exact.target()) < 2e-14
    # Two independent pointer input columns verify the coherent controlled
    # isometry, not only one selected state or one output probability.
    initial = np.column_stack((blank()[:, 0], np.roll(blank()[:, 0], 4096)))
    controlled = execute.complete(c.preparation_word(), 13, initial)
    target = np.column_stack((e.ideal_wait(F())[:, 0],
                             np.r_[-blank(12)[:, 0], exact.target()]/np.sqrt(2)))
    assert np.linalg.norm(controlled.success-target) < 8e-14
    np.testing.assert_allclose(controlled.success.conj().T@controlled.success+controlled.failure_gram,
                               np.eye(2), atol=1e-12, rtol=0)
    event = int(np.argmax(controlled.failure_weights.sum(axis=1)))
    branches = controlled.failure(event)
    expected = [np.vdot(branch, branch).real for branch in branches]
    np.testing.assert_allclose(expected, controlled.failure_weights[event, :len(branches)], atol=1e-42, rtol=1e-12)


def test_complete_executor_keeps_a_genuinely_nonzero_source_abort(monkeypatch):
    matrix = np.eye(6, dtype=complex)
    angle = .2
    matrix[np.ix_([0, 2], [0, 2])] = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    monkeypatch.setattr(compiler, 'native_pulse', lambda kind, angle: matrix.copy())
    run = execute.complete([compiler.pulse('I', 0., [0])], 13, blank())
    assert run.failure_probability == pytest.approx(np.sin(.2)**2)
    assert np.linalg.norm(run.success)**2+run.failure_probability == pytest.approx(1)
    assert len(run.failure(0)) == 4
    assert np.vdot(run.failure(0)[0], run.failure(0)[0]).real == pytest.approx(run.failure_probability)


def test_scalar_phase_omission_and_pointer_dephasing_destroy_the_ramsey_claim():
    plan = list(c.preparation_plan())
    pauli, _ = plan[0]; plan[0] = pauli, F()  # retain pulse but omit the required scalar phase
    wrong_word = [compiler.pulse('Y', -np.pi/4, [12])]+c.emit(c.controlled(plan))+c.analyzer_word(F())
    wrong = execute.complete(wrong_word, 13, blank())
    p = [np.vdot(x, x).real for x in execute.pointer_outputs(wrong.success)]
    assert p[0] == pytest.approx(.5, abs=1e-12)
    expected = e.ideal_wait(F())
    dephased = np.zeros((8192, 2), complex)
    dephased[:4096, 0] = expected[:4096, 0]
    dephased[4096:, 1] = expected[4096:, 0]
    run = execute.complete(c.analyzer_word(F()), 13, dephased)
    assert np.vdot(execute.pointer_outputs(run.success)[0], execute.pointer_outputs(run.success)[0]).real == pytest.approx(.5, abs=1e-12)


def test_exact_literal_source_certificate_rejects_changed_control_or_target():
    word = c.preparation_word()
    value = c.certify_word(word, c.preparation_template())
    assert float(value['complete_source_operator_error']) < 2e-13
    word[4] = {**word[4], 'targets': [11]}
    with pytest.raises(ArithmeticError):
        c.certify_word(word, c.preparation_template())
    word = c.preparation_word()
    word[0] = {**word[0], 'angle': word[0]['angle']+.1}
    assert float(c.certify_word(word, c.preparation_template())['complete_source_operator_error']) > .09


@pytest.mark.parametrize('bits,initial', [(11, np.zeros(2048)), (14, np.zeros(16384)),
                                        (13, np.zeros(4096)), (13, np.zeros((8192, 5))),
                                        (13, np.full(8192, np.nan))])
def test_executor_invalid_complete_domains_reject(bits, initial):
    with pytest.raises(ValueError):
        execute.complete([], bits, initial)
