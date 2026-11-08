from dataclasses import replace
from fractions import Fraction as F
import pytest
from . import chsh


def test_complete_unconditioned_score_includes_abort_sector():
    result = chsh.run()
    assert F(result['complete_CHSH_score']) > F(12, 5)
    assert result['conditional_acceptance_normalization_used'] is False
    assert result['abort_output_pair'] == [1, 1]
    assert len(result['abort_herald_records']) == 23
    assert all(len(setting['accepted_herald_histories']) == 2 for setting in result['settings'])
    assert all(len(history['complete_native_readout_law']) == 625
               for setting in result['settings'] for history in setting['accepted_herald_histories'])
    # Removing the scored abort sector gives the familiar sub-threshold
    # weighted success score; the complete convention is materially used.
    assert F(result['complete_CHSH_score']) - 2 * F(result['abort_mass']) < 2


def test_dephasing_the_actual_native_herald_fails_the_quantum_score(monkeypatch):
    actual = chsh.native.execute
    calls = 0
    def dephased(program, density, **kwargs):
        nonlocal calls
        history = actual(program, density, **kwargs)
        calls += 1
        if calls == 1:
            # Keep every native record and its probability exactly. Remove
            # only source coherences in the herald output. A classical
            # replacement must not pass this quantum correlation test.
            history = replace(history, branches=tuple(replace(row, density={
                key: value for key, value in row.density.items() if key[0] == key[2]
            }) for row in history.branches))
        return history
    monkeypatch.setattr(chsh.native, 'execute', dephased)
    with pytest.raises(ArithmeticError, match='stored complete CHSH score exceeds'):
        chsh.run()


def test_discarding_zero_native_records_fails_the_complete_protocol(monkeypatch):
    actual = chsh.native.execute
    def incomplete(program, density, **kwargs):
        history = actual(program, density, **kwargs)
        return replace(history, branches=tuple(row for row in history.branches if row.density))
    monkeypatch.setattr(chsh.native, 'execute', incomplete)
    with pytest.raises(ArithmeticError, match='retain all 25 raw slots'):
        chsh.run()
