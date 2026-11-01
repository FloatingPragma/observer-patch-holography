from fractions import Fraction as F

import numpy as np
import pytest

from ..observer.source_scalar import compiler
from ..reference_recovery.recovery import execute_complete
from . import circuits as c, execute


@pytest.mark.parametrize('targets', [(4, 0), (1, 4), (3, 0, 4), (4, 2, 0)])
def test_local_block_contraction_matches_every_direct_source_prefix(targets):
    rng = np.random.default_rng(41)
    bits = len(targets)
    word = (compiler.cnot(0, 1)+[compiler.pulse('Y', .31, [0]), compiler.pulse('Z', -.27, [bits-1])])
    if bits == 3: word += compiler.cnot(2, 0)
    direct_word = [{**event, 'targets': [targets[j] for j in event['targets']]} for event in word]
    initial = rng.normal(size=(32, 3))+1j*rng.normal(size=(32, 3)); initial /= np.linalg.norm(initial)
    direct = execute_complete(direct_word, 5, initial)
    block = c.block(word, bits, targets)
    contracted = execute.complete([block], 5, initial)
    np.testing.assert_allclose(contracted.success, direct.success, atol=4e-16, rtol=0)
    np.testing.assert_allclose(contracted.failure_gram, direct.failure_gram, atol=1e-45, rtol=2e-12)
    np.testing.assert_allclose(contracted.failure_weights, direct.failure_weights, atol=1e-45, rtol=2e-12)
    for event, (_, _, local_targets) in enumerate(block[0].tape):
        for slot in range(6-(1 << len(local_targets))):
            reconstructed = contracted.failure(0, event, slot)
            expected = direct.failure(event)[slot]
            np.testing.assert_allclose(reconstructed.conj().T@reconstructed, expected.conj().T@expected,
                                       atol=1e-45, rtol=2e-12)


def test_generic_order_one_export_failure_remains_complete(monkeypatch):
    rng = np.random.default_rng(73)
    a = rng.normal(size=(6, 6))+1j*rng.normal(size=(6, 6))
    u, _ = np.linalg.qr(a)
    execute.compile_block.cache_clear()
    monkeypatch.setattr(compiler, 'native_pulse', lambda *_: u)
    try:
        word = [compiler.pulse('Y', .731, [0]), compiler.pulse('J', -.893, [2, 1])]
        initial = rng.normal(size=(32, 2))+1j*rng.normal(size=(32, 2)); initial /= np.linalg.norm(initial)
        run = execute.complete([c.block(word, 3, (3, 0, 4))], 5, initial)
        routed = [{**event, 'targets': [(3, 0, 4)[j] for j in event['targets']]} for event in word]
        direct = execute_complete(routed, 5, initial)
        assert run.failure_probability > .1
        np.testing.assert_allclose(run.success, direct.success, atol=3e-16, rtol=0)
        for event, (_, _, targets) in enumerate(run.blocks[0][0].tape):
            for slot in range(6-(1 << len(targets))):
                x, y = run.failure(0, event, slot), direct.failure(event)[slot]
                np.testing.assert_allclose(x.conj().T@x, y.conj().T@y, atol=3e-16, rtol=0)
        np.testing.assert_allclose(run.success.conj().T@run.success+run.failure_gram,
                                   initial.conj().T@initial, atol=1e-14, rtol=0)
        assert run.failure_sum_error <= run.failure_sum_admission
    finally: execute.compile_block.cache_clear()


def test_later_block_failure_reconstructs_all_earlier_successes():
    word = [compiler.pulse('Y', .39, [0]), compiler.pulse('Z', -.82, [1])]
    first = c.block(word, 2, (2, 0)); second = c.block(word, 2, (1, 2))
    state = np.eye(8, dtype=complex)[:, :3]/np.sqrt(3)
    run = execute.complete([first, second], 3, state)
    native = [dict(kind=k, angle=a, targets=[targets[j] for j in local])
              for block, targets in (first, second) for k, a, local in block.tape]
    direct = execute_complete(native, 3, state)
    for slot in range(4):
        x, y = run.failure(1, 1, slot), direct.failure(3)[slot]
        np.testing.assert_allclose(x.conj().T@x, y.conj().T@y, atol=1e-45, rtol=2e-12)


def test_pointer_control_is_exact_local_number_phase_and_analyzer_reverses_it():
    control = c.block(c.control_word(), 3, (0, 1, 2))[0]
    expected = np.diag([np.exp(1j*np.pi*(mask & 3).bit_count()/2) if mask & 4 else 1 for mask in range(8)])
    np.testing.assert_allclose(control.success, expected, atol=2e-14, rtol=0)
    # Both pointer phases and all four local occupation inputs are retained.
    blank_pointer = np.eye(8, dtype=complex)[:, :4]
    prepared = execute.complete([c.block(c.pointer_prepare(), 3, (0, 1, 2))], 3, blank_pointer)
    reversed_run = execute.complete([c.block(c.pointer_analyzer(F()), 3, (0, 1, 2))], 3, prepared.success)
    np.testing.assert_allclose(reversed_run.success, blank_pointer, atol=3e-14, rtol=0)


def test_field_blocks_are_the_identical_original_source_word():
    blocks = c.field_binding(F(1, 2), 1)
    assert len(blocks) == 416
    assert len(c.flatten(blocks)) == 10560


def test_padding_is_not_an_extra_failure_record():
    initial = np.zeros(8); initial[0] = 1
    run = execute.complete([c.block([compiler.pulse('J', .1, [0, 1])], 2, (2, 0))], 3, initial)
    assert np.array_equal(run.failure_weights[:, 2:], np.zeros((1, 2)))
    with pytest.raises(ValueError): run.failure(0, 0, 2)
