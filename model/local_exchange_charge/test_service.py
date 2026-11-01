"""Independent complete-column, preparation and native-history checks."""
from itertools import product
from dataclasses import replace
import warnings

import numpy as np
import pytest

from . import service as s


def record(mask):
    return tuple((mask >> j) & 1 for j in range(4))


@pytest.mark.parametrize('law', ('C', 'E'))
def test_independent_target_retains_all_inputs_records_and_collisions(law):
    target = s.independent_target(law)
    assert target.shape == (512, 16)
    assert np.allclose(target.conj().T@target, np.eye(16), atol=1e-15)
    for mask, coin in product(range(16), range(2)):
        bits = [(mask >> j) & 1 for j in range(4)]
        counts = (sum(bits[:2]), sum(bits[2:]))
        found = np.flatnonzero(target[:, mask])
        destination = next(int(x) for x in found if (x >> 4) & 1 == coin)
        raw, output = destination >> 5, destination & 15
        if max(counts) == 2:
            assert raw == 8 and output == mask
        elif sum(counts) == 1:
            flavor = next(j % 2 for j, value in enumerate(bits) if value)
            endpoint = 0 if coin else 1
            assert output == 1 << (2*endpoint+flavor)
            current = counts[0]-int(endpoint == 0)
            assert raw == {0: 0, 1: 1, -1: 2}[current]+(4 if law == 'E' else 0)
        else:
            assert output == mask and raw == 0
    # A full site SWAP would incorrectly exchange these occupied endpoints.
    for mask in (s.AB, s.BA, 5, 10):
        assert target[mask, mask] == 1/np.sqrt(2)


@pytest.mark.parametrize('law', ('C', 'E'))
def test_literal_complete_source_and_all_arbitrary_reference_columns(law):
    run = s.repair_run(law)
    target = s.independent_target(law)
    assert np.linalg.norm(run.success-target)+np.sqrt(run.failure_probability) < s.ADMISSION
    assert len(s.repair_kraus(law)) == 32
    complete = sum(k.conj().T@k for k in s.repair_kraus(law).values())+run.failure_gram
    assert np.linalg.norm(complete-np.eye(16)) < s.ADMISSION
    rng = np.random.default_rng(715)
    psi = rng.normal(size=(16, 3))+1j*rng.normal(size=(16, 3))
    psi /= np.linalg.norm(psi)
    after = sum((k@psi).conj().T@(k@psi) for k in s.repair_kraus(law).values())
    after += psi.conj().T@run.failure_gram@psi
    assert np.linalg.norm(after-psi.conj().T@psi) < s.ADMISSION
    # Actual prefix complement amplitudes are addressable, including zeros.
    for event in (0, len(run.word)//2, len(run.word)-1):
        failures = run.failure(event)
        assert len(failures) in (2, 4)
        assert all(branch.shape[-1] == 16 for branch in failures)


def test_preparation_is_nonunital_local_complete_channel_not_unitary_conjugation():
    for reference in (False, True):
        effects = np.zeros((16, 16), complex)
        for raw, selector in product(s.READOUT_RECORDS, range(2)):
            k = s.preparation_kraus(raw, selector, reference=reference)
            effects += k.conj().T@k
        effects += s.preparation_failure_effect(reference=reference)
        assert np.linalg.norm(effects-np.eye(16)) < s.ADMISSION
    prepared = s.prepared_density()
    expected = s.local_reference().astype(complex)
    expected[s.AB, s.BA] = expected[s.BA, s.AB] = 1/264
    assert np.linalg.norm(prepared-expected) < s.ADMISSION
    assert np.linalg.norm(s.prepared_density(reference=True)-s.local_reference()) < s.ADMISSION
    assert abs(np.trace(prepared).real-1) < s.ADMISSION
    assert np.linalg.norm(prepared-s.local_reference()) > 1/300
    assert not np.any(s.preparation_kraus((2, 0, 0, 0), 0))
    # Preparation selector records really differ; only their traced outputs
    # belong to the subsequent native-tape equality statement.
    assert s._selector_weights(False) != s._selector_weights(True)


def test_pair_charge_analyzer_has_three_complete_outcomes_and_keeps_noncode_slots():
    q = np.zeros((16, 16)); q[s.AB, s.BA] = q[s.BA, s.AB] = 1
    effects = np.zeros((16, 16), complex)
    measured = np.zeros_like(effects)
    for raw in s.READOUT_RECORDS:
        k = s.readout_kraus(raw)
        effects += k.conj().T@k
        measured += s.readout_value(raw)*(k.conj().T@k)
    assert np.linalg.norm(effects+s.analyzer_run().failure_gram-np.eye(16)) < s.ADMISSION
    assert np.linalg.norm(measured-q) < s.ADMISSION
    assert s.measurement_certificate()['complete_raw_slots'] == 1296
    assert np.count_nonzero(s.readout_kraus((5, 0, 1, 1))) == 0
    assert np.linalg.norm(s.analyzer_run().success-s.ideal_basis().T) < s.ADMISSION


def test_complete_preparation_native_C_E_charge_readout_without_postselection():
    defect, reference = s.end_to_end(), s.end_to_end(reference=True)
    assert len(defect['rows']) == 32*32
    assert abs(defect['complete_mass']-1) < 4*s.ADMISSION
    assert np.max(abs(defect['charge_outcomes_plus_minus_zero']-[1/88, 1/264, 65/66])) < 4*s.ADMISSION
    assert np.max(abs(reference['charge_outcomes_plus_minus_zero']-[1/132, 1/132, 65/66])) < 4*s.ADMISSION
    # Subsequent native records are label-blind after the prep archive is traced.
    for changed, unchanged in zip(defect['rows'], reference['rows'], strict=True):
        assert changed['C'] == unchanged['C'] and changed['E'] == unchanged['E']
        assert abs(sum(changed['all16_code_probabilities'])-sum(unchanged['all16_code_probabilities'])) < 4*s.ADMISSION
    psi = np.zeros((16, 2), complex); psi[s.AB] = [1, 1j]
    output = s.history_amplitude(record(s.AB), 0, (0, 0), (1, 0), record(s.BA), psi)
    expected = np.zeros_like(psi); expected[s.BA] = np.sqrt(.75)*psi[s.AB]/2
    assert np.linalg.norm(output-expected) < 4*s.ADMISSION


@pytest.mark.parametrize('law', ('C', 'E'))
def test_actual_three_port_motion_retains_ambient_rows_and_all_records(law):
    changed = s.moving_history(law)
    reference = s.moving_history(law, reference=True)
    assert changed['prepared_density'].shape == (64, 64)
    assert len(changed['rows']) == 32
    assert np.max(abs(changed['charge_outcomes_plus_minus_zero']-[5/528, 1/176, 65/66])) < s.HISTORY_ADMISSION
    assert abs(changed['complete_mass']-1) < s.HISTORY_ADMISSION
    assert abs(changed['charge_outcomes_plus_minus_zero'][0]-changed['charge_outcomes_plus_minus_zero'][1]-1/264) < s.HISTORY_ADMISSION
    for row, before in zip(changed['rows'], reference['rows'], strict=True):
        assert row['native_record'] == before['native_record']
        assert abs(sum(row['all16_code_probabilities'])-sum(before['all16_code_probabilities'])) < s.HISTORY_ADMISSION
    # A source particle at port0 moves past the untouched port1 into port5.
    # Tensor-local label transport must not add a spectator-dependent CAR sign.
    moved_raw = 1 if law == 'C' else 5
    local = s.repair_kraus(law)[0, moved_raw]
    lifted = s.on_three_ports(local, (0, 2))
    input_ab = (1 << 0)+(1 << 3)
    output_ab = (1 << 4)+(1 << 3)
    assert abs(lifted[output_ab, input_ab]-1/np.sqrt(2)) < s.ADMISSION
    complete = sum(s.on_three_ports(k, (0, 2)).conj().T@s.on_three_ports(k, (0, 2))
                   for k in s.repair_kraus(law).values())
    complete += s.on_three_ports(s.repair_run(law).failure_gram, (0, 2))
    assert np.linalg.norm(complete-np.eye(64)) < 4*s.ADMISSION
    psi = np.zeros((64, 3), complex); psi[s.AB] = [1, 1j, -1]
    output = s.moving_amplitude(law, record(s.AB), 0, (0, moved_raw), record(s.BA), psi)
    expected = np.zeros_like(psi); expected[24] = np.sqrt(.75/2)*psi[s.AB]
    assert np.linalg.norm(output-expected) < s.HISTORY_ADMISSION


def test_terminal_demolition_can_destroy_global_charge_on_a_split_pair():
    state = np.zeros(64); state[s.AB] = state[s.BA] = 1/np.sqrt(2)
    density = np.outer(state, state)
    swap = np.zeros((64, 64))
    for mask in range(64):
        target = sum((((mask >> (2*j)) & 1) << (2*j+1))+
                     (((mask >> (2*j+1)) & 1) << (2*j)) for j in range(3))
        swap[target, mask] = 1
    observed = np.zeros((64, 64), complex)
    for mask in range(16):
        k = s.on_three_ports(s.readout_kraus(record(mask)), (1, 2))
        observed += k@density@k.conj().T
    assert abs(np.trace(swap@density)-1) < 1e-14
    assert abs(np.trace(swap@observed)) < s.ADMISSION


def test_public_run_snapshots_cannot_poison_cached_array_metadata():
    run = s.preparation_run(s.AB, 1)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        run.success.shape = (16,)
        run.initial.shape = (16,)
    assert s.preparation_run(s.AB, 1).success.shape == (16, 1)
    assert s.preparation_run(s.AB, 1).initial.shape == (16, 1)
    for wrong in (True, 1.0, np.int64(1)):
        with pytest.raises(ValueError):
            s.preparation_run(s.AB, wrong)
    with pytest.raises(ValueError):
        s.preparation_run(s.AB, 1).success.setflags(write=True)


def test_physical_native_readout_includes_every_noncode_zero_slot():
    assert s.measurement_certificate(5)['complete_raw_slots'] == 7776
    assert np.array_equal(s.repair_record_kraus('C', (0, 1, 0, 0, 0)), s.repair_kraus('C')[0, 1])
    assert not np.any(s.repair_record_kraus('C', (4, 1, 0, 0, 0)))


@pytest.mark.parametrize('law,raw', [('C', 2), ('E', 6)])
def test_moving_poststate_preserves_unread_unknown_label_and_reference(law, raw):
    actual = s.moving_amplitude(law, record(0), 0, (1, raw), record(0), np.eye(64))
    target = np.zeros((64, 64), complex)
    target[1, 16] = target[2, 32] = np.sqrt(.75/2)
    assert np.linalg.norm(actual-target) < s.HISTORY_ADMISSION
    entangled = np.zeros((64, 2), complex)
    entangled[16, 0] = 1/np.sqrt(2); entangled[32, 1] = 1j/np.sqrt(2)
    assert np.linalg.norm(actual@entangled-target@entangled) < s.HISTORY_ADMISSION
    wrong = actual.copy(); wrong[2] *= -1
    # A label-Z on the unread port changes the actual output despite matching
    # every computational output probability of this input.
    assert np.allclose(abs(wrong@entangled)**2, abs(actual@entangled)**2)
    assert np.linalg.norm(wrong@entangled-target@entangled) > .5


def test_nonfinite_complete_source_output_cannot_pass_a_ceiling(monkeypatch):
    good = s.analyzer_run()
    s._analyzer_run.cache_clear()
    bad = replace(good, success=np.full((16, 16), np.nan, complex))
    with monkeypatch.context() as patch:
        patch.setattr(s, 'execute_complete', lambda *args: bad)
        with pytest.raises(ArithmeticError):
            s.analyzer_run()
    s._analyzer_run.cache_clear()


@pytest.mark.parametrize('deleted', [8, 15])
def test_dropping_collision_or_padded_zero_native_slots_is_rejected(monkeypatch, deleted):
    original = s.repair_kraus
    monkeypatch.setattr(s, 'repair_kraus', lambda law: {key: value for key, value in original(law).items() if key[1] != deleted})
    with pytest.raises(ArithmeticError):
        s.end_to_end()


def test_probability_invisible_unread_label_phase_is_rejected(monkeypatch):
    original = s.on_three_ports
    def wrong(matrix, sites):
        value = original(matrix, sites)
        if sites == (1, 2):
            value[[mask for mask in range(64) if mask & 3 == 2]] *= -1
        return value
    monkeypatch.setattr(s, 'on_three_ports', wrong)
    with pytest.raises(ArithmeticError):
        s.moving_history()


@pytest.mark.parametrize('call', [lambda: s.independent_target('bad'),
                                  lambda: s.preparation_kraus((0, 0, 0, 0), True),
                                  lambda: s.preparation_run(0, 0),
                                  lambda: s.readout_kraus((0, 0, 0, 6)),
                                  lambda: s.prepared_density(reference=1),
                                  lambda: s.repair_run(np.str_('C')),
                                  lambda: s.history_amplitude(record(s.AB), 0, (False, 0), (0, 0), record(s.BA), np.ones((16, 1))),
                                  lambda: s.moving_amplitude('C', record(s.AB), 0, (0, 1.), record(s.BA), np.ones((64, 1))),
                                  lambda: s.on_three_ports(np.ma.array(np.eye(16)), (0, 1)),
                                  lambda: s.history_amplitude(record(s.AB), 0, (0, 0), (0, 0), record(s.BA), np.ma.array(np.ones((16, 1))))])
def test_malformed_local_service_inputs_fail_closed(call):
    with pytest.raises(ValueError):
        call()
