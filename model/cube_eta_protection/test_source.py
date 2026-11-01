"""Native fixed-number target, retained aborts and independent source error."""
from fractions import Fraction as F

import numpy as np
import pytest

from ..observer.source_scalar import compiler
from . import source


def test_compact_native_tree_copy_and_local_phase_with_all_failures():
    value = source.report()
    # Compare the diagnostic with the independently proved complete-source
    # ceiling, rather than a tighter platform-specific observed residual.
    proved = F(source.certificate()['complete_source_state_error_bound'])
    assert proved < F(848, 10**15)
    assert value['native_preparation_error'] < float(proved)
    assert value['complete_preparation_and_phase_error'] < float(proved)
    assert [r['pulses'] for r in value['source_records']] == [5379, 24, 3]
    assert len(value['all130_following_repair_records']) == 130
    assert value['complete_output_probability'] == pytest.approx(1, abs=1e-12)
    tree, copied, phase = source.runs()
    assert tree.bits == 8 and copied.bits == phase.bits == 16
    for run in (tree, copied, phase):
        event = int(np.argmax(run.failure_weights.sum(axis=1)))
        branches = run.failure(event)
        np.testing.assert_allclose([np.vdot(x, x).real for x in branches],
                                   run.failure_weights[event, :len(branches)], rtol=1e-11, atol=1e-42)


def test_exact_source_precision_includes_tree_and_actual_literal_controls():
    value = source.certificate()
    assert F(value['complete_source_state_error_bound']) < F(848, 10**15)
    assert value['signed_amplitude_tree_nodes_checked'] == 255
    assert value['native_field_wait_or_full_detector_executed_here'] is False
    assert value['charge_postselection_used'] is False


def test_complete_cube_executor_retains_a_nonzero_source_failure(monkeypatch):
    full = np.eye(6, dtype=complex)
    full[np.ix_([0, 2], [0, 2])] = [[np.cos(.2), -np.sin(.2)], [np.sin(.2), np.cos(.2)]]
    monkeypatch.setattr(compiler, 'native_pulse', lambda *_: full.copy())
    blank = np.zeros(65536); blank[0] = 1
    value = source.complete16([compiler.pulse('I', 0., [0])], blank)
    assert value.failure_probability == pytest.approx(np.sin(.2)**2)
    assert len(value.failure(0)) == 4
    assert np.linalg.norm(value.success)**2+value.failure_probability == pytest.approx(1)


@pytest.mark.parametrize('bad', [np.zeros(32768), np.zeros((65536, 5)), np.full(65536, np.nan)])
def test_short_cube_source_rejects_an_incomplete_input(bad):
    with pytest.raises(ValueError): source.complete16([], bad)


def test_reference_phase_signs_are_essential_for_protection():
    from . import cube
    unsigned = np.abs(source.ideal_reference())
    residual = cube.operators().H.real@unsigned-48*unsigned
    assert np.linalg.norm(residual) > 1
