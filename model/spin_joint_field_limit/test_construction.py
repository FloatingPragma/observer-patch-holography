from fractions import Fraction as F

import pytest

from model.native_spin_observer import bounds as native_bounds
from . import construction


def test_cutoff_safe_against_exact_inherited_binary_bound():
    for radius in (3,10,30):
        for time in (F(),F(1,100),F(1,8)):
            actual=native_bounds.history_cutoff(radius,1,2,F(1),time,2)
            upper=construction.history_cutoff(radius,1,2,F(1),time,2)
            assert actual<=upper<=2


def test_finite_prepared_complete_HB_history_budget():
    row=construction.probe_history(2,(('X',(0,0,0)),('Z',(1,0,0))),
                                   (F(1,32),F(1,16),F()),F(1,10))
    assert row['complete_output_half_trace_upper']<=F(1,10)
    assert row['ideal_pointer_record_count']==4
    assert row['all_native_failure_and_pointer_records_retained']
    assert row['retained_field_clock']['selection_field']==row['field']
    assert row['prepared_source']['field']==row['field']
    assert not row['effective_continuum_history_accuracy_included']
    assert not row['large_history_executed']


def test_zero_duration_history_needs_no_clock_pulses():
    row=construction.probe_history(2,(('Z',(0,0,0)),),(F(),F()),F(1,10))
    assert row['retained_field_clock']['native_pulse_count_upper']==0
    assert row['two_clock_cutoffs_combined_half_trace_upper']==0
    assert row['ideal_pointer_record_count']==2


@pytest.mark.parametrize('probes,times',[
    ((),(0,)),((('bad',(0,0,0)),),(0,0)),
    ((('X',(0,0)),),(0,0)),((('X',(0,0,0)),),(0,)),
    ((('X',(0,0,0)),),(0,0.0))])
def test_incomplete_or_inexact_history_rejected(probes,times):
    with pytest.raises(ValueError):construction.probe_history(2,probes,times,F(1,10))
