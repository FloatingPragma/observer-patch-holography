from fractions import Fraction as F
import math

import numpy as np

from model.cube_eta_protection import cube
from . import experiment as e


def test_full_native_ground_experiment_matches_the_source_spectral_intervals(ground_readout):
    value,_=ground_readout
    result=value['executed_native_readout']
    for row in result['settings']:
        low,high=map(F,row['unconditioned_pointer_zero_interval'])
        assert float(low)<=row['all_six_pointer_probabilities'][0]<=float(high)
        assert row['complete_output_error_diagnostic']<.0006
        assert abs(row['complete_probability']-1)<result['host_numerical_admission']
        assert len(row['all_six_pointer_probabilities'])==6
    assert result['settings'][1]['all_six_pointer_probabilities'][0]>.7
    assert not result['all78_lines_native_resolved']
    assert not result['success_or_charge_postselection_used']


def test_final_analysis_retains_the_physical_N6_component(ground_readout):
    _,arrays=ground_readout
    masks=cube.weights(6)
    for j in (0,1):
        state=arrays[f'analyzer_{j}_full_success_state'].reshape(2,65536)
        mass=float(np.vdot(state[:,masks],state[:,masks]).real)
        assert mass>1e-6


def test_complete_repair_records_and_individually_addressable_states(ground_readout):
    value,_=ground_readout
    records=value['executed_native_readout']['all130_repair_outcome_slots']
    assert len(records)==130
    assert sum(row['record']==['outside'] for row in records)==10
    for law,label in [('C',('native',0,0)),('E',('native',1,7)),('C',('outside',))]:
        state=e.repair_output(law,0,label)
        actual=float(np.vdot(state,state).real)
        row=next(r for r in records if r['law']==law and r['seam']==records[0]['seam'] and r['record']==list(label))
        assert math.isclose(actual,row['probability'],rel_tol=1e-8,abs_tol=1e-42)


def test_zero_time_native_analyzers_have_expected_marginals(ground_readout):
    value,_=ground_readout
    rows=value['executed_native_readout']['zero_time_native_control']
    assert abs(rows[0]['all_six_pointer_probabilities'][0]-1)<1e-10
    assert abs(rows[1]['all_six_pointer_probabilities'][0]-.5)<1e-10
    assert all(len(row['all_six_pointer_probabilities'])==6 for row in rows)
