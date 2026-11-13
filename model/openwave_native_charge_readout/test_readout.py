"""Independent controls on full padding, phase observations and clock recovery."""
from fractions import Fraction as F
import numpy as np
from . import source as s


def test_full_feedback_padding_and_compiled_inverse():
    permutation=s.shift_permutation((0,5),1)
    program=s.architecture.compile_permutation(permutation)
    actual=s.clean_circuit_columns(program)
    assert actual==permutation
    inverse={**program,'gates':list(reversed(program['gates']))}
    restored=s.architecture.execute_basis(inverse,np.array(actual,dtype=np.int64))
    assert np.array_equal(restored,np.arange(s.SIZE))
    active={s.source_address(pair)for pair in s.charge.PAIRS}
    for x,y in enumerate(actual):
        if x%s.SOURCE_DIM not in active or x//s.SOURCE_DIM>=6:assert x==y


def test_actual_source_preparation_quadratures_and_altered_voltage():
    coupled=s.quadrature_records(True);control=s.quadrature_records(False)
    assert coupled['unconditional_label_record_covariance_exact']=='3/8'
    assert control['unconditional_label_record_covariance_exact']=='0'
    for p in (coupled,control):
        for row in p['complete_six_character_experiments']:
            assert len(row['all_125_raw_native_readout_slots'])==125
            assert abs(sum(x['joint_probability']for x in row['all_125_raw_native_readout_slots'])+row['initial_admission_failure_probability']-row['total_record_weight'])<2e-14
            assert not row['numerical_state_renormalized']
            assert row['moving_output_scratch_or_memory_failure_probability']>=0
            assert row['character_preparation']['all_nonblank_inputs_retained_as_failure']


def test_clock_population_preservation_does_not_imply_coherence_preservation():
    p=s.uncompensated_clock_leakage();w=p['nonzero_charge_witness']
    assert w['coherence_multiplier_real_exact']=='85/96'
    assert w['logical_labels']==[0,2]and w['charge_character']==1
    assert w['matching_code_population_leakage_exact']=='0'
    for row in p['complete_uncompensated_one_event_transfer']:
        assert all(F(x)==0 for r in row['leakage_probability_by_charge_and_logical_label_exact']for x in r)
        for r in row['logical_coherence_multipliers_real_plus_i_sqrt3_exact'][0]:
            assert all(x==['1','0']for x in r)
    assert not p['uncompensated_clock_theorem_transferred']


def test_full_actual_witness_execution_retains_dirty_scratch():
    actual,report=s.executed_witness_unitary(True)
    assert actual.shape==(32,32)
    assert report['complete_clean_and_dirty_basis_columns']==32
    assert report['raw_source_execution_error_frobenius']<2e-8
    assert np.linalg.norm(actual.conj().T@actual-np.eye(32))<2e-8


def test_numerical_preparation_and_padding_amplitudes_enter_failure_records(monkeypatch):
    source=s.prepared_source_and_export()
    changed={**source,'prep_actual':source['prep_actual'].copy()}
    outside=next(i for i in range(512)if i not in [s.source_address(p)for p in ((0,1),(0,6))])
    changed['prep_actual'][outside]+=1e-6
    original_character=s.prepared_character
    original_witness=s.executed_witness_unitary
    def character(k):
        value,report=original_character(k);value=value.copy();value[6]+=1e-6
        return value,report
    def witness(extra=True):
        U,report=original_witness(extra);U=U.copy()
        theta=1e-6;old=U[[0,12]].copy()
        U[0]=np.cos(theta)*old[0]-np.sin(theta)*old[1]
        U[12]=np.sin(theta)*old[0]+np.cos(theta)*old[1]
        return U,report
    monkeypatch.setattr(s,'prepared_source_and_export',lambda:changed)
    monkeypatch.setattr(s,'prepared_character',character)
    monkeypatch.setattr(s,'executed_witness_unitary',witness)
    result=s.quadrature_records(True)
    for row in result['complete_six_character_experiments']:
        assert row['initial_admission_failure_probability']>1e-12
        assert row['moving_output_scratch_or_memory_failure_probability']>1e-15
        recorded=sum(x['joint_probability']for x in row['all_125_raw_native_readout_slots'])
        assert abs(recorded+row['initial_admission_failure_probability']-row['total_record_weight'])<2e-14
