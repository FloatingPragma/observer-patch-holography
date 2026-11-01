from fractions import Fraction as F
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import modular as m
from . import controls


@pytest.mark.parametrize('beta',[0,-1,True,'1',float('nan'),float('inf'),9])
def test_invalid_reference_calibrations_rejected(beta):
    with pytest.raises(ValueError):m.gibbs(np.eye(16),beta)


def test_gibbs_offset_is_not_a_reference_or_modular_flow_choice():
    h=np.diag(np.linspace(-1,1,16))
    a=m.gibbs(h);b=m.gibbs(h+3*np.eye(16))
    np.testing.assert_allclose(a['rho'],b['rho'],atol=1e-14)
    np.testing.assert_allclose(a['K'],b['K'],atol=1e-14)
    assert abs(b['log_partition']-a['log_partition']+3)<1e-14


def test_singular_reference_is_rejected_without_support_clipping():
    with pytest.raises(ValueError,match='faithfulness'):m.faithful_log(np.diag([1.,0.]))


def test_quantum_relative_entropy_retains_offdiagonal_coherence():
    sigma=np.diag([.7,.3]);angle=.6
    U=np.array([[math.cos(angle),-math.sin(angle)],[math.sin(angle),math.cos(angle)]])
    rho=U@np.diag([.9,.1])@U.T
    quantum=m.relative_entropy(rho,sigma)
    populations=np.diag(rho)
    classical=float(np.sum(populations*np.log(populations/np.diag(sigma))))
    assert quantum-classical>.1


def test_bregman_remainder_removes_linear_trace_drift():
    rho=np.diag([.6,.4]);direction=.3*rho;epsilon=.02
    raw=m.relative_entropy(rho+epsilon*direction,rho)
    corrected=m.entropy_bregman_remainder(rho,direction,epsilon)
    exact=(1+.3*epsilon)*math.log(1+.3*epsilon)-.3*epsilon
    assert abs(corrected-exact)<1e-14
    assert raw>100*corrected>0


def test_exact_full_fock_boundary_obstruction():
    result=m.exact_boundary_obstruction()
    assert F(result['cross_term_Hilbert_Schmidt_norm_squared_exact'])==4
    assert result['both_partial_traces_exactly_zero']


def test_actual_reference_is_same_action_faithful_and_not_forced_to_factor():
    reference=m.reference();region=m.regional_certificate()
    assert reference['minimum_eigenvalue']>0
    np.testing.assert_allclose(reference['K'],float(m.BETA)*reference['H']+reference['log_partition']*np.eye(16),atol=1e-12)
    assert region['nested_all_matrix_units_error']<1e-12
    assert region['regional_minimum_eigenvalue']>=4*reference['minimum_eigenvalue']-1e-12
    assert region['same_beta_bare_energy_mismatch']>1e-7
    assert region['retraction_control_error']>1e-3
    assert region['modular_covariance_control_error']>1e-3
    assert abs(region['nonadditive_H_HS_squared']-4)<1e-12


def test_same_reference_observations_preserve_arbitrary_ancilla_operators():
    a,b,d=m.observation_tower(m.reference()['rho'])
    x=np.arange(32**2,dtype=float).reshape(32,32)
    x=x+1j*x.T
    np.testing.assert_allclose(b.observe(a.observe(x,2),2),d.observe(x,2),atol=5e-11,rtol=0)
    recovered=a.recover(a.target_reference)
    np.testing.assert_allclose(recovered,m.reference()['rho'],atol=1e-13)
    wrong=m.ReferenceObservation(np.eye(4)/4,m.car_restriction(2,(0,)))
    assert np.linalg.norm(wrong.observe(a.observe(x,2),2)-d.observe(x,2))>1e-4


def test_source_preparation_retains_all_purification_and_failure_outputs():
    result=m.prepare_reference();columns=result['columns'];run=result['run']
    assert columns.shape==(16,16)
    np.testing.assert_allclose(columns@columns.conj().T,m.reference()['rho'],atol=m.ADMISSION)
    Q=m.joint.charge()
    np.testing.assert_allclose(Q@columns-columns@Q,0,atol=m.ADMISSION)
    assert run.failure_weights.shape==(len(result['source_word']),4)
    assert abs(np.linalg.norm(columns)**2+run.failure_probability-1)<m.ADMISSION
    assert result['preparation_certificate']['preparation_work_bits']==8


def test_native_petz_recovery_is_complete_on_all_regional_reference_columns():
    result=m.source_recovery();run=result['run']
    assert run.success.shape==(64,4)
    np.testing.assert_allclose(run.success.conj().T@run.success+run.failure_gram,np.eye(4),atol=m.ADMISSION)
    np.testing.assert_allclose(run.success,result['target'],atol=m.ADMISSION)
    region=result['observation'].target_reference
    kraus=run.success.reshape(4,16,4)
    restored=sum(k@region@k.conj().T for k in kraus)
    np.testing.assert_allclose(restored,m.reference()['rho'],atol=m.ADMISSION)
    assert run.failure_weights.shape==(len(result['word']),4)


@pytest.mark.parametrize('law,policy',[('C','coherent'),('C','flagged'),('E','coherent'),('E','flagged')])
def test_complete_actual_histories_obey_quantum_and_record_entropy_balance(law,policy):
    r=m.history_entropy(law,policy)
    assert r['outputs'].shape==(1024,16,16)
    assert r['probabilities'].shape==(1024,) and r['probabilities'].min()>=0
    assert abs(r['total_probability']-1)<m.ADMISSION
    assert max(abs(r['unread_identity_error']),abs(r['cq_identity_error']),abs(r['firstlaw_error']))<m.ADMISSION
    assert r['reference_stationarity_error']>1e-5
    assert r['noncommuting_variation_HS']>1e-5
    assert r['cq_entropy_change']>r['unread_entropy_change']
    assert -1e-10<r['relative_entropy_remainder']<=r['analytic_remainder_upper_numerical']+1e-10


def test_frozen_controls_are_only_actual_native_words_and_reject_bad_indices():
    words=controls.load()
    for word in words.values():
        rebuilt=controls.decode(controls.encode(word))
        assert rebuilt==word
    bad=controls.encode(words['preparation'][:1]);bad[0,0]=.5
    with pytest.raises(ValueError):controls.decode(bad)


def test_bad_preparation_control_cannot_preserve_the_claimed_reference(monkeypatch):
    source=controls.load();bad=[dict(e,targets=list(e['targets'])) for e in source['preparation']]
    index=next(i for i in range(len(bad)-1,-1,-1) if bad[i]['kind']=='Z' and bad[i]['targets']==[0])
    bad[index]['angle']+=.1
    monkeypatch.setattr(controls,'load',lambda:{**source,'preparation':bad})
    m.prepare_reference.cache_clear()
    try:
        with pytest.raises(ArithmeticError):m.prepare_reference()
    finally:
        m.prepare_reference.cache_clear()
