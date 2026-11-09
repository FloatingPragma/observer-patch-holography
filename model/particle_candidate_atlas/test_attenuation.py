from fractions import Fraction as F
import copy
import numpy as np
import pytest
import sympy as s
from . import attenuation as a

@pytest.mark.parametrize('r',['0','3/5','1'])
@pytest.mark.parametrize('z',[s.Integer(1),s.Rational(99,101)+s.I*s.Rational(20,101)])
def test_complete_blank_instrument_and_original_weights(r,z):
    row=a.family(r,z);ks=a.kraus(r,z);rr=s.Rational(r)
    assert a.clean(sum((k.H*k for k in ks),s.zeros(3)))==s.eye(3)
    assert a.clean(ks[0].H*ks[0]-row['Y'])==s.zeros(3)
    assert row['positive_roots']
    assert row['signed_Q']==(1+rr*rr)/3
    assert sum(row['original_tracial_records'])==1
    assert sum(row['original_tracial_records'][8:])==5*(1-rr*rr)/36
    assert row['original_tracial_records'][3:8]==(s.Rational(1,36),s.Rational(25,36),0,0,0)
    assert s.trace(row['Y'])==5*(1+rr*rr)/36
    for flag,k in enumerate(ks):assert a.analyzer_block(flag,0,r,z)[:3,:3]==k


def test_loss_flags_are_necessary_for_full_input_gram():
    ks=a.kraus('3/5')
    assert a.clean(sum((k.H*k for k in ks[:8]),s.zeros(3)))!=s.eye(3)
    # A chosen three-dimensional reference is retained, with no normalization.
    psi=s.Matrix([1,0,0,0,s.I,0,0,0,2])/s.sqrt(6)
    out=[s.kronecker_product(k,s.eye(3))*psi for k in ks]
    assert s.simplify(sum((v.H*v)[0] for v in out))==1
    assert s.simplify(sum((v.H*v)[0] for v in out[8:]))>0


def test_all_64_input_columns_factor_through_unitaries():
    r=s.Rational(3,5);z=s.Rational(99,101)+s.I*s.Rational(20,101)
    prep,f,blocks,loss=a._factors(r,z)
    assert a.clean(prep.H*prep)==s.eye(8)
    assert a.clean(f.H*f)==s.eye(8)
    assert a.clean(loss.H*loss)==s.eye(16)
    assert all(a.clean(b.H*b)==s.eye(4)for b in blocks)
    # Independent numeric multiplication checks the public block indexing on every column.
    control=s.diag(*(blocks+blocks))
    whole=s.kronecker_product(s.eye(2),f,s.eye(4))*s.kronecker_product(loss,s.eye(4))*control*s.kronecker_product(s.eye(2),prep,s.eye(4))
    numeric=np.asarray(whole,complex)
    assert np.linalg.norm(numeric.conj().T@numeric-np.eye(64),2)<1e-12
    for out,inp in [(0,0),(8,0),(1,9),(9,1),(15,15)]:
        assert np.max(abs(numeric[4*out:4*out+4,4*inp:4*inp+4]-np.asarray(a.analyzer_block(out,inp,str(r),z),complex)))<1e-12


def test_native_loss_word_all_columns_and_failure_slots():
    word=a.attenuator_word('3/5')
    actual=a.compiler.execute(word,4)
    _,_,_,loss=a._factors(s.Rational(3,5),s.Integer(1))
    assert np.linalg.norm(actual-np.asarray(loss,complex),2)<2e-11
    assert np.linalg.norm(actual.conj().T@actual-np.eye(16),2)<2e-11
    assert sum(6-(4 if len(event['targets'])==2 else 2)for event in word)>0
    assert len(word)==len(a.attenuator_word('1'))==len(a.attenuator_word('0'))
    # Compiler contracts only exactly invariant proper-code columns; every export
    # complement remains part of its instrument, with zero weight on this code.


def test_complete_weak_projection_has_no_loss_postselection():
    ks=a.weak_kraus('3/5',1,'1/4')
    assert len(ks)==16
    assert a.clean(sum((k.H*k for k in ks),s.zeros(6)))==s.eye(6)
    row=a.family('3/5');y=row['Y'];co=s.Rational(15,17);si=s.Rational(8,17)
    x=s.Matrix([[0,1],[1,0]])
    assert ks[0]==a.clean(s.eye(6)+(co-1)*s.kronecker_product(s.eye(2),y)-s.I*si*s.kronecker_product(x,y))
    nonzero=a.clean(sum((k.H*k for k in ks[1:]),s.zeros(6)))
    assert nonzero==a.clean(2*(1-co)*s.kronecker_product(s.eye(2),y-y*y))


def test_finite_controller_counterfamily_preserves_weights_and_changes_shape():
    balanced=a.family('1');damped=a.family('3/5')
    assert balanced['signed_Q']==s.Rational(2,3)
    assert damped['signed_Q']==s.Rational(34,75)
    assert balanced['signed_Q']!=damped['signed_Q']
    assert sum(balanced['original_tracial_records'])==sum(damped['original_tracial_records'])==1
    assert a.chi_from_declared_u('1/25')==F(77,24627200)


@pytest.mark.parametrize('r',[True,0.5,'-1/10','11/10'])
def test_invalid_pass_amplitude_is_rejected(r):
    with pytest.raises(ValueError):a.family(r)


@pytest.mark.parametrize('u',['-1','77',True])
def test_invalid_declared_response_domain_is_rejected(u):
    with pytest.raises(ValueError):a.chi_from_declared_u(u)


def test_archive_bounds_and_mode_counts_remain_scoped():
    b=a.budget(4096,'1/16',round_isometry_error='1/1000000')
    assert F(b['full_CAR_complete_unhalved_error'])==108*F(1,16)**2/4096+2*4096*F(1,1000000)
    assert b['zero_outcome_postselection']is False
    assert b['fixed_native_service_clock_equals_weak_clock']is False
    d=a.certificate();before=copy.deepcopy(d)
    d['rows'][0]['Q']='forged'
    assert a.certificate()==before
    assert before['full_CAR_with_fresh_flags_modes']==1440
    assert before['chi_selected_from_primitive_source']is False


def test_complete_native_extended_analyzer():
    receipt=a.native_execution('3/5',s.Rational(99,101)+s.I*s.Rational(20,101))
    assert receipt['full_input_columns']==64
    assert receipt['native_full_analyzer_execution_passed']is True
    assert receipt['indexed_failure_slots']>receipt['pulses']
    assert receipt['global_phase_fitted']is False
    assert receipt['certified_round_dilation_error']is None


def test_native_round_has_all_128_columns_and_16_return_records():
    d=a.native_round_execution('3/5',1,'1/64')
    assert d['full_input_columns']==128 and d['retained_return_flags']==16
    assert sum(s.sympify(r['exact_tracial_probability'])for r in d['original_return_records'])==1
    assert sum(s.sympify(r['exact_tracial_probability'])for r in d['original_return_records'][1:])>0
    assert d['indexed_failure_slots']>2*d['pulses']
    assert d['certified_round_dilation_error']is None
    assert d['global_phase_fitted']is False
