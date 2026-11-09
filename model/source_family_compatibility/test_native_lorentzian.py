"""Full-source and complete-record controls for the native stopped clock."""
from fractions import Fraction as F

import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from model.native_spin_observer import execute
from . import native_lorentzian as n


@pytest.fixture(scope='module')
def receipt():
    return n.certificate()


def test_actual_full_source_and_prepared_direction():
    h,rho,a,terms,reference,radius=n._source()
    assert h.shape==rho.shape==a.shape==(16,16)
    assert sum((s.Rational(c.numerator,c.denominator)*n.tensor(p)for p,c in terms),s.zeros(16))==h
    assert a*a.H+a.H*a==s.eye(16)
    exact=s.zeros(16,1)
    exact[3]=exact[12]=1;exact[6]=4+s.sqrt(17);exact[9]=-exact[6]
    assert (h*exact-(4-s.sqrt(17))*exact/8).applyfunc(s.simplify)==s.zeros(16,1)
    assert (exact*exact.H/(exact.H*exact)[0]-rho).applyfunc(s.simplify)==s.zeros(16)
    lo,hi=n.root_interval()
    assert lo*lo<17<hi*hi
    assert reference==(F(1,2)-(lo+hi)/16)
    assert radius==(hi-lo)/16


def test_native_preparation_matches_actual_ground_without_replacing_it(receipt):
    run=n._execute();_,rho,_,_,_,_=n._source()
    vector=run['preparation'].success[:,0]
    state=np.outer(vector,vector.conj())
    assert np.linalg.norm(state-np.array(rho,complex),ord=1)<1e-10
    assert run['preparation'].bits==4
    assert run['preparation'].initial.shape==(16,1)
    assert np.count_nonzero(run['preparation'].initial)==1
    assert run['preparation'].native_failure_slots==200
    assert receipt['postselection_or_output_renormalization_used'] is False


def test_complete_native_charge_map_on_every_source_column():
    c=n._controls();a=np.array(n._source()[2],complex)
    initial=np.zeros((32,16),complex);initial[:16]=np.eye(16)
    run=execute.complete(c['charge'],5,initial)
    ideal=np.vstack((a.conj().T,a))
    assert np.linalg.norm(run.success-ideal,2)<2e-12
    embed=np.zeros((6,2));embed[:2]=np.eye(2)
    raw_outputs=[np.einsum('c,cir->ir',row,run.success.reshape(2,16,16))for row in embed]
    assert all(not x.any()for x in raw_outputs[2:])
    assert np.array_equal(raw_outputs[0],run.success[:16])
    assert np.array_equal(raw_outputs[1],run.success[16:])
    gram=run.success.conj().T@run.success
    for j,count in enumerate(n._physical_slots(c['charge'])):
        for slot in range(count):
            failed=execute.failure(c['charge'],run,j,slot)
            gram+=failed.conj().T@failed
    assert np.linalg.norm(gram-np.eye(16),2)<2e-12
    assert run.initial.shape[1]==16


def test_native_controlled_clock_preserves_relative_ground_phase():
    c=n._controls();h,_,_,_,reference,_=n._source()
    full=np.array(h,complex)-float(reference)*np.eye(16)
    ideal=np.eye(32,dtype=complex)
    ideal[16:,16:]=expm(-1j*float(n.STEP)*full)
    run=execute.complete(c['clock'],5,np.eye(32,dtype=complex))
    assert np.linalg.norm(run.success-ideal,2)<float(c['clock_step_error'])+2e-12
    # The split controlled terms need not cancel exactly on clock zero;
    # the same full product-formula bound covers that block.
    assert np.linalg.norm(run.success[:16,:16]-np.eye(16),2)<float(c['clock_step_error'])+2e-12
    wrong=np.eye(32,dtype=complex)
    wrong[16:,16:]=expm(-1j*float(n.STEP)*np.array(h,complex))
    assert np.linalg.norm(run.success-wrong,2)>1e-5
    assert run.native_failure_slots==1408


def _failure_gram(word,run):
    gram=np.zeros((run.initial.shape[1],run.initial.shape[1]),complex)
    for j,count in enumerate(n._physical_slots(word)):
        for slot in range(count):
            output=execute.failure(word,run,j,slot)
            gram+=output.conj().T@output
    return gram


def test_three_stage_global_instrument_on_all_input_and_reference_columns():
    """Complete Gram equality includes stopped branches, exports and dummy."""
    c=n._controls();initial=np.zeros((32,16),complex);initial[:16]=np.eye(16)
    charge=execute.complete(c['charge'],5,initial)
    gram=_failure_gram(c['charge'],charge)
    for q in range(2):
        tagged=np.zeros((32,16),complex);tagged[:16]=charge.success[16*q:16*(q+1)]
        initialize=execute.complete(c['hadamard'],5,tagged)
        gram+=_failure_gram(c['hadamard'],initialize)
        current=initialize.success
        for _ in range(3):
            selector_failure=sum(float(np.linalg.norm(x[:,0])**2)for x in c['coin'].failures[0])
            gram+=selector_failure*(current.conj().T@current)
            read=execute.complete(c['hadamard'],5,c['coin'].success[0,0]*current)
            gram+=_failure_gram(c['hadamard'],read)
            for b in range(6):
                # The complete service absorbs every code-complement export.
                # Its surviving pointer is embedded in the first two levels.
                embed=np.zeros((6,2));embed[:2]=np.eye(2)
                output=np.einsum('c,cir->ir',embed[b],read.success.reshape(2,16,16))
                if b>=2:assert not output.any()
                gram+=output.conj().T@output
            step=execute.complete(c['clock'],5,c['coin'].success[1,0]*current)
            gram+=_failure_gram(c['clock'],step);current=step.success
        gram+=current.conj().T@current
    assert np.linalg.norm(gram-np.eye(16),2)<2e-10
    assert np.linalg.eigvalsh(gram).min()>1-2e-10


def test_exact_geometric_error_uses_original_reach_weights():
    errors=n.complete_error();c=n._controls();tail=n.RATIO**n.COUNT
    selectors=(1-tail)/(1-n.RATIO);clocks=n.RATIO*selectors
    _,prep,_=n._preparation()
    expected=prep+c['charge_error']+c['hadamard_error']+selectors*c['coin_error']+clocks*c['clock_step_error']+(1-tail)*c['hadamard_error']
    assert expected==errors['actual_rational_error']<F(1,100000)
    assert clocks<256 and selectors<257
    assert n.COUNT*c['clock_step_error']>6*expected
    assert errors['actual_failure_error']<F(1,10**20)


def test_rational_phase_intervals_enclose_independent_high_precision_score():
    rows,score,dummy=n._ideal_records()
    with mp.workdps(90):
        root=mp.sqrt(17);energies=(root/8-mp.mpf(1)/16,root/8+mp.mpf(1)/16)
        weights=(mp.mpf(1)/4+1/(4*root),mp.mpf(1)/4-1/(4*root))
        ratio=mp.mpf(256)/257;step=mp.mpf(1)/256
        total=mp.mpf(0)
        for j in range(n.COUNT):
            char=sum(w*mp.cos(j*step*e)for w,e in zip(weights,energies,strict=True))
            weight=(1-ratio)*ratio**j
            total+=2*weight*char
            if j in(0,1,77,1000,2047):
                for sign,interval in zip((1,-1),rows[j],strict=True):
                    value=weight*(mp.mpf(1)/4+sign*char/2)
                    assert mp.mpf(interval[0].numerator)/interval[0].denominator<=value<=mp.mpf(interval[1].numerator)/interval[1].denominator
        assert mp.mpf(score[0].numerator)/score[0].denominator<=total<=mp.mpf(score[1].numerator)/score[1].denominator
    assert dummy==n.RATIO**n.COUNT/2


def test_executed_complete_original_records_and_dummy_environment(receipt):
    run=n._execute();records=receipt['all_proper_terminal_records']
    assert len(records)==8194
    assert len({(x['time_index'],x['charge'],x['clock_raw'])for x in records})==8194
    assert len([x for x in records if x['signed_score']==0])==2
    assert run['proper'].shape==(2048,2,2)
    assert run['prefixes'].shape==(2049,32,2)
    assert run['prefixes'].nbytes<3*1024*1024
    assert np.array_equal(run['dummy'],np.sum(abs(run['prefixes'][-1])**2,axis=0))
    assert float(run['proper'].sum()+run['dummy'].sum()+sum(run['failures'].values()))==pytest.approx(1,abs=2e-8)
    for q in range(2):
        assert run['proper'][:,q].sum()+run['dummy'][q]==pytest.approx(.5,abs=2e-8)
    assert run['dummy'].sum()>0
    assert receipt['dummy_clock_environment_retained'] is True


def test_all_zero_raw_labels_and_failure_slots_are_retained(receipt):
    assert receipt['raw_pointer_terminal_record_count']==36*2048+6
    assert receipt['physical_first_failure_records']==17449226
    assert receipt['first_failure_records_on_two_nonzero_charge_prefixes']==5816554
    for q in range(6):
        for b in range(6):
            interval=n.terminal_record_interval(0,q,b)
            if q>=2 or b>=2:assert interval==(0,0)
            else:assert interval[1]>0
    for q in range(2,6):
        assert n.terminal_record_interval(None,q,None)==(0,0)
        for stage,j,event,slot in [('initialize',None,0,0),('selector',2047,0,3),('clock',2047,0,0),('readout',0,1,3)]:
            output=n.native_failure_record(stage,j,q,event,slot)
            assert output.ndim==2 and output.shape[1]==1
            assert not output.any()


@pytest.mark.parametrize('stage,j,q,event,slot',[
    ('preparation',None,None,0,0),('charge',None,None,0,0),
    ('initialize',None,1,1,3),('selector',2047,0,0,3),
    ('clock',0,1,0,0),('clock',2047,0,403,3),('readout',1024,1,0,2),
])
def test_failure_vector_reconstructs_retained_unnormalized_source(stage,j,q,event,slot):
    output=n.native_failure_record(stage,j,q,event,slot)
    assert output.ndim==2 and output.shape[1]==1
    assert np.isfinite(output).all()
    assert np.linalg.norm(output)**2<1e-20
    if stage in('preparation','charge'):
        assert np.linalg.norm(output)**2==pytest.approx(n._execute()[stage].failure_weights[event,slot],abs=1e-35)


@pytest.mark.parametrize('arguments',[
    ('initialize',None,True,0,0),('selector',2048,0,0,0),
    ('clock',0,6,0,0),('clock',0,0,404,0),
    ('readout',0,3,0,4),('unknown',0,0,0,0),
    ('preparation',None,None,False,0),
])
def test_invalid_record_indices_are_not_silent_zero_outcomes(arguments):
    with pytest.raises(ValueError):n.native_failure_record(*arguments)


def test_padding_is_not_a_physical_failure():
    word=n._controls()['clock']
    event=next(j for j,x in enumerate(word)if x['kind']=='J')
    with pytest.raises(ValueError):n.native_failure_record('clock',0,0,event,2)
    with pytest.raises(ValueError):n.native_failure_record('clock',0,5,event,2)


def test_certified_score_and_original_window_do_not_use_float_residual(receipt):
    assert receipt['native_signed_score_interval'][0]<F(7993,10000)<receipt['native_signed_score_interval'][1]
    assert receipt['original_total_weight_in_0_to_1_lower']>F(59,100)
    assert receipt['original_total_weight_in_0_to_1_lower']<F(3,5)
    assert receipt['floating_contraction_diagnostics_used_as_interval_premise'] is False
    assert receipt['hardware_shots_sampled'] is False
    assert receipt['large_infinite_ground_collar_executed'] is False
    assert receipt['physical_species_or_mass_assigned'] is False
    actual=n._execute()['signed_score'];interval=receipt['native_signed_score_interval']
    assert float(interval[0])<actual<float(interval[1])


def test_public_result_is_detached_even_for_mutable_fraction_internals():
    first=n.certificate();saved=repr(first)
    first['native_signed_score_interval'][0]._numerator=0
    first['all_proper_terminal_records'][0]['original_probability_interval'][0]._numerator=0
    first['native_programs']['clock']['pulse_count']=0
    second=n.certificate()
    assert repr(second)==saved
    interval=n.terminal_record_interval(2,0,0);saved_interval=repr(interval)
    interval[0]._numerator=0
    assert repr(n.terminal_record_interval(2,0,0))==saved_interval
