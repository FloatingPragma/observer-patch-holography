from fractions import Fraction as F

import mpmath as mp
import numpy as np
import pytest

from . import instrument as i, composition as c, sensitivity as s
from .certificate import report, symbolic_checks
from tools import strict_json


def trace_distance(a,b):
    return np.linalg.svd(a-b,compute_uv=False).sum()/2


def entangled_state():
    # Register/reference entanglement includes the full current gap two.
    vector=np.zeros(6,dtype=complex);vector[0]=1/np.sqrt(2);vector[5]=1j/np.sqrt(2)
    return np.outer(vector,vector.conj())


@pytest.fixture(scope='module')
def instruments():
    return {'weak':i.weak_instrument(.17),'mass':i.mass_instrument(2/9,.17)}


def test_receipt_and_symbolic_operator_identities():
    assert all(symbolic_checks().values())
    from .certificate import HERE
    assert report()==strict_json.load(HERE/'results.json')


def test_four_source_projections_reconstruct_and_are_orthogonal(instruments):
    for parts in instruments.values():
        total=sum(parts.values())
        for a in i.LABELS:
            for b in i.LABELS:
                overlap=np.sum(parts[a].conj()*parts[b],axis=(0,1))
                if a!=b:assert np.max(abs(overlap))<2e-12
        assert np.max(abs(np.sum(abs(total)**2,axis=(0,1))-1))<2e-12


def test_all_outcomes_preserve_entangled_reference_without_postselection(instruments):
    rho=entangled_state()
    for parts in instruments.values():
        outputs=i.apply_with_reference(rho,i.gram_parts(parts))
        total=sum(outputs.values())
        assert abs(np.trace(total)-1)<2e-12
        assert all(np.linalg.eigvalsh((out+out.conj().T)/2).min()>-2e-12 for out in outputs.values())
        assert np.max(abs(np.trace(total.reshape(3,2,3,2),axis1=0,axis2=2)-np.eye(2)/2))<2e-12
        # A success-only map loses original mass and alters the reference marginal.
        assert np.trace(outputs['yes_yes']).real<1-1e-7


def test_fresh_swap_exports_the_unmodified_failure_and_its_correlations(instruments):
    for parts in instruments.values():
        for label,source in parts.items():
            branch=i.reset_branch(label,source)
            assert np.array_equal(branch.archived_source_register,source)
            assert not np.shares_memory(branch.archived_source_register,source)
            assert np.array_equal(branch.fresh_source,i.preparation())
            assert branch.record==label


def test_nonreturn_environment_has_rank_three_not_one():
    x,e,_=i.analyzer_data(2/9)
    parts=i.return_parts(x)
    gram=sum(i.gram_parts(parts)[name] for name in i.LABELS if name!='yes_yes')
    assert np.max(abs(gram-(np.diag(e)-np.outer(e,e))))<2e-12
    assert np.linalg.eigvalsh(gram).min()>1e-6
    false_single_kraus=np.outer(np.sqrt(np.diag(gram).real),np.sqrt(np.diag(gram).real))
    assert np.max(abs(false_single_kraus-gram))>1e-3


def test_complete_histories_compose_and_omitting_failures_changes_law(instruments):
    for parts in instruments.values():
        grams=i.gram_parts(parts)
        histories=c.history_multipliers(grams,3)
        assert len(histories)==64
        total=c.complete_multiplier(grams,3)
        assert np.max(abs(sum(histories.values())-total))<2e-12
        assert np.max(abs(grams['yes_yes']**3-total))>1e-6


def test_preparation_gate_is_an_actual_unitary_map():
    H=i.address_preparation_unitary()
    assert np.max(abs(H.T@H-np.eye(60)))<2e-13
    assert np.max(abs(H[:,0]-np.ones(60)/np.sqrt(60)))<2e-13


def test_full_current_span_is_included_in_the_reference_stable_error_bound():
    rho=entangled_state()
    current=np.array([-1,0,1]);delta=current[:,None]-current[None,:]
    ideal=np.exp(-1j*float(i.P)*delta)
    for n in (1,10,100,1000):
        actual=i.unconditional_weak_multiplier(1/n)**n
        actual_state=(actual[:,None,:,None]*rho.reshape(3,2,3,2)).reshape(6,6)
        ideal_state=(ideal[:,None,:,None]*rho.reshape(3,2,3,2)).reshape(6,6)
        bound=float(F(c.error_budget(n)['ideal_half_diamond_bound']))
        error=trace_distance(actual_state,ideal_state)
        assert error<=bound+1e-12
        if n==1000:assert error/bound>.999


def test_mass_limit_bound_on_entangled_inputs_and_stationary_family():
    rho=entangled_state();_,e,_=i.analyzer_data(2/9)
    ideal=np.exp(-1j*(e[:,None]-e[None,:]))
    for n in (1,10,100):
        multiplier=i.unconditional_mass_multiplier(2/9,1/n)**n
        actual=(multiplier[:,None,:,None]*rho.reshape(3,2,3,2)).reshape(6,6)
        target=(ideal[:,None,:,None]*rho.reshape(3,2,3,2)).reshape(6,6)
        assert trace_distance(actual,target)<=float(F(c.error_budget(n,cycle='mass')['ideal_half_diamond_bound']))+1e-12
        assert np.array_equal(c.fixed_entry_mask(i.unconditional_mass_multiplier(2/9,1/n)),np.eye(3,dtype=bool))
    for probabilities in ([1,0,0],[1/3]*3,[.2,.3,.5]):
        state=np.diag(probabilities)
        assert np.max(abs(i.unconditional_mass_multiplier(2/9,.1)*state-state))<1e-13


def test_resources_and_implementation_allowance_are_explicit():
    row=c.resources(1000)
    assert row['fresh_prepared_pairs_including_terminal_active_pair']==1001
    assert row['classical_record_bits']==2000 and row['exported_address_pairs']==1000
    assert row['formal_complete_history_count']=='4^1000'
    for kind,n in [('weak',172840),('mass',269291)]:
        assert c.required_rounds(F(1,10**6),cycle=kind)==n
        design=c.implementation_design(F(1,10**6),cycle=kind)
        assert F(design['certified_total_error'])<=F(1,10**6)
        assert F(design['maximum_combined_channel_and_fresh_pair_error_per_cycle'])>0
    assert F(c.error_budget(100,per_cycle_channel_error=F(1,10**5))['implemented_half_diamond_bound'])>F(c.error_budget(100)['ideal_half_diamond_bound'])


def test_sensitivity_matches_independent_differentiation_and_scale_cancels():
    with mp.workdps(50):
        theta=mp.mpf(2)/9;balance=mp.sqrt(2);jac=s.logarithmic_jacobian(theta,balance)
        for k in (0,1):
            assert abs(jac[k,0]-mp.diff(lambda t:mp.log(s.ratios(t,balance)[k]),theta))<mp.mpf('1e-40')
            assert abs(jac[k,1]-mp.diff(lambda b:mp.log(s.ratios(theta,b)[k]),balance))<mp.mpf('1e-40')
        assert abs(mp.det(jac))>200
        for normalization in (mp.mpf(1),mp.mpf(2),mp.mpf('0.01')):
            assert max(abs(a-b) for a,b in zip(s.ratios(theta,balance),s.ratios(theta,balance,normalization)))<mp.mpf('1e-40')


def test_fixed_marginals_do_not_determine_independent_exclusive_count():
    forward,written=F(1,12),F(1,6)
    assert forward+written-2*forward*written==F(2,9)
    for joint,expected in [(F(0),F(1,4)),(F(1,12),F(1,12))]:
        probabilities=(joint,forward-joint,written-joint,1-forward-written+joint)
        assert min(probabilities)>=0 and sum(probabilities)==1
        assert forward+written-2*joint==expected!=F(2,9)


def test_finite_complete_channel_has_no_single_phase_for_both_current_gaps():
    with mp.workdps(40):
        assert s.finite_phase(100,1)!=s.finite_phase(100,2)
        assert abs(s.finite_phase(10000,1)-mp.mpf(2)/9)<abs(s.finite_phase(100,1)-mp.mpf(2)/9)
    comparison=s.report()['comparison_only']
    assert 9.83<float(comparison['relative_residual_ppm'])<9.84
    assert comparison['phase_displacement_applied'] is False
    assert comparison['reference_used_to_select_parameters'] is False


@pytest.mark.parametrize('bad', [0,-1,True,1.5])
def test_invalid_round_count_rejected(bad):
    with pytest.raises(ValueError):c.resources(bad)


@pytest.mark.parametrize('bad', [float('nan'),float('inf'),True])
def test_invalid_kick_rejected(bad):
    with pytest.raises(ValueError):i.weak_instrument(bad)
