"""Complete-carrier, substitution and source-compilation negative controls."""
from copy import deepcopy
from fractions import Fraction as F
import numpy as np
import pytest
from . import contract as c,source as s,reflection as r


def test_manifest_rejects_hardware_clock_and_record_substitution():
    original=c.manifest()
    for path,value in [(('carrier','load'),1),(('carrier','event_register'),60),
                       (('policies','fixed_compensated_clock','time_per_block'),'1'),
                       (('policies','coherent_spatial','raw_records'),'discard')]:
        changed=deepcopy(original)
        target=changed
        for key in path[:-1]:target=target[key]
        target[path[-1]]=value
        with pytest.raises(ValueError,match='substitution'):c.validate_manifest(changed)
    assert c.validate_manifest(original)==original


def test_every_padded_swap_column_and_arbitrary_reference_are_retained():
    # Integer address labels provide a simultaneous independent check of all
    #262144 input columns, including invalid source/memory/event addresses.
    labels=np.arange(512*8*64,dtype=float).reshape(512,8,64,1)
    moved=s.event_swap(labels)
    assert np.array_equal(s.event_swap(moved),labels)
    assert len(np.unique(moved))==labels.size
    codes=set(c.source_addresses());changed=0
    for address in range(512):
        if address not in codes:assert np.array_equal(moved[address],labels[address])
        assert np.array_equal(moved[address,6:],labels[address,6:])
        assert np.array_equal(moved[address,:,60:],labels[address,:,60:])
    # The untouched reference is carried by the same column permutation.
    both=np.concatenate((labels,1j*(labels+1)),axis=3)
    assert np.array_equal(s.event_swap(both),np.concatenate((moved,1j*(moved+1)),axis=3))


def test_complete_terminal_failure_counts_source_memory_and_event_padding():
    state=np.zeros((512,8,64,1),complex);address=c.source_addresses()[0]
    state[address,:6,0,0]=c.memory_character(1)
    state[511,0,0,0]=1e-6;state[address,7,0,0]=2e-6;state[address,0,63,0]=3e-6
    result=s.terminal_records(state)
    assert abs(result['complement_probability']-14e-12)<1e-25
    assert abs(sum(map(sum,result['source_pair_character_probabilities']))+result['complement_probability']-result['input_norm_squared'])<2e-15
    assert not result['normalization_performed']


def test_compiled_reflection_complete_columns_sign_and_padding():
    U,R=s.compiled_event_matrices();u=np.zeros(64);u[:60]=1/np.sqrt(60)
    target=2*np.outer(u,u)-np.eye(64)
    assert np.linalg.norm(U[:,0]-u)<1e-10 and np.linalg.norm(R-target)<1e-9
    assert np.linalg.norm(R@u-u)<1e-10
    assert np.linalg.norm((-R)@u-u)>1.999
    assert np.linalg.norm(R[:,60:]+np.eye(64)[:,60:])<1e-10
    bound=r.complete_bound(len(r.plan(60)['reflection']),F(13,r.iv.AQ))
    assert np.linalg.norm(R-target)<float(F(bound['full_word_operator_error_upper']))


def test_full_state_inverse_and_nonunit_initial_norm():
    state=s.prepare();state*=.7
    final=s.spatial_step(s.spatial_step(state),inverse=True)
    assert np.linalg.norm(final-state)<2e-12
    result=s.terminal_records(final)
    assert abs(result['input_norm_squared']-.49)<1e-10
    assert not result['normalization_performed']


def test_complete_clock_failure_is_retained_without_repreparation():
    initial=s.prepare();initial[511,7,63,0]=1e-5
    record=s.observe_fixed_clock(initial,'Y')
    assert record['initial_complement_probability']>=.99999e-10
    assert abs(sum(map(sum,record['joint_probabilities']))+record['initial_complement_probability']-record['input_norm_squared'])<2e-10
    assert not record['source_repreparation_between_blocks']and not record['normalization_performed']


def test_exact_walsh_roundtrip_and_enlarged_full_columns():
    p=r.plan(120);U=r.factor_matrix(120);u=np.zeros(128);u[:120]=1/np.sqrt(120)
    assert np.linalg.norm(U[:,0]-u)<2e-11
    assert np.linalg.norm(U.T@U-np.eye(128))<2e-12
    assert r._walsh(r._walsh([2,-3,4,5]))==[8,-12,16,20]
    with pytest.raises(ValueError):r.plan(1261)


def test_actual_moving_matching_code_failure_and_dirty_workspace():
    state=np.zeros((512,8,64,1),complex)
    state[:,:6,:60,0]=c.matching_vector(0)[:,None,None]*c.memory_character(1)[None,:,None]/np.sqrt(60)
    moved=s.spatial_step(state)
    code=np.stack([c.matching_vector(j)for j in range(6)])
    projected=np.einsum('ls,smer->lmer',code.conj(),moved)
    retained=float(np.vdot(projected,projected).real)
    assert abs(retained-67/90)<2e-13 and abs(1-retained-23/90)<2e-13
    control=r.workspace_admission(np.array([[1,0],[0,1j/4]]))
    assert control['blank_mass']==1 and control['dirty_failure_mass']==1/16
    assert not control['normalization_performed']


def test_actual_original_coin_policy_and_extension_reject_substitutions():
    initial=s.prepare();rows=s.observe_committed_event(initial)
    assert rows['raw_records']==750 and len(rows['raw_record_probabilities'])==750
    assert abs(sum(rows['raw_record_probabilities'])+rows['initial_complement_probability']-rows['input_norm_squared'])<2e-10
    electric=s.electric_extension_protocol()
    assert len(electric['steps'])==3 and not any(electric['final_gauss_residual'])
    assert all(len(row['all_fifty_raw_probabilities'])==50 for row in electric['steps'])
    altered=deepcopy(electric['extension_descriptor']);altered['source_embedding']=66
    with pytest.raises(ValueError,match='substitution'):s.electric_extension_protocol(descriptor=altered)


def test_explicit_idle_mode_is_distinct_preparation_and_full_walk_eigenmode():
    result=s.idle_mode_control()
    assert result['full_W_eigenvalue_target']=='-1'
    assert result['stored_full_step_residual_diagnostic']<2e-10
    assert result['uniform_A_overlap_diagnostic']<2e-11
    assert not result['physical_particle_or_rest_energy_identification']


def test_actual_inverse_native_word_on_all_event_columns():
    _,forward=s.compiled_event_matrices();inverse=s.compiled_inverse_reflection()
    assert np.linalg.norm(inverse@forward-np.eye(64))<2e-10
    state=s.prepare();moved=s.spatial_step(state,stored=True)
    returned=s.spatial_step(moved,inverse=True,stored=True)
    assert np.linalg.norm(returned-state)<2e-10
    assert abs(np.vdot(returned,returned).real-np.vdot(state,state).real)<2e-10
