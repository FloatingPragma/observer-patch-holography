"""Independent fair-source, complete-filter and reference controls."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import measurement as m, spectral as s, packet, reference


def test_original_fair_source_basis_and_exact_spectral_certificate():
    counts=np.array(s.source_counts(),dtype=object)
    assert counts.shape==(66,66)
    assert np.array_equal(counts,counts.T)
    assert all(sum(row)==60 for row in counts)
    assert s.certify()['counted_energy_operator_error_upper']=='27/100000000'


@pytest.mark.parametrize('mutation',['eigenvalue','column','drop'])
def test_bad_fair_spectral_proposals_rejected(mutation):
    columns,values=s.load()
    if mutation=='eigenvalue':values=(F(1),)+values[1:]
    elif mutation=='column':columns=((F(1),)+columns[0][1:],)+columns[1:]
    else:columns=columns[:-1]
    with pytest.raises((ValueError,ArithmeticError)):s.certify(columns,values)


def test_ideal_filters_are_complete_with_arbitrary_untouched_reference():
    rng=np.random.default_rng(1140)
    columns=rng.normal(size=(66,3))+1j*rng.normal(size=(66,3))
    columns/=np.linalg.norm(columns)
    effects,residual=m.ideal_probabilities(columns)
    expected=columns.T@columns.conj()
    observed=sum(effects)+residual.T@residual.conj()
    assert np.linalg.norm(expected-observed)<2e-14
    assert all(np.linalg.eigvalsh(block).min()>-1e-14 for block in effects)
    assert np.linalg.norm(residual)<1e-9


def test_one_word_needs_no_survival_zero_correction():
    rng=np.random.default_rng(66)
    ray=rng.normal(size=66);ray/=np.linalg.norm(ray)
    embedding=np.zeros((512,66));embedding[list(m.syn.tomo.codes())]=np.eye(66)
    padding=next(i for i in range(512) if i not in m.syn.tomo.codes())
    z=embedding@ray;z[padding]=1;z/=np.sqrt(2)
    operator=np.eye(512)-2*np.outer(z,z)
    residual=embedding.T@operator@embedding
    assert np.linalg.norm(residual-(np.eye(66)-np.outer(ray,ray)))<1e-14
    assert np.linalg.norm((operator@embedding)[padding]+ray)<1e-14
    assert np.linalg.norm(residual.T@residual+np.outer(ray,ray)-np.eye(66))<1e-14


def test_true_generator_energy_resolution_without_eigenvalue_gaps():
    columns,values=s.proposed_arrays()
    vals,basis=np.linalg.eigh(np.array(s.source_counts(),float)/60)
    gramvals,gramvecs=np.linalg.eigh(columns.T@columns)
    polar=columns@(gramvecs/np.sqrt(gramvals))@gramvecs.T
    measured=-60*np.log(values);actual=-60*np.log(vals)
    overlap=abs(polar.T@basis)**2
    error=sum(overlap[j]*(measured[j]-actual)**2 for j in range(66))
    assert max(error)<float(F(27,10**8))**2
    assert m.resolution(F(1,1000))['nearby_spectral_bad_label_probability_upper']=='729/10000000000'


def test_fair_preparation_uses_original_records_and_weights():
    for initial in range(66):
        run=m.bounds.detector.run(initial,tuple(i%2 for i in range(11)))
        weight,mapping,records=m.preparation_branch(run['slots'])
        assert weight==F(1,2**11)
        assert mapping=={initial:run['final']}
        assert len(records)==11
    _,mapping,records=m.preparation_branch((0,)*11)
    assert mapping=={} and len(records)==11


def test_conditional_reset_maps_every_terminal_to_same_source_configuration():
    for terminal in range(66):
        difference=m.syn.tomo.codes()[terminal]^m.syn.tomo.codes()[0]
        operations=tuple(('X',bit) for bit in range(9) if difference>>bit&1)
        assert m.syn.tomo.permutation(operations)[m.syn.tomo.codes()[terminal]]==m.syn.tomo.codes()[0]
        assert len(m.reset_word(terminal))<=18


def test_exact_noncentral_preparation_join_keeps_rejection_and_residual():
    value=m.prepared_reference_histogram()
    probabilities=tuple(F(x) for x in value['complete_exact_ideal_probabilities'])
    assert len(probabilities)==68 and sum(probabilities)==1 and min(probabilities)>=0
    assert probabilities[0]==F(57765,65536)
    assert value['preparation_acceptance']=='7771/65536'
    assert value['reference_partition']==409
    assert value['repair_resources']['attempts_upper']==71
    weights=np.array(value['reference_weights_in_source_order'],float)/409
    blocks,residual=m.ideal_probabilities(np.diag(np.sqrt(weights)))
    expected=[float(F(7771,65536))*np.trace(block).real for block in blocks]
    assert np.max(abs(np.array(list(map(float,probabilities[1:-1])))-expected))<1e-13
    assert abs(float(probabilities[-1])-float(F(7771,65536))*np.linalg.norm(residual)**2)<1e-20


def test_histogram_sampling_and_resources_have_separate_scopes():
    budget=m.budget();sampling=m.sampling_budget(preparation_rejection=True)
    assert budget['native_pulses']==3853872
    assert budget['native_words']==66
    assert budget['maximum_literal_nodes_per_word']<=budget['word_admission_ceiling']
    assert budget['total_streamed_word_nodes']>budget['word_admission_ceiling']
    assert not budget['single_monolithic_repair_AST_admitted']
    assert sampling['independent_ideal_repetitions']==17250000
    assert sampling['coarse_histogram_labels']==69


def test_stored_first_projective_word_including_reference_and_padding():
    word=m.word(0)
    assert len(word)==58392
    m.syn.certify_word(word)
    initial=np.zeros((66,2),complex);initial[0,0]=initial[1,1]=1/np.sqrt(2)
    current,padding=m.menus.successful_word(word,initial)
    ray=np.array([float(x) for x,_ in m.vector(0)]);ray/=np.linalg.norm(ray)
    row=ray@initial
    assert np.linalg.norm(current-(initial-np.outer(ray,row)))<1e-8
    assert np.linalg.norm(padding-np.outer(row,row.conj()))<1e-8
    admitted=m.syn.interpreter.source_word(word,9)
    stopped=m.syn.interpreter.execute_word(admitted,np.zeros((1,1)))
    assert stopped['status']=='RESOURCE_STOP' and stopped['pulses_executed']==0


@pytest.mark.parametrize('value',[-1,66,True,1.0])
def test_non_original_stage_address_rejected(value):
    with pytest.raises(ValueError):m.vector(value)


def test_zero_energy_margin_and_invalid_raw_word_rejected():
    with pytest.raises(ValueError):m.resolution(F(0))
    with pytest.raises(ValueError):m.preparation_branch((0,)*10)
    with pytest.raises(ValueError):m.preparation_branch((True,)*11)


def test_retained_receipt_source_and_exact_scientific_payload():
    """Fast publication check; the full emitted-word replay is the CLI owner."""
    value=packet.strict_json.load(packet.OUTPUT)
    assert value['source_bindings']==packet.source_bindings()
    assert value['spectral_certificate']==s.certify()
    assert value['complete_instrument_budget']==m.budget()
    assert packet.canonical(value['native_prepared_complete_energy_law'])==packet.canonical(m.prepared_reference_histogram())
    assert value['counted_energy_label_intervals']==[[str(x) for x in pair] for pair in s.energy_labels()]
    words=value['all_native_words']
    assert [row['stage'] for row in words]==list(range(66))
    assert all(row['resources']['native_pulses']==58392 for row in words)
    assert sum(row['resources']['native_pulses'] for row in words)==3853872
    assert all(row['generic_host_status']=='RESOURCE_STOP' and row['generic_host_pulses_executed']==0 for row in words)
    assert value['scope']['all66_spectral_labels_and_final_residual_retained']
    assert not value['scope']['full_expanded_generic_instrument_executed']
    assert not value['scope']['rounded_instrument_is_exactly_trace_preserving']


def test_existing_source_preparation_is_same_energy_reference_and_clock():
    value=reference.certificate()
    assert value['trace_P']==56 and value['faithful_reference_floor']=='1/84'
    assert value['owning_independent_numerical_verifier_passed']
    counts=np.array(s.source_counts(),float)
    eigenvalues,basis=np.linalg.eigh(counts/60)
    rho=counts/(60*56)
    modular=(basis*(-np.log(eigenvalues/56)))@basis.T
    generator=(basis*(-60*np.log(eigenvalues)))@basis.T
    assert np.linalg.norm(modular-generator/60-np.log(56)*np.eye(66))<1e-12
    assert np.linalg.norm(generator@rho-rho@generator)<1e-13
    assert value['preparation_ideal_timeout_probability']==str(F(23,128)**16)
    assert value['diagonal_tree_reference_separation']['half_trace_distance_lower']=='1/336'
