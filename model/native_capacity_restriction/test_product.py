from fractions import Fraction as F

import pytest

from . import product as p


@pytest.mark.parametrize('shell',('H66','H78'))
@pytest.mark.parametrize('cells',(2,3,64))
def test_adaptive_complete_history_on_arbitrary_correlated_operator_parts(shell,cells):
    value = p._program(shell)
    for density in p._state(shell,cells):
        small,large = p.compare(value,cells,density)
        assert small.status == 'COMPLETE' and small.attempted_steps == 4
        assert len(small.branches) == 625
        assert any(not branch.density for branch in small.branches)
        assert small.program_digest != large.program_digest
        assert all(len(branch.archive) == 3 for branch in small.branches)


def test_more_than_one_retained_active_cell_and_exact_nested_restriction():
    value = p.lift_program(p._program(),2)
    value['commands'][0]['cell'] = 1
    value['declared_resources'] = p.source._native_cost(value)
    value = p.source.admit(value)
    for density in p._state('H66',4,retained=2):
        small,_ = p.compare(value,4,density)
        assert small.status == 'COMPLETE'
        once = p.restrict_operator(density,'H66',4,2)
        nested = p.restrict_operator(p.restrict_operator(density,'H66',4,3),'H66',3,2)
        assert once == nested


@pytest.mark.parametrize('selector',(
    {'kind':'uniform'},
    {'kind':'observer','positive':[1,7],'otherwise':[2,3]}))
def test_uniform_and_exploratory_source_selectors_keep_their_full_alphabets(selector):
    value = p.source.program(commands=[{
        'op':'repair','selector':selector,'law':'mixture','repeat':1}])
    for density in p._state('H66',2):
        small,_ = p.compare(value,2,density)
        assert len(small.branches) == (721 if selector['kind'] == 'uniform' else 745)


@pytest.mark.parametrize('case',('FIFO','attempt','branch'))
def test_stops_are_retained_and_synchronized(case):
    value = p._program(overflow='halt' if case == 'FIFO' else 'export_oldest')
    limits = {'max_steps':1} if case == 'attempt' else {'max_branches':1} if case == 'branch' else {}
    small,_ = p.compare(value,2,p._state('H66',2)[0],**limits)
    assert small.status == ('FIFO_TIMEOUT' if case == 'FIFO' else 'RESOURCE_STOP')
    assert all(branch.halted for branch in small.branches)


def test_off_diagonal_appended_factor_has_zero_trace_without_losing_raw_slots():
    value = p._program()
    density = {(0,0,66,1):F(1)}
    assert p.restrict_operator(density,'H66',2,1) == {}
    small,large = p.compare(value,2,density)
    assert all(not branch.density for branch in small.branches)
    assert any(branch.density for branch in large.branches)
    assert len(small.branches) == len(large.branches) == 625


def test_one_factor_lift_changes_only_typed_carrier_and_identifier():
    value = p._program(); lifted = p.lift_program(value,1)
    density = {(0,0,0,0):F(1)}
    small,large = p.compare(value,1,density)
    assert small.branches == large.branches and small.program_digest != large.program_digest
    assert lifted['parameters'] == value['parameters']
    assert lifted['memory'] == value['memory']


def test_resource_and_specification_reads_are_not_identified():
    value = p._program(); lifted = p.lift_program(value,64)
    assert value['carrier'] != lifted['carrier']
    assert value['declared_resources']['record_bits_per_slot_upper'] < lifted['declared_resources']['record_bits_per_slot_upper']
    assert value['declared_resources']['program_description_bytes'] != lifted['declared_resources']['program_description_bytes']


def test_early_carrier_guard_is_outside_the_execution_equivalence():
    density = p._state('H66',65)[0]
    with pytest.raises(ArithmeticError,match='does not intertwine'):
        p.compare(p._program(),65,density)


def test_direct_fair_federation_comparison_fails_on_actual_coin_records():
    result = p.fair_direct_countercontrol()
    assert result['complete_output_TV_lower_bound'] == '17/66'


def test_lift_rejects_compiled_coherent_service_words_and_federation_members():
    word = p.source.source_word([p.source.compiler.pulse('Z',0.,[0])],1)
    with pytest.raises(ValueError):
        p.lift_program(word,2)
    federation = p.source.program(federation={'level':0,'total':2})
    with pytest.raises(ValueError,match='only for local'):
        p.lift_program(federation,2)


@pytest.mark.parametrize('cells',(0,-1,True,1.5))
def test_malformed_product_sizes_are_rejected(cells):
    with pytest.raises(ValueError):
        p.lift_program(p._program(),cells)


def test_addressed_factor_cannot_be_removed_or_relabelled_as_a_different_target():
    value = p.lift_program(p._program(),2)
    with pytest.raises(ValueError):
        p.lift_program(value,1)
    history = p.source.execute(value,p._state('H66',2)[0])
    with pytest.raises(ValueError,match='target identifier'):
        p.restrict_history(history,'H66',2,1,target_program=value)
    with pytest.raises(ValueError):
        p.restrict_operator({},'H66',2,3)


def test_mutating_an_appended_cell_breaks_the_intertwining_hypothesis():
    value = p._program(); lifted = p.lift_program(value,2)
    lifted['commands'][0]['cell'] = 1
    lifted['declared_resources'] = p.source._native_cost(lifted)
    density = p._state('H66',2)[0]
    large = p.source.execute(lifted,density)
    with pytest.raises(ValueError,match='exact lifted target'):
        p.restrict_history(large,'H66',2,1,target_program=value)


def test_known_identifier_relabel_rejects_a_different_same_carrier_program():
    value = p._program(); density = p._state('H66',2)[0]
    large = p.source.execute(p.lift_program(value,2),density)
    other = p.source.admit(value)
    other['commands'][0]['law'] = 'E'
    other['declared_resources'] = p.source._native_cost(other)
    with pytest.raises(ValueError,match='exact lifted target'):
        p.restrict_history(large,'H66',2,1,target_program=other)
