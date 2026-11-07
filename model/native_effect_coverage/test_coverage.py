from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import json
import numpy as np
import pytest
from . import detector as d, coverage as c, packet


@pytest.fixture(scope='module')
def coherent():
    return c.coherent_certificate()


def test_literal_program_parameters_and_all_native_columns():
    cert=d.certificate()
    assert cert['source_kernel_columns_checked']==1452
    assert cert['native_record_storage_bits_with_terminal_status']==480
    assert cert['nonzero_rank_one_words']==135168
    assert cert['complete_raw_words']==25**11
    for edge in d.tree():
        rows=d.kernel(edge)
        assert len(rows)==25
        assert sum(not row.pairs for row in rows)==19


@pytest.mark.parametrize('pattern',[0,1,2,3,5,7,13,2047])
def test_all_source_basis_states_under_independent_coin_patterns(pattern):
    coins=tuple((pattern>>i)&1 for i in range(11))
    for initial in range(66):
        value=d.run(initial,coins)
        weight,mapping,records=d.branch(value['slots'])
        assert mapping=={initial:value['final']}
        assert weight==value['weight']
        assert records==value['records']


def test_current_correction_is_necessary():
    found=False
    for initial in range(66):
        run=d.run(initial,(1,)*11)
        naive=tuple(record.label for record in run['records'])
        truth=tuple(d.local.basis()[initial][a]-d.local.basis()[initial][b] for a,b in d.tree())
        found |= naive!=truth
    assert found


def test_zero_address_and_coin_mutations_fail(monkeypatch):
    original=d.kernel
    edge=d.tree()[0];rows=original(edge)
    dropped=next(i for i,row in enumerate(rows) if not row.pairs)
    monkeypatch.setattr(d,'kernel',lambda e:tuple(row for i,row in enumerate(rows) if i!=dropped)
                        if e==edge else original(e))
    with pytest.raises(ArithmeticError,match='raw address'):d.certificate()
    bad=list(rows);live=next(i for i,row in enumerate(rows) if row.pairs)
    bad[live]=replace(bad[live],weight=bad[live].weight+F(1,100))
    monkeypatch.setattr(d,'kernel',lambda e:tuple(bad) if e==edge else original(e))
    with pytest.raises(ArithmeticError,match='endpoint'):d.certificate()


def test_disconnected_schedule_and_corrupt_record_rejected():
    run=d.run(0,(0,)*11)
    bad=list(run['records']);bad[0]=replace(bad[0],current=99)
    with pytest.raises(ValueError):d.decode(tuple(bad))
    with pytest.raises(ValueError):d.decode(run['records'],((0,1),)*11)
    # Two unmeasured vertices can exchange occupancy without any contrast.
    measured=((0,1),(0,5),(0,7),(0,11),(1,8),(1,9),(4,5),(6,7),(3,8))
    x=tuple(int(i in (0,2)) for i in range(12));y=tuple(int(i in (0,10)) for i in range(12))
    assert x!=y
    assert tuple(x[a]-x[b] for a,b in measured)==tuple(y[a]-y[b] for a,b in measured)


def test_complete_two_step_superoperator_against_expanded_interpreter():
    edges=d.tree()[:2]
    rho={(0,0,0,0):F(2,5),(2,1,2,1):F(3,5),
         (0,0,2,1):F(1,7),(2,1,0,0):F(1,7)}
    actual=d.interpreter.execute(d.program(edges),rho)
    assert actual.status=='COMPLETE'
    assert len(actual.branches)==25**2
    for slots,outcome in zip(product(range(25),repeat=2),actual.branches,strict=True):
        assert outcome.density==d.apply_path(slots,rho,edges)


def test_complete_reference_word_and_zero_slots():
    run=d.run(12,tuple(i%2 for i in range(11)))
    rho={(i,a,j,b):F(1,1+i+a+j+b) for i,j in ((12,12),(12,4),(4,12),(4,4))
         for a,b in product(range(2),repeat=2)}
    output=d.continue_raw(run['slots'],rho)
    assert output['density']==d.apply_path(run['slots'],rho)
    assert len(output['archive'])==10 and len(output['memory'])==1
    for position in range(11):
        zero=next(k for k,row in enumerate(d.kernel(d.tree()[position])) if not row.pairs)
        slots=list(run['slots']);slots[position]=zero
        assert d.branch(tuple(slots))[1]=={}
        assert d.continue_raw(tuple(slots),rho)['density']=={}


def test_stock_evaluator_stop_is_preserved():
    value=packet.execution_certificate()
    assert value['expanded_host_status']=='RESOURCE_STOP'
    assert value['expanded_attempts_executed']==3
    assert value['expanded_retained_blocks']==15625
    assert not value['host_stop_relabelled_as_success']


def test_reference_safe_reset_and_restoration_aggregates():
    rho={(i,a,j,b):F(1,1+i+a+j+b) for i,j in ((0,0),(0,65),(65,0),(65,65))
         for a,b in product(range(2),repeat=2)}
    outputs=d.aggregated_output(rho,'prepare',target=17)
    total={}
    for block in outputs.values():
        for key,v in block.items():total[key]=total.get(key,F())+v
    assert total=={(17,a,17,b):sum(rho.get((i,a,i,b),F()) for i in range(66))
                  for a,b in product(range(2),repeat=2)}
    restored=d.aggregated_output(rho,'restore')
    for (i,j),block in restored.items():
        assert all(x==y==i for x,a,y,b in block)
    raw=d.aggregated_output(rho,'raw')
    assert any(raw[key] and raw[key]!=restored[key] for key in raw)


def test_actual_xor_words_and_padding_preserved():
    for initial,target in ((0,65),(20,3),(17,17)):
        word=c.reset_word(initial,target)
        vector=np.zeros(512,complex);vector[c.t.codes()[initial]]=1
        actual=c.compiler.execute(word,9,vector) if word else vector
        ideal=np.zeros(512,complex);ideal[c.t.codes()[target]]=1
        assert np.linalg.norm(actual-ideal)<1e-10
    assert c.reset_word(0,0)==[]
    with pytest.raises(ValueError):c.reset_word(True,4)


def test_inverse_y_preparation_not_forward_y():
    a,b=0,65
    initial=np.zeros(512,complex);initial[c.t.codes()[a]]=1
    inverse=c.compiler.execute(c.preparation_word(a,b,'y_analyzer'),9,initial)[list(c.t.codes())]
    forward=c.compiler.execute(c.t.word(a,b,'y_analyzer'),9,initial)[list(c.t.codes())]
    expected=c.ideal_preparation(a,b,'y_analyzer')
    assert np.linalg.norm(np.outer(inverse,inverse.conj())-expected)<1e-10
    assert np.linalg.norm(np.outer(forward,forward.conj())-expected)>.5
    assert c.preparation_certificate()['spanning_preparations']==4356


def test_complete_compiler_error_and_actual_span(coherent):
    assert coherent['raw_M6_primitive_kinds_checked']==16
    assert coherent['inherited_pair_routes_rechecked']==2145
    assert F(coherent['uniform_words']['analyzer']['complete_halfdiamond_upper'])<F(621,10**12)
    assert F(coherent['uniform_words']['coherent_preparation_after_tree']['complete_halfdiamond_upper'])<F(801,10**12)
    stable=c.stable_span_certificate(coherent)
    assert F(stable['effect_analysis_smallest_singular_value_lower'])>F(999999,1000000)
    assert F(stable['actual_preparation_separating_lower_factor'])>F(999999,1000000)
    corrupted=deepcopy(coherent)
    corrupted['uniform_words']['analyzer']['complete_halfdiamond_upper']='1/2'
    with pytest.raises(ArithmeticError):c.stable_span_certificate(corrupted)


def test_full_execution_and_structural_padding(coherent):
    report=packet.compiled_execution_certificate(coherent)
    assert report['structural_ideal_padding_abort_observed']
    assert report['analyzer_and_inverse_preparation_pairs_checked']==6


def test_arithmetic_guard_and_native_primitive_mutation(monkeypatch):
    assert c.block_arithmetic_bound()<F(1,10**11)
    with pytest.raises(ValueError):c.block_arithmetic_bound(36)
    with pytest.raises(ValueError):c.service_bound(41,131)
    original=c.clock._processor_bounds
    monkeypatch.setattr(c.clock,'_processor_bounds',lambda array,dimension:(F(1,100),F()))
    with pytest.raises(ArithmeticError,match='primitive'):c.coherent_certificate()
    monkeypatch.setattr(c.clock,'_processor_bounds',original)


def test_all_effect_coordinates_and_separating_preparations():
    rng=np.random.default_rng(191)
    matrix=rng.normal(size=(66,66))+1j*rng.normal(size=(66,66));matrix+=matrix.conj().T
    diagonal=np.diag(matrix).real
    for a,b in ((0,1),(0,65),(15,37),(64,65)):
        for kind in ('mixer','y_analyzer'):
            actual=np.trace(c.ideal_preparation(a,b,kind)@matrix).real
            wanted=(diagonal[a]+diagonal[b])/2+(matrix[a,b].real if kind=='mixer' else -matrix[a,b].imag)
            assert abs(actual-wanted)<1e-12
    assert c.coverage_certificate()['real_effect_span_dimension']==4356
    assert 66+8*66*65==34386


def test_rank_one_positive_menu_obstruction():
    v=np.array([.6,.8],complex);target=np.outer(v,v.conj())
    assert np.linalg.matrix_rank(target)==1
    atoms=[np.array([1,0]),np.array([0,1])]
    atoms += [np.array([1,z])/np.sqrt(2) for z in (1,-1,1j,-1j)]
    assert all(np.linalg.norm((np.eye(2)-target)@a)>.1 for a in atoms)
    result=c.coverage_certificate()
    assert not result['universal_all_effect_realizability_or_coexistent_refinement']


def test_native_diagonal_refinement_counts_and_coexistence_guard():
    f=tuple(F(i%5,8) for i in range(66));g=tuple(F((i+1)%3,8) for i in range(66))
    output=c.dyadic_diagonal_refinement(f,g,3)
    assert all(a+b+r==8 for a,b,r in zip(output['accept_f_counts'],output['accept_g_counts'],output['reject_counts'],strict=True))
    with pytest.raises(ValueError):c.dyadic_diagonal_refinement((F(1),)*66,(F(1),)*66,3)
    with pytest.raises(ValueError):c.dyadic_diagonal_refinement((F(1,3),)*66,(F(),)*66,3)
    with pytest.raises(ValueError):c.dyadic_diagonal_refinement(f,g,54)


@pytest.mark.parametrize('bad',[(True,)*11,(25,)*11,(-1,)*11,(0,)*10])
def test_raw_word_domain_guard(bad):
    with pytest.raises(ValueError):d.branch(bad)


def test_receipt_pins_no_documentation_and_duplicate_json_rejection(tmp_path):
    bindings=packet.source_bindings()
    assert not any(p.startswith('docs/') for p in bindings)
    assert 'model/native_program/interpreter.py' in bindings
    assert all(len(h)==64 for h in bindings.values())
    path=tmp_path/'duplicate.json';path.write_text('{"schema":1,"schema":2}')
    with pytest.raises(ValueError):packet.strict_json.load(path)


def test_retained_receipt_and_mutation(tmp_path):
    packet.verify()
    value=packet.strict_json.load(packet.OUTPUT)
    value['native_detector']['complete_basis_projectors']=65
    path=tmp_path/'bad.json';path.write_text(json.dumps(value))
    with pytest.raises(ValueError,match='receipt differs'):packet.verify(path)
