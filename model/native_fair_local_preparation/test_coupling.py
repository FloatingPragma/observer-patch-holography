"""Complete occupation-chart and service controls."""
from fractions import Fraction as F
import numpy as np
import pytest
from . import coupling as c, preparation as p


@pytest.mark.parametrize('cell',[0,1,17])
def test_chart_is_exact_for_all_local_numbers_and_reference(cell):
    for q in range(13):
        local=(1<<q)-1;outside=(1<<(12*(cell+1)))
        full=c.join_mask(local,outside,cell)
        assert c.split_mask(full,cell)==(local,outside)
        rho={(full,2,full,4):F(3,7)}
        assert c.join_operator(c.split_operator(rho,cell,5),cell,5)==rho


def test_complete_native_intertwining_not_population_only():
    assert c.intertwining_certificate()['complete_raw_blocks']==3000


def test_wrong_chart_does_not_intertwine_a_fixed_address():
    rho={(1,0,1,0):F(1)}
    # A swap of labels without transporting the seam changes the operation.
    edge=(0,1);wrong={(4,0,4,0):F(1)}
    actual=p.unread(rho,(edge,));changed=p.unread(wrong,(edge,))
    def swap_chart(mask):
        return mask^5 if ((mask>>0)^(mask>>2))&1 else mask
    transported={(swap_chart(i),a,swap_chart(j),b):v for (i,a,j,b),v in actual.items()}
    keys=set(transported)|set(changed)
    assert sum(abs(transported.get(k,F())-changed.get(k,F())) for k in keys)/2==F(1,2)


def test_existing_macros_preserve_each_ideal_number_block():
    for name in c.words():
        target=c.ideal_matrix(name);numbers=np.array([x.bit_count() for x in range(len(target))])
        assert not np.any(target[numbers[:,None]!=numbers[None,:]])
        np.testing.assert_allclose(target.conj().T@target,np.eye(len(target)),atol=1e-15,rtol=0)


def test_source_service_bounds_and_failures_are_complete():
    value=c.service_certificate()
    assert [row['source_resources']['native_pulses'] for row in value['services']]==[22,2]
    for row in value['services']:
        assert F(row['error']['complete_halfdiamond_upper'])<F(1,10**8)
        assert row['all_first_failure_outputs_retained']
        assert not row['rounded_off_sector_work_output_discarded']


def test_uniform_number_reference_cannot_give_phase_response():
    for name in c.words():
        u=c.ideal_matrix(name);n=len(u)
        state=np.diag([F(x.bit_count()+1,100) for x in range(n)]).astype(complex)
        np.testing.assert_allclose(u@state@u.conj().T,state,atol=1e-15,rtol=0)


def test_input_chart_and_cost_validation_fail_closed():
    with pytest.raises(ValueError):c.join_mask(1,1)
    with pytest.raises(ValueError):c.split_mask(0,True)
    with pytest.raises(ValueError):c.complete_bound(True)
    with pytest.raises(ValueError):c.split_operator({(0,1,0,0):F(1)},reference_dimension=1)


def test_conditional_witness_covers_every_nonconstant_terminal():
    from . import witness
    assert witness.selected_edge(0) is None and witness.selected_edge(4095) is None
    for x in range(1,4095):
        edge=witness.selected_edge(x)
        assert edge in p.SEAMS and ((x>>edge[0])^(x>>edge[1]))&1


def test_full_compiled_phase_witness_and_inconclusive_budget():
    from . import witness
    value=witness.certificate()
    assert value['total_source_pulses']==305
    assert value['active_terminal_configurations']==4094
    assert value['inconclusive_constant_configurations']==2
    assert not value['constant_branches_postselected']
    assert not value['rounded_work_returned_as_exact_native_state']


def test_bad_conditional_word_parameters_fail_closed():
    from . import witness
    with pytest.raises(ValueError):witness.compact_word(1)
    with pytest.raises(ValueError):witness.declared_word((False,1),True)


def test_compiled_pointer_keeps_zero_complement_outcomes():
    from . import witness
    value=witness.pointer_certificate()
    assert value['complete_label_words']==1296
    assert value['zero_complement_label_words']==1280
    assert value['raw_bit_readout_calls']==4


def test_record_derived_recipe_and_exact_reference_valued_witness():
    from . import witness
    for q in range(1,12):
        x=(1<<q)-1;trace=p.local_trace(x,0,(0,)*11)
        rho={(x,3,x,8):F(2,7)}
        status,prepared=witness.ideal_prepared_branch(trace['slots'],rho)
        assert status=='ACTIVE' and len(prepared)==4
        for enabled in (False,True):
            status,rows=witness.ideal_detector_branches(trace['slots'],rho,enabled)
            assert len(rows)==25
            score=sum(score*sum(v for (i,a,j,b),v in block.items() if i==j) for _,score,block in rows)
            assert score==F(2,7*2**11)*enabled


def test_constant_recipe_preserves_coherence_and_zero_words():
    from . import witness
    trace=p.local_trace(0,0,(0,)*11)
    rho={(0,1,p.FULL,2):F(1,3)}
    status,prepared=witness.ideal_prepared_branch(trace['slots'],rho)
    assert status=='INCONCLUSIVE_CONSTANT' and prepared=={(0,1,p.FULL,2):F(1,3*2**11)}
    assert witness.recipe_from_tree((0,)*11,False)['status']=='ZERO_RAW_WORD'
