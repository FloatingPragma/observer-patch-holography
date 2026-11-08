"""Independent controls for exact all-sector preparation."""
from fractions import Fraction as F
from math import comb
import numpy as np
import pytest
from . import preparation as p


@pytest.mark.parametrize('q',range(13))
def test_integer_minorization_and_native_column(q):
    matrix=p.block_numerators(q);states=p.masks(q);index={x:i for i,x in enumerate(states)}
    assert matrix.dtype==np.dtype('int64') and not matrix.flags.writeable
    assert 2*sum(map(int,matrix.min(axis=1)))>=2**60
    # Independent population recurrence follows the actual endpoint coins.
    values={states[len(states)//2]:1}
    for edge in p.SEAMS*2:
        following={}
        for x,count in values.items():
            for coin in (0,1):
                y,_,_=p.readout.native.image(x,edge,coin,12)
                following[y]=following.get(y,0)+count
        values=following
    assert [values.get(x,0) for x in states]==list(matrix[:,len(states)//2])
    assert int(matrix.sum(axis=0)[0])==2**60


def test_reference_channel_preserves_extreme_entanglement_and_number_weights():
    rho={(0,0,0,0):F(1,4),(p.FULL,1,p.FULL,1):F(1,4),
         (0,0,p.FULL,1):F(1,4),(p.FULL,1,0,0):F(1,4),(3,2,3,2):F(1,2)}
    target=p.reference_channel(rho)
    assert target[0,0,p.FULL,1]==F(1,4)
    assert all(target[x,2,x,2]==F(1,132) for x in p.masks(2))
    assert p.reference_channel(target)==target
    assert sum(v for (i,a,j,b),v in target.items() if i==j and a==b)==1


def test_tree_kills_nonconstant_coherences_with_reference():
    cross={(1,3,2,7):F(2,5)}
    assert p.unread(cross,p.TREE)=={}
    corner={(0,3,p.FULL,7):F(2,5)}
    assert p.unread(corner,p.edges(1))==corner


def test_retained_tree_reference_keeps_its_actual_record_weight():
    trace=p.local_trace(3,0,(0,)*11)
    rho={(3,4,3,8):F(3,7)}
    target=p.retained_tree_reference(trace['slots'],rho)
    assert target=={(z,4,z,8):F(3,7*66*2**11) for z in p.masks(2)}


def test_one_block_actual_tv_bound_and_column_stationarity():
    for q in (1,2,6):
        value=p.block_numerators(q);d=comb(12,q)
        for column in (0,d//2,d-1):
            # All numerator arithmetic, no eigensolver or floating tolerance.
            tv_numerator=sum(abs(int(x)*d-2**60) for x in value[:,column])
            assert tv_numerator<=d*2**60
        assert all(int(x)==2**60 for x in value.sum(axis=1))


def test_full_archive_factorization_control():
    assert p.archive_separation(6)==F(923,924)
    assert p.archive_separation(2)==F(65,66)
    assert p.archive_separation(2)>F(1,2**32)


@pytest.mark.parametrize('bad',[True,-1,69,1.0,'1'])
def test_budget_fails_closed(bad):
    with pytest.raises(ValueError):p.budget(bad)


def test_source_budget_and_no_new_operation():
    value=p.program(blocks=32)
    assert len(value['commands'])==1931
    assert value['declared_resources']['attempts_upper']==1931
    assert all(c['op']=='repair' and c['law']=='E' and c['repeat']==1 for c in value['commands'])
    assert p.budget(68)['AST_nodes']==4091
    with pytest.raises(ValueError):p.program(cell=True)
    with pytest.raises(ValueError):p.program(cell=20)


def test_canonical_and_bernoulli_restrictions_are_fixed_without_selecting_weights():
    rho={}
    for q in range(13):
        for x in p.masks(q):rho[x,0,x,0]=F(2**q,3**12)
    assert sum(rho.values())==1
    assert p.reference_channel(rho)==rho
    point={(0,0,0,0):F(1)}
    assert p.reference_channel(point)==point
    assert p.reference_channel(point)!=rho


def test_phase_record_obstruction_uses_all_raw_labels():
    assert p.phase_blindness_control()['actual_raw_effect_controls']==1500


def test_native_source_raw_word_with_spectator():
    value=p.source_query(total=3,mask=3,outside=1<<20,blocks=1)
    assert value['attempts']==71 and value['raw_weight']==str(F(1,2**71))
    assert value['exported_records']==70


def test_integer_certificate_arrays_cannot_be_mutated():
    value=p.block_numerators(2)
    with pytest.raises(ValueError):value.flags.writeable=True
    value.shape=(66*66,)
    assert p.block_numerators(2).shape==(66,66)
