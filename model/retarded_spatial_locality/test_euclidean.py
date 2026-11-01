from fractions import Fraction as F

import numpy as np
import pytest

from model.field_symmetry import symmetry
from . import euclidean as e


def test_exact_reflection_matrix_and_all_principal_minors():
    value=e.certificate()
    assert value['all_principal_minors_checked']==31
    assert value['signed_polynomial_norm_squared']=='209/256'
    assert value['exact_reflection_Gram']==[
        ['1/4','1/64','0','0','0'],['1/64','1/1024','0','0','0'],
        ['0','0','1/16','0','0'],['0','0','0','1/64','-1/64'],
        ['0','0','0','-1/64','1/64']]


def test_rational_CAR_fixture_agrees_with_independent_source_operator():
    for mode in (0,1):
        source=symmetry._annihilator(2,mode).matrix().toarray()
        assert np.array_equal(source,np.asarray(e._annihilator(mode),dtype=int))


def test_equal_time_odd_order_has_a_sign_and_cannot_be_sorted():
    words=e.witness_words();energies=(0,1,2,3)
    first=e.half_word(*words[3],energies)
    second=e.half_word(*words[4],energies)
    assert first==tuple(-x for x in second)
    assert e.reflected_entry(words[3],words[4],energies)==-F(1,64)
    # Sorting the two odd operators would flip a genuinely negative cross term.
    assert e.reflected_entry(words[3],words[3],energies)==F(1,64)


def test_reflection_requires_reversing_the_adjoint_word():
    operators,times=e.witness_words()[3];energies=(0,1,2,3)
    wrong=tuple(e.adjoint(a) for a in operators)+operators
    wrong_times=(-1,-1,1,1)
    assert e.ordered_expectation(wrong,wrong_times,energies)==-F(1,64)
    assert e.reflected_entry((operators,times),(operators,times),energies)==F(1,64)


def test_operator_gaps_are_not_interchangeable():
    a=e.adjoint(e._annihilator(0));b=e.adjoint(e._annihilator(1))
    assert e.half_word((a,b),(1,2),(0,1,2,3))[3]==F(1,32)
    assert e.half_word((b,a),(1,2),(0,1,2,3))[3]==-F(1,16)


def test_negative_generator_and_noncontractive_data_are_rejected():
    for bad in ((0,-1),(1,2),(0,F(1,2)),(0,True)):
        with pytest.raises(ValueError):e.damping(bad,1)
    for bad in (0,-1,F(3,2),0.5):
        with pytest.raises(ValueError):e.damping((0,1),1,bad)
    with pytest.raises(ValueError):e.damping((0,1),-1)
    a=e.adjoint(e._annihilator(0))
    with pytest.raises(ValueError):e.half_word((a,a),(2,1),(0,1,2,3))
    with pytest.raises(ValueError):e.half_word((a,),(-1,),(0,1,2,3))


def test_Poisson_tail_bound_handles_zero_gaps_and_multiple_integrals():
    assert e.poisson_tail_bound((0,), (1,))==0
    assert e.poisson_tail_bound((F(1,2),F(3,2)),(100,100))==F(1,75)
    assert e.poisson_tail_bound((100,),(1,))==1
    assert e.poisson_tail_bound((1,),(20,))==e.poisson_tail_bound((1,),(10,))/2
    for gaps,cutoffs in (((-1,),(1,)),((1,),(0,)),((1,),()),((0.5,),(1,))):
        with pytest.raises(ValueError):e.poisson_tail_bound(gaps,cutoffs)


def test_fixed_spacing_scope_does_not_add_rotational_or_continuum_axioms():
    value=e.certificate()
    assert value['all_ordered_local_Euclidean_npoint_limits_exist']
    for name in ('bosonic_permutation_invariance_inferred','Euclidean_O4_covariance_inferred',
                 'continuum_OS_axioms_completed','limiting_hierarchy_values_computed'):
        assert value[name] is False
