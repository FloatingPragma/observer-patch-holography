from fractions import Fraction as F
import numpy as np
import pytest
from . import tree,execute
from ..observer.source_scalar import compiler
from ..certified_reference import source


@pytest.mark.parametrize('values',[(1,-2,0,3),(0,0,-3,4),(1,-3,2,0,-5,6,0,1)])
def test_signed_gray_tree_native_full_output(values):
    bits=(len(values)-1).bit_length();levels=tree.design(values,bits)
    word=[e for bit,raw,c in levels for e in tree.level_word(bit,c,bits)]
    initial=np.zeros(1<<bits);initial[0]=1
    run=execute.complete(word,bits,initial)
    expected=np.array(values)/np.linalg.norm(values)
    assert np.linalg.norm(run.success[:,0]-expected)<2e-13
    assert abs(np.linalg.norm(run.success)**2+run.failure_probability-1)<2e-13
    for bit,raw,c in levels:
        a,b=tree.level_error(values,bit,raw,c,bits)
        assert a+b<F(1,10**12)
    actual=execute.failure(word,run,len(word)-1,0)
    assert abs(np.vdot(actual,actual).real-run.failure_weights[-1,0])<1e-40


def test_gray_cycle_final_cnot_is_necessary_on_general_control_input():
    bit=0;coeff=np.array([.2,.3,-.1,.4]);word=list(tree.level_word(bit,coeff,3))
    exact=compiler.execute(word,3);broken=compiler.execute(word[:-3],3)
    assert np.linalg.norm(exact-broken)>1
    # On every control configuration the target block is the Walsh angle.
    angles=tree.walsh(coeff)
    for control,angle in enumerate(angles):
        rows=[2*control,2*control+1]
        expected=np.array([[np.cos(angle),np.sin(angle)],[-np.sin(angle),np.cos(angle)]])
        assert np.linalg.norm(exact[np.ix_(rows,rows)]-expected)<2e-14


def test_scalar_port_error_dominates_independent_full_generator_check():
    for kind,angle in [('Y',.123),('Y',-3.1),('H',np.pi/2),('J',-np.pi)]:
        actual=source.generator_error(kind,F(angle),source.literal_controls(kind,angle))
        assert actual<=tree.port_bound({(kind,angle):1})


def test_gray_level_integer_inputs_and_dimensions():
    with pytest.raises(ValueError):tree.design((1,2,3),2)
    with pytest.raises(ValueError):tree.design((1.,2.),1)
    with pytest.raises(ValueError):list(tree.level_word(0,[1],3))


def test_failure_reconstruction_rejects_a_changed_program():
    initial=np.array([1.,0.]);word=[compiler.pulse('Y',.1,[0])]
    run=execute.complete(word,1,initial)
    with pytest.raises(ValueError):execute.failure([compiler.pulse('Y',.2,[0])],run,0,0)
