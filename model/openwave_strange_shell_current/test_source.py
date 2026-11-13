from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np
import pytest
import sympy as sy
from model.openwave_strange_shell_current import source as s
from model.openwave_strange_shell_current.native import contraction,frame,epsilon_column,vertex_words
from model.openwave_strange_composites.nonleptonic import word


def test_exact_current_gram_for_every_compact_lambda_vector():
    data=contraction();_,P=frame();G=P.T*P;M=data['seven_Gram_without_Cabibbo']
    assert all(x>=0 for x in (M-6*G).eigenvals())
    assert all(x>=0 for x in (7*G-M).eigenvals())
    assert data['matrix'].rows==1656 and data['matrix'].cols==56
    # Direct occupation output and reduced Gram must agree for a superposition
    # with cross terms, not only a source ground-vector fixture.
    v=sy.Matrix([1,-2,3,-4,5,-6,7]);out=data['matrix']*P*v
    assert (out.T*out)[0]/6==(v.T*M*v)[0]


def test_literal_onsite_singlet_has_exact_nonzero_current_and_pauli_signs():
    shapes=[((0,0),(0,1),(0,0)),((0,1),(0,0),(0,0))]
    initial={}
    for shape,coefficient in zip(shapes,(1,-1)):
        for mask,v in epsilon_column(shape).items():initial[mask]=initial.get(mask,0)+coefficient*v
    output={}
    for coefficient,ops in vertex_words():
        for mask,v in initial.items():
            final,sign=word(mask,ops)
            if sign:output[final]=output.get(final,0)+coefficient*sign*v
    assert sum(v*v for v in output.values())==6*sum(v*v for v in initial.values())
    # Independent elementary CAR ordering and Pauli blocking control.
    assert word(0,[(True,0),(True,1)])==(3,1)
    assert word(0,[(True,1),(True,0)])==(3,-1)
    assert word(1,[(True,0)])==(None,0)


def test_full_source_pole_transfer_and_all_probe_outcomes():
    r=s.run();response=r['full_pole_to_complete_compact_N5_response'];lo,hi=map(F,response['amplitude_norm_squared_interval'])
    assert F(28,100)<lo<=hi<F(34,100)
    p=r['complete_component_probe'];a,b=map(F,p['source_pole_compact_N5_probability_interval'])
    c,d=map(F,p['source_pole_noncompact_N5_probability_interval']);e,f=map(F,p['source_pole_probe_failure_probability_interval'])
    assert 0<a<=b and c==0 and d>=0 and 0<e<=f<1
    assert a+c+e<=1<=b+d+f
    assert r['physical_on_shell_decay_amplitude_proved']is False
    assert r['new_mass_or_current_calibration']is False
    shell=r['same_source_energy_shell'];assert F(shell['below_margin_GeV'])>0 and F(shell['above_margin_GeV'])>0


def test_trial_corruption_rejects_and_probe_completeness_on_independent_matrix():
    record=json.loads((Path(s.__file__).with_name('lambda_trial.json')).read_text());record['real_integer_coefficients'][0]='01'
    with pytest.raises(AssertionError):s.trial(record)
    J=np.array([[1.,2.,0.],[0.,3.,1.],[1.,0.,1.]])/36
    values,U=np.linalg.eigh(np.eye(3)-J.T@J);K=(U*np.sqrt(values))@U.T
    assert np.max(abs(J.T@J+K.T@K-np.eye(3)))<1e-14
    assert np.max(abs(np.kron(J.T@J+K.T@K,np.eye(2))-np.eye(6)))<1e-14
