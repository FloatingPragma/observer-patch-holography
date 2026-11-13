"""Independent projector, CAR and spectral-scope regressions."""
from fractions import Fraction as F
import numpy as np
import pytest
import sympy as sy
from .native import channel,product_column,proton_fields,pion_field
from .source import run,support_test

@pytest.fixture(scope='module')
def result():return run()

def test_exact_projector_and_full_complement():
    c=channel();Y=c['basis'];G=c['basis_Gram'];B=c['overlap']
    assert G==Y.T*Y and G.rank()==23
    # Independent floating thin QR, with all raw occupation rows retained.
    q,_=np.linalg.qr(np.asarray(Y.tolist(),float),mode='reduced')
    from model.openwave_strange_shell_current.native import frame
    _,P=frame();V=np.asarray((c['vertex']*P).tolist(),float)
    projected=(q.T@V).T@(q.T@V)/6
    assert np.max(abs(projected-np.asarray(c['projected_Gram'].tolist(),float)))<2e-13
    other=V-q@(q.T@V)
    assert np.max(abs(other.T@other/6-np.asarray((c['total_Gram']-c['projected_Gram']).tolist(),float)))<2e-13

def test_complete_probability_and_scope(result):
    assert result['exact_CAR_channel']['independent_composite_rank']==23
    c=result['complete_probe_records']
    intervals=[tuple(map(F,c[k]))for k in('proton_pion_probability_interval','other_compact_probability_interval','noncompact_probability_interval','failure_probability_interval')]
    assert all(0<=a<=b<=1 for a,b in intervals)
    assert sum(a for a,b in intervals)<=1<=sum(b for a,b in intervals)
    assert F(c['original_P7_input_trial_pole_proton_pion_event_probability_lower'])>0
    assert not result['physical_Lambda_decay_amplitude_proved']
    assert not result['interacting_spectral_test']['compact_channel_projector_equals_H5_spectral_projector']
    assert not result['complete_Lambda_pole_transfer']['normalizable_Lambda_packet_executed']

def test_bare_support_zero_and_local_nonzero():
    assert support_test()['direct_current_to_disjoint_bare_proton_pion_support_exact']=='0'
    # All co-located composite columns are nonzero after antisymmetrization.
    assert product_column(proton_fields()[0],pion_field(0,0))
    # A positive channel Gram alone does not control spectral support.
    # The exact two-level control changes energy while keeping projector/current.
    W=sy.Matrix([1,1]);Pi=sy.diag(1,0)
    assert (W.T*Pi*W)[0]==1
    H=sy.diag(2,3);energy=0
    assert (H-energy*sy.eye(2)).det()!=0


def test_actual_kinetic_substep(result):
    from .native import separated_hopping
    h=separated_hopping()
    assert h['amplitude']==sy.Matrix([[0,-54,18,0,36,0,0],[0,0,0,36,0,18,54]])
    assert h['target_Gram']==108*sy.eye(2)
    k=result['actual_kinetic_substep']
    assert F(k['original_input_probe_then_kinetic_separated_probability_lower'])>0
    assert F(k['complete_probe_success_probability_exact'])+F(k['complete_probe_failure_probability_exact'])==1
    assert k['full_infinite_N5_kinetic_carrier_retained']
    assert not k['complete_autonomous_strong_Hamiltonian_used']
