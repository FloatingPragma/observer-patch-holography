from fractions import Fraction as F
import numpy as np
from scipy.linalg import expm
from . import algebra as a, circuits as c, execute


def test_exact_spin_and_complete_repair_algebra():
    assert a.certificate()['adjacent_bond_commutator_nonzero'] is True
    assert c.pulse_counts()['exchange_bond']==62


def test_native_encoding_and_exchange_preserve_quantum_reference():
    from .experiment import initial_and_target
    for name,template in [('encoded_spin_choi',c.encode(2)),('exchange_choi',c.exchange(F(1,3)))]:
        bits,initial,target=initial_and_target(name)
        run=execute.complete(c.emit(template),bits,initial)
        assert np.linalg.norm(run.success-target)<2e-12
        assert run.failure_probability<1e-24
        assert np.linalg.norm(run.success.conj().T@run.success-np.eye(4)/4)<2e-12


def test_fermionic_swap_retains_double_occupation_sign():
    run=execute.complete(c.emit(c.fswap(0)),2,np.eye(4)/2)
    target=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,-1]],complex)/2
    assert np.linalg.norm(run.success-target)<1e-12
    ordinary=target.copy();ordinary[3,3]=.5
    assert np.linalg.norm(run.success-ordinary)>.99


def test_native_pauli_probe_in_all_three_axes():
    p=a.inclusion(1)
    for axis in 'XYZ':
        initial=np.zeros((8,2),complex);initial[:4]=p/np.sqrt(2)
        run=execute.complete(c.emit(c.probe(axis,1,0,2)),3,initial)
        expected=np.vstack((p@(np.eye(2)+a.PAULI[axis])/2,p@(np.eye(2)-a.PAULI[axis])/2))/np.sqrt(2)
        assert np.linalg.norm(run.success-expected)<2e-12


def test_failure_reconstruction_binds_the_literal_program():
    word=list(c.emit(c.exchange(F(1,3))))
    run=execute.complete(word,4,a.inclusion(2)/2)
    branch=execute.failure(word,run,0,0)
    assert abs(np.vdot(branch,branch).real-run.failure_weights[0,0])<1e-40
    wrong=[dict(x) for x in word];wrong[0]['angle']+=.1
    import pytest
    with pytest.raises(ValueError):execute.failure(wrong,run,0,0)
