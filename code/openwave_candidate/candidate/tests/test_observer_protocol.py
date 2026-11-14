import numpy as np
import pytest
from oph_openwave.construction import Construction
from oph_openwave.quantum_records import JointRecords
from oph_openwave.observer_protocol import ObserverProtocol,RecordWindow,run
from oph_openwave.geometry import carrier


def test_complete_readout_matches_actual_instrument_ledger():
    j=JointRecords(Construction(3,((0,1),(1,2),(2,0)),.23))
    p=ObserverProtocol(j)
    psi=np.arange(1,len(j.states)+1,dtype=complex); psi/=np.linalg.norm(psi)
    rho=np.outer(psi,psi.conj())
    _,ledger=j.commit(rho)
    readout=p.readout(psi)
    for status,probability in readout['monitored_attempt_outcomes'].items():
        assert np.isclose(probability,sum(r['probability'] for r in ledger if r['status']==status))
    assert np.isclose(readout['port_probability']+readout['transit_probability'],1)
    assert np.isclose(sum(readout['applied_current_probabilities'].values()),1)


def test_sparse_trajectory_ensemble_agrees_with_exact_density():
    j=JointRecords(Construction(2,((0,1),),.23)); p=ObserverProtocol(j)
    psi=np.zeros(len(j.states),complex); psi[0]=1
    exact,_=j.step(j.initial_density(),.37,.29)
    rng=np.random.default_rng(9317); counts=np.zeros(len(psi))
    for _ in range(500):
        evolved,record=p.trajectory_step(psi,.37,.29,rng)
        counts+=np.abs(evolved)**2
    assert np.max(np.abs(counts/500-np.diag(exact).real))<.065


def test_carrier_window_sparse_gauss_and_probability():
    j=RecordWindow(carrier(.23),radius=4,max_states=4096)
    p=ObserverProtocol(j)
    assert (p.h-p.h.T).nnz==0
    assert all(j.gauss(s)==j.gauss(j.states[0]) for s in j.states)
    psi=np.zeros(len(j.states),complex); psi[0]=1
    evolved=p.evolve(psi,.37)
    assert np.isclose(np.vdot(evolved,evolved),1)
    assert np.isclose(sum(p.readout(evolved)['position_probabilities']),1)
    assert len(j.states)<4096


def test_monitoring_controls_and_invalid_inputs():
    result=run()
    assert all(result['checks'].values()) and all(result['mutations'].values())
    j=JointRecords(Construction(2,((0,1),))); p=ObserverProtocol(j)
    psi=np.zeros(len(j.states),complex); psi[0]=1
    for value in (float('nan'),float('inf'),-1,2):
        with pytest.raises(ValueError):p.trajectory_step(psi,.1,value,np.random.default_rng(0))
    with pytest.raises(ValueError):p.evolve(psi,float('nan'))
    with pytest.raises(ValueError):p.readout(psi*2)


def test_window_bound_against_complete_triangle_sector():
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    full=JointRecords(c)
    for radius in (0,1,2,3):
        window=RecordWindow(c,radius=radius)
        psi=np.zeros(len(window.states),complex); psi[0]=1
        projected=ObserverProtocol(window).evolve(psi,.37)
        embedded=np.zeros(len(full.states),complex)
        for k,state in enumerate(window.states): embedded[full.index[state]]=projected[k]
        initial=np.zeros(len(full.states),complex); initial[0]=1
        exact=ObserverProtocol(full).evolve(initial,.37)
        assert np.linalg.norm(exact-embedded)<=window.initial_state_error_bound(.37)+1e-12


def test_bare_compact_circulation_does_not_transfer_to_integer_records():
    from oph_openwave.observer_protocol import circulation_control
    result=circulation_control()
    assert result['bare_incidence_shape']==[3,3]
    assert result['bare_exact_rank']==2
    assert result['joint_incidence_shape']==[7,6]
    assert result['joint_exact_rank']==6
    assert result['bare_port_leakage_probability']<1e-20
    assert result['joint_port_leakage_probability']>.01
    assert result['joint_prepared_support_survival_probability']<.99
    # Independent exact graph check: the lifted bipartite state graph is a
    # connected path, with seven port vertices and six transit vertices.
    j=JointRecords(Construction(3,((0,1),(1,2),(2,0)),.23))
    h=ObserverProtocol(j).h.copy(); h.setdiag(0); h.eliminate_zeros()
    degrees=np.diff(h.indptr)
    assert sorted(degrees.tolist())==[1,1]+[2]*11
