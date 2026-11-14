import numpy as np
from oph_openwave.construction import Construction
from oph_openwave.quantum_records import JointRecords, run


def test_integer_joint_generator_and_monitored_native_instrument():
    result = run()
    assert all(result['checks'].values())
    assert all(result['mutations'].values())


def test_kraus_completeness_gauss_and_native_probability():
    j = JointRecords(Construction(3,((0,1),(1,2),(2,0)), .11))
    before = j.gauss(j.states[0])
    total = np.zeros(len(j.states))
    for source,state in enumerate(j.states):
        for e in range(3):
            for q in (0,1):
                target,status,current = j.native_target(state,e,q)
                assert j.gauss(target)==before
                total[source] += 1/6
    assert np.allclose(total,1)
    # Explicit Kraus matrices on a coherent input independently reproduce commit.
    psi = np.arange(1,len(j.states)+1,dtype=complex)*(1+1j)
    psi /= np.linalg.norm(psi)
    rho = np.outer(psi,psi.conj())
    expected = np.zeros_like(rho)
    for source,state in enumerate(j.states):
        for e in range(3):
            for q in (0,1):
                target,_,_ = j.native_target(state,e,q)
                k = np.zeros_like(rho)
                k[j.index[target],source]=1/np.sqrt(6)
                expected += k@rho@k.conj().T
    actual,ledger = j.commit(rho)
    assert np.allclose(actual,expected)
    assert np.isclose(sum(r['probability'] for r in ledger),1)


def test_unregulated_edge_sector_matches_bare_generator_and_protocol_limits():
    c = Construction(2,((0,1),), .27)
    j = JointRecords(c)
    assert len(j.states)==3  # tree sector has no independent record cycles
    reorder = [j.index[next(s for s in j.states if s[0]==p)] for p in range(3)]
    assert np.allclose(j.hamiltonian[np.ix_(reorder,reorder)],c.hamiltonian)
    rho = j.initial_density()
    coherent,ledger = j.step(rho,.31,0)
    psi = c.evolve(np.array([1,0,0]),.31)
    assert np.allclose(coherent[np.ix_(reorder,reorder)],np.outer(psi,psi.conj()))
    committed,_ = j.step(rho,0,1)
    assert np.isclose(committed[0,0],.5)
    assert np.isclose(np.trace(committed),1)
    assert np.linalg.eigvalsh(committed).min()>-1e-12


def test_complete_ce_native_records_include_waits_and_contrast():
    j = JointRecords(Construction(2,((0,1),)))
    c, lc = j.commit(j.initial_density(),law='C')
    e, le = j.commit(j.initial_density(),law='E')
    assert np.allclose(c,e)
    assert {tuple(r['native_record']) for r in lc if r['source']==0} == {(0,1),(1,0)}
    assert {tuple(r['native_record']) for r in le if r['source']==0} == {(0,1),(1,1)}
    _, ledger = j.commit(j.initial_density())
    for law in ('C','E'):
        assert np.isclose(sum(r['probability'] for r in ledger if r['law']==law), .5)


def test_protocol_rejects_invalid_parameters_and_boundary_is_explicit():
    import pytest
    c = Construction(3,((0,1),(1,2),(2,0)))
    for kwargs in ({'cutoff':1.5},{'cutoff':0},{'initial_port':.5},{'max_states':0}):
        with pytest.raises(ValueError):
            JointRecords(c,**kwargs)
    j = JointRecords(c)
    for t,p in ((float('nan'),.5),(.1,float('nan')),(float('inf'),.5),(.1,-.1),(.1,1.1)):
        with pytest.raises(ValueError):
            j.step(j.initial_density(),t,p)
    for e,q in ((.5,0),(0,.5),(-1,0),(0,2)):
        with pytest.raises(ValueError):
            j.native_target(j.states[0],e,q)
    rejected = []
    for state in j.states:
        for e in range(3):
            for q in (0,1):
                target,status,current = j.native_target(state,e,q)
                if status=='boundary_rejected':
                    assert target == state
                    rejected.append(state)
    assert rejected


def test_joint_mass_square_uses_lifted_incidence_not_bare_graph():
    j = JointRecords(Construction(3,((0,1),(1,2),(2,0)),.37))
    gamma = np.diag([1 if s[0]<3 else -1 for s in j.states])
    d = j.hamiltonian-.37*gamma
    assert np.allclose(d@gamma+gamma@d,0)
    assert np.allclose(j.hamiltonian@j.hamiltonian,d@d+.37**2*np.eye(len(j.states)))


def test_rejected_native_event_has_zero_applied_current_and_no_native_record():
    j = JointRecords(Construction(3,((0,1),(1,2),(2,0))))
    rho = np.eye(len(j.states))/len(j.states)
    _,ledger = j.commit(rho)
    found = False
    for row in ledger:
        before = j.states[row['source']][1]
        after = j.states[row['target']][1]
        delta = np.array(after)-before
        expected = np.zeros(6,dtype=int)
        expected[2*row['seam']] = row['applied_current']
        expected[2*row['seam']+1] = -row['applied_current']
        assert np.array_equal(delta,expected)
        if row['status']=='boundary_rejected':
            found = True
            assert row['applied_current']==row['current']==0
            assert row['native_event_executed'] is False
            assert row['native_record'] is None
    assert found
