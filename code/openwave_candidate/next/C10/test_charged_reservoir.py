import numpy as np
from charged_reservoir import BoundarySource,FermionicSharedRecords,Construction,run


def test_complete_charged_instrument_and_energy_outcomes():
    r=run();assert all(r['checks'].values())
    p=r['final_source_number_probabilities'];assert p['2']>0 and p['3']>0 and p['4']>0
    for row in r['records']:
        assert abs(sum(x['probability'] for x in row['complete_number_energy_eigenspace_records'])-row['probability'])<1e-12
        assert abs(sum(x['energy']*x['probability'] for x in row['complete_number_energy_eigenspace_records'])-row['unnormalized_source_energy'])<1e-12
    assert not r['full_criterion_met']


def test_same_N3_source_and_every_native_K_exactly():
    full=BoundarySource(3);old=FermionicSharedRecords(Construction(3,((0,1),(1,2),(2,0)),.23),7,3)
    order=[full.index[state] for state in old.states]
    assert np.array_equal(full.hamiltonian[np.ix_(order,order)],old.hamiltonian)
    for law in old.operators:
        for ops,other in zip(full.operators[law],old.operators[law]):
            keyed={tuple(m.values()):k for m,k in ops}
            for meta,k in other:assert np.array_equal(keyed[tuple(meta.values())][np.ix_(order,order)],k)


def test_all_Gauss_native_and_charged_maps_with_idle_reference():
    m=BoundarySource(3);g=(1,1,1,0,0,0)
    assert all(m.gauss(s)==g for s in m.states)
    for law in m.operators:
        for ops in m.operators[law]:assert np.allclose(sum(k.T@k for _,k in ops),np.eye(m.n))
    # Kraus completeness gives identity dual action on every idle reference;
    # check two-dimensional entangled input with no density positivity shortcut.
    rng=np.random.default_rng(123);psi=rng.normal(size=(m.n,2))+1j*rng.normal(size=(m.n,2));psi/=np.linalg.norm(psi)
    before=psi.conj().T@psi;after=np.zeros((2,2),complex)
    for meta,k in m.reservoir_instrument():
        branch=k@psi;after+=branch.conj().T@branch
    assert np.max(abs(before-after))<1e-12
