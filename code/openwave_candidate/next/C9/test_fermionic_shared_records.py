import numpy as np
from fermionic_shared_records import FermionicSharedRecords,Construction,run


def test_complete_actual_fixedN_channel_recovery_and_work():
    r=run();assert all(r['checks'].values());assert r['dimension']==60
    assert not r['full_criterion_met']
    for output in r['outputs'].values():
        assert abs(output['trace']-1)<1e-12
        assert output['minimum_eigenvalue']>-1e-12
        assert abs(output['energy_after']-output['energy_before']-output['mean_instrument_work'])<1e-12


def test_actual_ground_positive_excitation_and_monitoring_work():
    m=FermionicSharedRecords(Construction(3,((0,1),(1,2),(2,0)),.23),7,3)
    energies,v=np.linalg.eigh(m.hamiltonian);ground=np.outer(v[:,0],v[:,0])
    for law in ('C','E','MIX'):
        out,ledger=m.commit(ground,law)
        assert ledger['mean_instrument_work']>1e-4
        assert np.trace((m.hamiltonian-energies[0]*np.eye(m.n))@out).real>=-1e-12


def test_external_reference_and_per_outcome_injective_maps():
    m=FermionicSharedRecords(Construction(3,((0,1),(1,2),(2,0)),.23),7,3)
    rng=np.random.default_rng(2026);psi=rng.normal(size=2*m.n)+1j*rng.normal(size=2*m.n);psi/=np.linalg.norm(psi)
    rho=np.outer(psi,psi.conj());output=np.zeros_like(rho)
    for meta,k in m.operators['MIX'][0]:
        nz=abs(k)>1e-12
        assert max(nz.sum(axis=0))<=1 and max(nz.sum(axis=1))<=1
        lifted=np.kron(k,np.eye(2));output+=lifted@rho@lifted.T
    assert abs(np.trace(output)-1)<1e-12
    assert np.linalg.eigvalsh(output).min()>-1e-12
    before=rho.reshape(m.n,2,m.n,2).trace(axis1=0,axis2=2)
    after=output.reshape(m.n,2,m.n,2).trace(axis1=0,axis2=2)
    assert np.max(abs(before-after))<1e-12
