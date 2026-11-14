import numpy as np
import pytest
from oph_openwave.modular_native import LocalCoarseModular,Construction,run


def model():return LocalCoarseModular(Construction(3,((0,1),(1,2),(2,0)),.23),6)


def test_complete_local_kraus_and_gauss_blocks():
    m=model();gauss=m.gauss_residue(m.states[0])
    for seam,operators in enumerate(m.operators):
        assert np.allclose(sum(k.T@k for _,k in operators),np.eye(m.n))
        for record,k in operators:
            assert set(record)=={'seam','law','coin','record_value'}
            for target,source in zip(*np.nonzero(k)):
                assert m.gauss_residue(m.states[target])==m.gauss_residue(m.states[source])==gauss
    assert m.n==36


def test_complete_branches_and_transit_coherence_no_global_measurement():
    m=model();transit=[i for i,s in enumerate(m.states) if s[0]>=3]
    psi=np.zeros(m.n,complex);psi[transit[0]]=psi[transit[-1]]=1/np.sqrt(2)
    rho=np.outer(psi,psi.conj());actual,ledger=m.commit(rho,seam=0)
    assert np.allclose(actual,rho) # untouched coherent transit subspace
    assert np.isclose(sum(r['probability'] for r in ledger),1)
    measured,_=m.commit(rho,seam=0,fully_monitored_control=True)
    assert not np.allclose(actual,measured)
    assert np.linalg.eigvalsh(actual).min()>-1e-12


def test_generator_trace_complement_and_zero_rate():
    m=model();rng=np.random.default_rng(9317)
    psi=rng.normal(size=m.n)+1j*rng.normal(size=m.n);psi/=np.linalg.norm(psi)
    rho=np.outer(psi,psi.conj());rate=.29/(3*.37)
    derivative=(m.per_seam_generator(rate)@rho.reshape(-1,order='F')).reshape((m.n,m.n),order='F')
    expected=-1j*(m.hamiltonian@rho-rho@m.hamiltonian)
    for seam in range(3):expected+=rate*(m.commit(rho,seam)[0]-rho)
    assert np.allclose(derivative,expected)
    assert abs(np.trace(derivative))<1e-12
    zero=m.evolve_per_seam(rho,.37,0)
    step,_=m.step(rho,.37,0)
    assert np.allclose(zero,step)
    for value in (-1,float('nan'),float('inf')):
        with pytest.raises(ValueError):m.per_seam_generator(value)


def test_declared_run_distinguishes_instruments_and_regulators():
    result=run()
    assert all(result['checks'].values()) and all(result['mutations'].values())
    assert result['full_criterion_met'] is False


def test_rate_per_seam_is_unchanged_by_unoccupied_remote_component():
    # No two-load construction: the one prepared load occupies component0.
    local=LocalCoarseModular(Construction(2,((0,1),),.23),6)
    remote=LocalCoarseModular(Construction(4,((0,1),(2,3)),.23),6)
    assert local.n==remote.n==3
    assert np.array_equal(local.hamiltonian,remote.hamiltonian)
    rate=.29/(3*.37)
    assert np.allclose(local.per_seam_generator(rate).toarray(),remote.per_seam_generator(rate).toarray())
    a=local.evolve_per_seam(local.initial_density(),.37,rate)
    b=remote.evolve_per_seam(remote.initial_density(),.37,rate)
    assert np.allclose(a,b)
    assert np.allclose(remote.commit(remote.initial_density(),seam=1)[0],remote.initial_density())


def test_local_kraus_channel_preserves_external_reference_and_strict_schema():
    # Explicit extension C tensor I on an entangled state, without discarding R.
    m=LocalCoarseModular(Construction(2,((0,1),),.23),6)
    psi=np.zeros(m.n*2,complex);psi[0]=psi[3]=1/np.sqrt(2)
    rho=np.outer(psi,psi.conj());out=np.zeros_like(rho)
    for _,k in m.operators[0]:
        extended=np.kron(k,np.eye(2));out+=extended@rho@extended.T
    reference_before=np.einsum('aiaj->ij',rho.reshape(m.n,2,m.n,2))
    reference_after=np.einsum('aiaj->ij',out.reshape(m.n,2,m.n,2))
    assert np.allclose(reference_before,reference_after)
    assert np.isclose(np.trace(out),1) and np.linalg.eigvalsh(out).min()>-1e-12
    from oph_openwave.runner import validate_result
    validate_result('modular_native',run())
