import numpy as np
from oph_openwave.uniform_native import UniformFiberReduction,tiny_intertwiner_audit,Construction,LocalCoarseModular


def test_entire_local_instrument_intertwines_for_composite_modulus():
    for n in (3,6):
        result=tiny_intertwiner_audit(n)
        assert result['fiber_counts']==[n]*6
        assert result['H_intertwining_error']<1e-12
        assert result['every_coarse_Kraus_intertwining_error']<1e-12
        assert result['eight_step_lifted_density_error']<1e-12
        assert result['half_record_readout_error']<1e-12
        assert result['global_source_monitor_code_leakage']>.1
        assert result['definite_record_substitution_density_error']>.01


def test_external_quantum_reference_and_per_seam_generator():
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    full=LocalCoarseModular(c,6);reduced=UniformFiberReduction(c,6)
    v=np.zeros((full.n,reduced.n))
    for pos in range(reduced.n):
        indices=[i for i,s in enumerate(full.states) if s[0]==pos];v[indices,pos]=1/np.sqrt(6)
    psi=np.zeros(reduced.n*2,complex);psi[0]=psi[3]=1/np.sqrt(2)
    rho=np.outer(psi,psi.conj());isometry=np.kron(v,np.eye(2));lift=isometry@rho@isometry.T
    a=np.zeros_like(rho);b=np.zeros_like(lift)
    for (_,ka),(_,kb) in zip(reduced.kraus_for(0),full.operators[0]):
        ext_a=np.kron(ka,np.eye(2));ext_b=np.kron(kb,np.eye(2))
        a+=ext_a@rho@ext_a.T;b+=ext_b@lift@ext_b.T
    assert np.allclose(b,isometry@a@isometry.T)
    rho=reduced.initial_density();lift=v@rho@v.T;rate=.29/(3*.37)
    reference=(full.per_seam_generator(rate)@lift.reshape(-1,order='F')).reshape((full.n,full.n),order='F')
    assert np.allclose(reference,v@reduced.generator_action(rho,rate)@v.T)


def test_bridge_record_marginal_is_cut_charge_not_uniform_history():
    c=Construction(2,((0,1),),.23);model=UniformFiberReduction(c,6)
    assert model.fiber_cardinality==1
    rho,_=model.step(model.initial_density(),.37,.29)
    marginals=model.half_record_marginals(rho)
    assert all(r['type']=='bridge_cut_charge' for r in marginals)
    full=LocalCoarseModular(c,6);evolved,_=full.step(full.initial_density(),.37,.29)
    for h,row in enumerate(marginals):
        expected=[sum(evolved[i,i].real for i,s in enumerate(full.states) if s[1][h]==residue) for residue in range(6)]
        assert np.allclose(expected,row['residue_probabilities'])


def test_live_matrix_free_scheduler_adjoint_full_sector_and_schema():
    from oph_openwave.uniform_native import run
    from oph_openwave.runner import validate_result
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    model=UniformFiberReduction(c,6);full=LocalCoarseModular(c,6)
    op,trace=model.per_seam_operator(.2);rng=np.random.default_rng(65)
    x=rng.normal(size=model.n**2)+1j*rng.normal(size=model.n**2)
    y=rng.normal(size=model.n**2)+1j*rng.normal(size=model.n**2)
    assert abs(np.vdot(y,op@x)-np.vdot(op.rmatvec(y),x))<1e-11
    dense=np.column_stack([op@np.eye(model.n**2)[:,i] for i in range(model.n**2)])
    assert abs(np.trace(dense)-trace)<1e-11
    v=np.zeros((full.n,model.n))
    for i,(pos,_) in enumerate(full.states):v[i,pos]=1/np.sqrt(6)
    rho=model.initial_density();rate=.29/(3*.37)
    actual=model.evolve_per_seam(rho,.37,rate)
    expected=full.evolve_per_seam(v@rho@v.T,.37,rate)
    assert np.allclose(v@actual@v.T,expected,atol=1e-12)
    result=run();validate_result('uniform_native',result)
    assert 'complete_step_marginal_ledgers' not in result
