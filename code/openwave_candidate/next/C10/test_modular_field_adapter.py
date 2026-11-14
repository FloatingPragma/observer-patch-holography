import numpy as np
from modular_field_adapter import run,ModularField,LocalCoarseModular,Construction


def test_actual_fullstate_reciprocal_field_and_native_work():
    r=run()
    for experiment in r['experiments']:
        assert all(experiment['checks'].values())
        state=np.array(experiment['final_state_real'])+1j*np.array(experiment['final_state_imag'])
        assert abs(np.trace(state)-1)<1e-12
        assert np.linalg.eigvalsh(state).min()>-1e-12
        energy=experiment['energy']
        assert abs(energy['source_after']-energy['source_before']+energy['field_gain'])<1e-12
        assert abs(energy['held_source_field_energy_gain'])<1e-12
        for law,readout in experiment['all_native_law_work'].items():
            assert abs(sum(x['probability'] for x in readout['complete_outcomes'])-1)<1e-12
            assert abs(sum(x['unnormalized_total_energy'] for x in readout['complete_outcomes'])-readout['after_energy'])<1e-12
    assert not r['full_criterion_met']


def test_zero_field_and_noncompressible_periodic_energy():
    m=LocalCoarseModular(Construction(3,((0,1),(1,2),(2,0)),.23),3)
    zero=ModularField(m,0,0);assert np.array_equal(zero.hamiltonian,m.hamiltonian)
    f=ModularField(m);n=m.n
    v=np.zeros((n,6))
    for i,(position,_) in enumerate(m.states):v[i,position]=1/np.sqrt(3)
    assert np.linalg.norm((np.eye(n)-v@v.T)@f.electric@v)>.1
    # Scalar compression cannot reproduce actual state evolution after adding HE.
    assert np.linalg.norm(f.hamiltonian@v-v@(v.T@f.hamiltonian@v))>.1
