import numpy as np
from character_reference import run,blocks,FermionicSharedRecords,Construction


def test_all_characters_global_reference_complete_work():
    r=run();assert all(r['checks'].values())
    assert np.allclose(r['global_ground_character_weights'],[0,.5,.5],atol=1e-12)
    for out in r['outputs'].values():
        assert abs(sum(row['probability'] for row in out['outcomes'])-1)<1e-12
        assert abs(out['excitation_energy']-out['mean_instrument_work'])<1e-12
        assert abs(sum(row['unnormalized_final_energy'] for row in out['outcomes'])-out['energy_after'])<1e-12


def test_individual_record_measurement_mixes_characters():
    m=FermionicSharedRecords(Construction(3,((0,1),(1,2),(2,0)),.23),7,3)
    _,sectors=blocks(m);v=sectors[0][1];psi=v[:,0];rho=np.outer(psi,psi.conj());post=np.zeros_like(rho)
    for residue in range(3):
        mask=np.array([state[1][4]==residue for state in m.states]);post+=rho*mask[:,None]*mask[None,:]
    assert abs(np.trace(v@v.conj().T@post).real-1/3)<1e-12
